!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_rotframe.f90                                              #
!#                                                                           #
!# Copyright (C) 2010-2021                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
!#                                                                           #
!# This program is free software; you can redistribute it and/or modify      #
!# it under the terms of the GNU General Public License as published by      #
!# the Free Software Foundation; either version 2 of the License, or (at     #
!# your option) any later version.                                           #
!#                                                                           #
!# This program is distributed in the hope that it will be useful, but       #
!# WITHOUT ANY WARRANTY; without even the implied warranty of                #
!# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
!# NON INFRINGEMENT.  See the GNU General Public License for more            #
!# details.                                                                  #
!#                                                                           #
!# You should have received a copy of the GNU General Public License         #
!# along with this program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!----------------------------------------------------------------------------!
!> \addtogroup sources
!! - parameters of \link sources_rotframe_mod sources_rotframe \endlink as key-values
!!  \key{disable_centaccel,INTEGER,enable/disable (0|1) centrifugal acceleration,0)}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \attention Axis of rotation is assumed to point in z-direction; the angular
!!            is provided in Mesh%Omega
!!
!! \brief source terms module for inertial forces caused by a rotating grid
!----------------------------------------------------------------------------!
MODULE sources_rotframe_mod
  USE sources_c_accel_mod
  USE sources_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "inertial forces"
  TYPE, EXTENDS(sources_c_accel) :: sources_rotframe
    TYPE(marray_base), POINTER :: caccel => null(),  & !< rot. frame centrifugal accel
                                  vphi => null(), &    !< rot. frame azimuthal velocity
                                  twoOmega
    LOGICAL   :: disable_centaccel
  CONTAINS
    PROCEDURE :: InitSources_rotframe
    PROCEDURE :: InfoSources
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: Convert2RotatingFrame
    FINAL :: Finalize
  END TYPE

  PUBLIC :: &
       ! classes
       sources_rotframe
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of the rotating reference frame module
  SUBROUTINE InitSources_rotframe(this,Mesh,Physics,Fluxes,config,IO)
    USE physics_eulerisotherm_mod, ONLY: physics_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe)          :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),          POINTER :: config, IO
    !------------------------------------------------------------------------!
    TYPE(marray_base), POINTER :: caccel3D => null(), vphi3D => null(), Omez => null()
    INTEGER           :: stype,disable_centaccel
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitLogging(stype,source_name)

    CALL GetAttr(config,"disable_centaccel",disable_centaccel,0)
    IF (disable_centaccel.GT.0) THEN
      this%disable_centaccel = .TRUE.
    ELSE
      this%disable_centaccel = .FALSE.
    END IF

    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    SELECT TYPE(Physics)
    CLASS IS(physics_eulerisotherm)
      ! do nothing
    CLASS DEFAULT
       CALL this%Error("sources_rotframe::InitSources","physics not supported")
    END SELECT


    ALLOCATE(this%accel,caccel3D,vphi3D,Omez)
    this%accel = marray_base(Physics%VDIM)
    this%accel%data1d(:) = 0.

    ! compute centrifugal acceleration for rotation around vertical direction
    ! with center of rotation in r_0, i.e. -Omega**2 * ez x (ez x (r-r_0))
    caccel3D = marray_base(3)   ! 3D centrifugal acceleration
    vphi3D   = marray_base(3)   ! 3D azimuthal velocity caused by rotating frame
    Omez     = marray_base(3)   ! 3D angular velocity of rotating frame (= Omega*ez)
    ! set cartesian components of shifted center of rotation, i.e. r_0x, r_0y, r_0z,
    ! using this%cent as temporary storage
    caccel3D%data2d(:,1) = Mesh%rotcent(1)
    caccel3D%data2d(:,2) = Mesh%rotcent(2)
    caccel3D%data2d(:,3) = Mesh%rotcent(3)
    ! compute curvilinear components of shift vector
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,caccel3D%data4d,caccel3D%data4d)
    ! subtract the result from the position vector:
    ! this gives you the curvilinear components of all vectors pointing
    ! from the center of rotation to the bary center of any cell on the mesh,
    ! i.e. r-r_0
    caccel3D%data4d(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - caccel3D%data4d(:,:,:,:)
    ! set cartesian components of Omega*ez within computational domain
    Omez%data2d(:,1:2) = 0.0  ! no x and y component
    WHERE (Mesh%without_ghost_zones%mask1d(:))
      Omez%data2d(:,3) = Mesh%Omega
    ELSEWHERE
      Omez%data2d(:,3) = 0.0
    END WHERE
    ! compute curvinlinear components of Omega*ez with respect to the given geometry of the mesh
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,Omez%data4d,Omez%data4d)
    ! compute curvilinear components of the local azimuthal velocity caused by the rotating frame
    ! and the centrifual acceleration
    vphi3D = Omez.x.caccel3D ! = Omega*ez x (r-r0)
    caccel3D = vphi3D.x.Omez ! = -Omega**2 * ez x (ez x (r-r0)) = (Omega*ez x (r-r0)) x Omega*ez

    ! set the projected components of the centrifugal acceleration
    IF (Physics%VDIM.LT.3) THEN
      ALLOCATE(this%caccel,this%vphi)
      this%caccel = marray_base(Physics%VDIM) ! dimensionality of velocity vector
      this%vphi   = marray_base(Physics%VDIM)
      IF(Physics%VDIM.EQ.2) THEN
        ALLOCATE(this%twoOmega)
        this%twoOmega = marray_base(1)  ! = 2*Omega*(ez*e_perp)) = 2*projection of Omega onto 3rd suppressed dimension
      END IF
    END IF
    SELECT CASE(Mesh%VECTOR_COMPONENTS)
    CASE(VECTOR_X) ! 1D momentum in x-direction (1st curvilinear direction)
      ! vy = vz = my = mz = 0
      this%caccel%data2d(:,1) = caccel3D%data2d(:,1)
      this%vphi%data2d(:,1)   = vphi3D%data2d(:,1)
    CASE(VECTOR_Y) ! 1D momentum in y-direction (2nd curvilinear direction)
      ! vx = vz = mx = mz = 0
      this%caccel%data2d(:,1) = caccel3D%data2d(:,2)
      this%vphi%data2d(:,1)   = vphi3D%data2d(:,2)
    CASE(VECTOR_Z) ! 1D momentum in z-direction (3rd curvilinear direction)
      ! vx = vy = mx = my = 0
      this%caccel%data2d(:,1) = caccel3D%data2d(:,3)
      this%vphi%data2d(:,1)   = vphi3D%data2d(:,3)
    CASE(IOR(VECTOR_X,VECTOR_Y)) ! 2D momentum in x-y-plane
      ! vz = mz = 0
      this%caccel%data2d(:,1:2) = caccel3D%data2d(:,1:2)
      this%vphi%data2d(:,1:2)   = vphi3D%data2d(:,1:2)
      this%twoOmega%data1d(:)   = 2*Omez%data2d(:,3)
    CASE(IOR(VECTOR_X,VECTOR_Z)) ! 2D momentum in x-z-plane
      ! vy = my = 0
      this%caccel%data2d(:,1) = caccel3D%data2d(:,1)
      this%vphi%data2d(:,1)   = vphi3D%data2d(:,1)
      this%caccel%data2d(:,2) = caccel3D%data2d(:,3)
      this%vphi%data2d(:,2)   = vphi3D%data2d(:,3)
      ! the minus is because e1 and e3 become the new 2D
      ! curvilinear basis vectors and we must flip the direction of e2
      ! in order to retain a right handed triad and hence
      ! we must project Omez onto -e2 instead of e2
      this%twoOmega%data1d(:)   = -2*Omez%data2d(:,2)
    CASE(IOR(VECTOR_Y,VECTOR_Z)) ! 2D momentum in y-z-plane
      ! vx = mx = 0
      this%caccel%data2d(:,1:2) = caccel3D%data2d(:,2:3)
      this%vphi%data2d(:,1:2)   = vphi3D%data2d(:,2:3)
      this%twoOmega%data1d(:)   = 2*Omez%data2d(:,1)
    CASE DEFAULT
      ! full 3D case: caccel = caccel3D
      this%caccel => caccel3D
      this%vphi => vphi3D
      this%twoOmega => Omez
      this%twoOmega%data1d(:) = 2.0*this%twoOmega%data1d(:)
    END SELECT
    IF(Physics%VDIM.LT.3) DEALLOCATE(caccel3D,vphi3D,Omez)
    IF(this%disable_centaccel) DEALLOCATE(this%caccel)
  END SUBROUTINE InitSources_rotframe


  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(IN) :: this
    CLASS(mesh_base),        INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: omega_str
    !------------------------------------------------------------------------!
    WRITE (omega_str,'(ES9.2)') Mesh%OMEGA
    CALL this%Info("            angular velocity:  " // TRIM(omega_str))
    IF (this%disable_centaccel) THEN
      CALL this%Info("            centrifugal accel: " // "disabled")
    END IF
  END SUBROUTINE InfoSources


  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    USE physics_eulerisotherm_mod, ONLY : statevector_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    CLASS(physics_base),     INTENT(INOUT) :: Physics
    CLASS(fluxes_base),      INTENT(IN)    :: Fluxes
    CLASS(sources_base),     INTENT(INOUT) :: Sources
    REAL,                    INTENT(IN)    :: time, dt
    CLASS(marray_compound),  INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_eulerisotherm)
      SELECT CASE(Physics%VDIM)
      CASE(1) ! 1D acceleration
        ! no coriolis acceleration, because it is always perpendicular to the
        ! velocity vector and hence causes no acceleration along this direction
        IF (.NOT.this%disable_centaccel) THEN
          this%accel = this%caccel
        END IF
        ! otherwise do nothing, this%accel should be zero (see initialization)
      CASE(2) ! 2D acceleration
        IF (this%disable_centaccel) THEN
          this%accel%data2d(:,1) = this%twoOmega%data1d(:)*p%velocity%data2d(:,2)
          this%accel%data2d(:,2) = -this%twoOmega%data1d(:)*p%velocity%data2d(:,1)
        ELSE
          this%accel%data2d(:,1) = this%caccel%data2d(:,1) + this%twoOmega%data1d(:)*p%velocity%data2d(:,2)
          this%accel%data2d(:,2) = this%caccel%data2d(:,2) - this%twoOmega%data1d(:)*p%velocity%data2d(:,1)
        END IF
      CASE(3) ! 3D acceleration
        IF (this%disable_centaccel) THEN
          this%accel = p%velocity.x.this%twoOmega ! = -2*Omega x v
        ELSE
          this%accel = this%caccel + (p%velocity.x.this%twoOmega) ! = caccel - 2*Omega x v
        END IF
      CASE DEFAULT
        ! this should not happen
      END SELECT
    CLASS DEFAULT
      ! do nothing
    END SELECT

    ! inertial forces source terms
    CALL Physics%ExternalSources(this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_single

  SUBROUTINE Convert2RotatingFrame(this,Mesh,Physics,pvar)
    USE physics_eulerisotherm_mod, ONLY : statevector_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(IN)    :: this
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    CLASS(physics_base),     INTENT(IN)    :: Physics
    CLASS(marray_compound),  INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    ! Convert velocities to the rotating frame
    SELECT TYPE (p => pvar)
    TYPE IS(statevector_eulerisotherm)
      p%velocity%data1d = p%velocity%data1d - this%vphi%data1d
    CLASS DEFAULT
      ! nothing happens
    END SELECT
  END SUBROUTINE Convert2RotatingFrame

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_rotframe), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%vphi)
    IF (ASSOCIATED(this%caccel)) DEALLOCATE(this%caccel)
    IF (ASSOCIATED(this%twoOmega)) DEALLOCATE(this%twoOmega)
    NULLIFY(this%vphi,this%caccel,this%twoOmega)
    ! deallocation of this%accel is done in inherited destructor
    ! which is called automatically
  END SUBROUTINE Finalize

END MODULE sources_rotframe_mod
