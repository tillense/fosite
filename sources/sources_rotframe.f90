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
!!  \key{gparam,REAL,geometry parameter (needed for bianglspherical
!!       geometry - corresponds to the radius)}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \attention This module works currently only in 2D. Mesh%OMEGA is only
!!            one value and assumened to point in z-direction.
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
    TYPE(marray_base), ALLOCATABLE :: cent, &     !< rot. frame centrifugal
                                      centproj, & !< rot.frame centr.3d->2d
                                      cos1, &
                                      sin1
    INTEGER           :: issphere
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
    USE physics_euler_mod, ONLY: physics_euler
    USE physics_eulerisotherm_mod, ONLY: physics_eulerisotherm
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    USE geometry_spherical_mod, ONLY: geometry_spherical
    USE geometry_spherical_planet_mod, ONLY: geometry_spherical_planet
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe)          :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),          POINTER :: config, IO
    !------------------------------------------------------------------------!
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitLogging(stype,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    SELECT TYPE(Physics)
    TYPE IS(physics_euler)
      ! do nothing
    TYPE IS(physics_eulerisotherm)
      ! do nothing
    CLASS DEFAULT
       CALL this%Error("sources_rotframe::InitSources","physics not supported")
    END SELECT

    IF (Mesh%FARGO.GT.0) &
      CALL this%Error("sources_rotframe::InitSources", &
         "rotating frame with FARGO seems broken, don't use it together")

!    IF (Physics%VDIM.NE.2) &
!      CALL this%Error("sources_rotframe::InitSources", &
!         "Only 2D simulations working with rotating frame at the moment." )

    ALLOCATE(this%accel,this%cent,this%centproj,this%cos1,this%sin1)
    this%accel = marray_base(Physics%VDIM)
    this%accel%data1d(:) = 0.
    this%cent = marray_base(3)
    this%cent%data1d(:) = 0.
    CALL GetAttr(config, "gparam", this%gparam, 1.0)
    CALL GetAttr(config, "issphere", this%issphere, 0)

    ! define position vectors
    ! for planetery atmosphere no shift of axis possible
    ! (for a planet reasonable)
    SELECT TYPE(geo => Mesh%geometry)
    TYPE IS(geometry_spherical_planet)
      IF (this%issphere.EQ.1) THEN
        this%centproj = marray_base(3)
        this%centproj%data1d(:) = 0.
        this%cos1 = marray_base(3)
        this%cos1%data1d(:) = 0.
 
        this%centproj%data4d(:,:,:,1) = this%gparam*SIN(Mesh%bcenter(:,:,:,1))*&
                                COS(Mesh%bcenter(:,:,:,1))
        ! for better performance
        this%cos1%data4d(:,:,:,1)  = COS(Mesh%bcenter(:,:,:,1))
        this%cos1%data4d(:,:,:,2)  = COS(Mesh%bcenter(:,:,:,2))
      END IF
    TYPE IS(geometry_cylindrical)
      IF (ABS(Mesh%rotcent(1)).LE.TINY(Mesh%rotcent(1)) &
       .AND.ABS(Mesh%rotcent(2)).LE.TINY(Mesh%rotcent(2))) THEN
        ! no shift of point mass: set position vector to Mesh defaults
        this%cent%data4d(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:)
      ELSE
        ! shifted center of rotation:
        ! compute curvilinear components of shift vector
        this%cent%data4d(:,:,:,1) = Mesh%rotcent(1)
        this%cent%data4d(:,:,:,2) = Mesh%rotcent(2)
        this%cent%data4d(:,:,:,3) = Mesh%rotcent(3)
        CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,this%cent%data4d,this%cent%data4d)
        ! subtract the result from the position vector:
        ! this gives you the curvilinear components of all vectors pointing
        ! from the center of rotation to the bary center of any cell on the mesh
        this%cent%data4d(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - this%cent%data4d(:,:,:,:)
      END IF
    TYPE IS(geometry_spherical)
        this%cos1 = marray_base()
        this%cos1%data3d(:,:,:) = COS(Mesh%bcenter(:,:,:,2))
        this%sin1 = marray_base()
        this%sin1%data3d(:,:,:) = SIN(Mesh%bcenter(:,:,:,2))
        this%centproj = marray_base()
        this%centproj%data4d(:,:,:,1) = Mesh%radius%bcenter(:,:,:)*&
                                this%Sin1%data3d(:,:,:)
      CLASS DEFAULT
        CALL this%Error('Sources:Sources_rotframe','Geometry not supported at the moment')    
    END SELECT
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
  END SUBROUTINE InfoSources


  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    USE geometry_spherical_mod, ONLY: geometry_spherical
    USE geometry_spherical_planet_mod, ONLY: geometry_spherical_planet
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
    SELECT TYPE(geo => Mesh%Geometry)
    TYPE IS(geometry_spherical_planet)
      ! Two cases due to different angles between angular velocity to mesh
      ! 1. only a projected part plays role for bianglespherical geometry
      IF (this%issphere.EQ.1) THEN
!NEC$ outerloop_unroll(8)
        DO k=Mesh%KMIN,Mesh%KMAX
         DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ ivdep
            DO i=Mesh%IMIN,Mesh%IMAX
              this%accel%data4d(i,j,k,1) = Mesh%OMEGA*(Mesh%OMEGA*&
                     this%centproj%data4d(i,j,k,1) + 2.0*this%cos1%data4d(i,j,k,1)*&
                     pvar%data4d(i,j,k,Physics%YVELOCITY))
              this%accel%data4d(i,j,k,2) = -Mesh%OMEGA*2.0*this%cos1%data4d(i,j,k,1)*&
                     pvar%data4d(i,j,k,Physics%XVELOCITY)
            END DO
          END DO
        END DO
      END IF
    TYPE IS(geometry_cylindrical)
    ! 2. OMEGA is directing in z direction
!NEC$ outerloop_unroll(8)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ ivdep
            DO i=Mesh%IMIN,Mesh%IMAX
              ! components of centrifugal and coriolis acceleration
              this%accel%data4d(i,j,k,1) = Mesh%OMEGA*(Mesh%OMEGA*this%cent%data4d(i,j,k,1) &
                   + 2.0*pvar%data4d(i,j,k,Physics%YVELOCITY))
              this%accel%data4d(i,j,k,2) = -Mesh%OMEGA*2.0*pvar%data4d(i,j,k,Physics%XVELOCITY)
              this%accel%data4d(i,j,k,3) = 0.0
            END DO
          END DO
        END DO
    TYPE IS(geometry_spherical)
      ! 2. OMEGA is directing in z direction
!NEC$ outerloop_unroll(8)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ ivdep
            DO i=Mesh%IMIN,Mesh%IMAX
              this%accel%data4d(i,j,k,1) = Mesh%Omega*(Mesh%Omega* &
                     this%sin1%data3d(i,j,k)*this%centproj%data4d(i,j,k,1) &
                     + 2.0*pvar%data4d(i,j,k,Physics%ZVELOCITY)*this%Sin1%data3d(i,j,k))
              this%accel%data4d(i,j,k,2) = Mesh%Omega*(Mesh%Omega* &
                     this%cos1%data3d(i,j,k)*this%centproj%data4d(i,j,k,1) &
                     + 2.0*pvar%data4d(i,j,k,Physics%ZVELOCITY)*this%Cos1%data3d(i,j,k))
              this%accel%data4d(i,j,k,3) = -2.0*Mesh%Omega*( &
                          pvar%data4d(i,j,k,Physics%XVELOCITY)*this%Sin1%data3d(i,j,k) &
                        + pvar%data4d(i,j,k,Physics%YVELOCITY)*this%Cos1%data3d(i,j,k))
            END DO
          END DO
        END DO
    END SELECT

    ! inertial forces source terms
    CALL Physics%ExternalSources(this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_single

  SUBROUTINE Convert2RotatingFrame(this,Mesh,Physics,pvar)
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    USE geometry_spherical_mod, ONLY: geometry_spherical
    USE geometry_spherical_planet_mod, ONLY: geometry_spherical_planet
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(IN)    :: this
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    CLASS(physics_base),     INTENT(IN)    :: Physics
    CLASS(marray_compound),  INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    ! Convert velocities to the rotating frame
    SELECT TYPE(geo => Mesh%Geometry)
    TYPE IS(geometry_spherical_planet)
      IF (this%issphere.EQ.1) THEN
        ! no change in pvar(:,:,:,Physics%XVELOCITY)
        pvar%data4d(:,:,:,Physics%YVELOCITY) = pvar%data4d(:,:,:,Physics%YVELOCITY) - &
                              Mesh%OMEGA*SIN(Mesh%bcenter(:,:,:,1))*this%gparam
      END IF
    TYPE IS(geometry_cylindrical)
      pvar%data4d(:,:,:,Physics%XVELOCITY) = pvar%data4d(:,:,:,Physics%XVELOCITY) + &
                                      Mesh%OMEGA * this%cent%data4d(:,:,:,2)
      pvar%data4d(:,:,:,Physics%YVELOCITY) = pvar%data4d(:,:,:,Physics%YVELOCITY) - &
                                      Mesh%OMEGA * this%cent%data4d(:,:,:,1)
      pvar%data4d(:,:,:,Physics%ZVELOCITY) = pvar%data4d(:,:,:,Physics%ZVELOCITY)
    TYPE IS(geometry_spherical)
      pvar%data4d(:,:,:,Physics%XVELOCITY) = pvar%data4d(:,:,:,Physics%XVELOCITY)
      pvar%data4d(:,:,:,Physics%YVELOCITY) = pvar%data4d(:,:,:,Physics%YVELOCITY)
      pvar%data4d(:,:,:,Physics%ZVELOCITY) = pvar%data4d(:,:,:,Physics%ZVELOCITY) - &
                                      Mesh%OMEGA * this%centproj%data4d(:,:,:,1)
    CLASS DEFAULT
      ! this should not happen (see InitSources_rotframe)
    END SELECT
  END SUBROUTINE Convert2RotatingFrame

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_rotframe), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%cent,this%centproj,this%cos1,this%sin1)
    ! deallocation of this%accel is done in inherited destructor
    ! which is called automatically
  END SUBROUTINE Finalize

END MODULE sources_rotframe_mod
