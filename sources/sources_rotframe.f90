!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_rotframe.f03                                              #
!#                                                                           #
!# Copyright (C) 2010-2018                                                   #
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
!! - parameters of \link sources_shearing \endlink as key-values
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
  USE marray_compound_mod
  USE common_dict
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "inertial forces"
  TYPE, EXTENDS(sources_c_accel) :: sources_rotframe
    REAL, DIMENSION(:,:,:,:), POINTER :: cent         !< rot. frame centrifugal
    REAL, DIMENSION(:,:,:,:), POINTER :: centproj     !< rot.frame centr.3d->2d
    REAL, DIMENSION(:,:,:,:), POINTER :: cos1
  CONTAINS
    PROCEDURE :: InitSources_rotframe
    PROCEDURE :: InfoSources
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: Convert2RotatingFrame
    PROCEDURE :: Finalize
  END TYPE

  PUBLIC :: &
       ! classes
       sources_rotframe
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of the rotating reference frame module
  SUBROUTINE InitSources_rotframe(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe)          :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),          POINTER :: config, IO
    !------------------------------------------------------------------------!
    INTEGER           :: stype
    INTEGER           :: err
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitLogging(stype,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    SELECT CASE(Physics%GetType())
    CASE(EULER2D,EULER2D_ISOTHERM,EULER2D_IAMT,EULER3D,EULER3D_ISOTHERM)
       ! do nothing
    CASE DEFAULT
       CALL this%Error("ExternalSources_rotframe","physics not supported")
    END SELECT

    ALLOCATE( &
         this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VDIM), &
         this%cent(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
         this%centproj(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
         this%cos1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
         STAT = err)
    IF (err.NE.0) CALL this%Error("InitSources_rotframe", "Unable allocate memory!")

    CALL GetAttr(config, "gparam", this%gparam, 1.0)

    ! define position vectors
    ! for bianglespherical no shift of axis possible
    ! (for a planet reasonable)
    IF ((Mesh%Geometry%GetType()).EQ.BIANGLESPHERICAL) THEN

      this%centproj(:,:,:,1) = this%gparam*SIN(Mesh%bcenter(:,:,:,1))*&
                              COS(Mesh%bcenter(:,:,:,1))
      ! for better performance
      this%cos1(:,:,:,1)  = COS(Mesh%bcenter(:,:,:,1))
      this%cos1(:,:,:,2)  = COS(Mesh%bcenter(:,:,:,2))
    ELSE
      IF (ABS(Mesh%rotcent(1)).LE.TINY(Mesh%rotcent(1)) &
       .AND.ABS(Mesh%rotcent(2)).LE.TINY(Mesh%rotcent(2))) THEN
        ! no shift of point mass: set position vector to Mesh defaults
        this%cent(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:)
      ELSE
        ! shifted center of rotation:
        ! compute curvilinear components of shift vector
        this%cent(:,:,:,1) = Mesh%rotcent(1)
        this%cent(:,:,:,2) = Mesh%rotcent(2)
        this%cent(:,:,:,3) = Mesh%rotcent(3)
        CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,this%cent,this%cent)
        ! subtract the result from the position vector:
        ! this gives you the curvilinear components of all vectors pointing
        ! from the center of rotation to the bary center of any cell on the mesh
        this%cent(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - this%cent(:,:,:,:)
      END IF
    END IF

    ! reset acceleration term
    this%accel(:,:,:,:) = 0.0
  END SUBROUTINE InitSources_rotframe


  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(IN) :: this
    CLASS(mesh_base),        INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: omega_str
    !------------------------------------------------------------------------!
    WRITE (omega_str,'(ES8.2)') Mesh%OMEGA
    CALL this%Info("            angular velocity:  " // TRIM(omega_str))
  END SUBROUTINE InfoSources


  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
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
    ! Two cases due to different angles between angular velocity to mesh
    ! 1. only a projected part plays role for bianglespherical geometry
    IF ((Mesh%Geometry%GetType()).EQ.BIANGLESPHERICAL) THEN
!CDIR OUTERUNROLL=8
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
            this%accel(i,j,k,1) = Mesh%OMEGA*(Mesh%OMEGA*&
                   this%centproj(i,j,k,1) + 2.0*this%cos1(i,j,k,1)*&
                   pvar%data4d(i,j,k,Physics%YVELOCITY))
            this%accel(i,j,k,2) = -Mesh%OMEGA*2.0*this%cos1(i,j,k,1)*&
                   pvar%data4d(i,j,k,Physics%XVELOCITY)
          END DO
        END DO
      END DO
    ! 2. OMEGA is always perpendicular to other curvilinear coordinates
    ELSE
!CDIR OUTERUNROLL=8
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
  !CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
            ! components of centrifugal and coriolis acceleration
            this%accel(i,j,k,1) = Mesh%OMEGA*(Mesh%OMEGA*this%cent(i,j,k,1) &
                 + 2.0*pvar%data4d(i,j,k,Physics%YVELOCITY))
            this%accel(i,j,k,2) = Mesh%OMEGA*(Mesh%OMEGA*this%cent(i,j,k,2) &
                 - 2.0*pvar%data4d(i,j,k,Physics%XVELOCITY))
          END DO
        END DO
      END DO
    END IF

    ! inertial forces source terms
    CALL Physics%ExternalSources(Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_single

  SUBROUTINE Convert2RotatingFrame(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(IN)    :: this
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    CLASS(physics_base),     INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                             INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    ! Convert velocities to the rotating frame
    IF ((Mesh%Geometry%GetType()).EQ.BIANGLESPHERICAL) THEN
      pvar(:,:,:,Physics%XVELOCITY) = pvar(:,:,:,Physics%XVELOCITY)
      pvar(:,:,:,Physics%YVELOCITY) = pvar(:,:,:,Physics%YVELOCITY) - &
                            Mesh%OMEGA*SIN(Mesh%bcenter(:,:,:,1))*this%gparam
    ELSE
      pvar(:,:,:,Physics%XVELOCITY) = pvar(:,:,:,Physics%XVELOCITY) + &
                                      Mesh%OMEGA * this%cent(:,:,:,2)
      pvar(:,:,:,Physics%YVELOCITY) = pvar(:,:,:,Physics%YVELOCITY) - &
                                      Mesh%OMEGA * this%cent(:,:,:,1)
    END IF
  END SUBROUTINE Convert2RotatingFrame

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_rotframe), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel,this%cent,this%cos1,this%centproj)
    CALL this%next%Finalize()
  END SUBROUTINE Finalize

END MODULE sources_rotframe_mod
