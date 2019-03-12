!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_shearbox.f90                                              #
!#                                                                           #
!# Copyright (C) 2010-2018                                                   #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!! - parameters of \link sources_shearbox_mod sources_shearbox \endlink as key-values
!!
!! \brief source term module for shearing forces in a shearingsheet
!!
!! \todo This module can only be used with fixed 2D plane in x-y direction.
!!       It should be no problem to make it also run in 3D.
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!!
!! \brief Source terms module for fictious forces in a shearingsheet.
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_shearbox_mod
  USE sources_c_accel_mod
  USE sources_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
    CHARACTER(LEN=32) :: source_name = "forces in shearing-box"
  !--------------------------------------------------------------------------!

  TYPE, EXTENDS(sources_c_accel) :: sources_shearbox
    REAL      :: SIGN1, SIGN2
    INTEGER   :: VEL1, VEL2
    INTEGER   :: I1, I2
    INTEGER   :: MOMENTUM1, MOMENTUM2
  CONTAINS
    PROCEDURE :: InitSources_shearbox
    PROCEDURE :: InfoSources
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: Finalize
  END TYPE

  PUBLIC :: &
       ! classes
       sources_shearbox

CONTAINS

  !> \public Constructor of sources shearbox module.
  !!
  !! This subroutine reads the necessary config data for the source module
  !! for a shearingsheet. It initializes the sources type and the
  !! accelerations array.
  SUBROUTINE InitSources_shearbox(this,Mesh,Physics,Fluxes,config,IO)
    USE physics_euler_mod, ONLY: physics_euler
    USE physics_eulerisotherm_mod, ONLY: physics_eulerisotherm
    USE geometry_cartesian_mod, ONLY: geometry_cartesian
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox)          :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),          POINTER :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                          :: stype
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"stype",stype)
    CALL this%InitLogging(stype,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    SELECT TYPE(Physics)
    TYPE IS(physics_euler)
      ! do nothing
    TYPE IS(physics_eulerisotherm)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("InitSources_shearbox","physics not supported")
    END SELECT

    SELECT TYPE(geo=>Mesh%Geometry)
    TYPE IS(geometry_cartesian)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("InitSources_shearbox","mesh not supported")
    END SELECT

    this%accel = marray_base(Physics%VDIM)
    this%accel%data1d(:) = 0.

    IF(Mesh%shear_dir.EQ.2) THEN
      this%Vel1=Physics%YVELOCITY
      this%Vel2=Physics%XVELOCITY
      this%SIGN1 = 1.0
      this%SIGN2 = -1.0
      this%I1 = 1
      this%I2 = 2
      this%MOMENTUM1 = Physics%XMOMENTUM
      this%MOMENTUM2 = Physics%YMOMENTUM
    ELSE IF(Mesh%shear_dir.EQ.1) THEN
      this%Vel1=Physics%XVELOCITY
      this%Vel2=Physics%YVELOCITY
      this%SIGN1 = -1.0
      this%SIGN2 = 1.0
      this%I1 = 2
      this%I2 = 1
      this%MOMENTUM1 = Physics%YMOMENTUM
      this%MOMENTUM2 = Physics%XMOMENTUM
    END IF

    ! set vertical gravitational acceleration for 3D shearing
    ! box simulations using bary center cartesian z-coordinate
    IF (Mesh%KNUM.GT.1) THEN
      this%accel%data2d(:,3) = -Mesh%OMEGA**2 * Mesh%cart%data3d(:,2,3)
    END IF
  END SUBROUTINE InitSources_shearbox

  !> \public Write shearbox fictious forces parameters to screen.
  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox), INTENT(IN) :: this
    CLASS(mesh_base),        INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)       :: omega_str,q_str
    !------------------------------------------------------------------------!
    WRITE (omega_str,'(ES8.2)') Mesh%OMEGA
    WRITE (q_str,'(ES8.2)')     Mesh%Q
    CALL this%Info("            angular velocity:  " // TRIM(omega_str))
    CALL this%Info("            shearing parameter:" // TRIM(q_str))
  END SUBROUTINE InfoSources

  !> \public Compute fictious forces within a shearingsheet.
  !!
  !! Computation is done via
  !! \f[
  !!    \mathbf{a} = - 2 \mathbf{\Omega} \times \mathbf{v} + 2 q \Omega^2 x
  !!                    \mathbf{\hat{e}_x} - \Omega^2 z \mathbf{\hat{e}_z}
  !! \f]
  !! See for example \cite gammie2001 or \cite hawley1995 .
  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    CLASS(physics_base),     INTENT(INOUT) :: Physics
    CLASS(fluxes_base),      INTENT(IN)    :: Fluxes
    CLASS(sources_base),     INTENT(INOUT) :: Sources
    REAL,                    INTENT(IN)    :: time, dt
    CLASS(marray_compound),INTENT(INOUT)   :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    SELECT CASE(Mesh%FARGO)
    CASE(0)
      ! fargo transport disabled
!NEC$ IVDEP
      this%accel%data4d(:,:,:,this%I1) = 2*Mesh%OMEGA &
              * (Mesh%Q*Mesh%OMEGA*Mesh%bcenter(:,:,:,this%I1) &
               + this%SIGN1*pvar%data4d(:,:,:,this%VEL1))
      this%accel%data4d(:,:,:,this%I2) = 2*Mesh%OMEGA*this%SIGN2 &
              * pvar%data4d(:,:,:,this%VEL2)
      ! shearingsheet inertial forces source terms
      CALL Physics%ExternalSources(this%accel,pvar,cvar,sterm)
    CASE(3)
      ! fargo transport type 3 enabled
      sterm%data2d(:,Physics%DENSITY) = 0.0
!NEC$ IVDEP
      sterm%data2d(:,this%MOMENTUM1) = pvar%data2d(:,Physics%DENSITY) &
        *Mesh%OMEGA*2.0*this%SIGN1*pvar%data2d(:,this%VEL1)
      sterm%data2d(:,this%MOMENTUM2) = pvar%data2d(:,Physics%DENSITY) &
        *Mesh%OMEGA*(2.0-Mesh%Q)*this%SIGN2*pvar%data2d(:,this%VEL2)
      IF (Mesh%KNUM.GT.1) &
        sterm%data2d(:,Physics%ZMOMENTUM) = pvar%data2d(:,Physics%DENSITY) &
          *this%accel%data2d(:,3)
      IF (Physics%PRESSURE .GT. 0) THEN
        sterm%data2d(:,Physics%ENERGY) = &
             this%SIGN1*pvar%data2d(:,Physics%DENSITY)*Mesh%Q*Mesh%OMEGA* &
             pvar%data2d(:,this%VEL2)*pvar%data2d(:,this%VEL1)
        IF (Mesh%KNUM.GT.1) &
          sterm%data2d(:,Physics%ENERGY) = sterm%data2d(:,Physics%ENERGY) &
            + cvar%data2d(:,Physics%ZMOMENTUM)*this%accel%data2d(:,3)
      END IF
    CASE DEFAULT
      ! other fargo transport schemes are not supported
      CALL this%Error("sources_shearbox::ExternalSources_single", &
                      "currently only Fargo transport type 3 is supported")
    END SELECT
  END SUBROUTINE ExternalSources_single

  !> \public Closes the shearingsheet source term.
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%accel%Destroy()
  END SUBROUTINE Finalize

END MODULE sources_shearbox_mod
