!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_c_accel.f90                                               #
!#                                                                           #
!# Copyright (C) 2009-2024                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <tillense@astrophysik.uni-kiel.de>                       #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief source terms module for constant acceleration
!!
!! \extends sources_common
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_c_accel_mod
  USE sources_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
    CHARACTER(LEN=32) :: source_name = "constant acceleration"

  TYPE, EXTENDS(sources_base) :: sources_c_accel
    TYPE(marray_base), ALLOCATABLE :: accel          !< acceleration      !
  CONTAINS
    PROCEDURE :: InitSources_c_accel
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcTimestep_single
    PROCEDURE :: InfoSources
    FINAL :: Finalize
  END TYPE

  PUBLIC :: &
       ! classes
       sources_c_accel

CONTAINS

  SUBROUTINE InitSources_c_accel(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_c_accel), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(fluxes_base),     INTENT(IN)    :: Fluxes
    CLASS(physics_base),    INTENT(IN)    :: Physics
    TYPE(Dict_TYP),           POINTER     :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: k,stype
    REAL    :: accel(3)
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"stype",stype)
    CALL this%InitLogging(stype,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    ALLOCATE(this%accel)
    this%accel = marray_base(Physics%VDIM)
    this%accel%data1d(:) = 0.0

    ! initialize constant acceleration
    CALL GetAttr(config, "xaccel", accel(1), 0.0)
    CALL GetAttr(config, "yaccel", accel(2), 0.0)
    CALL GetAttr(config, "zaccel", accel(3), 0.0)
    DO k=1,Physics%VDIM
      this%accel%data2d(:,k) = accel(k)
    END DO

  END SUBROUTINE InitSources_c_accel

  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_c_accel),INTENT(INOUT):: this
    CLASS(mesh_base),INTENT(IN)         :: Mesh
    CLASS(physics_base),INTENT(INOUT)   :: Physics
    CLASS(fluxes_base),INTENT(IN)       :: Fluxes
    CLASS(sources_base), INTENT(INOUT)  :: Sources
    REAL,INTENT(IN)                     :: time, dt
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    ! compute source terms due to constant acceleration
    CALL Physics%ExternalSources(this%accel,pvar,cvar,sterm)
  END SUBROUTINE

  SUBROUTINE CalcTimestep_single(this,Mesh,Physics,Fluxes,pvar,cvar,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_c_accel), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(INOUT) :: Physics
    CLASS(fluxes_base),     INTENT(IN)    :: Fluxes
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    REAL,                   INTENT(IN)    :: time
    REAL,                   INTENT(OUT)   :: dt
    !------------------------------------------------------------------------!
    dt = HUGE(dt)
  END SUBROUTINE CalcTimestep_single

  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_c_accel), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
  END SUBROUTINE InfoSources

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_c_accel), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (ALLOCATED(this%accel)) DEALLOCATE(this%accel)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE sources_c_accel_mod
