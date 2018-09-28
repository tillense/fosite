!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_c_accel.f03                                               #
!#                                                                           #
!# Copyright (C) 2009,2011,2018                                              #
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
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
    CHARACTER(LEN=32) :: source_name = "constant acceleration"

  TYPE, EXTENDS(sources_base) :: sources_c_accel
    REAL, DIMENSION(:,:,:,:), POINTER :: accel          !< acceleration      !
 CONTAINS
    PROCEDURE :: InitSources_c_accel
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcTimestep_single
    PROCEDURE :: InfoSources
    PROCEDURE :: Finalize
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
    REAL    :: xaccel, yaccel, zaccel
    INTEGER :: err
    !------------------------------------------------------------------------!
    CALL this%InitLogging(C_ACCEL,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    ALLOCATE(&
      this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,&
                 Mesh%KGMIN:Mesh%KGMAX,Physics%DIM), &
         STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitSources_c_accel","memory allocation failed")

    ! initialize constant acceleration
    CALL GetAttr(config, "xaccel", xaccel, 0.0)
    this%accel(:,:,:,1) = xaccel

    CALL GetAttr(config, "yaccel", yaccel, 0.0)
    this%accel(:,:,:,2) = yaccel

    IF (Physics%DIM .GE. 3) THEN
      CALL GetAttr(config, "zaccel", zaccel, 0.0)
      this%accel(:,:,:,3) = zaccel
    END IF
  END SUBROUTINE InitSources_c_accel

  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,time,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_c_accel), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(INOUT) :: Physics
    CLASS(fluxes_base),     INTENT(IN)    :: Fluxes
    REAL,                   INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(IN)    :: cvar,pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(OUT)   :: sterm
    !------------------------------------------------------------------------!
    ! compute source terms due to constant acceleration
    CALL Physics%ExternalSources(Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_single

  SUBROUTINE CalcTimestep_single(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_c_accel), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(INOUT) :: Physics
    CLASS(fluxes_base),     INTENT(IN)    :: Fluxes
    REAL,                   INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(IN)    :: pvar,cvar
    REAL,                   INTENT(OUT)   :: dt
    !------------------------------------------------------------------------!
    !\todo How should be handled routines where no CalcTimestep normally
    !!     included. Is this handling ok with just setting the timestep huge?
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
    CLASS(sources_c_accel), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel)
    CALL this%FinalizeSources()
  END SUBROUTINE Finalize

END MODULE sources_c_accel_mod
