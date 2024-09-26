!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2021                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!!
!! \brief Includes all source terms
!!
!! - general parameters of sources group as key-values
!! \key{stype,INTEGER,Type of source}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief generic source terms module providing functionaly common
!! to all source terms
!!
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE marray_compound_mod
  USE marray_base_mod
  USE gravity_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE common_dict
#if defined(HAVE_FFTW)
  USE fftw
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, ABSTRACT, EXTENDS(logging_base) :: sources_base
     !> \name Variables
     CLASS(sources_base), POINTER    :: next => null() !< next source in list
     REAL                            :: time           !< simulation time
     REAL                            :: cvis
!      REAL                            :: gparam         !< geometry parameter
!      REAL, DIMENSION(:,:,:), POINTER :: invr         !< 1./radius
!      INTEGER                         :: use_envelope !< enable vicosity envelope

  CONTAINS

    PROCEDURE :: InitSources_base
    PROCEDURE (InitSources),     DEFERRED :: InitSources
    PROCEDURE (ExternalSources), DEFERRED :: ExternalSources
    PROCEDURE (CalcTimestep),    DEFERRED :: CalcTimestep
    PROCEDURE :: GetSourcesPointer
  END TYPE sources_base
  ABSTRACT INTERFACE
    SUBROUTINE InitSources(this,Mesh,Physics,Fluxes,config,IO)
      IMPORT Sources_base, Mesh_base, Physics_base, Fluxes_base, marray_compound, Dict_TYP
      IMPLICIT NONE
      CLASS(sources_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
      TYPE(Dict_TYP), POINTER            :: config, IO
    END SUBROUTINE
    SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes,pvar,cvar,time,dt,dtcause)
      IMPORT Sources_base, Mesh_base, Physics_base, Fluxes_base, marray_compound
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(sources_base),INTENT(INOUT)   :: this
      CLASS(mesh_base),INTENT(IN)         :: Mesh
      CLASS(physics_base),INTENT(INOUT)   :: Physics
      CLASS(fluxes_base),INTENT(IN)       :: Fluxes
      CLASS(marray_compound),INTENT(INOUT):: pvar,cvar
      REAL,INTENT(IN)                     :: time
      REAL, INTENT(INOUT)                 :: dt
      INTEGER, INTENT(OUT)                :: dtcause
    END SUBROUTINE
    SUBROUTINE ExternalSources(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
      IMPORT sources_base, mesh_base, physics_base, fluxes_base, marray_compound
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(sources_base),INTENT(INOUT)  :: this
      CLASS(mesh_base),INTENT(IN)        :: Mesh
      CLASS(physics_base),INTENT(INOUT)  :: Physics
      CLASS(fluxes_base),INTENT(IN)      :: Fluxes
      CLASS(sources_base), INTENT(INOUT) :: Sources
      REAL,INTENT(IN)                    :: time, dt
      CLASS(marray_compound),INTENT(INOUT):: pvar,cvar,sterm
    END SUBROUTINE
    SUBROUTINE Finalize(this)
      IMPORT sources_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(sources_base),INTENT(INOUT) :: this
    END SUBROUTINE
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       sources_base
  !--------------------------------------------------------------------------!

CONTAINS

  !> Initialize data in sources
  SUBROUTINE InitSources_base(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), INTENT(INOUT) :: this
    INTEGER,             INTENT(IN)    :: stype
    CHARACTER(LEN=*),    INTENT(IN)    :: sname
    !------------------------------------------------------------------------!
    IF (this%Initialized()) &
      CALL this%Error("sources_base::InitSources_base","source term already initialized")
    CALL this%InitLogging(stype,sname)
    CALL this%Info(" SOURCES--> source term:       " // this%GetName())
!     CALL this%InfoSources(Mesh)
  END SUBROUTINE InitSources_base


  !> \public return pointer to requested source term
  FUNCTION GetSourcesPointer(this,stype) RESULT(sp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), TARGET, INTENT(IN) :: this
    CLASS(sources_base), POINTER    :: sp
    INTEGER, INTENT(IN)             :: stype
    !------------------------------------------------------------------------!
    sp => this
    DO WHILE (ASSOCIATED(sp))
!CDIR IEXPAND
      IF (sp%GetType().EQ.stype) EXIT
      sp => sp%next
    END DO
  END FUNCTION GetSourcesPointer


!   !> Destructor
!   SUBROUTINE Finalize_base(this)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(sources_base), INTENT(INOUT) :: this
!     !------------------------------------------------------------------------!
!     IF (.NOT.this%Initialized()) &
!       CALL this%Error("sources_base::Finalize_base"," called for uninitialized sources")
!     IF(ASSOCIATED(temp_sterm)) DEALLOCATE(temp_sterm)
!   END SUBROUTINE Finalize_base

END MODULE sources_base_mod
