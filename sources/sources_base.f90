!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2019                                                   #
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
     REAL                            :: gparam         !< geometry parameter
     LOGICAL                         :: update_disk_height !< enable/disable computation of disk scale height
     TYPE(marray_base)               :: invheight2   !< energy sink due to cooling
     TYPE(marray_base)               :: height       !< energy sink due to cooling
     TYPE(marray_base)               :: h_ext        !< energy sink due to cooling
     TYPE(marray_base)               :: pot          !< gravitational potential
     REAL, DIMENSION(:,:,:), POINTER :: invr         !< 1./radius
     INTEGER                         :: use_envelope !< enable vicosity envelope

  CONTAINS

    PROCEDURE :: InitSources
    PROCEDURE (InfoSources),            DEFERRED :: InfoSources
    PROCEDURE :: ExternalSources
    PROCEDURE (ExternalSources_single), DEFERRED :: ExternalSources_single
    PROCEDURE :: CalcTimestep
    PROCEDURE (CalcTimestep_single),    DEFERRED :: CalcTimestep_single
    PROCEDURE :: GetSourcesPointer
    PROCEDURE (Finalize),               DEFERRED :: Finalize
    PROCEDURE :: Finalize_base
  END TYPE sources_base
  ABSTRACT INTERFACE
    SUBROUTINE InfoSources(this,Mesh)
      IMPORT sources_base, Mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(Sources_base),INTENT(IN) :: this
      CLASS(Mesh_base),INTENT(IN)    :: Mesh
    END SUBROUTINE
    SUBROUTINE CalcTimestep_single(this,Mesh,Physics,Fluxes,pvar,cvar,time,dt)
      IMPORT Sources_base, Mesh_base, Physics_base, Fluxes_base, marray_compound
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(sources_base),INTENT(INOUT) :: this
      CLASS(mesh_base),INTENT(IN)         :: Mesh
      CLASS(physics_base),INTENT(INOUT)   :: Physics
      CLASS(fluxes_base),INTENT(IN)       :: Fluxes
      CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
      REAL,INTENT(IN)                       :: time
      REAL, INTENT(OUT)                     :: dt
    END SUBROUTINE
    SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
      IMPORT sources_base, mesh_base, physics_base, fluxes_base, marray_compound
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(sources_base),INTENT(INOUT)  :: this
      CLASS(mesh_base),INTENT(IN)        :: Mesh
      CLASS(physics_base),INTENT(INOUT)  :: Physics
      CLASS(sources_base), INTENT(INOUT) :: Sources
      CLASS(fluxes_base),INTENT(IN)      :: Fluxes
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
  ! tempory storage for source terms
  CLASS(marray_compound), POINTER, SAVE :: temp_sterm => null()
  ! flags for source terms
  INTEGER, PARAMETER :: GRAVITY          = 1
!  INTEGER, PARAMETER :: DISK_THOMSON     = 2
  INTEGER, PARAMETER :: VISCOSITY        = 3
  INTEGER, PARAMETER :: C_ACCEL          = 4
!  INTEGER, PARAMETER :: COOLING          = 5
  INTEGER, PARAMETER :: ROTATING_FRAME   = 20
!  INTEGER, PARAMETER :: SGS              = 23
  INTEGER, PARAMETER :: DISK_COOLING     = 24
!  INTEGER, PARAMETER :: WAVE_DAMPING     = 25
!  INTEGER, PARAMETER :: FORCING          = 26
!  INTEGER, PARAMETER :: PLANET_HEATING   = 27
!  INTEGER, PARAMETER :: PLANET_COOLING   = 28
!  INTEGER, PARAMETER :: STELLAR_HEATING  = 29
  INTEGER, PARAMETER :: SHEARBOX         = 30
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       sources_base, &
       ! constants
       VISCOSITY, C_ACCEL, SHEARBOX, GRAVITY, DISK_COOLING, ROTATING_FRAME
  !--------------------------------------------------------------------------!

CONTAINS

  !> Initialize data in sources
  SUBROUTINE InitSources(this,Mesh,Fluxes,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), INTENT(IN)    :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    CLASS(physics_base), INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER            :: config, IO
    !------------------------------------------------------------------------!
    CALL Physics%new_statevector(temp_sterm,CONSERVATIVE)

    CALL this%Info(" SOURCES--> source term:       " // this%GetName())
    CALL this%InfoSources(Mesh)
  END SUBROUTINE InitSources


  SUBROUTINE ExternalSources(this,Mesh,Fluxes,Physics,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), TARGET, INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)         :: Mesh
    CLASS(fluxes_base),  INTENT(IN)         :: Fluxes
    CLASS(physics_base), INTENT(INOUT)      :: Physics
    REAL,                INTENT(IN)         :: time,dt
    CLASS(marray_compound), INTENT(INOUT)   :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: srcptr
    !------------------------------------------------------------------------!
    ! reset sterm
    sterm%data1d(:) = 0.0
    ! go through all source terms in the list
    srcptr => this
    DO WHILE (ASSOCIATED(srcptr))

      CALL srcptr%ExternalSources_single(Mesh,Physics,Fluxes,this,time,dt,pvar,cvar,temp_sterm)

      ! add to the sources
      sterm%data1d(:) = sterm%data1d(:) + temp_sterm%data1d(:)

      ! next source term
      srcptr => srcptr%next

    END DO
    ! reset ghost cell data
    IF (Mesh%GINUM.GT.0) THEN
      sterm%data4d(Mesh%IGMIN:Mesh%IMIN-Mesh%IP1,:,:,:) = 0.0
      sterm%data4d(Mesh%IMAX+Mesh%IP1:Mesh%IGMAX,:,:,:) = 0.0
    END IF
    IF (Mesh%GJNUM.GT.0) THEN
      sterm%data4d(:,Mesh%JGMIN:Mesh%JMIN-Mesh%JP1,:,:) = 0.0
      sterm%data4d(:,Mesh%JMAX+Mesh%JP1:Mesh%JGMAX,:,:) = 0.0
    END IF
    IF (Mesh%GKNUM.GT.0) THEN
      sterm%data4d(:,:,Mesh%KGMIN:Mesh%KMIN-Mesh%KP1,:) = 0.0
      sterm%data4d(:,:,Mesh%KMAX+Mesh%KP1:Mesh%KGMAX,:) = 0.0
    END IF
  END SUBROUTINE ExternalSources


  SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes,pvar,cvar,time,dt,dtcause)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), TARGET, INTENT(IN)    :: this
    CLASS(mesh_base),            INTENT(IN)    :: Mesh
    CLASS(physics_base),         INTENT(INOUT) :: Physics
    CLASS(fluxes_base),          INTENT(IN)    :: Fluxes
    CLASS(marray_compound),      INTENT(INOUT) :: pvar,cvar
    REAL, INTENT(IN)              :: time
    REAL, INTENT(OUT)             :: dt
    INTEGER, INTENT(OUT)          :: dtcause
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: srcptr
    REAL              :: dt_new
    !------------------------------------------------------------------------!
    dt_new = dt

    ! go through all source terms in the list
    srcptr => this
    DO WHILE(ASSOCIATED(srcptr))

       CALL srcptr%CalcTimestep_single(Mesh,Physics,Fluxes,pvar,cvar,time,dt_new)


       IF (dt_new .LT. dt) dtcause=srcptr%GetType()
       dt = MIN(dt,dt_new)
       ! next source term
       srcptr => srcptr%next
    END DO
  END SUBROUTINE CalcTimestep

  !> \public
  FUNCTION GetSourcesPointer(list,stype) RESULT(sp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), TARGET, INTENT(IN) :: list
    CLASS(sources_base), POINTER :: sp
    INTEGER, INTENT(IN)          :: stype
    !------------------------------------------------------------------------!
    sp => list
    DO
       IF (ASSOCIATED(sp).EQV..FALSE.) EXIT
!CDIR IEXPAND
       IF (sp%GetType().EQ.stype) RETURN
       sp => sp%next
    END DO
  END FUNCTION GetSourcesPointer

  !> Destructor
  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("sources_base::Finalize_base","not initialized")

    IF(ASSOCIATED(temp_sterm)) THEN
      IF (ASSOCIATED(temp_sterm%data1d)) CALL temp_sterm%Destroy()
      DEALLOCATE(temp_sterm)
      NULLIFY(temp_sterm)
    END IF
  END SUBROUTINE Finalize_base

END MODULE sources_base_mod
