!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2016-2024                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Manuel Jung
!! \author Lars Boesch
!! \author Jannes Klee
!! \author Tobias Illenseer
!!
!! \brief module to manage list of source terms
!!
!! This module allocates the sources class and decides which specific
!! source to use from the config.
!----------------------------------------------------------------------------!
MODULE sources_generic_mod
  USE sources_base_mod
  USE sources_c_accel_mod
  USE sources_cooling_mod
  USE sources_diskcooling_mod
  USE sources_gravity_mod
  USE sources_planetheating_mod
  USE sources_planetcooling_mod
  USE sources_rotframe_mod
  USE sources_shearbox_mod
  USE sources_viscosity_mod
  USE marray_compound_mod
  USE marray_base_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> container class to manage the list of source terms
  TYPE, EXTENDS(sources_base) :: sources_list
    !> \name Variables
    CLASS(marray_compound), POINTER :: sterm => null() !< tempory storage for source terms

  CONTAINS
    !> \name Methods
    PROCEDURE :: InitSources
    PROCEDURE :: ExternalSources
    PROCEDURE :: CalcTimestep
    FINAL     :: Finalize
  END TYPE sources_list
  !--------------------------------------------------------------------------!
  ! flags for source terms
  INTEGER, PARAMETER :: LIST             = 0
  INTEGER, PARAMETER :: GRAVITY          = 1
!  INTEGER, PARAMETER :: DISK_THOMSON     = 2
  INTEGER, PARAMETER :: VISCOSITY        = 3
  INTEGER, PARAMETER :: C_ACCEL          = 4
  INTEGER, PARAMETER :: COOLING          = 5
  INTEGER, PARAMETER :: ROTATING_FRAME   = 20
!  INTEGER, PARAMETER :: SGS              = 23
  INTEGER, PARAMETER :: DISK_COOLING     = 24
!  INTEGER, PARAMETER :: WAVE_DAMPING     = 25
!  INTEGER, PARAMETER :: FORCING          = 26
  INTEGER, PARAMETER :: PLANET_HEATING   = 27
  INTEGER, PARAMETER :: PLANET_COOLING   = 28
!  INTEGER, PARAMETER :: STELLAR_HEATING  = 29
  INTEGER, PARAMETER :: SHEARBOX         = 30
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    ! types
    sources_list, &

    ! constants
    VISCOSITY, C_ACCEL, SHEARBOX, GRAVITY, DISK_COOLING, ROTATING_FRAME, &
    PLANET_HEATING, PLANET_COOLING, COOLING, &
    MOLECULAR,ALPHA,BETA,POWERLAW,ALPHA_ALT,viscosity_name, &
    GRAY, GAMMIE, GAMMIE_SB


CONTAINS

  SUBROUTINE InitSources(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_list), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    TYPE(Dict_TYP), POINTER            :: config, IO
    !------------------------------------------------------------------------!
    CLASS(sources_base), POINTER :: newsrc => null(), tmpsrc => null()
    TYPE(Dict_TYP),      POINTER :: dir,next,src,IOsrc,gsrc => null(),gdir => null()
    INTEGER                      :: update_disk_height = 0
    INTEGER                      :: stype
    !------------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("sources_generic::InitSources","physics and/or mesh module uninitialized")
    IF (this%Initialized()) &
      CALL this%Error("sources_generic::InitSources","list of source term already initialized")

    CALL this%InitLogging(LIST,"sources list")

    ALLOCATE(this%sterm,STAT=this%err)
    IF (this%err.NE.0) &
      CALL this%Error("sources_generic::InitSources","memory allocation failed")

    CALL Physics%new_statevector(this%sterm,CONSERVATIVE)

    dir => config
    DO WHILE(ASSOCIATED(dir))
      NULLIFY(IOsrc)
      IF(HasChild(dir)) THEN
        src => GetChild(dir)
        CALL GetAttr(src, "stype", stype)

        ! object creation
        SELECT CASE(stype)
        CASE(GRAVITY)
           ! skip initialization of gravity modules here and initialize them
           ! at the end to make sure gravity is the first source term in the list
           gsrc => src
           gdir => dir
           NULLIFY(newsrc)
        CASE(C_ACCEL)
          ALLOCATE(sources_c_accel::newsrc)
        CASE(COOLING)
          ALLOCATE(sources_cooling::newsrc)
        CASE(DISK_COOLING)
          ALLOCATE(sources_diskcooling::newsrc)
        CASE(PLANET_HEATING)
          ALLOCATE(sources_planetheating::newsrc)
        CASE(PLANET_COOLING)
          ALLOCATE(sources_planetcooling::newsrc)
        CASE(ROTATING_FRAME)
          tmpsrc => this%GetSourcesPointer(ROTATING_FRAME)
          IF (ASSOCIATED(tmpsrc)) &
            CALL this%Error("sources_generic::InitSources","only one rotating frame source term allowed")
          ALLOCATE(sources_rotframe::newsrc)
        CASE(SHEARBOX)
          tmpsrc => this%GetSourcesPointer(SHEARBOX)
          IF (ASSOCIATED(tmpsrc)) &
            CALL this%Error("sources_generic::InitSources","only one shearing box source term allowed")
          ALLOCATE(sources_shearbox::newsrc)
        CASE(VISCOSITY)
          ALLOCATE(sources_viscosity::newsrc)
        CASE DEFAULT
          CALL this%Error("sources_generic::InitSources","Unknown source type")
        END SELECT

        ! basic initialization of all source terms except gravity
        IF (ASSOCIATED(newsrc)) THEN
          CALL newsrc%InitSources(Mesh,Physics,Fluxes,src,IOsrc)
          SELECT TYPE(obj => newsrc)
          TYPE IS (sources_diskcooling)
            IF (obj%cooling%GetType().EQ.GRAY) update_disk_height = 1
          TYPE IS (sources_viscosity)
            IF (obj%viscosity%GetType().EQ.ALPHA_ALT) update_disk_height = 1
          END SELECT

          ! prepend new source term to list
          newsrc%next => this%next
          this%next => newsrc
        END IF

        ! process the output dictionary
        IF(ASSOCIATED(IOsrc)) CALL SetAttr(IO, GetKey(dir), IOsrc)

      END IF
      dir => GetNext(dir)
    END DO

    ! finally initialize gravity
    IF(ASSOCIATED(gsrc)) THEN
      NULLIFY(IOsrc)

      ALLOCATE(sources_gravity::newsrc)
      IF (ASSOCIATED(newsrc)) THEN
        CALL SetAttr(gsrc,"update_disk_height", update_disk_height)
        CALL newsrc%InitSources(Mesh,Physics,Fluxes,gsrc,IOsrc)
        ! prepend new source term to list
        newsrc%next => this%next
        this%next => newsrc
        ! check for disk cooling source term and set pointer to gravity source
        tmpsrc => this%GetSourcesPointer(DISK_COOLING)
        IF (ASSOCIATED(tmpsrc)) THEN
          SELECT TYPE(sp => tmpsrc)
          CLASS IS(sources_diskcooling)
            SELECT TYPE(grav => newsrc)
            CLASS IS(sources_gravity)
              sp%grav => grav
            END SELECT
          END SELECT
        END IF
      END IF

      IF(ASSOCIATED(IOsrc)) CALL SetAttr(IO, GetKey(gdir), IOsrc)
    END IF

  END SUBROUTINE InitSources


  SUBROUTINE ExternalSources(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_list), INTENT(INOUT)      :: this
    CLASS(mesh_base),    INTENT(IN)         :: Mesh
    CLASS(physics_base), INTENT(INOUT)      :: Physics
    CLASS(fluxes_base),  INTENT(IN)         :: Fluxes
    CLASS(sources_base), INTENT(INOUT)      :: Sources
    REAL,                INTENT(IN)         :: time,dt
    CLASS(marray_compound), INTENT(INOUT)   :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    CLASS(sources_base), POINTER :: srcptr
    !------------------------------------------------------------------------!
    ! reset sterm
    sterm%data1d(:) = 0.0
    ! go through all source terms in the list
    srcptr => this%next
    DO WHILE (ASSOCIATED(srcptr))

      CALL srcptr%ExternalSources(Mesh,Physics,Fluxes,this,time,dt,pvar,cvar,this%sterm)

      ! add to the sources
      sterm%data1d(:) = sterm%data1d(:) + this%sterm%data1d(:)

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
    CLASS(sources_list),         INTENT(INOUT) :: this
    CLASS(mesh_base),            INTENT(IN)    :: Mesh
    CLASS(physics_base),         INTENT(INOUT) :: Physics
    CLASS(fluxes_base),          INTENT(IN)    :: Fluxes
    CLASS(marray_compound),      INTENT(INOUT) :: pvar,cvar
    REAL, INTENT(IN)              :: time
    REAL, INTENT(INOUT)           :: dt
    INTEGER, INTENT(OUT)          :: dtcause
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: srcptr
    REAL              :: dt_new
    !------------------------------------------------------------------------!
    dt_new = dt

    ! go through all source terms in the list
    srcptr => this%next
    DO WHILE(ASSOCIATED(srcptr))

       CALL srcptr%CalcTimestep(Mesh,Physics,Fluxes,pvar,cvar,time,dt_new,dtcause)


       IF (dt_new .LT. dt) dtcause=srcptr%GetType()
       dt = MIN(dt,dt_new)
       ! next source term
       srcptr => srcptr%next
    END DO
  END SUBROUTINE CalcTimestep

  !> \private destructor for the sources list
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_list) :: this
    !------------------------------------------------------------------------!
    CLASS(sources_base), POINTER :: srcptr,next
    !------------------------------------------------------------------------!
    ! loop over all source terms and call finalization subroutines
    srcptr => this%next
    DO WHILE (ASSOCIATED(srcptr))
      next => srcptr%next
      DEALLOCATE(srcptr)
      srcptr => next
    END DO
    IF (ASSOCIATED(this%sterm)) DEALLOCATE(this%sterm)
  END SUBROUTINE Finalize

END MODULE sources_generic_mod
