!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2016-2018                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!!
!! \brief constructor for sources class
!!
!! This module allocates the sources class and decides which specific
!! source to use from the config.
!----------------------------------------------------------------------------!
MODULE sources_generic_mod
  USE sources_base_mod
  USE sources_c_accel_mod
  USE sources_diskcooling_mod
  USE sources_gravity_mod
  USE sources_rotframe_mod
  USE sources_shearbox_mod
  USE sources_viscosity_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE common_dict

CONTAINS

  SUBROUTINE new_sources(this,Mesh,Fluxes,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), POINTER :: this
    CLASS(mesh_base)             :: Mesh
    CLASS(fluxes_base)           :: Fluxes
    CLASS(physics_base)          :: Physics
    TYPE(DICT_TYP),      POINTER :: config, IO
    !------------------------------------------------------------------------!
    CLASS(sources_base), POINTER :: newsrc => null(), tmpsrc => null()
    TYPE(Dict_TYP),      POINTER :: dir,src,IOsrc,gsrc => null(),gdir => null()
    INTEGER                      :: update_disk_height = 0
    INTEGER                      :: stype
    !------------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("InitGravity","physics and/or mesh module uninitialized")

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
           newsrc => null()
        CASE(C_ACCEL)
          ALLOCATE(sources_c_accel::newsrc)
        CASE(DISK_COOLING)
          ALLOCATE(sources_diskcooling::newsrc)
        CASE(ROTATING_FRAME)
          ALLOCATE(sources_rotframe::newsrc)
        CASE(SHEARBOX)
          ALLOCATE(sources_shearbox::newsrc)
        CASE(VISCOSITY)
          ALLOCATE(sources_viscosity::newsrc)
        CASE DEFAULT
          CALL this%Error("new_sources","Unknown source type")
        END SELECT

        ! basic initialization of all source terms except gravity
        IF (ASSOCIATED(newsrc)) THEN
          SELECT TYPE(obj => newsrc)
          TYPE IS (sources_c_accel)
            CALL obj%InitSources_c_accel(Mesh,Physics,Fluxes,src,IOsrc)
          TYPE IS (sources_diskcooling)
            CALL obj%InitSources_diskcooling(Mesh,Physics,Fluxes,src,IOsrc)
            IF (obj%cooling%GetType().EQ.GRAY) update_disk_height = 1
          TYPE IS (sources_rotframe)
            CALL obj%InitSources_rotframe(Mesh,Physics,Fluxes,src,IOsrc)
          TYPE IS (sources_shearbox)
            CALL obj%InitSources_shearbox(Mesh,Physics,Fluxes,src,IOsrc)
          TYPE IS (sources_viscosity)
            CALL obj%InitSources_viscosity(Mesh,Physics,Fluxes,src,IOsrc)
            IF (obj%viscosity%GetType().EQ.ALPHA_ALT) update_disk_height = 1
          END SELECT
          ! check if list of source terms is empty
          IF (.NOT.ASSOCIATED(this)) THEN
            ! first entry is the new source term
            this => newsrc
            NULLIFY(this%next)
          ELSE
            ! prepend new source term to the list of source terms
            tmpsrc => this
            this => newsrc
            this%next => tmpsrc
          END IF
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

      tmpsrc => this
      this => newsrc
      this%next => tmpsrc

      CALL SetAttr(gsrc,"update_disk_height", update_disk_height)
      SELECT TYPE(obj => newsrc)
      TYPE IS(sources_gravity)
        CALL obj%InitSources_gravity(Mesh,Physics,Fluxes,gsrc,IOsrc)
      END SELECT

      IF(ASSOCIATED(IOsrc)) CALL SetAttr(IO, GetKey(gdir), IOsrc)
    END IF

  END SUBROUTINE New_Sources

  SUBROUTINE close_sources()
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! add new source term to beginning of
    ! list of source terms
  END SUBROUTINE

END MODULE sources_generic_mod
