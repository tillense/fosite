!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2016-2018                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!! \author Jannes Klee
!!
!! \brief constructor for gravity class
!!
!! This module allocates the gravity class and decides which specific
!! source to use from the config.
!----------------------------------------------------------------------------!
MODULE gravity_generic_mod
  USE gravity_base_mod
  USE gravity_binary_mod
  USE gravity_pointmass_mod
  USE gravity_sboxspectral_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE common_dict

CONTAINS

  SUBROUTINE new_gravity(this,Mesh,Fluxes,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER :: this
    CLASS(mesh_base)             :: Mesh
    CLASS(fluxes_base)           :: Fluxes
    CLASS(physics_base)          :: Physics
    TYPE(DICT_TYP),      POINTER :: config, IO
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER :: newgrav, tmpgrav
    TYPE(Dict_TYP),      POINTER :: dir,grav,IOgrav
    INTEGER                      :: gtype
    !------------------------------------------------------------------------!
    dir => config
    DO WHILE(ASSOCIATED(dir))
      NULLIFY(IOgrav)
      IF (HasChild(dir) .AND. (TRIM(GetKey(dir)).NE."output")) THEN
        grav => GetChild(dir)
        CALL GetAttr(grav, "gtype", gtype)

        ! object creation
        SELECT CASE(gtype)
        CASE(POINTMASS)
          ALLOCATE(gravity_pointmass::newgrav)
        CASE(POINTMASS_BINARY)
          ALLOCATE(gravity_binary::newgrav)
        CASE(SBOXSPECTRAL)
          ALLOCATE(gravity_sboxspectral::newgrav)
        CASE DEFAULT
          CALL this%Error("new_gravity","Unknown gravity type")
        END SELECT

        IF (.NOT.ASSOCIATED(this)) THEN
          this => newgrav
          NULLIFY(this%next)
        ELSE
          tmpgrav => this
          this => newgrav
          this%next => tmpgrav
        END IF

        SELECT TYPE(obj => newgrav)
        TYPE IS (gravity_pointmass)
           ! gravitational acceleration due to point mass
           CALL obj%InitGravity_pointmass(Mesh,Physics,grav,IOgrav)
        TYPE IS (gravity_binary)
           ! gravitational acceleration due to two circling point masses
           CALL obj%InitGravity_binary(Mesh,Physics,grav,IOgrav)
        TYPE IS (gravity_sboxspectral)
           ! self-gravitation in flat geometries periodic in both dimensions
           CALL obj%InitGravity_sboxspectral(Mesh,Physics,grav,IOgrav)
        END SELECT

        ! print some information
!        IF (ASSOCIATED(this)) THEN
!           CALL this%nextg%Info(" GRAVITY--> Gravity term:      " // GetName(this%glist))
!           ! print setup information of the individual Gravity terms
!           CALL InfoGravity(this%nextg)
!        END IF

        IF(ASSOCIATED(IOgrav)) CALL SetAttr(IO, GetKey(dir), IOgrav)
      END IF
      dir => GetNext(dir)
    END DO

  END SUBROUTINE new_gravity

  SUBROUTINE CloseGravity()
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! add new source term to beginning of
    ! list of source terms
  END SUBROUTINE CloseGravity

END MODULE gravity_generic_mod
