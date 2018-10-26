!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2016                                                        #
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
!!
!! \brief constructor for physics class
!!
!! This module allocates the physics class and decides which specific
!! physics to use from the config.
!----------------------------------------------------------------------------!
MODULE physics_generic_mod
  USE physics_base_mod
  USE physics_eulerisotherm_mod
  USE physics_euler3Dit_mod
  USE physics_euler3D_mod
  USE physics_euler2Dit_mod
  USE physics_euler2D_mod
  USE mesh_base_mod
  USE common_dict

!  INTERFACE physics_base
!    MODULE PROCEDURE new_physics
!  END INTERFACE

CONTAINS

  SUBROUTINE new_physics(Physics,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), ALLOCATABLE :: Physics
    CLASS(mesh_base), INTENT(IN)     :: Mesh
    TYPE(DICT_TYP), POINTER          :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                          :: problem
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"problem",problem)

    ! allocate data
    SELECT CASE(problem)
    CASE(EULER_ISOTHERM)
      ALLOCATE(physics_eulerisotherm::Physics)
    CASE(EULER3D_ISOTHERM)
      ALLOCATE(physics_euler3Dit::Physics)
    CASE(EULER3D)
      ALLOCATE(physics_euler3D::Physics)
    CASE(EULER2D_ISOTHERM)
      ALLOCATE(physics_euler2Dit::Physics)
    CASE(EULER2D)
      ALLOCATE(physics_euler2D::Physics)
    CASE DEFAULT
      CALL Physics%Error("new_physics","Unknown physics type.")
    END SELECT

    ! call initialization
    SELECT TYPE(obj => Physics)
    TYPE IS (physics_eulerisotherm)
      CALL obj%InitPhysics_eulerisotherm(Mesh,config,IO)
    TYPE IS (physics_euler3Dit)
      CALL obj%InitPhysics_euler3Dit(Mesh,config,IO)
    TYPE IS (physics_euler3D)
      CALL obj%InitPhysics_euler3D(Mesh,config,IO)
    TYPE IS (physics_euler2Dit)
      CALL obj%InitPhysics_euler2Dit(Mesh,config,IO)
    TYPE IS (physics_euler2D)
      CALL obj%InitPhysics_euler2D(Mesh,config,IO)
    END SELECT
  END SUBROUTINE
END MODULE physics_generic_mod
