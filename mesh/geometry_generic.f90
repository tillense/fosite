!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2016 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
!#                    Jannes Klee      <jklee@astrophysik.uni-kiel.de>       #
!#                    Manuel Jung      <mjung@astrophysik.uni-kiel.de>       #
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
!> \author Jannes Klee
!! \author Manuel Jung
!!
!! \brief constructor for geometry class
!!
!! This module allocates the geometry class and decides which specific
!! geometry to use from the config.
!----------------------------------------------------------------------------!
MODULE geometry_generic_mod
  USE geometry_base_mod
  USE geometry_cartesian_mod
  USE geometry_cylindrical_mod
  USE geometry_logcylindrical_mod
  USE geometry_spherical_mod
  USE common_dict

!  INTERFACE geometry_base
!    MODULE PROCEDURE new_geometry
!  END INTERFACE

CONTAINS

  SUBROUTINE new_geometry(Geometry,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), ALLOCATABLE :: Geometry
    TYPE(DICT_TYP),       POINTER     :: config
    !------------------------------------------------------------------------!
    INTEGER :: geometry_type
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"geometry",geometry_type)
    ! allocate data
    SELECT CASE(geometry_type)
    CASE(CARTESIAN)
      ALLOCATE(geometry_cartesian::Geometry)
    CASE(CYLINDRICAL)
      ALLOCATE(geometry_cylindrical::Geometry)
    CASE(LOGCYLINDRICAL)
      ALLOCATE(geometry_logcylindrical::Geometry)
    CASE(SPHERICAL)
      ALLOCATE(geometry_spherical::Geometry)
    CASE DEFAULT
      CALL Geometry%Error("new_geometry","Unknown geometry")
    END SELECT

    ! call initialization
    SELECT TYPE(geometry_child => Geometry)
    TYPE IS (geometry_cartesian)
      CALL geometry_child%InitGeometry_cartesian(config)
    TYPE IS (geometry_cylindrical)
      CALL geometry_child%InitGeometry_cylindrical(config)
    TYPE IS (geometry_logcylindrical)
      CALL geometry_child%InitGeometry_logcylindrical(config)
    TYPE IS (geometry_spherical)
      CALL geometry_child%InitGeometry_spherical(config)
    END SELECT
  END SUBROUTINE
END MODULE geometry_generic_mod
