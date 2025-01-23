!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2016 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!! \brief constructor for constants class
!!
!! This module allocates the mesh class and decides which specific
!! mesh to use from the config.
!----------------------------------------------------------------------------!
MODULE constants_generic_mod
  USE constants_base_mod
  USE constants_cgs_mod
  USE constants_geometrical_mod
  USE constants_SI_mod

CONTAINS

  SUBROUTINE new_constants(Constants,units)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(constants_base), ALLOCATABLE :: Constants
    !------------------------------------------------------------------------!
    INTEGER                            :: units
    !------------------------------------------------------------------------!

    ! allocate data
    SELECT CASE(units)
    CASE(CGS)
      ALLOCATE(constants_cgs::Constants)
    CASE(GEOMETRICAL)
      ALLOCATE(constants_geometrical::Constants)
    CASE(SI)
      ALLOCATE(constants_SI::Constants)
    END SELECT

    ! call initialization
    SELECT TYPE(obj => Constants)
    TYPE IS (constants_cgs)
      CALL obj%InitConstants_cgs()
    TYPE IS (constants_geometrical)
      CALL obj%InitConstants_geometrical()
    TYPE IS (constants_SI)
      CALL obj%InitConstants_SI()
    END SELECT
  END SUBROUTINE
END MODULE constants_generic_mod
