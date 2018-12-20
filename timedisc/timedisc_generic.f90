!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
!! \brief constructor for timedisc class
!!
!! This module allocates the mesh class and decides which specific
!! mesh to use from the config.
!----------------------------------------------------------------------------!
MODULE timedisc_generic_mod
  USE timedisc_base_mod
  USE timedisc_modeuler_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict

!  INTERFACE timedisc_base
!    MODULE PROCEDURE new_timedisc
!  END INTERFACE

CONTAINS

  SUBROUTINE new_timedisc(Timedisc,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), ALLOCATABLE :: Timedisc
    CLASS(mesh_base), INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    TYPE(DICT_TYP), POINTER       :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: method
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"method",method)

    ! allocate data
    SELECT CASE(method)
    CASE(MODIFIED_EULER)
      ALLOCATE(timedisc_modeuler::Timedisc)
    CASE DEFAULT
      ALLOCATE(timedisc_modeuler::Timedisc)
      CALL Timedisc%Error("new_timedisc","Unknown timedisc integration scheme")
    END SELECT
    
    ! call initialization
    SELECT TYPE(obj => Timedisc)
    TYPE IS (timedisc_modeuler)
      CALL obj%InitTimedisc_modeuler(Mesh,Physics,config,IO)
    END SELECT
  END SUBROUTINE
END MODULE timedisc_generic_mod
