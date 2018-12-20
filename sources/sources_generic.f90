!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2016                                                        #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!!
!! \brief constructor for mesh class
!!
!! This module allocates the mesh class and decides which specific
!! mesh to use from the config.
!----------------------------------------------------------------------------!
MODULE sources_generic_mod
  USE sources_base_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE common_dict

!  INTERFACE sources_base
!    MODULE PROCEDURE new_sources
!  END INTERFACE

CONTAINS

  SUBROUTINE new_sources(Sources,Mesh,Fluxes,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), POINTER :: Sources
    CLASS(mesh_base), INTENT(IN) :: Mesh
    CLASS(fluxes_base), INTENT(IN) :: Fluxes
    CLASS(physics_base), INTENT(IN) :: Physics
!    CLASS(timedisc_base), INTENT(IN) :: Timedisc
    TYPE(DICT_TYP), POINTER       :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: stype
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"stype",stype)

!    ! allocate data
!    SELECT CASE(meshtype)
!    CASE(MIDPOINT)
!      ALLOCATE(mesh_midpoint::new_mesh)
!    CASE(TRAPEZOIDAL)
!      ALLOCATE(mesh_trapezoidal::new_mesh)
!    END SELECT
!
!    ! call initialization
!    SELECT TYPE(mesh_child => new_mesh)
!    TYPE IS (mesh_midpoint)
!      CALL mesh_child%InitMesh_midpoint(config,IO)
!    TYPE IS (mesh_trapezoidal)
!      CALL mesh_child%InitMesh_trapezoidal(config,IO)
!    END SELECT
  END SUBROUTINE
END MODULE sources_generic_mod
