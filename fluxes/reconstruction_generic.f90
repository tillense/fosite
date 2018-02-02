!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: reconstruction_generic.f90                                        #
!#                                                                           #
!# Copyright (C) 2016                                                        #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief constructor for reconstruction class
!!
!! This module allocates the reconstruction class and decides which specific
!! fluxtype to use from the config.
!----------------------------------------------------------------------------!
MODULE reconstruction_generic_mod
  USE reconstruction_base_mod
  USE reconstruction_constant_mod
  USE reconstruction_linear_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict

!  INTERFACE reconstruction_base
!    MODULE PROCEDURE new_reconstruction
!  END INTERFACE

CONTAINS

  SUBROUTINE new_reconstruction(Reconstruction,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_base), ALLOCATABLE :: Reconstruction
    CLASS(mesh_base), INTENT(IN)            :: Mesh
    CLASS(physics_base), INTENT(IN)         :: Physics
    TYPE(DICT_TYP),     POINTER             :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                                 :: order
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"order",order)
    ! allocate data
    SELECT CASE(order)
    CASE(CONSTANT)
      ALLOCATE(reconstruction_constant::Reconstruction)
    CASE(LINEAR)
      ALLOCATE(reconstruction_linear::Reconstruction)
    END SELECT

    ! call initialization
    SELECT TYPE(reconstruction_child => Reconstruction)
    TYPE IS (reconstruction_constant)
      CALL reconstruction_child%InitReconstruction_constant(Mesh,Physics,config,IO)
    TYPE IS (reconstruction_linear)
      CALL reconstruction_child%InitReconstruction_linear(Mesh,Physics,config,IO)
    END SELECT
  END SUBROUTINE
END MODULE reconstruction_generic_mod
