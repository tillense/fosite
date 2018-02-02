!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fluxes_generic.f90                                                #
!#                                                                           #
!# Copyright (C) 2016 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
!#                    Jannes Klee      <jklee@astrophysik.uni-kiel.de>       #
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
!! \author Jannes Klee
!!
!! \brief constructor for fluxes class
!!
!! This module allocates the fluxes class and decides which specific
!! fluxtype to use from the config.
!----------------------------------------------------------------------------!
MODULE fluxes_generic_mod
  USE fluxes_base_mod
  USE fluxes_kt_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: new_fluxes

CONTAINS

  SUBROUTINE new_fluxes(Fluxes,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_base), ALLOCATABLE :: Fluxes
    CLASS(mesh_base), INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    TYPE(DICT_TYP),     POINTER     :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                         :: fluxtype
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"fluxtype",fluxtype)
    ! allocate data
    SELECT CASE(fluxtype)
    CASE(KT)
      ALLOCATE(fluxes_kt::Fluxes)
    END SELECT

    ! call initialization
    SELECT TYPE(fluxes_child => Fluxes)
    TYPE IS (fluxes_kt)
      CALL fluxes_child%InitFluxes_KT(Mesh,Physics,config,IO)
    END SELECT
  END SUBROUTINE
END MODULE fluxes_generic_mod
