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
!> \author Lars Boesch
!!
!! \brief constructor for sources class
!!
!! This module allocates the sources class and decides which specific
!! source to use from the config.
!----------------------------------------------------------------------------!
MODULE sources_generic_mod
  USE sources_base_mod
  USE sources_viscosity_mod
  USE mesh_base_mod
  USE Fluxes_base_mod
  USE Physics_base_mod
  USE common_dict

CONTAINS


  SUBROUTINE New_Sources(this,Mesh,Fluxes,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),POINTER :: this
    CLASS(Mesh_base)    :: Mesh
    CLASS(Fluxes_base)  :: Fluxes
    CLASS(Physics_base) :: Physics
    TYPE(DICT_TYP), POINTER   :: config, IO

    !------------------------------------------------------------------------!
    CLASS(Sources_base),POINTER :: newsrc, tmpsrc
    !CLASS(Sources_base),POINTER :: errsrc      ! we need this only for error reporting !
    INTEGER           :: err
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    ! allocate memory for new source term
    CALL GetAttr(config,"vis/stype",stype)

    SELECT CASE(stype)
    CASE(Viscosity)
      ALLOCATE(sources_viscosity::newsrc)
      IF (err.NE.0) CALL newsrc%Error("New_Sources", "Unable allocate memory!")
    CASE DEFAULT
      CALL this%Error("New_Sources","Unknown source type")
    END SELECT
    
    ! basic initialization
    SELECT TYPE(obj => newsrc)
    TYPE IS (sources_viscosity)
      CALL obj%InitSources_all(Mesh,Fluxes,Physics,config,IO)
    END SELECT
    ! add new source term to beginning of
    ! list of source terms
    IF (.NOT.ASSOCIATED(this)) THEN
       this => newsrc
       NULLIFY(this%next)
    ELSE
       tmpsrc => this
       this => newsrc
       this%next => tmpsrc
    END IF
  END SUBROUTINE New_Sources




END MODULE sources_generic_mod
