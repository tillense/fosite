!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_axis.f90                                                 #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!!
!! \brief Boundary module for axis boundaries
!!
!! \extends boundary_reflecting
!! \ingroup boundary
!> \todo compiles, but strange behavior at axis. Be carefull
!----------------------------------------------------------------------------!
MODULE boundary_axis_mod
  USE mesh_base_mod
  USE boundary_base_mod
  USE boundary_reflecting_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_reflecting) :: boundary_axis
  CONTAINS
    PROCEDURE                :: InitBoundary_axis
  END TYPE
  CHARACTER(LEN=32), PARAMETER :: boundcond_name = "axis"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       boundary_axis
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for the axis boundary condition
  SUBROUTINE InitBoundary_axis(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_axis),    INTENT(INOUT) :: this
    CLASS(physics_base),     INTENT(IN)    :: Physics
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    TYPE(Dict_TYP), POINTER, INTENT(IN)    :: config
    INTEGER                                :: dir
    !------------------------------------------------------------------------!
    INTEGER                                :: err
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,AXIS,boundcond_name,dir,config)

    ALLOCATE(this%reflX(Physics%vnum), &
         this%reflY(Physics%vnum),     &
         this%reflZ(Physics%vnum),     &
         STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_axis", "Unable to allocate memory.")
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IMPORTANT:
    ! I don't know if this works for all geometries and all physics
    ! along all coordinate directions;
    ! check the underlying physics/geometry  modules
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL Physics%AxisMasks(Mesh,this%reflX,this%reflY,this%reflZ)
  END SUBROUTINE InitBoundary_axis

END MODULE boundary_axis_mod
