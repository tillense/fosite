!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_cellvector.f90                                             #
!#                                                                           #
!# Copyright (C) 2018                                                        #
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
!! \brief derived mesh array class for vector cell data
!!
!! \extends marray_base
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_cellvector_mod
  USE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> data types and methods
  TYPE, EXTENDS(marray_base) :: marray_cellvector
    REAL, DIMENSION(:,:,:,:), POINTER   :: center, &     !< geometric center
                                         bcenter         !< bary center

    REAL, DIMENSION(:,:,:,:,:), POINTER :: faces, &      !< cell face centers
                                         corners         !< cell corners
    CONTAINS
    PROCEDURE :: AssignPointers
  END TYPE
  INTERFACE marray_cellvector
    MODULE PROCEDURE CreateMArray_cellvector
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_cellvector

CONTAINS

  FUNCTION CreateMArray_cellvector() RESULT(new_cv)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_cellvector) :: new_cv
    !-------------------------------------------------------------------!
    ! create new rank 1 mesh array
    new_cv = marray_base(1+1+6+8,3)
  END FUNCTION CreateMArray_cellvector
  
  SUBROUTINE AssignPointers(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_cellvector),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%marray_base%AssignPointers()
    ! assign array pointers
    this%center  => this%RemapBounds(this%data5d(:,:,:,1,:))
    this%bcenter => this%RemapBounds(this%data5d(:,:,:,2,:))
    this%faces   => this%RemapBounds(this%data5d(:,:,:,3:8,:))
    this%corners => this%RemapBounds(this%data5d(:,:,:,9:16,:))    
  END SUBROUTINE AssignPointers
  
END MODULE marray_cellvector_mod
