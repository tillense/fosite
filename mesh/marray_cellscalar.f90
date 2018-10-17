!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_cellscalar.f90                                             #
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
!! \brief derived mesh array class for scalar cell data
!!
!! \extends marray_base
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_cellscalar_mod
  USE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> data types and methods
  TYPE, EXTENDS(marray_base) :: marray_cellscalar
    REAL, DIMENSION(:,:,:), POINTER   :: center, &       !< geometric center
                                         bcenter         !< bary center

    REAL, DIMENSION(:,:,:,:), POINTER :: faces, &        !< cell face centers
                                         corners         !< cell corners
    CONTAINS
    PROCEDURE :: AssignPointers
  END TYPE
  INTERFACE marray_cellscalar
    MODULE PROCEDURE CreateMArray_cellscalar
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_cellscalar

CONTAINS

  FUNCTION CreateMArray_cellscalar() RESULT(new_cs)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_cellscalar) :: new_cs
    !-------------------------------------------------------------------!
    ! create new rank 1 mesh array
    new_cs = marray_base(1+1+6+8)
    ! assign array pointers
    new_cs%center  => new_cs%RemapBounds(new_cs%data4d(:,:,:,1))
    new_cs%bcenter => new_cs%RemapBounds(new_cs%data4d(:,:,:,2))
    new_cs%faces   => new_cs%RemapBounds(new_cs%data4d(:,:,:,3:8))
    new_cs%corners => new_cs%RemapBounds(new_cs%data4d(:,:,:,9:16))
  END FUNCTION CreateMArray_cellscalar
  
  SUBROUTINE AssignPointers(this)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_cellscalar),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%marray_base%AssignPointers()
    ! assign array pointers
    this%center  => this%RemapBounds(this%data4d(:,:,:,1))
    this%bcenter => this%RemapBounds(this%data4d(:,:,:,2))
    this%faces   => this%RemapBounds(this%data4d(:,:,:,3:8))
    this%corners => this%RemapBounds(this%data4d(:,:,:,9:16))    
  END SUBROUTINE AssignPointers
  
END MODULE marray_cellscalar_mod
