!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_compound.f90                                               #
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
!! \brief derived class for compound of mesh arrays
!!
!! \extends marray_base
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_compound_mod
  USE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE item_properties
    INTEGER :: RANK, DIMS(2)
  END TYPE
  !> data types and methods
  TYPE, EXTENDS(marray_base) :: marray_compound
    TYPE(item_properties), ALLOCATABLE :: item(:)
    CONTAINS
    PROCEDURE :: AppendMArray
  END TYPE
  INTERFACE marray_compound
    MODULE PROCEDURE CreateMArray_compound
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_compound

CONTAINS

  FUNCTION CreateMArray_compound() RESULT(new_cp)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_base) :: ma
    TYPE(marray_compound) :: new_cp
    !-------------------------------------------------------------------!
    ! create new empty mesh array of rank 2 with dims=0
    new_cp = marray_base(0,0)
    ! be sure rank is 2, even if the compound consists only of scalar data
    this%RANK = 2
  END FUNCTION CreateMArray_cellscalar
  
  SUBROUTINE AppendMArray(this,ma)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(marray_base), INTENT(INOUT) :: ma
    !-------------------------------------------------------------------!
    REAL, DIMENSION(:),POINTER :: data1d
    TYPE(item_properties), ALLOCATABLE :: item(:)
    !-------------------------------------------------------------------!
    IF (.NOT.ASSOCIATED(ma%data1d)) THEN
      !> \todo this should not happen; think of how to implement error
      !! handling in mesh arrays
    ELSE
      this%DIMS(1) = ma%DIMS(1)*ma%DIMS(2) ! 1st dimension (excluding mesh
          ! dimensions) is collapsed into one index; if ma is a scalar, then this is 1;
          ! if ma is rank 1, then this is the vector component index;
          ! if ma is rank 2, the two indices are collapsed into one index
      this%DIMS(2) = this%DIMS(2)+1 ! 2nd dimension (excluding mesh dimensions)
          ! counts the items, i.e., the index of the mesh array in the compound;
          ! thus this%data5d(:,:,:,:,n) selects the nth mesh array in the compound

      IF (ASSOCIATED(this%data1d)) THEN
        ! allocate new contiguous data block for old+new data and for meta data
        ALLOCATE(data1d(SIZE(this%data1d)+SIZE(ma%data1d)),item(SIZE(this%item)+1))
        ! copy the meta data and move the allocation
        item(1:SIZE(this%item)) = this%item(:)
        item(SIZE(this%item)+1)%RANK = ma%RANK
        item(SIZE(this%item)+1)%DIMS(:) = ma%DIMS(:)
        CALL MOVE_ALLOC(item,this%item)
        ! copy the data
        data1d(1:SIZE(this%data1d)) = this%data1d(:)
        data1d(SIZE(this%data1d)+1:SIZE(this%data1d)+SIZE(ma%data1d)) = ma%data1d(:)
        ! assign 1D pointers
        ma%data1d => data1d(SIZE(this%data1d)+1:)
        this%data1d => data1d
        ! reassign multi-dim pointers in old mesh array
        CALL ma%AssignPointers()
      ELSE ! i.e., .NOT.ASSOCIATED(this%data1d)
        ! append ma to an empty compound -> just assign the pointer
        this%data1d => ma%data1d
        ! allocate memory for the meta data
        ALLOCATE(this%item(1))
        this%item(1)%RANK = ma%RANK
        this%item(1)%DIMS(:) = ma%DIMS(:)
      END IF
      ! assign new multi-dim pointers in the compound
      CALL this%AssignPointers()
    END IF
  END SUBROUTINE AppendMArray
  
END MODULE marray_compound_mod
