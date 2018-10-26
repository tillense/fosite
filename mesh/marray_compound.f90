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
    PROCEDURE :: AssignMArray_0
    PROCEDURE :: AppendMArray
    PROCEDURE :: Destroy
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
    TYPE(marray_compound) :: new_cp
    !-------------------------------------------------------------------!
    ! create new empty mesh array of rank 2 with dims=0
    new_cp = marray_base(0,0)
    ALLOCATE(new_cp%item(0))
    ! be sure rank is 2, even if the compound consists only of scalar data
    new_cp%RANK = 2
  END FUNCTION CreateMArray_compound
  
  !> assigns one compound of mesh arrays to another compound of mesh arrays
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_compound),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    CALL this%marray_base%AssignMArray_0(ma)
    SELECT TYPE(src => ma)
    CLASS IS(marray_compound)
      IF (.NOT.ALLOCATED(src%item)) THEN
        ! error
      ELSE
        ALLOCATE(this%item(SIZE(src%item)),SOURCE=src%item)
      END IF
    CLASS DEFAULT
      ! error
    END SELECT
  END SUBROUTINE AssignMArray_0
  
  
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

      IF (ASSOCIATED(this%data1d).AND.SIZE(this%data1d).GT.0) THEN
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
        IF (ALLOCATED(this%item)) DEALLOCATE(this%item) ! maybe allocated before with wrong size
        ALLOCATE(this%item(1))
        this%item(1)%RANK = ma%RANK
        this%item(1)%DIMS(:) = ma%DIMS(:)
      END IF
      ! print debug info
!       PRINT '(A,I2,A)',"item no. ",SIZE(this%item)," appended to compound"
!       PRINT '(2(A,I6),A,2(I3))', "  size", SIZE(this%data1d), &
!                         "  rank", this%item(SIZE(this%item))%RANK, &
!                         "  dims", this%item(SIZE(this%item))%DIMS(:)
      ! assign new multi-dim pointers in the compound
      CALL this%AssignPointers()
    END IF
  END SUBROUTINE AppendMArray
  
  !> deconstructor of the compound
  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    !-------------------------------------------------------------------!
    ! free meta data array
    DEALLOCATE(this%item)
    ! call deconstructor of the base class
    CALL this%marray_base%Destroy()
  END SUBROUTINE Destroy

END MODULE marray_compound_mod
