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
!! \todo improve error handling for compound arrays
!!
!! \extends marray_base
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_compound_mod
  USE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE compound_item
    TYPE(marray_base), POINTER :: item => null()
    TYPE(compound_item), POINTER :: next => null()
    INTEGER :: extent,entry_num
  END TYPE
  !> data types and methods
  TYPE, EXTENDS(marray_base) :: marray_compound
    PRIVATE
    TYPE(compound_item), POINTER :: list => null()
    INTEGER :: num_entries = 0
  CONTAINS
    PROCEDURE :: AssignPointers
    PROCEDURE :: AssignMArray_0
    PROCEDURE :: AppendMArray
    PROCEDURE :: LastItem
    PROCEDURE :: FirstItem
    PROCEDURE :: NextItem
    PROCEDURE :: GetItem
    PROCEDURE :: AppendItem
    PROCEDURE :: Destroy
  END TYPE
  INTERFACE marray_compound
    MODULE PROCEDURE CreateMArray_compound
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_compound

CONTAINS

  !> constructor for compound of mesh arrays
  !!
  !! There exist two types of compounds differing by their rank:
  !! - rank 0 compounds are associated with rank 1 mesh arrays;
  !!   hence they basically store one vector of arbitrary dimension
  !!   in each cell and are mostly used for cell center data
  !! - rank 1 compounds are associated with rank 2 mesh arrays;
  !!   they store an arbitrary number of vectors of arbitrary
  !!   dimension in each cell and are mostly used for cell face data
  FUNCTION CreateMArray_compound(n) RESULT(new_cp)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_compound) :: new_cp
    INTEGER, OPTIONAL, INTENT(IN) :: n
    !-------------------------------------------------------------------!
    IF (PRESENT(n)) THEN
      IF (n.LE.0) THEN
        PRINT *,"ERROR in CreateMArray_compound: n should be larger than 0"
        STOP 1
      ELSE
        ! create empty mesh array of rank 2 with first dimension set to n
        new_cp = marray_base(n,0)
        new_cp%RANK = 2
      END IF
    ELSE
      ! default is rank 1 mesh array
      new_cp = marray_base(0)
      new_cp%RANK = 1
    END IF
  END FUNCTION CreateMArray_compound
  
  !> \public assign pointers of different shapes to the 1D data
  SUBROUTINE AssignPointers(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_compound),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p
    INTEGER :: m,n
    !------------------------------------------------------------------------!
    ! set the multi-dim. pointers of the compound
    CALL this%marray_base%AssignPointers()
    ! go through the list of items and set the 1D
    ! and multi-dim. pointers for each item
    p => this%FirstItem()
    m = 0
    n = 0
    DO WHILE (ASSOCIATED(p))
      m = m + n
      n = p%extent
      ! point to the 1D data segment in the new array
      p%item%data1d(1:n) => this%data1d(m+1:m+n)
      ! restore all multi-dim. pointers
      CALL p%item%AssignPointers()
      p => this%NextItem(p)
    END DO
  END SUBROUTINE AssignPointers

  !> assigns one compound of mesh arrays to another compound of mesh arrays
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_compound),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p
    !------------------------------------------------------------------------!
    CALL this%marray_base%AssignMArray_0(ma)
    SELECT TYPE(src => ma)
    CLASS IS(marray_compound)
      this%num_entries = src%num_entries
      p => src%FirstItem()
      IF (.NOT.ASSOCIATED(p).OR.src%num_entries.LT.1) THEN
        this%num_entries = 0
        NULLIFY(this%list)
      ELSE
        ! p is associated => generate new this%list with data from p, i.e. src
        DO
          CALL this%AppendItem(p%item)
          IF (ASSOCIATED(p%next)) THEN
            p => p%next
          ELSE
            EXIT
          END IF
        END DO
        ! finally reassign all pointers
        CALL this%AssignPointers()
      END IF
! IF (ASSOCIATED(src%data1d).AND.SIZE(src%data1d).GT.0) &
! PRINT *,"marray_compound::AssignMArray_0 ma",src%data1d(LBOUND(src%data1d,1)),src%data1d(UBOUND(src%data1d,1))
! IF (ASSOCIATED(this%data1d).AND.SIZE(this%data1d).GT.0) &
! PRINT *,"marray_compound::AssignMArray_0 this",this%data1d(LBOUND(this%data1d,1)),this%data1d(UBOUND(this%data1d,1))
    CLASS DEFAULT
      ! do nothing, is this ok ???
    END SELECT
  END SUBROUTINE AssignMArray_0
  
  
  SUBROUTINE AppendMArray(this,ma)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound), INTENT(INOUT) :: this
    TYPE(marray_base), POINTER :: ma
    !-------------------------------------------------------------------!
    INTEGER :: m,n
    REAL, DIMENSION(:),POINTER :: data1d
    !-------------------------------------------------------------------!
    IF (.NOT.(ASSOCIATED(ma%data1d).AND.ASSOCIATED(this%data1d))) THEN
      PRINT *,"ERROR in marray_compound::AppendMArray: at least one of this%data1d,ma%data1d is not associated"
      STOP 1
    ELSE
      ! check the rank of the compound (either 1 or 2)
      SELECT CASE(this%RANK)
      CASE(1) ! rank 1 compound
        this%DIMS(1) = this%DIMS(1) + ma%DIMS(1)*ma%DIMS(2) ! 1st dimension (excluding mesh
            ! dimensions) is collapsed into one index; if ma is a scalar, then this is 1;
            ! if ma is rank 1, then this is the vector component index;
            ! if ma is rank 2, the two indices are collapsed into one index
      CASE(2) ! rank 2 compound
        ! check if the 1st dimension (excluding mesh dimensions) are the same for
        ! both compound and mesh array added to the compound
        IF (this%DIMS(1).NE.ma%DIMS(1)) THEN
          PRINT *,"ERROR in marray_compound::AppendMArray: this%dims(1) != ma%dims(1)"
          STOP 1
        ELSE
          ! this%DIMS(1) remains the same
          this%DIMS(2) = this%DIMS(2) + ma%DIMS(2) ! extend the last dimension of the compound;
            ! so far only scalars and vectors can be added to rank 2 compounds
        END IF
      CASE DEFAULT
        ! this should not happen
      END SELECT
      IF (ASSOCIATED(this%data1d)) THEN
        m = SIZE(this%data1d(:))
        IF (m.EQ.0) DEALLOCATE(this%data1d) ! associated but size 0
      ELSE
        m = 0
      END IF
      IF (m.GT.0) THEN
        n = SIZE(ma%data1d(:))
        ! allocate new contiguous data block for old+new data and for meta data
        ALLOCATE(data1d(m+n))
        ! copy the data, new data is appended to the old compound
        data1d(1:m) = this%data1d(1:m)
        data1d(m+1:m+n) = ma%data1d(1:n)
        ! append ma to the list of items
        CALL this%AppendItem(ma)
        ! free old data
        DEALLOCATE(ma%data1d,this%data1d)
        ! set 1D compound pointer to the new data
        this%data1d(1:) => data1d(1:)
      ELSE ! i.e., .NOT.ASSOCIATED(this%data1d)
        ! append ma to an empty compound => just assign the data pointer ...
        this%data1d(1:) => ma%data1d(1:)
        ! ... and append ma to the list of compound items
        CALL this%AppendItem(ma)
      END IF
      ! go through the whole compound and restore all 1D and multi-dim. pointers
      CALL this%AssignPointers()
      ! print debug info
! PRINT '(3(A,I2),I2)',"item no. ",this%num_entries, " appended to rank ",this%RANK," compound with dims ",this%DIMS(1:2)
! PRINT '(A,I6,A,I2,A,2(I3))', "  size", SIZE(this%data1d), &
!                   "  item%rank", ma%RANK, &
!                   "  item%dims", ma%DIMS(:)
! SELECT CASE(this%RANK)
! CASE(1)
!   PRINT '(A,4(I6))', "  shape4d", SHAPE(this%data4d)
! CASE(2)
!   PRINT '(A,5(I6))', "  shape5d", SHAPE(this%data5d)
! END SELECT
    END IF
  END SUBROUTINE AppendMArray
  
  !> get the first item from the list of compound elements
  FUNCTION FirstItem(this) RESULT(item)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(compound_item), POINTER :: item
    !-------------------------------------------------------------------!
    item => this%list
  END FUNCTION FirstItem

  !> get the last item from the list of compound elements
  FUNCTION LastItem(this) RESULT(item)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(compound_item), POINTER :: item
    !-------------------------------------------------------------------!
    item => this%list
    IF (.NOT.ASSOCIATED(item)) RETURN
    DO WHILE (ASSOCIATED(item%next))
      item => item%next
    END DO
  END FUNCTION LastItem

  !> get the next item from the list of compound elements
  FUNCTION NextItem(this,item) RESULT(next)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(compound_item), POINTER :: item,next
    !-------------------------------------------------------------------!
    IF (ASSOCIATED(item)) THEN
      next => item%next
    ELSE
      next => null()
    END IF
  END FUNCTION NextItem

  !> get pointer to the mesh array from a given item
  FUNCTION GetItem(this,list_entry) RESULT(ma)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(compound_item), POINTER :: list_entry
    TYPE(marray_base), POINTER :: ma
    !-------------------------------------------------------------------!
    ma => list_entry%item
  END FUNCTION GetItem

  !> append an item to the list of compound elements
  SUBROUTINE AppendItem(this,ma)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(marray_base), POINTER :: ma
    !-------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p
    !-------------------------------------------------------------------!
    IF (ASSOCIATED(ma).AND.ASSOCIATED(ma%data1d).AND.SIZE(ma%data1d).GT.0) THEN
      p => this%LastItem()
      IF (.NOT.ASSOCIATED(p)) THEN
        ! list is empty => create it
        ALLOCATE(this%list)
        p => this%list
        this%num_entries = 1 ! first entry in the list
      ELSE
        ! list exists => create a new entry at the end
        ALLOCATE(p%next)
        p => p%next
        this%num_entries = this%num_entries + 1 ! one entry added
      END IF
      ! assign item pointer and store size of 1D data
      p%item => ma
      p%extent = SIZE(ma%data1d)
      p%entry_num = this%num_entries
    END IF
  END SUBROUTINE AppendItem

  !> deconstructor of the compound
  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    !-------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p,q
    !-------------------------------------------------------------------!
    ! free list memory
    p => this%FirstItem()
    DO WHILE (ASSOCIATED(p))
      q => p
      p => this%NextItem(p)
      DEALLOCATE(q)
    END DO
    ! call deconstructor of the base class
    CALL this%marray_base%Destroy()
  END SUBROUTINE Destroy

END MODULE marray_compound_mod
