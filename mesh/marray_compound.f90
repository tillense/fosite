!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_compound.f90                                               #
!#                                                                           #
!# Copyright (C) 2018-2021                                                   #
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
    TYPE(compound_item), POINTER :: next => null()
    TYPE(marray_base), POINTER :: item => null()
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
    PROCEDURE :: AssignItemPointers
    FINAL     :: Finalize
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
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_compound::CreateMArray_compound: creating new compound"
#endif
    IF (PRESENT(n)) THEN
      ! create empty mesh array of rank 2 with first dimension set to n
      IF (new_cp%Init(n,0)) return ! immediately return if successful
    ELSE
      ! default is rank 1 mesh array
      IF (new_cp%Init(0)) return ! immediately return if successful
    END IF
#ifdef DEBUG
    PRINT *,"ERROR in marray_compound::CreateMArray: compound initialization failed"
    STOP 1
#endif
  END FUNCTION CreateMArray_compound
  
  !> \public assign pointers of different shapes to the 1D data
  FUNCTION AssignPointers(this) RESULT(success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_compound),INTENT(INOUT) :: this
    LOGICAL :: success
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_compound::AssignPointers: restoring compound pointers"
#endif
    ! set the multi-dim. pointers of the compound
    success = this%marray_base%AssignPointers()
    IF (success) success = this%AssignItemPointers()
  END FUNCTION AssignPointers

  !> assigns one compound of mesh arrays to another compound of mesh arrays
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_compound),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p,q
#ifndef __GFORTRAN__
    TYPE(marray_base), POINTER :: new_ma
    INTEGER :: err
#endif
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_compound::AssignMArray_0: compound assignment called"
#endif
    ! do basic marray assignment
    CALL this%marray_base%AssignMArray_0(ma)
    ! check if rhs is of compound class
    SELECT TYPE(src => ma)
    CLASS IS(marray_compound)
      ! check if rhs is empty
      p => src%FirstItem()
      IF (.NOT.ASSOCIATED(p)) THEN
#if DEBUG > 2
        PRINT *,"DEBUG INFO in marray_compound::AssignMArray_0: empty compound on rhs"
#endif
#ifdef DEBUG
        IF (src%num_entries.GT.0.OR.SIZE(src%data1d).GT.0) THEN
          PRINT *, "ERROR in marray_compound::AssignMArray_0: unassigned item list on rhs but compound not empty"
          STOP 1
        END IF
#endif
        this%num_entries = 0
        NULLIFY(this%list)
      ELSE ! rhs is not empty
        q => this%FirstItem()
        ! check if lhs is an empty compound
        IF (.NOT.ASSOCIATED(q)) THEN
          ! p is associated but q isn't -> lhs of assignment should be an empty compound
#if DEBUG > 2
          PRINT *,"DEBUG INFO in marray_compound::AssignMArray_0: empty compound on lhs"
#endif
#ifdef DEBUG
          IF (this%num_entries.GT.0) THEN
            PRINT *, "ERROR in marray_compound::AssignMArray_0: empty compound on lhs expected"
            STOP 1
          END IF
#endif
          ! ATTENTION: this part depends on the compiler, see comment in marray_base::AssignMArray_0
#ifdef __GFORTRAN__
          this%num_entries = src%num_entries
          this%list => src%list
#else
          DO WHILE (ASSOCIATED(p%item))
            ALLOCATE(new_ma,SOURCE=p%item,STAT=err)
            IF (err.NE.0) THEN
#ifdef DEBUG
              PRINT *,"ERROR in marray_compound::AssignMArray_0: memory allocation failed for new_ma"
              STOP 1
#else
              return
#endif
            END IF
            CALL this%AppendItem(new_ma)
            p => p%next
            IF (.NOT.ASSOCIATED(p)) EXIT
          END DO
          ! assign all item pointers
          IF (.NOT.this%AssignItemPointers()) THEN
#ifdef DEBUG
            PRINT *,"ERROR in marray_compound::AssignMArray_0: assignment of item pointers failed"
            STOP 1
#else
            RETURN
#endif
          END IF
#endif
        ELSE
          ! p and q are both associated
          IF (.NOT.(this.MATCH.ma)) THEN
#ifdef DEBUG
            PRINT *,"ERROR in marray_compound::AssignMArray_0: shape mismatch"
            STOP 1
#else
            return
#endif
          END IF
        END IF
      END IF
    CLASS DEFAULT
#ifdef DEBUG
      PRINT *, "ERROR in marray_compound::AssignMArray_0: rhs must be of class marray_compound"
      STOP 1
#else
      RETURN
#endif
    END SELECT
  END SUBROUTINE AssignMArray_0
  
  
  FUNCTION ShapesMatch(this,ma) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_compound),INTENT(IN) :: this
    CLASS(marray_base), INTENT(IN) :: ma
    !------------------------------------------------------------------------!
    TYPE(compound_item),POINTER   :: p,q
    LOGICAL                       :: res
    !------------------------------------------------------------------------!
    res = this%marray_base%ShapesMatch(ma)
    IF (.NOT.res) THEN
#if DEBUG > 1
      PRINT *,"WARNING in marray_compound::ShapesMatch: mismatch in marray_base"
#endif
      RETURN
    END IF
    SELECT TYPE(that => ma)
    CLASS IS(marray_compound)
      res = this%num_entries.EQ.that%num_entries
      IF (.NOT.res) THEN
#if DEBUG > 2
        PRINT *,"DEBUG INFO in marray_compound::ShapesMatch: number of entries do not match"
#endif
        RETURN
      END IF
      ! get pointers to the first items of the compound element lists
      p => this%FirstItem()
      q => that%FirstItem()
      ! if both compound element lists are associated
      DO WHILE (ASSOCIATED(p))
        res = res.AND.ASSOCIATED(q)
        IF (.NOT.res) RETURN ! p is associated, but q is not or extent mismatch
        res = res.AND.(p%extent.EQ.q%extent)
        p => p%next
        q => q%next
      END DO
      ! if p isn't associated then q should not be associated as well
      res = res.AND..NOT.ASSOCIATED(q)
    CLASS DEFAULT
      res = .FALSE.
    END SELECT
  END FUNCTION ShapesMatch

  FUNCTION AppendMArray(this,ma) RESULT(success)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    TYPE(marray_base), POINTER :: ma
    LOGICAL :: success
    !-------------------------------------------------------------------!
    INTEGER :: m,n,i,err
    REAL, DIMENSION(:),POINTER,CONTIGUOUS :: data1d
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_compound::AppendMArray: appending marray to compound"
#endif
    success = .FALSE.
    ! some sanity checks
    IF (.NOT.(ASSOCIATED(ma%data1d).OR.SIZE(ma%data1d).EQ.0)) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_compound::AppendMArray: input marray not associated of empty"
#else
      RETURN ! immediately return if data1d of input marray is not associated or empty
#endif
    END IF
    IF (.NOT.ASSOCIATED(this%data1d)) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_compound::AppendMArray: compound uninitialized"
#else
      RETURN ! immediately return if data1d of compound is not associated
#endif
    END IF

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
#ifdef DEBUG
        PRINT *,"ERROR in marray_compound::AppendMArray: this%dims(1) != ma%dims(1)"
#endif
        RETURN
      ELSE
        ! this%DIMS(1) remains the same
        this%DIMS(2) = this%DIMS(2) + ma%DIMS(2) ! extend the last dimension of the compound;
          ! so far only scalars and vectors can be added to rank 2 compounds
      END IF
    CASE DEFAULT
      ! this should not happen
#ifdef DEBUG
      PRINT *,"ERROR in marray_compound::AppendMArray: rank of compound marrays should be either 1 or 2"
#endif
      RETURN
    END SELECT

    ! determine size of compound data1d array
    m = SIZE(this%data1d(:))
    IF (m.EQ.0) THEN
      ! initialized but empty compound
      DEALLOCATE(this%data1d)
      ! append ma to an empty compound => just assign the data pointer ...
      this%data1d(1:) => ma%data1d(1:)
      ! ... and append ma to the list of compound items
      CALL this%AppendItem(ma)
    ELSE ! m > 0 -> compound has at least one component
      ! determine size of new component
      n = SIZE(ma%data1d(:))
      ! allocate new contiguous data block for old+new data and for meta data
      ALLOCATE(data1d(m+n),STAT=err)
      IF (err.NE.0) THEN
#ifdef DEBUG
        PRINT *,"ERROR in marray_compound::AppendMArray: memory allocation failed"
#else
        RETURN ! immediatly return on failure
#endif
      END IF
      ! copy the data, new data is appended to the old compound
      DO CONCURRENT (i=1:m)
        data1d(i) = this%data1d(i)
      END DO
      DO CONCURRENT (i=1:n)
        data1d(m+i) = ma%data1d(i)
      END DO
      ! append ma to the list of items
      CALL this%AppendItem(ma)
      ! free old data
      DEALLOCATE(ma%data1d,this%data1d)
      ! set 1D compound pointer to the new data
      this%data1d(1:) => data1d(1:)
    END IF
    ! go through the whole compound and restore all 1D and multi-dim. pointers
    IF (.NOT.this%AssignPointers()) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_compound::AppendMArray: pointer assignment failed"
#else
      RETURN ! immediatly return on failure
#endif
    END IF
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_compound::AppendMArray"
    PRINT '(3(A,I2),I2)',"  item no. ",this%num_entries, " appended to rank ",this%RANK, &
                        " compound with dims ",this%DIMS(1:2)
    PRINT '(A,I10,A,I2,A,2(I3))', "    size ", SIZE(this%data1d), &
                      "    item%rank", ma%RANK, &
                      "    item%dims", ma%DIMS(:)
    SELECT CASE(this%RANK)
    CASE(1)
      PRINT '(A,4(I6))', "    shape4d", SHAPE(this%data4d)
    CASE(2)
      PRINT '(A,5(I6))', "    shape5d", SHAPE(this%data5d)
    END SELECT
#endif
    success = .TRUE.
  END FUNCTION AppendMArray
  
  !> get the first item from the list of compound elements
  FUNCTION FirstItem(this) RESULT(item)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound), INTENT(IN) :: this
    TYPE(compound_item), POINTER :: item
    !-------------------------------------------------------------------!
    item => this%list
  END FUNCTION FirstItem

  !> get the last item from the list of compound elements
  FUNCTION LastItem(this) RESULT(item)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound), INTENT(IN) :: this
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
    CLASS(marray_compound), INTENT(IN) :: this
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
    INTEGER :: err
    !-------------------------------------------------------------------!
    IF (ASSOCIATED(ma).AND.ASSOCIATED(ma%data1d).AND.SIZE(ma%data1d).GT.0) THEN
      p => this%LastItem()
      IF (.NOT.ASSOCIATED(p)) THEN
        ! list is empty => create it
#if DEBUG > 2
        PRINT *,"DEBUG INFO in marray_compound::AppendItem: creating new list"
#endif
        ALLOCATE(this%list,STAT=err)
        IF (err.NE.0) THEN
#ifdef DEBUG
          PRINT *,"ERROR in marray_compound::AppendItem: memory allocation failed for new list"
#endif
          return
        END IF
        p => this%list
        this%num_entries = 1 ! first entry in the list
      ELSE
        ! list exists => create a new entry at the end
#if DEBUG > 2
        PRINT *,"DEBUG INFO in marray_compound::AppendItem: appending item to list of elements"
#endif
        ALLOCATE(p%next,STAT=err)
        IF (err.NE.0) THEN
#ifdef DEBUG
          PRINT *,"ERROR in marray_compound::AppendItem: memory allocation failed for new entry"
#endif
          return
        END IF
        p => p%next
        this%num_entries = this%num_entries + 1 ! one entry added
      END IF
      ! assign item pointer and store size of 1D data
      p%item => ma
      p%extent = SIZE(ma%data1d)
      p%entry_num = this%num_entries
      NULLIFY(p%next)
    END IF
  END SUBROUTINE AppendItem

  FUNCTION AssignItemPointers(this) RESULT(success)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_compound) :: this
    LOGICAL :: success
    !-------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p
    INTEGER :: m,n
    !------------------------------------------------------------------------!
    ! go through the list of items and set the 1D
    ! and multi-dim. pointers for each item
    p => this%FirstItem()
    m = 0
    n = 0
    success = .TRUE.
    DO WHILE (success.AND.ASSOCIATED(p))
#if DEBUG > 2
      PRINT '(A,I2)'," DEBUG INFO in marray_compound::AssignItemPointers: restoring entry no. ", p%entry_num
#endif
      m = m + n
      n = p%extent
      ! point to the 1D data segment in the new array
      p%item%data1d(1:n) => this%data1d(m+1:m+n)
      ! restore all multi-dim. pointers
      success = p%item%AssignPointers()
      p => this%NextItem(p)
    END DO
  END FUNCTION AssignItemPointers

  !> destructor of compounds
  !! ATTENTION: the data array itself is deallocated by the inherited finalizer
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_compound), INTENT(INOUT) :: this
    !-------------------------------------------------------------------!
    TYPE(compound_item), POINTER :: p,q
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_compound::Finalize: deallocating compound components"
#endif
    ! free list memory
    p => this%FirstItem()
    DO WHILE (ASSOCIATED(p))
#if DEBUG > 2
      PRINT '(A,I2)'," DEBUG INFO in marray_compound::Finalize: deleting entry no. ", p%entry_num
#endif
      q => p%next
      IF (ASSOCIATED(p%item)) THEN
#if DEBUG > 2
        PRINT '(A,I2)'," DEBUG INFO in marray_compound::Finalize: deallocating item data"
#endif
        ! skip deallocation of the compound element data when invoking the finalizer of p%item
        IF (ASSOCIATED(p%item%data1d)) NULLIFY(p%item%data1d)
        DEALLOCATE(p%item)
      END IF
      DEALLOCATE(p)
      p => q
    END DO
    this%num_entries = 0
    NULLIFY(this%list)
  END SUBROUTINE Finalize

END MODULE marray_compound_mod
