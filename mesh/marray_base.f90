!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_base.f90                                                   #
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
!> \defgroup marray Marray
!! \{
!! \brief Module class for mesh arrays
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief base class for mesh arrays
!!
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> \name Public Attributes
  !! #### mesh indices
  INTEGER, SAVE          :: IGMIN,IGMAX                   !< 1st dim
  INTEGER, SAVE          :: JGMIN,JGMAX                   !< 2nd dim
  INTEGER, SAVE          :: KGMIN,KGMAX                   !< 3rd dim
  INTEGER, SAVE          :: INUM,JNUM,KNUM                !< array sizes
  LOGICAL, SAVE          :: IDX_INIT = .FALSE.            !< init status
  !> type for selecting parts of an marray
  TYPE selection_base
     !> \name Variables
     INTEGER           :: imin,imax                  !< selection in x-direction
     INTEGER           :: jmin,jmax                  !< selection in y-direction
     INTEGER           :: kmin,kmax                  !< selection in z-direction
     LOGICAL, POINTER, CONTIGUOUS :: mask1d(:), &    !< 1d selection mask
                                     mask2d(:,:), &  !< 2d selection mask
                                     mask3d(:,:,:)   !< 3d selection mask
  CONTAINS
     !> \name Methods
     PROCEDURE :: Init => Init_selection
     PROCEDURE :: AssignSelection
     PROCEDURE :: Cuboid
     PROCEDURE :: Everything
     PROCEDURE :: Destroy_selection
     GENERIC :: ASSIGNMENT (=) => AssignSelection
     GENERIC :: Destroy => Destroy_selection
     FINAL   :: Destructor_selection
  END TYPE selection_base
  !> basic mesh array class
  TYPE :: marray_base
    INTEGER       :: RANK = -1,DIMS(2) = 0
    REAL, POINTER, CONTIGUOUS :: data1d(:) => null()
    REAL, POINTER, CONTIGUOUS :: data2d(:,:) => null()
    REAL, POINTER, CONTIGUOUS :: data3d(:,:,:) => null()
    REAL, POINTER, CONTIGUOUS :: data4d(:,:,:,:) => null()
    REAL, POINTER, CONTIGUOUS :: data5d(:,:,:,:,:) => null()
    CONTAINS
    PROCEDURE :: Init
    PROCEDURE :: AssignPointers
    PROCEDURE :: RemapBounds_0
    PROCEDURE :: RemapBounds_1
    PROCEDURE :: RemapBounds_2
    GENERIC   :: RemapBounds => RemapBounds_0, RemapBounds_1, RemapBounds_2
    PROCEDURE :: AssignMArray_0
    PROCEDURE :: AssignMArray_1
    PROCEDURE :: AssignMArray_2
    PROCEDURE :: AssignMArray_3
    PROCEDURE :: AssignMArray_4
    PROCEDURE :: AssignMArray_5
    GENERIC   :: ASSIGNMENT (=) => AssignMArray_0 , AssignMArray_1, AssignMArray_2, &
                                AssignMArray_3 , AssignMArray_4, AssignMArray_5
    PROCEDURE :: ShapesMatch
    GENERIC   :: OPERATOR(.MATCH.) => ShapesMatch
    PROCEDURE :: AddMArray_0
    PROCEDURE :: AddMArray_1
    PROCEDURE :: AddMArray_2
    PROCEDURE :: AddMArray_3
    PROCEDURE :: AddMArray_4
    PROCEDURE :: AddMArray_5
    GENERIC   :: OPERATOR (+) => AddMArray_0, AddMArray_1, AddMArray_2, &
                                 AddMArray_3, AddMArray_4, AddMArray_5
    PROCEDURE :: MultMArray_0
    PROCEDURE :: MultMArray_1
    PROCEDURE :: MultMArray_2
    PROCEDURE :: MultMArray_3
    PROCEDURE :: MultMArray_4
    PROCEDURE :: MultMArray_5
    GENERIC   :: OPERATOR (*) => MultMArray_0, MultMArray_1, MultMArray_2, &
                                 MultMArray_3, MultMArray_4, MultMArray_5
    PROCEDURE :: Destroy
    FINAL     :: Destructor
  END TYPE
  INTERFACE marray_base
    MODULE PROCEDURE CreateMArray
  END INTERFACE
  INTERFACE selection_base
    MODULE PROCEDURE CreateSelection
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_base, &
    selection_base, &
    InitMeshProperties, &
    CloseMeshProperties

  CONTAINS

  !> constructor for mesh arrays
  FUNCTION CreateMArray(m,n) RESULT(new_ma)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_base) :: new_ma
    INTEGER, OPTIONAL, INTENT(IN) :: m,n
    !-------------------------------------------------------------------!
    INTEGER :: err
    !-------------------------------------------------------------------!
    IF (new_ma%Init(m,n)) return ! immediately return if successful
    ! something went wrong
#ifdef DEBUG > 0
    PRINT *,"ERROR in marray_base::CreateMArray: marray initialization failed"
    STOP 1
#endif
  END FUNCTION CreateMArray

  !> basic initialization of mesh array class
  FUNCTION Init(this,m,n) RESULT(success)
    !-------------------------------------------------------------------!
    CLASS(marray_base), INTENT(INOUT) :: this
    INTEGER, OPTIONAL, INTENT(IN) :: m,n
    !-------------------------------------------------------------------!
    LOGICAL :: success
    INTEGER :: err
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::Init: marray initialization"
#endif
    success = .FALSE.
    IF (.NOT.IDX_INIT) return ! with success == .false.
    ! check input parameter
    ! ATTENTION: m=0 and n=0 is permitted to allow for creation of empty
    ! mesh_arrays which is required, e.g., during initialization of compounds
    IF (PRESENT(m)) THEN
      IF (m.LT.0) THEN
#ifdef DEBUG
        PRINT *,"ERROR in marray_base::Init: 1st dimension of mesh array should be >= 0"
#endif
        return ! with success == .false.
      END IF
      this%DIMS(1) = m
      IF (PRESENT(n)) THEN
        IF (n.LT.0) THEN
#ifdef DEBUG
          PRINT *,"ERROR in marray_base::Init: 2nd dimension of mesh array should be >= 0"
#endif
          return ! with success == .false.
        END IF
        this%DIMS(2) = n
        this%RANK = 2
      ELSE
        this%DIMS(2) = 1
        this%RANK = 1
      END IF
    ELSE
      this%DIMS(:) = 1
      this%RANK = 0
    END IF

    IF (.NOT.ASSOCIATED(this%data1d)) THEN
#if DEBUG > 2
      PRINT '(A,I2,A,2(I4))',"  creating marray with rank ",this%RANK," and dimensions ",this%DIMS(1:2)
#endif
      ! allocate memory for data1d array
      ALLOCATE(this%data1d(INUM*JNUM*KNUM*this%DIMS(1)*this%DIMS(2)),STAT=err)
      IF (err.NE.0) THEN
#ifdef DEBUG
        PRINT *,"ERROR in marray_base::Init: memory allocation failed for data1d array"
#endif
        return ! with success == .false.
      END IF
    ELSE
      IF (SIZE(this%data1d).NE.INUM*JNUM*KNUM*this%DIMS(1)*this%DIMS(2)) THEN
#ifdef DEBUG
        PRINT *,"ERROR in marray_base::Init: data1d array size mismatch"
#endif
        return ! with success == .false.
      END IF
    END IF
    ! assign the 2d,3d,... pointers
    IF (SIZE(this%data1d).GT.0) THEN ! do not try pointer assignment with zero-size data1d array
      IF (.NOT.this%AssignPointers()) return ! pointer assignment failed -> return immediately
    END IF
    ! report success
    success=.TRUE.
  END FUNCTION Init

  !> sets global mesh properties
  !!
  !! This subroutine should be called only once in MeshInit.
  SUBROUTINE InitMeshProperties(igmin_,igmax_,jgmin_,jgmax_,kgmin_,kgmax_)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: igmin_,igmax_,jgmin_,jgmax_,kgmin_,kgmax_
    !-------------------------------------------------------------------!
    IF (.NOT.IDX_INIT) THEN
      IGMIN = igmin_
      IGMAX = igmax_
      JGMIN = jgmin_
      JGMAX = jgmax_
      KGMIN = kgmin_
      KGMAX = kgmax_
      INUM = IGMAX-IGMIN+1
      JNUM = JGMAX-JGMIN+1
      KNUM = KGMAX-KGMIN+1
      IDX_INIT = .TRUE.
    END IF
  END SUBROUTINE InitMeshProperties
  
  !> unsets global mesh properties
  !!
  !! This subroutine should be called only once in MeshClose.
  SUBROUTINE CloseMeshProperties
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    IF (IDX_INIT) THEN
      IDX_INIT = .FALSE.
    END IF
  END SUBROUTINE CloseMeshProperties
  
  !> \public assign pointers of different shapes to the 1D data
  FUNCTION AssignPointers(this) RESULT(success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    LOGICAL :: success
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignPointers: assign 2d,3d,... pointers"
#endif
    success=.FALSE.
    IF (.NOT.ASSOCIATED(this%data1d).OR.SIZE(this%data1d).EQ.0) return
    ! assign pointers depending on rank
    SELECT CASE(this%RANK)
    CASE(0)
      this%data2d(1:INUM*JNUM,KGMIN:KGMAX) => this%data1d
      this%data3d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX) => this%data1d
      this%data4d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX,1:1) => this%data1d
      this%data5d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX,1:1,1:1) => this%data1d
    CASE(1)
      this%data2d(1:INUM*JNUM*KNUM,1:this%DIMS(1)) => this%data1d
      this%data3d(1:INUM*JNUM,KGMIN:KGMAX,1:this%DIMS(1)) => this%data1d
      this%data4d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX,1:this%DIMS(1)) => this%data1d
      this%data5d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX,1:this%DIMS(1),1:1) => this%data1d
    CASE(2)
      this%data2d(1:INUM*JNUM*KNUM*this%DIMS(1),1:this%DIMS(2)) => this%data1d
      this%data3d(1:INUM*JNUM*KNUM,1:this%DIMS(1),1:this%DIMS(2)) => this%data1d
      this%data4d(1:INUM*JNUM,KGMIN:KGMAX,1:this%DIMS(1),1:this%DIMS(2)) => this%data1d
      this%data5d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX,1:this%DIMS(1),1:this%DIMS(2)) => this%data1d
    CASE DEFAULT
      return ! wrong rank
    END SELECT
    ! report success
    success=.TRUE.
  END FUNCTION AssignPointers

  !> \public remap lower bounds in the first 3 dimensions of rank 0 mesh arrays
  !!
  !! This is a short hack to obviate a restriction in the generation of
  !! subarray pointers. The indices of subarrays usually start with a lower
  !! bound of 1, but Fosite requires that all mesh data arrays start with
  !! lower bounds of IGMIN, JGMIN and KGMIN, which are not equal to 1 in general.
  FUNCTION RemapBounds_0(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base)                :: this !< \param [in,out] this
    REAL, DIMENSION(IGMIN:,JGMIN:,KGMIN:), TARGET :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_0

  !> \public remap lower bounds in the first 3 dimensions of rank 1 mesh arrays
  FUNCTION RemapBounds_1(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base)                  :: this !< \param [in,out] this
    REAL, DIMENSION(IGMIN:,JGMIN:,KGMIN:,:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                        :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_1

  !> \public remap lower bounds in the first 3 dimensions of rank 2 mesh arrays
  FUNCTION RemapBounds_2(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base) :: this      !< \param [in,out] 
    REAL, DIMENSION(IGMIN:,JGMIN:,KGMIN:,:,:), TARGET &
                       :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_2

  !> assigns one mesh array to another mesh array
#ifndef DEBUG
  PURE &
#endif
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    LOGICAL :: success = .FALSE.
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignMArray_0: marray assignment"
#endif
    IF (.NOT.ASSOCIATED(ma%data1d).OR.ma%rank.LT.0) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::AssignMArray_0: rhs of assignment not initialized"
#endif
      return
    END IF

    IF (.NOT.ASSOCIATED(this%data1d)) THEN
      ! lhs of assignment uninitialized -> initialize new mesh array
      ! ATTENTION: finalization of derived types works different for
      !   GNU Fortran, hence to prevent memory leaks, one has to point
      !   the data1d array of the lhs (this%data1d) to the already associated
      !   data1d array of the rhs (ma%data1d).
      !   Other compilers, e.g., ifort (intel) & nfort (NEC) require generation
      !   of a new marray with data1d array which is destroyed on exit.
#ifdef __GFORTRAN__
      this%RANK = ma%RANK
      this%DIMS(:) = ma%DIMS(:)
      this%data1d => ma%data1d
      IF (this%AssignPointers()) return ! pointer assignment ok -> return immediately
#else
      SELECT CASE(ma%RANK)
      CASE(0)
        success = this%Init()
      CASE(1)
        success = this%Init(m=ma%DIMS(1))
      CASE(2)
        success = this%Init(m=ma%DIMS(1),n=ma%DIMS(2))
      CASE DEFAULT
#ifdef DEBUG
        PRINT *,"ERROR in marray_base::AssignMArray_0: rank > 2 currently not supported"
#endif
        return
      END SELECT
      IF (.NOT.success) THEN
#ifdef DEBUG
        PRINT *,"ERROR in marray_base::AssignMArray_0: marray initialization failed"
#endif
        return
      END IF
#endif
    END IF

    IF (.NOT.(this.MATCH.ma)) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::AssignMArray_0: shape mismatch"
#endif
      return
    END IF
    IF (SIZE(this%data1d).NE.SIZE(ma%data1d)) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::AssignMArray_0: size mismatch of data1d array"
#endif
      return
    END IF
    ! copy data
    this%data1d(:) = ma%data1d(:)
    IF (this%AssignPointers()) return ! pointer assignment ok -> return immediately
    ! something bad happend
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::AssignMArray_0: final pointer reassignment failed for lhs"
#endif
  END SUBROUTINE AssignMArray_0
  
  !> assign 1D fortran array to mesh array
#ifndef DEBUG
  PURE &
#endif
  SUBROUTINE AssignMArray_1(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(INUM*JNUM*KNUM*this%DIMS(1)*this%DIMS(2)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignMArray_1: assigning 1D Fortran array"
#endif
    this%data1d(:) = a(:)
  END SUBROUTINE AssignMArray_1

  !> assign 2D fortran array to mesh array
#ifndef DEBUG
  PURE &
#endif
  SUBROUTINE AssignMArray_2(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data2d,1),SIZE(this%data2d,2)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignMArray_2: assigning 2D Fortran array"
#endif
    this%data2d(:,:) = a(:,:)
  END SUBROUTINE AssignMArray_2

  !> assign 3D fortran array to mesh array
#ifndef DEBUG
  PURE &
#endif
  SUBROUTINE AssignMArray_3(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data3d,1),SIZE(this%data3d,2),SIZE(this%data3d,3)), &
                          INTENT(IN) :: a
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignMArray_3: assigning 3D Fortran array"
#endif
    this%data3d(:,:,:) = a(:,:,:)
  END SUBROUTINE AssignMArray_3

  !> assign 4D fortran array to mesh array
#ifndef DEBUG
  PURE &
#endif
  SUBROUTINE AssignMArray_4(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data4d,1),SIZE(this%data4d,2),SIZE(this%data4d,3), &
                SIZE(this%data4d,4)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignMArray_1: assigning 4D Fortran array"
#endif
    this%data4d(:,:,:,:) = a(:,:,:,:)
  END SUBROUTINE AssignMArray_4

  !> assign 5D fortran array to mesh array
#ifndef DEBUG
  PURE &
#endif
  SUBROUTINE AssignMArray_5(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data5d,1),SIZE(this%data5d,2),SIZE(this%data5d,3), &
                SIZE(this%data5d,4),SIZE(this%data5d,5)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
    this%data5d(:,:,:,:,:) = a(:,:,:,:,:)
  END SUBROUTINE AssignMArray_5

  PURE FUNCTION ShapesMatch(this,ma) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: this, ma
    !------------------------------------------------------------------------!
    LOGICAL                       :: res
    !------------------------------------------------------------------------!
    res = .FALSE.
    IF (this%rank.NE.ma%rank) return
    IF (SIZE(this%dims).NE.SIZE(ma%dims)) return
    IF (.NOT.ALL(this%dims(:).EQ.ma%dims(:))) return
    res = .TRUE.
  END FUNCTION ShapesMatch

  !> add 2 mesh arrays
  PURE FUNCTION AddMArray_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data1d)) :: c
    !------------------------------------------------------------------------!
    IF (SIZE(a%data1d).EQ.SIZE(b%data1d)) &
        c(:) = a%data1d(:) + b%data1d(:)
  END FUNCTION AddMArray_0

  !> add 1D fortran array and mesh array
  PURE FUNCTION AddMArray_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data1d)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data1d)) :: c
    !------------------------------------------------------------------------!
    c(:) = a%data1d(:) + b(:)
  END FUNCTION AddMArray_1

  !> add 2D fortran array and mesh array
  PURE FUNCTION AddMArray_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data2d,1),SIZE(a%data2d,2)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data2d,1),SIZE(a%data2d,2)) :: c
    !------------------------------------------------------------------------!
    c(:,:) = a%data2d(:,:) + b(:,:)
  END FUNCTION AddMArray_2

  !> add 3D fortran array and mesh array
  PURE FUNCTION AddMArray_3(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data3d,1),SIZE(a%data3d,2),SIZE(a%data3d,3)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data3d,1),SIZE(a%data3d,2),SIZE(a%data3d,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data3d(:,:,:) + b(:,:,:)
  END FUNCTION AddMArray_3

  !> add 4D fortran array and mesh array
  PURE FUNCTION AddMArray_4(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data4d,1),SIZE(a%data4d,2),SIZE(a%data4d,3), &
          SIZE(a%data4d,4)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data4d,1),SIZE(a%data4d,2),SIZE(a%data4d,3), &
          SIZE(a%data4d,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data4d(:,:,:,:) + b(:,:,:,:)
  END FUNCTION AddMArray_4

  !> add 5D fortran array and mesh array
  PURE FUNCTION AddMArray_5(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data5d,1),SIZE(a%data5d,2),SIZE(a%data5d,3), &
          SIZE(a%data5d,4),SIZE(a%data5d,5)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data5d,1),SIZE(a%data5d,2),SIZE(a%data5d,3), &
          SIZE(a%data5d,4),SIZE(a%data5d,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a%data5d(:,:,:,:,:) + b(:,:,:,:,:)
  END FUNCTION AddMArray_5

  !> multiply 2 mesh arrays
#ifndef DEBUG
  PURE &
#endif
  FUNCTION MultMArray_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data1d)) :: c
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::MultMArray_0: multiply 2 marrays"
#endif
#ifdef DEBUG
    IF (.NOT.(a.MATCH.b)) THEN
      PRINT *,"ERROR in marray_base::MultMArray_0: shape mismatch"
    END IF
#endif
    c(:) = a%data1d(:) * b%data1d(:)
  END FUNCTION MultMArray_0

  !> multiply 1D fortran array and mesh arrays
#ifndef DEBUG
  PURE &
#endif
  FUNCTION MultMArray_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data1d)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data1d)) :: c
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::MultMArray_0: multiply marray with 1d Fortran array"
#endif
    c(:) = a%data1d(:) * b(:)
  END FUNCTION MultMArray_1

  !> multiply 2D fortran array and mesh arrays
#ifndef DEBUG
  PURE &
#endif
  FUNCTION MultMArray_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data2d,1),SIZE(a%data2d,2)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data2d,1),SIZE(a%data2d,2)) :: c
    !------------------------------------------------------------------------!
    c(:,:) = a%data2d(:,:) * b(:,:)
  END FUNCTION MultMArray_2

  !> multiply 3D fortran array and mesh arrays
#ifndef DEBUG
  PURE &
#endif
  FUNCTION MultMArray_3(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data3d,1),SIZE(a%data3d,2),SIZE(a%data3d,3)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data3d,1),SIZE(a%data3d,2),SIZE(a%data3d,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data3d(:,:,:) * b(:,:,:)
  END FUNCTION MultMArray_3

  !> multiply 4D fortran array and mesh arrays
#ifndef DEBUG
  PURE &
#endif
  FUNCTION MultMArray_4(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data4d,1),SIZE(a%data4d,2),SIZE(a%data4d,3), &
          SIZE(a%data4d,4)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data4d,1),SIZE(a%data4d,2),SIZE(a%data4d,3), &
          SIZE(a%data4d,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data4d(:,:,:,:) * b(:,:,:,:)
  END FUNCTION MultMArray_4

  !> multiply 5D fortran array and mesh arrays
#ifndef DEBUG
  PURE &
#endif
  FUNCTION MultMArray_5(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data5d,1),SIZE(a%data5d,2),SIZE(a%data5d,3), &
          SIZE(a%data5d,4),SIZE(a%data5d,5)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data5d,1),SIZE(a%data5d,2),SIZE(a%data5d,3), &
          SIZE(a%data5d,4),SIZE(a%data5d,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a%data5d(:,:,:,:,:) * b(:,:,:,:,:)
  END FUNCTION MultMArray_5

  !> polymorphic destructor of all mesh_array classes
  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_base) :: this
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::Destroy: deallocating data1d"
#endif
    IF (ASSOCIATED(this%data1d)) DEALLOCATE(this%data1d)
    NULLIFY(this%data1d,this%data2d,this%data3d,this%data4d)
    this%rank    =-1
    this%dims(:) = 0
  END SUBROUTINE Destroy

  !> actual destructor of mesh_array - this is called automatically if
  !! deallocate is invoked
  SUBROUTINE Destructor(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_base) :: this
    !-------------------------------------------------------------------!
    CALL this%Destroy()
  END SUBROUTINE Destructor

  FUNCTION CreateSelection(idx) RESULT(new_sel)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(selection_base) :: new_sel
    INTEGER, OPTIONAL    :: idx(6)
    !-------------------------------------------------------------------!
    INTEGER :: err
    !-------------------------------------------------------------------!
    IF (new_sel%Init(idx)) return ! immediately return if successful
#ifdef DEBUG
    PRINT *,"ERROR in selection_base::CreateSelection: initialization failed"
    STOP 1
#endif
  END FUNCTION CreateSelection

  !> basic initialization of selection
  FUNCTION Init_selection(this,idx) RESULT(success)
    !-------------------------------------------------------------------!
    CLASS(selection_base), INTENT(INOUT) :: this
    INTEGER, OPTIONAL    :: idx(6)
    !-------------------------------------------------------------------!
    LOGICAL :: success
    INTEGER :: err
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::Init_selection: selection initialization"
#endif
    success = .FALSE.
    IF (.NOT.IDX_INIT) return ! with success == .false.
    ! allocate 1D mask array
    ALLOCATE(this%mask1d(1:INUM*JNUM*KNUM),STAT=err)
    IF (err.NE.0) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::Init_selection: memory allocation failed"
#endif
      return ! with success == .false.
    END IF
    ! set 2D & 3D pointers
    this%mask2d(1:INUM*JNUM,KGMIN:KGMAX) => this%mask1d
    this%mask3d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX) => this%mask1d
    IF (PRESENT(idx)) THEN
      CALL this%Cuboid(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
    ELSE
      CALL this%Everything()
    END IF
    ! report success
    success=.TRUE.
  END FUNCTION Init_selection

  !> assigns one selection to another selection
  SUBROUTINE AssignSelection(this,sel)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(selection_base),INTENT(INOUT) :: this
    CLASS(selection_base),INTENT(IN)    :: sel
    !------------------------------------------------------------------------!
    LOGICAL :: success = .FALSE.
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_base::AssignSelection: selection assignment"
#endif
    IF (.NOT.ASSOCIATED(sel%mask1d)) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::AssignSelection: rhs of assignment not initialized"
#endif
      return
    END IF

    IF (.NOT.ASSOCIATED(this%mask1d)) THEN
      ! lhs of assignment uninitialized -> initialize new selection
      ! ATTENTION: finalization of derived types works different for
      !   GNU Fortran, hence to prevent memory leaks, one has to point
      !   the mask1d array of the lhs (this%mask1d) to the already associated
      !   mask1d array of the rhs (ma%mask1d).
      !   Other compilers, e.g., ifort (intel) & nfort (NEC) require generation
      !   of a new selection with mask1d array which is destroyed on exit.
#ifdef __GFORTRAN__
      this%imin = sel%imin
      this%imax = sel%imax
      this%jmin = sel%jmin
      this%jmax = sel%jmax
      this%kmin = sel%kmin
      this%kmax = sel%kmax
      this%mask1d => sel%mask1d
      ! set 2D & 3D pointers
      this%mask2d(1:INUM*JNUM,KGMIN:KGMAX) => this%mask1d
      this%mask3d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX) => this%mask1d
      ! immediately return
      return
#else
      IF (.NOT.this%Init((/sel%imin,sel%imax,sel%jmin,sel%jmax,sel%kmin,sel%kmax/))) THEN
#ifdef DEBUG
        PRINT *,"ERROR in marray_base::AssignSelection: initialization failed"
#endif
        return
      END IF
#endif
    ELSE
      this%imin = sel%imin
      this%imax = sel%imax
      this%jmin = sel%jmin
      this%jmax = sel%jmax
      this%kmin = sel%kmin
      this%kmax = sel%kmax
    END IF
    IF (.NOT.SIZE(this%mask1d).EQ.SIZE(sel%mask1d)) THEN
#ifdef DEBUG
      PRINT *,"ERROR in marray_base::AssignSelection: size mismatch"
#endif
      return
    END IF
    ! copy data
    this%mask1d(:) = sel%mask1d(:)
  END SUBROUTINE AssignSelection

  SUBROUTINE Cuboid(this,imin,imax,jmin,jmax,kmin,kmax)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(selection_base) :: this
    INTEGER, INTENT(IN)   :: imin,imax,jmin,jmax,kmin,kmax
    !-------------------------------------------------------------------!
    INTEGER :: i,j,k
    !-------------------------------------------------------------------!
    this%imin = MIN(MAX(imin,IGMIN),IGMAX)
    this%imax = MAX(MIN(imax,IGMAX),IGMIN)
    this%jmin = MIN(MAX(jmin,JGMIN),JGMAX)
    this%jmax = MAX(MIN(jmax,JGMAX),JGMIN)
    this%kmin = MIN(MAX(kmin,KGMIN),KGMAX)
    this%kmax = MAX(MIN(kmax,KGMAX),KGMIN)
    DO k=KGMIN,KGMAX
      DO j=JGMIN,JGMAX
        DO i=IGMIN,IGMAX
          IF ((i.GE.this%imin.AND.i.LE.this%imax).AND. &
              (j.GE.this%jmin.AND.j.LE.this%jmax).AND. &
              (k.GE.this%kmin.AND.k.LE.this%kmax)) THEN
            this%mask3d(i,j,k) = .TRUE.
          ELSE
            this%mask3d(i,j,k) = .FALSE.
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE Cuboid

  SUBROUTINE Everything(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(selection_base) :: this
    !-------------------------------------------------------------------!
    this%imin = IGMIN
    this%imax = IGMAX
    this%jmin = JGMIN
    this%jmax = JGMAX
    this%kmin = KGMIN
    this%imax = KGMAX
    this%mask1d(:) = .TRUE.
  END SUBROUTINE Everything

  !> destructor of all selection classes
  SUBROUTINE Destroy_selection(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(selection_base) :: this
    !-------------------------------------------------------------------!
    IF (ASSOCIATED(this%mask1d)) DEALLOCATE(this%mask1d)
    NULLIFY(this%mask2d,this%mask3d)
  END SUBROUTINE Destroy_selection

  !> actual destructor of selection_base
  SUBROUTINE Destructor_selection(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(selection_base) :: this
    !-------------------------------------------------------------------!
    CALL this%Destroy_selection()
  END SUBROUTINE Destructor_selection

END MODULE marray_base_mod
