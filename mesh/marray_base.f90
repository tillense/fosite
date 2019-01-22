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
!> \addtogroup marray
!!
!! \author Tobias Illenseer
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
     LOGICAL, POINTER  :: mask1d(:)                  !< 1d selection mask
     LOGICAL, POINTER  :: mask2d(:,:)                !< 2d selection mask
     LOGICAL, POINTER, CONTIGUOUS :: mask3d(:,:,:)   !< 3d selection mask
  CONTAINS
     !> \name Methods
     PROCEDURE :: Cuboid
     PROCEDURE :: Everything
     PROCEDURE :: Destroy_selection
     GENERIC :: Destroy => Destroy_selection
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
    IF (IDX_INIT) THEN
      IF (PRESENT(m)) THEN
        new_ma%DIMS(1) = m
        IF (PRESENT(n)) THEN
          new_ma%DIMS(2) = n
          new_ma%RANK = 2
        ELSE
          new_ma%DIMS(2) = 1
          new_ma%RANK = 1
        END IF
      ELSE
        new_ma%DIMS(:) = 1
        new_ma%RANK = 0
      END IF
    END IF
  END FUNCTION CreateMArray
  
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
  SUBROUTINE AssignPointers(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(this%data1d).AND.SIZE(this%data1d).GT.0) THEN
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
      END SELECT
    END IF
  END SUBROUTINE AssignPointers

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
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    INTEGER :: err
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(ma%data1d)) THEN
      ! copy meta data
      this%RANK    = ma%RANK
      this%DIMS(:) = ma%DIMS(:)
      IF (ASSOCIATED(this%data1d)) THEN
        IF (SIZE(this%data1d).EQ.0) THEN
          ! this%data1d allocated with size zero
          DEALLOCATE(this%data1d)
          ! allocate and assign data
          ALLOCATE(this%data1d(SIZE(ma%data1d,1)),SOURCE=ma%data1d,STAT=err)
          IF (err.NE.0) THEN
            PRINT *,"ERROR in marray_base::AssignMArray_0: memory allocation failed"
            STOP 1
          END IF
        ELSE IF (SIZE(this%data1d).EQ.SIZE(ma%data1d)) THEN
          ! just copy the data
          this%data1d(:) = ma%data1d(:)
        ELSE
          !> \todo improve error handling for mesh arrays
          PRINT *,"ERROR in marray_base::AssignMArray_0: size of input and output do not match"
          STOP 1
        END IF
      ELSE
        ! allocate and copy data
        ALLOCATE(this%data1d(SIZE(ma%data1d,1)),SOURCE=ma%data1d,STAT=err)
        IF (err.NE.0) THEN
          PRINT *,"ERROR in marray_base::AssignMArray_0: memory allocation failed"
          STOP 1
        END IF
      END IF
      CALL this%AssignPointers()
    ELSE IF ((.NOT.ASSOCIATED(ma%data1d)).AND.(ma%RANK.GE.0)) THEN
      this%RANK    = ma%RANK
      this%DIMS(:) = ma%DIMS(:)
      ! just allocate data of newly created marray without copying
      ALLOCATE(this%data1d(INUM*JNUM*KNUM*this%DIMS(1)*this%DIMS(2)),STAT=err)
      IF (err.NE.0) THEN
        PRINT *,"ERROR in marray_base::AssignMArray_0: memory allocation failed"
        STOP 1
      END IF
      CALL this%AssignPointers()
    ELSE IF (ASSOCIATED(this%data1d)) THEN
      ! both of same type ma%data1d not associated, but this%data1d associated
      PRINT *,"ERROR in marray_base::AssignMArray_0: ma%data1d not associated but this%data1d is"
      STOP 1
    ELSE
      ! both of same type, but none of the data1d arrays allocated
      PRINT *,"ERROR in marray_base::AssignMArray_0: ma%data1d and this%data1d not associated"
      STOP 1
    END IF
  END SUBROUTINE AssignMArray_0
  
  !> assign 1D fortran array to mesh array
  PURE SUBROUTINE AssignMArray_1(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(INUM*JNUM*KNUM*this%DIMS(1)*this%DIMS(2)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
    this%data1d(:) = a(:)
  END SUBROUTINE AssignMArray_1

  !> assign 2D fortran array to mesh array
  PURE SUBROUTINE AssignMArray_2(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data2d,1),SIZE(this%data2d,2)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
    this%data2d(:,:) = a(:,:)
  END SUBROUTINE AssignMArray_2

  !> assign 3D fortran array to mesh array
  PURE SUBROUTINE AssignMArray_3(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data3d,1),SIZE(this%data3d,2),SIZE(this%data3d,3)), &
                          INTENT(IN) :: a
    !------------------------------------------------------------------------!
    this%data3d(:,:,:) = a(:,:,:)
  END SUBROUTINE AssignMArray_3

  !> assign 4D fortran array to mesh array
  PURE SUBROUTINE AssignMArray_4(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data4d,1),SIZE(this%data4d,2),SIZE(this%data4d,3), &
                SIZE(this%data4d,4)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
    this%data4d(:,:,:,:) = a(:,:,:,:)
  END SUBROUTINE AssignMArray_4

  !> assign 5D fortran array to mesh array
  PURE SUBROUTINE AssignMArray_5(this,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    REAL, DIMENSION(SIZE(this%data5d,1),SIZE(this%data5d,2),SIZE(this%data5d,3), &
                SIZE(this%data5d,4),SIZE(this%data5d,5)), INTENT(IN) :: a
    !------------------------------------------------------------------------!
    this%data5d(:,:,:,:,:) = a(:,:,:,:,:)
  END SUBROUTINE AssignMArray_5

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
  PURE FUNCTION MultMArray_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data1d)) :: c
    !------------------------------------------------------------------------!
    IF (SIZE(a%data1d).EQ.SIZE(b%data1d)) &
        c(:) = a%data1d(:) * b%data1d(:)
  END FUNCTION MultMArray_0

  !> multiply 1D fortran array and mesh arrays
  PURE FUNCTION MultMArray_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data1d)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data1d)) :: c
    !------------------------------------------------------------------------!
    c(:) = a%data1d(:) * b(:)
  END FUNCTION MultMArray_1

  !> multiply 2D fortran array and mesh arrays
  PURE FUNCTION MultMArray_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data2d,1),SIZE(a%data2d,2)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data2d,1),SIZE(a%data2d,2)) :: c
    !------------------------------------------------------------------------!
    c(:,:) = a%data2d(:,:) * b(:,:)
  END FUNCTION MultMArray_2

  !> multiply 3D fortran array and mesh arrays
  PURE FUNCTION MultMArray_3(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data3d,1),SIZE(a%data3d,2),SIZE(a%data3d,3)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data3d,1),SIZE(a%data3d,2),SIZE(a%data3d,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data3d(:,:,:) * b(:,:,:)
  END FUNCTION MultMArray_3

  !> multiply 4D fortran array and mesh arrays
  PURE FUNCTION MultMArray_4(a,b) RESULT(c)
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
  PURE FUNCTION MultMArray_5(a,b) RESULT(c)
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

  !> deconstructor of the mesh array
  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_base) :: this
    !-------------------------------------------------------------------!
    IF (ASSOCIATED(this%data1d)) DEALLOCATE(this%data1d)
    NULLIFY(this%data1d,this%data2d,this%data3d,this%data4d,this%data5d)
    this%RANK = 0
    this%DIMS = 0
  END SUBROUTINE Destroy

  FUNCTION CreateSelection(idx) RESULT(new_sel)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(selection_base) :: new_sel
    INTEGER, OPTIONAL    :: idx(6)
    !-------------------------------------------------------------------!
    INTEGER :: err
    !-------------------------------------------------------------------!
    IF (IDX_INIT) THEN
      ! allocate 1D mask array
      ALLOCATE(new_sel%mask3d(IGMIN:IGMAX,JGMIN:JGMAX,KGMIN:KGMAX),STAT=err)
      IF (err.NE.0) THEN
        PRINT *,"ERROR in marray_base::CreateSelection: memory allocation failed"
        STOP 1
      END IF
      ! create 2D & 3D pointers into the 1D mask array
      new_sel%mask1d(1:INUM*JNUM*KNUM) => new_sel%mask3d
      new_sel%mask2d(1:INUM*JNUM,KGMIN:KGMAX) => new_sel%mask3d
      IF (PRESENT(idx)) THEN
        CALL new_sel%Cuboid(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
      ELSE
        CALL new_sel%Everything()
      END IF
    END IF
  END FUNCTION CreateSelection

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

  !> deconstructor of the mesh array
  SUBROUTINE Destroy_selection(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(selection_base) :: this
    !-------------------------------------------------------------------!
    IF (ASSOCIATED(this%mask3d)) DEALLOCATE(this%mask3d)
    NULLIFY(this%mask1d,this%mask2d)
  END SUBROUTINE Destroy_selection

END MODULE marray_base_mod
