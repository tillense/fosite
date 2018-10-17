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
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief base class for mesh arrays
!!
!! \ingroup mesh
!----------------------------------------------------------------------------!
MODULE marray_base_mod
  USE logging_base_mod
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
  !> basic mesh array class
  TYPE :: marray_base
!     INTEGER, PRIVATE, POINTER :: IGMIN => null(),IGMAX => null()
!     INTEGER, PRIVATE, POINTER :: JGMIN => null(),JGMAX => null()
!     INTEGER, PRIVATE, POINTER :: KGMIN => null(),KGMAX => null()
    INTEGER       :: RANK = 0,DIMS(2) = 1
    REAL, POINTER :: data1d(:) => null()
    REAL, POINTER :: data2d(:,:) => null()
    REAL, POINTER :: data3d(:,:,:) => null()
    REAL, POINTER :: data4d(:,:,:,:) => null()
    REAL, POINTER :: data5d(:,:,:,:,:) => null()
    CONTAINS
    PROCEDURE, PRIVATE :: AssignPointers
    PROCEDURE :: AssignMArray_0 !, AssignMArray_1, AssignMArray_2
    GENERIC   :: ASSIGNMENT (=) => AssignMArray_0 !, AssignMArray_1, AssignMArray_2
!     PROCEDURE :: AddMArray_0, AddMArray_1, AddMArray_2
!     GENERIC   :: OPERATOR (+) => AddMArray_0, AddMArray_1, AddMArray_2
!     PROCEDURE :: MultMArray_0, MultMArray_1, MultMArray_2
!     GENERIC   :: OPERATOR (*) => MultMArray_0, MultMArray_1, MultMArray_2
    PROCEDURE :: Destroy
  END TYPE
  INTERFACE marray_base
    MODULE PROCEDURE CreateMArray
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_base, &
    InitMeshProperties, &
    CloseMeshProperties

  CONTAINS

  FUNCTION CreateMArray(m,n) RESULT(new_ma)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_base) :: new_ma
    INTEGER, OPTIONAL, INTENT(IN) :: m,n
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
      ALLOCATE(new_ma%data1d(INUM*JNUM*KNUM*new_ma%DIMS(1)*new_ma%DIMS(2)))
      CALL new_ma%AssignPointers()
    END IF
  END FUNCTION CreateMArray
  
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
  
  SUBROUTINE CloseMeshProperties
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    IF (IDX_INIT) THEN
      IDX_INIT = .FALSE.
    END IF
  END SUBROUTINE CloseMeshProperties
  
  SUBROUTINE AssignPointers(this)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(this%data1d)) THEN
      SELECT CASE(this%RANK)
      CASE(0)
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data2d,[inum*jnum,knum])
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data3d,[inum,jnum,knum])
      CASE(1)
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data2d,[inum*jnum*knum,this%DIMS(1)])
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data3d,[inum*jnum,knum,this%DIMS(1)])
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data4d,[inum,jnum,knum,this%DIMS(1)])
      CASE(2)
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data2d,[inum*jnum*knum*this%DIMS(1),this%DIMS(2)])
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data3d,[inum*jnum*knum,this%DIMS(1),this%DIMS(2)])
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data4d,[inum,jnum,knum,this%DIMS(1)*this%DIMS(2)])
        CALL C_F_POINTER(C_LOC(this%data1d(LBOUND(this%data1d,1))), &
                        this%data5d,[inum,jnum,knum,this%DIMS(1),this%DIMS(2)])
      END SELECT
    END IF
  END SUBROUTINE AssignPointers

  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_base),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    this%RANK    = ma%RANK
    this%DIMS(:) = ma%DIMS(:)
    IF (.NOT.ASSOCIATED(this%data1d)) THEN
        ALLOCATE(this%data1d(SIZE(ma%data1d,1)),SOURCE=ma%data1d)
    ELSE
        this%data1d(:) = ma%data1d(:)
    END IF
    CALL this%AssignPointers()
  END SUBROUTINE AssignMArray_0
  
  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_base) :: this
    !-------------------------------------------------------------------!
    DEALLOCATE(this%data1d)
  END SUBROUTINE Destroy


END MODULE marray_base_mod
