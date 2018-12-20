!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: field_base.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006-2017                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!> \defgroup field
!! \{
!! \brief Family of field modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief basic field module
!!
!----------------------------------------------------------------------------!
MODULE field_base_mod
  IMPLICIT NONE

  !--------------------------------------------------------------------------!
  !> abstract class for a set of scalar fields
  TYPE, ABSTRACT :: FieldSet
    CHARACTER(LEN=16) :: name
  END TYPE

  !> type for scalar (rank 0) field data on the mesh
  TYPE FieldS
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: data
    CHARACTER :: name
  CONTAINS
    PROCEDURE :: AssignFieldS_0
    PROCEDURE :: AssignFieldS_1
    PROCEDURE :: AddFieldS_0
    PROCEDURE :: AddFieldS_1
    PROCEDURE :: MultFieldS_0
    PROCEDURE :: MultFieldS_1
    GENERIC :: ASSIGNMENT(=) => AssignFieldS_0, AssignFieldS_1
    GENERIC :: OPERATOR(+) => AddFieldS_0, AddFieldS_1
    GENERIC :: OPERATOR(*) => MultFieldS_0, MultFieldS_1
  END TYPE

  !> type for vector (rank 1) array data on the mesh
  TYPE FieldV
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: data
    CHARACTER :: name
  CONTAINS
    PROCEDURE :: AssignFieldV_0
    PROCEDURE :: AssignFieldV_1
    PROCEDURE :: AddFieldV_0
    PROCEDURE :: AddFieldV_1
    GENERIC :: ASSIGNMENT(=) => AssignFieldV_0, AssignFieldV_1
    GENERIC :: OPERATOR(+) => AddFieldV_0, AddFieldV_1
  END TYPE

  !> type for tensor (rank 2) array data on the mesh
  TYPE FieldT
    REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: data
    CHARACTER :: name
  CONTAINS
    PROCEDURE :: AssignFieldT_0
    PROCEDURE :: AssignFieldT_1
    PROCEDURE :: AddFieldT_0
    PROCEDURE :: AddFieldT_1
    GENERIC :: ASSIGNMENT(=) => AssignFieldT_0, AssignFieldT_1
    GENERIC :: OPERATOR(+) => AddFieldT_0, AddFieldT_1
  END TYPE
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public copy data field of scalar mesh array
  SUBROUTINE AssignFieldS_0(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN) :: in
    CLASS(FieldS),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:) = in%data(:,:,:)
  END SUBROUTINE AssignFieldS_0

  !> \public copy rank 3 array to data field of scalar mesh array
  SUBROUTINE AssignFieldS_1(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:), INTENT(IN) :: in
    CLASS(FieldS),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:) = in(:,:,:)
  END SUBROUTINE AssignFieldS_1

  !> \public copy data field of vector mesh array
  SUBROUTINE AssignFieldV_0(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldV),INTENT(IN) :: in
    CLASS(FieldV),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:) = in%data(:,:,:,:)
  END SUBROUTINE AssignFieldV_0

  !> \public copy rank 4 array to data field of vector mesh array
  SUBROUTINE AssignFieldV_1(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: in
    CLASS(FieldV),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:) = in(:,:,:,:)
  END SUBROUTINE AssignFieldV_1

  !> \public copy data field of tensor mesh array
  SUBROUTINE AssignFieldT_0(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldT),INTENT(IN) :: in
    CLASS(FieldT),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:,:) = in%data(:,:,:,:,:)
  END SUBROUTINE AssignFieldT_0

  !> \public copy rank 5 array to data field of tensor mesh array
  SUBROUTINE AssignFieldT_1(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: in
    CLASS(FieldT),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:,:) = in(:,:,:,:,:)
  END SUBROUTINE AssignFieldT_1

  !> \public add two scalar mesh arrays
  PURE FUNCTION AddFieldS_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data(:,:,:)+b%data(:,:,:)
  END FUNCTION AddFieldS_0

  !> \public add scalar mesh arrays and rank 3 data array
  PURE FUNCTION AddFieldS_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data(:,:,:)+b(:,:,:)
  END FUNCTION AddFieldS_1

  !> \public add rank 3 data array and scalar mesh arrays
  PURE FUNCTION AddFieldS_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a(:,:,:)+b%data(:,:,:)
  END FUNCTION AddFieldS_2

  !> \public add two vector mesh arrays
  PURE FUNCTION AddFieldV_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldV),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data(:,:,:,:)+b%data(:,:,:,:)
  END FUNCTION AddFieldV_0

  !> \public add vector mesh array and rank 4 data array
  PURE FUNCTION AddFieldV_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldV),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data(:,:,:,:)+b(:,:,:,:)
  END FUNCTION AddFieldV_1

  !> \public add rank 4 data array and vector mesh array
  PURE FUNCTION AddFieldV_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldV),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a(:,:,:,:)+b%data(:,:,:,:)
  END FUNCTION AddFieldV_2

  !> \public add two tensor mesh arrays
  PURE FUNCTION AddFieldT_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldT),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a%data(:,:,:,:,:)+b%data(:,:,:,:,:)
  END FUNCTION AddFieldT_0

  !> \public add tensor mesh array and rank 5 data array
  PURE FUNCTION AddFieldT_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldT),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a%data(:,:,:,:,:)+b(:,:,:,:,:)
  END FUNCTION AddFieldT_1

  !> \public add rank 5 data array and tensor mesh array
  PURE FUNCTION AddFieldT_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldT),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4),SIZE(b%data,5)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4),SIZE(b%data,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a(:,:,:,:,:)+b%data(:,:,:,:,:)
  END FUNCTION AddFieldT_2

  !> \public multiply two scalar mesh arrays
  PURE FUNCTION MultFieldS_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN)  :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data(:,:,:)*b%data(:,:,:)
  END FUNCTION MultFieldS_0

  !> \public multiply scalar mesh array and rank 3 data array
  PURE FUNCTION MultFieldS_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN)  :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a%data(:,:,:)*b(:,:,:)
  END FUNCTION MultFieldS_1

  !> \public multiply rank 3 data array and scalar mesh array
  PURE FUNCTION MultFieldS_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(FieldS),INTENT(IN)  :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:) = a(:,:,:)*b%data(:,:,:)
  END FUNCTION MultFieldS_2

END MODULE field_base_mod
