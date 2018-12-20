!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_common.f90                                                   #
!#                                                                           #
!# Copyright (C) 2006-2016                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \defgroup array
!! \{
!! \brief Family of array moduless
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief basic array module
!!
!----------------------------------------------------------------------------!
MODULE array
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE AssignMesharrayS_0, AssignMesharrayS_1, &
                      AssignMesharrayV_0, AssignMesharrayV_1, &
                      AssignMesharrayT_0, AssignMesharrayT_1
  END INTERFACE
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE AddMesharrayS_0, AddMesharrayS_1, AddMesharrayS_2, &
                      AddMesharrayV_0, AddMesharrayV_1, AddMesharrayV_2, &
                      AddMesharrayT_0, AddMesharrayT_1, AddMesharrayT_2
  END INTERFACE
  INTERFACE OPERATOR (*)
     MODULE PROCEDURE MultMesharrayS_0, MultMesharrayS_1, MultMesharrayS_2
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> type for scalar (rank 0) array data on the mesh
  TYPE MArrayS_TYP
     REAL, DIMENSION(:,:,:,:), POINTER :: data
     REAL, DIMENSION(:,:,:), POINTER :: &
                          center, &       !< geometric center
                          bcenter         !< bary center

     REAL, DIMENSION(:,:,:,:), POINTER :: &
                          faces, &        !< cell face centers
                          corners         !< cell corners
  END TYPE
  !> type for vector (rank 1) array data on the mesh
  TYPE MArrayV_TYP
     REAL, DIMENSION(:,:,:,:,:), POINTER :: data
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                          center, &       !< geometric center
                          bcenter         !< bary center

     REAL, DIMENSION(:,:,:,:,:), POINTER :: &
                          faces, &        !< cell face centers
                          corners         !< cell corners
  END TYPE
  !> type for tensor (rank 2) array data on the mesh
  TYPE MArrayT_TYP
     REAL, DIMENSION(:,:,:,:,:,:), POINTER :: data
     REAL, DIMENSION(:,:,:,:,:), POINTER :: &
                          center, &       !< geometric center
                          bcenter         !< bary center

     REAL, DIMENSION(:,:,:,:,:,:), POINTER :: &
                          faces, &        !< cell face centers
                          corners         !< cell corners
  END TYPE
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       MArrayS_TYP, MArrayV_TYP, MArrayT_TYP, &
       OPERATOR(+), &
       OPERATOR(*), &
       ASSIGNMENT(=)
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public copy data field of scalar mesh array
  SUBROUTINE AssignMesharrayS_0(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN) :: in
    TYPE(MArrayS_TYP),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:) = in%data(:,:,:,:)
  END SUBROUTINE AssignMesharrayS_0

  !> \public copy rank 3 array to data field of scalar mesh array
  SUBROUTINE AssignMesharrayS_1(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: in
    TYPE(MArrayS_TYP),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:) = in(:,:,:,:)
  END SUBROUTINE AssignMesharrayS_1

  !> \public copy data field of vector mesh array
  SUBROUTINE AssignMesharrayV_0(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayV_TYP),INTENT(IN) :: in
    TYPE(MArrayV_TYP),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:,:) = in%data(:,:,:,:,:)
  END SUBROUTINE AssignMesharrayV_0

  !> \public copy rank 4 array to data field of vector mesh array
  SUBROUTINE AssignMesharrayV_1(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: in
    TYPE(MArrayV_TYP),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:,:) = in(:,:,:,:,:)
  END SUBROUTINE AssignMesharrayV_1

  !> \public copy data field of tensor mesh array
  SUBROUTINE AssignMesharrayT_0(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayT_TYP),INTENT(IN) :: in
    TYPE(MArrayT_TYP),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:,:,:) = in%data(:,:,:,:,:,:)
  END SUBROUTINE AssignMesharrayT_0

  !> \public copy rank 5 array to data field of tensor mesh array
  SUBROUTINE AssignMesharrayT_1(out,in)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: in
    TYPE(MArrayT_TYP),INTENT(OUT) :: out
    !------------------------------------------------------------------------!
    out%data(:,:,:,:,:,:) = in(:,:,:,:,:,:)
  END SUBROUTINE AssignMesharrayT_1

  !> \public add two scalar mesh arrays
  PURE FUNCTION AddMesharrayS_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data(:,:,:,:)+b%data(:,:,:,:)
  END FUNCTION AddMesharrayS_0

  !> \public add scalar mesh arrays and rank 3 data array
  PURE FUNCTION AddMesharrayS_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data(:,:,:,:)+b(:,:,:,:)
  END FUNCTION AddMesharrayS_1

  !> \public add rank 3 data array and scalar mesh arrays
  PURE FUNCTION AddMesharrayS_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a(:,:,:,:)+b%data(:,:,:,:)
  END FUNCTION AddMesharrayS_2

  !> \public add two vector mesh arrays
  PURE FUNCTION AddMesharrayV_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayV_TYP),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a%data(:,:,:,:,:)+b%data(:,:,:,:,:)
  END FUNCTION AddMesharrayV_0

  !> \public add vector mesh array and rank 4 data array
  PURE FUNCTION AddMesharrayV_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayV_TYP),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a%data(:,:,:,:,:)+b(:,:,:,:,:)
  END FUNCTION AddMesharrayV_1

  !> \public add rank 4 data array and vector mesh array
  PURE FUNCTION AddMesharrayV_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayV_TYP),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4),SIZE(b%data,5)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4),SIZE(b%data,5)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:) = a(:,:,:,:,:)+b%data(:,:,:,:,:)
  END FUNCTION AddMesharrayV_2

  !> \public add two tensor mesh arrays
  PURE FUNCTION AddMesharrayT_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayT_TYP),INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5),SIZE(a%data,6)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:,:) = a%data(:,:,:,:,:,:)+b%data(:,:,:,:,:,:)
  END FUNCTION AddMesharrayT_0

  !> \public add tensor mesh array and rank 5 data array
  PURE FUNCTION AddMesharrayT_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayT_TYP),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5),SIZE(a%data,6)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4),SIZE(a%data,5),SIZE(a%data,6)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:,:) = a%data(:,:,:,:,:,:)+b(:,:,:,:,:,:)
  END FUNCTION AddMesharrayT_1

  !> \public add rank 5 data array and tensor mesh array
  PURE FUNCTION AddMesharrayT_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayT_TYP),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4),SIZE(b%data,5),SIZE(b%data,6)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4),SIZE(b%data,5),SIZE(b%data,6)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:,:,:) = a(:,:,:,:,:,:)+b%data(:,:,:,:,:,:)
  END FUNCTION AddMesharrayT_2

  !> \public multiply two scalar mesh arrays
  PURE FUNCTION MultMesharrayS_0(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN)  :: a,b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data(:,:,:,:)*b%data(:,:,:,:)
  END FUNCTION MultMesharrayS_0

  !> \public multiply scalar mesh array and rank 3 data array
  PURE FUNCTION MultMesharrayS_1(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN)  :: a
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)),INTENT(IN) :: b
    REAL, DIMENSION(SIZE(a%data,1),SIZE(a%data,2),SIZE(a%data,3),SIZE(a%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a%data(:,:,:,:)*b(:,:,:,:)
  END FUNCTION MultMesharrayS_1

  !> \public multiply rank 3 data array and scalar mesh array
  PURE FUNCTION MultMesharrayS_2(a,b) RESULT(c)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(MArrayS_TYP),INTENT(IN)  :: b
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4)),INTENT(IN) :: a
    REAL, DIMENSION(SIZE(b%data,1),SIZE(b%data,2),SIZE(b%data,3),SIZE(b%data,4)) :: c
    !------------------------------------------------------------------------!
    c(:,:,:,:) = a(:,:,:,:)*b%data(:,:,:,:)
  END FUNCTION MultMesharrayS_2

END MODULE array
