!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_cellvector.f90                                             #
!#                                                                           #
!# Copyright (C) 2018,2021                                                   #
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
!! \brief derived mesh array class for vector cell data
!!
!! \extends marray_base
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_cellvector_mod
  USE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> data types and methods
  TYPE, EXTENDS(marray_base) :: marray_cellvector
    REAL, DIMENSION(:,:,:,:), POINTER &
                             :: center => null(), &      !< geometric center
                                bcenter => null()        !< bary center

    REAL, DIMENSION(:,:,:,:,:), POINTER &
                             :: faces => null(), &       !< cell face centers
                                corners => null()        !< cell corners
    CONTAINS
    PROCEDURE :: AssignPointers
    PROCEDURE :: Destroy
    FINAL     :: Finalize
  END TYPE
  INTERFACE marray_cellvector
    MODULE PROCEDURE CreateMArray_cellvector
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_cellvector

CONTAINS

  FUNCTION CreateMArray_cellvector() RESULT(new_cv)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_cellvector) :: new_cv
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_cellvector::CreateMArray_cellvector: new cellvector"
#endif
    ! 1 center + 1 bcenter + 6 faces + 8 corners = 16
    IF (new_cv%Init(16,3)) return ! immediately return if successful
#ifdef DEBUG
    PRINT *,"ERROR in marray_cellvector::CreateMArray: cellvector initialization failed"
    STOP 1
#endif
  END FUNCTION CreateMArray_cellvector
  
  FUNCTION AssignPointers(this) RESULT(success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_cellvector),INTENT(INOUT) :: this
    LOGICAL :: success
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_cellvector::AssignPointers: assigning pointers"
#endif
    success = this%marray_base%AssignPointers()
    ! assign array pointers
    IF (success) THEN
      this%center  => this%RemapBounds(this%data5d(:,:,:,1,:))
      this%bcenter => this%RemapBounds(this%data5d(:,:,:,2,:))
      this%faces   => this%RemapBounds(this%data5d(:,:,:,3:8,:))
      this%corners => this%RemapBounds(this%data5d(:,:,:,9:16,:))
#ifdef DEBUG
    ELSE
      PRINT *,"ERROR in marray_cellvector::AssignPointers: pointer assignment failed"
#endif
    END IF
  END FUNCTION AssignPointers

  !> polymorphic destructor of all mesh_cellvector classes
  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_cellvector) :: this
    !-------------------------------------------------------------------!
    CALL this%marray_base%Destroy() ! call inherited destructor
    NULLIFY(this%center,this%bcenter,this%faces,this%corners)
  END SUBROUTINE Destroy

  !> actual destructor of mesh_cellvector - this is called automatically if
  !! deallocate is invoked
  SUBROUTINE Destructor(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_cellvector) :: this
    !-------------------------------------------------------------------!
    CALL this%Destroy() ! call inherited marray destructor
  END SUBROUTINE Destructor

END MODULE marray_cellvector_mod
