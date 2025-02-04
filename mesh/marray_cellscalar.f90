!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: marray_cellscalar.f90                                             #
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
!! \brief derived mesh array class for scalar cell data
!!
!! \extends marray_base
!! \ingroup marray
!----------------------------------------------------------------------------!
MODULE marray_cellscalar_mod
  USE marray_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> data types and methods
  TYPE, EXTENDS(marray_base) :: marray_cellscalar
    REAL, DIMENSION(:,:,:), POINTER :: &
                                    center => null(), &  !< geometric center
                                    bcenter=> null()     !< bary center

    REAL, DIMENSION(:,:,:,:), POINTER :: &
                                    faces => null(), &   !< cell face centers
                                    corners => null()    !< cell corners
    CONTAINS
    PROCEDURE :: AssignPointers
    PROCEDURE :: Destroy
    FINAL     :: Finalize
  END TYPE
  INTERFACE marray_cellscalar
    MODULE PROCEDURE CreateMArray_cellscalar
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: marray_cellscalar

CONTAINS

  FUNCTION CreateMArray_cellscalar() RESULT(new_cs)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_cellscalar) :: new_cs
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_cellscalar::CreateMArray_cellscalar: creating new cellscalar"
#endif
    ! 1 center + 1 bcenter + 6 faces + 8 corners = 16
    IF (new_cs%Init(16)) return ! immediately return if successful
#ifdef DEBUG
    PRINT *,"ERROR in marray_cellscalar::CreateMArray: cellscalar initialization failed"
    STOP 1
#endif
  END FUNCTION CreateMArray_cellscalar

  FUNCTION AssignPointers(this) RESULT(success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(marray_cellscalar),INTENT(INOUT) :: this
    LOGICAL :: success
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_cellscalar::AssignPointers: assigning pointers"
#endif
    success = this%marray_base%AssignPointers()
    IF (success) THEN
      ! assign array pointers
      this%center  => this%RemapBounds(this%data4d(:,:,:,1))
      this%bcenter => this%RemapBounds(this%data4d(:,:,:,2))
      this%faces   => this%RemapBounds(this%data4d(:,:,:,3:8))
      this%corners => this%RemapBounds(this%data4d(:,:,:,9:16))
#ifdef DEBUG
    ELSE
      PRINT *,"ERROR in marray_cellscalar::AssignPointers: pointer assignment failed"
#endif
    END IF
  END FUNCTION AssignPointers

  !> destructor of cellscalar - manual deallocation / reset of pointers
  SUBROUTINE Destroy(this,called_from_finalize)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(marray_cellscalar) :: this
    LOGICAL, OPTIONAL :: called_from_finalize
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_cellscalar::Destroy: nullify pointers"
#endif
    NULLIFY(this%center,this%bcenter,this%faces,this%corners)
    ! only call inherited destructor if not called from Finalize
    IF (PRESENT(called_from_finalize)) THEN
       IF (called_from_finalize) RETURN
    END IF
    CALL this%marray_base%Destroy()
  END SUBROUTINE Destroy

  !> destructor of cellscalar - this is called automatically on deallocation
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(marray_cellscalar) :: this
    !-------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in marray_cellscalar::Finalize called"
#endif
    CALL this%Destroy(.TRUE.)
  END SUBROUTINE Finalize

END MODULE marray_cellscalar_mod
