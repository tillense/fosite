!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_inner.f90                                                #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> \author Jannes Klee
!!
!! \brief Boundary module for the inner boundaries (only necessary in parallel
!!        runs)
!!
!! \attention It is not possible to set the boundaries with these boundary
!!            conditions. It is a placeholder class, in order to state that
!!            there is an inner boundary, which needs to be handled by MPI.
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_inner_mod
  USE boundary_base_mod
  USE marray_compound_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_inner
  CONTAINS
    PROCEDURE :: InitBoundary_inner
    PROCEDURE :: Finalize
    PROCEDURE :: SetBoundaryData
  END TYPE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "none"
  !--------------------------------------------------------------------------!
  PUBLIC      :: boundary_inner
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for inner boundary conditions
  SUBROUTINE InitBoundary_inner(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_inner), INTENT(INOUT) :: this
    CLASS(physics_base),   INTENT(IN)    :: Physics
    CLASS(mesh_base),      INTENT(IN)    :: Mesh
    TYPE(Dict_TYP), POINTER              :: config
    INTEGER                              :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)                           :: dir
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL this%InitBoundary(Mesh,Physics,NONE,boundcond_name,dir,config)
#endif
  END SUBROUTINE InitBoundary_inner

  !> \public Applies the inner boundary condition
  !!
  !! \attention This routine does nothing.
  SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_inner), INTENT(INOUT) :: this
    CLASS(mesh_base),      INTENT(IN) :: Mesh
    CLASS(physics_base),   INTENT(IN) :: Physics
    REAL,                  INTENT(IN) :: time
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    ! nothing to do
  END SUBROUTINE SetBoundaryData

  !> \public Destructor for inner boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_inner), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize


END MODULE boundary_inner_mod
