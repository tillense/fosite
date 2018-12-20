!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_periodic.f90                                             #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!! \brief Boundary module for periodic boundary conditions
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_periodic_mod
  USE boundary_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_periodic
  CONTAINS
    PROCEDURE :: InitBoundary_periodic
    FINAL :: Finalize
    PROCEDURE :: SetBoundaryData
  END TYPE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "periodic"  
  !--------------------------------------------------------------------------!
  PUBLIC :: boundary_periodic
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for periodic boundary conditions
  SUBROUTINE InitBoundary_periodic(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_periodic), INTENT(INOUT) :: this
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(mesh_base), INTENT(IN) :: Mesh
    TYPE(Dict_TYP),POINTER &
                       :: config
    INTEGER            :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: dir
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,PERIODIC,boundcond_name,dir,config)
  END SUBROUTINE InitBoundary_periodic

  !> \public Applies the periodic boundary condition
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_periodic), INTENT(IN) :: this
    CLASS(mesh_base), INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER            :: i,j
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: pvar   
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(this%Direction%GetType())
    CASE(WEST)
!CDIR NODEP
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMAX-i+1,:,:)
       END DO
    CASE(EAST)
!CDIR NODEP
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMIN+i-1,:,:)
       END DO
    CASE(SOUTH)
!CDIR NODEP
       DO j=1,Mesh%GNUM
          pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMAX-j+1,:)
       END DO
    CASE(NORTH)
!CDIR NODEP
       DO j=1,Mesh%GNUM
          pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMIN+j-1,:)
       END DO
    END SELECT
  END SUBROUTINE SetBoundaryData


  !> \public Destructor for periodic boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(boundary_periodic), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizeBoundary()
  END SUBROUTINE Finalize


END MODULE boundary_periodic_mod
