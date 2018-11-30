!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_periodic.f90                                             #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jubin Lirawi     <jlirawi@astrophysik.uni-kiel.de>                        #
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
!> \author Jubin Lirawi
!!
!! \brief Boundary module for periodic boundary conditions
!!
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
    PROCEDURE :: Finalize
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
    CLASS(physics_base),      INTENT(IN)    :: Physics
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    TYPE(Dict_TYP),POINTER                  :: config
    INTEGER                                 :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: dir
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,PERIODIC,boundcond_name,dir,config)
  END SUBROUTINE InitBoundary_periodic

  !> \public Applies the periodic boundary condition
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_periodic), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN) :: Mesh
    CLASS(physics_base),      INTENT(IN) :: Physics
    REAL,                     INTENT(IN) :: time
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER                              :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(INOUT)                        :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(this%Direction%GetType())
    CASE(WEST)
!NEC$ IVDEP
       DO i=1,Mesh%GINUM
         pvar(Mesh%IMIN-i,:,:,:) = pvar(Mesh%IMAX-i+1,:,:,:)
       END DO
    CASE(EAST)
!NEC$ IVDEP
       DO i=1,Mesh%GINUM
         pvar(Mesh%IMAX+i,:,:,:) = pvar(Mesh%IMIN+i-1,:,:,:)
       END DO
    CASE(SOUTH)
!NEC$ IVDEP
       DO j=1,Mesh%GJNUM
         pvar(:,Mesh%JMIN-j,:,:) = pvar(:,Mesh%JMAX-j+1,:,:)
       END DO
    CASE(NORTH)
!NEC$ IVDEP
       DO j=1,Mesh%GJNUM
         pvar(:,Mesh%JMAX+j,:,:) = pvar(:,Mesh%JMIN+j-1,:,:)
       END DO
    CASE(BOTTOM)
!NEC$ IVDEP
      DO k=1,Mesh%GKNUM
        pvar(:,:,Mesh%KMIN-k,:) = pvar(:,:,Mesh%KMAX-k+1,:)
      END DO
    CASE(TOP)
!NEC$ IVDEP
      DO k=1,Mesh%GKNUM
        pvar(:,:,Mesh%KMAX+k,:) = pvar(:,:,Mesh%KMIN+k-1,:)
      END DO
    END SELECT
  END SUBROUTINE SetBoundaryData


  !> \public Destructor for periodic boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_periodic), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize


END MODULE boundary_periodic_mod
