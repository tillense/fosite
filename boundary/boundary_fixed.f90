!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_fixed.f90                                                #
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
!! \brief Boundary module for fixed in/outflow conditions
!!
!! Implementation of sub/supersonic in/outflow conditions with
!! user defined fixed data.
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_fixed_mod
  USE marray_compound_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE boundary_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base)  :: boundary_fixed
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: data   !< boundary data array
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: fixed  !< mask array for fixed bc
  CONTAINS
    PROCEDURE :: InitBoundary_fixed
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE boundary_fixed
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "fixed in/outflow"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
      boundary_fixed
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for fixed boundary conditions
  SUBROUTINE InitBoundary_fixed(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_fixed) :: this
    CLASS(Mesh_base)      :: Mesh
    CLASS(Physics_base)   :: Physics
    TYPE(Dict_TYP),POINTER:: config
    INTEGER               :: dir
    !------------------------------------------------------------------------!
    INTEGER            :: err = 0
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,config
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,FIXED,boundcond_name,dir,config)

    ! allocate memory for boundary data and mask
    SELECT CASE(this%direction%GetType())
    CASE(WEST,EAST)
       ALLOCATE(this%data(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            this%fixed(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            STAT=err)
     CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            this%fixed(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,Physics%VNUM), &
            STAT=err)
     CASE(Bottom,Top)
       ALLOCATE(this%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
            this%fixed(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_fixed", "Unable to allocate memory.")
    END IF
    ! fixed(:,:,:) defaults to EXTRAPOLATION everywhere
    this%fixed(:,:,:) = .FALSE.
    this%data(:,:,:,:) = 0.0
  END SUBROUTINE InitBoundary_fixed


  !> \public Applies the fixed boundary condition
  SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_Fixed),INTENT(INOUT) :: this
    CLASS(Mesh_base),INTENT(IN)      :: Mesh
    CLASS(Physics_base),INTENT(IN)   :: Physics
    REAL,INTENT(IN) :: time
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    SELECT CASE(this%direction%GetType())
    CASE(WEST)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!NEC$ UNROLL(4)
       DO i=1,Mesh%GINUM
          WHERE(this%fixed(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:))
             ! set fixed boundary data
             pvar%data4d(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
          ELSEWHERE
             ! first order extrapolation
             pvar%data4d(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) = &
             (i+1)*pvar%data4d(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
             - i*pvar%data4d(Mesh%IMIN+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
          END WHERE

       END DO
    CASE(EAST)
!NEC$ UNROLL(4)
       DO i=1,Mesh%GINUM
          WHERE(this%fixed(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:))
             ! set fixed boundary data
             pvar%data4d(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) = &
                  this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
          ELSEWHERE
             ! first order extrapolation
             pvar%data4d(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) = &
                  (i+1)*pvar%data4d(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
                  - i*pvar%data4d(Mesh%IMAX-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
          END WHERE
        END DO
    CASE(SOUTH)
!NEC$ UNROLL(4)
       DO j=1,Mesh%GJNUM
          WHERE(this%fixed(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,:))
             ! set fixed boundary data
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,:) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:)
          ELSEWHERE
             ! first order extrapolation
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,:) = &
               (j+1)*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Mesh%KMIN:Mesh%KMAX,:) &
               - j*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+1,Mesh%KMIN:Mesh%KMAX,:)
          END WHERE
        END DO
    CASE(NORTH)
!NEC$ UNROLL(4)
       DO j=1,Mesh%GJNUM
          WHERE(this%fixed(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,:))
             ! set fixed boundary data
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,:) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:)
          ELSEWHERE
             ! first order extrapolation
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,:) = &
                  (j+1)*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
                  - j*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-1,Mesh%KMIN:Mesh%KMAX,:)
          END WHERE

       END DO
    CASE(BOTTOM)
!NEC$ UNROLL(4)
          DO k=1,Mesh%GKNUM
            WHERE(this%fixed(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:))
             ! set fixed boundary data
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,:) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:)
          ELSEWHERE
             ! first order extrapolation
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,:) = &
                  (k+1)*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN,:) &
                  - k*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN+1,:)
          END WHERE
          END DO
    CASE(TOP)
!NEC$ UNROLL(4)
       DO k=1,Mesh%GKNUM
          WHERE(this%fixed(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:))
             ! set fixed boundary data
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,:) = &
                  this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:)
          ELSEWHERE
             ! first order extrapolation
             pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,:) = &
                  (k+1)*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX,:) &
                  - k*pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX-1,:)
          END WHERE
      END DO
    END SELECT
  END SUBROUTINE SetBoundaryData

  !> \public Destructor for fixed boundary conditions
  SUBROUTINE FINALIZE(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_fixed),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data,this%fixed)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE boundary_fixed_mod
