!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_axis.f90                                                 #
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
!! \brief Boundary module for axis boundaries
!!
!! \extends boundary_reflecting
!! \ingroup boundary
!> \todo compiles, but strange behavior at axis. Be carefull
!----------------------------------------------------------------------------!
MODULE boundary_axis_mod
  USE mesh_base_mod
  USE boundary_base_mod 
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_axis
    CONTAINS
      PROCEDURE :: InitBoundary_axis
      PROCEDURE :: SetBoundaryData
      FINAL     :: Finalize
  END TYPE
  CHARACTER(LEN=32), PARAMETER :: boundcond_name = "axis"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       boundary_axis, &
       WEST, EAST, SOUTH, NORTH, BOTTOM, TOP
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for the axis boundary condition
  SUBROUTINE InitBoundary_axis(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_axis),     INTENT(INOUT) :: this
    CLASS(physics_base),      INTENT(IN)    :: Physics
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    TYPE(Dict_TYP), POINTER,  INTENT(IN)    :: config
    INTEGER            :: dir
    !------------------------------------------------------------------------!
    INTEGER       :: err
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,AXIS,boundcond_name,dir,config)

    ALLOCATE(this%reflX(Physics%vnum), &
         this%reflY(Physics%vnum), &
         this%reflZ(Physics%vnum), &
         STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_axis", "Unable to allocate memory.")
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IMPORTANT:
    ! I don't know if this works for all geometries and all physics
    ! along all coordinate directions;
    ! check the underlying physics/geometry  modules
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL Physics%AxisMasks(this%reflX,this%reflY,this%reflZ)
  END SUBROUTINE InitBoundary_axis
 
  !TODO: NOT VERIFIED
  !> \public Applies the axis/reflecting boundary condition
 PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_axis),       INTENT(IN)    :: this
    CLASS(mesh_base),           INTENT(IN)    :: Mesh
    CLASS(physics_base),        INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                                INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    SELECT CASE(this%direction%GetType())
    CASE(WEST)
!CDIR UNROLL=8
       DO k=Mesh%KMIN,Mesh%KMAX
         DO j=Mesh%JMIN, Mesh%JMAX
!CDIR NODEP
          DO i=1,Mesh%GINUM
             WHERE (this%reflX)
                pvar(Mesh%IMIN-i,j,k,:) = -pvar(Mesh%IMIN+i-1,j,k,:)
             ELSEWHERE
                pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMIN+i-1,j,k,:)
             END WHERE
          END DO
        END DO
       END DO
    CASE(EAST)
!CDIR UNROLL=8
       DO k=Mesh%KMIN,Mesh%KMAX
         DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=1,Mesh%GINUM
             WHERE (this%reflX)
                pvar(Mesh%IMAX+i,j,k,:) = -pvar(Mesh%IMAX-i+1,j,k,:)
             ELSEWHERE
                pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMAX-i+1,j,k,:)
             END WHERE
          END DO
        END DO
       END DO
    CASE(SOUTH)
!CDIR UNROLL=4
      DO k=Mesh%KMIN,Mesh%KMAX
       DO j=1,Mesh%GJNUM
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%reflY)
                pvar(i,Mesh%JMIN-j,k,:) = -pvar(i,Mesh%JMIN+j-1,k,:)
             ELSEWHERE
                pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMIN+j-1,k,:)
             END WHERE
          END DO
       END DO
     END DO
    CASE(NORTH)
!CDIR UNROLL=4
      DO k=Mesh%KMIN,Mesh%KMAX
       DO j=1,Mesh%GJNUM
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%reflY)
                pvar(i,Mesh%JMAX+j,k,:) = -pvar(i,Mesh%JMAX-j+1,k,:)
             ELSEWHERE
                pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMAX-j+1,k,:)
             END WHERE
          END DO
       END DO
     END DO

    CASE(BOTTOM)
!CDIR UNROLL=4
      DO k=1,Mesh%GKNUM
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%reflZ)
                pvar(i,j,Mesh%KMIN-k,:) = -pvar(i,j,Mesh%KMIN+k-1,:)
             ELSEWHERE
                pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMIN+k-1,:)
             END WHERE
          END DO
       END DO
     END DO
    CASE(TOP)
!CDIR UNROLL=4
      DO k=1,Mesh%GKNUM
       DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
             WHERE (this%reflZ)
                pvar(i,j,Mesh%KMAX+k,:) = -pvar(i,j,Mesh%KMAX-k+1,:)
             ELSEWHERE
                pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMAX-k+1,:)
             END WHERE
          END DO
       END DO
     END DO

    END SELECT
  END SUBROUTINE SetBoundaryData

  !> \public Destructor for axis boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(boundary_axis), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%reflX,this%reflY,this%reflZ)
    CALL this%FinalizeBoundary()
  END SUBROUTINE Finalize



END MODULE boundary_axis_mod
