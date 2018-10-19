!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_noslip.f03                                               #
!#                                                                           #
!# Copyright (C) 2010-2018                                                   #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
!> \author Tobias Illenseer, Bjoern Sperling
!!
!!
!! \brief Boundary module for noslip conditions
!!    (see [wikipedia](https://en.wikipedia.org/wiki/No-slip_condition) )
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_noslip_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE boundary_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base)  :: boundary_noslip
  CONTAINS
    PROCEDURE :: InitBoundary_noslip
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE boundary_noslip
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "noslip"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       Boundary_noslip, &
       WEST, EAST, SOUTH, NORTH, BOTTOM, TOP
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for noslip boundary conditions
  SUBROUTINE InitBoundary_noslip(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_noslip),INTENT(INOUT) :: this
    CLASS(Mesh_base)      ,INTENT(IN) :: Mesh
    CLASS(Physics_base)   ,INTENT(IN) :: Physics
    TYPE(Dict_TYP),POINTER,INTENT(IN) :: config
    INTEGER               ,INTENT(IN) :: dir
    !------------------------------------------------------------------------!
    INTEGER       :: err = 0
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,NOSLIP,boundcond_name,dir,config)

    ! allocate memory for boundary data
    SELECT CASE(dir)
    CASE(WEST,EAST)
       ALLOCATE(this%data(Mesh%GINUM,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            STAT=err)
       this%data(:,:,:,:)=0.0
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            STAT=err)
       this%data(:,:,:,:)=0.0
    CASE(BOTTOM,TOP)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%GKNUM, Physics%VNUM), &
             STAT=err)
       this%data(:,:,:,:)=0.0
    END SELECT

    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_noslip", "Unable to allocate memory.")
    END IF
  END SUBROUTINE InitBoundary_noslip

  !> \public Applies the noslip boundary condition
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_noslip),INTENT(INOUT) :: this
    CLASS(Mesh_base)      ,INTENT(IN) :: Mesh
    CLASS(Physics_base)   ,INTENT(IN) :: Physics
    REAL                  ,INTENT(IN) :: time
    REAL,INTENT(INOUT) :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum)
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
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY) &
               = pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY) &
               = -pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY)
          IF (Physics%ZVELOCITY.GT.0) THEN
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY)
          END IF
          ! vanishing pressure gradient at the boundary
          IF (Physics%PRESSURE.GT.0) THEN
             pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE) &
                  = pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE)
          END IF
       END DO
    CASE(EAST)
!NEC$ UNROLL(4)
       DO i=1,Mesh%GNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY) &
               = pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY) &
               = -pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY)
          IF (Physics%ZVELOCITY.GT.0) THEN
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY) &
               = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY)
          END IF
          ! vanishing pressure gradient at the boundary
          IF (Physics%PRESSURE.GT.0) THEN
             pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE) &
                  = pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE)
          END IF
       END DO
    CASE(SOUTH)
!NEC$ UNROLL(4)
       DO j=1,Mesh%GJNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY) &
               = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY) &
               = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY)
          IF (Physics%ZVELOCITY.GT.0) THEN
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY)
          END IF
          ! vanishing pressure gradient at the boundary
          IF (Physics%PRESSURE.GT.0) THEN
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE) &
                  = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE)
          END IF
       END DO
    CASE(NORTH)
!NEC$ UNROLL(4)
       DO j=1,Mesh%GNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY) &
               = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY)
          ! normal velocity
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY) &
               = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY)
          ! tangential velocities
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,Physics%XVELOCITY)
          IF (Physics%ZVELOCITY.GT.0) THEN
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,Physics%ZVELOCITY)
          END IF
          ! vanishing pressure gradient at the boundary
          IF (Physics%PRESSURE.GT.0) THEN
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE) &
                  = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE)
          END IF
       END DO
    CASE(BOTTOM)
!NEC$ UNROLL(4)
       DO k=1,Mesh%GKNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,Physics%DENSITY) &
               = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%kMIN+k-1,Physics%DENSITY)
          ! normal velocity
          IF (Physics%ZVELOCITY.GT.0) THEN
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,Physics%ZVELOCITY) &
               = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN+k-1,Physics%ZVELOCITY)
          END IF
          ! tangential velocities
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,Physics%XVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,Physics%XVELOCITY)
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,Physics%YVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,Physics%YVELOCITY)
          ! vanishing pressure gradient at the boundary
          IF (Physics%PRESSURE.GT.0) THEN
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,Physics%PRESSURE) &
                  = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN+k-1,Physics%PRESSURE)
          END IF
       END DO
    CASE(TOP)
!NEC$ UNROLL(4)
       DO k=1,Mesh%GKNUM
          ! vanishing density gradient at the boundary
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,Physics%DENSITY) &
               = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX-k+1,Physics%DENSITY)
          ! normal velocity
          IF (Physics%ZVELOCITY.GT.0) THEN
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,Physics%ZVELOCITY) &
               = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX-k+1,Physics%ZVELOCITY)
          END IF
          ! tangential velocities
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,Physics%XVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,Physics%XVELOCITY)
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,Physics%YVELOCITY) &
               = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,Physics%YVELOCITY)
          ! vanishing pressure gradient at the boundary
          IF (Physics%PRESSURE.GT.0) THEN
             pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,Physics%PRESSURE) &
                  = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX-k+1,Physics%PRESSURE)
          END IF
       END DO

    END SELECT
  END SUBROUTINE SetBoundaryData

  !> \public Destructor for noslip boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_noslip), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%data)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE boundary_noslip_mod
