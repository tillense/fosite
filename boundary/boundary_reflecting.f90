!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_reflecting.f90                                           #
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
!! \brief Boundary module for reflecting boundaries
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_reflecting_mod
  USE mesh_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_reflecting
    !> mask arrays to specify reflection, i.e. change of sign
    LOGICAL, DIMENSION(:), ALLOCATABLE     :: reflX,reflY,reflZ
  CONTAINS
    PROCEDURE :: InitBoundary_reflecting
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE boundary_reflecting
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "reflecting"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
      boundary_reflecting
  !--------------------------------------------------------------------------!

CONTAINS
  !TODO: NOT VERIFIED
  !> \public Constructor for reflecting boundary conditions
  SUBROUTINE InitBoundary_reflecting(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_reflecting), INTENT(INOUT) :: this
    CLASS(physics_base),        INTENT(IN)    :: Physics
    CLASS(mesh_base),           INTENT(IN)    :: Mesh
    TYPE(Dict_TYP), POINTER,    INTENT(IN)    :: config
    INTEGER,                    INTENT(IN)    :: dir
    !------------------------------------------------------------------------!
    INTEGER                                   :: err
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,REFLECTING,boundcond_name,dir,config)

    ALLOCATE(                      &
         this%reflX(Physics%VNUM), &
         this%reflY(Physics%VNUM), &
         this%reflZ(Physics%VNUM), &
         STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_reflecting", "Unable to allocate memory.")
    END IF
    ! this tells us which vars get the opposite sign/vanish at cell faces;
    ! e.g. vertical velocities (depends on the underlying physics)
    CALL Physics%ReflectionMasks(Mesh,this%reflX,this%reflY,this%reflZ)
  END SUBROUTINE InitBoundary_reflecting

  !> \public Applies the reflecting boundary condition
  !!
  !! \warning Not rigorously tested!
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_reflecting), INTENT(INOUT)    :: this
    CLASS(mesh_base),           INTENT(IN)    :: Mesh
    CLASS(physics_base),        INTENT(IN)    :: Physics
    REAL,                       INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                                INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER                                   :: i,j,k,m
    !------------------------------------------------------------------------!
    SELECT CASE(this%direction%GetType())
    CASE(WEST)
!NEC$ SHORTLOOP
      DO m=1,Physics%VNUM
        IF (this%reflX(m)) THEN
!NEC$ SHORTLOOP
          DO i=1,Mesh%GINUM
            pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m) &
              = -pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        ELSE
!NEC$ SHORTLOOP
          DO i=1,Mesh%GINUM
            pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m) &
              = pvar(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        END IF
      END DO
    CASE(EAST)
!NEC$ SHORTLOOP
      DO m=1,Physics%VNUM
        IF (this%reflX(m)) THEN
!NEC$ SHORTLOOP
          DO i=1,Mesh%GINUM
            pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m) &
              = -pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        ELSE
!NEC$ SHORTLOOP
          DO i=1,Mesh%GINUM
            pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m) &
              = pvar(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        END IF
      END DO
    CASE(SOUTH)
!NEC$ SHORTLOOP
      DO m=1,Physics%VNUM
        IF (this%reflY(m)) THEN
!NEC$ SHORTLOOP
          DO j=1,Mesh%GJNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,m) &
              = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        ELSE
!NEC$ SHORTLOOP
          DO j=1,Mesh%GJNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,m) &
              = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        END IF
      END DO
    CASE(NORTH)
!NEC$ SHORTLOOP
      DO m=1,Physics%VNUM
        IF (this%reflY(m)) THEN
!NEC$ SHORTLOOP
          DO j=1,Mesh%GJNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,m) &
              = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        ELSE
!NEC$ SHORTLOOP
          DO j=1,Mesh%GJNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,m) &
              = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX,m)
          END DO
        END IF
      END DO
    CASE(BOTTOM)
!NEC$ SHORTLOOP
      DO m=1,Physics%VNUM
        IF (this%reflZ(m)) THEN
!NEC$ SHORTLOOP
          DO j=1,Mesh%GKNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,m) &
              = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN+k-1,m)
          END DO
        ELSE
!NEC$ SHORTLOOP
          DO j=1,Mesh%GKNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,m) &
              = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN+k-1,m)
          END DO
        END IF
      END DO
    CASE(TOP)
!NEC$ SHORTLOOP
      DO m=1,Physics%VNUM
        IF (this%reflZ(m)) THEN
!NEC$ SHORTLOOP
          DO j=1,Mesh%GKNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,m) &
              = -pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX-k+1,m)
          END DO
        ELSE
!NEC$ SHORTLOOP
          DO j=1,Mesh%GKNUM
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,m) &
              = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX-k+1,m)
          END DO
        END IF
      END DO
    END SELECT
  END SUBROUTINE SetBoundaryData

  !> \public Destructor for reflecting boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_reflecting), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%reflX,this%reflY,this%reflZ)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

  END MODULE boundary_reflecting_mod
