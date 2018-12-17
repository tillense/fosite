!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_custom.f03                                               #
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
!! \author Jannes Klee
!!
!! \brief Boundary module for custom conditions
!!
!! \attention Not well tested, use with care, since the transition to
!!            Fosite 3D. Use with care.
!----------------------------------------------------------------------------!
MODULE boundary_custom_mod
  USE boundary_base_mod
  USE boundary_fixed_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE cbspec_typ
    INTEGER, DIMENSION(:), ALLOCATABLE   :: bc
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: region
  END TYPE
  TYPE, EXTENDS(boundary_fixed) :: boundary_custom
    TYPE(cbspec_typ), DIMENSION(:), ALLOCATABLE :: cbspec 
  CONTAINS
    PROCEDURE :: InitBoundary_custom
    PROCEDURE :: Finalize
    PROCEDURE :: SetBoundaryData
  END TYPE
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "custom"
  !--------------------------------------------------------------------------!
  ! create bit masks for custom boundary conditions
  INTEGER, PARAMETER :: CUSTOM_NOGRAD   = INT(b'000000000000000') ! no gradients (default)
  INTEGER, PARAMETER :: CUSTOM_PERIOD   = INT(b'000000000000001') ! periodic
  INTEGER, PARAMETER :: CUSTOM_REFLECT  = INT(b'000000000000010') ! reflecting
  INTEGER, PARAMETER :: CUSTOM_REFLNEG  = INT(b'000000000000100') ! reflect and change sign
  INTEGER, PARAMETER :: CUSTOM_EXTRAPOL = INT(b'000000000001000') ! linear extrapolation
  INTEGER, PARAMETER :: CUSTOM_FIXED    = INT(b'000000000010000') ! set fixed boundary data
  INTEGER, PARAMETER :: CUSTOM_LOGEXPOL = INT(b'000000000100000') ! extrapolation of log values
  INTEGER, PARAMETER :: CUSTOM_OUTFLOW  = INT(b'000000001000000') ! nograd/reflect depending on flow direction
  INTEGER, PARAMETER :: CUSTOM_KEPLER   = INT(b'000000010000000') ! extrapolate according to Keplers law
  INTEGER, PARAMETER :: CUSTOM_ANGKEPLER= INT(b'000000100000000') ! same but for angular momentum instead of velocity
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       boundary_custom, &
       ! constants
       CUSTOM_NOGRAD, CUSTOM_PERIOD, CUSTOM_REFLECT, CUSTOM_REFLNEG, &
       CUSTOM_EXTRAPOL, CUSTOM_FIXED, CUSTOM_LOGEXPOL, &
       CUSTOM_OUTFLOW, CUSTOM_KEPLER, CUSTOM_ANGKEPLER
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for custom boundary conditions
  SUBROUTINE InitBoundary_custom(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_custom), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    TYPE(Dict_TYP),       POINTER       :: config
    INTEGER,              INTENT(IN)    :: dir
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,err = 0
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,CUSTOM,boundcond_name,dir,config)


    ! allocate memory for boundary data and mask
    SELECT CASE(this%GetDirection())
    CASE(WEST,EAST)
       ALLOCATE(this%data(Mesh%GINUM,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%cbtype(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%Rscale(Mesh%GINUM,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
            this%invRscale(Mesh%GINUM,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%cbtype(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%Rscale(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX), &
            this%invRscale(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX), &
            STAT=err)
    CASE(BOTTOM,TOP)
       ALLOCATE(this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%GKNUM,Physics%VNUM), &
            this%cbtype(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%Rscale(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%GKNUM), &
            this%invRscale(Mesh%IMIN:Mesh%IMAX,Mesh%KGMIN:Mesh%JMAX,Mesh%GKNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_custom", "Unable to allocate memory.")
    END IF

    SELECT CASE(this%GetDirection())
    CASE(WEST)
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GINUM
            this%Rscale(i,j,k) = Mesh%radius%bcenter(Mesh%IMIN-i,j,k) / Mesh%radius%bcenter(Mesh%IMIN,j,k)
          END DO
        END DO
      END DO
    CASE(EAST)
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GINUM
             this%Rscale(i,j,k) = Mesh%radius%bcenter(Mesh%IMAX+i,j,k) / Mesh%radius%bcenter(Mesh%IMAX,j,k)
          END DO
        END DO
      END DO
    CASE(SOUTH)
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=1,Mesh%GJNUM
          DO i= Mesh%IMIN,Mesh%IMAX
            this%Rscale(i,j,k) = Mesh%radius%bcenter(i,Mesh%JMIN-j,k) / Mesh%radius%bcenter(i,Mesh%JMIN,k)
          END DO
        END DO
      END DO
    CASE(NORTH)
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=1,Mesh%GJNUM
          DO i= Mesh%IMIN,Mesh%IMAX
            this%Rscale(i,j,k) = Mesh%radius%bcenter(i,Mesh%JMAX+j,k) / Mesh%radius%bcenter(i,Mesh%JMAX,k)
          END DO
        END DO
      END DO
    CASE(BOTTOM)
      DO k=1,Mesh%GKNUM
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i= Mesh%IMIN,Mesh%IMAX
            this%Rscale(i,j,k) = Mesh%radius%bcenter(i,j,Mesh%KMIN-k) / Mesh%radius%bcenter(i,j,Mesh%KMIN)
          END DO
        END DO
      END DO
    CASE(TOP)
      DO k=1,Mesh%GKNUM
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i= Mesh%IMIN,Mesh%IMAX
            this%Rscale(i,j,k) = Mesh%radius%bcenter(i,j,Mesh%KMAX+k) / Mesh%radius%bcenter(i,j,Mesh%KMAX)
          END DO
        END DO
      END DO
    END SELECT
    this%Rscale(:,:,:) = SQRT(this%Rscale(:,:,:))
    this%invRscale(:,:,:) = 1.0 / (this%Rscale(:,:,:) + TINY(1.0))
    ! this array contains the boundary condition for each primitive variable;
    ! the default setting is NO_GRADIENTS; the user has to assign reasonable
    ! values after initialization of the boundary module, e.g. set
    ! this%cbtype(:,:,1..Physics%VNUM) = {CUSTOM_NOGRAD | CUSTOM_PERIOD | ...}
    ! for each physical variable at each custom boundary
    this%cbtype(:,:,:) = CUSTOM_NOGRAD
    this%data(:,:,:,:) = 0.0
  END SUBROUTINE InitBoundary_custom


  !> \public Applies the custom boundary conditions
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_custom), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    REAL,                 INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,MESH%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                          INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    ! ATTENTION: If this%bctype(:) contains bogus values, this boundary module
    !            behaves like the NO_GRADIENTS boundary condition.
    SELECT CASE(this%GetDirection())
    CASE(WEST)
      DO i=1,Mesh%GNUM
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            WHERE (this%cbtype(j,k,:).EQ.CUSTOM_PERIOD)
               pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMAX-i+1,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_REFLECT)
               pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMIN+i-1,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_REFLNEG)
               pvar(Mesh%IMIN-i,j,k,:) = -pvar(Mesh%IMIN+i-1,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_EXTRAPOL)
               pvar(Mesh%IMIN-i,j,k,:) = (i+1)*pvar(Mesh%IMIN,j,k,:) - i*pvar(Mesh%IMIN+1,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_FIXED)
               pvar(Mesh%IMIN-i,j,k,:) = this%data(i,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_LOGEXPOL)
               pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMIN,j,k,:) &
                  * ABS(pvar(Mesh%IMIN,j,k,:) / pvar(Mesh%IMIN+1,j,k,:))**i
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_OUTFLOW &
                    .AND. pvar(Mesh%IMIN,j,k,:) .GE. 0.0 )
               !REFLNEG, else default (NO_GRADIENTS)
               pvar(Mesh%IMIN-i,j,k,:) = -pvar(Mesh%IMIN+i-1,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_KEPLER)
               pvar(Mesh%IMIN-i,j,k,:) = (pvar(Mesh%IMIN,j,k,:) + Mesh%radius%bcenter(Mesh%IMIN,j,k)*Mesh%OMEGA)&
                  * this%invRscale(i,j,k) - Mesh%radius%bcenter(Mesh%IMIN-i,j,k)*Mesh%OMEGA
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_ANGKEPLER)
               pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMIN,j,k,:) &
                  * this%Rscale(i,j,k)
!            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_POISSON.AND.this%accel(Mesh%IMIN,j,k,1).LT.TINY(1.e0))
!                pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMIN,j,k,:) &
!                  + (Mesh%bhy(Mesh%IMIN-i,j,k)-Mesh%bhy(Mesh%IMIN,j,k)) &
!                    * 0.5*(SQRT(Mesh%bhy(Mesh%IMIN+1,j,k)*(-this%accel(Mesh%IMIN+1,j,k,1))) &
!                       -SQRT(Mesh%bhy(Mesh%IMIN-1,j,k)*(-this%accel(Mesh%IMIN-1,j,k,1))))/Mesh%dlx(Mesh%IMIN,j,k)
            ELSEWHERE
               ! defaults to NO_GRADIENTS
               pvar(Mesh%IMIN-i,j,k,:) = pvar(Mesh%IMIN,j,k,:)
            END WHERE
          END DO
        END DO
      END DO
    CASE(EAST)
      DO i=1,Mesh%GNUM
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            WHERE (this%cbtype(j,k,:).EQ.CUSTOM_PERIOD)
               pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMIN+i-1,j,k,:)
            ELSEWHERE (this%cbtype(j,k,:).EQ.CUSTOM_REFLECT)
               pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMAX-i+1,j,k,:)
            ELSEWHERE (this%cbtype(j,k,:).EQ.CUSTOM_REFLNEG)
               pvar(Mesh%IMAX+i,j,k,:) = -pvar(Mesh%IMAX-i+1,j,k,:)
            ELSEWHERE (this%cbtype(j,k,:).EQ.CUSTOM_EXTRAPOL)
               pvar(Mesh%IMAX+i,j,k,:) = (i+1)*pvar(Mesh%IMAX,j,k,:) - i*pvar(Mesh%IMAX-1,j,k,:)
            ELSEWHERE (this%cbtype(j,k,:).EQ.CUSTOM_FIXED)
               pvar(Mesh%IMAX+i,j,k,:) = this%data(i,j,k,:)
            ELSEWHERE (this%cbtype(j,k,:).EQ.CUSTOM_LOGEXPOL)
               pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMAX,j,k,:) &
                  * ABS(pvar(Mesh%IMAX,j,k,:) / pvar(Mesh%IMAX-1,j,k,:))**i
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_OUTFLOW &
                    .AND. pvar(Mesh%IMAX,j,k,:) .LE. 0.0 )
               !REFLNEG, else default (NO_GRADIENTS)
               pvar(Mesh%IMAX+i,j,k,:) = -pvar(Mesh%IMAX-i+1,j,k,:)
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_KEPLER)
               pvar(Mesh%IMAX+i,j,k,:) = (pvar(Mesh%IMAX,j,k,:) + Mesh%radius%bcenter(Mesh%IMAX,j,k)*Mesh%OMEGA)&
                   * this%invRscale(i,j,k) - Mesh%radius%bcenter(Mesh%IMAX+i,j,k)*Mesh%OMEGA
            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_ANGKEPLER)
               pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMAX,j,k,:) &
                  * this%Rscale(i,j,k)
!            ELSEWHERE(this%cbtype(j,k,:).EQ.CUSTOM_POISSON.AND.this%accel(Mesh%IMAX,j,k,1).LT.TINY(1.e0))
!               pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMAX,j,k,:) &
!                  + (Mesh%bhy(Mesh%IMAX+i,j,k,k)-Mesh%bhy(Mesh%IMAX,j,k,k)) &
!                    * 0.5*(SQRT(Mesh%bhy(Mesh%IMAX+1,j,k,k)*(-this%accel(Mesh%IMAX+1,j,k,1))) &
!                       -SQRT(Mesh%bhy(Mesh%IMAX-1,j,k,k)*(-this%accel(Mesh%IMAX-1,j,k,1))))/Mesh%dlx(Mesh%IMAX,j,k,k)
            ELSEWHERE
               ! defaults to NO_GRADIENTS
               pvar(Mesh%IMAX+i,j,k,:) = pvar(Mesh%IMAX,j,k,:)
            END WHERE
          END DO
        END DO
      END DO
    CASE(SOUTH)
      DO j=1,Mesh%GNUM
        DO k=Mesh%KMIN,Mesh%KMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            WHERE (this%cbtype(i,k,:).EQ.CUSTOM_PERIOD)
               pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMAX-j+1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_REFLECT)
               pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMIN+j-1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_REFLNEG)
               pvar(i,Mesh%JMIN-j,k,:) = -pvar(i,Mesh%JMIN+j-1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_EXTRAPOL)
               pvar(i,Mesh%JMIN-j,k,:) = (j+1)*pvar(i,Mesh%JMIN,k,:) - j*pvar(i,Mesh%JMIN+1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_FIXED)
               pvar(i,Mesh%JMIN-j,k,:) = this%data(i,j,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_LOGEXPOL)
               pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMIN,k,:) &
                  * ABS(pvar(i,Mesh%JMIN,k,:) / pvar(i,Mesh%JMIN+1,k,:))**j
            ELSEWHERE(this%cbtype(i,k,:).EQ.CUSTOM_OUTFLOW &
                    .AND. pvar(i,Mesh%JMIN,k,:) .GE. 0.0 )
               !REFLNEG, else default (NO_GRADIENTS)
               pvar(i,Mesh%JMIN-j,k,:) = -pvar(i,Mesh%JMIN+j-1,k,:)
            ELSEWHERE(this%cbtype(i,k,:).EQ.CUSTOM_KEPLER)
               pvar(i,Mesh%JMIN-j,k,:) = (pvar(i,Mesh%JMIN,k,:) + Physics%bcradius(i,Mesh%JMIN,k)*Mesh%OMEGA)&
                  * this%invRscale(i,j,k) - Physics%bcradius(i,Mesh%JMIN-j,k)*Mesh%OMEGA
!                pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMIN,k,:) &
!                   * this%invRscale(i,j,k)
            ELSEWHERE(this%cbtype(i,k,:).EQ.CUSTOM_ANGKEPLER)
               pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMIN,k,:) &
                  * this%Rscale(i,j,k)
            ELSEWHERE
               ! defaults to NO_GRADIENTS
               pvar(i,Mesh%JMIN-j,k,:) = pvar(i,Mesh%JMIN,k,:)
            END WHERE
          END DO
        END DO
      END DO
    CASE(NORTH)
      DO j=1,Mesh%GNUM
        DO k=Mesh%KMIN,Mesh%KMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            WHERE (this%cbtype(i,k,:).EQ.CUSTOM_PERIOD)
                pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMIN+j-1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_REFLECT)
                pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMAX-j+1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_REFLNEG)
                pvar(i,Mesh%JMAX+j,k,:) = -pvar(i,Mesh%JMAX-j+1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_EXTRAPOL)
                pvar(i,Mesh%JMAX+j,k,:) = (j+1)*pvar(i,Mesh%JMAX,k,:) - j*pvar(i,Mesh%JMAX-1,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_FIXED)
                pvar(i,Mesh%JMAX+j,k,:) = this%data(i,j,k,:)
            ELSEWHERE (this%cbtype(i,k,:).EQ.CUSTOM_LOGEXPOL)
                pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMAX,k,:) &
                   * ABS(pvar(i,Mesh%JMAX,k,:) / pvar(i,Mesh%JMAX-1,k,:))**j
            ELSEWHERE(this%cbtype(i,k,:).EQ.CUSTOM_OUTFLOW &
                    .AND. pvar(i,Mesh%JMAX,k,:) .LE. 0.0 )
               !REFLNEG, else default (NO_GRADIENTS)
               pvar(i,Mesh%JMAX+j,k,:) = -pvar(i,Mesh%JMAX-j+1,k,:)
            ELSEWHERE(this%cbtype(i,k,:).EQ.CUSTOM_KEPLER)
               pvar(i,Mesh%JMAX+j,k,:) = (pvar(i,Mesh%JMAX,k,:) + Physics%bcradius(i,Mesh%JMAX,k)*Mesh%OMEGA)&
                  * this%invRscale(i,j,k) - Physics%bcradius(i,Mesh%JMAX+j,k)*Mesh%OMEGA
!                pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMAX,k,:) &
!                   * this%invRscale(i,j,k)
            ELSEWHERE(this%cbtype(i,k,:).EQ.CUSTOM_ANGKEPLER)
               pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMAX,k,:) &
                  * this%Rscale(i,j,k)
            ELSEWHERE
                ! defaults to NO_GRADIENTS
                pvar(i,Mesh%JMAX+j,k,:) = pvar(i,Mesh%JMAX,k,:)
            END WHERE
          END DO
        END DO
      END DO
    CASE(BOTTOM)
      DO k=1,Mesh%GNUM
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            WHERE (this%cbtype(i,k,:).EQ.CUSTOM_PERIOD)
               pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMAX-k+1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_REFLECT)
               pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMIN+k-1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_REFLNEG)
               pvar(i,j,Mesh%KMIN-k,:) = -pvar(i,j,Mesh%KMIN+k-1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_EXTRAPOL)
               pvar(i,j,Mesh%KMIN-k,:) = (k+1)*pvar(i,j,Mesh%KMIN,:) - k*pvar(i,j,Mesh%KMIN+1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_FIXED)
               pvar(i,j,Mesh%KMIN-k,:) = this%data(i,j,k,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_LOGEXPOL)
               pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMIN,:) &
                  * ABS(pvar(i,j,Mesh%KMIN,:) / pvar(i,j,Mesh%KMIN+1,:))**k
            ELSEWHERE(this%cbtype(i,j,:).EQ.CUSTOM_OUTFLOW &
                    .AND. pvar(i,j,Mesh%KMIN,:) .GE. 0.0 )
               !REFLNEG, else default (NO_GRADIENTS)
               pvar(i,j,Mesh%KMIN-k,:) = -pvar(i,j,Mesh%KMIN+k-1,:)
            ELSEWHERE(this%cbtype(i,j,:).EQ.CUSTOM_KEPLER)
               pvar(i,j,Mesh%KMIN-k,:) = (pvar(i,j,Mesh%KMIN,:) + Physics%bcradius(i,j,Mesh%KMIN)*Mesh%OMEGA)&
                  * this%invRscale(i,j,k) - Physics%bcradius(i,j,Mesh%KMIN-k)*Mesh%OMEGA
!                pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMIN,:) &
!                   * this%invRscale(i,j,k)
            ELSEWHERE(this%cbtype(i,j,:).EQ.CUSTOM_ANGKEPLER)
               pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMIN,:) &
                  * this%Rscale(i,j,k)
            ELSEWHERE
               ! defaults to NO_GRADIENTS
               pvar(i,j,Mesh%KMIN-k,:) = pvar(i,j,Mesh%KMIN,:)
            END WHERE
          END DO
        END DO
      END DO
    CASE(TOP)
      DO k=1,Mesh%GNUM
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            WHERE (this%cbtype(i,j,:).EQ.CUSTOM_PERIOD)
                pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMIN+k-1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_REFLECT)
                pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMAX-k+1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_REFLNEG)
                pvar(i,j,Mesh%KMAX+k,:) = -pvar(i,j,Mesh%KMAX-k+1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_EXTRAPOL)
                pvar(i,j,Mesh%KMAX+k,:) = (k+1)*pvar(i,j,Mesh%KMAX,:) - k*pvar(i,j,Mesh%KMAX-1,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_FIXED)
                pvar(i,j,Mesh%KMAX+k,:) = this%data(i,j,k,:)
            ELSEWHERE (this%cbtype(i,j,:).EQ.CUSTOM_LOGEXPOL)
                pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMAX,:) &
                   * ABS(pvar(i,j,Mesh%KMAX,:) / pvar(i,j,Mesh%KMAX-1,:))**k
            ELSEWHERE(this%cbtype(i,j,:).EQ.CUSTOM_OUTFLOW &
                    .AND. pvar(i,j,Mesh%KMAX,:) .LE. 0.0 )
               !REFLNEG, else default (NO_GRADIENTS)
               pvar(i,j,Mesh%KMAX+k,:) = -pvar(i,j,Mesh%KMAX-k+1,:)
            ELSEWHERE(this%cbtype(i,j,:).EQ.CUSTOM_KEPLER)
               pvar(i,j,Mesh%KMAX+k,:) = (pvar(i,j,Mesh%KMAX,:) + Physics%bcradius(i,j,Mesh%KMAX)*Mesh%OMEGA)&
                  * this%invRscale(i,j,k) - Physics%bcradius(i,j,Mesh%KMAX+k)*Mesh%OMEGA
!                pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMAX,:) &
!                   * this%invRscale(i,j,k)
            ELSEWHERE(this%cbtype(i,j,:).EQ.CUSTOM_ANGKEPLER)
               pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMAX,:) &
                  * this%Rscale(i,j,k)
            ELSEWHERE
                ! defaults to NO_GRADIENTS
                pvar(i,j,Mesh%KMAX+k,:) = pvar(i,j,Mesh%KMAX,:)
            END WHERE
          END DO
        END DO
      END DO
    END SELECT
  END SUBROUTINE SetBoundaryData

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_custom), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(this%data)) DEALLOCATE(this%data)
    IF (ASSOCIATED(this%cbtype)) DEALLOCATE(this%cbtype)
    IF (ASSOCIATED(this%Rscale)) DEALLOCATE(this%Rscale)
    IF (ASSOCIATED(this%invRscale)) DEALLOCATE(this%invRscale)

    CALL this%Finalize_base()
  END SUBROUTINE Finalize
END MODULE boundary_custom_mod
