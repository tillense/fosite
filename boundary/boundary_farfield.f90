!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_farfield.f90                                             #
!#                                                                           #
!# Copyright (C) 2006-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \author Lars Boesch
!!
!! \brief Boundary module for far field conditions
!! 
!! Implementation of inflow/outflow boundary conditions using Riemann invariants.
!!
!! \extends boundary_fixed 
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_farfield_mod
  USE marray_compound_mod
  USE mesh_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_farfield
    REAL, DIMENSION(:,:,:),   ALLOCATABLE :: Rinv, lambda
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: data
    LOGICAL                               :: first_call = .TRUE.
  CONTAINS
    PROCEDURE :: InitBoundary_farfield
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE boundary_farfield
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "far-field in-/ouflow"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    boundary_farfield
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for farfield boundary conditions
  SUBROUTINE InitBoundary_farfield(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_farfield), INTENT(INOUT) :: this
    CLASS(Mesh_base),         INTENT(IN)    :: Mesh
    CLASS(Physics_base),      INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER,  INTENT(IN)    :: config
    INTEGER,                  INTENT(IN)    :: dir
    !------------------------------------------------------------------------!
    INTEGER            :: err = 0
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,FARFIELD,boundcond_name,dir,config)
    ! check if physics supports absorbing boundary conditions
    IF (.NOT.Physics%supports_farfield) &
       CALL this%Error("InitBoundary_farfield", &
                  "boundary condition not supported for this type of physics")
    
    ! allocate memory for boundary data and mask
!CDIR IEXPAND
    SELECT CASE(this%direction%GetType())
    CASE(WEST,EAST)
       ALLOCATE(this%Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%data(Mesh%GINUM,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%data(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            STAT=err)
     CASE(BOTTOM,TOP)
       ALLOCATE(this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%GKNUM,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_farfield", "Unable to allocate memory.")
    END IF
    ! this ensures that Rinv is computed once for the data provided in this%data
    ! by the user during initialization
    this%first_call = .TRUE.
    ! reset Rinv and lambda
    this%Rinv(:,:,:) = 0.0
    this%lambda(:,:,:) = 0.0
    this%data(:,:,:,:) = 0.0
  END SUBROUTINE InitBoundary_farfield


  !> \public Applies the farfield boundary condition
  SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_farfield), INTENT(INOUT) :: this
    CLASS(Mesh_base),         INTENT(IN)    :: Mesh
    CLASS(Physics_base),      INTENT(IN)    :: Physics
    REAL,                     INTENT(IN)    :: time
    CLASS(marray_compound),   INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    !------------------------------------------------------------------------!
   SELECT CASE(this%direction%GetType())
   CASE(WEST)
     ! this must be done only once, but after the general boundary initialization
     IF (this%first_call) THEN
        ! compute Riemann invariants for the data (given in primitive variables)
        ! provided by the user in the data array
        DO i=1,Mesh%GINUM
           ! temporaryly store boundary data in ghost cells
           pvar%data4d(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
             = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
           ! compute Riemann invariants
           CALL Physics%CalculatePrim2RiemannX(Mesh,Mesh%IMIN-i,&
                                       pvar,this%lambda,this%Rinv)
           ! store Riemann invariants in data array
           this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
             = this%Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
        END DO
        ! skip the above computation for subsequent calls
        this%first_call = .FALSE.
     END IF

     DO i=1,Mesh%GINUM
       ! compute Riemann invariants
       CALL Physics%CalculatePrim2RiemannX(Mesh,Mesh%IMIN-i+1,&
                                   pvar,this%lambda,this%Rinv)
       ! set infinity Riemann invariants for inflow
       WHERE (this%lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:).GE.0.0)
         this%Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
           = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
       END WHERE
       ! transform back to primitive variables in ghost cells
       CALL Physics%CalculateRiemann2PrimX(Mesh,Mesh%IMIN-i,this%Rinv,pvar)
     END DO
   CASE(EAST)
     ! this must be done only once, but after the general boundary initialization
     IF (this%first_call) THEN
        ! compute Riemann invariants for the data (given in primitive variables)
        ! provided by the user in the data array
        DO i=1,Mesh%GINUM
           ! temporaryly store boundary data in ghost cells
           pvar%data4d(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
             = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
           ! compute Riemann invariants
           CALL Physics%CalculatePrim2RiemannX(Mesh,Mesh%IMAX+i,&
                                       pvar,this%lambda,this%Rinv)
           ! store Riemann invariants in data array
           this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
             = this%Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
        END DO
        ! skip the above computation for subsequent calls
        this%first_call = .FALSE.
     END IF

     DO i=1,Mesh%GINUM
       ! compute Riemann invariants
       CALL Physics%CalculatePrim2RiemannX(Mesh,Mesh%IMAX+i-1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:).LE.0.0)
         this%Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
           = this%data(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL Physics%CalculateRiemann2PrimX(Mesh,Mesh%IMAX+i,this%Rinv,pvar)
     END DO
   CASE(SOUTH)
     ! this must be done only once, but after the general boundary initialization
     IF (this%first_call) THEN
        ! compute Riemann invariants for the data (given in primitive variables)
        ! provided by the user in the data array
        DO j=1,Mesh%GJNUM
           ! temporaryly store boundary data in ghost cells
           pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,:) &
             = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:)
           ! compute Riemann invariants 
           CALL Physics%CalculatePrim2RiemannY(Mesh,Mesh%JMIN-j,&
                                       pvar,this%lambda,this%Rinv)
           ! store Riemann invariants in data array
           this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:) &
             = this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,:)
        END DO
        ! skip the above computation for subsequent calls
        this%first_call = .FALSE.
     END IF

     DO j=1,Mesh%GJNUM
       ! compute Riemann invariants
       CALL Physics%CalculatePrim2RiemannY(Mesh,Mesh%JMIN-j+1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow
       WHERE (this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,:).GE.0.0)
             this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,:) &
               = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL Physics%CalculateRiemann2PrimY(Mesh,Mesh%JMIN-j,this%Rinv,pvar) 
     END DO
   CASE(NORTH)
     ! this must be done only once, but after the general boundary initialization
     IF (this%first_call) THEN
        ! compute Riemann invariants for the data (given in primitive variables)
        ! provided by the user in the data array
        DO j=1,Mesh%GJNUM
           ! temporaryly store boundary data in ghost cells
           pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,:) &
             = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:)
           ! compute Riemann invariants 
           CALL Physics%CalculatePrim2RiemannY(Mesh,Mesh%JMAX+j,&
                                       pvar,this%lambda,this%Rinv)
           ! store Riemann invariants in data array
           this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:) &
             = this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,:)
        END DO
        ! skip the above computation for subsequent calls
        this%first_call = .FALSE.
     END IF

     DO j=1,Mesh%GJNUM
       ! compute Riemann invariants
       CALL Physics%CalculatePrim2RiemannY(Mesh,Mesh%JMAX+j-1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,:).LE.0.0)
         this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,:) &
           = this%data(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL Physics%CalculateRiemann2PrimY(Mesh,Mesh%JMAX+j,this%Rinv,pvar)
     END DO
   CASE(BOTTOM)
     ! this must be done only once, but after the general boundary initialization
     IF (this%first_call) THEN
        ! compute Riemann invariants for the data (given in primitive variables)
        ! provided by the user in the data array
        DO k=1,Mesh%GKNUM
           ! temporaryly store boundary data in ghost cells
           pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,:) &
             = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:)

           ! compute Riemann invariants
           CALL Physics%CalculatePrim2RiemannZ(Mesh,Mesh%KMIN-k,&
                                       pvar,this%lambda,this%Rinv)

           ! store Riemann invariants in data array
           this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:) &
             = this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:)
        END DO
        ! skip the above computation for subsequent calls
        this%first_call = .FALSE.
     END IF

     DO k=1,Mesh%GKNUM
       ! compute Riemann invariants
       CALL Physics%CalculatePrim2RiemannZ(Mesh,Mesh%KMIN-k+1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow
       WHERE (this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:).GE.0.0)
         this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:) &
           = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL Physics%CalculateRiemann2PrimZ(Mesh,Mesh%KMIN-k,this%Rinv,pvar)
     END DO
   CASE(TOP)
     ! this must be done only once, but after the general boundary initialization
     IF (this%first_call) THEN
        ! compute Riemann invariants for the data (given in primitive variables)
        ! provided by the user in the data array
        DO k=1,Mesh%GKNUM
           ! temporaryly store boundary data in ghost cells
           pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,:) &
             = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:)
           ! compute Riemann invariants
           CALL Physics%CalculatePrim2RiemannZ(Mesh,Mesh%KMAX+k,&
                                       pvar,this%lambda,this%Rinv)
           ! store Riemann invariants in data array
           this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:) &
             = this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:)
        END DO
        ! skip the above computation for subsequent calls
        this%first_call = .FALSE.
     END IF

     DO k=1,Mesh%GKNUM
       ! compute Riemann invariants
       CALL Physics%CalculatePrim2RiemannZ(Mesh,Mesh%KMAX+k-1,&
                                  pvar,this%lambda,this%Rinv)
       ! set infinity Riemanns for inflow 
       WHERE (this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:).LE.0.0)
         this%Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,:) &
           = this%data(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,:)
       END WHERE
       ! transform back to primitive variables at the boundary
       CALL Physics%CalculateRiemann2PrimZ(Mesh,Mesh%KMAX+k,this%Rinv,pvar)
     END DO
    END SELECT
  END SUBROUTINE SetBoundaryData

  !> \public Destructor for farfield boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_farfield), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%Rinv,this%lambda,this%data)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE boundary_farfield_mod
