!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_absorbing.f90                                            #
!#                                                                           #
!# Copyright (C) 2009-2014                                                   #
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
!! \brief Boundary module for absorbing (non-reflecting) conditions
!!
!! This module uses characteristic variables and wave speeds at the
!! boundary to determine the state of the flow. Depending on this it
!! damps oszillations by setting the characteristic variables to zero
!! for imcomming waves.
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_absorbing_mod
  USE marray_compound_mod
  USE mesh_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_absorbing
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: xvar, & !< characteristic variables for absorbing bc
                                         lambda    !< eigenvalues for absorbing bc
  CONTAINS
    PROCEDURE :: InitBoundary_absorbing
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE boundary_absorbing
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "absorbing"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    boundary_absorbing
    !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for absorbing boundary conditions
  !!
  !! Initilizes the boundary condition type and direction and allocates
  !! memory for characteristic variables and wave speeds.
  SUBROUTINE InitBoundary_absorbing(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_absorbing), INTENT(INOUT) :: this
    CLASS(Mesh_base),          INTENT(IN)    :: Mesh
    CLASS(Physics_base),       INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER,   INTENT(IN)    :: config
    INTEGER,                   INTENT(IN)    :: dir
    !------------------------------------------------------------------------!
    INTEGER                                  :: err
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,ABSORBING,boundcond_name,dir,config)
    ! check if physics supports absorbing boundary conditions
    IF (.NOT.Physics%supports_absorbing) &
       CALL this%Error("InitBoundary_absorbing", &
                  "boundary condition not supported for this type of physics")
    SELECT CASE(this%direction%GetType())
    CASE(WEST,EAST)
       ALLOCATE(this%xvar(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            this%lambda(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%xvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            this%lambda(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
            STAT=err)
    CASE(BOTTOM,TOP)
       ALLOCATE(this%xvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            this%lambda(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM), &
            STAT=err)
    END SELECT
    IF (err.NE.0) &
         CALL this%Error("InitBoundary_absorbing", "Unable to allocate memory.")
    ! initialize the data
    this%xvar(:,:,:) = 0.0
    this%lambda(:,:,:) = 0.0
  END SUBROUTINE InitBoundary_absorbing


  !> \public Applies the absorbing boundary condition
  !!
  !! This is an implementation of characteristic variable extrapolation.
  !! The algorithm first computes the characteristic (pseudo-) variables at
  !! the boundary and then sets them to zero for incomming waves, i. e.
  !! for those associates with positive (western/southern) or negative
  !! (eastern/northern) wave speeds. After that it transforms the
  !! new set of characteristic variables back to primitive variables.
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_absorbing),INTENT(INOUT) :: this
    CLASS(Mesh_base),         INTENT(IN) :: Mesh
    CLASS(Physics_base),      INTENT(IN) :: Physics
    REAL,                     INTENT(IN) :: time
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    SELECT CASE(this%direction%GetType())
    CASE(WEST)
       DO i=1,Mesh%GINUM
          CALL Physics%CalculateCharSystemX(Mesh,Mesh%IMIN-i+1,+1,pvar%data4d,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:,:).GE.0.0)
             this%xvar(:,:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL Physics%CalculateBoundaryDataX(Mesh,Mesh%IMIN-i+1,-1,this%xvar,pvar%data4d)
       END DO
    CASE(EAST)
       ! characteristic variables
       DO i=1,Mesh%GINUM
          CALL Physics%CalculateCharSystemX(Mesh,Mesh%IMAX+i-1,-1,pvar%data4d,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:,:).LE.0.0)
             this%xvar(:,:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL Physics%CalculateBoundaryDataX(Mesh,Mesh%IMAX+i-1,+1,this%xvar,pvar%data4d)
       END DO
     CASE(SOUTH)
       DO j=1,Mesh%GJNUM
          CALL Physics%CalculateCharSystemY(Mesh,Mesh%JMIN-j+1,+1,pvar%data4d,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:,:).GE.0.0)
             this%xvar(:,:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL Physics%CalculateBoundaryDataY(Mesh,Mesh%JMIN-j+1,-1,this%xvar,pvar%data4d)
       END DO
    CASE(NORTH)
       DO j=1,Mesh%GJNUM
          CALL Physics%CalculateCharSystemY(Mesh,Mesh%JMAX+j-1,-1,pvar%data4d,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:,:).LE.0.0)
             this%xvar(:,:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL Physics%CalculateBoundaryDataY(Mesh,Mesh%JMAX+j-1,+1,this%xvar,pvar%data4d)
       END DO   
     CASE(BOTTOM)
       DO k=1,Mesh%GKNUM
          CALL Physics%CalculateCharSystemZ(Mesh,Mesh%KMIN-k+1,+1,pvar%data4d,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:,:).GE.0.0)
             this%xvar(:,:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL Physics%CalculateBoundaryDataZ(Mesh,Mesh%KMIN-k+1,-1,this%xvar,pvar%data4d)
       END DO
    CASE(TOP)
       DO k=1,Mesh%GKNUM
          CALL Physics%CalculateCharSystemZ(Mesh,Mesh%KMAX+k-1,-1,pvar%data4d,this%lambda,this%xvar)
          ! set characteristic variables to zero for all incomming waves
          WHERE (this%lambda(:,:,:).LE.0.0)
             this%xvar(:,:,:) = 0.0
          END WHERE
          ! transform back to primitive variables at the boundary
          CALL Physics%CalculateBoundaryDataZ(Mesh,Mesh%KMAX+k-1,+1,this%xvar,pvar%data4d)
       END DO
     END SELECT
  END SUBROUTINE SetBoundaryData


  !> \public Destructor for absorbing boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_absorbing), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%xvar,this%lambda)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE boundary_absorbing_mod
