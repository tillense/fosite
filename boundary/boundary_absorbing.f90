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
       ALLOCATE(this%xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            STAT=err)
    CASE(SOUTH,NORTH)
       ALLOCATE(this%xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
            STAT=err)
    CASE(BOTTOM,TOP)
       ALLOCATE(this%xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
            this%lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%VNUM), &
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
      ! get characteristic variables
      CALL Physics%CalculateCharSystemX(Mesh,Mesh%IMIN,Mesh%IMIN+1,pvar,this%lambda,this%xvar)
      ! set characteristic variables to zero for all incomming waves
      WHERE (this%lambda(:,:,:).GE.0.0)
        this%xvar(:,:,:) = 0.0
      END WHERE
      ! transform back to primitive variables at the boundary
      CALL Physics%CalculateBoundaryDataX(Mesh,Mesh%IMIN,Mesh%IMIN-1,this%xvar,pvar)
      ! copy data to outer ghost cells
      DO i=2,Mesh%GINUM
        pvar%data4d(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
          = pvar%data4d(Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
      END DO
    CASE(EAST)
      ! get characteristic variables
      CALL Physics%CalculateCharSystemX(Mesh,Mesh%IMAX,Mesh%IMAX-1,pvar,this%lambda,this%xvar)
      ! set characteristic variables to zero for all incomming waves
      WHERE (this%lambda(:,:,:).LE.0.0)
        this%xvar(:,:,:) = 0.0
      END WHERE
      ! transform back to primitive variables at the boundary
      CALL Physics%CalculateBoundaryDataX(Mesh,Mesh%IMAX,Mesh%IMAX+1,this%xvar,pvar)
      ! copy data to outer ghost cells
      DO i=2,Mesh%GINUM
        pvar%data4d(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) &
          = pvar%data4d(Mesh%IMAX+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
      END DO
    CASE(SOUTH)
      ! get characteristic variables
      CALL Physics%CalculateCharSystemY(Mesh,Mesh%JMIN,Mesh%JMIN+1,pvar,this%lambda,this%xvar)
      ! set characteristic variables to zero for all incomming waves
      WHERE (this%lambda(:,:,:).GE.0.0)
        this%xvar(:,:,:) = 0.0
      END WHERE
      ! transform back to primitive variables at the boundary
      CALL Physics%CalculateBoundaryDataY(Mesh,Mesh%JMIN,Mesh%JMIN-1,this%xvar,pvar)
      ! copy data to outer ghost cells
      DO j=2,Mesh%GJNUM
        pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,:) &
          = pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-1,Mesh%KMIN:Mesh%KMAX,:)
      END DO
    CASE(NORTH)
      ! get characteristic variables
      CALL Physics%CalculateCharSystemY(Mesh,Mesh%JMAX,Mesh%JMAX-1,pvar,this%lambda,this%xvar)
      ! set characteristic variables to zero for all incomming waves
      WHERE (this%lambda(:,:,:).LE.0.0)
        this%xvar(:,:,:) = 0.0
      END WHERE
      ! transform back to primitive variables at the boundary
      CALL Physics%CalculateBoundaryDataY(Mesh,Mesh%JMAX,Mesh%JMAX+1,this%xvar,pvar)
      ! copy data to outer ghost cells
      DO j=2,Mesh%GJNUM
        pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,:) &
          = pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1,Mesh%KMIN:Mesh%KMAX,:)
      END DO
    CASE(BOTTOM)
      ! get characteristic variables
      CALL Physics%CalculateCharSystemZ(Mesh,Mesh%KMIN,Mesh%KMIN+1,pvar,this%lambda,this%xvar)
      ! set characteristic variables to zero for all incomming waves
      WHERE (this%lambda(:,:,:).GE.0.0)
        this%xvar(:,:,:) = 0.0
      END WHERE
      ! transform back to primitive variables at the boundary
      CALL Physics%CalculateBoundaryDataZ(Mesh,Mesh%KMIN,Mesh%KMIN-1,this%xvar,pvar)
      ! copy data to outer ghost cells
      DO k=2,Mesh%GKNUM
        pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,:) &
          = pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-1,:)
      END DO
    CASE(TOP)
      ! get characteristic variables
      CALL Physics%CalculateCharSystemZ(Mesh,Mesh%KMAX,Mesh%KMAX-1,pvar,this%lambda,this%xvar)
      ! set characteristic variables to zero for all incomming waves
      WHERE (this%lambda(:,:,:).LE.0.0)
        this%xvar(:,:,:) = 0.0
      END WHERE
      ! transform back to primitive variables at the boundary
      CALL Physics%CalculateBoundaryDataZ(Mesh,Mesh%KMAX,Mesh%KMAX+1,this%xvar,pvar)
      ! copy data to outer ghost cells
      DO k=2,Mesh%GKNUM
        pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,:) &
          = pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+1,:)
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
