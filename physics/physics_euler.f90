!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_euler.f90                                                 #
!#                                                                           #
!# Copyright (C) 2007-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup physics
!! - non-isothermal gas dynamics
!!   \key{gamma,REAL,ratio of specific heats (default is for diatomic
!!      molecular gas),1.4}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Jannes Klee
!!
!! \brief physics module for 1D,2D and 3D non-isothermal Euler equations
!!
!! \extends physics_eulerisotherm
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler_mod
  USE logging_base_mod
  USE physics_base_mod
  USE physics_eulerisotherm_mod, ONLY: physics_eulerisotherm, statevector_eulerisotherm
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler"
  !--------------------------------------------------------------------------!
  TYPE,  EXTENDS(physics_eulerisotherm) :: physics_euler
    REAL                :: gamma                 !< ratio of spec. heats
  CONTAINS
    PROCEDURE :: InitPhysics             !< constructor
    PROCEDURE :: new_statevector
    !------Convert2Primitve--------!
    PROCEDURE :: Convert2Primitive_all
    PROCEDURE :: Convert2Primitive_subset
    !------Convert2Conservative----!
    PROCEDURE :: Convert2Conservative_all
    PROCEDURE :: Convert2Conservative_subset
    !------soundspeed routines-----!
    PROCEDURE :: UpdateSoundSpeed
    !------flux routines-----------!
    PROCEDURE :: CalcFluxesX
    PROCEDURE :: CalcFluxesY
    PROCEDURE :: CalcFluxesZ
    !------fargo routines----------!
    PROCEDURE :: AddBackgroundVelocityX
    PROCEDURE :: AddBackgroundVelocityY
    PROCEDURE :: AddBackgroundVelocityZ
    PROCEDURE :: SubtractBackgroundVelocityX
    PROCEDURE :: SubtractBackgroundVelocityY
    PROCEDURE :: SubtractBackgroundVelocityZ
    PROCEDURE :: AddFargoSourcesX
    PROCEDURE :: AddFargoSourcesY
    PROCEDURE :: AddFargoSourcesZ
    !------HLLC routines-----------!
!    PROCEDURE :: CalcIntermediateStateX
!    PROCEDURE :: CalcIntermediateStateY

    PROCEDURE :: ExternalSources
    PROCEDURE :: GeometricalSources
    PROCEDURE :: ViscositySources

    ! boundarie routines
    PROCEDURE :: CalculateCharSystemX          ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemY          ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemZ          ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataX        ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataY        ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataZ        ! for absorbing boundaries
    PROCEDURE :: CalculatePrim2RiemannX        ! for farfield boundaries
    PROCEDURE :: CalculatePrim2RiemannY        ! for farfield boundaries
    PROCEDURE :: CalculatePrim2RiemannZ        ! for farfield boundaries
    PROCEDURE :: CalculateRiemann2PrimX        ! for farfield boundaries
    PROCEDURE :: CalculateRiemann2PrimY        ! for farfield boundaries
    PROCEDURE :: CalculateRiemann2PrimZ        ! for farfield boundaries

    PROCEDURE :: Finalize
  END TYPE
  TYPE, EXTENDS(statevector_eulerisotherm) :: statevector_euler
    TYPE(marray_base), POINTER :: &
                               pressure => null(), &
                               energy => null()
    CONTAINS
    PROCEDURE :: AssignMArray_0
    FINAL     :: Finalize_statevector
  END TYPE
  INTERFACE statevector_euler
    MODULE PROCEDURE CreateStateVector
  END INTERFACE
  INTERFACE SetFlux
    MODULE PROCEDURE SetFlux1d, SetFlux2d, SetFlux3d
  END INTERFACE
  INTERFACE Cons2Prim
    MODULE PROCEDURE Cons2Prim1d, Cons2Prim2d, Cons2Prim3d
  END INTERFACE
  INTERFACE Prim2Cons
    MODULE PROCEDURE Prim2Cons1d, Prim2Cons2d, Prim2Cons3d
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       physics_euler, &
       statevector_euler
  !--------------------------------------------------------------------------!

CONTAINS

  !> constructor of physics_euler class
  SUBROUTINE InitPhysics(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN) :: Mesh
    TYPE(Dict_TYP), POINTER             :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                             :: err,next_idx
    CHARACTER(LEN=64)                   :: info_str
    !------------------------------------------------------------------------!
    ! call InitPhysics from base class
    CALL this%InitPhysics_base(Mesh,config,IO,EULER,problem_name)

    ! set the total number of variables in a state vector
    this%VNUM = this%VDIM + 2

    ! ratio of specific heats
    CALL GetAttr(config, "gamma", this%gamma, 1.4)
    
    ! allocate memory for arrays used in euler
    ALLOCATE(this%pvarname(this%VNUM),this%cvarname(this%VNUM),this%bccsound, &
             this%fcsound, &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitPhysics_euler", "Unable to allocate memory.")

    !> \todo remove / improve in future version
    !! set array indices for 1st,2nd,3rd non-vanishing velocities
    !! this may actually not coincide with the x,y and z-velocities
    this%DENSITY   = 1                                 ! mass density        !
    this%pvarname(this%DENSITY)   = "density"
    this%cvarname(this%DENSITY)   = "density"
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    IF (this%VDIM.GE.2) THEN
      this%YVELOCITY = 3                               ! y-velocity          !
      this%YMOMENTUM = 3                               ! y-momentum          !
    ELSE
      this%YVELOCITY = 0                               ! no y-velocity       !
      this%YMOMENTUM = 0                               ! no y-momentum       !
    END IF
    IF (this%VDIM.EQ.3) THEN
      this%ZVELOCITY = 4                               ! z-velocity          !
      this%ZMOMENTUM = 4                               ! z-momentum          !
    ELSE
      this%ZVELOCITY = 0                               ! no z-velocity       !
      this%ZMOMENTUM = 0                               ! no z-momentum       !
    END IF
    this%PRESSURE  = this%VNUM                         ! pressure            !
    this%ENERGY    = this%VNUM                         ! total energy        !
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%ENERGY)    = "energy"

    ! check which vector components are available and
    ! set names shown in the data file
    next_idx = 2
    IF (BTEST(Mesh%VECTOR_COMPONENTS,0)) THEN
      this%pvarname(next_idx) = "xvelocity"
      this%cvarname(next_idx) = "xmomentum"
      next_idx = next_idx + 1
    END IF
    IF (BTEST(Mesh%VECTOR_COMPONENTS,1)) THEN
      this%pvarname(next_idx) = "yvelocity"
      this%cvarname(next_idx) = "ymomentum"
      next_idx = next_idx + 1
    END IF
    IF (BTEST(Mesh%VECTOR_COMPONENTS,2)) THEN
      this%pvarname(next_idx) = "zvelocity"
      this%cvarname(next_idx) = "zmomentum"
    END IF

    ! not used in non-isotherml physics
    this%csiso = 0.0

    ! create new mesh arrays for sound speeds
    this%bccsound = marray_base()
    this%fcsound = marray_base(Mesh%NFACES)

    ! enable support for absorbing and farfield boundary conditions
    this%supports_absorbing = .TRUE.
    this%supports_farfield  = .TRUE.

    CALL this%SetOutput(Mesh,config,IO)

    ! print some information
    WRITE(info_str,'("spec. heat ratio:  ",F4.2)') this%gamma
    CALL this%Info(REPEAT(" ",12) // TRIM(info_str))
  END SUBROUTINE InitPhysics

  !> \public allocate and initialize new non-isothermal state vector
  SUBROUTINE new_statevector(this,new_sv,flavour,num)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(marray_compound), POINTER :: new_sv
    INTEGER, OPTIONAL, INTENT(IN)   :: flavour,num
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(new_sv)) THEN
      DEALLOCATE(new_sv)
#ifdef DEBUG
      CALL this%Warning("physics_euler::new_statevector","new statevector already associated")
#endif
    END IF
    ALLOCATE(statevector_euler::new_sv)
    new_sv = statevector_euler(this,flavour,num)
  END SUBROUTINE new_statevector

  !> Converts conservative to primitive variables on the whole mesh
  PURE SUBROUTINE Convert2Primitive_all(this,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)      :: this
    CLASS(marray_compound), INTENT(INOUT) :: cvar,pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    TYPE IS (statevector_euler)
      SELECT TYPE(p => pvar)
      TYPE IS (statevector_euler)
        IF (c%flavour.EQ.CONSERVATIVE.AND.p%flavour.EQ.PRIMITIVE) THEN
          ! conservative -> primitive
          SELECT CASE(this%VDIM)
          CASE(1)
            CALL Cons2Prim(this%gamma,c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%energy%data1d(:),p%density%data1d(:), &
                          p%velocity%data2d(:,1),p%pressure%data1d(:))
          CASE(2)
            CALL Cons2Prim(this%gamma,c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%momentum%data2d(:,2),c%energy%data1d(:), &
                          p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%velocity%data2d(:,2),p%pressure%data1d(:))
          CASE(3)
            CALL Cons2Prim(this%gamma,c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%momentum%data2d(:,2),c%momentum%data2d(:,3), &
                          c%energy%data1d(:),p%density%data1d(:), &
                          p%velocity%data2d(:,1),p%velocity%data2d(:,2),&
                          p%velocity%data2d(:,3),p%pressure%data1d(:))
          END SELECT
          p%fargo_transformation_applied = c%fargo_transformation_applied
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Primitive_all

  !> Converts conservative to primitive variables on a subset of the data
  PURE SUBROUTINE Convert2Primitive_subset(this,i1,i2,j1,j2,k1,k2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),      INTENT(IN) :: this
    INTEGER,                   INTENT(IN) :: i1,i2,j1,j2,k1,k2
    CLASS(marray_compound), INTENT(INOUT) :: cvar,pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    TYPE IS (statevector_euler)
      SELECT TYPE(p => pvar)
      TYPE IS (statevector_euler)
        IF (c%flavour.EQ.CONSERVATIVE.AND.p%flavour.EQ.PRIMITIVE) THEN
          SELECT CASE (c%density%RANK)
          CASE(0) ! state vector contains cell center values
            ! perform the transformation depending on dimensionality
            SELECT CASE(this%VDIM)
            CASE(1)
              CALL Cons2Prim(this%gamma,c%density%data3d(i1:i2,j1:j2,k1:k2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,1), &
                            c%energy%data3d(i1:i2,j1:j2,k1:k2), &
                            p%density%data3d(i1:i2,j1:j2,k1:k2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,1), &
                            p%pressure%data3d(i1:i2,j1:j2,k1:k2))
            CASE(2)
              CALL Cons2Prim(this%gamma,c%density%data3d(i1:i2,j1:j2,k1:k2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,1), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,2), &
                            c%energy%data3d(i1:i2,j1:j2,k1:k2), &
                            p%density%data3d(i1:i2,j1:j2,k1:k2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,1), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,2), &
                            p%pressure%data3d(i1:i2,j1:j2,k1:k2))
            CASE(3)
              CALL Cons2Prim(this%gamma,c%density%data3d(i1:i2,j1:j2,k1:k2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,1), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,3), &
                            c%energy%data3d(i1:i2,j1:j2,k1:k2), &
                            p%density%data3d(i1:i2,j1:j2,k1:k2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,1), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,2),&
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,3), &
                            p%pressure%data3d(i1:i2,j1:j2,k1:k2))
            END SELECT
          CASE(1) ! state vector contains cell face / corner values
            ! perform the transformation depending on dimensionality
            SELECT CASE(this%VDIM)
            CASE(1)
              CALL Cons2Prim(this%gamma,c%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            c%energy%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            p%pressure%data4d(i1:i2,j1:j2,k1:k2,:))
            CASE(2)
              CALL Cons2Prim(this%gamma,c%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            c%energy%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            p%pressure%data4d(i1:i2,j1:j2,k1:k2,:))
            CASE(3)
              CALL Cons2Prim(this%gamma,c%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,3), &
                            c%energy%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,2),&
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,3), &
                            p%pressure%data4d(i1:i2,j1:j2,k1:k2,:))
            END SELECT
          CASE DEFAULT
            ! do nothing
          END SELECT
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Primitive_subset

  !> Converts primitive to conservative variables on the whole mesh
  PURE SUBROUTINE Convert2Conservative_all(this,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)      :: this
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS (statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS (statevector_euler)
        IF (p%flavour.EQ.PRIMITIVE.AND.c%flavour.EQ.CONSERVATIVE) THEN
          ! perform the transformation depending on dimensionality
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D velocity / momentum
            CALL Prim2Cons(this%gamma,p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%pressure%data1d(:),c%density%data1d(:), &
                          c%momentum%data2d(:,1),c%energy%data1d(:))
          CASE(2) ! 2D velocity / momentum
            CALL Prim2Cons(this%gamma,p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%velocity%data2d(:,2),p%pressure%data1d(:), &
                          c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%momentum%data2d(:,2),c%energy%data1d(:))
          CASE(3) ! 3D velocity / momentum
            CALL Prim2Cons(this%gamma,p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%velocity%data2d(:,2),p%velocity%data2d(:,3), &
                          p%pressure%data1d(:),c%density%data1d(:), &
                          c%momentum%data2d(:,1),c%momentum%data2d(:,2), &
                          c%momentum%data2d(:,3),c%energy%data1d(:))
          END SELECT
          c%fargo_transformation_applied = p%fargo_transformation_applied
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Conservative_all

  !> Converts primitive to conservative variables on a subset of the data
  PURE SUBROUTINE Convert2Conservative_subset(this,i1,i2,j1,j2,k1,k2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),      INTENT(IN) :: this
    INTEGER,                   INTENT(IN) :: i1,i2,j1,j2,k1,k2
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS (statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS (statevector_euler)
        IF (p%flavour.EQ.PRIMITIVE.AND.c%flavour.EQ.CONSERVATIVE) THEN
          SELECT CASE (p%density%RANK)
          CASE(0) ! state vector contains cell center values
            ! perform the transformation depending on dimensionality
            SELECT CASE(this%VDIM)
            CASE(1) ! 1D velocity / momentum
              CALL Prim2Cons(this%gamma,p%density%data3d(i1:i2,j1:j2,k1:k2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,1), &
                            p%pressure%data3d(i1:i2,j1:j2,k1:k2), &
                            c%density%data3d(i1:i2,j1:j2,k1:k2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,1), &
                            c%energy%data3d(i1:i2,j1:j2,k1:k2))
            CASE(2) ! 2D velocity / momentum
              CALL Prim2Cons(this%gamma,p%density%data3d(i1:i2,j1:j2,k1:k2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,1), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,2), &
                            p%pressure%data3d(i1:i2,j1:j2,k1:k2), &
                            c%density%data3d(i1:i2,j1:j2,k1:k2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,1), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,2), &
                            c%energy%data3d(i1:i2,j1:j2,k1:k2))
            CASE(3) ! 3D velocity / momentum
              CALL Prim2Cons(this%gamma,p%density%data3d(i1:i2,j1:j2,k1:k2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,1), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,2), &
                            p%velocity%data4d(i1:i2,j1:j2,k1:k2,3), &
                            p%pressure%data3d(i1:i2,j1:j2,k1:k2), &
                            c%density%data3d(i1:i2,j1:j2,k1:k2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,1), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,2), &
                            c%momentum%data4d(i1:i2,j1:j2,k1:k2,3), &
                            c%energy%data3d(i1:i2,j1:j2,k1:k2))
            END SELECT
          CASE(1) ! state vector contains cell face / corner values
            ! perform the transformation depending on dimensionality
            SELECT CASE(this%VDIM)
            CASE(1) ! 1D velocity / momentum
              CALL Prim2Cons(this%gamma,p%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            p%pressure%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            c%energy%data4d(i1:i2,j1:j2,k1:k2,:))
            CASE(2) ! 2D velocity / momentum
              CALL Prim2Cons(this%gamma,p%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            p%pressure%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            c%energy%data4d(i1:i2,j1:j2,k1:k2,:))
            CASE(3) ! 3D velocity / momentum
              CALL Prim2Cons(this%gamma,p%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            p%velocity%data5d(i1:i2,j1:j2,k1:k2,:,3), &
                            p%pressure%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%density%data4d(i1:i2,j1:j2,k1:k2,:), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,1), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,2), &
                            c%momentum%data5d(i1:i2,j1:j2,k1:k2,:,3), &
                            c%energy%data4d(i1:i2,j1:j2,k1:k2,:))
            END SELECT
          CASE DEFAULT
            ! do nothing
          END SELECT
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Conservative_subset

  !> Calculate Fluxes in x-direction
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalcFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: nmin,nmax
    CLASS(marray_compound), INTENT(INOUT) :: prim,cons,xfluxes
    !------------------------------------------------------------------------!
    SELECT TYPE(p => prim)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cons)
      TYPE IS(statevector_euler)
        SELECT TYPE(f => xfluxes)
        TYPE IS(statevector_euler)
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%energy%data2d(:,nmin:nmax))
          CASE(2) ! 2D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%energy%data2d(:,nmin:nmax))
          CASE(3) ! 3D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,3),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,3),    &
                      f%energy%data2d(:,nmin:nmax))
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcFluxesX

  !> Calculate Fluxes in y-direction
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalcFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: nmin,nmax
    CLASS(marray_compound), INTENT(INOUT) :: prim,cons,yfluxes
    !------------------------------------------------------------------------!
    SELECT TYPE(p => prim)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cons)
      TYPE IS(statevector_euler)
        SELECT TYPE(f => yfluxes)
        TYPE IS(statevector_euler)
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%energy%data2d(:,nmin:nmax))
          CASE(2) ! 2D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,2),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%energy%data2d(:,nmin:nmax))
          CASE(3) ! 3D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,2),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,3),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,3),    &
                      f%energy%data2d(:,nmin:nmax))
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcFluxesY

  !> Calculate Fluxes in z-direction
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalcFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: nmin,nmax
    CLASS(marray_compound), INTENT(INOUT) :: prim,cons,zfluxes
    !------------------------------------------------------------------------!
    SELECT TYPE(p => prim)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cons)
      TYPE IS(statevector_euler)
        SELECT TYPE(f => zfluxes)
        TYPE IS(statevector_euler)
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%energy%data2d(:,nmin:nmax))
          CASE(2) ! 2D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,2),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%energy%data2d(:,nmin:nmax))
          CASE(3) ! 3D flux
            CALL SetFlux(p%density%data2d(:,nmin:nmax),    &
                      p%velocity%data3d(:,nmin:nmax,3),    &
                      p%pressure%data2d(:,nmin:nmax),      &
                      c%momentum%data3d(:,nmin:nmax,3),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%energy%data2d(:,nmin:nmax),        &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,3),    &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%energy%data2d(:,nmin:nmax))
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcFluxesZ

!  !> Reconstruction of the intermediate state for HLLC
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcIntermediateStateX(this,Mesh,prim,cons, &
!                  amin,amax,cstar,astar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM) &
!                                       :: prim,cons
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: cstar
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
!                                       :: amin,amax
!    REAL,                  INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
!                                       :: astar
!    !------------------------------------------------------------------------!
!    INTEGER                            :: i,j,k
!    !------------------------------------------------------------------------!
!    DO k=Mesh%KMIN,Mesh%KMAX
!      DO j=Mesh%JMIN,Mesh%JMAX
!!NEC$ IVDEP
!        DO i=Mesh%IMIN-1,Mesh%IMAX
!          CALL SetIntermediateState( &
!               prim(i,j,k,2,this%DENSITY),prim(i+1,j,k,1,this%DENSITY), &
!               prim(i,j,k,2,this%XVELOCITY),prim(i+1,j,k,1,this%XVELOCITY), &
!               prim(i,j,k,2,this%YVELOCITY),prim(i+1,j,k,1,this%YVELOCITY), &
!               prim(i,j,k,2,this%PRESSURE),prim(i+1,j,k,1,this%PRESSURE), &
!               cons(i,j,k,2,this%ENERGY),cons(i+1,j,k,1,this%ENERGY), &
!               amin(i,j,k),amax(i,j,k), &
!               cstar(i,j,k,this%DENSITY),cstar(i,j,k,this%XMOMENTUM), &
!               cstar(i,j,k,this%YMOMENTUM),cstar(i,j,k,this%ENERGY), &
!               astar(i,j,k))
!         END DO
!       END DO
!    END DO
!  END SUBROUTINE CalcIntermediateStateX
!
!
!  !> Reconstruction of the intermediate state for HLLC
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcIntermediateStateY(this,Mesh,prim,cons, &
!                  bmin,bmax,cstar,bstar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(INOUT) :: this
!    CLASS(mesh_base),       INTENT(IN)    :: Mesh
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM) &
!                                          :: prim,cons
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                          :: cstar
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
!                                          :: bmin,bmax
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
!                                          :: bstar
!    !------------------------------------------------------------------------!
!    INTEGER                               :: i,j,k
!    !------------------------------------------------------------------------!
!    DO k=Mesh%KGMIN,Mesh%KGMAX
!      DO j=Mesh%JMIN-1,Mesh%JMAX
!!NEC$ IVDEP
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          CALL SetIntermediateState( &
!               prim(i,j,k,4,this%DENSITY),prim(i,j+1,k,3,this%DENSITY),     &
!               prim(i,j,k,4,this%YVELOCITY),prim(i,j+1,k,3,this%YVELOCITY), &
!               prim(i,j,k,4,this%XVELOCITY),prim(i,j+1,k,3,this%XVELOCITY), &
!               prim(i,j,k,4,this%PRESSURE),prim(i,j+1,k,3,this%PRESSURE),   &
!               cons(i,j,k,4,this%ENERGY),cons(i,j+1,k,3,this%ENERGY),       &
!               bmin(i,j,k),bmax(i,j,k),                                     &
!               cstar(i,j,k,this%DENSITY),cstar(i,j,k,this%YMOMENTUM),       &
!               cstar(i,j,k,this%XMOMENTUM),cstar(i,j,k,this%ENERGY),        &
!               bstar(i,j,k))
!        END DO
!      END DO
!    END DO
!  END SUBROUTINE CalcIntermediateStateY

  PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i1,i2,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),         INTENT(IN)    :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    INTEGER,                      INTENT(IN)    :: i1,i2
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM), &
                                  INTENT(OUT)   :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: iL,iR
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      iL = MIN(i1,i2)
      iR = MAX(i1,i2)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        ! compute eigenvalues at i
        CALL SetEigenValues1d(this%gamma, &
              p%density%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars1d(this%gamma, &
              p%density%data3d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%density%data3d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%pressure%data3d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
      CASE(2) ! 2D
        ! compute eigenvalues at i1
        CALL SetEigenValues2d(this%gamma, &
              p%density%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars2d(this%gamma, &
              p%density%data3d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%density%data3d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%pressure%data3d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%pressure%data3d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
      CASE(3) ! 3D
        ! compute eigenvalues at i1
        CALL SetEigenValues3d(this%gamma, &
              p%density%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars3d(this%gamma, &
              p%density%data3d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%density%data3d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              p%velocity%data4d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              p%pressure%data3d(iL,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%pressure%data3d(iR,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateCharSystemX

  PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j1,j2,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),         INTENT(IN)    :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    INTEGER,                      INTENT(IN)    :: j1,j2
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM), &
                                  INTENT(OUT)   :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: jL,jR,vn,vt
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      jL = MIN(j1,j2)
      jR = MAX(j1,j2)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        ! compute eigenvalues at j
        CALL SetEigenValues1d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,2), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars1d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
      CASE(2) ! 2D
        ! check which velocity component is along y-direction
        SELECT CASE(Mesh%VECTOR_COMPONENTS)
        CASE(IOR(VECTOR_X,VECTOR_Y)) ! 2D velocities in x-y-plane
          vt = 1
          vn = 2
        CASE(IOR(VECTOR_Y,VECTOR_Z)) ! 2D velocities in y-z-plane
          vt = 2
          vn = 1
        CASE DEFAULT
          ! this should not happen
          RETURN
        END SELECT
        ! compute eigenvalues at j
        CALL SetEigenValues2d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,vn), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars2d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX,vn), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX,vn), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX,vt), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX,vt), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4))
      CASE(3) ! 3D
        ! compute eigenvalues at j
        CALL SetEigenValues3d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,2), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars3d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX,3), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX,3), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,jL,Mesh%KMIN:Mesh%KMAX), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,jR,Mesh%KMIN:Mesh%KMAX), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateCharSystemY


  PURE SUBROUTINE CalculateCharSystemZ(this,Mesh,k1,k2,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),         INTENT(IN)    :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    INTEGER,                      INTENT(IN)    :: k1,k2
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,this%VNUM), &
                                  INTENT(OUT)   :: lambda,xvar
    !------------------------------------------------------------------------!
    INTEGER           :: kL,kR
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      kL = MIN(k1,k2)
      kR = MAX(k1,k2)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        ! compute eigenvalues at k
        CALL SetEigenValues1d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,3), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars1d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
      CASE(2) ! 2D
        ! compute eigenvalues at k
        CALL SetEigenValues2d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,3), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars2d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL,2), & ! 2nd component is vz
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL,1), & ! 1st component: vx or vy
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4))
      CASE(3) ! 3D
        ! compute eigenvalues at k
        CALL SetEigenValues3d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,3), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))
        ! compute characteristic variables using cell mean values of adjacent
        ! cells to calculate derivatives and the isothermal speed of sound
        ! at the intermediate cell face
        CALL SetCharVars3d(this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL,3), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR,3), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR,2), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kL), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,kR), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateCharSystemZ

  !> extrapolate pvar using characteristic pseudo variables (absorbing boundaries)
  PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i1,i2,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),         INTENT(IN)    :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    INTEGER,                      INTENT(IN)    :: i1,i2
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM), &
                                  INTENT(IN)    :: xvar
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL SetBoundaryData1d(i2-i1,this%gamma, &
              p%density%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              p%density%data3d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CASE(2) ! 2D
        CALL SetBoundaryData2d(i2-i1,this%gamma, &
              p%density%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%pressure%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
              p%density%data3d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%pressure%data3d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CASE(3) ! 3D
        CALL SetBoundaryData3d(i2-i1,this%gamma, &
              p%density%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              p%pressure%data3d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
              xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5), &
              p%density%data3d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
              p%pressure%data3d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateBoundaryDataX

  PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j1,j2,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),         INTENT(IN)    :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    INTEGER,                      INTENT(IN)    :: j1,j2
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM), &
                                  INTENT(IN)    :: xvar
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: vt,vn
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL SetBoundaryData1d(j2-j1,this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX))
      CASE(2) ! 2D
        ! check which velocity component is along y-direction
        SELECT CASE(Mesh%VECTOR_COMPONENTS)
        CASE(IOR(VECTOR_X,VECTOR_Y)) ! 2D velocities in x-y-plane
          vt = 1
          vn = 2
        CASE(IOR(VECTOR_Y,VECTOR_Z)) ! 2D velocities in y-z-plane
          vt = 2
          vn = 1
        CASE DEFAULT
          ! this should not happen
          RETURN
        END SELECT
        CALL SetBoundaryData2d(j2-j1,this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,vn), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,vt), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,vn), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,vt), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX))
      CASE(3) ! 3D
        CALL SetBoundaryData3d(j2-j1,this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,3), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,3), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateBoundaryDataY

  PURE SUBROUTINE CalculateBoundaryDataZ(this,Mesh,k1,k2,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler),         INTENT(IN)    :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    INTEGER,                      INTENT(IN)    :: k1,k2
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,this%VNUM), &
                                  INTENT(IN)    :: xvar
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL SetBoundaryData1d(k2-k1,this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2))
      CASE(2) ! 2D
        CALL SetBoundaryData2d(k2-k1,this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,1), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2))
      CASE(3) ! 3D
        CALL SetBoundaryData3d(k2-k1,this%gamma, &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,3), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,2), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
              xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
              p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,3), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,1), &
              p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,2), &
              p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateBoundaryDataZ


  !> Conversion from primitive to riemann invariants for farfield boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculatePrim2RiemannX(this,Mesh,i,pvar,lambda,Rinv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: i
    CLASS(marray_compound), INTENT(IN) :: pvar
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM) &
                                       :: lambda
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM) &
                                       :: Rinv
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL SetEigenValues1d(this%gamma, &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),&
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
        ! compute Riemann invariants
        CALL Prim2Riemann1d(this%gamma, &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
      CASE(2) ! 2D
        CALL SetEigenValues2d(this%gamma, &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),&
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
        ! compute Riemann invariants
        CALL Prim2Riemann2d(this%gamma, &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
      CASE(3) ! 3D
        ! compute eigenvalues at i
        CALL SetEigenValues3d(this%gamma, &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),&
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5))
        ! compute Riemann invariants
        CALL Prim2Riemann3d(this%gamma, &
          p%density%data3d( i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5))
      END SELECT
    END SELECT
  END SUBROUTINE CalculatePrim2RiemannX


  !> Conversion from primitive to riemann invariants for farfield boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculatePrim2RiemannY(this,Mesh,j,pvar,lambda,Rinv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: j
    CLASS(marray_compound), INTENT(IN) :: pvar
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM) &
                                       :: lambda
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM) &
                                       :: Rinv
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
         ! compute eigenvalues at j
        CALL SetEigenValues1d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX),&
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
        ! compute Riemann invariants
        CALL Prim2Riemann1d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
      CASE(2) ! 2D
        ! compute eigenvalues at j
        CALL SetEigenValues2d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,2), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX),&
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4))
        ! compute Riemann invariants
        CALL Prim2Riemann2d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,2), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4))
      CASE(3) ! 3D
        ! compute eigenvalues at j
        CALL SetEigenValues3d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,2), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX),&
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5))
        ! compute Riemann invariants
        CALL Prim2Riemann3d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,2), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,3), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5))
      END SELECT
    END SELECT
  END SUBROUTINE CalculatePrim2RiemannY


  !> Conversion from primitive to riemann invariants for farfield boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculatePrim2RiemannZ(this,Mesh,k,pvar,lambda,Rinv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: k
    CLASS(marray_compound), INTENT(IN) :: pvar
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,this%VNUM) &
                                       :: lambda
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,this%VNUM) &
                                       :: Rinv
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        ! compute eigenvalues at k
        CALL SetEigenValues1d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k),&
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3))
        ! compute Riemann invariants
        CALL Prim2Riemann1d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3))
      CASE(2) ! 2D
        ! compute eigenvalues at k
        CALL SetEigenValues2d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,2), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k),&
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4))
        ! compute Riemann invariants
        CALL Prim2Riemann2d(this%gamma, &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,2), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4))

      CASE(3) ! 3D
        ! compute eigenvalues at k
        CALL SetEigenValues3d(this%gamma, &
          p%density%data3d( Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,3), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k),&
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))

        ! compute Riemann invariants
        CALL Prim2Riemann3d(this%gamma, &
          p%density%data3d( Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,3), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,2), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))

      END SELECT
    END SELECT
  END SUBROUTINE CalculatePrim2RiemannZ


  !> Convert Riemann invariants to primitives for farfield boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateRiemann2PrimX(this,Mesh,i,Rinv,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: i
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM) &
                                       :: Rinv
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL Riemann2Prim1d(this%gamma, &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CASE(2) ! 2D
        CALL Riemann2Prim2d(this%gamma, &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CASE(3) ! 3D
        CALL Riemann2Prim3d(this%gamma, &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5), &
          p%density%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          p%pressure%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateRiemann2PrimX


  !> Convert Riemann invariants to primitives for farfield boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateRiemann2PrimY(this,Mesh,j,Rinv,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: j
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,this%VNUM) &
                                       :: Rinv
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL Riemann2Prim1d(this%gamma, &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX))
      CASE(2) ! 2D
        CALL Riemann2Prim2d(this%gamma, &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,2), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX))
      CASE(3) ! 3D
        CALL Riemann2Prim3d(this%gamma, &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,5), &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,2), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,3), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateRiemann2PrimY


  !> Convert Riemann invariants to primitives for farfield boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateRiemann2PrimZ(this,Mesh,k,Rinv,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: k
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,this%VNUM) &
                                       :: Rinv
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D
        CALL Riemann2Prim1d(this%gamma, &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k))
      CASE(2) ! 2D
        CALL Riemann2Prim2d(this%gamma, &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,2), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k))
      CASE(3) ! 3D
        CALL Riemann2Prim3d(this%gamma, &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
          p%density%data3d( Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,3), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,1), &
          p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,2), &
          p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k))
      END SELECT
    END SELECT
  END SUBROUTINE CalculateRiemann2PrimZ

  !> Calculates geometrical sources
  PURE SUBROUTINE GeometricalSources(this,Mesh,pvar,cvar,sterm)
    USE geometry_generic_mod
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    ! compute geometrical source only for non-cartesian mesh
    SELECT TYPE(mgeo => Mesh%geometry)
    TYPE IS(geometry_cartesian)
      ! do nothing
    CLASS DEFAULT
      ! curvilinear geometry
    SELECT TYPE(p => pvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_euler)
          SELECT TYPE(s => sterm)
          TYPE IS(statevector_euler)
            ! no source terms
            s%density%data1d(:) = 0.0
            s%energy%data1d(:) = 0.0
            SELECT CASE(Mesh%VECTOR_COMPONENTS)
            CASE(VECTOR_X) ! 1D momentum in x-direction
              ! vy = vz = my = mz = 0
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),0.0,0.0, &
                  p%pressure%data1d(:), &
                  0.0,0.0)
            CASE(VECTOR_Y) ! 1D momentum in y-direction
              ! vx = vz = mx = mz = 0
              s%momentum%data2d(:,1) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,p%velocity%data2d(:,1),0.0, &
                  p%pressure%data1d(:), &
                  0.0,0.0)
            CASE(VECTOR_Z) ! 1D momentum in z-direction
              ! vx = vy = mx = my = 0
              s%momentum%data2d(:,1) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,0.0,p%velocity%data2d(:,1), &
                  p%pressure%data1d(:), &
                  0.0,0.0)
            CASE(IOR(VECTOR_X,VECTOR_Y)) ! 2D momentum in x-y-plane
              ! vz = mz = 0
              ! x-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2),0.0, &
                  p%pressure%data1d(:), &
                  c%momentum%data2d(:,2),0.0)
              ! y-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2),0.0, &
                  p%pressure%data1d(:), &
                  c%momentum%data2d(:,1),0.0)
            CASE(IOR(VECTOR_X,VECTOR_Z)) ! 2D momentum in x-z-plane
              ! vy = my = 0
              ! x-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),0.0,p%velocity%data2d(:,2), &
                  p%pressure%data1d(:), &
                  0.0,c%momentum%data2d(:,2))
              ! z-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),0.0,p%velocity%data2d(:,2), &
                  p%pressure%data1d(:), &
                  c%momentum%data2d(:,1),0.0)
            CASE(IOR(VECTOR_Y,VECTOR_Z)) ! 2D momentum in y-z-plane
              ! vx = mx = 0
              ! y-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%pressure%data1d(:), &
                  0.0,c%momentum%data2d(:,2))
              ! z-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%pressure%data1d(:), &
                  0.0,c%momentum%data2d(:,1))
            CASE(IOR(IOR(VECTOR_X,VECTOR_Y),VECTOR_Z)) ! 3D momentum
              ! x-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%velocity%data2d(:,3),p%pressure%data1d(:), &
                  c%momentum%data2d(:,2),c%momentum%data2d(:,3))
              ! y-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%velocity%data2d(:,3),p%pressure%data1d(:), &
                  c%momentum%data2d(:,1),c%momentum%data2d(:,3))
              ! z-momentum
              s%momentum%data2d(:,3) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%velocity%data2d(:,3),p%pressure%data1d(:), &
                  c%momentum%data2d(:,1),c%momentum%data2d(:,2))
            CASE DEFAULT
              ! return NaN
              s%momentum%data1d(:) = 0.0 !NAN_DEFAULT_REAL
            END SELECT
          END SELECT
        END SELECT
      END SELECT
      ! reset ghost cell data
      IF (Mesh%INUM.GT.1) THEN
        sterm%data4d(Mesh%IGMIN:Mesh%IMIN+Mesh%IM1,:,:,:) = 0.0
        sterm%data4d(Mesh%IMAX+Mesh%IP1:Mesh%IGMAX,:,:,:) = 0.0
      END IF
      IF (Mesh%JNUM.GT.1) THEN
        sterm%data4d(:,Mesh%JGMIN:Mesh%JMIN+Mesh%JM1,:,:) = 0.0
        sterm%data4d(:,Mesh%JMAX+Mesh%JP1:Mesh%JGMAX,:,:) = 0.0
      END IF
      IF (Mesh%KNUM.GT.1) THEN
        ! collapse the first 2 dimensions
        sterm%data3d(:,Mesh%KGMIN:Mesh%KMIN+Mesh%KM1,:) = 0.0
        sterm%data3d(:,Mesh%KMAX+Mesh%KP1:Mesh%KGMAX,:) = 0.0
      END IF
    END SELECT
  END SUBROUTINE GeometricalSources

  !> compute momentum and energy sources given an external force
  SUBROUTINE ExternalSources(this,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_euler),   INTENT(IN)  :: this
    CLASS(marray_base),     INTENT(IN)  :: accel
    CLASS(marray_compound), INTENT(IN)  :: pvar,cvar
    CLASS(marray_compound), INTENT(INOUT) :: sterm
    !------------------------------------------------------------------------!
    INTEGER :: m
    !------------------------------------------------------------------------!
    CALL this%physics_eulerisotherm%ExternalSources(accel,pvar,cvar,sterm)
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(s => sterm)
        TYPE IS(statevector_euler)
          SELECT CASE(this%VDIM)
          CASE(1)
            DO m=1,SIZE(s%energy%data1d)
              s%energy%data1d(m) = c%momentum%data2d(m,1) * accel%data2d(m,1)
            END DO
          CASE(2)
            DO m=1,SIZE(s%energy%data1d)
              s%energy%data1d(m) = c%momentum%data2d(m,1) * accel%data2d(m,1) &
                                 + c%momentum%data2d(m,2) * accel%data2d(m,2)
            END DO
          CASE(3)
            DO m=1,SIZE(s%energy%data1d)
              s%energy%data1d(m) = c%momentum%data2d(m,1) * accel%data2d(m,1) &
                                 + c%momentum%data2d(m,2) * accel%data2d(m,2) &
                                 + c%momentum%data2d(m,3) * accel%data2d(m,3)
            END DO
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE ExternalSources

  SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(marray_compound), INTENT(INOUT) :: pvar,sterm
    REAL,                   INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                                           :: btxx,btxy,btxz,btyy,btyz,btzz
   !------------------------------------------------------------------------!
   CALL this%physics_eulerisotherm%ViscositySources(Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
   SELECT TYPE(p => pvar)
   CLASS IS(statevector_euler)
     SELECT TYPE(s => sterm)
     CLASS IS(statevector_euler)
        SELECT CASE(this%VDIM)
      !   CASE(1)
      !     this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:)
      !
      !     CALL Mesh%Divergence(this%tmp(:,:,:),sterm(:,:,:,this%ENERGY))
        CASE(2)
          !compute scalar product of v and tau (x-component)
          this%tmp(:,:,:) = p%velocity%data4d(:,:,:,1)*btxx(:,:,:) &
                          + p%velocity%data4d(:,:,:,2)*btxy(:,:,:)

          !compute scalar product of v and tau (y-component)
          this%tmp1(:,:,:) = p%velocity%data4d(:,:,:,1)*btxy(:,:,:) &
                           + p%velocity%data4d(:,:,:,2)*btyy(:,:,:)

          ! compute vector divergence of scalar product v and tau
          CALL Mesh%Divergence(this%tmp(:,:,:),this%tmp1(:,:,:), &
                  s%energy%data3d(:,:,:))
        CASE(3)
          !compute scalar product of v and tau (x-component)
          this%tmp(:,:,:) = p%velocity%data4d(:,:,:,1)*btxx(:,:,:) &
                          + p%velocity%data4d(:,:,:,2)*btxy(:,:,:) &
                          + p%velocity%data4d(:,:,:,3)*btxz(:,:,:)

          !compute scalar product of v and tau (y-component)
          this%tmp1(:,:,:) = p%velocity%data4d(:,:,:,1)*btxy(:,:,:) &
                           + p%velocity%data4d(:,:,:,2)*btyy(:,:,:) &
                           + p%velocity%data4d(:,:,:,3)*btyz(:,:,:)

          !compute scalar product of v and tau (z-component)
          this%tmp2(:,:,:) = p%velocity%data4d(:,:,:,1)*btxz(:,:,:) &
                           + p%velocity%data4d(:,:,:,2)*btyz(:,:,:) &
                           + p%velocity%data4d(:,:,:,3)*btzz(:,:,:)
          ! compute vector divergence of scalar product v and tau
          CALL Mesh%Divergence(this%tmp(:,:,:),this%tmp1(:,:,:),this%tmp2(:,:,:), &
                  s%energy%data3d(:,:,:))
        CASE DEFAULT
          ! return NaN
          s%data1d(:) = NAN_DEFAULT_REAL
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE ViscositySources

  !> Adds a background velocity field for fargo routines
  !!
  !! Calculates
  !! \f{eqnarray*}{
  !!    E   &=& E' + m_x' w + \frac{1}{2}\varrho w^2 \\
  !!    v_x &=& v_x' +  w \\
  !!    m_x &=& m_x' +  \varrho w,
  !! \f}
  !! with \f$ E, v_y, m_y \f$ the total energy, velocity and momentum. The
  !! \f$ ' \f$ denotes the residual part. \f$ w \f$ is the velocity shift.
  PURE SUBROUTINE AddBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                            INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        IF (c%fargo_transformation_applied) THEN
          IF (.NOT.p%fargo_transformation_applied) RETURN ! should not happen
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                ! ATTENTION: don't change the order; on the RHS of the first
                !            assignment there must be the old momentum
                c%energy%data3d(i,j,k) = c%energy%data3d(i,j,k) &
                    + w(j,k)*(c%momentum%data4d(i,j,k,1) &
                    + 0.5*c%density%data3d(i,j,k)*w(j,k))
                p%velocity%data4d(i,j,k,1) = p%velocity%data4d(i,j,k,1) + w(j,k)
                c%momentum%data4d(i,j,k,1) = c%momentum%data4d(i,j,k,1) &
                    + c%density%data3d(i,j,k)*w(j,k)
              END DO
            END DO
          END DO
          c%fargo_transformation_applied = .FALSE.
          p%fargo_transformation_applied = .FALSE.
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE AddBackgroundVelocityX

  !> Adds a background velocity field for fargo routines
  !!
  !! Calculates
  !! \f{eqnarray*}{
  !!    E   &=& E' + m_y' w + \frac{1}{2}\varrho w^2 \\
  !!    v_y &=& v_y' +  w \\
  !!    m_y &=& m_y' +  \varrho w,
  !! \f}
  !! with \f$ E, v_y, m_y \f$ the total energy, velocity and momentum. The
  !! \f$ ' \f$ denotes the residual part. \f$ w \f$ is the velocity shift.
  PURE SUBROUTINE AddBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                            INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k,v_idx
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        IF (c%fargo_transformation_applied) THEN
          IF (.NOT.p%fargo_transformation_applied) RETURN ! should not happen
          ! check if x-component of velocity vector is available
          IF (BTEST(Mesh%VECTOR_COMPONENTS,0)) THEN
            ! y-component is the 2nd component
            v_idx = 2
          ELSE
            ! no x-component -> y-component is the first component
            v_idx = 1
          END IF
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                ! ATTENTION: don't change the order; on the RHS of the first
                !            assignment there must be the old momentum
                c%energy%data3d(i,j,k) = c%energy%data3d(i,j,k) &
                    + w(i,k)*(c%momentum%data4d(i,j,k,v_idx) &
                    + 0.5*c%density%data3d(i,j,k)*w(i,k))
                p%velocity%data4d(i,j,k,v_idx) = p%velocity%data4d(i,j,k,v_idx) + w(i,k)
                c%momentum%data4d(i,j,k,v_idx) = c%momentum%data4d(i,j,k,v_idx) &
                    + c%density%data3d(i,j,k)*w(i,k)
              END DO
            END DO
          END DO
          c%fargo_transformation_applied = .FALSE.
          p%fargo_transformation_applied = .FALSE.
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE AddBackgroundVelocityY

  !> Adds a background velocity field for fargo routines
  !!
  !! Calculates
  !! \f{eqnarray*}{
  !!    E   &=& E' + m_z' w + \frac{1}{2}\varrho w^2 \\
  !!    v_z &=& v_z' +  w \\
  !!    m_z &=& m_z' +  \varrho w,
  !! \f}
  !! with \f$ E, v_z, m_z \f$ the total energy, velocity and momentum. The
  !! \f$ ' \f$ denotes the residual part. \f$ w \f$ is the velocity shift.
  PURE SUBROUTINE AddBackgroundVelocityZ(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                            INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        IF (c%fargo_transformation_applied) THEN
          IF (.NOT.p%fargo_transformation_applied) RETURN ! should not happen
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                ! ATTENTION: don't change the order; on the RHS of the first
                !            assignment there must be the old momentum
                c%energy%data3d(i,j,k) = c%energy%data3d(i,j,k) &
                      + w(i,j)*(c%momentum%data4d(i,j,k,this%VDIM) &
                      + 0.5*c%density%data3d(i,j,k)*w(i,j))
                p%velocity%data4d(i,j,k,this%VDIM) = p%velocity%data4d(i,j,k,this%VDIM) + w(i,j)
                c%momentum%data4d(i,j,k,this%VDIM) = c%momentum%data4d(i,j,k,this%VDIM) &
                      + c%density%data3d(i,j,k)*w(i,j)
              END DO
            END DO
          END DO
          c%fargo_transformation_applied = .FALSE.
          p%fargo_transformation_applied = .FALSE.
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE AddBackgroundVelocityZ

  !> Substracts a background velocity field for fargo routines
  !!
  !! Calculates
  !! \f{eqnarray*}{
  !!    E'   &=& E - m_x w + \frac{1}{2}\varrho w^2 \\
  !!    v_x' &=& v_x -  w \\
  !!    m_x' &=& m_x -  \varrho w,
  !! \f}
  !! with \f$ E, v_x, m_x \f$ the total energy, velocity and momentum. The
  !! \f$ ' \f$ denotes the residual part. \f$ w \f$ is the velocity shift.
  !! ATTENTION: the "+" before the 3rd term in the energy transformation is
  !!            correct!
  PURE SUBROUTINE SubtractBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                            INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        IF (.NOT.c%fargo_transformation_applied) THEN
          IF (p%fargo_transformation_applied) RETURN ! should not happen
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                ! ATTENTION: don't change the order; on the RHS of the first
                !            assignment there must be the old momentum
                c%energy%data3d(i,j,k) = c%energy%data3d(i,j,k) &
                    - w(j,k)*(c%momentum%data4d(i,j,k,1) &
                    - 0.5*c%density%data3d(i,j,k)*w(j,k))
                p%velocity%data4d(i,j,k,1) = p%velocity%data4d(i,j,k,1) - w(j,k)
                c%momentum%data4d(i,j,k,1) = c%momentum%data4d(i,j,k,1) &
                    - c%density%data3d(i,j,k)*w(j,k)
              END DO
            END DO
          END DO
          c%fargo_transformation_applied = .TRUE.
          p%fargo_transformation_applied = .TRUE.
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE SubtractBackgroundVelocityX

  !> Substracts a background velocity field for fargo routines
  !!
  !! Calculates
  !! \f{eqnarray*}{
  !!    E'   &=& E - m_y w + \frac{1}{2}\varrho w^2 \\
  !!    v_y' &=& v_y -  w \\
  !!    m_y' &=& m_y -  \varrho w,
  !! \f}
  !! with \f$ E, v_y, m_y \f$ the total energy, velocity and momentum. The
  !! \f$ ' \f$ denotes the residual part. \f$ w \f$ is the velocity shift.
  !! ATTENTION: the "+" before the 3rd term in the energy transformation is
  !!            correct!
  PURE SUBROUTINE SubtractBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                            INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k,v_idx
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        IF (.NOT.c%fargo_transformation_applied) THEN
          IF (p%fargo_transformation_applied) RETURN ! should not happen
          ! check if x-component of velocity vector is available
          IF (BTEST(Mesh%VECTOR_COMPONENTS,0)) THEN
            ! y-component is the 2nd component
            v_idx = 2
          ELSE
            ! no x-component -> y-component is the first component
            v_idx = 1
          END IF
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                ! ATTENTION: don't change the order; on the RHS of the first
                !            assignment there must be the old momentum
                c%energy%data3d(i,j,k) = c%energy%data3d(i,j,k) &
                    - w(i,k)*(c%momentum%data4d(i,j,k,v_idx) &
                    - 0.5*c%density%data3d(i,j,k)*w(i,k))
                p%velocity%data4d(i,j,k,v_idx) = p%velocity%data4d(i,j,k,v_idx) - w(i,k)
                c%momentum%data4d(i,j,k,v_idx) = c%momentum%data4d(i,j,k,v_idx) &
                    - c%density%data3d(i,j,k)*w(i,k)
              END DO
            END DO
          END DO
          c%fargo_transformation_applied = .TRUE.
          p%fargo_transformation_applied = .TRUE.
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE SubtractBackgroundVelocityY

  !> Substracts a background velocity field for fargo routines
  !!
  !! Calculates
  !! \f{eqnarray*}{
  !!    E'   &=& E - m_z w + \frac{1}{2}\varrho w^2 \\
  !!    v_z' &=& v_z -  w \\
  !!    m_z' &=& m_z -  \varrho w,
  !! \f}
  !! with \f$ E, v_z, m_z \f$ the total energy, velocity and momentum. The
  !! \f$ ' \f$ denotes the residual part. \f$ w \f$ is the velocity shift.
  !! ATTENTION: the "+" before the 3rd term in the energy transformation is
  !!            correct!
  PURE SUBROUTINE SubtractBackgroundVelocityZ(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                            INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        IF (.NOT.c%fargo_transformation_applied) THEN
          IF (p%fargo_transformation_applied) RETURN ! should not happen
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                ! ATTENTION: don't change the order; on the RHS of the first
                !            assignment there must be the old momentum
                c%energy%data3d(i,j,k) = c%energy%data3d(i,j,k) &
                      - w(i,j)*(c%momentum%data4d(i,j,k,this%VDIM) &
                      - 0.5*c%density%data3d(i,j,k)*w(i,j))
                p%velocity%data4d(i,j,k,this%VDIM) = p%velocity%data4d(i,j,k,this%VDIM) - w(i,j)
                c%momentum%data4d(i,j,k,this%VDIM) = c%momentum%data4d(i,j,k,this%VDIM) &
                      - c%density%data3d(i,j,k)*w(i,j)
              END DO
            END DO
          END DO
          c%fargo_transformation_applied = .TRUE.
          p%fargo_transformation_applied = .TRUE.
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE SubtractBackgroundVelocityZ


  !> \public sources terms for fargo advection along x-direction
  !!
  !! If the background velocity \f$\vec{w}=w\,\hat{e}_\xi\f$ with
  !! \f$w\f$ independent of \f$\xi\f$ and \f$t\f$ is subtracted from
  !! the overall velocity of the flow, an additional source term occurs
  !! in the \f$\xi\f$-momentum equation:
  !! \f[
  !!     S_\mathrm{Fargo} = -\varrho \vec{v} \cdot \nabla \vec{w}
  !! \f]
  !! where \f$\vec{v}\f$ is the real velocity (including the background
  !! velocity \f$\vec{w}\f$ ).
  PURE SUBROUTINE AddFargoSourcesX(this,Mesh,w,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base), INTENT(IN)             :: Mesh
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), INTENT(IN) :: w
    CLASS(marray_compound), INTENT(INOUT)    :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    CLASS IS(statevector_euler)
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(s => sterm)
        CLASS IS(statevector_euler)
          ! ATTENTION: fargo sources are added to the given sterm
          SELECT CASE(Mesh%VECTOR_COMPONENTS)
          CASE(IOR(VECTOR_X,VECTOR_Y))
            ! 2D (x,y) momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,2) * 0.5 * (w(j+1,k)-w(j-1,k)) &
                    / Mesh%dly%data3d(i,j,k)
                  ! x-momentum source
                  s%momentum%data4d(i,j,k,1) = s%momentum%data4d(i,j,k,1) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceX( &
                    Mesh%cxyx%data4d(i,j,k,2),Mesh%cxzx%data4d(i,j,k,2), &
                    Mesh%cyxy%data4d(i,j,k,2),Mesh%czxz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2),0.0, &
                    0.0,c%momentum%data4d(i,j,k,2),0.0)
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                    - p%velocity%data4d(i,j,k,1)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE(IOR(VECTOR_X,VECTOR_Z))
            ! 2D (x,z) momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,2) * 0.5 * (w(j,k+1)-w(j,k-1)) &
                    / Mesh%dlz%data3d(i,j,k)
                  ! x-momentum source
                  s%momentum%data4d(i,j,k,1) = s%momentum%data4d(i,j,k,1) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceX( &
                    Mesh%cxyx%data4d(i,j,k,2),Mesh%cxzx%data4d(i,j,k,2), &
                    Mesh%cyxy%data4d(i,j,k,2),Mesh%czxz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),0.0,p%velocity%data4d(i,j,k,2), &
                    0.0,0.0,c%momentum%data4d(i,j,k,2))
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                    - p%velocity%data4d(i,j,k,1)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE(IOR(IOR(VECTOR_X,VECTOR_Y),VECTOR_Z))
            ! 3D momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,2) * 0.5 * (w(j+1,k)-w(j-1,k)) / Mesh%dly%data3d(i,j,k) &
                    + c%momentum%data4d(i,j,k,3) * 0.5 * (w(j,k+1)-w(j,k-1)) / Mesh%dlz%data3d(i,j,k)
                  ! x-momentum source
                  s%momentum%data4d(i,j,k,1) = s%momentum%data4d(i,j,k,1) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceX( &
                    Mesh%cxyx%data4d(i,j,k,2),Mesh%cxzx%data4d(i,j,k,2), &
                    Mesh%cyxy%data4d(i,j,k,2),Mesh%czxz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,3), &
                    0.0,c%momentum%data4d(i,j,k,2),c%momentum%data4d(i,j,k,3))
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                    - p%velocity%data4d(i,j,k,1)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE DEFAULT
            ! do nothing in any other case
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE AddFargoSourcesX

  !> \public sources terms for fargo advection along y-direction
  !!
  !! If the background velocity \f$\vec{w}=w\,\hat{e}_\eta\f$ with
  !! \f$w\f$ independent of \f$\eta\f$ and \f$t\f$ is subtracted from
  !! the overall velocity of the flow, an additional source term occurs
  !! in the \f$\eta\f$-momentum equation:
  !! \f[
  !!     S_\mathrm{Fargo} = -\varrho \vec{v} \cdot \nabla \vec{w}
  !! \f]
  !! where \f$\vec{v}\f$ is the real velocity (including the background
  !! velocity \f$\vec{w}\f$ ).
  PURE SUBROUTINE AddFargoSourcesY(this,Mesh,w,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base), INTENT(IN)             :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), INTENT(IN) :: w
    CLASS(marray_compound), INTENT(INOUT)    :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    CLASS IS(statevector_euler)
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(s => sterm)
        CLASS IS(statevector_euler)
          SELECT CASE(Mesh%VECTOR_COMPONENTS)
          CASE(IOR(VECTOR_X,VECTOR_Y))
            ! 2D (x,y) momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,1) * 0.5 * (w(i+1,k)-w(i-1,k)) &
                      / Mesh%dlx%data3d(i,j,k)
                  ! y-momentum source
                  s%momentum%data4d(i,j,k,2) = s%momentum%data4d(i,j,k,2) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceY( &
                    Mesh%cxyx%data4d(i,j,k,2),Mesh%cyxy%data4d(i,j,k,2), &
                    Mesh%cyzy%data4d(i,j,k,2),Mesh%czyz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2),0.0, &
                    0.0,c%momentum%data4d(i,j,k,1),0.0)
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                    - p%velocity%data4d(i,j,k,2)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE(IOR(VECTOR_Y,VECTOR_Z))
            ! 2D (y,z) momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,2) * 0.5 * (w(i,k+1)-w(i,k-1)) &
                    / Mesh%dlz%data3d(i,j,k)
                  ! y-momentum source
                  s%momentum%data4d(i,j,k,1) = s%momentum%data4d(i,j,k,1) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceY( &
                    Mesh%cxyx%data4d(i,j,k,2),Mesh%cyxy%data4d(i,j,k,2), &
                    Mesh%cyzy%data4d(i,j,k,2),Mesh%czyz%data4d(i,j,k,2), &
                    0.0,p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2), &
                    0.0,0.0,c%momentum%data4d(i,j,k,2))
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                    - p%velocity%data4d(i,j,k,1)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE(IOR(IOR(VECTOR_X,VECTOR_Y),VECTOR_Z))
            ! 3D momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,1) * 0.5 * (w(i+1,k)-w(i-1,k)) / Mesh%dlx%data3d(i,j,k) &
                    + c%momentum%data4d(i,j,k,3) * 0.5 * (w(j,k+1)-w(j,k-1)) / Mesh%dlz%data3d(i,j,k)
                  ! y-momentum source
                  s%momentum%data4d(i,j,k,2) = s%momentum%data4d(i,j,k,2) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceY( &
                    Mesh%cxyx%data4d(i,j,k,2),Mesh%cyxy%data4d(i,j,k,2), &
                    Mesh%cyzy%data4d(i,j,k,2),Mesh%czyz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,3),0.0, &
                    c%momentum%data4d(i,j,k,1),c%momentum%data4d(i,j,k,3))
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                      - p%velocity%data4d(i,j,k,2)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE DEFAULT
            ! do nothing in any other case
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE AddFargoSourcesY

  !> \public sources terms for fargo advection along z-direction
  !!
  !! If the background velocity \f$\vec{w}=w\,\hat{e}_phi\f$ with
  !! \f$w\f$ independent of \f$\phi\f$ and \f$t\f$ is subtracted from
  !! the overall velocity of the flow, an additional source term occurs
  !! in the \f$\phi\f$-momentum equation:
  !! \f[
  !!     S_\mathrm{Fargo} = -\varrho \vec{v} \cdot \nabla \vec{w}
  !! \f]
  !! where \f$\vec{v}\f$ is the real velocity (including the background
  !! velocity \f$\vec{w}\f$ ).
  PURE SUBROUTINE AddFargoSourcesZ(this,Mesh,w,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT)         :: this
    CLASS(mesh_base), INTENT(IN)             :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), INTENT(IN) :: w
    CLASS(marray_compound), INTENT(INOUT)    :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    CLASS IS(statevector_euler)
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(s => sterm)
        CLASS IS(statevector_euler)
          SELECT CASE(Mesh%VECTOR_COMPONENTS)
          CASE(IOR(VECTOR_X,VECTOR_Z))
            ! 2D (x,z) momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,1) * 0.5 * (w(i+1,j)-w(i-1,j)) / Mesh%dlx%data3d(i,j,k)
                  ! z-momentum source
                  s%momentum%data4d(i,j,k,2) = s%momentum%data4d(i,j,k,2) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceZ( &
                    Mesh%cxzx%data4d(i,j,k,2),Mesh%cyzy%data4d(i,j,k,2), &
                    Mesh%czxz%data4d(i,j,k,2),Mesh%czyz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),0.0,p%velocity%data4d(i,j,k,2), &
                    0.0,c%momentum%data4d(i,j,k,1),0.0)
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                      - p%velocity%data4d(i,j,k,2)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE(IOR(VECTOR_Y,VECTOR_Z))
            ! 2D (y,z) momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,1) * 0.5 * (w(i,j+1)-w(i,j-1)) / Mesh%dly%data3d(i,j,k)
                  ! z-momentum source
                  s%momentum%data4d(i,j,k,2) = s%momentum%data4d(i,j,k,2) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceZ( &
                    Mesh%cxzx%data4d(i,j,k,2),Mesh%cyzy%data4d(i,j,k,2), &
                    Mesh%czxz%data4d(i,j,k,2),Mesh%czyz%data4d(i,j,k,2), &
                    0.0,p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2), &
                    0.0,0.0,c%momentum%data4d(i,j,k,1))
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                      - p%velocity%data4d(i,j,k,2)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE(IOR(IOR(VECTOR_X,VECTOR_Y),VECTOR_Z))
            ! 3D momentum vector
            DO k=Mesh%KMIN,Mesh%KMAX
              DO j=Mesh%JMIN,Mesh%JMAX
                DO i=Mesh%IMIN,Mesh%IMAX
                  this%tmp(i,j,k) = c%momentum%data4d(i,j,k,1) * 0.5 * (w(i+1,j)-w(i-1,j)) / Mesh%dlx%data3d(i,j,k) &
                    + c%momentum%data4d(i,j,k,2) * 0.5 * (w(j,j+1)-w(j,j-1)) / Mesh%dly%data3d(i,j,k)
                  ! z-momentum source
                  s%momentum%data4d(i,j,k,3) = s%momentum%data4d(i,j,k,3) - this%tmp(i,j,k)
                  ! add geometrical terms
                  this%tmp(i,j,k) = this%tmp(i,j,k) + GetGeometricalSourceZ( &
                    Mesh%cxzx%data4d(i,j,k,2),Mesh%cyzy%data4d(i,j,k,2), &
                    Mesh%czxz%data4d(i,j,k,2),Mesh%czyz%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,1),p%velocity%data4d(i,j,k,2), &
                    p%velocity%data4d(i,j,k,3), &
                    0.0,c%momentum%data4d(i,j,k,1),c%momentum%data4d(i,j,k,2))
                  ! energy source
                  s%energy%data3d(i,j,k) = s%energy%data3d(i,j,k) &
                      - p%velocity%data4d(i,j,k,3)*this%tmp(i,j,k)
                END DO
              END DO
            END DO
          CASE DEFAULT
            ! do nothing in any other case
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE AddFargoSourcesZ

  PURE SUBROUTINE UpdateSoundSpeed(this,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(statevector_euler),INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i
    !------------------------------------------------------------------------!
    IF (pvar%density%RANK.EQ.0) THEN
      DO CONCURRENT (i=1:SIZE(this%bccsound%data1d))
        this%bccsound%data1d(i) = GetSoundSpeed(this%gamma, &
            pvar%density%data1d(i),pvar%pressure%data1d(i))
      END DO
    ELSE
      DO CONCURRENT (i=1:SIZE(this%fcsound%data1d))
        this%fcsound%data1d(i) = GetSoundSpeed(this%gamma, &
            pvar%density%data1d(i),pvar%pressure%data1d(i))
      END DO
    END IF
  END SUBROUTINE UpdateSoundSpeed

  !> \public Destructor of the physics_euler class
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%physics_eulerisotherm%Finalize()
  END SUBROUTINE Finalize

!----------------------------------------------------------------------------!
!> \par methods of class statevector_euler

  !> \public Constructor of statevector_euler
  !!
  !! This is not a class member itself, instead its an ordinary
  !! module procedure. The function name is overloaded with
  !! the class name.
  FUNCTION CreateStateVector(Physics,flavour,num) RESULT(new_sv)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: Physics
    INTEGER, OPTIONAL, INTENT(IN) :: flavour,num
    TYPE(statevector_euler) :: new_sv
    !-------------------------------------------------------------------!
    LOGICAL :: success = .FALSE.
    INTEGER :: err
    !-------------------------------------------------------------------!
    ! call inherited function
    new_sv = statevector_eulerisotherm(Physics,flavour,num)
    ! add entries specific for euler physics
    ALLOCATE(new_sv%pressure,STAT=err)
    IF (err.EQ.0) THEN
      IF (PRESENT(num)) THEN
        new_sv%pressure = marray_base(num)      ! num scalars
      ELSE
        new_sv%pressure = marray_base()         ! one scalar
      END IF
      success = new_sv%AppendMArray(new_sv%pressure)
      SELECT CASE(new_sv%flavour)
      CASE(PRIMITIVE)
        new_sv%energy => null()
      CASE(CONSERVATIVE)
        new_sv%energy => new_sv%pressure
        new_sv%pressure => null()
      CASE DEFAULT
#ifdef DEBUG
        PRINT *,"ERROR in physics_euler::CreateStateVector: unknown flavour"
#endif
      END SELECT
    END IF
    IF (.NOT.success) &
      CALL Physics%Error("physics_euler::CreateStateVector", &
                         "state vector initialization failed")
  END FUNCTION CreateStateVector

  !> assigns one state vector to another state vector
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(statevector_euler),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in physics_euler::AssignMArray_0: assigning 2 state vectors"
#endif
    CALL this%statevector_eulerisotherm%AssignMArray_0(ma)
    IF (SIZE(this%data1d).LE.0) RETURN ! empty compound
    SELECT TYPE(src => ma)
    CLASS IS(statevector_euler)
#if DEBUG > 2
      PRINT *,"DEBUG INFO in physics_euler::AssignMArray_0: restoring pressure/energy pointers"
#endif
      SELECT CASE(this%flavour)
      CASE(PRIMITIVE)
        ! pressure is the third item
        this%pressure => this%GetItem(this%NextItem(this%NextItem(this%FirstItem())))
        this%energy => null()
      CASE(CONSERVATIVE)
        ! energy is the third item
        this%energy => this%GetItem(this%NextItem(this%NextItem(this%FirstItem())))
        this%pressure => null()
      CASE DEFAULT
        ! error, this should not happen
      END SELECT
    CLASS DEFAULT
      ! do nothing: ma may not be of type euler, i.e. during initialization (see CreateStatevector)
    END SELECT
  END SUBROUTINE AssignMArray_0

    !> actual destructor of the statevector_eulerisotherm type
  SUBROUTINE Finalize_statevector(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(statevector_euler),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
#if DEBUG > 2
    PRINT *,"DEBUG INFO in physics_euler::Finalize: cleanup statevector"
#endif
    NULLIFY(this%pressure,this%energy)
  END SUBROUTINE Finalize_statevector


!----------------------------------------------------------------------------!
!> \par elemental non-class subroutines / functions


  ELEMENTAL FUNCTION GetSoundSpeed(gamma,density,pressure) RESULT(cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,density,pressure
    REAL              :: cs
    !------------------------------------------------------------------------!
    cs = SQRT(MAX(2.0*TINY(cs),gamma*pressure/density))
  END FUNCTION GetSoundSpeed

  !> \todo NOT VERIFIED
  !! only for advanced wavespeeds
  ELEMENTAL SUBROUTINE SetRoeAverages(gamma,rhoL,rhoR,ul,uR,vL,vR,pL,pR,eL,eR,u,cs)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rhoL,rhoR,uL,uR,vL,vR,pL,pR,eL,eR
    REAL, INTENT(OUT) :: u,cs
    !------------------------------------------------------------------------!
    REAL              :: sqrtrhoL,sqrtrhoR,invsqrtrho,v,hL,hR,h
    !------------------------------------------------------------------------!
    sqrtrhoL = SQRT(rhoL)
    sqrtrhoR = SQRT(rhoR)
    ! density
    invsqrtrho = 1./ (sqrtrhoL + sqrtrhoR)
    ! velocities
    u = (sqrtrhoL*uL + sqrtrhoR*uR) * invsqrtrho
    v = (sqrtrhoL*vL + sqrtrhoR*vR) * invsqrtrho
    ! enthalpy
    hL = (eL + pL) / rhoL
    hR = (eR + pR) / rhoR
    h  = (sqrtrhoL * hL + sqrtrhoR * hR) * invsqrtrho
    ! sound speed
    cs = SQRT((gamma-1.)*(h-0.5*(u**2+v**2)))
  END SUBROUTINE SetRoeAverages

  !> set minimal and maximal wave speeds
  ELEMENTAL SUBROUTINE SetWaveSpeeds(cs,v,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: minwav,maxwav
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    minwav = MIN(0.,v-cs)
    maxwav = MAX(0.,v+cs)
  END SUBROUTINE SetWaveSpeeds

  !> set all eigenvalues for 1D transport (used in absorbing boundary conditions)
  ELEMENTAL SUBROUTINE SetEigenValues1d(gamma,rho,v,P,l1,l2,l3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    ! adiabatic sound speed
    cs = GetSoundSpeed(gamma,rho,P)
    ! call subroutine for isothermal case with the adiabatic sound speed
    l1 = v - cs
    l2 = v
    l3 = v + cs
  END SUBROUTINE SetEigenValues1d

  !> set all eigenvalues for 2D transport (used in absorbing boundary conditions)
  ELEMENTAL SUBROUTINE SetEigenValues2d(gamma,rho,v,P,l1,l2,l3,l4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4
    !------------------------------------------------------------------------!
    CALL SetEigenValues1d(gamma,rho,v,P,l1,l2,l4)
    l3 = v
  END SUBROUTINE SetEigenValues2d

  !> set all eigenvalues for 3D transport (used in absorbing boundary conditions)
  ELEMENTAL SUBROUTINE SetEigenValues3d(gamma,rho,v,P,l1,l2,l3,l4,l5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4,l5
    !------------------------------------------------------------------------!
    CALL SetEigenValues2d(gamma,rho,v,P,l1,l2,l3,l5)
    l4 = v
  END SUBROUTINE SetEigenValues3d

  ! \todo NOT VERIFIED
  !! only for HLLC fluxes
  ELEMENTAL SUBROUTINE SetIntermediateState(rhoL,rhoR,uL,uR,vl,vR,&
       pL,pR,eL,eR,amin,amax,rho,mu,mv,e,a)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rhoL,rhoR,uL,uR,vl,vR,pL,pR,eL,eR,amin,amax
    REAL, INTENT(OUT) :: rho,mu,mv,e,a
    !------------------------------------------------------------------------!
    REAL              :: qL,qR
    !------------------------------------------------------------------------!
    qL = rhoL * (uL-amin)
    qR = rhoR * (uR-amax)
    ! wave speed in intermediate region
    a = (pR - pL + qR*uR - qL*uL) / (qR - qL)
    ! set left or right state vector depending on wave speed "a"
    IF (a.GT.0.0) THEN
       ! left state
       rho = qL / (a - amin)
       mu  = rho * a
       mv  = rho * vL
       e   = rho * (eL/rhoL + (a - uL) * (a - pL/qL))
    ELSE
       ! right state
       rho = qR / (a - amax)
       mu  = rho * a
       mv  = rho * vR
       e   = rho * (eR/rhoR + (a - uR) * (a - pR/qR))
    END IF
  END SUBROUTINE SetIntermediateState

  !> \private compute characteristic variables for 1D transport
  ELEMENTAL SUBROUTINE SetCharVars1d(gamma,rho1,rho2,u1,u2,P1,P2,l1,l3, &
                                     xvar1,xvar2,xvar3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,P1,P2,l1,l3
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3
    !------------------------------------------------------------------------!
    REAL              :: gammadu_cs,dlnP
    !------------------------------------------------------------------------!
    dlnP = LOG(P2/P1)      ! = LOG(P2)-LOG(P1)
    ! gamma/cs = 2*gamma / (v+cs - (v-cs))
    gammadu_cs = 2*gamma*(u2-u1) / (l3-l1)
    ! characteristic variables
    xvar1 = dlnP - gammadu_cs
    xvar2 = dlnP - gamma * LOG(rho2/rho1)
    xvar3 = dlnP + gammadu_cs
  END SUBROUTINE SetCharVars1d

  !> \private compute characteristic variables for 2D transport
  ELEMENTAL SUBROUTINE SetCharVars2d(gamma,rho1,rho2,u1,u2,v1,v2,P1,P2,l1,l4, &
                                     xvar1,xvar2,xvar3,xvar4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,P1,P2,l1,l4
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4
    !------------------------------------------------------------------------!
    CALL SetCharVars1d(gamma,rho1,rho2,u1,u2,P1,P2,l1,l4,xvar1,xvar2,xvar4)
    xvar3 = (v2-v1)
  END SUBROUTINE SetCharVars2d

  !> \private compute characteristic variables for 3D transport
  ELEMENTAL SUBROUTINE SetCharVars3d(gamma,rho1,rho2,u1,u2,v1,v2,w1,w2,P1,P2, &
                                     l1,l5,xvar1,xvar2,xvar3,xvar4,xvar5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,w1,w2,P1,P2,l1,l5
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4,xvar5
    !------------------------------------------------------------------------!
    CALL SetCharVars2d(gamma,rho1,rho2,u1,u2,v1,v2,P1,P2,l1,l5,xvar1,xvar2,xvar3,xvar5)
    xvar4 = (w2-w1)
  END SUBROUTINE SetCharVars3d

  !> \private extrapolate primitive variables using characteristic pseudo pevariables
  !! 1D transport
  ELEMENTAL SUBROUTINE SetBoundaryData1d(delta,gamma,rho1,u1,P1,xvar1, &
       xvar2,xvar3,rho2,u2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: delta
    REAL, INTENT(IN)    :: gamma,rho1,u1,P1,xvar1,xvar2,xvar3
    REAL, INTENT(OUT)   :: rho2,u2,P2
    !------------------------------------------------------------------------!
    REAL                :: dlnP,cs
    !------------------------------------------------------------------------!
    dlnP = 0.5*(xvar3+xvar1)
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(delta/gamma*(dlnP-xvar2))
    P2   = P1 * EXP(delta*dlnP)
    ! compute sound speed with arithmetic mean of density and pressure
    ! (the factor 0.5 cancels out, because cs depends on the quotient P/rho)
    cs = GetSoundSpeed(gamma,rho1+rho2,P1+P2)
    u2 = u1 + 0.5*delta*cs/gamma * (xvar3-xvar1)
  END SUBROUTINE SetBoundaryData1d

  !> \private extrapolate primitive variables using characteristic pseudo pevariables
  !! 2D transport
  ELEMENTAL SUBROUTINE SetBoundaryData2d(delta,gamma,rho1,u1,v1,P1,xvar1, &
       xvar2,xvar3,xvar4,rho2,u2,v2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: delta
    REAL, INTENT(IN)    :: gamma,rho1,u1,v1,P1,xvar1,xvar2,xvar3,xvar4
    REAL, INTENT(OUT)   :: rho2,u2,v2,P2
    !------------------------------------------------------------------------!
    CALL SetBoundaryData1d(delta,gamma,rho1,u1,P1,xvar1,xvar2,xvar4,rho2,u2,P2)
    v2 = v1 + delta*xvar3
  END SUBROUTINE SetBoundaryData2d

  !> \private extrapolate primitive variables using characteristic pseudo pevariables
  !! 2D transport
  ELEMENTAL SUBROUTINE SetBoundaryData3d(delta,gamma,rho1,u1,v1,w1,P1,xvar1, &
       xvar2,xvar3,xvar4,xvar5,rho2,u2,v2,w2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: delta
    REAL, INTENT(IN)    :: gamma,rho1,u1,v1,w1,P1,xvar1,xvar2,xvar3,xvar4,xvar5
    REAL, INTENT(OUT)   :: rho2,u2,v2,w2,P2
    !------------------------------------------------------------------------!
    CALL SetBoundaryData2d(delta,gamma,rho1,u1,v1,P1,xvar1,xvar2,xvar3,xvar5,rho2,u2,v2,P2)
    w2 = w1 + delta*xvar4
  END SUBROUTINE SetBoundaryData3d

  ! \todo NOT VERIFIED
  !! only for farfield boundary conditions
  ELEMENTAL SUBROUTINE Prim2Riemann1d(gamma,rho,vx,p,l1,l2,l3,Rminus,Rs,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,p,l1,l2,l3
    REAL, INTENT(OUT) :: Rminus,Rs,Rplus
    !------------------------------------------------------------------------!
    REAL              :: cs
    !------------------------------------------------------------------------!
    cs = l3-l2 ! l2 = v, l3 = v+cs
    ! compute 1st Riemann invariant (R+)
    Rplus = vx + 2./(gamma-1.0) * cs
    ! compute 2st Riemann invariant (R-)
    Rminus = vx - 2./(gamma-1.0) * cs
    ! compute entropy
    Rs = p/rho**gamma
  END SUBROUTINE Prim2Riemann1d


  ! \todo NOT VERIFIED
  !! only for farfield boundary conditions
  ELEMENTAL SUBROUTINE Prim2Riemann2d(gamma,rho,vx,vy,p,&
                                        l1,l2,l3,l4,Rminus,Rs,Rvt,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,vy,p,l1,l2,l3,l4
    REAL, INTENT(OUT) :: Rminus,Rs,Rvt,Rplus
    !------------------------------------------------------------------------!
    REAL              :: cs
    !------------------------------------------------------------------------!
    CALL Prim2Riemann1d(gamma,rho,vx,p,l1,l2,l4,Rminus,Rs,Rplus)
    ! tangential velocities
    Rvt = vy
  END SUBROUTINE Prim2Riemann2d


  ! \todo NOT VERIFIED
  !! only for farfield boundary conditions
  ELEMENTAL SUBROUTINE Prim2Riemann3d(gamma,rho,vx,vy,vz,p,&
                                        l1,l2,l3,l4,l5,Rminus,Rs,Rvt,Rwt,Rplus)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,vx,vy,vz,p,l1,l2,l3,l4,l5
    REAL, INTENT(OUT) :: Rminus,Rs,Rvt,Rwt,Rplus
    !------------------------------------------------------------------------!
    REAL              :: cs
    !------------------------------------------------------------------------!
    CALL Prim2Riemann2d(gamma,rho,vx,vy,p,l1,l2,l3,l5,Rminus,Rs,Rvt,Rplus)
    Rwt = vz
  END SUBROUTINE Prim2Riemann3d

  ! \todo NOT VERIFIED
  !! only for farfield boundary conditions
  ELEMENTAL SUBROUTINE Riemann2Prim1d(gamma,Rminus,Rs,Rplus,rho,vx,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rplus
    REAL, INTENT(OUT) :: rho,vx,p
    !------------------------------------------------------------------------!
    REAL              :: cs2gam
    !------------------------------------------------------------------------!
    ! normal velocity
    vx = 0.5*(Rplus+Rminus)
    ! cs**2 / gamma
    cs2gam = (0.25*(gamma-1.0)*(Rplus-Rminus))**2 / gamma
    ! density
    rho = (cs2gam/Rs)**(1./(gamma-1.0))
    ! pressure
    p = cs2gam * rho
  END SUBROUTINE Riemann2Prim1d

  ! \todo NOT VERIFIED
  !! only for farfield boundary conditions
  ELEMENTAL SUBROUTINE Riemann2Prim2d(gamma,Rminus,Rs,Rvt,Rplus,&
                                            rho,vx,vy,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rvt,Rplus
    REAL, INTENT(OUT) :: rho,vx,vy,p
    !------------------------------------------------------------------------!
    REAL              :: cs2gam
    !------------------------------------------------------------------------!
    ! tangential velocity
    vy = Rvt

    CALL Riemann2Prim1d(gamma,Rminus,Rs,Rplus,rho,vx,p)
  END SUBROUTINE Riemann2Prim2d


  ! \todo NOT VERIFIED
  !! only for farfield boundary conditions
  ELEMENTAL SUBROUTINE Riemann2Prim3d(gamma,Rminus,Rs,Rvt,Rwt,Rplus,&
       rho,vx,vy,vz,p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rvt,Rwt,Rplus
    REAL, INTENT(OUT) :: rho,vx,vy,vz,p
    !------------------------------------------------------------------------!
    REAL              :: cs2gam
    !------------------------------------------------------------------------!
    ! tangential velocity
    vz = Rwt

    CALL Riemann2Prim2d(gamma,Rminus,Rs,Rvt,Rplus,rho,vx,vy,p)
  END SUBROUTINE Riemann2Prim3d

  !> \private set mass, 1D momentum and energy flux for transport along the 1st dimension
  ELEMENTAL SUBROUTINE SetFlux1d(rho,u,P,mu,E,f1,f2,f3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,u,P,mu,E
    REAL, INTENT(OUT) :: f1,f2,f3
    !------------------------------------------------------------------------!
    f1 = rho*u
    f2 = mu*u + P
    f3 = (E+P)*u
  END SUBROUTINE SetFlux1d

  !> \private set mass, 2D momentum and energy flux for transport along the 1st dimension
  ELEMENTAL SUBROUTINE SetFlux2d(rho,u,P,mu,mv,E,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,u,P,mu,mv,E
    REAL, INTENT(OUT) :: f1,f2,f3,f4
    !------------------------------------------------------------------------!
    CALL SetFlux1d(rho,u,P,mu,E,f1,f2,f4)
    f3 = mv*u
  END SUBROUTINE SetFlux2d

  !> \private set mass, 3D momentum and energy flux for transport along the 1st dimension
  ELEMENTAL SUBROUTINE SetFlux3d(rho,u,P,mu,mv,mw,E,f1,f2,f3,f4,f5)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,u,P,mu,mv,mw,E
    REAL, INTENT(OUT) :: f1,f2,f3,f4,f5
    !------------------------------------------------------------------------!
    CALL SetFlux2d(rho,u,P,mu,mv,E,f1,f2,f3,f5)
    f4 = mw*u
  END SUBROUTINE SetFlux3d

  !> \private Convert from 1D conservative to primitive variables
  ELEMENTAL SUBROUTINE Cons2Prim1d(gamma,rho_in,mu,E,rho_out,u,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,E
    REAL, INTENT(OUT) :: rho_out,u,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * mu*mu)
  END SUBROUTINE Cons2Prim1d

  !> \private Convert from 2D conservative to primitive variables
  ELEMENTAL SUBROUTINE Cons2Prim2d(gamma,rho_in,mu,mv,E,rho_out,u,v,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,E
    REAL, INTENT(OUT) :: rho_out,u,v,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu + mv*mv))
  END SUBROUTINE Cons2Prim2d

  !> \private Convert from 3D conservative to primitive variables
  ELEMENTAL SUBROUTINE Cons2Prim3d(gamma,rho_in,mu,mv,mw,E,rho_out,u,v,w,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,mw,E
    REAL, INTENT(OUT) :: rho_out,u,v,w,P
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    w = mw * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu + mv*mv + mw*mw))
  END SUBROUTINE Cons2Prim3d

  !> \private Convert from 1D primitive to conservative variables
  ELEMENTAL SUBROUTINE Prim2Cons1d(gamma,rho_in,u,P,rho_out,mu,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,P
    REAL, INTENT(OUT) :: rho_out,mu,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    E = P/(gamma-1.) + 0.5 * rho_in * u*u
  END SUBROUTINE Prim2Cons1d

  !> \private Convert from 2D primitive to conservative variables
  ELEMENTAL SUBROUTINE Prim2Cons2d(gamma,rho_in,u,v,P,rho_out,mu,mv,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    E = P/(gamma-1.) + 0.5 * rho_in * (u*u + v*v)
  END SUBROUTINE Prim2Cons2d

  !> \private Convert from 3D primitive to conservative variables
  ELEMENTAL SUBROUTINE Prim2Cons3d(gamma,rho_in,u,v,w,P,rho_out,mu,mv,mw,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,w,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,mw,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    mw = rho_in * w
    E = P/(gamma-1.) + 0.5 * rho_in * (u*u + v*v + w*w)
  END SUBROUTINE Prim2Cons3d

  !> geometrical momentum source terms
  !! P is the either isothermal pressure rho*cs**2 or the real pressure.
  !!
  !! \attention These elemental functions exist multiple times for performance
  !!  reasons (inlining). Please keep this in mind for changes.
  !!  Other modules with this function:
  !!      - physics_eulerisotherm_mod
  !!
  !! x-momentum geometrical source term
  ELEMENTAL FUNCTION GetGeometricalSourceX(cxyx,cxzx,cyxy,czxz,vx,vy,vz,P,my,mz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cxyx,cxzx,cyxy,czxz,vx,vy,vz,P,my,mz
    REAL :: GetGeometricalSourceX
    !------------------------------------------------------------------------!
    GetGeometricalSourceX = my * (cyxy*vy - cxyx*vx) + mz * (czxz*vz - cxzx*vx) &
                            + (cyxy + czxz) * P
  END FUNCTION GetGeometricalSourceX

  !> y-momentum geometrical source term
  ELEMENTAL FUNCTION GetGeometricalSourceY(cxyx,cyxy,cyzy,czyz,vx,vy,vz,P,mx,mz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cxyx,cyxy,cyzy,czyz,vx,vy,vz,P,mx,mz
    REAL :: GetGeometricalSourceY
    !------------------------------------------------------------------------!
    GetGeometricalSourceY = mz * (czyz*vz - cyzy*vy) + mx * (cxyx*vx - cyxy*vy) &
                            + (cxyx + czyz) * P
  END FUNCTION GetGeometricalSourceY

  !> z-momentum geometrical source term
  ELEMENTAL FUNCTION GetGeometricalSourceZ(cxzx,cyzy,czxz,czyz,vx,vy,vz,P,mx,my)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cxzx,cyzy,czxz,czyz,vx,vy,vz,P,mx,my
    REAL :: GetGeometricalSourceZ
    !------------------------------------------------------------------------!
    GetGeometricalSourceZ = mx * (cxzx*vx - czxz*vz) + my * (cyzy*vy - czyz*vz) &
                            + (cxzx + cyzy) * P
  END FUNCTION GetGeometricalSourceZ

END MODULE physics_euler_mod
