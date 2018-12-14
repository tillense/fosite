!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_euler.f90                                                 #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
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
    PROCEDURE :: InitPhysics_euler             !< constructor
    PROCEDURE :: PrintConfiguration_euler
    PROCEDURE :: new_statevector
    !------Convert2Primitve--------!
    PROCEDURE :: Convert2Primitive_new
    PROCEDURE :: Convert2Primitive_centsub
    PROCEDURE :: Convert2Primitive_facesub
    !------Convert2Conservative----!
    PROCEDURE :: Convert2Conservative_new
    PROCEDURE :: Convert2Conservative_centsub
    PROCEDURE :: Convert2Conservative_facesub
    !------soundspeed routines-----!
    PROCEDURE :: UpdateSoundSpeed
    !------flux routines-----------!
    PROCEDURE :: CalcFluxesX
    PROCEDURE :: CalcFluxesY
    PROCEDURE :: CalcFluxesZ
    !------fargo routines----------!
    PROCEDURE :: AddBackgroundVelocityX
    PROCEDURE :: SubtractBackgroundVelocityX
    PROCEDURE :: AddBackgroundVelocityY
    PROCEDURE :: SubtractBackgroundVelocityY
    !------HLLC routines-----------!
!    PROCEDURE :: CalcIntermediateStateX
!    PROCEDURE :: CalcIntermediateStateY

    PROCEDURE :: ExternalSources
    PROCEDURE :: GeometricalSources
    PROCEDURE :: ReflectionMasks                ! for reflecting boundaries

    ! boundarie routines
    PROCEDURE :: CalculateCharSystemX          ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemY          ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemZ          ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataX        ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataY        ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataZ        ! for absorbing boundaries
!    PROCEDURE :: CalcPrim2RiemannX        ! for farfield boundaries
!    PROCEDURE :: CalcPrim2RiemannY        ! for farfield boundaries
!    PROCEDURE :: CalcPrim2RiemannZ        ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimX        ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimY        ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimZ        ! for farfield boundaries
    PROCEDURE :: AxisMasks                 ! for axis boundaries
    PROCEDURE :: ViscositySources
    PROCEDURE :: CalcStresses_euler


    PROCEDURE :: Finalize
  END TYPE
  TYPE, EXTENDS(statevector_eulerisotherm) :: statevector_euler
    TYPE(marray_base), POINTER :: &
                               pressure => null(), &
                               energy => null()
    CONTAINS
    PROCEDURE :: AssignMArray_0
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
  SUBROUTINE InitPhysics_euler(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN)   :: Mesh
    TYPE(Dict_TYP), POINTER, INTENT(IN)   :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: err,next_idx
    !------------------------------------------------------------------------!
    ! call InitPhysics from base class
    CALL this%InitPhysics(Mesh,config,IO,EULER,problem_name)

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

    ! set names shown in the data file
    next_idx = 2
    IF (Mesh%INUM.GT.1.OR.Mesh%ROTSYM.EQ.1) THEN
      this%pvarname(next_idx) = "xvelocity"
      this%cvarname(next_idx) = "xmomentum"
      next_idx = next_idx + 1
    END IF
    IF (Mesh%JNUM.GT.1.OR.Mesh%ROTSYM.EQ.2) THEN
      this%pvarname(next_idx) = "yvelocity"
      this%cvarname(next_idx) = "ymomentum"
      next_idx = next_idx + 1
    END IF
    IF (Mesh%KNUM.GT.1.OR.Mesh%ROTSYM.EQ.3) THEN
      this%pvarname(next_idx) = "zvelocity"
      this%cvarname(next_idx) = "zmomentum"
    END IF

    ! not used in non-isotherml physics
    this%csiso = 0.0

    ! create new mesh arrays for sound speeds
    this%bccsound = marray_base()
    this%fcsound = marray_base(Mesh%NFACES)

    CALL this%EnableOutput(Mesh,config,IO)
  END SUBROUTINE InitPhysics_euler

  SUBROUTINE PrintConfiguration_euler(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%PrintConfiguration()
  END SUBROUTINE PrintConfiguration_euler

  !> \public allocate an initialize new non-isothermal state vector
  SUBROUTINE new_statevector(this,new_sv,flavour,num)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(marray_compound), POINTER :: new_sv
    INTEGER, OPTIONAL, INTENT(IN)   :: flavour,num
    !------------------------------------------------------------------------!
    ALLOCATE(statevector_euler::new_sv)
    SELECT TYPE(sv => new_sv)
    TYPE IS (statevector_euler)
      sv = statevector_euler(this,flavour,num)
    END SELECT
  END SUBROUTINE new_statevector

  PURE SUBROUTINE Convert2Primitive_new(this,cvar,pvar)
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
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Primitive_new

    !> Converts to conservative at cell centers using state vectors
  PURE SUBROUTINE Convert2Conservative_new(this,pvar,cvar)
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
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Conservative_new

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

!  !> Characteristic variables for absorbing boundary conditions
!  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: i,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i1,i2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at i
    CALL SetEigenValues(this%gamma, &
          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE),&
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
    CALL SetCharVars(this%gamma, &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
 
 END SUBROUTINE CalculateCharSystemX


  !> Characteristic variables for absorbing boundary conditions
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: j,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) ,INTENT(OUT):: lambda
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
    INTEGER           :: j1,j2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at j
    CALL SetEigenValues(this%gamma, &
          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
    CALL SetCharVars(this%gamma, &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4))

        END SUBROUTINE CalculateCharSystemY


  !> Characteristic variables for absorbing boundary conditions
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateCharSystemZ(this,Mesh,k,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: k,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
!     INTEGER           :: k1,k2
    !------------------------------------------------------------------------!
    !TODO Should not exist in 2D !
    !    ! compute eigenvalues at k
!    CALL SetEigenValues(this%gamma, &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%PRESSURE), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))
!    ! compute characteristic variables
!    k1 = k + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
!    k2 = MAX(k,k1)
!    k1 = MIN(k,k1)
!    CALL SetCharVars(this%gamma, &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%PRESSURE), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%PRESSURE), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))


  END SUBROUTINE CalculateCharSystemZ


  !> Calculate boundary data for absorbing boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: i1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER                            :: i2
    !------------------------------------------------------------------------!
    i2 = i1 + SIGN(1,dir)  ! i +/- 1 depending on the sign of dir
    CALL SetBoundaryData(this%gamma,1.0*SIGN(1,dir), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE))
  END SUBROUTINE CalculateBoundaryDataX


  !> Calculate boundary data for absorbing boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: j1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER                            :: j2
    !------------------------------------------------------------------------!
    j2 = j1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
    CALL SetBoundaryData(this%gamma,1.0*SIGN(1,dir), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,4), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%PRESSURE))
  END SUBROUTINE CalculateBoundaryDataY


  !> Calculate boundary data for absorbing boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateBoundaryDataZ(this,Mesh,k1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: k1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
!     INTEGER                            :: k2
    !------------------------------------------------------------------------!
    !TODO Should not exist in 2D !
    !    k2 = k1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
!    CALL SetBoundaryData(this%gamma,1.0*SIGN(1,dir), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%PRESSURE), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%PRESSURE))
  END SUBROUTINE CalculateBoundaryDataZ


!  !> Conversion from primitive to riemann invariants for farfield boundaries
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcPrim2RiemannX(this,Mesh,i,pvar,lambda,Rinv)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    INTEGER,                INTENT(IN) :: i
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: pvar
!    REAL,                   INTENT(INOUT), &
!      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: lambda
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: Rinv
!    !------------------------------------------------------------------------!
!    ! compute eigenvalues at i
!    CALL SetEigenValues(this%gamma, &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE),&
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5))
!    ! compute Riemann invariants
!    CALL Prim2Riemann(this%gamma, &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
!          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5))
!  END SUBROUTINE CalcPrim2RiemannX
!
!
!  !> Conversion from primitive to riemann invariants for farfield boundaries
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcPrim2RiemannY(this,Mesh,j,pvar,lambda,Rinv)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    INTEGER,                INTENT(IN) :: j
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: pvar
!    REAL,                   INTENT(INOUT), &
!      DIMENSION(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,this%VNUM) &
!                                       :: lambda
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,this%VNUM) &
!                                       :: Rinv
!    !------------------------------------------------------------------------!
!    ! compute eigenvalues at j
!    CALL SetEigenValues(this%gamma, &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%PRESSURE),&
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,5))
!    ! compute Riemann invariants
!    CALL Prim2Riemann(this%gamma, &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%PRESSURE), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4), &
!          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,5), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,5))
!  END SUBROUTINE CalcPrim2RiemannY
!
!
!  !> Conversion from primitive to riemann invariants for farfield boundaries
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcPrim2RiemannZ(this,Mesh,k,pvar,lambda,Rinv)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    INTEGER,                INTENT(IN) :: k
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: pvar
!    REAL,                   INTENT(INOUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
!                                       :: lambda
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
!                                       :: Rinv
!    !------------------------------------------------------------------------!
!    ! compute eigenvalues at k
!    CALL SetEigenValues(this%gamma, &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%PRESSURE),&
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))
!    ! compute Riemann invariants
!    CALL Prim2Riemann(this%gamma, &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%PRESSURE), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5))
!  END SUBROUTINE CalcPrim2RiemannZ
!
!
!  !> Convert Riemann invariants to primitives for farfield boundaries
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcRiemann2PrimX(this,Mesh,i,Rinv,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    INTEGER,                INTENT(IN) :: i
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: Rinv
!    REAL,                   INTENT(INOUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: pvar
!    !------------------------------------------------------------------------!
!    CALL Riemann2Prim(this%gamma, &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
!          Rinv(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,5), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
!          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%PRESSURE))
!  END SUBROUTINE CalcRiemann2PrimX
!
!
!  !> Convert Riemann invariants to primitives for farfield boundaries
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcRiemann2PrimY(this,Mesh,j,Rinv,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    INTEGER,                INTENT(IN) :: j
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,this%VNUM) &
!                                       :: Rinv
!    REAL,                   INTENT(INOUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: pvar
!    !------------------------------------------------------------------------!
!    CALL Riemann2Prim(this%gamma, &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4), &
!          Rinv(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,5), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%PRESSURE))
!  END SUBROUTINE CalcRiemann2PrimY
!
!
!  !> Convert Riemann invariants to primitives for farfield boundaries
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcRiemann2PrimZ(this,Mesh,k,Rinv,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    INTEGER,                INTENT(IN) :: k
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
!                                       :: Rinv
!    REAL,                   INTENT(INOUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: pvar
!    !------------------------------------------------------------------------!
!    CALL Riemann2Prim(this%gamma, &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
!          Rinv(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,5), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%ZVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%PRESSURE))
!  END SUBROUTINE CalcRiemann2PrimZ

  !> Calculates geometrical sources
  PURE SUBROUTINE GeometricalSources(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    ! compute geometrical source only for non-cartesian mesh
    IF (Mesh%Geometry%GetType().NE.CARTESIAN) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_euler)
          SELECT TYPE(s => sterm)
          TYPE IS(statevector_euler)
            ! no source terms
            s%density%data1d(:) = 0.0
            s%energy%data1d(:) = 0.0
            SELECT CASE(this%VDIM)
            CASE(1) ! 1D
              IF (Mesh%INUM.GT.1) THEN
                ! x-momentum
                ! vy = vz = my = mz = 0
                s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                   Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                   Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                   p%velocity%data2d(:,1),0.0,0.0, &
                   p%pressure%data1d(:), &
                   0.0,0.0)
              ELSE IF (Mesh%JNUM.GT.1) THEN
                ! y-momentum
                ! vx = vz = mx = mz = 0
                s%momentum%data2d(:,1) = GetGeometricalSourceY( &
                   Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                   Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                   0.0,p%velocity%data2d(:,1),0.0, &
                   p%pressure%data1d(:), &
                   0.0,0.0)
              ELSE IF (Mesh%KNUM.GT.1) THEN
                ! z-momentum
                ! vx = vy = mx = my = 0
                s%momentum%data2d(:,1) = GetGeometricalSourceZ( &
                   Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                   Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                   0.0,0.0,p%velocity%data2d(:,1), &
                   p%pressure%data1d(:), &
                   0.0,0.0)
              END IF
            CASE(2) ! 2D
              IF (Mesh%KNUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.3) THEN
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
              ELSE IF (Mesh%JNUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.2) THEN
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
                   0.0,c%momentum%data2d(:,2))
              ELSE IF (Mesh%INUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.1) THEN
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
                   0.0,c%momentum%data2d(:,2))
              END IF
            CASE(3) ! 3D
                ! x-momentum
                s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                   Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                   Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                   p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                   p%velocity%data2d(:,3), &
                   p%pressure%data1d(:), &
                   c%momentum%data2d(:,2),c%momentum%data2d(:,3))
                ! y-momentum
                s%momentum%data2d(:,2) = GetGeometricalSourceY( &
                   Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                   Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                   p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                   p%velocity%data2d(:,3), &
                   p%pressure%data1d(:), &
                   c%momentum%data2d(:,1),c%momentum%data2d(:,3))
                ! z-momentum
                s%momentum%data2d(:,3) = GetGeometricalSourceZ( &
                   Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                   Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                   p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                   p%velocity%data2d(:,3), &
                   p%pressure%data1d(:), &
                   c%momentum%data2d(:,1),c%momentum%data2d(:,2))
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
    END IF
  END SUBROUTINE GeometricalSources

  !> compute momentum and energy sources given an external force
  PURE SUBROUTINE ExternalSources(this,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(marray_base),       INTENT(IN)  :: accel
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    CALL this%physics_eulerisotherm%ExternalSources(accel,pvar,cvar,sterm)
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      SELECT TYPE(c => cvar)
      TYPE IS(statevector_euler)
        SELECT TYPE(s => sterm)
        TYPE IS(statevector_euler)
          s%energy%data1d(:) = SUM(c%momentum%data2d(:,:) * accel%data2d(:,:),DIM=2)
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE ExternalSources

  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL,                   INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                           :: pvar
    REAL,                   INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                                           :: btxx,btxy,btxz,btyy,btyz,btzz
    REAL,                   INTENT(OUT), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: sterm
   !------------------------------------------------------------------------!
   CALL this%physics_eulerisotherm%ViscositySources(Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
   SELECT CASE(this%VDIM)
!   CASE(1)
!     this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:)
!
!     CALL Mesh%Divergence(this%tmp(:,:,:),sterm(:,:,:,this%ENERGY))
   CASE(2)
     !compute scalar product of v and tau (x-component)
    this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btxy(:,:,:)

    !compute scalar product of v and tau (y-component)
    this%tmp1(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxy(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btyy(:,:,:)

    ! compute vector divergence of scalar product v and tau
    CALL Mesh%Divergence(this%tmp(:,:,:),this%tmp1(:,:,:), &
          sterm(:,:,:,this%ENERGY))
  CASE(3) 
    !compute scalar product of v and tau (x-component)
    this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btxy(:,:,:) & 
                    + pvar(:,:,:,this%ZVELOCITY)*btxz(:,:,:)

    !compute scalar product of v and tau (y-component)
    this%tmp1(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxy(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btyy(:,:,:) &
                    + pvar(:,:,:,this%ZVELOCITY)*btyz(:,:,:)

    !compute scalar product of v and tau (z-component)
    this%tmp2(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxz(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btyz(:,:,:) &
                    + pvar(:,:,:,this%ZVELOCITY)*btzz(:,:,:)
    ! compute vector divergence of scalar product v and tau
    CALL Mesh%Divergence(this%tmp(:,:,:),this%tmp1(:,:,:),this%tmp2(:,:,:), &
          sterm(:,:,:,this%ENERGY))
   END SELECT
 END SUBROUTINE ViscositySources




  ! identical to isothermal case. 
  PURE SUBROUTINE CalcStresses_euler(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_euler), INTENT(INOUT) :: this
    CLASS(Mesh_base), INTENT(IN)          :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,dynvis,bulkvis
    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    ! compute components of the stress tensor at cell bary centers
    ! inside the computational domain including one slice of ghost cells

    ! compute bulk viscosity first and store the result in this%tmp
    SELECT CASE (this%VDIM) 
!    CASE(1)
!      CALL Mesh%Divergence(pvar(:,:,:,this%XVELOCITY),this%tmp(:,:,:))
    CASE(2)
      CALL Mesh%Divergence(pvar(:,:,:,this%XVELOCITY),pvar(:,:,:,this%YVELOCITY),this%tmp(:,:,:))
    CASE(3)
      CALL Mesh%Divergence(pvar(:,:,:,this%XVELOCITY),pvar(:,:,:,this%YVELOCITY),pvar(:,:,:,this%ZVELOCITY),this%tmp(:,:,:))
    END SELECT
    this%tmp(:,:,:) = bulkvis(:,:,:)*this%tmp(:,:,:)

    SELECT CASE(this%VDIM)
!    CASE(1)
!      !NEC$ OUTERLOOP_UNROLL(8)
!      DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
!        DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!          !NEC$ IVDEP
!          DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
!            ! compute the diagonal elements of the stress tensor
!            btxx(i,j,k) = dynvis(i,j,k) * &
!                ((pvar(i+1,j,k,this%XVELOCITY) - pvar(i-1,j,k,this%XVELOCITY)) / Mesh%dlx%data3d(i,j,k) &
!               + this%tmp(i,j,k)
!          END DO
!        END DO
!      END DO
    CASE(2)
      !NEC$ OUTERLOOP_UNROLL(8)
      DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
        DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
          !NEC$ IVDEP
          DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
            ! compute the diagonal elements of the stress tensor
            btxx(i,j,k) = dynvis(i,j,k) * &
                ((pvar(i+1,j,k,this%XVELOCITY) - pvar(i-1,j,k,this%XVELOCITY)) / Mesh%dlx%data3d(i,j,k) &
               + 2.0 * Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY))  &
               + this%tmp(i,j,k)

            btyy(i,j,k) = dynvis(i,j,k) * &
               ( (pvar(i,j+1,k,this%YVELOCITY) - pvar(i,j-1,k,this%YVELOCITY)) / Mesh%dly%data3d(i,j,k) &
               + 2.0 * Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY))  &
               + this%tmp(i,j,k)

            ! compute the off-diagonal elements (no bulk viscosity)
            btxy(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i+1,j,k,this%YVELOCITY) - pvar(i-1,j,k,this%YVELOCITY)) / Mesh%dlx%data3d(i,j,k) &
               + (pvar(i,j+1,k,this%XVELOCITY) - pvar(i,j-1,k,this%XVELOCITY)) / Mesh%dly%data3d(i,j,k) ) &
               - Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               - Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

          END DO
        END DO
      END DO
    CASE(3)
      !NEC$ OUTERLOOP_UNROLL(8)
      DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
        DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
          !NEC$ IVDEP
          DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
            ! compute the diagonal elements of the stress tensor
            btxx(i,j,k) = dynvis(i,j,k) * &
                ((pvar(i+1,j,k,this%XVELOCITY) - pvar(i-1,j,k,this%XVELOCITY)) / Mesh%dlx%data3d(i,j,k) &
               + 2.0 * Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY)  &
               + 2.0 * Mesh%cxzx%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) ) &
               + this%tmp(i,j,k)

            btyy(i,j,k) = dynvis(i,j,k) * &
               ( (pvar(i,j+1,k,this%YVELOCITY) - pvar(i,j-1,k,this%YVELOCITY)) / Mesh%dly%data3d(i,j,k) &
               + 2.0 * Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY)  &
               + 2.0 * Mesh%cyzy%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) ) &
               + this%tmp(i,j,k)

            btzz(i,j,k) = dynvis(i,j,k) * &
               ( (pvar(i,j,k+1,this%ZVELOCITY) - pvar(i,j,k-1,this%ZVELOCITY)) / Mesh%dlz%data3d(i,j,k) &
               + 2.0 * Mesh%czxz%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               + 2.0 * Mesh%czyz%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) ) &
               + this%tmp(i,j,k)

            ! compute the off-diagonal elements (no bulk viscosity)
            btxy(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i+1,j,k,this%YVELOCITY) - pvar(i-1,j,k,this%YVELOCITY)) / Mesh%dlx%data3d(i,j,k) &
               + (pvar(i,j+1,k,this%XVELOCITY) - pvar(i,j-1,k,this%XVELOCITY)) / Mesh%dly%data3d(i,j,k) ) &
               - Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               - Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

            btxz(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i+1,j,k,this%ZVELOCITY) - pvar(i-1,j,k,this%ZVELOCITY)) / Mesh%dlx%data3d(i,j,k) &
               + (pvar(i,j,k+1,this%XVELOCITY) - pvar(i,j,k-1,this%XVELOCITY)) / Mesh%dlz%data3d(i,j,k) ) &
               - Mesh%czxz%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) &
               - Mesh%cxzx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) )

            btyz(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i,j,k+1,this%YVELOCITY) - pvar(i,j,k-1,this%YVELOCITY)) / Mesh%dlz%data3d(i,j,k) &
               + (pvar(i,j+1,k,this%ZVELOCITY) - pvar(i,j-1,k,this%ZVELOCITY)) / Mesh%dly%data3d(i,j,k) ) &
               - Mesh%czyz%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) &
               - Mesh%cyzy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

          END DO
        END DO
      END DO
    END SELECT
  END SUBROUTINE CalcStresses_euler


  !> Convert to from conservative to primitive variables at cell-centers
  PURE SUBROUTINE Convert2Primitive_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(IN)  :: cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Cons2Prim(this%gamma, &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ENERGY),    &
                     pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE))
    CASE(2) ! 2D velocity / momentum
      CALL Cons2Prim(this%gamma, &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ENERGY),    &
                     pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE))
    CASE(3) ! 3D velocity / momentum
      CALL Cons2Prim(this%gamma, &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ENERGY),    &
                     pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE))
    END SELECT
  END SUBROUTINE Convert2Primitive_centsub

  !> Convert to from conservative to primitive variables at cell-faces
  PURE SUBROUTINE Convert2Primitive_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Cons2Prim(this%gamma, &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY),    &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE))
    CASE(2) ! 2D velocity / momentum
      CALL Cons2Prim(this%gamma, &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY),    &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE))
    CASE(3) ! 3D velocity / momentum
      CALL Cons2Prim(this%gamma, &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY),    &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE))
    END SELECT
  END SUBROUTINE Convert2Primitive_facesub

  !> Convert to from primitve to conservative variables at cell-centers
  PURE SUBROUTINE Convert2Conservative_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(IN)  :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Prim2Cons(this%gamma, &
                     pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE),  &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ENERGY))
    CASE(2) ! 2D velocity / momentum
      CALL Prim2Cons(this%gamma, &
                     pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE),  &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ENERGY))
    CASE(3) ! 3D velocity / momentum
      CALL Prim2Cons(this%gamma, &
                     pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE),  &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ENERGY))
    END SELECT
  END SUBROUTINE Convert2Conservative_centsub

  !> Convert to from primitve to conservative variables at cell-faces
  PURE SUBROUTINE Convert2Conservative_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: prim
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Prim2Cons(this%gamma, &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE) , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY))
    CASE(2) ! 2D velocity / momentum
      CALL Prim2Cons(this%gamma, &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE) , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY))
    CASE(3) ! 3D velocity / momentum
      CALL Prim2Cons(this%gamma, &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE) , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY))
    END SELECT
  END SUBROUTINE Convert2Conservative_facesub


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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                            INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_yvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! ATTENTION: don't change the order; on the RHS of the first
            !            assignment there must be the old momentum
            cvar(i,j,k,this%ENERGY) = cvar(i,j,k,this%ENERGY) &
                                  + w(i,k)*(cvar(i,j,k,this%YMOMENTUM) &
                                  + 0.5*cvar(i,j,k,this%DENSITY)*w(i,k))
            pvar(i,j,k,this%YVELOCITY) = pvar(i,j,k,this%YVELOCITY) + w(i,k)
            cvar(i,j,k,this%YMOMENTUM) = cvar(i,j,k,this%YMOMENTUM) &
                                     + cvar(i,j,k,this%DENSITY)*w(i,k)
          END DO
        END DO
      END DO
      this%transformed_yvelocity = .FALSE.
    END IF
  END SUBROUTINE AddBackgroundVelocityY

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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                            INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_xvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! ATTENTION: don't change the order; on the RHS of the first
             !            assignment there must be the old momentum
             cvar(i,j,k,this%ENERGY) = cvar(i,j,k,this%ENERGY) &
                                     + w(j,k)*(cvar(i,j,k,this%XMOMENTUM) &
                                     + 0.5*cvar(i,j,k,this%DENSITY)*w(j,k))
             pvar(i,j,k,this%XVELOCITY) = pvar(i,j,k,this%XVELOCITY) + w(j,k)
             cvar(i,j,k,this%XMOMENTUM) = cvar(i,j,k,this%XMOMENTUM) &
                                      + cvar(i,j,k,this%DENSITY)*w(j,k)
          END DO
        END DO
      END DO
      this%transformed_xvelocity = .FALSE.
    END IF
  END SUBROUTINE AddBackgroundVelocityX


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
  PURE SUBROUTINE SubtractBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                            INTENT(IN)    :: w
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                            INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_yvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! ATTENTION: don't change the order; on the RHS of the first
            !            assignment there must be the old momentum
            cvar(i,j,k,this%ENERGY) = cvar(i,j,k,this%ENERGY) &
                                  - w(i,k)*(cvar(i,j,k,this%YMOMENTUM) &
                                  - 0.5*cvar(i,j,k,this%DENSITY)*w(i,k))
            pvar(i,j,k,this%YVELOCITY) = pvar(i,j,k,this%YVELOCITY) - w(i,k)
            cvar(i,j,k,this%YMOMENTUM) = cvar(i,j,k,this%YMOMENTUM) &
                                     - cvar(i,j,k,this%DENSITY)*w(i,k)
          END DO
        END DO
      END DO
      this%transformed_yvelocity = .TRUE.
    END IF
  END SUBROUTINE SubtractBackgroundVelocityY

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
  PURE SUBROUTINE SubtractBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                            INTENT(IN)    :: w
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                            INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_xvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! ATTENTION: don't change the order; on the RHS of the first
            !            assignment there must be the old momentum
            cvar(i,j,k,this%ENERGY) = cvar(i,j,k,this%ENERGY) &
                                  - w(j,k)*(cvar(i,j,k,this%XMOMENTUM) &
                                  - 0.5*cvar(i,j,k,this%DENSITY)*w(j,k))
            pvar(i,j,k,this%XVELOCITY) = pvar(i,j,k,this%XVELOCITY) - w(j,k)
            cvar(i,j,k,this%XMOMENTUM) = cvar(i,j,k,this%XMOMENTUM) &
                                     - cvar(i,j,k,this%DENSITY)*w(j,k)
          END DO
        END DO
      END DO
      this%transformed_xvelocity = .TRUE.
    END IF
  END SUBROUTINE SubtractBackgroundVelocityX


  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN)  :: this
    LOGICAL, DIMENSION(this%VNUM+this%PNUM), &
                            INTENT(OUT) :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .FALSE.
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
    ! bottomer / topper boundary
    reflZ(this%DENSITY)   = .FALSE.
    reflZ(this%XVELOCITY) = .FALSE.
    reflZ(this%YVELOCITY) = .FALSE.
    reflZ(this%PRESSURE)  = .FALSE.
  END SUBROUTINE ReflectionMasks

  ! \todo \warning not clear since 3D version if this is correct. Most probably
  ! axis boundaries can be applied always in two dimensions. Now only x-y plane
  PURE SUBROUTINE AxisMasks(this,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: this
    LOGICAL, DIMENSION(this%VNUM+this%PNUM), &
                            INTENT(OUT) :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .TRUE.
    reflX(this%PRESSURE)  = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.        !old: .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%PRESSURE)  = .FALSE.
    ! bottomer / topper boundary
    reflZ(this%DENSITY)   = .FALSE.
    reflZ(this%XVELOCITY) = .FALSE.
    reflZ(this%YVELOCITY) = .FALSE.
    reflZ(this%PRESSURE)  = .FALSE.
  END SUBROUTINE AxisMasks

  PURE SUBROUTINE UpdateSoundSpeed(this,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(INOUT) :: this
    CLASS(statevector_euler),INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    IF (pvar%density%RANK.EQ.0) THEN 
      this%bccsound%data1d(:) = GetSoundSpeed(this%gamma, &
            pvar%density%data1d(:),pvar%pressure%data1d(:))
    ELSE
      this%fcsound%data1d(:) = GetSoundSpeed(this%gamma, &
            pvar%density%data1d(:),pvar%pressure%data1d(:))
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
  !! \attention This is not a class member itself, instead its an ordinary
  !!            module procedure. The function name is overloaded with
  !!            the class name.
  FUNCTION CreateStateVector(Physics,flavour,num) RESULT(new_sv)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(physics_euler), INTENT(IN) :: Physics
    INTEGER, OPTIONAL, INTENT(IN) :: flavour,num
    TYPE(statevector_euler) :: new_sv
    !-------------------------------------------------------------------!
    ! call inherited function
    new_sv = statevector_eulerisotherm(Physics,flavour,num)
    ! add entries specific for euler physics
    SELECT CASE(flavour)
    CASE(PRIMITIVE)
      ! allocate memory for pressure mesh array
      ALLOCATE(new_sv%pressure)
      new_sv%pressure  = marray_base(num)           ! one/num scalars
      ! append to compound
      CALL new_sv%AppendMArray(new_sv%pressure)
    CASE(CONSERVATIVE)
      ! allocate memory for energy mesh array
      ALLOCATE(new_sv%energy)
      new_sv%energy  = marray_base(num)             ! one/num scalars
      ! append to compound
      CALL new_sv%AppendMArray(new_sv%energy)
    CASE DEFAULT
      CALL Physics%Warning("physics_euler::CreateStateVector", "incomplete state vector")
    END SELECT
  END FUNCTION CreateStateVector

  !> assigns one state vector to another state vector
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(statevector_euler),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    CALL this%statevector_eulerisotherm%AssignMArray_0(ma)
    IF (SIZE(this%data1d).GT.0) THEN
      SELECT TYPE(src => ma)
      CLASS IS(statevector_euler)
        SELECT CASE(this%flavour)
        CASE(PRIMITIVE)
          ! pressure is the third item
          this%pressure => this%GetItem(this%NextItem(this%NextItem(this%FirstItem())))
        CASE(CONSERVATIVE)
          ! energy is the third item
          this%energy => this%GetItem(this%NextItem(this%NextItem(this%FirstItem())))
        CASE DEFAULT
          ! error, this should not happen
        END SELECT
      CLASS DEFAULT
        ! error, this should not happen
      END SELECT
    END IF
  END SUBROUTINE AssignMArray_0


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

  !> \todo NOT VERIFIED
  !! only for boundary conditions - absorbing
  ELEMENTAL SUBROUTINE SetEigenValues(gamma,rho,v,P,l1,l2,l3,l4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho,v,P
    REAL, INTENT(OUT) :: l1,l2,l3,l4
    !------------------------------------------------------------------------!
    REAL :: cs
    !------------------------------------------------------------------------!
    ! adiabatic sound speed
    cs = GetSoundSpeed(gamma,rho,P)
    ! call subroutine for isothermal case with the adiabatic sound speed
    l1 = v - cs
    l2 = v
    l3 = v
    l4 = v + cs
  END SUBROUTINE SetEigenValues

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

  ! \todo NOT VERIFIED
  !! only for absorbing boundary conditions
  ELEMENTAL SUBROUTINE SetCharVars(gamma,rho1,rho2,u1,u2,v1,v2,P1,P2, &
       l1,l2,l3,l4,xvar1,xvar2,xvar3,xvar4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho1,rho2,u1,u2,v1,v2,P1,P2,l1,l2,l3,l4
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4
    !------------------------------------------------------------------------!
    REAL              :: gamcs,dlnP,du
    !------------------------------------------------------------------------!
    gamcs= gamma / (l4-l1) ! = 2*gamma/cs
    dlnP = LOG(P2/P1)         ! = LOG(P2)-LOG(P1)
    du   = u2-u1
    ! characteristic variables
    xvar1 = dlnP - gamcs * du
    xvar2 = -dlnP + gamma * LOG(rho2/rho1)
    xvar3 = (v2-v1)
    xvar4 = dlnP + gamcs * du
  END SUBROUTINE SetCharVars

  ! \todo NOT VERIFIED
  !! only for absorbing boundary conditions
  ELEMENTAL SUBROUTINE SetBoundaryData(gamma,dir,rho1,u1,v1,P1,xvar1, &
       xvar2,xvar3,xvar4,rho2,u2,v2,P2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,dir,rho1,u1,v1,P1,xvar1,xvar2,xvar3,xvar4
    REAL, INTENT(OUT) :: rho2,u2,v2,P2
    !------------------------------------------------------------------------!
    REAL              :: dlnP,csgam
    !------------------------------------------------------------------------!
    dlnP = 0.5 * (xvar4+xvar1)
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(dir*(dlnP-xvar2)/gamma)
    P2   = P1 * EXP(dir*dlnP)
    csgam= GetSoundSpeed(gamma,rho1+rho2,P1+P2) / gamma
    u2   = u1 + dir*csgam * 0.5*(xvar4-xvar1)
    v2   = v1 + dir*xvar3
  END SUBROUTINE SetBoundaryData

!  ! \todo NOT VERIFIED
!  !! only for farfield boundary conditions
!  ELEMENTAL SUBROUTINE Prim2Riemann(gamma,rho,vx,vy,vz,p,&
!                                        l1,l2,l3,l4,l5,Rminus,Rs,Rvt,Rwt,Rplus)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: gamma,rho,vx,vy,vz,p,l1,l2,l3,l4,l5
!    REAL, INTENT(OUT) :: Rminus,Rs,Rvt,Rwt,Rplus
!    !------------------------------------------------------------------------!
!    REAL              :: cs
!    !------------------------------------------------------------------------!
!    cs = l5-l2 ! l2 = v, l5 = v+cs
!    ! compute 1st Riemann invariant (R+)
!    Rplus = vx + 2./(gamma-1.0) * cs
!    ! compute 2st Riemann invariant (R-)
!    Rminus = vx - 2./(gamma-1.0) * cs
!    ! compute entropy
!    Rs = p/rho**gamma
!    ! tangential velocities
!    Rvt = vy
!    Rwt = vz
!  END SUBROUTINE Prim2Riemann

!  ! \todo NOT VERIFIED
!  !! only for farfield boundary conditions
!  ELEMENTAL SUBROUTINE Riemann2Prim(gamma,Rminus,Rs,Rvt,Rwt,Rplus,&
!       rho,vx,vy,vz,p)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: gamma,Rminus,Rs,Rvt,Rwt,Rplus
!    REAL, INTENT(OUT) :: rho,vx,vy,vz,p
!    !------------------------------------------------------------------------!
!    REAL              :: cs2gam
!    !------------------------------------------------------------------------!
!    ! tangential velocity
!    vy = Rvt
!    vz = Rwt
!    ! normal velocity
!    vx = 0.5*(Rplus+Rminus)
!    ! cs**2 / gamma
!    cs2gam = (0.25*(gamma-1.0)*(Rplus-Rminus))**2 / gamma
!    ! density
!    rho = (cs2gam/Rs)**(1./(gamma-1.0))
!    ! pressure
!    p = cs2gam * rho
!  END SUBROUTINE Riemann2Prim

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
