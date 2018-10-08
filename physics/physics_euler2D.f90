!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_euler2D.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
!! \brief basic module for 2D Euler equations
!----------------------------------------------------------------------------!
MODULE physics_euler2D_mod
  USE physics_base_mod
  USE physics_euler2Dit_mod, ONLY: physics_euler2Dit 
  USE mesh_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 4              ! number of variables
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D"
  !--------------------------------------------------------------------------!
  TYPE,  EXTENDS(physics_euler2Dit) :: physics_euler2D
  CONTAINS
    PROCEDURE :: InitPhysics_euler2D             !< constructor
    !------Convert2Primitve--------!
    PROCEDURE :: Convert2Primitive_centsub
    PROCEDURE :: Convert2Primitive_facesub
    !------Convert2Conservative----!
    PROCEDURE :: Convert2Conservative_centsub
    PROCEDURE :: Convert2Conservative_facesub
    !------soundspeed routines-----!
    PROCEDURE :: UpdateSoundSpeed_center
    PROCEDURE :: UpdateSoundSpeed_faces
    !------flux routines-----------!
    PROCEDURE :: CalcFluxesX
    PROCEDURE :: CalcFluxesY
    PROCEDURE :: CalcFluxesZ                    ! empty routine
    !------fargo routines----------!
    PROCEDURE :: AddBackgroundVelocity
    PROCEDURE :: SubtractBackgroundVelocity
    !------HLLC routines-----------!
!    PROCEDURE :: CalcIntermediateStateX
!    PROCEDURE :: CalcIntermediateStateY

    PROCEDURE :: ExternalSources
    PROCEDURE :: GeometricalSources_center
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
!    PROCEDURE :: GeometricalSources_faces
    PROCEDURE :: ViscositySources
    PROCEDURE :: ViscositySources_euler2D
    PROCEDURE :: CalcStresses_euler


    PROCEDURE :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       physics_euler2D, &
       ! elemental procedures
       GetSoundSpeed, &
       SetEigenValues
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2D(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN)   :: Mesh
    TYPE(Dict_TYP), POINTER, INTENT(IN)   :: config, IO

    !------------------------------------------------------------------------!
!    IF (PRESENT(pname)) THEN
!       CALL this%InitPhysics(problem,pname,num_var)
!    ELSE
!       CALL this%InitPhysics(problem,problem_name,num_var)
!    END IF
    CALL this%InitPhysics(Mesh,config,IO,EULER2D,problem_name,num_var)
    ! set array indices
    this%DENSITY   = 1                                 ! mass density        !
    this%PRESSURE  = num_var                           ! pressure            !
    this%ENERGY    = num_var                           ! total energy        !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%pvarname(this%PRESSURE)  = "pressure"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%cvarname(this%ENERGY)    = "energy"
    this%DIM = 2
  END SUBROUTINE InitPhysics_euler2D

  !> Calculate Fluxes in x-direction
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalcFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: xfluxes
    !------------------------------------------------------------------------!
    CALL SetFlux( &
         prim(:,:,:,nmin:nmax,this%DENSITY),prim(:,:,:,nmin:nmax,this%XVELOCITY), &
         prim(:,:,:,nmin:nmax,this%PRESSURE),cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
         cons(:,:,:,nmin:nmax,this%YMOMENTUM),cons(:,:,:,nmin:nmax,this%ENERGY), &
         xfluxes(:,:,:,nmin:nmax,this%DENSITY),xfluxes(:,:,:,nmin:nmax,this%XMOMENTUM), &
         xfluxes(:,:,:,nmin:nmax,this%YMOMENTUM),xfluxes(:,:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesX

  !> Calculate Fluxes in y-direction
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalcFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: yfluxes
    !------------------------------------------------------------------------!
    CALL SetFlux( &
         prim(:,:,:,nmin:nmax,this%DENSITY),prim(:,:,:,nmin:nmax,this%YVELOCITY), &
         prim(:,:,:,nmin:nmax,this%PRESSURE),cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
         cons(:,:,:,nmin:nmax,this%XMOMENTUM),cons(:,:,:,nmin:nmax,this%ENERGY), &
         yfluxes(:,:,:,nmin:nmax,this%DENSITY),yfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
         yfluxes(:,:,:,nmin:nmax,this%XMOMENTUM),yfluxes(:,:,:,nmin:nmax,this%ENERGY))
  END SUBROUTINE CalcFluxesY

  !> Calculate Fluxes in z-direction
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalcFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: zfluxes
    !------------------------------------------------------------------------!
    ! routine does not exist in 2D
  END SUBROUTINE CalcFluxesZ

!  !> Reconstruction of the intermediate state for HLLC
!  !\todo NOT VERIFIED
!  PURE SUBROUTINE CalcIntermediateStateX(this,Mesh,prim,cons, &
!                  amin,amax,cstar,astar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler2D), INTENT(IN) :: this
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
!    CLASS(physics_euler2D), INTENT(INOUT) :: this
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
    CLASS(physics_euler2D), INTENT(IN)    :: this
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
    CLASS(physics_euler2D), INTENT(IN)    :: this
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
    CLASS(physics_euler2D), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: k,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
    INTEGER           :: k1,k2
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
    CLASS(physics_euler2D), INTENT(IN) :: this
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
    CLASS(physics_euler2D), INTENT(IN) :: this
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
    CLASS(physics_euler2D), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: k1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER                            :: k2
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
!    CLASS(physics_euler2D), INTENT(IN) :: this
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
!    CLASS(physics_euler2D), INTENT(IN) :: this
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
!    CLASS(physics_euler2D), INTENT(IN) :: this
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
!    CLASS(physics_euler2D), INTENT(IN) :: this
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
!    CLASS(physics_euler2D), INTENT(IN) :: this
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
!    CLASS(physics_euler2D), INTENT(IN) :: this
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


  ! \todo HIER BLICKE ICH NICHT DURCH!!! HIER FEHLEN SICHERLICH FAKTOREN
  !> Calculate geometrical sources at the center
  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(IN)    :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT)   :: sterm
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    ! compute geometrical source only for non-cartesian mesh except for the
    ! EULER2D_IAMROT case for which geometrical sources are always necessary.
    IF ((Mesh%Geometry%GetType().NE.CARTESIAN)) THEN

      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
           DO i=Mesh%IGMIN,Mesh%IGMAX
              sterm(i,j,k,this%DENSITY) = 0.
              sterm(i,j,k,this%ENERGY)  = 0.
              !geometrical source terms in momentum equations
              sterm(i,j,k,this%XMOMENTUM) = MomentumSourcesX( &
                cvar(i,j,k,this%YMOMENTUM), &
                0., & !ZMOMENTUM for 2D is 0
                pvar(i,j,k,this%XVELOCITY), &
                pvar(i,j,k,this%YVELOCITY), &
                0., & !ZVELOCITY for 2D is 0
                pvar(i,j,k,this%PRESSURE), &
                Mesh%cxyx%bcenter(i,j,k), &
                Mesh%cyxy%bcenter(i,j,k), &
                Mesh%czxz%bcenter(i,j,k), &
                Mesh%cxzx%bcenter(i,j,k))

              sterm(i,j,k,this%YMOMENTUM) = MomentumSourcesY( &
                0., & !ZMOMENTUM for 2D is 0
                cvar(i,j,k,this%XMOMENTUM), &
                pvar(i,j,k,this%XVELOCITY), &
                pvar(i,j,k,this%YVELOCITY), &
                0., & !ZVELOCITY for 2D is 0
                pvar(i,j,k,this%PRESSURE), &
                Mesh%cxyx%bcenter(i,j,k), &
                Mesh%cyxy%bcenter(i,j,k), &
                Mesh%czyz%bcenter(i,j,k), &
                Mesh%cyzy%bcenter(i,j,k))

      !       CALL CalcGeometricalSources(cvar(i,j,k,this%XMOMENTUM),   &
      !                                    cvar(i,j,k,this%YMOMENTUM),   &
      !                                    pvar(i,j,k,this%XVELOCITY),   &
      !                                    pvar(i,j,k,this%YVELOCITY),   &
      !                                    pvar(i,j,k,this%PRESSURE),     &
      !                                    Mesh%cxyx%bcenter(i,j,k),     &
      !                                    Mesh%cyxy%bcenter(i,j,k),     &
      !                                    Mesh%czxz%bcenter(i,j,k),     &
      !                                    Mesh%czyz%bcenter(i,j,k),     &
      !                                    sterm(i,j,k,this%DENSITY),    &
      !                                    sterm(i,j,k,this%XMOMENTUM),  &
      !                                    sterm(i,j,k,this%YMOMENTUM)   &
      !                                    )
           END DO
        END DO
      END DO
      ! reset ghost cell data
      sterm(Mesh%IGMIN:Mesh%IMIN-Mesh%ip1,:,:,:) = 0.0
      sterm(Mesh%IMAX+Mesh%ip1:Mesh%IGMAX,:,:,:) = 0.0
      sterm(:,Mesh%JGMIN:Mesh%JMIN-Mesh%jp1,:,:) = 0.0
      sterm(:,Mesh%JMAX+Mesh%jp1:Mesh%JGMAX,:,:) = 0.0
      sterm(:,:,Mesh%KGMIN:Mesh%KMIN-Mesh%kp1,:) = 0.0
      sterm(:,:,Mesh%KMAX+Mesh%kp1:Mesh%KGMAX,:) = 0.0
    END IF
  END SUBROUTINE GeometricalSources_center

  !TODO GENAUSO WIE IN GEOMETRICALSOURCES_CENTER, cxyx,... not at all touched!
!  PURE SUBROUTINE GeometricalSources_faces(this,Mesh,prim,cons,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler2D), INTENT(IN) :: this
!    CLASS(mesh_base),       INTENT(IN) :: Mesh
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,6,this%VNUM) &
!                                       :: prim,cons
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                       :: sterm
!    !------------------------------------------------------------------------!
!    INTEGER                            :: i,j,k
!    !------------------------------------------------------------------------!
!    DO k=Mesh%KGMIN,Mesh%KGMAX
!      DO j=Mesh%JGMIN,Mesh%JGMAX
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          ! no geometrical density or energy sources
!          sterm(i,j,k,this%DENSITY)   = 0.
!          sterm(i,j,k,this%ENERGY)    = 0.
!          ! momentum sources (sum up corner values, don't use SUM function,
!          ! because it prevents COLLAPSING and causes poor vectorization
!          sterm(i,j,k,this%XMOMENTUM) = MomentumSourcesX( &
!              cons(i,j,k,1,this%YMOMENTUM),cons(i,j,k,1,this%ZMOMENTUM), &
!              prim(i,j,k,1,this%XVELOCITY),prim(i,j,k,1,this%YVELOCITY),&
!              prim(i,j,k,1,this%ZVELOCITY),prim(i,j,k,1,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,1),Mesh%cyxy%faces(i,j,k,1),Mesh%czxz%faces(i,j,k,1)) &
!            + MomentumSourcesX(&
!              cons(i,j,k,2,this%YMOMENTUM),cons(i,j,k,2,this%ZMOMENTUM), &
!              prim(i,j,k,2,this%XVELOCITY),prim(i,j,k,2,this%YVELOCITY), &
!              prim(i,j,k,2,this%ZVELOCITY),prim(i,j,k,2,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,2),Mesh%cyxy%faces(i,j,k,2),Mesh%czxz%faces(i,j,k,2)) &
!            + MomentumSourcesX(&
!              cons(i,j,k,3,this%YMOMENTUM),cons(i,j,k,3,this%ZMOMENTUM), &
!              prim(i,j,k,3,this%XVELOCITY),prim(i,j,k,3,this%YVELOCITY), &
!              prim(i,j,k,3,this%ZVELOCITY),prim(i,j,k,3,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,3),Mesh%cyxy%faces(i,j,k,3),Mesh%czxz%faces(i,j,k,3)) &
!            + MomentumSourcesX(&
!              cons(i,j,k,4,this%YMOMENTUM),cons(i,j,k,4,this%ZMOMENTUM), &
!              prim(i,j,k,4,this%XVELOCITY),prim(i,j,k,4,this%YVELOCITY), &
!              prim(i,j,k,4,this%ZVELOCITY),prim(i,j,k,4,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,4),Mesh%cyxy%faces(i,j,k,4),Mesh%czxz%faces(i,j,k,4)) &
!            + MomentumSourcesX(&
!              cons(i,j,k,5,this%YMOMENTUM),cons(i,j,k,5,this%ZMOMENTUM), &
!              prim(i,j,k,5,this%XVELOCITY),prim(i,j,k,5,this%YVELOCITY), &
!              prim(i,j,k,5,this%ZVELOCITY),prim(i,j,k,5,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,5),Mesh%cyxy%faces(i,j,k,5),Mesh%czxz%faces(i,j,k,5)) &
!            + MomentumSourcesX(&
!              cons(i,j,k,6,this%YMOMENTUM),cons(i,j,k,6,this%ZMOMENTUM), &
!              prim(i,j,k,6,this%XVELOCITY),prim(i,j,k,6,this%YVELOCITY), &
!              prim(i,j,k,6,this%ZVELOCITY),prim(i,j,k,6,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,6),Mesh%cyxy%faces(i,j,k,6),Mesh%czxz%faces(i,j,k,6))
!
!          sterm(i,j,k,this%YMOMENTUM) = &
!              MomentumSourcesY(&
!              cons(i,j,k,1,this%ZMOMENTUM),cons(i,j,k,1,this%XMOMENTUM), &
!              prim(i,j,k,1,this%XVELOCITY),prim(i,j,k,1,this%YVELOCITY), &
!              prim(i,j,k,1,this%ZVELOCITY),prim(i,j,k,1,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,1),Mesh%cyxy%faces(i,j,k,1),Mesh%czyz%faces(i,j,k,1)) &
!            + MomentumSourcesY(&
!              cons(i,j,k,2,this%ZMOMENTUM),cons(i,j,k,2,this%XMOMENTUM), &
!              prim(i,j,k,2,this%XVELOCITY),prim(i,j,k,2,this%YVELOCITY), &
!              prim(i,j,k,2,this%ZVELOCITY),prim(i,j,k,2,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,2),Mesh%cyxy%faces(i,j,k,2),Mesh%czyz%faces(i,j,k,2)) &
!            + MomentumSourcesY(&
!              cons(i,j,k,3,this%ZMOMENTUM),cons(i,j,k,3,this%XMOMENTUM), &
!              prim(i,j,k,3,this%XVELOCITY),prim(i,j,k,3,this%YVELOCITY), &
!              prim(i,j,k,3,this%ZVELOCITY),prim(i,j,k,3,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,3),Mesh%cyxy%faces(i,j,k,3),Mesh%czyz%faces(i,j,k,3)) &
!            + MomentumSourcesY(&
!              cons(i,j,k,4,this%ZMOMENTUM),cons(i,j,k,4,this%XMOMENTUM), &
!              prim(i,j,k,4,this%XVELOCITY),prim(i,j,k,4,this%YVELOCITY), &
!              prim(i,j,k,4,this%ZVELOCITY),prim(i,j,k,4,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,4),Mesh%cyxy%faces(i,j,k,4),Mesh%czyz%faces(i,j,k,4)) &
!            + MomentumSourcesY(&
!              cons(i,j,k,5,this%ZMOMENTUM),cons(i,j,k,5,this%XMOMENTUM), &
!              prim(i,j,k,5,this%XVELOCITY),prim(i,j,k,5,this%YVELOCITY), &
!              prim(i,j,k,5,this%ZVELOCITY),prim(i,j,k,5,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,5),Mesh%cyxy%faces(i,j,k,5),Mesh%czyz%faces(i,j,k,5)) &
!            + MomentumSourcesY(&
!              cons(i,j,k,6,this%ZMOMENTUM),cons(i,j,k,6,this%XMOMENTUM), &
!              prim(i,j,k,6,this%XVELOCITY),prim(i,j,k,6,this%YVELOCITY), &
!              prim(i,j,k,6,this%ZVELOCITY),prim(i,j,k,6,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,6),Mesh%cyxy%faces(i,j,k,6),Mesh%czyz%faces(i,j,k,6))
!
!          sterm(i,j,k,this%ZMOMENTUM) = &
!              MomentumSourcesZ(&
!              cons(i,j,k,1,this%XMOMENTUM),cons(i,j,k,1,this%YMOMENTUM), &
!              prim(i,j,k,1,this%XVELOCITY),prim(i,j,k,1,this%YVELOCITY), &
!              prim(i,j,k,1,this%ZVELOCITY),prim(i,j,k,1,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,1),Mesh%cyxy%faces(i,j,k,1),Mesh%czyz%faces(i,j,k,1)) &
!            + MomentumSourcesZ(&
!              cons(i,j,k,2,this%XMOMENTUM),cons(i,j,k,2,this%YMOMENTUM), &
!              prim(i,j,k,2,this%XVELOCITY),prim(i,j,k,2,this%YVELOCITY), &
!              prim(i,j,k,2,this%ZVELOCITY),prim(i,j,k,2,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,2),Mesh%cyxy%faces(i,j,k,2),Mesh%czyz%faces(i,j,k,2)) &
!            + MomentumSourcesZ(&
!              cons(i,j,k,3,this%XMOMENTUM),cons(i,j,k,3,this%YMOMENTUM), &
!              prim(i,j,k,3,this%XVELOCITY),prim(i,j,k,3,this%YVELOCITY), &
!              prim(i,j,k,3,this%ZVELOCITY),prim(i,j,k,3,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,3),Mesh%cyxy%faces(i,j,k,3),Mesh%czyz%faces(i,j,k,3)) &
!            + MomentumSourcesZ(&
!              cons(i,j,k,4,this%XMOMENTUM),cons(i,j,k,4,this%YMOMENTUM), &
!              prim(i,j,k,4,this%XVELOCITY),prim(i,j,k,4,this%YVELOCITY), &
!              prim(i,j,k,4,this%ZVELOCITY),prim(i,j,k,4,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,4),Mesh%cyxy%faces(i,j,k,4),Mesh%czyz%faces(i,j,k,4)) &
!            + MomentumSourcesZ(&
!              cons(i,j,k,5,this%XMOMENTUM),cons(i,j,k,5,this%YMOMENTUM), &
!              prim(i,j,k,5,this%XVELOCITY),prim(i,j,k,5,this%YVELOCITY), &
!              prim(i,j,k,5,this%ZVELOCITY),prim(i,j,k,5,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,5),Mesh%cyxy%faces(i,j,k,5),Mesh%czyz%faces(i,j,k,5)) &
!            + MomentumSourcesZ(&
!              cons(i,j,k,6,this%XMOMENTUM),cons(i,j,k,6,this%YMOMENTUM), &
!              prim(i,j,k,6,this%XVELOCITY),prim(i,j,k,6,this%YVELOCITY), &
!              prim(i,j,k,6,this%ZVELOCITY),prim(i,j,k,6,this%PRESSURE), &
!              Mesh%cxyx%faces(i,j,k,6),Mesh%cyxy%faces(i,j,k,6),Mesh%czyz%faces(i,j,k,6))
!        END DO
!      END DO
!    END DO
!  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS) &
                                       :: accel
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar,cvar
    REAL,                   INTENT(OUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: sterm
    !------------------------------------------------------------------------!
    INTEGER                            :: i,j,k
    !------------------------------------------------------------------------!
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          sterm(i,j,k,this%DENSITY)   = 0.
          sterm(i,j,k,this%XMOMENTUM) = pvar(i,j,k,this%DENSITY) * accel(i,j,k,1)
          sterm(i,j,k,this%YMOMENTUM) = pvar(i,j,k,this%DENSITY) * accel(i,j,k,2)
          sterm(i,j,k,this%ENERGY)    = &
               cvar(i,j,k,this%XMOMENTUM) * accel(i,j,k,1) + &
               cvar(i,j,k,this%YMOMENTUM) * accel(i,j,k,2)
        END DO
      END DO
    END DO
  END SUBROUTINE ExternalSources

  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(INOUT) :: this
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
   CALL this%ViscositySources_euler2D(Mesh,pvar,btxx,btxy,btyy,sterm)
 
   !compute scalar product of v and tau (x-component)
   this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btxy(:,:,:) 
 !                   + pvar(:,:,:,this%ZVELOCITY)*btxz(:,:,:)

   !compute scalar product of v and tau (y-component)
   this%tmp1(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxy(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btyy(:,:,:) 
  !                  + pvar(:,:,:,this%ZVELOCITY)*btyz(:,:,:)

   !compute scalar product of v and tau (z-component)
   !this%tmp2(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxz(:,:,:) &
   !                 + pvar(:,:,:,this%YVELOCITY)*btyz(:,:,:) &
   !                 + pvar(:,:,:,this%ZVELOCITY)*btzz(:,:,:)
   ! compute vector divergence of scalar product v and tau
   CALL Mesh%Divergence(this%tmp(:,:,:),this%tmp1(:,:,:), &
        sterm(:,:,:,this%ENERGY))
 END SUBROUTINE ViscositySources




  ! identical to isothermal case
  PURE SUBROUTINE ViscositySources_euler2D(this,Mesh,pvar,btxx,btxy,btyy,sterm)
    IMPLICIT NONE
   !------------------------------------------------------------------------!
    CLASS(Physics_euler2D),INTENT(IN)  :: this
    CLASS(Mesh_base),INTENT(IN)        :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: &
          pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: &
                    btxx,btxy,btyy
   !------------------------------------------------------------------------!
   !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,btxx,btxy,btyy
    INTENT(OUT)       :: sterm
   !------------------------------------------------------------------------!
   ! mean values of stress tensor components across the cell interfaces

   ! viscosity source terms
    sterm(:,:,:,this%DENSITY) = 0.0 

   ! compute viscous momentum sources
   ! divergence of stress tensor with symmetry btyx=btxy
    CALL Mesh%Divergence(btxx,btxy,btxy,btyy,sterm(:,:,:,this%XMOMENTUM), &
                         sterm(:,:,:,this%YMOMENTUM))
  END SUBROUTINE ViscositySources_euler2D

  ! identical to isothermal case. 
  PURE SUBROUTINE CalcStresses_euler(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_euler2D), INTENT(INOUT) :: this
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
    CALL Mesh%Divergence(pvar(:,:,:,this%XVELOCITY),pvar(:,:,:,this%YVELOCITY),this%tmp(:,:,:))
    this%tmp(:,:,:) = bulkvis(:,:,:)*this%tmp(:,:,:)

!NEC$ OUTERLOOP_UNROLL(8)
  DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
    DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!NEC$ IVDEP
       DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
          ! compute the diagonal elements of the stress tensor
          btxx(i,j,k) = dynvis(i,j,k) * &
                ((pvar(i+1,j,k,this%XVELOCITY) - pvar(i-1,j,k,this%XVELOCITY)) / Mesh%dlx(i,j,k) &
               + 2.0 * Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) ) &
 !              + 2.0 * Mesh%cxzx%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) ) &
               + this%tmp(i,j,k)

          btyy(i,j,k) = dynvis(i,j,k) * &
               ( (pvar(i,j+1,k,this%YVELOCITY) - pvar(i,j-1,k,this%YVELOCITY)) / Mesh%dly(i,j,k) &
               + 2.0 * Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) ) &
!               + 2.0 * Mesh%cyzy%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) ) &
               + this%tmp(i,j,k)

!          btzz(i,j,k) = dynvis(i,j,k) * &
!               ( (pvar(i,j,k+1,this%ZVELOCITY) - pvar(i,j,k-1,this%ZVELOCITY)) / Mesh%dlz(i,j,k) &
!               + 2.0 * Mesh%czxz%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
!               + 2.0 * Mesh%czyz%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) ) &
!               + this%tmp(i,j,k)

          ! compute the off-diagonal elements (no bulk viscosity)
          btxy(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i+1,j,k,this%YVELOCITY) - pvar(i-1,j,k,this%YVELOCITY)) / Mesh%dlx(i,j,k) &
               + (pvar(i,j+1,k,this%XVELOCITY) - pvar(i,j-1,k,this%XVELOCITY)) / Mesh%dly(i,j,k) ) &
               - Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               - Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

!          btxz(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
!               ( (pvar(i+1,j,k,this%ZVELOCITY) - pvar(i-1,j,k,this%ZVELOCITY)) / Mesh%dlx(i,j,k) &
!               + (pvar(i,j,k+1,this%XVELOCITY) - pvar(i,j,k-1,this%XVELOCITY)) / Mesh%dlz(i,j,k) ) &
!               - Mesh%czxz%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) &
!               - Mesh%cxzx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) )

!          btyz(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
!               ( (pvar(i,j,k+1,this%YVELOCITY) - pvar(i,j,k-1,this%YVELOCITY)) / Mesh%dlz(i,j,k) &
!               + (pvar(i,j+1,k,this%ZVELOCITY) - pvar(i,j-1,k,this%ZVELOCITY)) / Mesh%dly(i,j,k) ) &
!               - Mesh%czyz%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) &
!               - Mesh%cyzy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

       END DO
    END DO
  END DO
  END SUBROUTINE CalcStresses_euler



!  !> Calculate viscous forces
!  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btyy,btyz,btzz,btzx,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler2D), INTENT(INOUT) :: this
!    CLASS(mesh_base),       INTENT(IN)    :: Mesh
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                          :: pvar
!    REAL,                   INTENT(IN), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
!                                          :: btxx,btxy,btyy,btyz,btzz,btzx
!    REAL,                   INTENT(OUT), &
!      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
!                                          :: sterm
!    !------------------------------------------------------------------------!
!    CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!
!    !compute scalar product of v and tau (x-component)
!    this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:) &
!                    + pvar(:,:,:,this%YVELOCITY)*btxy(:,:,:) &
!                    + pvar(:,:,:,this%ZVELOCITY)*btzx(:,:,:)
!
!    !compute scalar product of v and tau (y-component)
!    this%tmp1(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxy(:,:,:) &
!                     + pvar(:,:,:,this%YVELOCITY)*btyy(:,:,:) &
!                     + pvar(:,:,:,this%ZVELOCITY)*btyz(:,:,:)
!
!    !compute scalar product of v and tau (z-component)
!    this%tmp2(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btzx(:,:,:) &
!                     + pvar(:,:,:,this%YVELOCITY)*btyz(:,:,:) &
!                     + pvar(:,:,:,this%ZVELOCITY)*btzz(:,:,:)
!    ! compute vector divergence of scalar product v and tau
!    CALL Divergence(Mesh,this%tmp(:,:,:),this%tmp1(:,:,:), &
!      this%tmp2(:,:,:),sterm(:,:,:,this%ENERGY))
!  END SUBROUTINE ViscositySources


  !> Convert to from conservative to primitive variables at cell-centers
  PURE SUBROUTINE Convert2Primitive_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(IN)  :: cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL Cons2Prim(this%gamma, &
                   cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                   cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                   cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                   cvar(i1:i2,j1:j2,k1:k2,this%ENERGY),    &
                   pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                   pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                   pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                   pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_centsub

  !> Convert to from conservative to primitive variables at cell-faces
  PURE SUBROUTINE Convert2Primitive_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim(this%gamma, &
                   cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY),   &
                   cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY),    &
                   prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY),   &
                   prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE))
  END SUBROUTINE Convert2Primitive_facesub

  !> Convert to from primitve to conservative variables at cell-centers
  PURE SUBROUTINE Convert2Conservative_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(IN)  :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons(this%gamma, &
         pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
         pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
         pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
         pvar(i1:i2,j1:j2,k1:k2,this%PRESSURE),  &
         cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
         cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
         cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
         cvar(i1:i2,j1:j2,k1:k2,this%ENERGY))
  END SUBROUTINE Convert2Conservative_centsub

  !> Convert to from primitve to conservative variables at cell-faces
  PURE SUBROUTINE Convert2Conservative_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    CLASS(mesh_base),       INTENT(IN)  :: Mesh
    INTEGER,                INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)  :: prim
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    !\todo Is direct call ok?
    CALL Prim2Cons(this%gamma, &
         prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY),   &
         prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
         prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
         prim(i1:i2,j1:j2,k1:k2,:,this%PRESSURE),  &
         cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY),   &
         cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
         cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
         cons(i1:i2,j1:j2,k1:k2,:,this%ENERGY))
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
  PURE SUBROUTINE AddBackgroundVelocity(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX), &
                            INTENT(IN)    :: w
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! ATTENTION: don't change the order; on the RHS of the first
          !            assignment there must be the old momentum
          cvar(i,j,k,this%ENERGY) = cvar(i,j,k,this%ENERGY) &
                                + w(i)*(cvar(i,j,k,this%YMOMENTUM) &
                                + 0.5*cvar(i,j,k,this%DENSITY)*w(i))
          pvar(i,j,k,this%YVELOCITY) = pvar(i,j,k,this%YVELOCITY) + w(i)
          cvar(i,j,k,this%YMOMENTUM) = cvar(i,j,k,this%YMOMENTUM) &
                                   + cvar(i,j,k,this%DENSITY)*w(i)
        END DO
      END DO
    END DO
  END SUBROUTINE AddBackgroundVelocity


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
  PURE SUBROUTINE SubtractBackgroundVelocity(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX), &
                            INTENT(IN)    :: w
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! ATTENTION: don't change the order; on the RHS of the first
          !            assignment there must be the old momentum
          cvar(i,j,k,this%ENERGY) = cvar(i,j,k,this%ENERGY) &
                                - w(i)*(cvar(i,j,k,this%YMOMENTUM) &
                                - 0.5*cvar(i,j,k,this%DENSITY)*w(i))
          pvar(i,j,k,this%YVELOCITY) = pvar(i,j,k,this%YVELOCITY) - w(i)
          cvar(i,j,k,this%YMOMENTUM) = cvar(i,j,k,this%YMOMENTUM) &
                                   - cvar(i,j,k,this%DENSITY)*w(i)
        END DO
      END DO
    END DO
  END SUBROUTINE SubtractBackgroundVelocity

  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(IN)  :: this
    LOGICAL, DIMENSION(this%VNUM), &
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
    CLASS(physics_euler2D), INTENT(IN) :: this
    LOGICAL, DIMENSION(this%VNUM), &
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

  PURE SUBROUTINE UpdateSoundSpeed_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
    !------------------------------------------------------------------------!
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          this%bccsound(i,j,k) = GetSoundSpeed( &
            this%gamma, &
            pvar(i,j,k,this%DENSITY), &
            pvar(i,j,k,this%PRESSURE))
        END DO
      END DO
    END DO
  END SUBROUTINE UpdateSoundSpeed_center

  PURE SUBROUTINE UpdateSoundSpeed_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                            INTENT(IN)    :: prim
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k,l
    !------------------------------------------------------------------------!
    DO l=1,Mesh%NFACES
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            this%fcsound(i,j,k,l) = GetSoundSpeed(&
              this%gamma,&
              prim(i,j,k,l,this%DENSITY),&
              prim(i,j,k,l,this%PRESSURE))
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE UpdateSoundSpeed_faces

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

  !> \todo NOT VERIFIED
  ELEMENTAL SUBROUTINE SetFlux(rho,v,P,m1,m2,E,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho,v,P,m1,m2,E
    REAL, INTENT(OUT) :: f1, f2, f3, f4
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + P
    f3 = m2*v
    f4 = (E+P)*v
  END SUBROUTINE SetFlux

  !> \todo NOT VERIFIED
  ELEMENTAL SUBROUTINE Cons2Prim(gamma,rho_in,mu,mv,E,rho_out,u,v,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,mu,mv,E
    REAL, INTENT(OUT) :: rho_out,u,v,P
    !------------------------------------------------------------------------!
    REAL              :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
    P = (gamma-1.)*(E - 0.5 * inv_rho * (mu*mu+mv*mv))
  END SUBROUTINE Cons2Prim

  !> \todo NOT VERIFIED
  ELEMENTAL SUBROUTINE Prim2Cons(gamma,rho_in,u,v,P,rho_out,mu,mv,E)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gamma,rho_in,u,v,P
    REAL, INTENT(OUT) :: rho_out,mu,mv,E
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    E = P/(gamma-1.) + 0.5 * rho_in * (u*u+v*v)
  END SUBROUTINE Prim2Cons

  !> Similar to isothermal version
!  ELEMENTAL SUBROUTINE CalcGeometricalSources(this,mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz,srho,smx,smy,smz)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_euler2D), INTENT(IN) :: this
!    REAL, INTENT(IN)  :: mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz
!    REAL, INTENT(OUT) :: srho, smx, smy, smz
!    !------------------------------------------------------------------------!
!    srho =  0.
!    smx  = -my * (cxyx * vx - cyxy * vy) + mz * (czxz * vz - cxzx * vx) + (cyxy + czxz) * P
!    smy  =  mx * (cxyx * vx - cyxy * vy) + mz * (czyz * vz - cyzy * vy) + (cxyx + czyz) * P
!    smz  =  mx * (cxzx * vx - czxz * vz) + my * (cyzy * vy - czyz * vz) + (cxzx + cyzy) * P
!  END SUBROUTINE CalcGeometricalSources

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2D), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE physics_euler2D_mod
