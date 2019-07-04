!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sedov3d.f90                                                       #
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
!> 3D Sedov explosion
!! \author Tobias Illenseer
!!
!! References:
!! \cite sedov1945
!! \cite sedov1959
!! \cite padmanabhan2000
!----------------------------------------------------------------------------!
PROGRAM sedov3d
  USE fosite_mod
  USE solutions
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 1.0         ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4         ! ratio of specific heats
  REAL, PARAMETER    :: R       = 1.033       ! scaling parameter: 1.0 for 2D, 1.033 for 3D
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHO0    = 1.0         ! ambient density
  REAL, PARAMETER    :: P0      = 1.0E-05     ! ambient pressure
  REAL, PARAMETER    :: E1      = 1.0         ! initial energy input
  ! Spatial with of the initial pulse should be at least 5 cells;
  ! if you wish to compare the results on different grids
  ! R0 should be of the same order
  REAL, PARAMETER    :: R0      = 5.0E-2
  ! mesh settings
!  INTEGER, PARAMETER :: MGEO    = LOGSPHERICAL   ! geometry
  INTEGER, PARAMETER :: MGEO    = SPHERICAL   ! geometry
!  INTEGER, PARAMETER :: MGEO    = CARTESIAN   ! geometry
  INTEGER, PARAMETER :: XRES    = 200          ! x-resolution
  INTEGER, PARAMETER :: YRES    = 6           ! y-resolution
  INTEGER, PARAMETER :: ZRES    = 1           ! z-resolution
  REAL, PARAMETER    :: RMAX    = 2.0         ! outer radius of comput. domain
  REAL, PARAMETER    :: GPAR    = 0.2         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME  = 'sedov3d'
  !-------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE  :: Sim
  REAL, DIMENSION(:), ALLOCATABLE :: sigma
  REAL :: sum_numer, sum_denom
  INTEGER :: n,DEN,VEL,PRE
  LOGICAL :: ok
  !-------------------------------------------------------------------------!
  ALLOCATE(Sim)
  CALL Sim%InitFosite()

#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
TAP_PLAN(4)
#ifdef PARALLEL
  END IF
#endif

  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()

  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  CALL Run(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  ok = .NOT.Sim%aborted
  ! compare with exact solution if requested
  IF (ASSOCIATED(Sim%Timedisc%solution)) THEN
    ALLOCATE(sigma(Sim%Physics%VNUM))
    DO n=1,Sim%Physics%VNUM
      ! use L1 norm to estimate the deviation from the exact solution:
      !   Σ |pvar - pvar_exact| / Σ |pvar_exact|
      sum_numer = SUM(ABS(Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n) &
                        -Sim%Timedisc%solution%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n)))
      sum_denom = SUM(ABS(Sim%Timedisc%solution%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n)))
#ifdef PARALLEL
      IF (Sim%GetRank().GT.0) THEN
        CALL MPI_Reduce(sum_numer,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
        CALL MPI_Reduce(sum_denom,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
      ELSE
        CALL MPI_Reduce(MPI_IN_PLACE,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
        CALL MPI_Reduce(MPI_IN_PLACE,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
#endif
      IF (ABS(sum_denom).LE.2*TINY(sum_denom)) THEN
        sigma(n) = 0.0
      ELSE
        sigma(n) = sum_numer / sum_denom
      END IF
#ifdef PARALLEL
      END IF
#endif
    END DO
  ELSE
    sigma(:) = 0.0
  END IF

  DEN = Sim%Physics%DENSITY
  VEL = Sim%Physics%XVELOCITY
  PRE = Sim%Physics%PRESSURE

!   PRINT *,sigma(:)
#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
TAP_CHECK(ok,"stoptime reached")
! These lines are very long if expanded. So we can't indent it or it will be cropped.
TAP_CHECK_SMALL(sigma(DEN),5.0E-02,"density deviation < 5%")
TAP_CHECK_SMALL(sigma(VEL),4.0E-02,"radial velocity deviation < 4%")
! skip azimuthal velocity deviation, because exact value is 0
TAP_CHECK_SMALL(sigma(PRE),5.0E-02,"pressure deviation < 5%")
TAP_DONE
#ifdef PARALLEL
  END IF
#endif

! #ifdef PARALLEL
!   ALLOCATE(pvar(Sim%Mesh%IMIN:Sim%Mesh%IMAX),pvar_all(Sim%GetNumProcs()*(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1)))
!   ALLOCATE(radius(Sim%Mesh%IMIN:Sim%Mesh%IMAX),radius_all(Sim%GetNumProcs()*(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1)))
!   pvar(:) = Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:SIM%Mesh%IMAX,1,1,1)
!   radius(:)= Sim%Mesh%radius%center(Sim%Mesh%IMIN:Sim%Mesh%IMAX,1,1)
!  
!   ! Compare results with analytical solution. Check only for correct shock velocity
!   ! even if the full analytical solution is implemented
!   IF (ASSOCIATED(Sim%Timedisc%solution)) THEN
!     CALL sedov(GAMMA,E1,RHO0,P0,TSIM,3, SIM%Mesh%radius%bcenter(1:XRES,1,1), &
!                Sim%Timedisc%solution)
!   ELSE
!     sigma(ic,sd) = 0.0
!   END IF
! 
!   
!   !The shock location is where the density gradient is greatest
!   CALL MPI_Gather(pvar,int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE, &
!     pvar_all, int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE,0,MPI_COMM_WORLD,err)
!   CALL MPI_Gather(radius,int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE, &
!     radius_all, int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE,0,MPI_COMM_WORLD,err)
! 
!  IF(Sim%GetRank().EQ.0) THEN 
!   DO i=1,XRES-1
!     pvar_diff(i) = pvar_all(i)-pvar_all(i+1)
!   END DO
!   rshock = radius_all(MAXLOC(pvar_diff,DIM=1))
!   
!   !analytical shock position
!   Rt = R*(E1*TSIM**2/RHO0)**0.2
! 
!   !Check whether analytical solution is within +/- one cell of simulated shock
!   TAP_CHECK((Rt.LT.rshock+Sim%Mesh%dx).AND.(Rt.GT.rshock-Sim%Mesh%dx),"Shock velocity correct")
! END IF
! #else
!   DO i=1,XRES
!       pvar_diff(i) = Sim%Timedisc%pvar%data4d(i,1,1,1)-Sim%Timedisc%pvar%data4d(i+1,1,1,1)
!   END DO
!   Rshock = Sim%Mesh%radius%center(MAXLOC(pvar_diff,DIM=1),1,1)
!   
!   !analytical shock position
!   Rt = R*(E1*TSIM**2/RHO0)**0.2
! 
!   TAP_CHECK((Rt.LT.rshock+Sim%Mesh%dx).AND.(Rt.GT.rshock-Sim%Mesh%dx),"Shock velocity correct")
! #endif
  CALL Sim%Finalize()
!   TAP_DONE
  DEALLOCATE(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh,Acosh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP),POINTER  :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                 :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes
    REAL                    :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 = -0.5
       x2 =  0.5
       y1 = -0.5
       y2 =  0.5
       z1 = -0.5
       z2 =  0.5
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
       bc(BOTTOM)= NO_GRADIENTS
       bc(TOP)   = NO_GRADIENTS
    CASE(SPHERICAL)
       x1 = 0.0
       x2 = RMAX
       y1 = 0.0
       y2 = PI
       z1 = 0.0
       z2 = 0.0 !2*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
       bc(BOTTOM)= PERIODIC
       bc(TOP)   = PERIODIC
    CASE(LOGSPHERICAL)
       x1 = LOG(1e-3/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0
       y2 = PI
       z1 = 0.0
       z2 = 2*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
       bc(BOTTOM)= PERIODIC
       bc(TOP)   = PERIODIC
    CASE(CYLINDRICAL)
       x1 = -RMAX
       x2 = RMAX
       y1 = 0.0
       y2 = 2.0*PI
       z1 = 0.0
       z2 = RMAX
       bc(WEST)  = NO_GRADIENTS    !default: ABSORBING
       bc(EAST)  = NO_GRADIENTS    !default: ABSORBING
       bc(SOUTH) = NO_GRADIENTS    !default: AXIS
       bc(NORTH) = NO_GRADIENTS    !default: ABSORBING
       bc(BOTTOM)= NO_GRADIENTS
       bc(TOP)   = NO_GRADIENTS
!     CASE(OBLATE_SPHEROIDAL)       !TODO anderen Meshs zu 3D konvertieren
!        x1 = 0.0
!        x2 = Acosh(RMAX/GPAR)
!        y1 = -0.5*PI
!        y2 = 0.5*PI
!        bc(WEST)  = FOLDED
!        bc(EAST)  = NO_GRADIENTS
!        bc(SOUTH) = AXIS
!        bc(NORTH) = AXIS
!     CASE(TANCYLINDRICAL)
!        x1 = ATAN(-RMAX/GPAR)
!        x2 = ATAN(RMAX/GPAR)
!        y1 = 0.0
!        y2 = RMAX
!        bc(WEST)  = NO_GRADIENTS
!        bc(EAST)  = NO_GRADIENTS
!        bc(SOUTH) = AXIS
!        bc(NORTH) = NO_GRADIENTS
!     CASE(SINHSPHERICAL)
!        x1 = 0.0
!        x2 = Asinh(RMAX/GPAR)
!        y1 = 0.0
!        y2 = PI
!        bc(WEST)  = REFLECTING
!        bc(EAST)  = NO_GRADIENTS
!        bc(SOUTH) = AXIS
!        bc(NORTH) = AXIS
    CASE DEFAULT
       CALL Sim%Physics%Error("InitProgram","geometry not supported for 3D Sedov explosion")
    END SELECT

    ! mesh settings
    mesh => Dict( &
              "meshtype"    / MIDPOINT, &
              "geometry"    /     MGEO, &
              "inum"        /     XRES, &
              "jnum"        /     YRES, &
              "knum"        /     ZRES, &
              "xmin"        /       x1, &
              "xmax"        /       x2, &
              "ymin"        /       y1, &
              "ymax"        /       y2, &
              "zmin"        /       z1, &
              "zmax"        /       z2, &
              "gparam"      /     GPAR  )

    ! boundary conditions
    boundary => Dict( &
              "western"     / bc(WEST),   &
              "eastern"     / bc(EAST),   &
              "southern"    / bc(SOUTH),  &
              "northern"    / bc(NORTH),  &
              "bottomer"    / bc(BOTTOM), &
              "topper"      / bc(TOP))

    ! physics settings
    physics => DICT( &
              "problem"     / EULER, &
              "gamma"       / GAMMA)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
                "order"       / LINEAR,    &
                "fluxtype"    / KT,        &
                "variables"   / PRIMITIVE, &
                "limiter"     / VANLEER)

    ! time discretization settings
    timedisc => Dict( &
      "method"      / DORMAND_PRINCE,               &
              "order"       / 5,                            &
              "cfl"         / 0.4,                          &
              "stoptime"    / TSIM,                         &
              "dtlimit"     / 1.0E-13,                      &
              "tol_rel"     / 0.01,                         &
!               "tol_abs"     / (/1e-5,1e-5,1e-5,1e-5,1e-5/), &
              "output/solution" / 1, &
              "maxiter"     / 1000000                       )

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"  / VTK,                          &
              "filename"    / (TRIM(ODIR) // TRIM(OFNAME)), &
              "count"       / ONUM)

    config => Dict( &
              "mesh"        / mesh,     &
              "physics"     / physics,  &
              "boundary"    / boundary, &
              "fluxes"      / fluxes,   &
              "timedisc"    / timedisc, &
              "datafile"    / datafile  )
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                             :: n
    REAL                                :: P1
    !------------------------------------------------------------------------!
    ! isothermal modules are excluded
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler)
      SELECT TYPE (pvar => Timedisc%pvar)
      CLASS IS(statevector_euler)
        ! uniform density everywhere
        pvar%density%data1d(:)  = RHO0
        ! vanishing initial velocities
        pvar%velocity%data1d(:) = 0.0
        ! set initial peak pressure P1 inside sphere with radius R0 centered on the origin
        n  = 3 ! 3 for 3D
        P1 = 3.*(phys%gamma - 1.0)*E1 / ((n + 1)*PI*R0**n)
        WHERE (Mesh%radius%bcenter(:,:,:).LE.R0)
          ! behind the shock front
          pvar%pressure%data3d(:,:,:)  = P1
        ELSEWHERE
          ! in front of the shock front (ambient medium)
          pvar%pressure%data3d(:,:,:)  = P0
        END WHERE
      CLASS DEFAULT
        ! abort
        CALL phys%Error("sedov3d:InitData","statevector must be of class euler")
      END SELECT
    CLASS DEFAULT
      ! abort
      CALL phys%Error("sedov3d:InitData","physics not supported")
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 3D Sedov explosion")

  END SUBROUTINE InitData

  SUBROUTINE Run(Mesh,Physics,Timedisc)
    USE sedov_mod
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    TYPE(sedov_typ), ALLOCATABLE :: sedov
    REAL :: T0,vr
    REAL, DIMENSION(:,:,:,:), POINTER :: er
    INTEGER :: n,m,i,j,k
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(Timedisc%solution)) THEN
      ! initialize sedov class
      ! parameters: ndim, gamma, w, rho0, E1, P0
      sedov = sedov_typ(3,GAMMA,0.0,RHO0,E1,P0)
#ifdef PARALLEL
      IF (Sim%GetRank().EQ.0) &
#endif
      CALL sedov%PrintConfiguration()
      ! store initial condition in exact solution array
      Timedisc%solution = Timedisc%pvar

      ! compute initial time using R0 as the initial position of the shock
      T0 = SQRT((E1/RHO0)*(R0/sedov%alpha)**5)

      ! set radial unit vector to compute full 3D radial velocities in any geometry
      ALLOCATE(er(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
      DO n=1,3
        er(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,n) = &
            Mesh%posvec%bcenter(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,n) &
          / Mesh%radius%bcenter(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
      END DO
      DO WHILE((Timedisc%maxiter.LE.0).OR.(Sim%iter.LE.Timedisc%maxiter))
        SELECT TYPE(solution => Timedisc%solution)
        TYPE IS (statevector_euler)
          ! compute exact solution of the spherical Sedov explosion problem
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
              DO i=Mesh%IMIN,Mesh%IMAX
                CALL sedov%ComputeSolution(Timedisc%time+T0,Mesh%radius%bcenter(i,j,k), &
                  solution%density%data3d(i,j,k),vr,solution%pressure%data3d(i,j,k))
                  ! set the velocities
                  m = 1
                  DO n=1,3
                    IF (BTEST(Mesh%VECTOR_COMPONENTS,n-1)) THEN
                      solution%velocity%data4d(i,j,k,m) = vr * er(i,j,k,n)
                      m = m + 1
                    END IF
                  END DO
              END DO
            END DO
          END DO
        END SELECT
        IF(Sim%Step()) EXIT
      END DO
      CALL sedov%Destroy()
    ELSE
      CALL Sim%Run()
    END IF
  END SUBROUTINE Run

END PROGRAM sedov3d
