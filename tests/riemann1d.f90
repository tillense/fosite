!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: riemann1d.f90                                                     #
!#                                                                           #
!# Copyright (C) 2006-2024                                                   #
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
!> 1D Riemann problems
!! \author Tobias Illenseer
!!
!! References:
!! - \cite toro1999 Toro, E. F.: Riemann Solvers and Numerical Methods for Fluid Dynamics,
!!     A Practical Introduction, Springer-Verlag 1999, 2nd ed., Chapter 4.3.3
!! - \cite sod1978 Sod, G. A.: A survey of several finite difference methods for systems of
!!     nonlinear hyperbolic conservation laws, J. Comput. Phys. 27 (1978), 1-31
!!     DOI: 10.1016/0021-9991(78)90023-2
!! - \cite noh1987 noh Noh, W. F.: Errors for calculations of strong shocks using an artificial
!!     viscosity and an artificial heat-flux, J. Comput. Phys. 72 (1987), 78-120
!!     DOI: 10.1016/0021-9991(87)90074-X
!!
!! \test Eight different tests, which show the numerical solution for 1D 
!!       isothermal and non-isothermal Euler equations with piecewise
!!       constant initial conditions.
!!
!! initial conditions, one of
!!  1. Sod problem, see ref. \cite toro1999, \cite sod1978
!!  2. Toro test no. 2, see ref. \cite toro1999
!!  3. Toro test no. 3, see ref. \cite toro1999
!!  4. Toro test no. 4, see ref. \cite toro1999
!!  5. Toro test no. 5, see ref. \cite toro1999
!!  6. Noh problem, see ref. \cite noh1987
!!  7. Isothermal shock tube
!!  8. simple isothermal Riemann problem
!----------------------------------------------------------------------------!
PROGRAM riemann1d
  USE fosite_mod
  USE solutions
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: TESTNUM = 8 ! number of available tests (do not modify)
  CHARACTER(LEN=32), PARAMETER :: TESTSTR(TESTNUM) = (/ &
                              "Sod shock tube                  ", &
                              "Toro test no. 2                 ", &
                              "Toro test no. 3                 ", &
                              "Toro test no. 4                 ", &
                              "Toro test no. 5                 ", &
                              "Noh problem                     ", &
                              "Isothermal shock tube           ", &
                              "Isothermal Noh problem          " /)
  !--------------------------------------------------------------------------!
  ! mesh settings
  CHARACTER(LEN=3), PARAMETER :: TESTDIR  = "all" !"x" ! direction: x,y,z or all
  INTEGER, PARAMETER :: RES  = 400       ! resolution
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 1         ! number of output data sets
  CHARACTER(LEN=8), PARAMETER :: OFNAME(TESTNUM) = (/ & ! output file name
    "sod     ","toro2   ","toro3   ","toro4   ","toro5   ","noh     ", &
    "sodiso  ","nohiso  " /)
  CHARACTER(LEN=256), PARAMETER &        ! output data dir
                     :: ODIR = './'
  ! some global constants
  REAL               :: GAMMA            ! ratio of specific heats
  REAL               :: TSIM             ! simulation time
  REAL               :: CSISO = 0.0      ! isothermal sound speed (test no. 6)
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE     :: Sim
  CHARACTER(LEN=128) :: verbose_tap_output
  INTEGER            :: ic,sd,dir_min,dir_max,n
  REAL, DIMENSION(:,:), ALLOCATABLE :: sigma
  REAL, PARAMETER    :: sigma_tol(TESTNUM) &
                        = (/ 3.0E-3, 1.5E-2, 5.0E-3, 7.0E-3, 4.0E-3, 5.0E-3, 1.0E+0, 1.0E+0 /)
  REAL               :: sum_numer, sum_denom
  !--------------------------------------------------------------------------!

  ! check whether we perform tests in all directions
  SELECT CASE (TRIM(TESTDIR))
  CASE("all")
    dir_min = 1
    dir_max = 3
  CASE("y")
    dir_min = 2
    dir_max = 2
  CASE("z")
    dir_min = 3
    dir_max = 3
  CASE DEFAULT
    ! x-direction is also the default
    dir_min = 1
    dir_max = 1
  END SELECT

  ALLOCATE(Sim,sigma(TESTNUM,dir_min:dir_max))

  ! initialize Fosite
  CALL Sim%InitFosite()

#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
TAP_PLAN(TESTNUM*(dir_max-dir_min+1))
#ifdef PARALLEL
  END IF
#endif

  ! loop over all tests
  tests: DO ic=1,TESTNUM
    ! loop over selected directions
    DO sd=dir_min,dir_max
      ! get configuration for the particular test and direction
      CALL MakeConfig(Sim, Sim%config,ic,sd)
!       CALL PrintDict(config)

      ! setup the simulation
      CALL Sim%Setup()

      ! set initial conditions and run the simulation
      CALL Run(Sim,ic,sd)

      IF (Sim%aborted) THEN
        sigma(ic,sd) = HUGE(1.0)
      ELSE
        sigma(ic,sd) = 0.0 ! default
        sum_numer = 0.0
        sum_denom = 0.0
        IF (ASSOCIATED(Sim%Timedisc%solution)) THEN
          ! use L1 norm to estimate the deviation from the exact solution:
          !   Σ |pvar - pvar_exact| / Σ |pvar_exact|
          DO n=1,Sim%Physics%VNUM
            sum_numer = SUM(ABS(Sim%Timedisc%pvar%data2d(:,n)-Sim%Timedisc%solution%data2d(:,n)), &
                                MASK=Sim%Mesh%without_ghost_zones%mask1d(:))
            sum_denom = SUM(ABS(Sim%Timedisc%solution%data2d(:,n)),MASK=Sim%Mesh%without_ghost_zones%mask1d(:))
          END DO
#ifdef PARALLEL
          CALL MPI_Allreduce(MPI_IN_PLACE,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,MPI_COMM_WORLD,Sim%ierror)
          CALL MPI_Allreduce(MPI_IN_PLACE,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,MPI_COMM_WORLD,Sim%ierror)
#endif
          IF (ABS(sum_denom).GT.4*TINY(sum_denom)) & ! avoid division by zero error
            sigma(ic,sd) = sum_numer / sum_denom
        END IF
      END IF

      ! reset fosite
      IF (.NOT.(ic.EQ.TESTNUM.AND.sd.EQ.dir_max)) THEN
         CALL Sim%Finalize(mpifinalize_=.FALSE.)
         DEALLOCATE(Sim)
         ALLOCATE(Sim)
         CALL Sim%InitFosite()
      END IF
    END DO
  END DO tests

  ! check results
  ! loop over all tests
#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
  DO ic=1,TESTNUM
    ! loop over selected directions
    DO sd=dir_min,dir_max
      WRITE (verbose_tap_output,'(A,ES10.2)') TESTSTR(ic) // " " // ACHAR(119+sd) &
        // "-direction: " // " sigma = ", sigma(ic,sd)
! This line is long if expanded. So we can't indent it or it will be cropped.
TAP_CHECK_SMALL(sigma(ic,sd),sigma_tol(ic),TRIM(verbose_tap_output))
    END DO
  END DO

TAP_DONE
#ifdef PARALLEL
  END IF
#endif

CALL Sim%Finalize()
DEALLOCATE(Sim,sigma)

CONTAINS

  SUBROUTINE MakeConfig(Sim,config,ic,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: ic,dir
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes
    INTEGER           :: xres,yres,zres
    REAL              :: xmax,ymax,zmax
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    INTENT(IN)        :: ic,dir
    !------------------------------------------------------------------------!
    ! mesh settings depends on the direction selected above
    SELECT CASE(DIR)
    CASE(1)
      xmax = 1.0
      ymax = 0.0
      zmax = 0.0
      xres = RES
      yres = 1
      zres = 1
    CASE(2)
      xmax = 0.0
      ymax = 1.0
      zmax = 0.0
      xres = 1
      yres = RES
      zres = 1
    CASE(3)
      xmax = 0.0
      ymax = 0.0
      zmax = 1.0
      xres = 1
      yres = 1
      zres = RES
    CASE DEFAULT
      CALL Sim%Error("riemann1d::InitData","direction must be one of 1,2,3")
    END SELECT

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
               "inum" / xres, &            ! resolution in x and            !
               "jnum" / yres, &            !   y direction                  !
               "knum" / zres, &            !   z direction                  !
               "xmin" / (0.0), &
               "xmax" / xmax, &
               "ymin" / (-0.0), &
               "ymax" / ymax, &
               "zmin" / (0.0), &
               "zmax" / zmax  &
               )
    
    SELECT CASE(ic)
    CASE(1)
       GAMMA  = 1.4
       TSIM   = 0.25
    CASE(2)
       GAMMA  = 1.4
       TSIM   = 0.15
    CASE(3)
       GAMMA  = 1.4
       TSIM   = 0.012
    CASE(4)
       GAMMA  = 1.4
       TSIM   = 0.035
    CASE(5)
       GAMMA  = 1.4
       TSIM   = 0.012
    CASE(6)
       GAMMA  = 5./3.
       TSIM   = 1.0
    CASE(7)
       CSISO  = 1.0
       TSIM   = 0.25
    CASE(8)
       CSISO  = 1.0
       TSIM   = 0.5
    CASE DEFAULT
       CALL Sim%Mesh%Error("riemann1d::InitData", "Test problem number should be 1,2,3,4,5,6,7 or 8")
    END SELECT


    ! physics settings
    IF (CSISO.NE.0.0) THEN
       physics => Dict("problem" / EULER_ISOTHERM, &
                 "cs"      / CSISO)
    ELSE   
       physics => Dict("problem" / EULER, &
                 "gamma"   / GAMMA)         ! ratio of specific heats        !
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict(&
             "fluxtype"  / KT, &
             "order"     / LINEAR, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2 &          ! optional parameter for limiter !
             )
    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.4, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-10, &
!           "output/xmomentum" / 1, &
!           "output/ymomentum" / 1, &
!           "output/zmomentum" / 1, &
!           "output/energy" / 1, &
!           "output/rhs" / 1 &
           "output/solution" / 0, &
           "maxiter"  / 100000 &
           )

    ! compute and output exact solution in non-isothermal simulations
    IF (CSISO.EQ.0.0) CALL SetAttr(timedisc,"output/solution",1)

    ! boundary conditions
    boundary => Dict( &
           "western"  / NO_GRADIENTS, &
           "eastern"  / NO_GRADIENTS, &
           "southern" / NO_GRADIENTS, &
           "northern" / NO_GRADIENTS, &
           "bottomer" / NO_GRADIENTS, &
           "topper"   / NO_GRADIENTS)

    ! initialize data input/output
    datafile => Dict( &
           "fileformat" / GNUPLOT, &
           "filename"   / (TRIM(ODIR) // TRIM(OFNAME(ic)) // "-" // ACHAR(119+dir)), &
           "filecycles" / 0, &
           "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE Run(this,ic,dir)
    USE solutions
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
!     CLASS(marray_compound), POINTER :: pvar
    INTEGER, INTENT(IN)          :: ic,dir
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL    :: x0, rho_l, rho_r, u_l, u_r, p_l, p_r
    REAL, DIMENSION(:), ALLOCATABLE :: x
    !------------------------------------------------------------------------!
    SELECT CASE(ic)
    CASE(1) ! Sod shock tube (Toro's test no. 1)
       x0    = 0.5
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1.0
       p_r   = 0.1
    CASE(2) ! Toro's test no. 2
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = -2.0
       u_r   = 2.0
       p_l   = 0.4
       p_r   = 0.4
    CASE(3) ! Toro's test no. 3
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 0.0
       u_r   = 0.0
       p_l   = 1000.
       p_r   = 0.01
    CASE(4) ! Toro's test no. 4
       x0    = 0.4
       rho_l = 5.99924
       rho_r = 5.99924
       u_l   = 19.5975
       u_r   = -6.19633
       p_l   = 460.894
       p_r   = 46.0950
    CASE(5) ! Toro's test no. 5
       x0    = 0.8
       rho_l = 1.0
       rho_r = 1.0
       u_l   = -19.5975
       u_r   = -19.5975
       p_l   = 1000.0
       p_r   = 0.01
    CASE(6) ! Noh problem
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 1.0
       u_r   = -1.0
       p_l   = 1.0E-05
       p_r   = 1.0E-05
    CASE(7) ! isothermal shock tube
       x0    = 0.5
       rho_l = 1.0
       rho_r = 0.125
       u_l   = 0.0
       u_r   = 0.0
    CASE(8) ! isothermal Noh problem
       x0    = 0.5
       rho_l = 1.0
       rho_r = 1.0
       u_l   = 1.0
       u_r   = -1.0
    CASE DEFAULT
       CALL this%Error("InitData", "Test problem should be 1,2,3,4,5,6,7 or 8")
    END SELECT

    ! initial condition
    ! always set density and first velocity
    SELECT TYPE(pvar => this%Timedisc%pvar)
    CLASS IS(statevector_eulerisotherm)
      ! set all velocities to 0
      pvar%velocity%data1d(:) = 0.0
      WHERE (this%Mesh%bcenter(:,:,:,dir).LT.x0)
        pvar%density%data3d(:,:,:)    = rho_l
        pvar%velocity%data4d(:,:,:,1) = u_l
      ELSEWHERE
        pvar%density%data3d(:,:,:)    = rho_r
        pvar%velocity%data4d(:,:,:,1) = u_r
      END WHERE
    END SELECT
    
    ! set pressure for non-isothermal HD
    SELECT TYPE(pvar => this%Timedisc%pvar)
    TYPE IS(statevector_euler)
      WHERE (this%Mesh%bcenter(:,:,:,dir).LT.x0)
        pvar%pressure%data3d(:,:,:)   = p_l
      ELSEWHERE
        pvar%pressure%data3d(:,:,:)   = p_r
      END WHERE
    END SELECT
   
    CALL this%Physics%Convert2Conservative(this%Timedisc%pvar,this%Timedisc%cvar)
    CALL this%Info(" DATA-----> initial condition: " &
      // "1D Riemann problem: " // TRIM(TESTSTR(ic)))

    IF (ASSOCIATED(this%Timedisc%solution)) THEN
      ! store initial condition in exact solution array
      this%Timedisc%solution = this%Timedisc%pvar

      ! set 1D coordinate array for exact solution
      ALLOCATE(x(SIZE(this%Mesh%cart%data3d(:,2,DIR),1)), &
        SOURCE=this%Mesh%cart%data3d(:,2,DIR))

      DO WHILE((this%Timedisc%maxiter.LE.0).OR.(this%iter.LE.this%Timedisc%maxiter))
        ! compute exact solution of the 1D Riemann problem (only non-isothermal HD)
        ! if output is written during next integration step
        IF (this%Datafile%time-this%Timedisc%time.LT.this%Timedisc%dt) THEN
          SELECT TYPE(phys => this%Physics)
          TYPE IS (physics_euler)
            CALL riemann(x0,phys%gamma,rho_l,u_l,p_l,rho_r,u_r,p_r, &
                        this%Datafile%time,x,this%Timedisc%solution%data2d(:,:))
          END SELECT
        END IF
        IF(this%Step()) EXIT
       END DO
    ELSE
       CALL this%Run()
    END IF
  END SUBROUTINE Run

END PROGRAM riemann1d
