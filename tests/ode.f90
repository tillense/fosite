!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: ode.f90                                                           #
!#                                                                           #
!# Copyright (C) 2006-2021                                                   #
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
!> \test time stepping methods using 2D sedov explosion simulation
!! \author Tobias Illenseer
!! \author Björn Sperling

!! References:
!! \cite sedov1945 Sedov, L. I.: Unsteady motions of compressible fluids,
!!     J. Appl. Math. Mech. 9 (1945)
!! \cite sedov1959  Sedov, L. I.: Similarity and Dimensional Methods in Mechanics
!!     Academic Press Ltd., New York (1959)
!----------------------------------------------------------------------------!
PROGRAM ode_test
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 0.05     ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHO0 = 1.0         ! ambient density
  REAL, PARAMETER    :: P0   = 1.0E-05     ! ambient pressure
  REAL, PARAMETER    :: E1   = 1.0         ! initial energy input
  ! Spatial with of the initial pulse should be at least 5 cells;
  ! if you wish to compare the results on different grids
  ! R0 should be of the same order
  REAL, PARAMETER    :: R0   = 3.0E-2
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry
  INTEGER, PARAMETER :: RES  = 64         ! x- & y-resolution
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 5          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'ode'
  ! ode solvers
  INTEGER, PARAMETER :: NUM_TESTS = 8
  !--------------------------------------------------------------------------!
  CLASS(fosite),ALLOCATABLE      :: Sim
  INTEGER               :: k
  CHARACTER(LEN=80)     :: testinfo(NUM_TESTS)
  !--------------------------------------------------------------------------!

  TAP_PLAN(NUM_TESTS)

  DO k=1,NUM_TESTS
    ALLOCATE(Sim)
    CALL Sim%InitFosite()

    CALL MakeConfig(Sim, Sim%config)

    ! specific test settings
    SELECT CASE(k)
    CASE(1) ! Modified Euler, Order 2
      CALL SetAttr(Sim%config, "/timedisc/method", MODIFIED_EULER)
      CALL SetAttr(Sim%config, "/timedisc/order", 2)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-mdeu2") )
    CASE(2) ! Modified Euler, Order 3
      CALL SetAttr(Sim%config, "/timedisc/method", MODIFIED_EULER)
      CALL SetAttr(Sim%config, "/timedisc/order", 3)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-mdeu3") )
    CASE(3) ! RK-Fehlberg, Order 3
      CALL SetAttr(Sim%config, "/timedisc/method", RK_FEHLBERG)
      CALL SetAttr(Sim%config, "/timedisc/order", 3)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-rkfeh3") )
    CASE(4) ! RK-Fehlberg, Order 5
      CALL SetAttr(Sim%config, "/timedisc/method", RK_FEHLBERG)
      CALL SetAttr(Sim%config, "/timedisc/order", 5)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-rkfeh5") )
    CASE(5) ! Cash-Karp, Order 5
      CALL SetAttr(Sim%config, "/timedisc/method", CASH_KARP)
      CALL SetAttr(Sim%config, "/timedisc/order", 5)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-cshka5") )
    CASE(6) ! Dormand-Prince, Order 5
      CALL SetAttr(Sim%config, "/timedisc/method", DORMAND_PRINCE)
      CALL SetAttr(Sim%config, "/timedisc/order", 5)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-drmpr5") )
    CASE(7) ! SSPRK, Order 5
      CALL SetAttr(Sim%config, "/timedisc/method", SSPRK)
      CALL SetAttr(Sim%config, "/timedisc/order", 3)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-ssprk3") )
    CASE(8) ! SSPRK, Order 5
      CALL SetAttr(Sim%config, "/timedisc/method", SSPRK)
      CALL SetAttr(Sim%config, "/timedisc/order", 5)
      CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "-ssprk5") )
    CASE DEFAULT
      CALL Sim%Error("ode","Unknown test number!")
    END SELECT

  !  CALL PrintDict(Sim%config)
    CALL Sim%Setup()
    ! set initial condition
    CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

    CALL Sim%Run()

    IF (.NOT. Sim%aborted) THEN
      TAP_CHECK(.TRUE.,"Simulation finished")
      CALL Sim%ComputeRunTime()
      WRITE(testinfo(k),'(A,I3,A38,A,I2,A11,F5.2,A2)') "#", k, &
        ". ODE solver: " // ADJUSTL(Sim%Timedisc%GetName()) // ", ", &
        "order", Sim%Timedisc%GetOrder(), ", runtime: ", Sim%run_time, " s"
    ELSE
      WRITE(testinfo(k),'(A)') "#", k, ". ODE solver:  " // ADJUSTL(Sim%Timedisc%GetName()) // "FAILED"
    END IF

    IF(k.LT.NUM_TESTS) THEN
      CALL Sim%Finalize(.FALSE.)
      DEALLOCATE(Sim)
    END IF
  END DO

  ! print test summary
  CALL Sim%Finalize()
  CALL Sim%Info(REPEAT("*",67))
  DO k=1,NUM_TESTS
    CALL Sim%Info(testinfo(k))
  END DO
  CALL Sim%Info(REPEAT("*",67))

  DEALLOCATE(Sim)

  TAP_DONE

CONTAINS

 SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
           "inum"     / RES, &
           "jnum"     / RES, &
           "knum"     / 1, &
           "xmin"     / (-0.3), &
           "xmax"     / 0.3, &
           "ymin"     / (-0.3), &
           "ymax"     / 0.3, &
           "zmin"     / 0.0, &
           "zmax"     / 0.0, &
           "gparam"   / 1.0)

    ! boundary conditions
    boundary => Dict("western" / NO_GRADIENTS, &
               "eastern" / NO_GRADIENTS, &
               "southern" / NO_GRADIENTS, &
               "northern" / NO_GRADIENTS, &
               "bottomer" / NO_GRADIENTS, &
               "topper"   / NO_GRADIENTS)

    ! physics settings
    physics => Dict("problem" / EULER, &
              "gamma"   / GAMMA)                 ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order" / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, & ! overwritten in main program
           "cfl"      / 0.3, &
           "ShowButcherTableau" / 1, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-15, &
           "tol_rel"  / 0.01, &
           "tol_abs"  / (/0.0,1e-3,1e-3,0.0/), &
           "output/error" / 1, &
           "maxiter"  / 1000000)

    ! initialize data input/output
    datafile => Dict("fileformat" / VTK, &
!    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), & ! overwritte in main program
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
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
        n  = 2 ! 2 for 2D
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
        CALL phys%Error("ode:InitData","statevector must be of class euler")
      END SELECT
    CLASS DEFAULT
      ! abort
      CALL phys%Error("ode:InitData","physics not supported")
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 2D Sedov explosion")

  END SUBROUTINE InitData


END PROGRAM ode_test
