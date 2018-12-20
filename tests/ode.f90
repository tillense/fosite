!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sedov2d.f90                                                       #
!#                                                                           #
!# Copyright (C) 2006-2012, 2013                                             #
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
!> ode test using a 2D sedov explosion
!! \author Tobias Illenseer
!! \author Björn Sperling

!! References:
!! [1] Sedov, L. I.: Unsteady motions of compressible fluids,
!!     J. Appl. Math. Mech. 9 (1945)
!! [2] Sedov, L. I.: Similarity and Dimensional Methods in Mechanics
!!     Academic Press Ltd., New York (1959)
!----------------------------------------------------------------------------!
PROGRAM ode_test
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
  USE common_dict
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
!!$  INTEGER, PARAMETER :: MGEO = POLAR
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES = 100         ! y-resolution
  REAL, PARAMETER    :: RMIN = 1.0E-3      ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX = 0.3         ! outer radius
  REAL, PARAMETER    :: GPAR = 0.2         ! geometry scaling parameter     !
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 5          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'ode'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)      :: Sim
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER    :: NUM = 6 ! last ODE-solver
  INTEGER               :: ode, i
  CHARACTER(LEN=256),DIMENSION(NUM) :: odename
  REAL,DIMENSION(NUM)   :: time
  CHARACTER(LEN=256)    :: PREFIX, buffer
  !--------------------------------------------------------------------------!

  TAP_PLAN(NUM)

DO ode=1,NUM
  CALL InitFosite(Sim)

  write(PREFIX,'((A),(I2.2),(A),(I2.2))') "_ode-",ode,"_geo-",MGEO
  CALL MakeConfig(Sim, Sim%config,ode)

!  CALL PrintDict(Sim%config)
  CALL SetupFosite(Sim)
  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL RunFosite(Sim)

  CALL ComputeRunTime(Sim)
  time(ode) = Sim%run_time
  odename(ode) = GetName(Sim%timedisc)

  TAP_CHECK(.TRUE.,"Simulation finished")
END DO

  CALL Info(Sim, repeat("*",64))
  DO i=1,NUM
    WRITE(buffer,'((a),(I2),(A),(F5.2))') "#", i, ". RUN, ODE solver:  " //TRIM(odename(i)) //", runtime: ", time(i)
    CALL Info(Sim, buffer)
  END DO
  CALL Info(Sim, repeat("*",64))

  CALL CloseFosite(Sim)

  TAP_DONE

CONTAINS

 SUBROUTINE MakeConfig(Sim, config, i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: i
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: i
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 =-RMAX
       x2 = RMAX
       y1 =-RMAX 
       y2 = RMAX
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0 
       y2 = 2*PI       
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
    CASE(SINHPOLAR)
       x1 = RMIN/GPAR
       x1 = LOG(x1+SQRT(1.0+x1*x1))  ! = ASINH(RMIN/GPAR))
       x2 = RMAX/GPAR
       x2 = LOG(x2+SQRT(1.0+x2*x2))  ! = ASINH(RMAX/GPAR))
       y1 = 0.0 
       y2 = 2*PI       
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for 2D Sedov explosion")
    END SELECT
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "gparam"   / GPAR)

    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
    CASE(POLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for 2D Sedov explosion")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    physics => Dict("problem" / EULER2D, &
              "gamma"   / GAMMA)                 ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! time discretization settings
    timedisc => Dict( &
           "method"   / i, &
           "cfl"      / 0.3, &
           "ShowButcherTableau" / 1, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-15, &
           "tol_rel"  / 0.01, &
           "tol_abs"  / (/0.0,1e-3,1e-3,0.0/), &
           "output/error" / 1, &
           "maxiter"  / 1000000)

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME) //TRIM(PREFIX)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: n
    REAL              :: dr,P1
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! compute peak pressure
    n  = 2 ! 2 for 2D
    P1 = 3.*(Physics%gamma-1.0) * E1 / ((n + 1) * PI * R0**n)

    ! uniform density
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
    ! vanishing initial velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    ! pressure
    WHERE ((Mesh%bccart(:,:,1)**2 + Mesh%bccart(:,:,2)**2).LE.R0**2)
       ! behind the shock front
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
    ELSEWHERE
       ! in front of the shock front (ambient medium)
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
    END WHERE
     
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: 2D Sedov explosion")

  END SUBROUTINE InitData


END PROGRAM ode_test
