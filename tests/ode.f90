!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
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
!> \test ode test using a 2D sedov explosion
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
!  INTEGER, PARAMETER :: MGEO = CYLINDRICAL   ! geometry
!!$  INTEGER, PARAMETER :: MGEO = POLAR
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES = 100         ! y-resolution
  INTEGER, PARAMETER :: ZRES = 1           ! z-resolution
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
  CLASS(fosite),ALLOCATABLE      :: Sim
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER    :: NUM = 5 ! last ODE-solver
  INTEGER               :: ode, i
  CHARACTER(LEN=256),DIMENSION(NUM) :: odename
  REAL,DIMENSION(NUM)   :: time
  CHARACTER(LEN=256)    :: PREFIX, buffer
  !--------------------------------------------------------------------------!

  TAP_PLAN(NUM)
DO ode=1,NUM
  ALLOCATE(SIM)
  CALL Sim%InitFosite()

  write(PREFIX,'((A),(I2.2),(A),(I2.2))') "_ode-",ode,"_geo-",MGEO
  CALL MakeConfig(Sim, Sim%config,ode)

!  CALL PrintDict(Sim%config)
  CALL Sim%Setup()
  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL Sim%Run()

  CALL Sim%ComputeRunTime()
  time(ode) = Sim%run_time
  odename(ode) = Sim%Timedisc%GetName()

  TAP_CHECK(.TRUE.,"Simulation finished")
  IF(ode.NE.NUM) THEN
    CALL Sim%Finalize()
    Deallocate(Sim)
  END IF
END DO
  CALL Sim%Info(repeat("*",64))
  DO i=1,NUM
    WRITE(buffer,'((a),(I2),(A),(F5.2))') "#", i, ". RUN, ODE solver:  " //TRIM(odename(i)) //", runtime: ", time(i)
    CALL Sim%Info(buffer)
  END DO
  CALL Sim%Info(repeat("*",64))

  CALL Sim%Finalize()
  Deallocate(Sim)

  TAP_DONE

CONTAINS

 SUBROUTINE MakeConfig(Sim, config, i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: i
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2,z1,z2
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
       z1 = 0.0
       z2 = 0.0
!     CASE(POLAR)
!       x1 = RMIN
!       x2 = RMAX
!       y1 = 0.0
!       y2 = 2*PI   
     CASE(CYLINDRICAL)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0
       y2 = 2*PI
       z1 = 0.0
       z2 = 0.0
!    CASE(LOGPOLAR)
!       x1 = LOG(RMIN/GPAR)
!       x2 = LOG(RMAX/GPAR)
!       y1 = 0.0 
!       y2 = 2*PI       
!    CASE(TANPOLAR)
!       x1 = ATAN(RMIN/GPAR)
!       x2 = ATAN(RMAX/GPAR)
!       y1 = 0.0 
!       y2 = 2*PI       
!    CASE(SINHPOLAR)
!       x1 = RMIN/GPAR
!       x1 = LOG(x1+SQRT(1.0+x1*x1))  ! = ASINH(RMIN/GPAR))
!       x2 = RMAX/GPAR
!       x2 = LOG(x2+SQRT(1.0+x2*x2))  ! = ASINH(RMAX/GPAR))
!       y1 = 0.0 
!       y2 = 2*PI       
    CASE DEFAULT
       CALL Sim%Error("InitProgram","mesh geometry not supported for 2D Sedov explosion")
    END SELECT
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "knum"     / ZRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "zmin"     / z1, &
           "zmax"     / z2, &
           "gparam"   / GPAR)

    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
       bc(BOTTOM)= NO_GRADIENTS
       bc(TOP)   = NO_GRADIENTS
    CASE(CYLINDRICAL)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
       bc(BOTTOM)= NO_GRADIENTS
       bc(TOP)   = NO_GRADIENTS
!     CASE(POLAR)
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
!    CASE(LOGPOLAR)
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
!    CASE(TANPOLAR)
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
!    CASE(SINHPOLAR)
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Sim%Error("InitProgram","mesh geometry not supported for 2D Sedov explosion")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH), &
               "bottomer" / bc(BOTTOM), &
               "topper"   / bc(TOP))

    ! physics settings
    physics => Dict("problem" / EULER, &
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
    datafile => Dict("fileformat" / VTK, &
!    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
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
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_base) :: Physics
    CLASS(Mesh_base)    :: Mesh
    CLASS(Timedisc_base):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: n
    REAL              :: P1
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! isothermal modules are excluded
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler)
      ! peak pressure
      n  = 2 ! 3 for 3D
      P1 = 3.*(phys%gamma - 1.0)*E1 / ((n + 1)*PI*R0**n)
    CLASS DEFAULT
      ! abort
      CALL phys%Error("InitData","physics not supported")
    END SELECT

    ! uniform density
    Timedisc%pvar%data4d(:,:,:,Physics%DENSITY)   = RHO0
    ! vanishing initial velocities
    Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
    ! pressure
    WHERE ((Mesh%bccart(:,:,:,1)**2 + Mesh%bccart(:,:,:,2)**2 + Mesh%bccart(:,:,:,3)**2).LE.R0**2)
       ! behind the shock front
       Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE)  = P1
    ELSEWHERE
       ! in front of the shock front (ambient medium)
       Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE)  = P0
    END WHERE
     
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 2D Sedov explosion")

  END SUBROUTINE InitData


END PROGRAM ode_test
