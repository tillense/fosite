!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sedov2d.f90                                                       #
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
!> 2D Sedov explosion
!! \author Tobias Illenseer
!!
!! References:
!! [1] Sedov, L. I.: Unsteady motions of compressible fluids,
!!     J. Appl. Math. Mech. 9 (1945)
!! [2] Sedov, L. I.: Similarity and Dimensional Methods in Mechanics
!!     Academic Press Ltd., New York (1959)
!----------------------------------------------------------------------------!
PROGRAM sedov2d
  USE fosite_mod
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
  REAL, PARAMETER    :: THETA0 = 0.
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
  INTEGER, PARAMETER :: ONUM = 1          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'sedov2d_'
  !--------------------------------------------------------------------------!
  TYPE(fosite)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config, MGEO)
!    CALL PrintDict(config)
  CALL Sim%Setup()

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,MGEO)

  CALL Sim%Run()
  TAP_CHECK(.TRUE.,"Simulation finished")

  TAP_DONE

CONTAINS

 SUBROUTINE MakeConfig(Sim, config,MGEO)
    USE functions, ONLY : Asinh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: MGEO
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2
    CHARACTER(LEN=16) :: str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: MGEO
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 =-RMAX
       x2 = RMAX
       y1 =-RMAX 
       y2 = RMAX
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
   CASE(SINHPOLAR)
       x1 = Asinh(RMIN/GPAR)
       x2 = Asinh(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2*PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(BIANGLESPHERICAL)
       x1 = 1.E-3
       x2 = PI-1.E-3
       y1 = 0.0
       y2 = 2.*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = REFLECTING
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
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
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.3, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-13, &
           "tol_rel"  / 5e-2, &
           "tol_abs"  / (/0.0,1e-3,1e-3,0.0/), &
           "maxiter"  / 1000000)

    WRITE(str,"(I4)") MGEO

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME) // TRIM(ADJUSTL(str))), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc,MGEO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: Mesh
    CLASS(physics_base),INTENT(IN) :: Physics
    CLASS(timedisc_base),INTENT(INOUT):: Timedisc
    INTEGER           :: MGEO
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: n
    REAL              :: P1
    !------------------------------------------------------------------------!
    INTENT(IN)        :: MGEO
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
    SELECT CASE(MGEO)
    CASE(BIANGLESPHERICAL)
      WHERE ((Mesh%bcenter(:,:,1).LE.(THETA0+R0/GPAR)))
        Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
      ELSEWHERE
       ! in front of the shock front (ambient medium)
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
      END WHERE
    CASE DEFAULT
      WHERE (Mesh%radius%bcenter(:,:).LE.R0)
        ! behind the shock front
        Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
      ELSEWHERE
        ! in front of the shock front (ambient medium)
        Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
      END WHERE
    END SELECT

    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 2D Sedov explosion")

  END SUBROUTINE InitData


END PROGRAM sedov2d
