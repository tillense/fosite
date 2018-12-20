!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
!! [1] Sedov, L. I.: Unsteady motions of compressible fluids,
!!     J. Appl. Math. Mech. 9 (1945)
!! [2] Sedov, L. I.: Similarity and Dimensional Methods in Mechanics
!!     Academic Press Ltd., New York (1959)
!! [3] Padmanabhan, T.:Theoretical Astrophysics, Vol. I: Astrophysical
!!     Processes, Cambridge University Press (2000), Chapter 8.12
!----------------------------------------------------------------------------!
PROGRAM sedov3d
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
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
  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry
!!$  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES  = 100        ! x-resolution
  INTEGER, PARAMETER :: YRES  = 1          ! y-resolution
  REAL, PARAMETER    :: RMAX  = 0.4        ! outer radius of comput. domain
  REAL, PARAMETER    :: GPAR  = 0.2        ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'sedov3d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Simulation finished")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh,Acosh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x1 = 0.0
       x2 = RMAX
       y1 = 0.0
       y2 = PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(CYLINDRICAL)
       x1 = -RMAX
       x2 = RMAX
       y1 = 0.0
       y2 = RMAX
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = ABSORBING
    CASE(OBLATE_SPHEROIDAL)
       x1 = 0.0
       x2 = Acosh(RMAX/GPAR)
       y1 = -0.5*PI
       y2 = 0.5*PI
       bc(WEST)  = FOLDED
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(TANCYLINDRICAL)
       x1 = ATAN(-RMAX/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0
       y2 = RMAX
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = ABSORBING
    CASE(SINHSPHERICAL)
       x1 = 0.0
       x2 = Asinh(RMAX/GPAR)
       y1 = 0.0
       y2 = PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","geometry not supported for 3D Sedov explosion")
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
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA)         ! ratio of specific heats        !

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
           "cfl"      / 0.4, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-13, &
           "tol_rel"  / 0.01, &
           "tol_abs"  / (/1e-5,1e-5,1e-5,1e-5,1e-5/), &
           "maxiter"  / 1000000)

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "decimals" / 6, "cartcoords" / 1, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)
 
    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
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
    REAL              :: P1
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! peak pressure
    n  = 3 ! 3 for 3D
    P1 = 3.*(Physics%gamma-1.0) * E1 / ((n + 1) * PI * R0**n)

    ! uniform density
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
    ! vanishing initial velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.
    ! pressure
    WHERE (Mesh%radius%bcenter(:,:).LE.R0)
       ! behind the shock front
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
    ELSEWHERE
       ! in front of the shock front (ambient medium)
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P0
    END WHERE
    
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 3D Sedov explosion")

  END SUBROUTINE InitData

END PROGRAM sedov3d
