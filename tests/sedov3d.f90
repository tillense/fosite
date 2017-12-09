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
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 0.05        ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4         ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHO0    = 1.0         ! ambient density
  REAL, PARAMETER    :: P0      = 1.0E-05     ! ambient pressure
  REAL, PARAMETER    :: E1      = 1.0         ! initial energy input
  ! Spatial with of the initial pulse should be at least 5 cells;
  ! if you wish to compare the results on different grids
  ! R0 should be of the same order
  REAL, PARAMETER    :: R0      = 0.1!3.0E-2
  ! mesh settings
  INTEGER, PARAMETER :: MGEO    = CARTESIAN   ! geometry
  INTEGER, PARAMETER :: XRES    = 40          ! x-resolution
  INTEGER, PARAMETER :: YRES    = 40          ! y-resolution
  INTEGER, PARAMETER :: ZRES    = 40          ! z-resolution
  REAL, PARAMETER    :: RMAX    = 0.4         ! outer radius of comput. domain
  REAL, PARAMETER    :: GPAR    = 0.2         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME  = 'sedov3d'
  !-------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE  :: Sim
  !-------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL Sim%Run()

  DEALLOCATE(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
  TAP_DONE


CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh,Acosh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)          :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 = -0.5
       x2 = 0.5
       y1 = -0.5
       y2 = 0.5
       z1 = -0.5
       z2 = 0.5
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
       bc(BOTTOM)= NO_GRADIENTS
       bc(TOP)   = NO_GRADIENTS
    CASE(SPHERICAL)
       x1 = 0.01
       x2 = RMAX
       y1 = 0.0
       y2 = PI
       z1 = 0.0
       z2 = 2*PI
       bc(WEST)  = REFLECTING       !default: REFLECTING
       bc(EAST)  = REFLECTING       !default: ABSORBING
       bc(SOUTH) = REFLECTING       !default: AXIS
       bc(NORTH) = REFLECTING       !default: AXIS
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
    CASE(OBLATE_SPHEROIDAL)       !TODO anderen Meshs zu 3D konvertieren
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
    mesh => Dict( &
              "meshtype"    / MIDPOINT, &
              "geometry"    / MGEO, &
              "inum"        / XRES, &
              "jnum"        / YRES, &
              "knum"        / ZRES, &
              "xmin"        / x1, &
              "xmax"        / x2, &
              "ymin"        / y1, &
              "ymax"        / y2, &
              "zmin"        / z1, &
              "zmax"        / z2, &
              "gparam"      / GPAR)

    ! boundary conditions
    boundary => Dict( &
              "western"     / bc(WEST), &
              "eastern"     / bc(EAST), &
              "southern"    / bc(SOUTH), &
              "northern"    / bc(NORTH), &
              "bottomer"    / bc(BOTTOM), &
              "topper"      / bc(TOP))

    ! physics settings
    physics => DICT( &
              "problem"     / EULER3D, &
              "gamma"       / GAMMA)

   ! flux calculation and reconstruction method
   fluxes => Dict( &
              "order"       / LINEAR, &
              "fluxtype"    / KT, &
              "variables"   / PRIMITIVE, &
              "limiter"     / VANLEER)

    ! time discretization settings
    timedisc => Dict( &
              "method"      / MODIFIED_EULER, &
              "order"       / 3, &
              "cfl"         / 0.4, &
              "stoptime"    / TSIM, &
              "dtlimit"     / 1.0E-13, &
              "tol_rel"     / 0.01, &
              "tol_abs"     / (/1e-5,1e-5,1e-5,1e-5,1e-5/), &
              "maxiter"     / 1000000)

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"  / VTK, &
              "filename"    / (TRIM(ODIR) // TRIM(OFNAME)), &
              "count"       / ONUM)

    config => Dict( &
              "mesh"        / mesh, &
              "physics"     / physics, &
              "boundary"    / boundary, &
              "fluxes"      / fluxes, &
              "timedisc"    / timedisc, &
              "datafile"    / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER :: n
    REAL    :: P1
    !------------------------------------------------------------------------!
    ! peak pressure
    n  = 3 ! 3 for 3D
    P1 = 3.*(Physics%gamma - 1.0)*E1 / ((n + 1)*PI*R0**n)

    ! uniform density
    Timedisc%pvar(:,:,:,Physics%DENSITY)   = RHO0
    ! vanishing initial velocities
    Timedisc%pvar(:,:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,:,Physics%ZVELOCITY) = 0.
    ! pressure
    WHERE (Mesh%radius%bcenter(:,:,:).LE.R0)
       ! behind the shock front
       Timedisc%pvar(:,:,:,Physics%PRESSURE)  = P1
    ELSEWHERE
       ! in front of the shock front (ambient medium)
       Timedisc%pvar(:,:,:,Physics%PRESSURE)  = P0
    END WHERE


!    DO k=Mesh%KGMIN,Mesh%KGMAX
!      DO j=Mesh%JGMIN,Mesh%JGMAX
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          Timedisc%pvar(i,j,k,Physics%DENSITY)  = RHO0 + RHO1*EXP(-LOG(2.0) &
!               * (radius(i,j,k)/RWIDTH)**2)
!        END DO
!      END DO
!    END DO



    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 3D Sedov explosion")

  END SUBROUTINE InitData

END PROGRAM sedov3d
