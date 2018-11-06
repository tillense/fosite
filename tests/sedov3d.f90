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
!! \cite sedoc1959
!! \cite padmanabhan2000
!----------------------------------------------------------------------------!
PROGRAM sedov3d
  USE fosite_mod
  USE solutions
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 0.05        ! simulation stop time
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
  INTEGER, PARAMETER :: XRES    = 50          ! x-resolution
  INTEGER, PARAMETER :: YRES    = 2           ! y-resolution
  INTEGER, PARAMETER :: ZRES    = 2           ! z-resolution
  REAL, PARAMETER    :: RMAX    = 0.4         ! outer radius of comput. domain
  REAL, PARAMETER    :: GPAR    = 0.2         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME  = 'sedov3d'
  !-------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE  :: Sim
  !-------------------------------------------------------------------------!
  REAL, DIMENSION(:),ALLOCATABLE  :: pvar_diff
#ifdef PARALLEL
  REAL, DIMENSION(:), POINTER     :: pvar,pvar_all,radius, radius_all
#endif
  INTEGER :: i,err
  REAL    :: Rt,rshock
  !-------------------------------------------------------------------------!
  ALLOCATE(Sim,pvar_diff(1:XRES))
  CALL Sim%InitFosite()
#ifdef PARALLEL
  IF(Sim%GetRank().EQ.0) &
#endif
  TAP_PLAN(1)

  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL Sim%Run()
  
#ifdef PARALLEL
  ALLOCATE(pvar(Sim%Mesh%IMIN:Sim%Mesh%IMAX),pvar_all(Sim%GetNumProcs()*(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1)))
  ALLOCATE(radius(Sim%Mesh%IMIN:Sim%Mesh%IMAX),radius_all(Sim%GetNumProcs()*(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1)))
  pvar(:) = Sim%Timedisc%pvar(Sim%Mesh%IMIN:SIM%Mesh%IMAX,1,1,1)
  radius(:)= Sim%Mesh%radius%center(Sim%Mesh%IMIN:Sim%Mesh%IMAX,1,1)
 
  ! Compare results with analytical solution. Check only for correct shock velocity
  ! even if the full analytical solution is implemented

!  CALL sedov(GAMMA,E1,RHO0,P0,TSIM,3, SIM%Mesh%radius%bcenter(1:XRES,1,1),pvar)
  
  !The shock location is where the density gradient is greatest
  CALL MPI_Gather(pvar,int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE, &
    pvar_all, int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE,0,MPI_COMM_WORLD,err)
  CALL MPI_Gather(radius,int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE, &
    radius_all, int(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1),MPI_DOUBLE,0,MPI_COMM_WORLD,err)

 IF(Sim%GetRank().EQ.0) THEN 
  DO i=1,XRES-1
    pvar_diff(i) = pvar_all(i)-pvar_all(i+1)
  END DO
  rshock = radius_all(MAXLOC(pvar_diff,DIM=1))
  
  !analytical shock position
  Rt = R*(E1*TSIM**2/RHO0)**0.2

  !Check whether analytical solution is within +/- one cell of simulated shock
  TAP_CHECK((Rt.LT.rshock+Sim%Mesh%dx).AND.(Rt.GT.rshock-Sim%Mesh%dx),"Shock velocity correct")
END IF
#else
  DO i=1,XRES
      pvar_diff(i) = Sim%Timedisc%pvar%data4d(i,1,1,1)-Sim%Timedisc%pvar%data4d(i+1,1,1,1)
  END DO
  Rshock = Sim%Mesh%radius%center(MAXLOC(pvar_diff,DIM=1),1,1)
  
  !analytical shock position
  Rt = R*(E1*TSIM**2/RHO0)**0.2

  TAP_CHECK((Rt.LT.rshock+Sim%Mesh%dx).AND.(Rt.GT.rshock-Sim%Mesh%dx),"Shock velocity correct")
#endif
  CALL Sim%Finalize()
  TAP_DONE
  DEALLOCATE(Sim,pvar_diff)

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
       z2 = 2*PI
       bc(WEST)  = REFLECTING
       bc(EAST)  = ABSORBING
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
       bc(EAST)  = ABSORBING
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
              "limiter"     / SUPERBEE)

    ! time discretization settings
    timedisc => Dict( &
      "method"      / DORMAND_PRINCE,               &
              "order"       / 5,                            &
              "cfl"         / 0.4,                          &
              "stoptime"    / TSIM,                         &
              "dtlimit"     / 1.0E-13,                      &
              "tol_rel"     / 0.01,                         &
              "tol_abs"     / (/1e-5,1e-5,1e-5,1e-5,1e-5/), &
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
      ! peak pressure
      n  = 3 ! 3 for 3D
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
    Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
    ! pressure
    WHERE (Mesh%radius%bcenter(:,:,:).LE.R0)
       ! behind the shock front
       Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE)  = P1
    ELSEWHERE
       ! in front of the shock front (ambient medium)
       Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE)  = P0
    END WHERE

    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar%data4d,Timedisc%cvar%data4d)
    CALL Mesh%Info(" DATA-----> initial condition: 3D Sedov explosion")

  END SUBROUTINE InitData

END PROGRAM sedov3d
