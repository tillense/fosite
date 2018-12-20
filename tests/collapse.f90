!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: collapse.f90                                                      #
!#                                                                           #
!# Copyright (C) 2008-2012                                                   #
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
!> collapse with/without angular momentum and self-gravity
!! \author Tobias Illenseer
!! \author Björn Sperling
!!
!! References:
!! [1] Norman, M. L.; Wilson, J. R.; Barton, R. T.
!!     "A new calculation on rotating protostar collapse"
!!     Astrophysical Journal, Part 1, vol. 239, Aug. 1, 1980, p. 968-981.
!!     DOI: 10.1086/158185
!! [2] Colgate, Stirling A.; White, Richard H.
!!     "The Hydrodynamic Behavior of Supernovae Explosions"
!!     Astrophysical Journal, vol. 143, p.626
!!     DOI: 10.1086/148549
!----------------------------------------------------------------------------!
PROGRAM collapse
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE sources_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM      = 1.0E-0 ! simulation time in terms of the
                                           !   free-fall time [TAU]
  REAL, PARAMETER    :: GAMMA     = 1.4    ! ratio of specific heats
  REAL, PARAMETER    :: MASS      = 1.0E+2 ! mass of spheroid
  REAL, PARAMETER    :: CENTMASS  = 0.0E+2 ! central pointmass (0.0 to disable)
  REAL, PARAMETER    :: RSPH      = 30.0   ! semi-minor axis of the spheroid
  REAL, PARAMETER    :: ECC       = 0.0    ! eccentricity (0.0 is sphere)
  REAL, PARAMETER    :: VOL0      = 4*PI/3 * RSPH*RSPH*RSPH / (1.-ECC*ECC)
                                           ! volume of the spheroid
  REAL, PARAMETER    :: OMEGA     = 0.0E-7 ! angular velocity (0.0 to disable)
  REAL, PARAMETER    :: ETA_P     = 100.0  ! ratio of p_(hydro_static) to p 
                                           ! (in case of self-gravity) approx
                                           ! 100 (free fall); approx 1 stable
                                           ! (without pointmass)
  REAL, PARAMETER    :: ETA_RHO   = 1.0E-6 ! density ratio rho / rho_inf 
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry of the mesh
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL
!!$  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES = 31          ! x-resolution
  INTEGER, PARAMETER :: YRES = 31          ! y-resolution
  REAL, PARAMETER    :: RMAX = 1.5         ! size of comput. domain in [RSPH]
  REAL, PARAMETER    :: GPAR = 0.5*RSPH    ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'collapse'
  !--------------------------------------------------------------------------!
  REAL               :: TAU                ! free-fall time scale
  REAL               :: RHO0               ! initial density of the sphere
  REAL               :: P0                 ! initial hydrostatic pressure
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

! multigrid not implemented for MPI runs
#ifndef PARALLEL
  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
#else
  TAP_CHECK(.TRUE.,"Nothing to do for MPI runs.")
#endif
  TAP_DONE

  CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4),sgbc
    DOUBLE PRECISION, PARAMETER :: GN = 6.6742E-11     ! [m^3/kg/s^2]
    TYPE(Dict_TYP), POINTER :: mesh,physics, boundary, datafile, logfile, &
                               grav,pmass, timedisc, fluxes, selfgravity
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x2 = RMAX*RSPH
       x1 = 2.*x2 / (XRES+2)        ! x_min = 2*dx
       y1 = 0.0
       y2 = PI
    CASE(CYLINDRICAL)
       x1 = -RMAX*RSPH
       x2 = RMAX*RSPH
       y1 = 0.0
       y2 = RMAX*RSPH
    CASE(OBLATE_SPHEROIDAL)
       x2 = RMAX*RSPH/GPAR
       x2 = LOG(x2+SQRT(x2**2-1.0)) ! = ACOSH(RMAX*RSPH/GPAR)
       x1 = 2.*x2/(XRES+2)          ! x_min = 2*dx
       y1 = -0.5*PI
       y2 = 0.5*PI
    CASE(TANCYLINDRICAL)
       x1 = ATAN(-RMAX*RSPH/GPAR)
       x2 = ATAN(RMAX*RSPH/GPAR)
       y1 = 0.0
       y2 = RMAX*RSPH
    CASE(SINHSPHERICAL)
       x2 = RMAX*RSPH/GPAR
       x2 = LOG(x2+SQRT(x2**2+1.0)) ! = ASINH(RMAX*RSPH/GPAR)
       x1 = 2.*x2/(XRES+2)          ! x_min = 2*dx
       y1 = 0.0
       y2 = PI
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
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
       sgbc = SPHERMULTEXPAN        ! use spherical multipole expansion for BC
                                    !   in the multigrid poisson solver
    CASE(CYLINDRICAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = ABSORBING
       sgbc = CYLINMULTEXPAN        ! cylindrical multipole expansion
    CASE(OBLATE_SPHEROIDAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
       sgbc = CYLINMULTEXPAN        ! cylindrical multipole expansion
    CASE(TANCYLINDRICAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = NO_GRADIENTS
       sgbc = CYLINMULTEXPAN        ! cylindrical multipole expansion
    CASE(SINHSPHERICAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
       sgbc = SPHERMULTEXPAN        ! spherical multipole expansion
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","geometry not supported for 3D Sedov explosion")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))
    
    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA)                 ! ratio of specific heats        !

     ! compute some derived simulation parameters
    RHO0 = MASS / VOL0                 ! initial density within the spheroid !
    ! "hydrostatic" pressure * ETA_P 
    !     => with ETA_P approx 100 => free-fall (in case of self-gravity)
    P0 = 4.0/3.0*PI*GN*RHO0**2*RSPH**2 / ETA_P
    ! free-fall time (at radius RSPH) with contributions from both
    ! the selfgravitating spheroid and the central point mass
    TAU = SQRT((RSPH**3)/GN/(4./3.*PI*RSPH**3*RHO0 + CENTMASS))
 
    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! source term due to a point mass
    pmass => Dict("gtype" / POINTMASS, &        ! grav. accel. of a point mass   !
            "mass"  / CENTMASS)          ! mass [kg]                      !

    ! source term due to self-gravity
    selfgravity => Dict( &
          "gtype"        / MULTIGRID, &   ! poisson solver for self-gravity!
!          "maxmult"      / 5, &           ! number of (spher.) multipol moments
!          "maxresidnorm" / 1.0E-7, &     ! accuracy of multigrid solver (max error)
          "relaxtype"    / BLOCK_GAUSS_SEIDEL, &      ! relaxation method
!          "relaxtype"    / RED_BLACK_GAUSS_SEIDEL, &
!          "relaxtype"    / GAUSS_SEIDEL, &
!          "npre"         / 1, &           ! number of pre smoothings
!          "npost"        / 1, &           ! and post smoothings
!          "minres"       / 3, &           ! resolution of coarsest grid
!          "nmaxcycle"    / 250, &         ! limit for iterations
          "bndrytype"    / sgbc)         ! multipole expansion (see above)


    IF (CENTMASS.GT.TINY(CENTMASS)) THEN
      grav => Dict("stype" / GRAVITY, & 
                   "self" / selfgravity, &
                   "point" / pmass)
    ELSE 
      grav => Dict("stype" / GRAVITY, & 
                   "self" / selfgravity)
    END IF

    ! time discretization settings
    timedisc => Dict( &
         "method"    / MODIFIED_EULER, &
         "order"     / 3, &
         "cfl"       / 0.4, &
         "stoptime"  / (TSIM * TAU), &
         "dtlimit"   / 1.0E-9, &
         "maxiter"   / 10000000)

    ! initialize data input/output
    datafile => Dict( &
!!$           "fileformat" / VTK, &
           "fileformat" / GNUPLOT, "filecycles" / 0, "cartcoords" / 1, &
           "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
           "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources/grav"  / grav, &
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
    CHARACTER(LEN=64) :: value
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = 0.
    ! angular velocity (angular velocity) * (axis distance)
    ! axis distance in 2.5D rotationally symmetric coordinates
    ! is always given by the scale factor hz
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = OMEGA * Mesh%hz%bcenter(:,:)
    
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       ! specific angular momentum
       Timedisc%pvar(:,:,Physics%ZVELOCITY) = &
            Timedisc%pvar(:,:,Physics%ZVELOCITY)*Mesh%hz%bcenter(:,:)
    END IF

    ! density
    WHERE ((1.0-ECC*ECC)*Mesh%cart%bcenter(:,:,1)**2+Mesh%cart%bcenter(:,:,2)**2.LE.RSPH*RSPH)
       Timedisc%pvar(:,:,Physics%DENSITY) = RHO0
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY) = RHO0 * ETA_RHO
    END WHERE

    ! pressure
    Timedisc%pvar(:,:,Physics%PRESSURE) = P0
   
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: " // &
         "Homogenious density with uniform angular frequency")
    IF (CENTMASS.GT.TINY(CENTMASS)) THEN
       WRITE (value,"(E12.4)") SQRT(RSPH**3/(CENTMASS*Physics%constants%GN))
       CALL Info(Mesh,"                               " // &
          "timescale of pointmass:    " //trim(value))
    END IF
    IF (MASS .GT. TINY(MASS)) THEN
       WRITE(value,"(E12.4)") SQRT(3.0*PI/(4.0*RHO0*Physics%constants%GN))
       CALL Info(Mesh,"                               " // &
          "timescale of self-gravity: " //trim(value))
    END IF
    IF (OMEGA .GT. TINY(OMEGA)) THEN
       WRITE (value,"(E12.4)") OMEGA
       CALL Info(Mesh,"                               " // &
          "angular velocity:          " // trim(value))
    END IF

    WRITE (value,"(E12.4)") TAU
    CALL Info(Mesh,"                               " // &
          "free-fall time:            " //trim(value))

  END SUBROUTINE InitData

END PROGRAM collapse
