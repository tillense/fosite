!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: collapse.f90                                                      #
!#                                                                           #
!# Copyright (C) 2008-2021                                                   #
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
!> \test collapse with/without angular momentum and self-gravity
!! \author Tobias Illenseer
!! \author Björn Sperling
!! \author Lars Bösch
!!
!! References:
!! - \cite norman1980 Norman, M. L.; Wilson, J. R.; Barton, R. T.
!!     "A new calculation on rotating protostar collapse"
!!     Astrophysical Journal, Part 1, vol. 239, Aug. 1, 1980, p. 968-981.
!!     DOI: 10.1086/158185
!! - \cite colgate1966  Colgate, Stirling A.; White, Richard H.
!!     "The Hydrodynamic Behavior of Supernovae Explosions"
!!     Astrophysical Journal, vol. 143, p.626
!!     DOI: 10.1086/148549
!----------------------------------------------------------------------------!
PROGRAM collapse
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  REAL, PARAMETER    :: GN        = 6.6742E-11     ! [m^3/kg/s^2]
  ! simulation parameters
  REAL, PARAMETER    :: TSIM      = 1.0E-0!/1000. ! simulation time in terms of the
                                           !   free-fall time [TAU]
  REAL, PARAMETER    :: GAMMA     = 1.4    ! ratio of specific heats
  REAL, PARAMETER    :: MASS      = 1.0E+2 ! mass of spheroid
  REAL, PARAMETER    :: CENTMASS  = 1.0E+2 ! central pointmass (0.0 to disable)
  REAL, PARAMETER    :: RSPH      = 30.0   ! semi-minor axis of the spheroid
  REAL, PARAMETER    :: ECC       = 0.0    ! eccentricity (0.0 is sphere)
  REAL, PARAMETER    :: VOL0      = 4*PI/3 * RSPH*RSPH*RSPH / (1.-ECC*ECC)
                                           ! volume of the spheroid
  REAL, PARAMETER    :: OMEGA     = 1.0E-7 ! angular velocity (0.0 to disable)
  REAL, PARAMETER    :: OMEGA_FRAME = 0.0*Omega ! angular velocity (0.0 to disable)
  REAL, PARAMETER    :: ETA_P     = 100.0  ! ratio of p_(hydro_static) to p 
                                           ! (in case of self-gravity) approx
                                           ! 100 (free fall); approx 1 stable
                                           ! (without pointmass)
  REAL, PARAMETER    :: ETA_RHO   = 1.0E-6 ! density ratio rho / rho_inf 
  ! mesh settings
!  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry of the mesh
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
  INTEGER, PARAMETER :: XRES = 150          ! x-resolution
  INTEGER, PARAMETER :: YRES = 1 !XRES          ! y-resolution
  INTEGER, PARAMETER :: ZRES = XRES           ! z-resolution
  REAL, PARAMETER    :: RMAX = 1.5         ! size of comput. domain in [RSPH]
  REAL, PARAMETER    :: GPAR = 0.5*RSPH    ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'collapse_cyl_rhs1_150'
  !--------------------------------------------------------------------------!
  REAL               :: TAU                ! free-fall time scale
  REAL               :: RHO0               ! initial density of the sphere
  REAL               :: P0                 ! initial hydrostatic pressure
  !--------------------------------------------------------------------------!
  CLASS(fosite),ALLOCATABLE   :: Sim
  LOGICAL  :: ok
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

! multigrid not implemented for MPI runs
#ifndef PARALLEL
  ALLOCATE(SIM)
  CALL Sim%InitFosite()

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL Sim%Setup()

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL Sim%Run()
  ok = .NOT.Sim%aborted

  CALL Sim%Finalize()

  TAP_CHECK(ok,"Finished simulation")
#else
  TAP_CHECK(.TRUE.,"Nothing to do for MPI runs.")
#endif
  TAP_DONE

  CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6),sgbc
    TYPE(Dict_TYP), POINTER :: mesh,physics, boundary, datafile, rotframe, &
                               grav,pmass, timedisc, fluxes, selfgravity
    REAL              :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x2 = RMAX*RSPH
       x1 = 2.*x2 / (XRES+2)        ! x_min = 2*dx
       y1 = 0.005
       y2 = PI-0.005
       z1 = 0.0
       z2 = 2.*PI
    CASE(CYLINDRICAL)
       x2 = RMAX*RSPH
       x1 = 2.*x2 / (XRES+2)        ! x_min = 2*dx
       y1 = 0.0
       y2 = 2.*PI
       z1 = -RMAX*RSPH
       z2 = RMAX*RSPH
!    CASE(OBLATE_SPHEROIDAL)
!       x2 = RMAX*RSPH/GPAR
!       x2 = LOG(x2+SQRT(x2**2-1.0)) ! = ACOSH(RMAX*RSPH/GPAR)
!       x1 = 2.*x2/(XRES+2)          ! x_min = 2*dx
!       y1 = -0.5*PI
!       y2 = 0.5*PI
    CASE(TANCYLINDRICAL)
       x1 = ATAN(-RMAX*RSPH/GPAR)
       x2 = ATAN(RMAX*RSPH/GPAR)
       y2 = RMAX*RSPH
       y1 = 2.*y2 / (YRES+2)        ! x_min = 2*dx
       z1 = 0.0
       z2 = 2.*PI
!    CASE(SINHSPHERICAL)
!       x2 = RMAX*RSPH/GPAR
!       x2 = LOG(x2+SQRT(x2**2+1.0)) ! = ASINH(RMAX*RSPH/GPAR)
!       x1 = 2.*x2/(XRES+2)          ! x_min = 2*dx
!       y1 = 0.0
!       y2 = PI
    CASE DEFAULT
       CALL Sim%Error("InitProgram","geometry not supported for collapse simulation")
    END SELECT

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "knum"     / ZRES, &
           "omega"    / Omega_Frame, &
           "use_fargo"/ 0,        &
           "fargo"    / 0,        &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "zmin"     / z1, &
           "zmax"     / z2, &
           "output/radius" / 1, &
           "output/volume" / 1, &
           "gparam"   / GPAR)
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       bc(WEST)  = REFLECTING !ABSORBING
       bc(EAST)  = REFLECTING !ABSORBING
       bc(SOUTH) = REFLECTING !AXIS
       bc(NORTH) = REFLECTING !AXIS
       bc(BOTTOM)= PERIODIC
       bc(TOP)   = PERIODIC
!       sgbc = SPHERMULTEXPAN        ! use spherical multipole expansion for BC
                                    !   in the multigrid poisson solver
    CASE(CYLINDRICAL)
       bc(WEST)  = REFLECTING !ABSORBING
       bc(EAST)  = REFLECTING !ABSORBING
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
       bc(BOTTOM)= REFLECTING !AXIS
       bc(TOP)   = REFLECTING !ABSORBING
!       sgbc = CYLINMULTEXPAN        ! cylindrical multipole expansion
!    CASE(OBLATE_SPHEROIDAL)
!       bc(WEST)  = ABSORBING
!       bc(EAST)  = ABSORBING
!       bc(SOUTH) = AXIS
!       bc(NORTH) = AXIS
!       sgbc = CYLINMULTEXPAN        ! cylindrical multipole expansion
    CASE(TANCYLINDRICAL)
       bc(WEST)  = REFLECTING !ABSORBING
       bc(EAST)  = REFLECTING !ABSORBING
       bc(SOUTH) = REFLECTING !AXIS
       bc(NORTH) = REFLECTING !NO_GRADIENTS
       bc(BOTTOM)= PERIODIC
       bc(TOP)   = PERIODIC
!       sgbc = CYLINMULTEXPAN        ! cylindrical multipole expansion
!    CASE(SINHSPHERICAL)
!       bc(WEST)  = ABSORBING
!       bc(EAST)  = ABSORBING
!       bc(SOUTH) = AXIS
!       bc(NORTH) = AXIS
!       sgbc = SPHERMULTEXPAN        ! spherical multipole expansion
    CASE DEFAULT
       CALL Sim%Error("InitProgram","geometry not supported for collapse simulation")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH), &
               "bottomer" / bc(BOTTOM), &
               "topper" / bc(TOP))
    
    ! physics settings
    physics => Dict("problem" / EULER, &
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
             "theta"     / 1.8)          ! optional parameter for limiter !

    ! source term due to a point mass
    pmass => Dict("gtype" / POINTMASS, &        ! grav. accel. of a point mass   !
            "mass"  / CENTMASS)          ! mass [kg]                      !

    ! source term due to self-gravity
    selfgravity => Dict( &
          "gtype"        / SPECTRAL)!, &   ! poisson solver for self-gravity!
!          "gtype"        / MULTIGRID, &   ! poisson solver for self-gravity!
!          "maxmult"      / 5, &           ! number of (spher.) multipol moments
!          "maxresidnorm" / 1.0E-7, &     ! accuracy of multigrid solver (max error)
!          "relaxtype"    / BLOCK_GAUSS_SEIDEL, &      ! relaxation method
!          "relaxtype"    / RED_BLACK_GAUSS_SEIDEL, &
!          "relaxtype"    / GAUSS_SEIDEL, &
!          "npre"         / 1, &           ! number of pre smoothings
!          "npost"        / 1, &           ! and post smoothings
!          "minres"       / 3, &           ! resolution of coarsest grid
!          "nmaxcycle"    / 250, &         ! limit for iterations
!          "bndrytype"    / sgbc)         ! multipole expansion (see above)
    IF (OMEGA_FRAME.GT.TINY(OMEGA_FRAME)) THEN
       rotframe => Dict("stype" / ROTATING_FRAME, &
               "x"     / 0.0, &
               "y"     / 0.0)
!       sources => Dict("rotframe" / rotframe)
    END IF


    IF (CENTMASS.GT.TINY(CENTMASS)) THEN
      grav => Dict("stype" / GRAVITY, & 
!                   "self" / selfgravity, &
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
         "rhstype"   / 1, &
         "maxiter"   / 10000000)

    ! initialize data input/output
    datafile => Dict( &
!           "fileformat" / VTK, &
           "fileformat" / BINARY, &
           "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
           "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources/grav"  / grav, &
!             "sources/rotframe" / rotframe, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Mesh_base)    :: Mesh
    CLASS(Physics_base) :: Physics
    CLASS(Timedisc_base):: Timedisc
    !------------------------------------------------------------------------!
    CHARACTER(LEN=64) :: value
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    IF(MGEO.EQ.SPHERICAL) THEN
      ! velocities
      Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
      Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
      ! angular velocity (angular velocity) * (axis distance)
      ! axis distance in 2.5D rotationally symmetric coordinates
      ! is always given by the scale factor hz
      Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = (OMEGA-OMEGA_FRAME) * Mesh%hz%bcenter(:,:,:)
    ELSE IF(MGEO.EQ.CYLINDRICAL) THEN
      Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
      Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
      Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = (OMEGA-OMEGA_FRAME) * Mesh%hy%bcenter(:,:,:)
    ELSE IF(MGEO.EQ.TANCYLINDRICAL) THEN
      Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
      Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
      Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = (OMEGA-OMEGA_FRAME) * Mesh%hz%bcenter(:,:,:)
    END IF
    
!    IF (Physics%GetType().EQ.EULER3D_ROTAMT) THEN
!       ! specific angular momentum
!       Timedisc%pvar(:,:,Physics%ZVELOCITY) = &
!            Timedisc%pvar(:,:,Physics%ZVELOCITY)*Mesh%hz%bcenter(:,:)
!    END IF

    ! density
!    DO k=Mesh%KGMIN, Mesh%KGMAX
!      DO j=Mesh%JGMIN, Mesh%JGMAX
!        DO i=Mesh%IGMIN, MEsh%IGMAX
!          IF(Mesh%radius%bcenter(i,j,k).LE.RSPH) THEN
!            Timedisc%pvar%data4d(i,j,k,Physics%DENSITY) = RHO0
!          ELSE
!            Timedisc%pvar%data4d(i,j,k,Physics%DENSITY) = RHO0 * ETA_RHO
!          END IF
!        END DO
!      END DO
!    END DO
    WHERE (Mesh%radius%bcenter(:,:,:).LE.RSPH)
       Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
    ELSEWHERE
       Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0 * ETA_RHO
    END WHERE

    ! pressure
    Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
   
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Homogenious density with uniform angular frequency")
    IF (CENTMASS.GT.TINY(CENTMASS)) THEN
       WRITE (value,"(E12.4)") SQRT(RSPH**3/(CENTMASS*Physics%constants%GN))
       CALL Mesh%Info("                               " // &
          "timescale of pointmass:    " //trim(value))
    END IF
    IF (MASS .GT. TINY(MASS)) THEN
       WRITE(value,"(E12.4)") SQRT(3.0*PI/(4.0*RHO0*Physics%constants%GN))
       CALL Mesh%Info("                               " // &
          "timescale of self-gravity: " //trim(value))
    END IF
    IF (OMEGA .GT. TINY(OMEGA)) THEN
       WRITE (value,"(E12.4)") OMEGA
       CALL Mesh%Info("                               " // &
          "angular velocity:          " // trim(value))
    END IF

    WRITE (value,"(E12.4)") TAU
    CALL Mesh%Info("                               " // &
          "free-fall time:            " //trim(value))
    SELECT TYPE(bwest => Timedisc%Boundary%Boundary(WEST)%p)
    CLASS IS(boundary_custom)
      IF(MGEO.EQ.SPHERICAL) THEN
        CALL bwest%SetCustomBoundaries(Mesh,Physics, &
          (/CUSTOM_NOGRAD,CUSTOM_REFLECT,CUSTOM_REFLECT,CUSTOM_FIXED,CUSTOM_REFLECT/))
        bwest%data(:,:,:,Physics%ZVELOCITY) = 0.0
      ELSE IF(MGEO.EQ.CYLINDRICAL) THEN
        CALL bwest%SetCustomBoundaries(Mesh,Physics, &
!          (/CUSTOM_NOGRAD,CUSTOM_REFLECT,CUSTOM_FIXED,CUSTOM_REFLECT/))
          (/CUSTOM_NOGRAD,CUSTOM_REFLECT,CUSTOM_FIXED,CUSTOM_NOGRAD,CUSTOM_REFLECT/))
        bwest%data(:,:,:,Physics%YVELOCITY) = 0.0
      END IF
    END SELECT

    SELECT TYPE(bbottom => Timedisc%Boundary%Boundary(BOTTOM)%p)
    CLASS IS(boundary_fixed)
      IF(MGEO.EQ.CYLINDRICAL) THEN
        bbottom%fixed(:,:,:) = .TRUE.
        bbottom%data(:,:,1,:) = Timedisc%pvar%data4d(:,:,Mesh%KMIN,:)
        bbottom%data(:,:,2,:) = Timedisc%pvar%data4d(:,:,Mesh%KMIN,:)
      END IF
    CLASS IS(boundary_farfield)
        bbottom%data(:,:,1,:) = Timedisc%pvar%data4d(:,:,Mesh%KMIN,:)
        bbottom%data(:,:,2,:) = Timedisc%pvar%data4d(:,:,Mesh%KMIN,:)
    END SELECT

    SELECT TYPE(btop => Timedisc%Boundary%Boundary(TOP)%p)
    CLASS IS(boundary_fixed)
      IF(MGEO.EQ.CYLINDRICAL) THEN
        btop%fixed(:,:,:) = .TRUE.
        btop%data(:,:,1,:) = Timedisc%pvar%data4d(:,:,Mesh%KMAX,:)
        btop%data(:,:,2,:) = Timedisc%pvar%data4d(:,:,Mesh%KMAX,:)
      END IF
    CLASS IS(boundary_farfield)
        btop%data(:,:,1,:) = Timedisc%pvar%data4d(:,:,Mesh%KMAX,:)
        btop%data(:,:,2,:) = Timedisc%pvar%data4d(:,:,Mesh%KMAX,:)
    END SELECT


  END SUBROUTINE InitData

END PROGRAM collapse

