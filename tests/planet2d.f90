!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: planet2d.f90                                                      #
!#                                                                           #
!# Copyright (C) 2006-2019                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee         <jklee@astrophysik.uni-kiel.de>                       #
!# Jubin Lirawi      <jlirawi@astrophysik.uni-kiel.de>                       #
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
!> \example planet2d.f90
!!
!! \author Tobias Illenseer
!! \author Jannes Klee
!! \author Jubin Lirawi
!!
!! \brief initialisation of a 2D planetary atmosphere simulation
!!
!! \warning use SI units
!----------------------------------------------------------------------------!
PROGRAM planet2d
  USE fosite_mod
  USE sources_base_mod, ONLY : sources_base
  USE sources_planetcooling_mod, ONLY : sources_planetcooling
  USE constants_SI_mod, ONLY : C, GN, KB, NA, SB, KE, AU, &
       MSUN, MJUPITER, MEARTH, RSUN, RJUPITER, REARTH, DAY
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants                               !                        !
!   REAL, PARAMETER    :: GN      = 6.6742D-11        ! Newtons grav. constant !
!   REAL, PARAMETER    :: CC      = 2.99792458D+8     ! speed of light   [m/s] !
!   REAL, PARAMETER    :: RG      = 8.31447           ! molar gas constant     !
!   REAL, PARAMETER    :: SBconst = 5.670367D-8       ! Stefan-Boltzmann const !
!   REAL, PARAMETER    :: AU      = 1.49597870691E+11 ! astronomical unit  [m] !
!   REAL, PARAMETER    :: YEAR    = 3.15576E+7        ! Julian year        [s] !
!   REAL, PARAMETER    :: DAY     = 8.6400E+4         ! Solar Day          [s] !
!   REAL, PARAMETER    :: REARTH  = 6.371E+6          ! radius of earth    [m] !
!   REAL, PARAMETER    :: RJUP    = 71.492E+6         ! radius of jupiter  [m] !
!   REAL, PARAMETER    :: RSUN    = 6.957E+8          ! radius of the sun  [m] !
!   REAL, PARAMETER    :: MEARTH  = 5.9723D+24        ! mass of the earth [kg] !
  ! planetary parameters                            !                        !
  REAL, PARAMETER    :: TYEAR   = 4.05 * DAY        ! tropical year      [s] !
  REAL, PARAMETER    :: RPLANET = 0.788*REARTH      ! planetary radius   [m] !
  REAL, PARAMETER    :: THETA0  = 0.0               ! axis-plane-angle [rad] !
  REAL, PARAMETER    :: PHI0    = 0.0               !                  [rad] !
  REAL, PARAMETER    :: OMEGA   = 2*PI/TYEAR/8      ! ang. rotation  [rad/s] !
  REAL, PARAMETER    :: FSUN    = OMEGA/(2*PI) - &  ! freq.(!) day-night[/s] !
                                     1.0/TYEAR      !                        !
  REAL, PARAMETER    :: MASS    = 0.297*MEARTH      ! mass of the planet [kg]!
  REAL, PARAMETER    :: GACC    = GN*MASS/RPLANET**2! grav. accel.  [m/s**2] !
  ! orbital parameters                              !                        !
  REAL, PARAMETER    :: SM_AXIS = 0.02219*AU        ! semi major axis    [m] !
  REAL, PARAMETER    :: EXCENT  = 0.0               ! numerical eccentricity !
  ! stellar parameters                              !                        !
  REAL, PARAMETER    :: RSTAR   = 0.121*RSUN        ! radius of the star [m] !
  REAL, PARAMETER    :: TSTAR   = 2511              ! eff. temperature   [K] !
  REAL, PARAMETER    :: FSTAR   = (RSTAR/AU)**2 &   ! rad flux density @ 1 AU!
                                  * SB*TSTAR**4     !   [W/m**2]             !
  ! simulation parameters                           !                        !
  REAL, PARAMETER    :: TSIM    = 10 * TYEAR        ! simulation stop time   !
  INTEGER, PARAMETER :: XRES    = 1                 ! radius-resolution [px] !
  INTEGER, PARAMETER :: YRES    = 16                ! theta-resolution  [px] !
  INTEGER, PARAMETER :: ZRES    = 32                ! phi-resolution    [px] !
  ! gas properties of the atmosphere
  REAL, PARAMETER    :: GAMMA   = 1.4               ! ratio of specific heats!
  REAL, PARAMETER    :: P0      = 1.0E7             ! surf. press.      [Pa] !
  REAL, PARAMETER    :: MU      = 2.8586E-2         ! molar mass    [kg/mol] !
  REAL, PARAMETER    :: T0      = 258.1             ! mean equil. temp.  [K] !
  ! optical properties of the atmosphere
  REAL, PARAMETER    :: ALBEDO  = 0.306             ! albedo of the planet   !
  ! output parameters                               !                        !
  INTEGER, PARAMETER :: ONUM    = 10                ! num. output data sets  !
  CHARACTER(LEN=256), PARAMETER &                   ! output data dir        !
                     :: ODIR    = './'              !                        !
  CHARACTER(LEN=256), PARAMETER &                   ! output data file name  !
                     :: OFNAME  = 'planet2d'        !                        !
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  CLASS(sources_base), POINTER :: sp => null()
  CLASS(sources_planetcooling), POINTER :: spcool => null()
  REAL :: Tmean
  LOGICAL            :: ok
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  ! get pointer on planetary cooling source term
  IF (ALLOCATED(Sim%Sources)) THEN
    sp => Sim%Sources%GetSourcesPointer(PLANET_COOLING)
    IF (ASSOCIATED(sp)) THEN
      SELECT TYPE(sp)
      CLASS IS(sources_planetcooling)
        spcool => sp
      CLASS DEFAULT
        NULLIFY(spcool)
      END SELECT
    END IF
  END IF
  IF (.NOT.ASSOCIATED(spcool)) &
    CALL Sim%Error("planet2d","planetary cooling sources not initialized")
  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
  ! run simulation
  CALL Sim%Run()
  ! compute final global mean surface temperature using local surface temperatures
  ! computed by the planetary cooling module
  Tmean = SUM(spcool%T_s%data1d(:)*Sim%Mesh%volume%data1d(:),MASK=Sim%Mesh%without_ghost_zones%mask1d(:)) &
    / SUM(Sim%Mesh%volume%data1d(:),MASK=Sim%Mesh%without_ghost_zones%mask1d(:))
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

  TAP_CHECK_CLOSE(T0,Tmean,T0*0.01,"deviation of final mean surface temperature from initial T0 < 1%")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP), POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile,  &
                               timedisc, fluxes, heating, sources, &
                               grav, pmass, cooling, rotframe
    !------------------------------------------------------------------------!
    !mesh settings
    mesh => Dict( &
         "meshtype"        / MIDPOINT, &
         "geometry"        / SPHERICAL_PLANET, &
         "decomposition"   / (/ 1, -1, 1/), &
         "omega"           / OMEGA, &
         "inum"            / YRES, &
         "jnum"            / ZRES, &
         "knum"            / XRES, &
         "xmin"            / (0.05), &
         "xmax"            / (PI-0.05), &
         "ymin"            / (0.0), &
         "ymax"            / (2.0*PI), &
         "zmin"            / RPLANET, &
         "zmax"            / RPLANET, &
         "gparam"          / RPLANET, &
         "output/rotation" / 0, &
         "output/volume"   / 1, &
         "output/dAz"      / 1)

    ! boundary conditions
    boundary => Dict( &
         "western"         / REFLECTING, &
         "eastern"         / REFLECTING, &
         "southern"        / PERIODIC, &
         "northern"        / PERIODIC, &
         "bottomer"        / REFLECTING, &
         "topper"          / REFLECTING)

    ! physics settings
    physics => Dict( &
         "problem"         / EULER, &
         "mu"              / MU, &
         "gamma"           / GAMMA)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
         "fluxtype"        / KT, &
         "order"           / LINEAR, &
!          "variables"       / CONSERVATIVE, &
         "variables"       / PRIMITIVE, &
         "limiter"         / VANLEER)

    ! rotating frame with disabled centrifugal forces
    ! Rotating planets are usually in an equilibrium state in which
    ! gravitational and centrifugal forces balance, so that the planet
    ! becomes a spheroid instead of a sphere. On the spheroid the
    ! tangential components of the centrifugal forces vanish and
    ! one can account for the vertical component defining an
    ! effective gravitational acceleration (in 3D simulations).
    rotframe => Dict( &
         "stype"           / ROTATING_FRAME, &
         "disable_centaccel" / 1)

    ! cooling in infrared
    cooling => Dict( &
         "stype"           / PLANET_COOLING, &
         "output/Qcool"    / 1, &
         "output/T_s"      / 1, &
         "output/RHO_s"    / 1, &
         "output/P_s"      / 1, &
         "radflux"         / FSTAR, &
         "albedo"          / ALBEDO, &
         "mean_distance"   / (SM_AXIS*(1. + 0.5*EXCENT**2)), &
         "T_0"             / T0, &
         "cvis"            / 0.1, &
         "gacc"            / GACC)

    ! heating by a star
    heating => Dict( &
         "stype"           / PLANET_HEATING, &
         "output/Qstar"    / 1, &
         "year"            / TYEAR, &
         "theta0"          / THETA0, &
         "phi0"            / PHI0, &
         "omegasun"        / FSUN, &
         "albedo"          / ALBEDO, &
         "radflux"         / FSTAR,&
         "semimajoraxis"   / SM_AXIS, &
         "excentricity"    / EXCENT, &
         "cvis"            / 0.1)

    sources => Dict( &
         "rotframe"        / rotframe, &
         "cooling"         / cooling, &
         "heating"         / heating)

    ! time discretization settings
    timedisc => Dict( &
         "method"          / SSPRK,   &
         "cfl"             / 0.4,     &
         "stoptime"        / TSIM,    &
         "tol_rel"         / 0.0095,  &
         "dtlimit"         / 1.0E-15, &
         "maxiter"         / 1000000000)

    datafile => Dict(&
         "fileformat"      / VTK, &
         "filename"        / (TRIM(ODIR) // TRIM(OFNAME)), &
         "count"           / ONUM)

    config => Dict( &
         "mesh"            / mesh,     &
         "physics"         / physics,  &
         "boundary"        / boundary, &
         "fluxes"          / fluxes,   &
         "timedisc"        / timedisc, &
         "sources"         / sources,  &
         "datafile"        / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),             INTENT(IN)    :: Physics
    CLASS(mesh_base),                INTENT(IN)    :: Mesh
    CLASS(marray_compound), POINTER, INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_euler)
      ! assuming hydrostatic equillibrium with constant (effective) gravitational
      ! acceleration yields constant column density, i.e. vertically integrated density,
      ! given by the ratio of mean surface pressure P0 and (effective) grav. acceleration GACC
      p%density%data1d(:) = P0/GACC

      ! assuming hydrostic equillibrium with mean surface
      ! temperature T0 yields the initial condition for the 2D vertically
      ! integrated pressure
      p%pressure%data1d(:) = GAMMA*Physics%Constants%RG / MU * T0 * p%density%data1d(:)

      ! vanishing velocities in comoving frame
      p%velocity%data1d(:) = 0.0
    CLASS DEFAULT
      CALL Mesh%Error("planet2d::InitData","state vector not supported")
    END SELECT

    CALL Physics%Convert2Conservative(pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 2D planetary atmosphere")

  END SUBROUTINE InitData

END PROGRAM planet2d
