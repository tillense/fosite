!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mestel.f90                                                        #
!#                                                                           #
!# Copyright (C) 2012-2019                                                   #
!# Manuel Jung    <mjung@astrophysik.uni-kiel.de>                            #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> self-gravitating accretion disk
!!
!! \author Manuel Jung
!! \author Björn Sperling
!! \author Tobias Illenseer
!!
!! 2D simulation of a geometrically thin, self-gravitating accretion disk
!! around a supermassive black hole in polar geometry with logarithmic
!! radial spacing. The standard setup solves the non-isothermal inviscid
!! Euler equations with thin-disk gray cooling (see \ref sources_diskcooling_mod.lambda_gray ).
!! Gravitational forces account for the central point mass as well as self-gravity
!! of the disk.
!!
!! The setup is based on those described in \cite britsch2006 .
!!
!! <div class="row"> <div class="col-md-6">
!!  Simulation parameters         | \f$ \quad \f$
!!  ------------------            | -----------------
!!  black hole mass               | \f$ 10^7\, \mathsf{M}_\odot \f$
!!  disk / black hole mass ratio  | \f$ 0.1 \f$
!!  mean molecular weight         | \f$ 6.02\cdot 10^{-4}\, \mathsf{kg/mol} \f$
!!  specific heat ratio           | \f$ 1.4 \f$
!!  inner radius                  | \f$ 0.05\, \mathsf{pc} \f$
!!  outer radius                  | \f$ 1\, \mathsf{pc} \f$
!!
!!  Initial condition                              | \f$ \quad \f$
!!  ------------------                             | -----------------
!!  power law surface density, i. e. Mestel's disk | \f$ \Sigma \propto 1/r + \mathsf{noise} \f$
!!  constant temperature                           | \f$ 100\, \mathsf{K} \f$
!!  centrifugal balance                            | \f$ v_\varphi^2 = -r \partial_r \Phi \f$
!!
!! </div> <div class="col-md-6">
!! <img src="http://www.astrophysik.uni-kiel.de/fosite/sgagndisk_mbh1e7_md1e6_0960.png" class="img-fluid img-thumbnail" alt="column density">
!! You can find a [time-lapse movie] (agndisk.html) showing the temporal evolution of the
!! column density in the [gallery] (gallery.html).
!! </div> </div>
!!
!! References:
!! - \cite britsch2006 M. Britsch, "Gravitational instability and fragmentation of self-gravitating accretion disks",
!!      PhD thesis, 2006
!!
!! \example mestel.f90
!----------------------------------------------------------------------------!
PROGRAM mestel
  USE fosite_mod
#ifdef NECSXAURORA
    USE asl_unified
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER :: GN = 6.67384E-11       ! Newtons grav. constant [m^3/kg/s^2]
  REAL, PARAMETER :: YEAR = 3.15576E+7      ! year [s]
  REAL, PARAMETER :: PARSEC = 3.0857E+16    ! parsec [m]
  REAL, PARAMETER :: AU = 1.49598E+11       ! astronomical unit [m]
  REAL, PARAMETER :: MSUN = 1.989E+30       ! solar mass [kg]
  REAL, PARAMETER :: RG = 8.31              ! gas constant
  ! simulation parameters
  REAL, PARAMETER :: SIMTIME = 1.0E+4*YEAR  ! simulation time [s]
  REAL, PARAMETER :: MBH = 1.0E+7*MSUN      ! initial black hole mass [kg]
  REAL, PARAMETER :: MRATIO = 0.1           ! initial mdisk/mbh mass ratio
  REAL, PARAMETER :: MDISK = MRATIO*MBH     ! initial disk mass [kg]
  REAL, PARAMETER :: TEMP = 100.0           ! initial temperature [K]
  REAL, PARAMETER :: NOISE = 0.3            ! initial noise level
  ! physics settings
  INTEGER, PARAMETER :: PHYS = EULER        ! transport model
  REAL, PARAMETER :: MU = 6.02E-4           ! mean molecular weight [kg/mol]
  REAL, PARAMETER :: GAMMA = 1.4            ! ratio of specific heats
  REAL, PARAMETER :: BETA_VIS = 1.0E-3      ! beta viscosity parameter
  ! mesh settings
  REAL, PARAMETER :: RMIN = 5.0E-2 * PARSEC ! inner radius [m]
  REAL, PARAMETER :: RMAX = 1.E+0 * PARSEC  ! outer radius [m]
  REAL, PARAMETER :: RGEO = 1.0 * PARSEC    ! geometry scaling constant
  INTEGER, PARAMETER :: MGEO = LOGCYLINDRICAL ! mesh geometry
  INTEGER, PARAMETER :: XRES = 64         ! mesh resolution (radial)
  INTEGER, PARAMETER :: YRES = 128         ! mesh resolution (azimuthal)
  INTEGER, PARAMETER :: ZRES = 1            ! mesh resolution (z)
  ! output settings
  INTEGER, PARAMETER :: ONUM = 100         ! number of output time steps
  CHARACTER(LEN=256), PARAMETER :: &
  OFNAME = 'markward', &                    ! data file name
  ODIR   = "./"                             ! output directory
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)

#ifdef NECSXAURORA
CALL asl_library_initialize()
#endif

CALL Sim%InitFosite()
CALL MakeConfig(Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics,Sim%Fluxes,Sim%Sources)
CALL Sim%Run()
CALL Sim%Finalize()

#ifdef NECSXAURORA
CALL asl_library_finalize()
#endif

DEALLOCATE(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              sources,fluxes,grav,physics,&
                              vis,cooling
    !------------------------------------------------------------------------!
    physics => Dict( &
        "problem"         / PHYS, &
        "output/bccsound" / 1, &
        "mu"              / MU, &
        "gamma"           / GAMMA, &
        "units"           / SI)

    fluxes => Dict( &
        "fluxtype"        / KT, &
        "order"           / LINEAR, &
        "variables"       / PRIMITIVE, &
        "limiter"         / VANLEER, &
        "theta"           / 1.2)

    mesh => Dict( &
        "meshtype"        / MIDPOINT, &
        "geometry"        / MGEO, &
        "inum"            / XRES, &
        "jnum"            / YRES, &
        "knum"            / ZRES, &
        "xmin"            / LOG(RMIN/PARSEC), &
        "xmax"            / LOG(RMAX/PARSEC), &
        "ymin"            / (-PI), &
        "ymax"            / ( PI), &
        "zmin"            / 0.0, &
        "zmax"            / 0.0, &
        "gparam"          / RGEO, &
        "use_fargo"       / 1, &
        "fargo"           / 1, &
        "decomposition"   / (/-1,1,1/), &
        "output/volume"   / 1 )

    boundary => Dict( &
        "western"         / CUSTOM, &
        "eastern"         / CUSTOM, &
        "southern"        / PERIODIC, &
        "northern"        / PERIODIC, &
        "bottomer"        / REFLECTING, &
        "topper"          / REFLECTING)

    grav => Dict( &
        "stype"           / GRAVITY, &
        "cvis"            / 0.9, &
        "energy"          / 0, &
        "output/height"   / 1, &
        "self/gtype"      / SPECTRAL, &
        "self/green"      / 1, &
        "pmass/gtype"     / POINTMASS, &
        "pmass/potential" / NEWTON, &
        "pmass/mass"      / MBH, &
        "pmass/outbound"  / 1)

    ! cooling model
    cooling => Dict( &
        "stype"           / DISK_COOLING, &
!         "method"          / GRAY, &
        "method"          / GAMMIE, &
        "b_cool"          / 10.0, &
        "Tmin"            / 3.0, &
        "rhomin"          / 1.0E-30, &
        "cvis"            / 0.1)

    ! viscosity model
    vis => Dict( &
        "stype"           / VISCOSITY, &
        "vismodel"        / BETA, &
        "dynconst"        / BETA_VIS, &
        "output/stress"   / 1, &
        "output/dynvis"   / 1, &
        "output/kinvis"   / 1, &
        "cvis"            / 0.5)

    ! collect source terms
    sources => Dict( &
        "diskcooling"     / cooling, &
!         "viscosity"       / vis, &
        "grav"            / grav)

    ! time discretization settings
    timedisc => Dict( &
        "method"          / SSPRK, &
        "tol_rel"         / 1.0E-3, &
        "cfl"             / 0.3, &
        "stoptime"        / SIMTIME, &
        "dtlimit"         / 1.0E-4, &
        "tmin"            / 3.0, &
        "rhstype"         / 1, &
        "maxiter"         / 2000000000)

    ! add absolute error bounds and output fields depending on physics
    SELECT CASE(PHYS)
    CASE(EULER)
       CALL SetAttr(timedisc, "tol_abs", (/1.0E-3, 1.0E-1, 1.0E-01, 1.0E+10/))
       CALL SetAttr(timedisc, "output/xmomentum", 0)
       CALL SetAttr(timedisc, "output/ymomentum", 0)
       CALL SetAttr(timedisc, "output/energy", 0)
       CALL SetAttr(timedisc, "output/rhs", 0)
       CALL SetAttr(timedisc, "output/external_sources", 0)
    CASE DEFAULT
       CALL Sim%Error("MakeConfig","Physics model not supported.")
    END SELECT

    ! data file settings
    datafile => Dict( &
        "fileformat"      / VTK, &
        "filename"        / "mestel", &
        "count"           / ONUM)

    ! create configuration
    config => Dict( &
        "physics"         / physics, &
        "fluxes"          / fluxes, &
        "mesh"            / mesh, &
        "boundary"        / boundary, &
        "sources"         / sources, &
        "timedisc"        / timedisc, &
        "datafile"        / datafile)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Timedisc,Mesh,Physics,Fluxes,Sources)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    USE sources_gravity_mod, ONLY : sources_gravity
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp => null()
    CLASS(sources_gravity), POINTER :: gp => null()
    CLASS(marray_base), ALLOCATABLE :: rands,Sigma
#ifdef PARALLEL
    INTEGER :: ierror
#endif
    REAL    :: mass
    CHARACTER(LEN=20) :: mdisk_str
#ifdef NECSXAURORA
    INTEGER :: rng, n
#endif
    !------------------------------------------------------------------------!
    ! get gravitational acceleration
    IF (ALLOCATED(Sources)) &
      sp => Sources%GetSourcesPointer(GRAVITY)
    IF (.NOT.ASSOCIATED(sp)) &
      CALL Physics%Error("mestel::InitData","no gravity term initialized")

    ! get random numbers for density noise
    ALLOCATE(rands,Sigma)
    rands = marray_base()
#ifndef NECSXAURORA
    CALL InitRandSeed(Physics)
    CALL RANDOM_NUMBER(rands%data1d)
#else
    CALL asl_random_create(rng, ASL_RANDOMMETHOD_MT19937_64)
    CALL asl_random_distribute_uniform(rng)
    n = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)
    CALL asl_random_generate_d(rng, n, rands%data1d)
    CALL asl_random_destroy(rng)
#endif
    rands%data1d(:) = rands%data1d(:) * NOISE * 2.0 + (1.0 - NOISE)

    ! set surface density using radial power law (1/r) with a little noise
    Sigma = marray_base()
    Sigma%data1d(:) = rands%data1d(:) * RMIN/Mesh%radius%data2d(:,2)

    ! determine disk mass
    mass = SUM(Mesh%volume%data1d(:)*Sigma%data1d(:),MASK=Mesh%without_ghost_zones%mask1d(:))
#ifdef PARALLEL
    CALL MPI_AllReduce(MPI_IN_PLACE,mass,1,DEFAULT_MPI_REAL,MPI_SUM, &
            Mesh%comm_cart,ierror)
#endif

    ! setting for custom boundary conditions (western boundary)
    SELECT TYPE(bwest => Timedisc%Boundary%boundary(WEST)%p)
    CLASS IS (boundary_custom)
      CALL bwest%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER,CUSTOM_NOGRAD/))
    END SELECT

    ! setting for custom boundary conditions (eastern boundary)
    SELECT TYPE(beast => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS (boundary_custom)
      CALL beast%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_REFLECT,CUSTOM_REFLECT,CUSTOM_LOGEXPOL,CUSTOM_REFLECT/))
    END SELECT

    SELECT TYPE (pvar => Timedisc%pvar)
    CLASS IS(statevector_euler)
      ! 1. set surface density by rescaling Sigma
      pvar%density%data1d(:) = Sigma%data1d(:) * MDISK / mass
      ! 2. initialize gravitational acceleration
      CALL gp%UpdateGravity(Mesh,Physics,Fluxes,pvar,0.0,0.0)
      ! 3. set azimuthal velocity: balance initial radial gravitational
      !    acceleration with centrifugal acceleration
      pvar%velocity%data4d(:,:,:,1:Physics%VDIM) = &
        Timedisc%GetCentrifugalVelocity(Mesh,Physics,Fluxes,Sources,(/0.,0.,1./),gp%accel%data4d)
      ! 4. transform velocities to rotating frame
      ! ATTENTION: this works only if the second velocity is the azimuthal velocity
      IF (Mesh%OMEGA.GT.0.0) &
        pvar%velocity%data2d(:,2) = pvar%velocity%data2d(:,2) - Mesh%OMEGA*Mesh%radius%data2d(:,2)
      ! 5. set pressure using surface density and initial temperature
      pvar%pressure%data1d(:) = Physics%constants%RG/Physics%MU * TEMP * pvar%density%data1d(:)
    CLASS DEFAULT
      CALL Physics%Error("mestel::InitData","unsupported state vector")
    END SELECT

    IF (Mesh%fargo%GetType().EQ.2) &
       Timedisc%w(:,:) = SQRT(Physics%constants%GN*(MBH/Mesh%radius%bcenter(:,Mesh%JMIN,:)))

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    ! print some information
    WRITE (mdisk_str, '(ES10.2)') mdisk/MSUN
    CALL Mesh%Info(" DATA-----> initial condition: Mestel's disk")
    CALL Mesh%Info("            disk mass:         " // TRIM(mdisk_str) // " M_sun")

    DEALLOCATE(rands,Sigma)
  END SUBROUTINE InitData

  SUBROUTINE InitRandSeed(Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: Physics
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    !------------------------------------------------------------------------!
    ! Initialize random number generator with a seed based on the systems time
    ! source: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
#ifdef PARALLEL
    seed = seed + Physics%GetRank()
#endif
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
  END SUBROUTINE InitRandSeed
END PROGRAM mestel
