!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: markward.f90                                                      #
!#                                                                           #
!# Copyright (C) 2012-2014                                                   #
!# Manuel Jung    <mjung@astrophysik.uni-kiel.de>                            #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> markward standard disk (self gravitating accretion disk around a SMBH)
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE gravity_generic
  USE timedisc_generic
  USE common_dict
  USE sources_rotframe, ONLY : Convert2RotatingFrame_rotframe
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER :: GN = 6.67384E-11       ! Newtons grav. constant [m^3/kg/s^2]
  REAL, PARAMETER :: YEAR = 3.15576E+7      ! year [s]
  REAL, PARAMETER :: PARSEC = 3.0857E+16    ! parsec [m]
  REAL, PARAMETER :: AU = 1.49598E+11       ! astronomical unit [m]
  REAL, PARAMETER :: MSUN = 1.989E+30       ! solar mass [kg]
  REAL, PARAMETER :: RG = 8.31              ! gas constant
  ! simulation parameters
  REAL, PARAMETER :: SIMTIME = 1.0E+2*YEAR  ! simulation time [s]
  REAL, PARAMETER :: MBH = 1.0E+7*MSUN      ! initial black hole mass [kg]
  REAL, PARAMETER :: MRATIO = 0.1           ! initial mbh/mdisk mass ratio
  REAL, PARAMETER :: MDISK = MRATIO*MBH     ! initial disk mass [kg]
  REAL, PARAMETER :: TEMP = 100.0           ! initial temperature [K]
  REAL, PARAMETER :: NOISE = 0.3            ! initial noise level
  ! physics settings
  INTEGER, PARAMETER :: PHYS = EULER2D    ! transport model 
!   INTEGER, PARAMETER :: PHYS = EULER2D_SGS
  REAL, PARAMETER :: MU = 6.02E-4           ! mean molecular weight [kg/mol]
  REAL, PARAMETER :: GAMMA = 1.4            ! ratio of specific heats
  REAL, PARAMETER :: BETA_VIS = 1.0E-3      ! beta viscosity parameter
  ! mesh settings
  REAL, PARAMETER :: RMIN = 5.0E-2 * PARSEC ! inner radius [m]
  REAL, PARAMETER :: RMAX = 1.E+0 * PARSEC  ! outer radius [m]
  REAL, PARAMETER :: RGEO = 0.1 * PARSEC    ! geometry scaling constant
  INTEGER, PARAMETER :: MGEO = LOGPOLAR     ! mesh geometry
!   INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 128          ! mesh resolution (radial)
  INTEGER, PARAMETER :: YRES = 256          ! mesh resolution (azimuthal)
  ! output settings
  INTEGER, PARAMETER :: ONUM = 100          ! number of output time steps
  CHARACTER(LEN=256), PARAMETER :: &
  OFNAME = 'markward', &                    ! data file name
  ODIR   = "./"                             ! output directory
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)          :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(Sim%config)

CALL SetupFosite(Sim)

! set initial condition
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics,Sim%Fluxes)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              sources,fluxes,grav,physics,rotframe,&
                              vis,cooling
    !------------------------------------------------------------------------!
    physics => Dict( &
        "problem" / PHYS, &
        "output/bccsound" / 1, &
        "mu" / MU, &
        "gamma" / GAMMA, &
        "units" / SI)

    fluxes => Dict( &
        "fluxtype"  / KT, &
        "order" / LINEAR, &
        "variables" / PRIMITIVE, &
        "limiter" / VANLEER, &
        "theta" / 1.2)

    mesh => Dict( &
        "meshtype" / MIDPOINT, &
        "geometry" / MGEO, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "xmin" / LOG(RMIN/PARSEC), &
        "xmax" / LOG(RMAX/PARSEC), &
        "ymin" / (-PI), &
        "ymax" / ( PI), &
        "gparam" / RGEO, &
        "decomposition" / (/-1,1/), &
        "output/volume" / 1 )

    boundary => Dict( &
        "western" / CUSTOM, &
        "eastern" / CUSTOM, &
!         "eastern" / NO_GRADIENTS, &
!        "eastern" / FARFIELD, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC)

    grav => Dict("stype" / GRAVITY, &
                 "cvis" / 0.9, &
        "output/height" / 1, &
        "pmass/gtype" / POINTMASS, & 
        "pmass/potential" / NEWTON, &
        "pmass/mass" / MBH, &
        "pmass/outbound" / 1 , &             ! enables accretion
        "self/gtype" / SPECTRAL, &
        "self/green" / 1)

    ! cooling model
    cooling => Dict("stype" / DISK_COOLING, &
!                     "method" / GRAY, &
                    "method" / GAMMIE, &
                    "b_cool" / 10.0, &
                    "Tmin" / 30.0, &        ! minimum midplane temperature
                    "rhomin" / 1.0E-30, &   ! minimum midplane density
                    "cvis" / 0.1)
    
    ! viscosity model
    IF (PHYS.EQ.EULER2D) THEN
       ! alpha/beta viscosity
       vis => Dict( &
          "stype" / VISCOSITY, &
          "vismodel" / BETA, &
          "dynconst" / BETA_VIS, &
          "output/stress" / 1, &
          "output/dynvis" / 1, &
          "output/kinvis" / 1, &
          "cvis" / 0.5)
    ELSE
       ! SGS viscosity source term
       vis => Dict(&
            "stype" / SGS,&
            "output/dynvis" / 1, &
            "output/kinvis" / 1, &
            "cvis" / 0.5)
    END IF

    ! collect source terms
    sources => Dict( &
       "diskcooling" / cooling, &
!        "viscosity" / vis, &
       "grav" / grav)

    ! time discretization settings
    timedisc => Dict( &
        "method" / SSPRK, &        ! time integration method
!         "fargo" / 0, &
!         "rhstype" / 1, &
        "tol_rel" / 1.0E-2, &      ! relative error
        "cfl" / 0.3, &             ! CFL number
        "stoptime" / SIMTIME, &
        "dtlimit" / 1.0E+2, &
        "maxiter" / 2000000000)

    ! add absolute error bounds and output fields depending on physics
    SELECT CASE(PHYS)
    CASE(EULER2D)
       CALL SetAttr(timedisc, "tol_abs", (/0.0, 1.0E-3, 0.0, 1.0E+06/))
       CALL SetAttr(timedisc, "output/xmomentum", 0)
       CALL SetAttr(timedisc, "output/ymomentum", 0)
       CALL SetAttr(timedisc, "output/energy", 0)
       CALL SetAttr(timedisc, "output/rhs", 0)
       CALL SetAttr(timedisc, "output/external_sources", 0)      
    CASE(EULER2D_SGS)
       CALL SetAttr(timedisc, "tol_abs", (/0.0, 1.0E-3, 0.0, 1.0E+03, 1.0E-3/))
       CALL SetAttr(timedisc, "output/xmomentum", 0)
       CALL SetAttr(timedisc, "output/ymomentum", 0)
       CALL SetAttr(timedisc, "output/energy", 0)
       CALL SetAttr(timedisc, "output/rhs", 0)
       CALL SetAttr(timedisc, "output/external_sources", 0)
    CASE DEFAULT
       CALL Error(Sim,"MakeConfig","Physics model not supported.")
    END SELECT

    ! data file settings
    datafile => Dict( &
        "fileformat" / XDMF, &
        "filename" / "markward", &
        "count" / ONUM)

    ! create configuration
    config => Dict( &
        "physics" / physics, &
        "fluxes" / fluxes, &
        "mesh" / mesh, &
        "boundary" / boundary, &
        "sources" / sources, &
        "timedisc" / timedisc, &
        "datafile" / datafile)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Timedisc,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Sources_TYP), POINTER :: sp
    INTEGER           :: i,j
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    REAL              :: mass,ephi(2),cs,vphi
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: rands
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: accel
    REAL,DIMENSION(:,:), POINTER :: r,Sigma
    REAL, DIMENSION(:,:,:), POINTER :: r_vec
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc, Physics,Fluxes
    !------------------------------------------------------------------------!
    ! distance from origin to cell bary centers and position vector
    r => RemapBounds(Mesh,Mesh%radius%bcenter(:,:))
    r_vec => RemapBounds(Mesh,Mesh%posvec%bcenter(:,:,:))
    ! pointer to density array
    Sigma => RemapBounds(Mesh,Timedisc%pvar(:,:,Physics%DENSITY))

    ! set surface density using radial power law (1/r) with a little noise
    CALL InitRandSeed(Timedisc)
    CALL RANDOM_NUMBER(rands)
    rands = rands * NOISE * 2.0 + (1.0 - NOISE)
    Sigma = rands*(RMIN/r(:,:))

    ! determine disk mass
    mass = SUM(Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) * Sigma(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
#ifdef PARALLEL
       CALL MPI_AllReduce(MPI_IN_PLACE,mass,1,DEFAULT_MPI_REAL,MPI_SUM, &
            Mesh%comm_cart,ierror)
#endif
    ! rescale disk mass
    Sigma(:,:) = Sigma(:,:) * MDISK / mass

    ! set pressure using surface density and initial temperature
    Timedisc%pvar(:,:,Physics%PRESSURE) = Physics%constants%RG/Physics%mu * TEMP * Sigma(:,:)

    ! set initial SGS pressure to very low value
    IF(GetType(Physics) .EQ. EULER2D_SGS)&
        Timedisc%pvar(:,:,Physics%SGSPRESSURE) = 1.0E-9*Timedisc%pvar(:,:,Physics%PRESSURE)

    ! reset velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = 0.0

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    ! set the velocity due to the centrifugal force
    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%XVELOCITY+Physics%DIM) = &
          GetCentrifugalVelocity(Timedisc,Mesh,Physics,Fluxes,(/0.,0.,1./))

    ! setting for custom boundary conditions (western boundary)
    IF(GetType(Timedisc%Boundary(WEST)).EQ.CUSTOM) THEN
      Timedisc%boundary(WEST)%cbtype(:,Physics%DENSITY)   = CUSTOM_NOGRAD
      Timedisc%boundary(WEST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_NOGRAD
      Timedisc%boundary(WEST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
      Timedisc%boundary(WEST)%cbtype(:,Physics%PRESSURE)  = CUSTOM_NOGRAD
      IF(GetType(Physics) .EQ. EULER2D_SGS)&
                Timedisc%boundary(WEST)%cbtype(:,Physics%SGSPRESSURE) = CUSTOM_NOGRAD
    END IF

    ! setting for custom boundary conditions (eastern boundary)
    IF(GetType(Timedisc%Boundary(EAST)).EQ.CUSTOM) THEN
      Timedisc%boundary(EAST)%cbtype(:,Physics%DENSITY)   = CUSTOM_NOGRAD
      Timedisc%boundary(EAST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_NOGRAD
      Timedisc%boundary(EAST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
      Timedisc%boundary(EAST)%cbtype(:,Physics%PRESSURE)  = CUSTOM_NOGRAD
      IF(GetType(Physics) .EQ. EULER2D_SGS)&
                Timedisc%boundary(EAST)%cbtype(:,Physics%SGSPRESSURE) = CUSTOM_NOGRAD
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE InitData

  SUBROUTINE InitRandSeed(Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP),INTENT(IN) :: Timedisc
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
    seed = seed + GetRank(Timedisc)
#endif
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
  END SUBROUTINE InitRandSeed
END PROGRAM Init
