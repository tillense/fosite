!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: shearingsheet.f90                                                 #
!#                                                                           #
!# Copyright (C) 2015-2018                                                   #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> Program and data initialization for a shearing-box simulation with
!! standard initialization
!!
!! \author Jannes Klee
!!
!! \example shearingsheet.f90
!!
!! This program runs the self-gravitating shearingsheet simulation first
!! done by Gammie (2001) \cite gammie2001 . For initialization a constant
!! background density \f$ \Sigma \f$ is used and the pressure is set
!! in a way to fullfill the Toomre-criterion \f$ Q = 1 \f$. The initial
!! velocities\f$ v_x, v_y \f$ are randomly distributed with sub-sonic values.
!! Whether you use a fast cooling with \f$ \beta = 2.0 \f$ or a slow cooling
!! \f$ \beta = 10.0 \f$, you see fragmentation or a settling into the
!! gravitoturbulent state.
!!
!! <div class="row"> <div class="col-md-6">
!!  Simulation parameters         ||
!!  ------------------            | -----------------
!!  cooling parameter \f$ \beta \f$| \f$ 10 \text{ (slow)}, 2 \text{ (fast)} \f$
!!  resolution \f$ N_x \times N_y \f$ | \f$ 1024 \times 1024 \f$
!!  flux solver                   | \f$ \mathtt{KT} \f$
!!  time-discretization           | \f$ \mathtt{DORMAND-PRINCE} \f$
!!  limiter                       | \f$ \mathtt{VANLEER} \f$
!!  boundaries (north, south)     | shearing
!!  boundaries (east, west)       | periodic
!!  heat capacity ratio \f$ \gamma \f$ | \f$ 2.0 \f$
!!
!!  Initial condition       | |
!!  ------------------      | -----------------
!!  density \f$ \Sigma \f$  | \f$ 1.0 \f$
!!  pressure \f$ P \f$      | \f$ \frac{2.5^2 \pi^2 G \Sigma^3}{\gamma \Omega^2} \f$
!!  velocity \f$ v_x \f$ | \f$ q \Omega y + \delta v_x \f$
!!  velocity \f$ v_y \f$ | \f$ \delta v_y \f$
!!  random velocities \f$ \delta v_x, \delta v_y \f$ | \f$ < c_{\mathrm{s}} \f$
!! </div> <div class="col-md-6">
!!   <img src="http://www.astrophysik.uni-kiel.de/fosite/1024_beta10_vanleer_gamma2.png" class="img-fluid img-thumbnail" alt="slow cooling">
!!   <img src="http://www.astrophysik.uni-kiel.de/fosite/1024_beta2_vanleer_gamma2.png" class="img-fluid img-thumbnail" alt="fastcooling">
!! </div> </div>
!! You can find some [movies] (shearingsheet.html) showing the temporal evolution of the
!! column density in the [gallery] (gallery.html).
!!
!! \n
!!
!! \attention Long runtime.
!!
!! References:
!! - \cite gammie2001 Charles F. Gammie (2001). "Nonlinear outcome of
!!     gravitational instability in cooling, gaseous disks"
!!     The Astrophysical Journal 553: 174-183
!----------------------------------------------------------------------------!
PROGRAM shearingsheet
  USE fosite_mod
#ifdef NECSXAURORA
  USE asl_unified
#endif
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPI_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: GN         = 1.0            ! grav. constant [GEOM]  !
  INTEGER, PARAMETER :: UNITS      = GEOMETRICAL    ! needs to be consistent !
  ! simulation parameter
  REAL, PARAMETER    :: OMEGA      = 1.0            ! rotation at fid. point !
  REAL, PARAMETER    :: SIGMA0     = 1.0            ! mean surf.dens.        !
  REAL, PARAMETER    :: TSIM       = 100./OMEGA     ! simulation time        !
  REAL, PARAMETER    :: GAMMA      = 2.0            ! dep. on vert. struct.  !
  REAL, PARAMETER    :: BETA_C     = 10.0           ! cooling parameter      !
!  REAL, PARAMETER    :: BETA_C     = 2.0           ! 2 -> collapse          !
  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter     !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN
  INTEGER, PARAMETER :: SHEAR_DIRECTION = 1         ! enables shearingsheet with direcion !
  INTEGER, PARAMETER :: XRES       = 1024           ! cells in x-direction   !
  INTEGER, PARAMETER :: YRES       = 1024           ! cells in y-direction   !
  INTEGER, PARAMETER :: ZRES       = 1              ! cells in z-direction   !
  REAL               :: DOMAINX    = 320.0          ! domain size [GEOM]     !
  REAL               :: DOMAINY    = 320.0          ! domain size [GEOM]     !
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM       = 100
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "shearingsheet"
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)

CALL asl_library_initialize()

CALL Sim%InitFosite()
CALL MakeConfig(Sim, Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
CALL Sim%Run()
CALL Sim%Finalize()


CALL asl_library_finalize()

DEALLOCATE(Sim)

CONTAINS
  SUBROUTINE MakeConfig(Sim,config)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP), POINTER :: config
    !--------------------------------------------------------------------------!
    ! local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh,physics,fluxes,grav,cooling, &
                               shearingbox,sources,timedisc,datafile
    REAL :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX, SOUNDSPEED
    !--------------------------------------------------------------------------!
    DOMAINX    = DOMAINX*GN*SIGMA0/(OMEGA*OMEGA)
    DOMAINY    = DOMAINY*GN*SIGMA0/(OMEGA*OMEGA)
    XMIN       = -0.5*DOMAINX
    XMAX       = +0.5*DOMAINX
    YMIN       = -0.5*DOMAINY
    YMAX       = +0.5*DOMAINY
    ZMIN       = 0.0
    ZMAX       = 0.0
    SOUNDSPEED = PI*GN*SIGMA0/OMEGA ! Toomre-criterion

    ! physics settings
    physics =>  Dict(&
                "problem"     / EULER, &
                "gamma"       / GAMMA, &
                "units"       / GEOMETRICAL &
                )

    ! mesh settings
    mesh =>     Dict(&
                "meshtype"    / MIDPOINT, &
                "geometry"    / MGEO, &
                "omega"       / OMEGA, &
                "shear_dir"   / SHEAR_DIRECTION, &
                "decomposition" / (/ 1, -1, 1/), &
                "inum"        / XRES, &
                "jnum"        / YRES, &
                "knum"        / ZRES, &
                "xmin"        / XMIN, &
                "xmax"        / XMAX, &
                "ymin"        / YMIN, &
                "ymax"        / YMAX, &
                "zmin"        / ZMIN, &
                "zmax"        / ZMAX, &
                "Q"           / Q &
                )

    ! fluxes settings
    fluxes =>   Dict(&
                "order"       / LINEAR, &
                "fluxtype"    / KT, &
                "variables"   / PRIMITIVE, &
                "limiter"     / VANLEER &
                )

    ! gravity settings (source term)
    grav =>     Dict(&
                "stype"               / GRAVITY, &
                "self/gtype"          / SBOXSPECTRAL, &
                "self/output/phi"     / 0, &
                "self/output/accel_x" / 0, &
                "self/output/Fmass2D" / 0, &
                "self/output/accel_y" / 0 &
                )

    ! parametrized cooling from Gammie (2001)
    cooling =>  Dict(&
                "stype"        / DISK_COOLING, &
                "method"       / GAMMIE_SB, &
                "output/Qcool" / 1, &
                "b_cool"       / BETA_C &
                )

    ! shearing box fictious forces
    shearingbox => Dict(&
                "stype"        / SHEARBOX &
                )

    ! sources settings (contains source terms)
    sources =>  Dict(&
                "cooling"     / cooling, &
                "shearing"    / shearingbox, &
                "grav"        / grav &
                )

    ! time discretization settings
    timedisc => Dict(&
                "method"      / DORMAND_PRINCE, &
                "cfl"         / 0.4, &
                "stoptime"    / TSIM, &
                "dtlimit"     / 1e-40, &
                "output/external_sources" / 0, &
                "output/density_cvar" / 0, &
                "output/energy" / 0, &
                "output/xmomentum" / 0, &
                "output/ymomentum" / 0, &
                "output/rhs" / 0, &
                "output/geometrical_sources" / 0, &
                "maxiter"     / 100000000, &
                "tol_rel"     / 1.0E-3 &
                )

    ! data i/o settings
    datafile => Dict(&
                "fileformat"  / XDMF, &
                "filepath"    / TRIM(ODIR), &
                "filename"    / TRIM(OFNAME), &
                "count"       / ONUM &
                )

    ! overall config settings
    config =>   Dict(&
                "mesh"        / mesh, &
                "physics"     / physics, &
                "fluxes"      / fluxes, &
                "sources"     / sources, &
                "timedisc"    / timedisc, &
                "datafile"    / datafile &
                )
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics, pvar, cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),    INTENT(IN)  :: Mesh
    CLASS(physics_base), INTENT(IN)  :: Physics
    CLASS(marray_compound), POINTER, INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    REAL              :: SOUNDSPEED
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                      :: rands
#ifdef NECSXAURORA
    INTEGER :: rng, n
#endif
    !------------------------ standard run ----------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      ! constant initial density
      p%density%data3d(:,:,:) = SIGMA0

      ! constant initial pressure determined by Q = 1
      SOUNDSPEED = PI*Physics%Constants%GN*SIGMA0/OMEGA ! Toomre-criterion
      p%pressure%data3d(:,:,:) = 2.5**2*PI*PI* &
                          Physics%Constants%GN**2.*SIGMA0**3./(GAMMA*OMEGA**2)

      ! create random numbers for setup of initial velocities
#ifndef NECSXAURORA
    CALL InitRandSeed(Physics)
#else
    CALL asl_random_create(rng, ASL_RANDOMMETHOD_MT19937_64)
    CALL asl_random_distribute_uniform(rng)
    n = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)
#endif
      IF (Mesh%shear_dir.EQ.2) THEN
#ifndef NECSXAURORA
        CALL RANDOM_NUMBER(rands)
#else
        CALL asl_random_generate_d(rng, n, rands)
#endif
        rands = (rands-0.5)*0.1*SOUNDSPEED
        p%velocity%data4d(:,:,:,1) = rands(:,:,:)

#ifndef NECSXAURORA
        CALL RANDOM_NUMBER(rands)
#else
        CALL asl_random_generate_d(rng, n, rands)
#endif
        rands = (rands-0.5)*0.1*SOUNDSPEED
        p%velocity%data4d(:,:,:,2)  = -Q*OMEGA*Mesh%bcenter(:,:,:,1) + rands(:,:,:)
      ELSE IF (Mesh%shear_dir.EQ.1) THEN
#ifndef NECSXAURORA
        CALL RANDOM_NUMBER(rands)
#else
        CALL asl_random_generate_d(rng, n, rands)
#endif
        rands = (rands-0.5)*0.1*SOUNDSPEED
        p%velocity%data4d(:,:,:,1) = Q*OMEGA*Mesh%bcenter(:,:,:,2) + rands(:,:,:)

#ifndef NECSXAURORA
        CALL RANDOM_NUMBER(rands)
#else
        CALL asl_random_generate_d(rng, n, rands)
#endif
        rands = (rands-0.5)*0.1*SOUNDSPEED
        p%velocity%data4d(:,:,:,2)  = rands(:,:,:)
      END IF
    CLASS DEFAULT
      CALL Physics%Error("shear::InitData","only non-isothermal HD supported")
    END SELECT
#ifdef NECSXAURORA
    CALL asl_random_destroy(rng)
#endif
    !------------------------------------------------------------------------!
    CALL Physics%Convert2Conservative(pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Standard run shearingsheet")
  END SUBROUTINE InitData

  !> random number generator
  SUBROUTINE InitRandSeed(Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),INTENT(IN) :: Physics
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
END PROGRAM shearingsheet
