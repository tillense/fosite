!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sgdisk.f90                                                        #
!#                                                                           #
!# Copyright (C) 2010-2019                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> Program and data initialization a self-gravitating disk (sgdisk)
!!
!! \author Jannes Klee
!!
!! \example sgdisk.f90
!----------------------------------------------------------------------------!
PROGRAM sgdisk
  USE fosite_mod
#ifdef NECSXAURORA
    USE asl_unified
#endif
#ifdef HAVE_KROME
  USE krome_main
  USE krome_user
#endif
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
  REAL,     PARAMETER :: GN      = 6.6742D-8      ! grav. constant (cgs)     !
  REAL,     PARAMETER :: RG      = 8.31447D+7     ! molar gas constant (cgs) !
  REAL,     PARAMETER :: MSUN    = 1.989D+33      ! solar mass [g]           !
  REAL,     PARAMETER :: AU      = 1.49597870691E+13 ! ast. unit [cm]        !
  REAL,     PARAMETER :: YEAR    = 3.15576E+7     ! Julian year [s]          !
  REAL,     PARAMETER :: MSCALE  = MSUN           ! mass scaling param.      !
  ! simulation parameters
  REAL,     PARAMETER :: TSIM    = 6.0            ! simulation time [ORP]    !
  REAL,     PARAMETER :: MU      = 2.35D+0        ! mean mol. mass [g/mol]  !
  REAL,     PARAMETER :: VALPHA  = 0.05           ! alpha visc. parameter    !
  ! disk & central object
  REAL,     PARAMETER :: MBH     = 1.0*MSCALE     ! central mass             !
  REAL,     PARAMETER :: MDISK   = 0.1*MSCALE     ! initial disk mass        !
  REAL,     PARAMETER :: T_INIT  = 120.0          ! initial temperatur       !
  REAL,     PARAMETER :: T_MIN   = 10.0            ! minimal temperatur       !
  ! cooling parameters
  REAL,     PARAMETER :: BETA_C  = 10.0           ! cooling parameter        !
  REAL,     PARAMETER :: GAMMA   = 1.6            ! adiabatic index          !
  ! mesh settings
  INTEGER,  PARAMETER :: MGEO    = LOGCYLINDRICAL ! geometry of the mesh     !
  REAL,     PARAMETER :: GPAR    = AU             ! geom. scal. parameter    !
  REAL,     PARAMETER :: RMIN    = 1.0*GPAR       ! inner radius of the disk !
  REAL,     PARAMETER :: RMAX    = 25.0*GPAR      ! outer radius of the grid !
  INTEGER,  PARAMETER :: XRES    = 256            ! x-resolution             !
  INTEGER,  PARAMETER :: YRES    = XRES*3          ! y-resolution             !
  INTEGER,  PARAMETER :: ZRES    = 1              ! z-resolution             !
  ! mestel
  LOGICAL             :: DAMP                     ! wave damping             !
  REAL, PARAMETER     :: NOISE   = 0.3            ! initial noise level      !
  ! output file parameter
  INTEGER, PARAMETER  :: ONUM    = 100            ! number of outputs        !
  CHARACTER(LEN=256), PARAMETER &
                      :: ODIR    = "./"           ! output directory         !
  CHARACTER(LEN=256)  :: OFNAME  = 'sgdisk'       ! output name              !
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  REAL                :: CSISO,OMEGA,ORP
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)

#ifdef NECSXAURORA
CALL asl_library_initialize()
#endif

CALL Sim%InitFosite()
CALL MakeConfig(Sim, Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc, Sim%Fluxes, Sim%Sources, &
              Sim%Timedisc%pvar, Sim%Timedisc%cvar)
CALL Sim%Run()
CALL Sim%Finalize()

#ifdef NECSXAURORA
    CALL asl_library_finalize()
#endif

DEALLOCATE(Sim)

CONTAINS

   SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: Sim
    TYPE(Dict_TYP), POINTER      :: config
    TYPE(Dict_TYP), POINTER      :: mesh, physics, fluxes, boundary, grav,&
                                    sources, timedisc, datafile, cooling, &
                                    rotframe, pmass, self, damping
    !------------------------------------------------------------------------!
    ! some derived simulation parameters
    CSISO   = SQRT(RG/MU*T_INIT)
    OMEGA   = 0.0
    ORP     = 2.*PI/(SQRT(GN*MBH/RMAX**3.))
    DAMP    = .FALSE.

    ! mesh settings
    ! stellar orbits must be inside the central hole of the mesh
    mesh => Dict(&
                "meshtype"        / MIDPOINT, &
                "geometry"        / MGEO, &
                "inum"            / XRES, &
                "jnum"            / YRES, &
                "knum"            / ZRES, &
                "xmin"            / LOG(RMIN/GPAR), &
                "xmax"            / LOG(RMAX/GPAR), &
                "ymin"            / 0.0, &
                "ymax"            / (2.0*PI), &
                "zmin"            / 0.0, &
                "zmax"            / 0.0, &
                "use_fargo"       / 1, &
                "fargo"           / 2, &
                "decomposition"   / (/ -1, 1, 1/), &
!                "output/radius"   / 1, &
!                "output/dl"       / 1, &
                "output/position_vector" / 1, &
                "gparam"          / GPAR)


    ! physics settings
    physics => Dict(&
                "problem"         / EULER, &
                "mu"              / MU, &
                "gamma"           / GAMMA, &
                "units"           / CGS)
!                "output/bccsound" / 1)

    ! boundary conditions
    boundary => Dict(&
                "western"         / CUSTOM,&
                "eastern"         / CUSTOM,&
                "southern"        / PERIODIC, &
                "northern"        / PERIODIC, &
                "bottomer"        / REFLECTING, &
                "topper"          / REFLECTING)


    ! numerical fluxes and reconstruction method
    fluxes => Dict(&
                "order"           / LINEAR, &
                "fluxtype"        / KT, &
                "variables"       / PRIMITIVE, &
                "passive_limiting" / .FALSE., &
                "limiter"         / VANLEER, &
                "theta"           / 1.2)


    ! viscosity source term
    rotframe => Dict(&
                "stype"           / ROTATING_FRAME)


    ! gravitational acceleration due to binary system
    pmass => Dict(&
                "gtype"           / POINTMASS, &
                "mass"            / MBH, &
                "output/position" / 1)

    ! cooling function
    cooling =>  Dict(&
                "stype"           / DISK_COOLING, &
                "output/Qcool"    / 1, &
                "rhomin"          / 1.0E-30, &
                "method"          / GAMMIE, &
                "b_cool"          / BETA_C)

    ! self-gravity
    self => Dict(&
                "gtype"           / SPECTRAL, &
                "output/potential" / 1, &
                "self/green"      / 1)

    ! collect all gravity-source terms
    grav => Dict(&
                "stype"             / GRAVITY, &
                "pmass"             / pmass, &
                "self"              / self, &
                "output/potential"  / 1, &
                "energy"            / 0, &
                "output/height"     / 1, &
                "output/accel"      / 1)

    ! collect all source terms
    sources => Dict(&
                "cooling"           / cooling, &
                "grav"              / grav)

    ! add wave damping if requested
    IF (DAMP) &
       CALL SetAttr(sources,"damping",damping)

    ! time discretization settings
    timedisc => Dict(&
                "method"            / MODIFIED_EULER, &
                "cfl"               / 0.3, &
                "dtlimit"           / 1.0E-50, &
                "tol_rel"           / 1.0E-3, &
                "rhstype"           / 1, &
                "tol_abs"           / (/ 1.0E-16, 1.0, 1.0E-16, 1.0 /), &
                "output/energy"     / 1, &
                "output/rhs"        / 1, &
                "stoptime"          / (ORP*TSIM), &
                "pmin"              /  1e-14, &
                "tmin"              / T_MIN, &
                "checkdata"         / IOR(CHECK_NOTHING,CHECK_TMIN), &
                "maxiter"           / 100000000)

    ! initialize data input/output
    datafile => Dict(&
               "fileformat"         / XDMF, &
               "filename"           / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"              / ONUM)

    config => Dict(&
                "physics"           / physics, &
                "fluxes"            / fluxes, &
                "mesh"              / mesh, &
                "boundary"          / boundary, &
                "sources"           / sources, &
                "timedisc"          / timedisc, &
                "datafile"          / datafile)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc,Fluxes,Sources,pvar,cvar)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(sources_base),  POINTER       :: Sources
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp => null()
    INTEGER           :: i,j,k
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    REAL              :: mass
    REAL, DIMENSION(:,:,:), POINTER :: r, Sigma
    REAL, DIMENSION(:,:,:,:), POINTER :: r_vec
    CHARACTER(LEN=32) :: mdisk_str
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                      :: rands
#ifdef NECSXAURORA
    INTEGER :: rng, n
#endif
    !------------------------------------------------------------------------!
    ! MESTEL'S DISK
    ! distance from origin to cell bary centers and position vector
    r => Mesh%RemapBounds(Mesh%radius%bcenter(:,:,:))
    r_vec => Mesh%RemapBounds(Mesh%posvec%bcenter(:,:,:,:))
    ! pointer to density array
    Sigma => Mesh%RemapBounds(pvar%data4d(:,:,:,Physics%DENSITY))

    ! set surface density using radial power law (1/r) with a little noise
#ifndef NECSXAURORA
    CALL InitRandSeed(Timedisc)
    CALL RANDOM_NUMBER(rands)
#else
    CALL asl_random_create(rng, ASL_RANDOMMETHOD_MT19937_64)
    CALL asl_random_distribute_uniform(rng)
    n = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)
    CALL asl_random_generate_d(rng, n, rands)
#endif
    rands = rands * NOISE * 2.0 + (1.0 - NOISE)
    Sigma = rands*(RMIN/r(:,:,:))

#ifdef NECSXAURORA
    CALL asl_random_destroy(rng)
#endif

    ! determine disk mass
    mass = SUM(Mesh%volume%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) * &
               Sigma(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
#ifdef PARALLEL
       CALL MPI_AllReduce(MPI_IN_PLACE,mass,1,DEFAULT_MPI_REAL,MPI_SUM, &
            Mesh%comm_cart,ierror)
#endif
    ! rescale disk mass
    Sigma(:,:,:) = Sigma(:,:,:) * MDISK / mass

    ! 2. azimuthal velocity: balance initial radial acceleration with centrifugal acceleration
    SELECT TYPE(phys => Physics)
    TYPE IS(physics_euler)
      pvar%data4d(:,:,:,Physics%PRESSURE) = &
      CSISO*CSISO*pvar%data4d(:,:,:,Physics%DENSITY)/GAMMA
      ! necessary in order to calculte height below in UpdateGravity
      SELECT TYPE(pvarious => pvar)
      CLASS IS(statevector_euler)
        CALL phys%UpdateSoundSpeed(pvarious)
      END SELECT
    END SELECT

    sp => Sources
    DO
      IF (.NOT.ASSOCIATED(sp)) EXIT 
      SELECT TYPE(sp)
      CLASS IS(sources_gravity)
        gp => sp
        EXIT
      END SELECT
      sp => sp%next
    END DO

    CALL gp%UpdateGravity(Mesh,Physics,Fluxes,pvar,0.0,0.0)

    ! reset velocities
    pvar%data2d(:,Physics%XVELOCITY:Physics%YVELOCITY) = 0.0
    ! ATTENTION: Don't use GetCentrifugalVelocity without the optional acceleration array!
    ! This would yield undefined data, because GetCentrifugalVelocity calls ComputeRHS
    ! which calls CenterBoundary. Since the FARFIELD boundary conditions are not
    ! initialized at this stage (see below), the result is undefined.
    pvar%data4d(:,:,:,Physics%XVELOCITY:Physics%XVELOCITY+Physics%VDIM-1) = &
        Timedisc%GetCentrifugalVelocity(Mesh,Physics,Fluxes,Sources,(/0.,0.,1./),gp%accel%data4d)

    IF (Mesh%fargo%GetType().EQ.2) &
       Timedisc%w(:,:) = SQRT(Physics%constants%GN*MBH/r(:,Mesh%JMIN,:))-Mesh%OMEGA*r(:,Mesh%JMIN,:)

    ! transform velocities to rotating frame
    pvar%data4d(:,:,:,Physics%YVELOCITY) = pvar%data4d(:,:,:,Physics%YVELOCITY) &
        - Mesh%OMEGA*r(:,:,:)

    ! get conservative variables
    CALL Physics%Convert2Conservative(pvar,cvar)

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
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER,CUSTOM_NOGRAD/))
    END SELECT

    CALL Physics%Convert2Conservative(pvar,cvar)
    ! print some information
    WRITE (mdisk_str, '(ES8.2)') mdisk
    CALL Mesh%Info(" DATA-----> initial condition: Mestel's disk")
    CALL Mesh%Info("            disk mass:         " // TRIM(mdisk_str))

  END SUBROUTINE InitData


  SUBROUTINE InitRandSeed(Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base),INTENT(IN) :: Timedisc
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
    seed = seed + Timedisc%GetRank()
#endif
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
  END SUBROUTINE InitRandSeed


END PROGRAM sgdisk
