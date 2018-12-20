!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: shearingbox.f90                                                   #
!#                                                                           #
!# Copyright (C) 2015                                                        #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!> Program and data initialization for a shearing-box simulation
!! References:
!! [1] Charles F. Gammie (2001). "Nonlinear outcome of gravitational
!!     instability in cooling, gaseous disks"
!!     The Astrophysical Journal 553: 174-183
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE mesh_generic
  USE physics_generic
  USE fluxes_generic
  USE reconstruction_generic
  USE boundary_generic
  USE sources_generic
  USE gravity_generic
  USE timedisc_generic
  USE fileio_generic
  USE common_dict
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
  ! simulation parameter
  REAL, PARAMETER    :: OMEGA      = 1.             ! rotation at fid. point !
  REAL, PARAMETER    :: TSIM       = 100./OMEGA     ! simulation time        !
  REAL, PARAMETER    :: GAMMA      = 2.0            ! ratio of spec. heats   !
!  REAL, PARAMETER    :: GAMMA      = 1.3            ! ratio of spec. heats   !
  ! not yet ordered parameter
  REAL, PARAMETER    :: GN         = 6.67384E-11    ! grav. constant [SI]    !
  REAL, PARAMETER    :: SIGMA0     = 1.             ! mean surf.dens.        !
  REAL, PARAMETER    :: DELSIGMA   = SIGMA0*5e-4    ! disturbance            !

!  REAL, PARAMETER    :: SIGMA0     = 1.0            ! mean surface density   !
!  REAL, PARAMETER    :: P0         = 1e-20          ! ambient pressure       !
!  REAL, PARAMETER    :: AMP        = 1e3            ! amplitude of the pulse !
!  REAL, PARAMETER    :: SIGMAWIDTH = 0.06           ! width of gauss         !

  INTEGER, PARAMETER :: GREEN      = 1              ! 1 = razor-sharp disk   !
  REAL, PARAMETER    :: TAU_C      = 10.0/OMEGA
  REAL, PARAMETER    :: B_COOL     = TAU_C

  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN

! Linear Theory Test
  INTEGER, PARAMETER :: XRES       = 1024           ! amount of cells in x-  !
  INTEGER, PARAMETER :: YRES       = 1024           ! y-direction (rho/phi)  !
  REAL, PARAMETER    :: DOMAINX    = 320*GN*SIGMA0/(OMEGA*OMEGA)! domain size !
  REAL, PARAMETER    :: DOMAINY    = 320*GN*SIGMA0/(OMEGA*OMEGA)! domain size !

  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter 1.5 !
                                                    ! for kep. disk          !
  REAL, PARAMETER    :: SOUNDSPEED = PI*GN*SIGMA0/OMEGA ! det. by. Toomre   !
  REAL, PARAMETER    :: XMIN       = -0.5*DOMAINX
  REAL, PARAMETER    :: XMAX       = +0.5*DOMAINX
  REAL, PARAMETER    :: YMIN       = -0.5*DOMAINY
  REAL, PARAMETER    :: YMAX       = +0.5*DOMAINY
  REAL, PARAMETER    :: SHEARVAL   = Q*OMEGA
                                                    ! shear value, necessary !
                                                    ! for boundaries         !
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM       = 100
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"!/sfs/fs3/work-sh1/supas295/fosite/tests/lintheory2/'
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sbox"
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(config)

CALL SetupFosite(Sim)

! set initial condition
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: config
    !--------------------------------------------------------------------------!
    ! local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh,physics,fluxes,boundary,&
                               grav,vis,cooling,shearingbox,sources,timedisc,&
                               datafile,&
                               logfile,reconstruction,selfgravity !TODO: these are not used until now.
    !--------------------------------------------------------------------------!
    ! mesh settings
    mesh =>     Dict(&
                "meshtype"    / MIDPOINT, &
                "geometry"    / MGEO, &
                "inum"        / XRES, &
                "jnum"        / YRES, &
                "xmin"        / XMIN, &
                "xmax"        / XMAX, &
                "ymin"        / YMIN, &
                "ymax"        / YMAX, &
                "output/rotation" / 1, &
                "output/volume"   / 1, &
                "output/bh"   / 1, &
                "output/dl"   / 1  &
                )

    ! physics settings
    physics =>  Dict(&
                "problem"     / EULER2D, &
                "gamma"       / GAMMA &
!                "cs"          / SOUNDSPEED &
                )

    ! fluxes settings
    fluxes =>   Dict(&
                "order"       / LINEAR, &
                "fluxtype"    / KT, &
                "variables"   / PRIMITIVE, &
                "limiter"     / VANLEER, &
                "theta"       / 2.0 &
                )

    ! boundary conditions
    boundary => Dict(&
                "western"     / SHEARING, &
                "eastern"     / SHEARING, &
                "southern"    / PERIODIC, &
                "northern"    / PERIODIC, &
                "shearval"    / SHEARVAL &
                )

    ! gravity settings (source term)
    grav =>     Dict(&
                "stype"               / GRAVITY, &
                "self/gtype"          / SBOXSPECTRAL, &
                "self/green"          / GREEN, &
                "self/omega_rot"      / OMEGA, &
                "output/accel"        / 1, &
                "self/output/phi"     / 1, &
                "self/output/accel_x" / 1, &
                "self/output/accel_y" / 1, &
                "self/Q"              / Q &
                )

    ! viscosity settings (source term)
    vis =>      Dict(&
                "stype"       / VISCOSITY, &
                "vismodel"    / ALPHA, &
                "cvis"        / 0.3 &
                )

    ! cooling settings (source term)
    cooling => Dict(&
               "stype"        / DISK_COOLING, &
               "method"       / GAMMIE_SB, &
               "b_cool"       / B_COOL &
               )


    ! shearing box fictious forces
    shearingbox => Dict(&
                "stype"           / SHEARBOX, &
                "output/accel_x"  / 1, &
                "output/accel_y"  / 1, &
                "omega"           / OMEGA &
                )

    ! sources settings (contains source terms)
    sources =>  Dict(&
                "grav"        / grav, &
!                "viscosity"   / vis, &
                "cooling"     / cooling, &
                "shearing"    / shearingbox &
                )

    ! time discretization settings
    timedisc => Dict(&
                "method"      / SSPRK, &
                "cfl"         / 0.4, &
                "stoptime"    / TSIM, &
                "dtlimit"     / 1e-30, &
                "maxiter"     / 1000000, &
                "tol_rel"     / 1.0E-3 &
                )

    ! data i/o settings
    datafile => Dict(&
                "fileformat"  / XDMF, &
                "unit"        / 5555, &
                "filepath"    / TRIM(ODIR), &
                "filename"    / TRIM(OFNAME), &
                "count"       / ONUM &
                )

    ! overall config settings
    config =>   Dict(&
                "mesh"        / mesh, &
                "physics"     / physics, &
                "fluxes"      / fluxes, &
                "boundary"    / boundary, &
                "sources"     / sources, &
                "timedisc"    / timedisc, &
                "datafile"    / datafile &
                )
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Timedisc,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! local variable declaration
    ! TODO: no variables yet
!    REAL, PARAMETER   :: K
    REAL              :: ylen, kx, ky, xlen
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX+200) :: rands2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: rands
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: K
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Physics,Timedisc
    !------------------------------------------------------------------------!
    !------------------------ standard run ----------------------------------!
    Timedisc%pvar(:,:,Physics%DENSITY)    = SIGMA0!*(1.+rands2)
    ! fullfills Q = 1
    IF (Physics%PRESSURE.GT.0) Timedisc%pvar(:,:,Physics%PRESSURE) = &
              2.5**2*PI*PI*GN*GN*SIGMA0**3./(GAMMA*OMEGA**2)

    ! create random numbers
    CALL InitRandSeed(Timedisc)
    CALL RANDOM_NUMBER(rands2)

    rands2 = (rands2-0.5)*0.05
    rands = SQRT(Timedisc%pvar(:,:,Physics%PRESSURE)/ &
                 Timedisc%pvar(:,:,Physics%DENSITY))*rands2(:,Mesh%JGMIN+200:Mesh%JGMAX)
!    rands = SOUNDSPEED*rands2(:,Mesh%JGMIN+200:Mesh%JGMAX)

    ! initial velocity
    Timedisc%pvar(:,:,Physics%XVELOCITY)  = rands(:,:)

    CALL RANDOM_NUMBER(rands2)
    rands2 = (rands2-0.5)*0.05
    rands = SQRT(Timedisc%pvar(:,:,Physics%PRESSURE)/ &
                 Timedisc%pvar(:,:,Physics%DENSITY))*rands2(:,Mesh%JGMIN+200:Mesh%JGMAX)
!    rands = SOUNDSPEED*rands2(:,Mesh%JGMIN+200:Mesh%JGMAX)
    Timedisc%pvar(:,:,Physics%YVELOCITY)  = -Q*OMEGA*Mesh%bcenter(:,:,1) + &
                                            rands(:,:)



    !------------------------------------------------------------------------!
!    !-------------------------- epicyclic motion ----------------------------!
!!   References:
!!   [2] Sijme-Jan Paardekooper (2012). "Numerical convergence in self-
!!       gravitating shearing sheet simulations and the sochastic nature of
!!       disc fragmentation"
!!       Mon.Not.R.Astron.Soc.421: 3286-3299
!
!    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.1*SOUNDSPEED
!    Timedisc%pvar(:,:,Physics%YVELOCITY) = -Q*Omega*Mesh%bcenter(:,:,1)
!    Timedisc%pvar(:,:,Physics%DENSITY)   = SIGMA0
!
!    !------------------------------------------------------------------------!
!    !------------------------ two blobs test --------------------------------!
!    ! Test for self-gravity
!    ! Set OMEGA = 0, use Euler2D and only(!) the gravity source module
!
!    ! initial velocity
!    Timedisc%pvar(:,:,Physics%XVELOCITY)  = 0.0
!    Timedisc%pvar(:,:,Physics%YVELOCITY)  = 0.0
!
!!     xlen = ABS(Mesh%xmax-Mesh%xmin)
!!     print *, Mesh%dx
!!     WHERE ((Mesh%bcenter(:,:,1).GT.(Mesh%xmin+0.49*xlen)).AND. &
!!          (Mesh%bcenter(:,:,1).LT.(Mesh%xmin+0.51*xlen)).AND. &
!!          (Mesh%bcenter(:,:,2).GT.(Mesh%xmin+0.49*xlen)).AND. &
!!          (Mesh%bcenter(:,:,2).LT.(Mesh%xmin+0.51*xlen)))
!!        Timedisc%pvar(:,:,Physics%DENSITY) = AMP
!!     END WHERE
!
!    ! initial density
!    FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
!       Timedisc%pvar(i,j,Physics%DENSITY) = SIGMA0 + &
!            AMP*EXP(-((SQRT((Mesh%bcenter(i,j,1))**2+Mesh%bcenter(i,j,2)**2))&
!            /SIGMAWIDTH)**2)
!!            AMP*EXP(-&
!!            ((SQRT((Mesh%bcenter(i,j,1)-0.1*DOMAINX)**2+Mesh%bcenter(i,j,2)**2))&
!!            /SIGMAWIDTH)**2) + &
!!            AMP*EXP(-&
!!            ((SQRT((Mesh%bcenter(i,j,1)+0.3*DOMAINX)**2+&
!!                   (Mesh%bcenter(i,j,2)+0.3*DOMAINY)**2))&
!!            /SIGMAWIDTH)**2)
!    END FORALL
!
!    ! initial pressure (-> very low value)
!    Timedisc%pvar(:,:,Physics%PRESSURE)   = P0
!    !------------------------------------------------------------------------!
!    !---------------------- linear theory test ------------------------------!
!!   References:
!!   [2] Charles F. Gammie (1996). "Linear Theory of Magnetized, Viscous, Self-
!!       gravitating Gas Disks"
!!       The Astrophysical Journal 462: 725-731
!!    There is also an python script solving the ode of the paper for the
!!    inviscid, unmagnetized case (!>\todo{not uploaded yet})
!
!    ! initial velocity
!    Timedisc%pvar(:,:,Physics%XVELOCITY)  = 0.0
!    Timedisc%pvar(:,:,Physics%YVELOCITY)  = -Q*OMEGA*Mesh%bcenter(:,:,1)
!
!    ! initial density
!    kx = -2*(2*PI/DOMAINX)
!    ky = 2*PI/DOMAINY
!    Timedisc%pvar(:,:,Physics%DENSITY)    = SIGMA0 + DELSIGMA*COS( &
!                                            kx*Mesh%bcenter(:,:,1) + &
!                                            ky*Mesh%bcenter(:,:,2))
!
!    K(:,:) = (PI*GN/OMEGA)**(2.)*Timedisc%pvar(:,:,Physics%DENSITY)**(3.-2.)
!
!    ! initial pressure (determined by Toomre-criterion)
!    IF (Physics%PRESSURE.GT.0) Timedisc%pvar(:,:,Physics%PRESSURE)   = &
!!               1e-20
!!              PI*PI*GN*GN*Timedisc%pvar(:,:,Physics%DENSITY)**3./(GAMMA*OMEGA**2) !> \todo{make clear which line should be used here}
!              PI*PI*GN*GN*SIGMA0**3./(GAMMA*OMEGA**2)
!!               K*Timedisc%pvar(:,:,Physics%DENSITY)**(GAMMA)
!    !------------------------------------------------------------------------!
!    !-------------------------- square test ---------------------------------!
!    ! Test for shearing and boundary module
!
!    Timedisc%pvar(:,:,Physics%XVELOCITY)  = 0.0
!    Timedisc%pvar(:,:,Physics%YVELOCITY)  = 0.0
!    Timedisc%pvar(:,:,Physics%PRESSURE)   = 1e-20
!
!    WHERE ((ABS(Mesh%bcenter(:,:,1)).LT.0.15*DOMAINX.AND. &
!            ABS(Mesh%bcenter(:,:,2)).LT.0.15*DOMAINY))
!      Timedisc%pvar(:,:,Physics%DENSITY)  = SIGMA0
!    ELSEWHERE
!      Timedisc%pvar(:,:,Physics%DENSITY)  = SIGMA0*1e-5
!    END WHERE
!    !------------------------------------------------------------------------!


    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Shearing-Box")
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
