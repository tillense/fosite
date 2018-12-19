!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: chemshock1d.f03                                                   #
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
!> 1D Chemical Shock
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
#ifdef HAVE_KROME
  USE krome_main
  USE krome_user
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  !! Some explanation to the setup. For a given central mass and a radius an
  !! angular velocity velocity results, which needs to be calculated by hand
  !! because SQRT is not available on some compilers on NEC machines.
  !!
  !! The Toomre criterion needs to be fullfilled with Q=1. This gives a
  !! limitation between density and pressure.
  ! general constants
  INTEGER, PARAMETER :: UNITS      = CGS            ! needs to be consistent !
  ! simulation parameter
  REAL, PARAMETER    :: MSUN       = 1.9884e33      ! mass of the sun [cgs]  !
  REAL, PARAMETER    :: PARSEC     = 3.0856776e18   ! parsec in [cgs]        !
  REAL, PARAMETER    :: GN         = 6.6742e-8      ! grav. constant [cgs]   !
  REAL, PARAMETER    :: AU         = 149.5978707e11 ! astron. units [cgs]    !
  ! simulation parameter
  REAL, PARAMETER    :: MBH        = 1.0e4*MSUN     ! m. of cent. obj. [cgs] !
!  REAL, PARAMETER    :: R0         = 100*AU         ! rad. dist. of SB [cgs] !
  REAL, PARAMETER    :: R0         = 300*AU         ! rad. dist. of SB [cgs] !
  ! ----- OMEGA = SQRT(GN*MBH/R0**3)                ! needs to be fullfilled !
!  REAL, PARAMETER    :: OMEGA0     = 1.99e-8        ! ang. velocity [cgs]    !
  REAL, PARAMETER    :: OMEGA0     = 3.832e-9        ! ang. velocity [cgs]    !
  REAL, PARAMETER    :: TSIM       = 100.0/OMEGA0   ! simulation time [cgs]  !
  REAL, PARAMETER    :: GAMMA      = 1.6            ! dep. on vert. struct.  !
  REAL, PARAMETER    :: GAMMA3D    = 1.6            ! ratio of spec. heats   !
  REAL, PARAMETER    :: BETA_C     = 10.0           ! cooling parameter      !
  REAL, PARAMETER    :: BETA_START = 10.0           ! cooling parameter      !
  REAL               :: T_START    = 0.0
  REAL               :: DT_BDEC    = 0.0/OMEGA0
  ! not yet ordered parameter
  REAL, PARAMETER    :: MU         = 2.3            ! mean mol. mass [cgs]   !
!  REAL, PARAMETER    :: TEMP       = 30.           ! Temp. in K. [cgs]      !
  ! ----- CS = SQRT(GAMMA3D*Physics%Constants%RG*TEMP/Physics%mu)
!  REAL, PARAMETER    :: CS         = 416.556       ! speed of sound [cgs]   !
!  REAL, PARAMETER    :: SIGMA0     = 15.0e3        ! surface density [cgs]  !
  REAL, PARAMETER    :: SIGMAMAX   = R0*OMEGA0*OMEGA0/(2.0*PI*GN) ! this value is only for comparison in order to ensure a KGS
  REAL, PARAMETER    :: SIGMA0     = SIGMAMAX*4e-4
!  REAL, PARAMETER    :: SIGMA0     = 50000.0
  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter     !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN
  INTEGER, PARAMETER :: XRES       = 32            ! amount of cells in x-  !
  INTEGER, PARAMETER :: YRES       = 32            ! y-direction (rho/phi)  !
  REAL               :: DOMAINX    = 320.0          ! domain size            !
  REAL               :: DOMAINY    = 320.0          ! domain size            !
  ! fargo 0=off, 3=on (for SB)
  INTEGER, PARAMETER :: FARGO      = 3              ! 3 = Shearingbox        !
  ! number of output time step
  INTEGER, PARAMETER :: ONUM       = 100
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sbox"
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!
  print *, SIGMA0

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
                               grav,chemics,cooling,shearingbox,sources,timedisc,&
                               datafile
    REAL :: XMIN,XMAX,YMIN,YMAX
    !--------------------------------------------------------------------------!
    DOMAINX    = DOMAINX*GN*SIGMA0/(OMEGA0*OMEGA0)
    DOMAINY    = DOMAINY*GN*SIGMA0/(OMEGA0*OMEGA0)
    XMIN       = -0.5*DOMAINX
    XMAX       = +0.5*DOMAINX
    YMIN       = -0.5*DOMAINY
    YMAX       = +0.5*DOMAINY

    ! physics settings
    physics =>  Dict(&
                "problem"     / EULER3D, &
                "gamma"       / GAMMA, &
                "mu"          / MU, &
                "units"       / UNITS &
                )

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
                "zmin"        / ZMIN, &
                "zmax"        / ZMAX, &
                "omega"       / OMEGA0, &
                "output/rotation" / 0, &
                "output/volume"   / 0, &
                "output/bh"   / 0, &
                "output/dl"   / 0  &
                )

    ! fluxes settings
    fluxes =>   Dict(&
                "order"       / LINEAR, &
                "fluxtype"    / KT, &
                "variables"   / PRIMITIVE, &
                "limiter"     / VANLEER &
                )

    ! boundary conditions
    boundary => Dict(&
                "western"     / PERIODIC, &
                "eastern"     / PERIODIC, &
                "southern"    / SHEARING, &
                "northern"    / SHEARING, &
                "fargo"       / FARGO &
                )

    ! gravity settings (source term)
    grav =>     Dict(&
                "stype"               / GRAVITY, &
                "self/gtype"          / SBOXSPECTRAL, &
                "self/fargo"          / FARGO, &
                "output/accel"        / 0, &
                "self/output/phi"     / 1, &
                "self/output/accel_x" / 0, &
                "self/output/accel_y" / 0, &
                "self/Q"              / Q, &
                "output/height"       / 0 &
                )

    ! cooling settings (source term)
    cooling =>  Dict(&
                "stype"        / DISK_COOLING, &
                "method"       / GAMMIE_SB, &
!                "Tmin"         / 10.0, &
                "b_cool"       / BETA_C, &
                "cvis"         / 0.01 &
                )

    ! shearing box fictious forces
    shearingbox => Dict(&
                "stype"           / SHEARBOX, &
                "fargo"           / FARGO, &
                "output/accel_x"  / 0, &
                "output/accel_y"  / 0 &
                )

    ! krome chemistry library
    chemics =>  Dict(&
                "stype"       / CHEMISTRY, &
                "num_species" / 12 &
                )

    ! sources settings (contains source terms)
    sources =>  Dict(&
                "grav"        / grav, &
                "chemics"     / chemics, &
                "cooling"     / cooling, &
                "shearing"    / shearingbox &
                )

    ! time discretization settings
    timedisc => Dict(&
                "method"      / DORMAND_PRINCE, &
                "cfl"         / 0.4, &
                "stoptime"    / TSIM, &
                "fargo"       / FARGO, &
                "dtlimit"     / 1e-40, &
                "maxiter"     / 100000000, &
                "tol_rel"     / 1.0E-2, &
!                "tol_abs"     / (/ 1.0E-16, 1., 1.0E-16 /) &
                "checkdata"   / IOR(CHECK_NOTHING,CHECK_TMIN), &
                "tmin"        / 2.7 &
!                "pmin"        / 1.E-50 &
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
    REAL              :: ylen, kx, ky, xlen
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: rands2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: rands
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: K
    INTEGER           :: i,j
    REAL              :: CS,CS2
    TYPE(Sources_TYP), POINTER :: chemics
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Physics,Timedisc
    !----------------------- realistic cooling-------------------------------!
    CS = PI*GN*SIGMA0/OMEGA0
    CS2 = CS*CS
    ! fullfills Toomre
    Timedisc%pvar(:,:,Physics%DENSITY) = SIGMA0

    ! with ideal gas law
    Timedisc%pvar(:,:,Physics%PRESSURE) = 2.5**2*PI*PI* &
                        Physics%Constants%GN**2.*SIGMA0**3./(GAMMA*OMEGA0**2)

    ! with ideal gas law
    Timedisc%pvar(:,:,Physics%PRESSURE) = &
          Timedisc%pvar(:,:,Physics%DENSITY)*CS2/GAMMA


    ! random velocities (slower than speed of sound)
    CALL InitRandSeed(Timedisc)
    CALL RANDOM_NUMBER(rands)
    rands = CS*(rands-0.5)
    CALL RANDOM_NUMBER(rands2)
    rands2 = CS*(rands2-0.5)

    IF (Mesh%WE_shear) THEN
      Timedisc%pvar(:,:,Physics%XVELOCITY)  = rands(:,:)
      Timedisc%pvar(:,:,Physics%YVELOCITY)  = -Q*OMEGA0*Mesh%bcenter(:,:,1) + &
                                              rands2(:,:)
    ELSE IF (Mesh%SN_shear) THEN
      Timedisc%pvar(:,:,Physics%XVELOCITY)  = Q*OMEGA0*Mesh%bcenter(:,:,2) + &
                                              rands(:,:)
      Timedisc%pvar(:,:,Physics%YVELOCITY)  = rands2(:,:)
    END IF

    ! initialize chemistry
    !------------------------------------------------------------------------!
#ifdef HAVE_KROME
    Timedisc%pvar(:,:,:,Physics%VNUM+1:Physics%VNUM+Physics%PNUM) = 0e0
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_H)     = 0.9225    !H
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_E)     = 1.0d-4    !E
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_Hj)    = 1.0d-4    !H+
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_D)     = 1.0d-20   !D
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_Dj)    = 1.0d-20   !D+
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_HE)    = 0.0972    !He
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_HEj)   = 1.0d-20   !He+
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_H2j)   = 1.0d-20   !H2+
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_H2)    = 1.0d-5    !H2
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_HD)    = 1.0d-8    !HD
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_Hk)    = 1.0d-20   !H-
    Timedisc%pvar(:,:,:,Physics%VNUM+krome_idx_HEjj)  = 1.0d-20   !He++

    DO k=Mesh%KMIN,Mesh%KMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          Timedisc%pvar(i,j,k,Physics%VNUM+1:Physics%VNUM+Physics%PNUM) = &
            Timedisc%pvar(i,j,k,Physics%VNUM+1:Physics%VNUM+Physics%PNUM) / &
            sum(Timedisc%pvar(i,j,k,Physics%VNUM+1:Physics%VNUM+Physics%PNUM))
        END DO
      END DO
    END DO
#endif

    !------------------------------------------------------------------------!

    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
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
