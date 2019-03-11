!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: pringle3d.f90                                                     #
!#                                                                           #
!# Copyright (C) 2008-2019                                                   #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
!> Jim Pringles famous viscous spreading ring solution in 3D
!!
!! \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \example pringle3d.f90
!!
!! This variant simulates the ring with vertical extent with and without
!! rotational symmetry.
!! The equations are solved in non-dimensional units, using the Keplerian
!! velocity at the location of the initial ring at R=1 as the velocity scale.
!!
!! \remark (1) This test is very sensitive to mesh settings. If there are
!!         not enough mesh points at small radii instabilities grow
!!         starting at the inner boundary.
!! \remark (2) There is an unsolved stability issue with the viscous
!!         source term. The stability condition doesn't seem to be
!!         sufficient. Thus in case of high Reynolds numbers one has
!!         to adjust the viscous CFL number cvis to enforce smaller time
!!         steps and preserve stability.
!!
!! References:
!!  - \cite lynden1974 Lynden-Bell, D. & Pringle, J. E. Evolution of Viscous Disks and Origin
!!      of Nebular Variables, M.N.R.A.S vol. 168(3) (1974) pp. 603-637
!!      http://adsabs.harvard.edu/abs/1974MNRAS.168..603L
!!  - \cite pringle1981 Pringle, J. E. Accretion discs in astrophysics
!!      Annual review of astronomy and astrophysics. vol. 19 (1981) pp. 137-162
!!      DOI: 10.1146/annurev.aa.19.090181.001033
!!  - \cite speith2003 Speith, R. and Kley, W. Stability of the viscously spreading ring
!!      Astron. Astrophys. vol. 399 (2003), pp. 395-407
!!      DOI: 10.1051/0004-6361:20021783
!!
!----------------------------------------------------------------------------!
PROGRAM pringle3d
  USE fosite_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 5.0E-3     ! simulation time [TVIS] see below
                                             ! should be < 1 ~ viscous time scale
  ! these are the basic parameters; for stable solutions one requires
  ! MA >> 1 and RE > MA
  REAL, PARAMETER    :: RE      = 1.0E+4     ! Reynolds number (at R=1)
  REAL, PARAMETER    :: MA      = 1.0E+2     ! Mach number (at R=1)
  REAL, PARAMETER    :: MDISK   = 1.0E-2     ! disk mass << 1 = central mass
  ! lower limit for nitial density
  REAL, PARAMETER    :: RHOMIN  = 1.0E-8    ! minimal initial density
  ! viscosity prescription
  INTEGER, PARAMETER :: VISTYPE = BETA
!   INTEGER, PARAMETER :: VISTYPE = POWERLAW
  REAL, PARAMETER    :: PL_EXP  = 0.0        ! exponent for power law viscosity
  REAL, PARAMETER    :: TVIS    = 4./3.*RE   ! viscous time scale
  REAL, PARAMETER    :: TINIT   = 1.0E-3     ! time for initial condition [TVIS]
  ! mesh settings
!   INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!   INTEGER, PARAMETER :: MGEO = LOGCYLINDRICAL
  INTEGER, PARAMETER :: MGEO = SPHERICAL !!! ATTENTION: not applicable in 1D
  INTEGER, PARAMETER :: XRES = 200        ! x-resolution
  INTEGER, PARAMETER :: YRES = 16           ! y-resolution
  INTEGER, PARAMETER :: ZRES = 6           ! z-resolution
  REAL, PARAMETER    :: RMIN = 0.1        ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 2.0         ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'pringle3d'
  ! derived parameters
  REAL, PARAMETER    :: CSISO = 1./MA      ! isothermal speed of sound
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim,Sim%Mesh, Sim%Physics, Sim%Timedisc)
  CALL Sim%Run()
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)          :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               grav, vis, timedisc, fluxes, sources
    REAL              :: x1,x2,y1,y2,z1,z2,cvis
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CYLINDRICAL)
      x1 = RMIN
      x2 = RMAX
      y1 = -PI
      y2 = PI
!       z1 = 0.0
!       z2 = 0.0
      z1 = -0.1
      z2 = 0.1
!       bc(WEST)  = AXIS
      bc(WEST)  = NO_GRADIENTS
      bc(EAST)  = NO_GRADIENTS
!       bc(WEST)  = ABSORBING
!       bc(EAST)  = ABSORBING
!       bc(WEST)  = CUSTOM
!       bc(EAST)  = CUSTOM
      bc(SOUTH) = PERIODIC
      bc(NORTH) = PERIODIC
!       bc(BOTTOM) = ABSORBING
!       bc(TOP)    = ABSORBING
      bc(BOTTOM) = REFLECTING
      bc(TOP)    = REFLECTING
      !!! ATTENTION:  empirical formula found for standard parameters
      !!!             seems to work for RMIN=0.05 and XRES <= 1000
      cvis = (XRES/1207.0)**1.9
    CASE(LOGCYLINDRICAL)
      x1 = LOG(RMIN/GPAR)
      x2 = LOG(RMAX/GPAR)
      y1 = -PI
      y2 = PI
      z1 = 0.0
      z2 = 0.0
      bc(WEST)  = NO_GRADIENTS
      bc(EAST)  = NO_GRADIENTS
      bc(SOUTH) = PERIODIC
      bc(NORTH) = PERIODIC
      bc(BOTTOM) = NO_GRADIENTS
      bc(TOP)    = NO_GRADIENTS
      cvis = 0.1
    CASE(SPHERICAL)
      x1 = RMIN
      x2 = RMAX
      y1 = PI/2 - 0.02*PI
      y2 = PI/2 + 0.02*PI
      z1 = -PI
      z2 = PI
!       bc(WEST)  = ABSORBING
      bc(EAST)  = ABSORBING
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
      bc(WEST)  = CUSTOM
!       bc(EAST)  = CUSTOM
!       bc(SOUTH) = NO_GRADIENTS
!       bc(NORTH) = NO_GRADIENTS
      bc(SOUTH) = ABSORBING
      bc(NORTH) = ABSORBING
      bc(BOTTOM) = PERIODIC
      bc(TOP)    = PERIODIC
      !!! ATTENTION:  empirical formula found for standard parameters
      !!!             seems to work for RMIN=0.05 and XRES <= 1000
      cvis = (XRES/1207.0)**1.9
    CASE DEFAULT
      CALL Sim%Error("pringle::MakeConfig","mesh geometry not supported for Pringle disk")
    END SELECT

    ! mesh settings
    mesh => Dict( &
           "meshtype"   / MIDPOINT, &
           "geometry"   / MGEO, &
           "inum"       / XRES, &
           "jnum"       / YRES, &
           "knum"       / ZRES, &
           "xmin"       / x1, &
           "xmax"       / x2, &
           "ymin"       / y1, &
           "ymax"       / y2, &
           "zmin"       / z1, &
           "zmax"       / z2, &
!            "use_fargo"  / 1, &
!            "fargo"      / 2, &
           "gparam"     / GPAR)

    ! boundary conditions
    boundary => Dict( &
            "western"   / bc(WEST), &
            "eastern"   / bc(EAST), &
            "southern"  / bc(SOUTH), &
            "northern"  / bc(NORTH), &
            "bottomer"  / bc(BOTTOM), &
            "topper"    / bc(TOP))


    ! physics settings
    physics => Dict( &
            "problem"   / EULER_ISOTHERM, &
            "units"     / GEOMETRICAL, &
            "cs"        / CSISO)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
            "order"     / LINEAR, &
            "fluxtype"  / KT, &
            "variables" / PRIMITIVE, &
            "limiter"   / VANLEER, &
            "theta"     / 1.2)

    ! source term due to viscosity
    vis => Dict( &
            "stype"     / VISCOSITY, &
            "vismodel"  / VISTYPE, &
            "dynconst"  / (1./RE), &
            "exponent"  / PL_EXP, &
            "cvis"      / cvis)

    ! source term due to a point mass
    grav => Dict( &
            "stype"     / GRAVITY, &
            "pmass/gtype"/ POINTMASS, &    ! grav. accel. of a point mass
            "pmass/mass" / 1.0, &          ! non-dim. mass of the accreting object
            "pmass/outbound" / 0)          ! disable accretion

    ! combine all source terms
    sources => Dict( &
            "viscosity" / vis, &
            "grav"      / grav)

    ! time discretization settings
    timedisc => Dict( &
            "method"    / SSPRK, &
            "stoptime"  / (TSIM * TVIS), &
            "cfl"       / 0.4, &
            "dtlimit"   / (1E-10 * TVIS), &
!             "rhstype"   / 1, &
            "maxiter"   / 100000000, &
            "tol_rel"   / 0.01, &
            "tol_abs"   / (/1e-40,1e-5,1e-8,1e-5/))

!     ! enable output of analytical solution
! #ifdef HAVE_FGSL
!     IF (MGEO.EQ.CYLINDRICAL) &
!       CALL SetAttr(timedisc,"output/solution",1)
! #endif

    ! initialize data input/output
    datafile => Dict( &
            "fileformat"/ VTK, &
            "filename"  / (TRIM(ODIR) // TRIM(OFNAME)), &
            "count"     / ONUM)

    config => Dict( &
            "mesh"      / mesh, &
            "physics"   / physics, &
            "boundary"  / boundary, &
            "fluxes"    / fluxes, &
            "sources"   / sources, &
            "timedisc"  / timedisc, &
            "datafile"  / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Sim,Mesh,Physics,Timedisc)
    USE geometry_generic_mod
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite),   INTENT(INOUT) :: Sim
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp
    CLASS(geometry_base), ALLOCATABLE :: geo_cyl
    TYPE(Dict_TYP), POINTER :: geo_config
    INTEGER           :: i,j,k
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                      :: bccyl
    REAL              :: r,z,vr,vphi,sden,mu,height,mass
    CHARACTER(LEN=64) :: value
    !------------------------------------------------------------------------!
    geo_config => Dict("geometry" / CYLINDRICAL)
    CALL new_geometry(geo_cyl,geo_config)
    CALL DeleteDict(geo_config)

    ! compute cylindrical coordinates for each cell bary center
    CALL geo_cyl%Convert2Curvilinear(Mesh%bccart,bccyl)

    ! custom boundary conditions if requested
    SELECT TYPE(bwest => Timedisc%Boundary%boundary(WEST)%p)
    CLASS IS (boundary_custom)
      CALL bwest%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER,CUSTOM_NOGRAD/))
    END SELECT
    SELECT TYPE(beast => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS (boundary_custom)
      CALL beast%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER,CUSTOM_NOGRAD/))
    END SELECT

    ! get gravitational acceleration
    sp => Sim%Sources
    DO
      IF (ASSOCIATED(sp).EQV..FALSE.) RETURN
      SELECT TYPE(sp)
      CLASS IS(sources_gravity)
        gp => sp
        EXIT
      END SELECT
      sp => sp%next
    END DO
    IF (.NOT.ASSOCIATED(sp)) CALL Physics%Error("pringle::InitData","no gravity term initialized")

    SELECT CASE(VISTYPE)
    CASE(BETA)
      mu=1.0
    CASE(POWERLAW)
      mu=PL_EXP
    CASE DEFAULT
      CALL Timedisc%Error("pringle::InitData","only pringle and beta viscosity possible")
    END SELECT

    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_eulerisotherm)
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IMIN,Mesh%IGMAX
            ! distance to axis (r-coordinate in cylindrical coordinates)
            r = bccyl(i,j,k,1)
            ! z-coordinate
            z = bccyl(i,j,k,3)
            ! Keplerian velocity
            vphi = r/(r**2+z**2)**0.75 
            ! compute pressure scale height
            height = CSISO * SQRT(r**3) ! h = cs / Omega = cs / SQRT(GM/r**3), with GM=1
#ifdef HAVE_FGSL
            ! use exact solution at time TINIT
            CALL ViscousRing_analytic(mu,TINIT,r,sden,vr)
#else
            ! use approximation for TINIT << 1
            ! no need for modified Bessel functions
            CALL ViscousRing_approx(mu,TINIT,r,sden,vr)
#endif
            pvar%density%data3d(i,j,k) = RHOMIN + sden / (SQRT(2*PI) * height) &
              * EXP(-0.5*(z/height)**2)
!             pvar%velocity%data4d(i,j,k,1) = vr
!             pvar%velocity%data4d(i,j,k,2) = vphi
!             pvar%velocity%data4d(i,j,k,3) = 0.0
          END DO
        END DO
      END DO
      ! determine disk mass
      mass = SUM(Mesh%volume%data1d(:)*pvar%density%data1d(:),MASK=Mesh%without_ghost_zones%mask1d)
#ifdef PARALLEL
      CALL MPI_AllReduce(MPI_IN_PLACE,mass,1,DEFAULT_MPI_REAL,MPI_SUM, &
                         Mesh%comm_cart,ierror)
#endif
      ! rescale disk mass
      pvar%density%data1d(:) = pvar%density%data1d(:) * MDISK / mass
      pvar%velocity%data1d(:) = 0.0 ! initialize all velocities with zero
      CALL gp%UpdateGravity(Mesh,Sim%Physics,Sim%Fluxes,pvar,0.0,0.0)
      pvar%velocity%data4d(:,:,:,1:Physics%VDIM) = &
        Timedisc%GetCentrifugalVelocity(Mesh,Sim%Physics,Sim%Fluxes,Sim%Sources,(/0.,0.,1./), &
          gp%accel%data4d)
    CLASS DEFAULT
      CALL Timedisc%Error("pringle::InitData","only isothermal physics possible")
    END SELECT

    IF (ASSOCIATED(Timedisc%solution)) THEN
      ! store initial condition in exact solution array
      Timedisc%solution = Timedisc%pvar
    END IF

    ! transform to conservative variables
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    CALL Mesh%Info(" DATA-----> initial condition: " // "pringle disk")
    WRITE(value,"(ES9.3)") TVIS
    CALL Mesh%Info("                               " // "viscous timescale:  " //TRIM(value))
    WRITE(value,"(ES9.3)") RE
    CALL Mesh%Info("                               " // "Reynolds number:    " //TRIM(value))
    WRITE(value,"(ES9.3)") MA
    CALL Mesh%Info("                               " // "Mach number:        " //TRIM(value))
    WRITE(value,"(ES9.3)") MDISK
    CALL Mesh%Info("                               " // "Disk mass:          " //TRIM(value))

    CALL geo_cyl%Finalize()
    DEALLOCATE(geo_cyl)
  END SUBROUTINE InitData


  SUBROUTINE ViscousRing_approx(mu,tau,r,Sigma,vr)
    IMPLICIT NONE
    !--------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mu,tau,r
    REAL, INTENT(OUT) :: Sigma,vr
    !--------------------------------------------------------------------!
    REAL :: lambda,x,r4l,invr4l
    !--------------------------------------------------------------------!
    lambda = 1./(4.0-mu)
    r4l    = r**(0.25/lambda)
    x      = 2*lambda**2 / tau * r**(0.25/lambda)
    invr4l = 1./r4l
    Sigma  = MDISK/SQRT((4*PI)**3 * tau) * r4l**(1.5-9*lambda) &
               * EXP(-lambda**2/tau*(1.0-r4l)**2)
    vr     = -1.5/RE * ((x*tau)/(2*lambda*lambda))**(4*lambda-2.) &
               * (1.0 - 0.5*x/lambda * (x*tau/(2*lambda*lambda) - 1.0))
  END SUBROUTINE ViscousRing_approx

#ifdef HAVE_FGSL
  SUBROUTINE ViscousRing_analytic(mu,tau,r,Sigma,vr)
    IMPLICIT NONE
    !--------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mu,tau,r
    REAL, INTENT(OUT) :: Sigma,vr
    !--------------------------------------------------------------------!
    REAL :: lambda,r4l,x
    !--------------------------------------------------------------------!
    lambda = 1./(4.0-mu)
    r4l    = r**(0.25/lambda)
    x      = 2*lambda**2 / tau * r4l
    Sigma  = MDISK/(4*PI) * (lambda/tau) * r4l**2 * r**(-2.25) &
               * EXP(-lambda**2/tau*(1.0+r4l**2)) &
               * modified_bessel_Inu(lambda,x)
    vr     = -1.5/RE * ((x*tau)/(2*lambda*lambda))**(4*lambda-2.) &
               * (1.0 - 0.5*x/lambda * (x*tau/(2*lambda*lambda) &
               - modified_bessel_Inu(1+lambda,x)/modified_bessel_Inu(lambda,x)))
  END SUBROUTINE ViscousRing_analytic

  FUNCTION modified_bessel_Inu(nu,x) RESULT(Inu)
    USE fgsl
    IMPLICIT NONE
    !--------------------------------------------------------------------!
    REAL, INTENT(IN) :: nu,x
    REAL             :: Inu
    !--------------------------------------------------------------------!
    INTEGER          :: n
    !--------------------------------------------------------------------!
    IF (AINT(nu).EQ.nu) THEN
      ! nu is an integer -> use modified bessel functions of integer order
      n = INT(nu)
      Inu = fgsl_sf_bessel_icn(n,x)
    ELSE
      ! nu is not an integer -> use representation of modified bessel functions
      ! through hypergeometric functions 0F1
      Inu = (0.5*x)**nu / fgsl_sf_gamma(1.0+nu) * fgsl_sf_hyperg_0f1(1.0+nu,(0.5*x)**2)
    END IF
  END FUNCTION modified_bessel_Inu
#endif

END PROGRAM pringle3d
