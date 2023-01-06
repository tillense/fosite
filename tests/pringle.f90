!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: pringle.f90                                                       #
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
!> \test Jim Pringles famous viscous spreading ring solution
!! \author Björn Sperling
!! \author Tobias Illenseer
!!
!! The equations are solved in non-dimensional units, using the Keplerian
!! velocity at the location of the initial ring at R=1 as the velocity scale.
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
!! \warning compile with autodouble
!----------------------------------------------------------------------------!
PROGRAM pringle_test
  USE fosite_mod
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 2.0E-2     ! simulation time [TVIS] see below
                                             ! should be < 1 ~ viscous time scale
  ! these are the basic parameters; for stable solutions one requires
  ! MA >> 1 and RE > MA
  REAL, PARAMETER    :: RE      = 1.0E+4     ! Reynolds number (at R=1)
  REAL, PARAMETER    :: MA      = 1.0E+2     ! Mach number (at R=1)
  REAL, PARAMETER    :: MDISK   = 1.0E-2     ! disk mass << 1 = central mass
  ! lower limit for nitial density
  REAL, PARAMETER    :: RHOMIN  = 1.0E-20    ! minimal initial density
  ! viscosity prescription
!   INTEGER, PARAMETER :: VISTYPE = BETA
  INTEGER, PARAMETER :: VISTYPE = POWERLAW
  REAL, PARAMETER    :: PL_EXP  = 0.0        ! exponent for power law viscosity
  REAL, PARAMETER    :: TVIS    = 4./3.*RE   ! viscous time scale
  REAL, PARAMETER    :: TINIT   = 1.0E-3     ! time for initial condition [TVIS]
  ! mesh settings
  !**************************************************************************!
  !> \remark #1: This test is very sensitive to mesh settings. If there are
  !!         not enough mesh points at small radii instabilities grow
  !!         starting at the inner boundary.
  !! \remark #2: The new implementation of the numerical stability criterion
  !!         for advection-diffusion problems seems to work far better than
  !!         before. However, it is sometimes necessary to adjust the viscous
  !!         Courant number especially in very high Mach number flows MA~1000.
  !**************************************************************************!
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!   INTEGER, PARAMETER :: MGEO = LOGCYLINDRICAL
!   INTEGER, PARAMETER :: MGEO = SPHERICAL !!! ATTENTION: not applicable in 1D
  INTEGER, PARAMETER :: XRES = 400        ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution
  INTEGER, PARAMETER :: ZRES = 1           ! z-resolution
  REAL, PARAMETER    :: RMIN = 0.05        ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 2.0         ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'pringle'
  ! derived parameters
  REAL, PARAMETER    :: CSISO = 1./MA      ! isothermal speed of sound
  REAL               :: X0 = 0.0           ! cartesian position of point
  REAL               :: Y0 = 0.0           !   mass (default: origin)
  REAL               :: Z0 = 0.0           !
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  REAL, DIMENSION(:), ALLOCATABLE :: sigma
  REAL :: sum_numer, sum_denom
  INTEGER :: n,DEN,VX,VY
  REAL    :: next_output_time
  INTEGER :: i,j,k
  LOGICAL :: ok
  !--------------------------------------------------------------------------!

  ALLOCATE(Sim)
  CALL Sim%InitFosite()

#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
#ifdef HAVE_FGSL
TAP_PLAN(4)
#else
TAP_PLAN(1)
#endif
#ifdef PARALLEL
  END IF
#endif

  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim,Sim%Mesh, Sim%Physics, Sim%Timedisc)
  next_output_time = Sim%Datafile%time
  IF (ASSOCIATED(Sim%Timedisc%solution)) THEN
#ifdef HAVE_FGSL
    DO WHILE((Sim%Timedisc%maxiter.LE.0).OR.(Sim%iter.LE.Sim%Timedisc%maxiter))
      ! advance numerical solution
      IF(Sim%Step()) EXIT
      IF (next_output_time.LT.Sim%Datafile%time) THEN
        ! compute exact solution for next output time step
        next_output_time = Sim%Datafile%time
        SELECT TYPE(pvar => Sim%Timedisc%solution)
        TYPE IS (statevector_eulerisotherm)
          k=Sim%Mesh%KMIN
          DO j=Sim%Mesh%JMIN,Sim%Mesh%JMAX
            DO i=Sim%Mesh%IMIN,Sim%Mesh%IMAX
              IF (VISTYPE.EQ.BETA) THEN
                CALL ViscousRing_analytic(1.0,next_output_time/TVIS+TINIT, &
                  Sim%Mesh%radius%bcenter(i,j,k), &
                  pvar%density%data3d(i,j,k),pvar%velocity%data4d(i,j,k,1))
              ELSE
                CALL ViscousRing_analytic(PL_EXP,next_output_time/TVIS+TINIT, &
                  Sim%Mesh%radius%bcenter(i,j,k), &
                  pvar%density%data3d(i,j,k),pvar%velocity%data4d(i,j,k,1))
              END IF
            END DO
          END DO
        END SELECT
      END IF
    END DO
#else
    CALL Sim%Error("pringle","computation of analytical solution requires GSL support")
#endif
  ELSE
    CALL Sim%Run()
  END IF
  ok = .NOT.Sim%aborted
  ALLOCATE(sigma(Sim%Physics%VNUM))
  ! compare with exact solution if requested
  IF (ASSOCIATED(Sim%Timedisc%solution)) THEN
    DO n=1,Sim%Physics%VNUM
      ! use L1 norm to estimate the deviation from the exact solution:
      !   Σ |pvar - pvar_exact| / Σ |pvar_exact|
      sum_numer = SUM(ABS(Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n) &
                        -Sim%Timedisc%solution%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n)))
      sum_denom = SUM(ABS(Sim%Timedisc%solution%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n)))
#ifdef PARALLEL
      IF (Sim%GetRank().GT.0) THEN
        CALL MPI_Reduce(sum_numer,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
        CALL MPI_Reduce(sum_denom,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
      ELSE
        CALL MPI_Reduce(MPI_IN_PLACE,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
        CALL MPI_Reduce(MPI_IN_PLACE,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
#endif
      sigma(n) = sum_numer / sum_denom
#ifdef PARALLEL
      END IF
#endif
    END DO
  ELSE
    sigma(:) = 0.0
  END IF

  DEN = Sim%Physics%DENSITY
  VX  = Sim%Physics%XVELOCITY
  VY  = Sim%Physics%YVELOCITY

#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
TAP_CHECK(ok,"stoptime reached")
! These lines are very long if expanded. So we can't indent it or it will be cropped.
#ifdef HAVE_FGSL
TAP_CHECK_SMALL(sigma(DEN),5.0E-02,"density deviation < 5%")
TAP_CHECK_SMALL(sigma(VX),8.0E-02,"radial velocity deviation < 8%")
! skip azimuthal velocity deviation, because exact value is 0
TAP_CHECK_SMALL(sigma(VY),2.0E-04,"azimuthal velocity deviation < 0.02%")
#endif
TAP_DONE
#ifdef PARALLEL
  END IF
#endif

  IF (ALLOCATED(sigma)) DEALLOCATE(sigma)
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
      z1 = 0.0
      z2 = 0.0
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(WEST)  = ABSORBING
!       bc(EAST)  = ABSORBING
      bc(WEST)  = CUSTOM
      bc(EAST)  = CUSTOM
      bc(SOUTH) = PERIODIC
      bc(NORTH) = PERIODIC
      bc(BOTTOM) = NO_GRADIENTS
      bc(TOP)    = NO_GRADIENTS
      cvis = 0.3
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
      cvis = 0.3
    CASE(SPHERICAL)
      x1 = RMIN
      x2 = RMAX
      y1 = PI/2 ! - 0.05*PI
      y2 = PI/2 ! + 0.05*PI
      z1 = -PI
      z2 = PI
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
      bc(WEST)  = CUSTOM
      bc(EAST)  = CUSTOM
      bc(SOUTH) = NO_GRADIENTS
      bc(NORTH) = NO_GRADIENTS
      bc(BOTTOM) = PERIODIC
      bc(TOP)    = PERIODIC
      cvis = 0.3
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
            "pmass/x"   / X0, &            ! cartesian position of point mass
            "pmass/y"   / Y0, &
            "pmass/outbound" / 0)          ! disable accretion

    ! combine all source terms
    sources => Dict( &
            "grav"      / grav, &
            "viscosity" / vis)

    ! time discretization settings
    timedisc => Dict( &
            "method"    / SSPRK, &
            "stoptime"  / (TSIM * TVIS), &
            "cfl"       / 0.4, &
            "dtlimit"   / (1E-10 * TVIS), &
!             "rhstype"   / 1, &
            "maxiter"   / 100000000, &
            "tol_rel"   / 0.01, &
            "tol_abs"   / (/0.0,1e-5,0.0/))

    ! enable output of analytical solution
#ifdef HAVE_FGSL
    IF (MGEO.EQ.CYLINDRICAL) &
      CALL SetAttr(timedisc,"output/solution",1)
#endif

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
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite),   INTENT(INOUT) :: Sim
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp => null()
    INTEGER           :: i,j,k
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                      :: posvec,ephi
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: radius
    REAL              :: r,vr,vphi,sden,mu
    CHARACTER(LEN=64) :: value
    !------------------------------------------------------------------------!
    IF (ABS(X0).LE.TINY(X0).AND.ABS(Y0).LE.TINY(Y0)) THEN
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:,:) = Mesh%radius%bcenter(:,:,:)
       posvec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,:,1) = X0
       posvec(:,:,:,2) = Y0
       posvec(:,:,:,3) = Z0
       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - posvec(:,:,:,:)
       ! compute its absolute value
       radius(:,:,:) = SQRT(posvec(:,:,:,1)**2+posvec(:,:,:,2)**2+posvec(:,:,:,3)**2)
    END IF

    ! curvilinear components of azimuthal unit vector
    ! (maybe with respect to shifted origin)
    ! from ephi = ez x er = ez x posvec/radius = ez x (rxi*exi + reta*eeta)/r
    !             = rxi/r*(ez x exi) + reta/r*(ez x eeta) = rxi/r*eeta - reta/r*exi
    ! because (ez,exi,eeta) is right handed orthonormal set of basis vectors
    ephi(:,:,:,1) = -posvec(:,:,:,2)/radius(:,:,:)
    ephi(:,:,:,2) = posvec(:,:,:,1)/radius(:,:,:)
    ephi(:,:,:,3) = 0.0

    ! custom boundary conditions if requested
    SELECT TYPE(bwest => Timedisc%Boundary%boundary(WEST)%p)
    CLASS IS (boundary_custom)
      CALL bwest%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER/))
    END SELECT
    SELECT TYPE(beast => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS (boundary_custom)
      CALL beast%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_LOGEXPOL,CUSTOM_OUTFLOW,CUSTOM_KEPLER/))
    END SELECT

    ! get gravitational acceleration
    sp => Sim%Sources
    DO
      IF (.NOT.ASSOCIATED(sp)) EXIT 
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
            ! distance to center of mass
            r = radius(i,j,k)
            ! Keplerian velocity
            vphi = SQRT(1.0/r)
#ifdef HAVE_FGSL
            ! use exact solution at time TINIT
            CALL ViscousRing_analytic(mu,TINIT,r,sden,vr)
#else
            ! use approximation for TINIT << 1
            ! no need for modified Bessel functions
            CALL ViscousRing_approx(mu,TINIT,r,sden,vr)
#endif
            pvar%density%data3d(i,j,k) = RHOMIN + sden
            ! curvilinear velocity components
            pvar%velocity%data4d(i,j,k,1:2) = &
                vr*posvec(i,j,k,1:2)/r + vphi*ephi(i,j,k,1:2)
!               pvar%velocity%data4d(i,j,k,1:2) = vr*posvec(i,j,k,1:2)/r
!               p%velocity%data4d(i,j,k,2) = vphi
          END DO
        END DO
      END DO
!         CALL gp%UpdateGravity(Mesh,Sim%Physics,Sim%Fluxes,pvar,0.0,0.0)
!         pvar%velocity%data4d(:,:,:,1:Physics%VDIM) = &
!           Timedisc%GetCentrifugalVelocity(Mesh,Sim%Physics,Sim%Fluxes,Sim%Sources,(/0.,0.,1./), &
!             gp%accel%data4d)
    CLASS DEFAULT
      CALL Timedisc%Error("pringle::InitData","only isothermal physics possible")
    END SELECT

    IF (ASSOCIATED(Timedisc%solution)) THEN
      ! store initial condition in exact solution array
      Timedisc%solution = Timedisc%pvar
    END IF

    ! check fargo and set background velocity field
    IF (Mesh%FARGO.EQ.2) &
       Timedisc%w(:,:) = SQRT(1.0/radius(:,Mesh%JMIN,:))

    ! transform to conservative variables
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    CALL Mesh%Info(" DATA-----> initial condition: " // "viscous spreading ring")
    WRITE(value,"(ES9.3)") TVIS
    CALL Mesh%Info("                               " // "viscous timescale:  " //TRIM(value))
    WRITE(value,"(ES9.3)") RE
    CALL Mesh%Info("                               " // "Reynolds number:    " //TRIM(value))
    WRITE(value,"(ES9.3)") MA
    CALL Mesh%Info("                               " // "Mach number:        " //TRIM(value))

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

END PROGRAM pringle_test
