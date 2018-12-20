!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: pringle.f90                                                       #
!#                                                                           #
!# Copyright (C) 2008-2012                                                   #
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
!> a thin viscous ring rotating around a central pointmass
!! \author Björn Sperling
!! \author Tobias Illenseer
!!
!! References:
!!  [1] Lynden-Bell, D. & Pringle, J. E. Evolution of Viscous Disks and Origin
!!      of Nebular Variables, M.N.R.A.S vol. 168(3) (1974) pp. 603-637
!!      http://adsabs.harvard.edu/abs/1974MNRAS.168..603L
!!  [2] Pringle, J. E. Accretion discs in astrophysics
!!      Annual review of astronomy and astrophysics. vol. 19 (1981) pp. 137-162
!!      DOI: 10.1146/annurev.aa.19.090181.001033
!!  [3] Speith, R. and Kley, W. Stability of the viscously spreading ring
!!      Astron. Astrophys. vol. 399 (2003), pp. 395-407
!!      DOI: 10.1051/0004-6361:20021783
!!
!! \warning compile with autodouble
!----------------------------------------------------------------------------!
PROGRAM pringle_test
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE geometry_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER :: GN      = 6.6742D-11  ! Newtons grav. constant [SI]
  ! simulation parameters
  REAL, PARAMETER :: TSIM     = 1.0E-1     ! simulation time [TAU] see below
  ! the equations are solved in non-dimensional units, using the Keplerian
  ! velocity at the location of the initial ring at R=1 as the velocity scale
  REAL, PARAMETER :: CENTMASS = 1./GN      ! fixed for non-dim. equations
  ! these are the basic parameters; for stable solutions one requires
  ! MA >> 1 and RE > MA
  REAL, PARAMETER :: RE       = 1.0E+4     ! Reynolds number (at R=1)
  REAL, PARAMETER :: MA       = 1.0E+2     ! Mach number (at R=1)
  ! lower limit for nitial density
  REAL, PARAMETER :: RHOMIN   = 1.0E-20    ! minimal initial density
  ! viscosity prescription
!!$  INTEGER, PARAMETER :: VISTYPE = BETA     
  INTEGER, PARAMETER :: VISTYPE = PRINGLE
  REAL, PARAMETER :: TAU0     = 0.01       ! time for initial condition [TAU]
  ! mesh settings
  !**************************************************************************!
  ! REMARK #1: This test is very sensitive to mesh settings. If there are
  !         not enough mesh points at small radii instabilities grow
  !         starting at the inner boundary.
  ! REMARK #2: There is an unsolved stability issue with the viscous
  !         source term. The stability condition doesn't seem to be
  !         sufficient. Thus in case of high Reynolds numbers one has
  !         to adjust the viscous CFL number cvis to enforce smaller time
  !         steps and preserve stability.
  !**************************************************************************!
  INTEGER, PARAMETER :: MGEO = POLAR       ! geometry
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!  INTEGER, PARAMETER :: MGEO = SINHPOLAR
!  INTEGER, PARAMETER :: MGEO = BIPOLAR
  INTEGER, PARAMETER :: XRES = 200         ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution
  REAL, PARAMETER    :: RMIN = 0.01        ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 2.0         ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 0.8         ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'pringle'
  ! derived parameters
  REAL, PARAMETER    :: CSISO = 1./MA      ! isothermal speed of sound
  REAL               :: TAU                ! viscous time scale
  REAL               :: OMEGA0             ! angular velocity at R0
  REAL               :: X0 = 0.0           ! cartesian position of point
  REAL               :: Y0 = 0.0           !   mass (default: origin)
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim, Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Simulation finished")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               grav, vis, timedisc, fluxes, sources
    REAL              :: x1,x2,y1,y2
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = -PI
       y2 = PI
       bc(WEST)  = ABSORBING
       bc(EAST)  = ABSORBING
!        bc(WEST)  = NO_GRADIENTS
!        bc(EAST)  = NO_GRADIENTS
!        bc(WEST)  = CUSTOM
!        bc(EAST)  = CUSTOM
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = -PI 
       y2 = PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = -PI 
       y2 = PI       
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       x1 = Asinh(RMIN/GPAR)
       x2 = Asinh(RMAX/GPAR)
       y1 = -PI
       y2 = PI       
!        bc(WEST)  = ABSORBING
!        bc(EAST)  = ABSORBING
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
!        bc(WEST)  = CUSTOM
!        bc(EAST)  = CUSTOM
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(BIPOLAR)
       x1 = 0.0
       x2 = 2*PI
       y2 = Asinh(GPAR/RMIN)
       y1 = -y2
       bc(WEST)  = PERIODIC
       bc(EAST)  = PERIODIC
!        bc(SOUTH) = NO_GRADIENTS
!        bc(NORTH) = FOLDED
!        bc(SOUTH) = ABSORBING
!        bc(NORTH) = ABSORBING
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
!        bc(SOUTH) = REFLECTING
!        bc(NORTH) = REFLECTING
       ! shift the point mass into the left hole at tau -> -inf
       X0 = GPAR*(1.0+2./(EXP(2.*y1)-1.0))
    CASE DEFAULT
       CALL Error(Sim,"InitProgram","mesh geometry not supported for Pringle disk")
    END SELECT

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "gparam"   / GPAR)

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    physics => Dict("problem" / EULER2D_ISOTHERM, &
              "gamma"   / 1.4, &                   ! ratio of specific heats        !
              "cs"      / CSISO)                   ! isothermal speed of sound      !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! compute viscous time scale
    SELECT CASE(VISTYPE)
    CASE(BETA)
       TAU = 4./27. * RE
    CASE(PRINGLE)
       TAU = RE / 12.0
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","viscosity type not supported")
    END SELECT

    ! source term due to viscosity
    vis => Dict("stype"    / VISCOSITY, &
          "vismodel" / VISTYPE, &
          "dynconst" / (1./RE), &
          "cvis"     / 0.003)

    ! source term due to a point mass
    grav => Dict("stype"  / GRAVITY, &
            "pmass/gtype" / POINTMASS, &    ! grav. accel. of a point mass     !
            "pmass/mass" / CENTMASS, &      ! mass of the accreting object[kg] !
            "pmass/x" / X0, &               ! cartesian position of point mass !
            "pmass/y" / Y0, &
            "pmass/outbound" / 0)           ! disable accretion

    ! combine all source terms
    sources => Dict("grav" / grav, &
                    "viscosity" / vis)

    ! time discretization settings
    timedisc => Dict( &
           "method"   / SSPRK, &
           "stoptime" / (TSIM * TAU), &
           "cfl"      / 0.4, &
           "dtlimit"  / (1E-10 * TAU), &
           "maxiter"  / 100000000, &
           "tol_rel"  / 0.01, &
           "tol_abs"  / (/0.0,1e-5,0.0/))

    ! initialize data input/output
    datafile => Dict( &
!           "fileformat" / NETCDF, &
!           "fileformat" / BINARY, &
           "fileformat" / GNUPLOT, "filecycles" / 0, &
           "decimals" / 6, "cartcoords" / 1, &
           "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
           "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Sim, Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: posvec,ephi
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: radius
    REAL              :: r,r34,vr,vphi
    CHARACTER(LEN=64) :: value
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Sim, Timedisc
    !------------------------------------------------------------------------!
    IF (ABS(X0).LE.TINY(X0).AND.ABS(Y0).LE.TINY(Y0)) THEN
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:) = Mesh%radius%bcenter(:,:)
       posvec(:,:,:) = Mesh%posvec%bcenter(:,:,:)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,1) = X0
       posvec(:,:,2) = Y0
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing 
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:) = Mesh%posvec%bcenter(:,:,:) - posvec(:,:,:)
       ! compute its absolute value
       radius(:,:) = SQRT(posvec(:,:,1)**2+posvec(:,:,2)**2)
    END IF

    ! curvilinear components of azimuthal unit vector
    ! (maybe with respect to shifted origin)
    ! from ephi = ez x er = ez x posvec/radius = ez x (rxi*exi + reta*eeta)/r
    !             = rxi/r*(ez x exi) + reta/r*(ez x eeta) = rxi/r*eeta - reta/r*exi
    ! because (ez,exi,eeta) is right handed orthonormal set of basis vectors
    ephi(:,:,1) = -posvec(:,:,2)/radius(:,:)
    ephi(:,:,2) = posvec(:,:,1)/radius(:,:)

    SELECT CASE(VISTYPE)
    CASE(BETA)
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! distance to center of mass
             r = radius(i,j)
             ! radius parameter r**(3/4)
             r34 = r**0.75
             ! Keplerian velocity
             vphi = SQRT(GN*CENTMASS/r)
             ! radial velocity (approximation for small TAU0)
             vr = 4.5/(RE*TAU0) * r34*(r34-1.0) * vphi
             ! surface density
             Timedisc%pvar(i,j,Physics%DENSITY) = RHOMIN &
                  + 3. / SQRT(TAU0*(4*PI*r34)**3) * EXP(-((1.0-r34)**2)/TAU0)
             ! curvilinear velocity components
             Timedisc%pvar(i,j,Physics%XVELOCITY:Physics%YVELOCITY) = &
                  vr*posvec(i,j,:)/r + vphi*ephi(i,j,:)
         END DO
      END DO
    CASE(PRINGLE)
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! distance to center of mass
             r = radius(i,j)
             ! Keplerian velocity
             vphi = SQRT(GN*CENTMASS/r)
             ! radial velocity
             IF (2*r/TAU0.LE.0.1) THEN
                ! use series expansion, see below
                vr = func_vr(TAU0,r)
             ELSE
                ! approximation for large values of 2*r/t
                vr = 3.0/RE * (0.25/r + 2.0/TAU0*(r-1.0))
             END IF
             ! surface density (approximate solution)
             Timedisc%pvar(i,j,Physics%DENSITY) = RHOMIN &
                  + 1. / SQRT(4*TAU0*(PI*SQRT(r))**3) * EXP(-(1.0-r)**2/TAU0)
             ! curvilinear velocity components
             Timedisc%pvar(i,j,Physics%XVELOCITY:Physics%YVELOCITY) = &
                  vr*posvec(i,j,:)/r + vphi*ephi(i,j,:)
          END DO
       END DO
    END SELECT

    ! set pressure to constant value if necessary
    IF (Physics%PRESSURE.NE.0) &
       Timedisc%pvar(:,:,Physics%PRESSURE) = 1.0E-15

    ! transform to conservative variables
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    ! custom boundary conditions if requested
    IF (GetType(Timedisc%boundary(WEST)).EQ.CUSTOM) THEN
       Timedisc%boundary(WEST)%cbtype(:,Physics%DENSITY) = CUSTOM_NOGRAD
       Timedisc%boundary(WEST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_EXTRAPOL
       Timedisc%boundary(WEST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
       IF (Physics%PRESSURE.NE.0) &
          Timedisc%boundary(WEST)%cbtype(:,Physics%PRESSURE) = CUSTOM_NOGRAD
    END IF
    IF (GetType(Timedisc%boundary(EAST)).EQ.CUSTOM) THEN
       Timedisc%boundary(EAST)%cbtype(:,Physics%DENSITY) = CUSTOM_NOGRAD
       Timedisc%boundary(EAST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_NOGRAD
       Timedisc%boundary(EAST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_NOGRAD
       IF (Physics%PRESSURE.NE.0) &
          Timedisc%boundary(EAST)%cbtype(:,Physics%PRESSURE) = CUSTOM_NOGRAD
    END IF

    CALL Info(Mesh, " DATA-----> initial condition: " // "pringle disk")
    WRITE(value,"(ES9.3)") TAU
    CALL Info(Mesh, "                               " // "viscous timescale:  " //TRIM(value))
    WRITE(value,"(ES9.3)") RE
    CALL Info(Mesh, "                               " // "Reynolds number:    " //TRIM(value))
    WRITE(value,"(ES9.3)") MA
    CALL Info(Mesh, "                               " // "Mach number:        " //TRIM(value))

  END SUBROUTINE InitData


  FUNCTION func_vr(t,r) RESULT(v)
    IMPLICIT NONE
    !--------------------------------------------------------------------!
    REAL, INTENT(IN) :: t,r
    REAL :: v
    !--------------------------------------------------------------------!
    INTEGER :: k
    REAL :: ak,qk,y,y2,y2k,sum1,sum2
    REAL, PARAMETER :: GAMMA_QUARTER = 3.62560990822191
    !--------------------------------------------------------------------!
    ak = 1./GAMMA_QUARTER ! 1./Gamma(0.25) with Gamma-function
    qk = ak
    y  = r/t
    y2 = y*y
    y2k = 1.0
    sum1 = qk
    sum2 = 4.*qk
    DO k=1,100
       ak = ak /(k*(k-0.75))
       y2k = y2k*y2
       qk = ak*y2k
       sum1 = sum1 + qk
       sum2 = sum2 + qk/(k+0.25)
    END DO
    ! sum/(y*sum2) = I_{-3/4}(2*r/t) / I_{1/4}(2*r/t) 
    ! where I_n are modified Bessel functions of the first kind
    v = 6./(RE*t)*(r-sum1/(y*sum2))
  END FUNCTION func_vr

END PROGRAM pringle_test
