!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: selfsimdisk1d.f90                                                 #
!#                                                                           #
!# Copyright (C) 2014-2015                                                   #
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
!> 1D self-similar accretion disk disk
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
  !--------------------------------------------------------------------------!
  ! the equations are solved in non-dimensional units;
  ! we use geometrical units (c=1,G=1)
  ! simulation parameters
  REAL, PARAMETER :: TSIM        = 2.0E+4        ! simulation time [dyn. time]
  REAL, PARAMETER :: BETA_VIS    = 1.0E-3        ! viscous coupling constant
  ! initial condition parameters
  REAL, PARAMETER :: MBH0        = 1.0           ! (non-dimensional) black hole mass
  REAL, PARAMETER :: KAPPA       = -1.0          ! power law index (-1.5<KAPPA<-0.75)
  REAL, PARAMETER :: CSISO       = 5E-2          ! isothermal speed of sound (<<1)
  REAL, PARAMETER :: Q           = 0.0           ! initial Toomre parameter (<1)
                                                 ! (if Q=0: globally constant speed of sound)
  ! mesh settings
  INTEGER, PARAMETER :: MGEO     = LOGPOLAR      ! geometry
!!$  INTEGER, PARAMETER :: MGEO     = TANPOLAR
!   INTEGER, PARAMETER :: MGEO     = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 40          ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution (1D!)
  REAL, PARAMETER    :: RMIN = 1e-1        ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 1e+3        ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 1e+0        ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'ssdisk1d'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Timedisc, Sim%Mesh, Sim%Physics, Sim%Fluxes)
  
  ! main loop
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4),sgbc
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, sources, &
                               grav, timedisc, fluxes, vis
    REAL              :: x1,x2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
    CASE(SINHPOLAR)
       x1 = LOG(RMIN/GPAR+SQRT(1.0+(RMIN/GPAR)**2))  ! = ASINH(RIN))
       x2 = LOG(RMAX/GPAR+SQRT(1.0+(RMAX/GPAR)**2))
    END SELECT
    

    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / 0.0, &
           "ymax"     / (2.0*PI), &
           "gparam"   / GPAR)

    ! physics settings
    physics => Dict("problem" / EULER2D_ISOTHERM, &
              "units"   / GEOMETRICAL, &        ! Newton constant = 1 !
              "output/bccsound" / 1)

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &   ! vars. to use for reconstruction!
             "limiter"   / OSPRE, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    boundary => Dict( &
!            "western"  / ABSORBING, &
!            "western"  / NO_GRADIENTS, &
           "western"  / CUSTOM, &          ! see below
!            "eastern"  / REFLECTING, &
!            "eastern"  / ABSORBING, &
           "eastern"  / CUSTOM, &          ! see below
!            "eastern"  / FARFIELD, &
!            "eastern"  / NO_GRADIENTS, &
           "southern" / PERIODIC, &
           "northern" / PERIODIC)

    ! source term due to a point mass
    grav => Dict("stype" / GRAVITY, &
            "output/height" / 1, &
            "pmass/gtype"     / POINTMASS, &                ! grav. accel. of a point mass    !
            "pmass/potential" / NEWTON, &                   ! type of gravitational potential !
            "pmass/mass"      / MBH0, & 
            "pmass/outbound" / 1, &           ! enables accretion onto central point mass
            "selfg/gtype" / MONOPOL, &        ! enables monopol approximation for self-gravity
            "selfg/output/enclosed_mass" / 1) ! output the enclosed mass; accounts only for the
                                              ! material in the disk; add the current point mass
                                              !  to get the total enclosed mass

    ! viscosity source term
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / BETA, &
          "dynconst"  / BETA_VIS, &
!           "bulkconst" / 0.0, &
          "cvis"      / 0.05)

    sources => Dict("grav"       / grav, &
              "vis"         / vis)

    ! time discretization settings
    timedisc => Dict( &
!            "method"   / MODIFIED_EULER, &
           "method"   / SSPRK, &
           "order"    / 5, &
           "cfl"      / 0.3, &
           "stoptime" / TSIM, &
           "dtlimit"  / (1.0E-16*TSIM), &
           "tol_rel" / 1.0E-2, &
           "tol_abs" / (/ 0.0, 1.0E-3, 1.0E-3 /), &
           "maxiter"  / 2000000000)

    ! initialize data input/output
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
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
    TYPE(Sources_TYP),POINTER :: sp
    CHARACTER(LEN=32) :: info_str
    INTEGER           :: i,j,k
    REAL              :: Sigma0,mass
    REAL              :: y0,x,y,z
    REAL              :: er(2),ephi(2),vphi,vr
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4) :: fcsound
    REAL, DIMENSION(:,:), POINTER :: r
    REAL, DIMENSION(:,:,:), POINTER :: r_vec
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Fluxes
    INTENT(INOUT)     :: Timedisc,Physics
    !------------------------------------------------------------------------!
    x(z) = -1.5 + (2./3.*KAPPA+1.0)*z/(1.+2./3.*z)
    y(z) = z/(1.+2./3.*z)**(2./3.*KAPPA+1.)            ! z = y0*xi**3/2

    ! check parameter
    IF (KAPPA.LT.-1.5.OR.KAPPA.GT.-0.5) &
       CALL Error(Timedisc,"InitData","parameter kappa out of range")

    ! set radius pointer
    r => RemapBounds(Mesh,Mesh%radius%bcenter)
    r_vec => RemapBounds(Mesh,Mesh%posvec%bcenter)

    ! set initial conditions using an approximate solution of the 
    ! self-similar disk evolution equation (roughly a broken power law)
    y0 = 1./SQRT(KAPPA**2*MBH0)
    Sigma0 = 0.8/(2*PI*KAPPA**2)
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=Mesh%IMIN,Mesh%IMAX
          ! dimensionless surface density
          Timedisc%pvar(i,j,Physics%DENSITY) = Sigma0*(2*x(y0*r(i,j)**1.5)+3.) &
                     *r(i,j)/(y(y0*r(i,j)**1.5))**2
       END DO
    END DO

    ! get gravitational acceleration
    sp => GetSourcesPointer(Physics%sources,GRAVITY)       
    IF (ASSOCIATED(sp)) THEN
       CALL UpdateGravity(sp,Mesh,Physics,Fluxes,0.0,0.0,Timedisc%pvar)
    ELSE
       CALL Error(Sim,"InitData","no gravity term initialized")
    END IF

    ! set isothermal speed of sound
    bccsound(:,:) = CSISO
    fcsound(:,:,:) = CSISO

    ! balance of gravitational and centrifugal forces
    DO j=Mesh%JMIN,Mesh%JMAX
       DO i=Mesh%IMIN,Mesh%IMAX
          ! compute centrifugal velocity
          vphi = SQRT(MAX(0.0,-r_vec(i,j,1)*sp%accel(i,j,1)-r_vec(i,j,2)*sp%accel(i,j,2)))
          ! radial velocity
!           vr = -2*BETA_VIS/KAPPA*r(i,j) &
!                      *(x(y0*r(i,j)**1.5)-KAPPA)/(2*x(y0*r(i,j)**1.5)+3.)
!           vr = -2* BETA_VIS * vphi
!                      *(x(y0*r(i,j)**1.5)-KAPPA)/(x(y0*r(i,j)**1.5)+1.5)
          vr = -1.56*BETA_VIS
          ! curvilinear components of azimuthal unit vector
          ephi(1) = -r_vec(i,j,2)/r(i,j)
          ephi(2) = r_vec(i,j,1)/r(i,j)
          ! curvilinear components of radial unit vector
          er(1) = r_vec(i,j,1)/r(i,j)
          er(2) = r_vec(i,j,2)/r(i,j)
          ! set curvilinear velocity field
          Timedisc%pvar(i,j,Physics%XVELOCITY:Physics%YVELOCITY) = vr * er(1:2) + vphi * ephi(1:2)
          IF (Q.GT.TINY(Q)) THEN
             ! compute local sound speed assuming a constant Toomre parameter Q
             bccsound(i,j) = Q*PI*Physics%Constants%GN*Timedisc%pvar(i,j,Physics%DENSITY) &
               * r(i,j) / (vphi * SQRT(2*x(y0*r(i,j)**1.5)+4.) )
             DO k=1,4
                fcsound(i,j,k) = Q*PI*Physics%Constants%GN*Timedisc%pvar(i,j,Physics%DENSITY) &
                   * Mesh%radius%faces(i,j,k) / (vphi * SQRT(2*x(y0*Mesh%radius%faces(i,j,k)**1.5)+4.) )
             END DO
          END IF
       END DO
    END DO

    CALL SetSoundSpeeds(Physics,Mesh,bccsound)
    CALL SetSoundSpeeds(Physics,Mesh,fcsound)

    ! custom boundary conditions at western boundary if requested
    IF (GetType(Timedisc%boundary(WEST)).EQ.CUSTOM) THEN
       Timedisc%boundary(WEST)%cbtype(:,Physics%DENSITY) = CUSTOM_KEPLER
       Timedisc%boundary(WEST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_EXTRAPOL
       Timedisc%boundary(WEST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
       IF (Physics%PRESSURE.NE.0) &
          Timedisc%boundary(WEST)%cbtype(:,Physics%PRESSURE) = CUSTOM_NOGRAD
    END IF

    ! custom boundary conditions at eastern boundary if requested
    IF (GetType(Timedisc%boundary(EAST)).EQ.CUSTOM) THEN
       Timedisc%boundary(EAST)%cbtype(:,Physics%DENSITY) = CUSTOM_KEPLER
       Timedisc%boundary(EAST)%cbtype(:,Physics%XVELOCITY) = CUSTOM_NOGRAD
       Timedisc%boundary(EAST)%cbtype(:,Physics%YVELOCITY) = CUSTOM_FIXED
       Timedisc%boundary(EAST)%invRscale(:,:) = Timedisc%boundary(EAST)%invRscale(:,:)**2
       IF (Physics%PRESSURE.NE.0) &
          Timedisc%boundary(EAST)%cbtype(:,Physics%PRESSURE) = CUSTOM_NOGRAD
       ! initial data for fixed boundary conditions
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             ! linear extrapolation of initial data
             Timedisc%boundary(EAST)%data(i,j,:) = (i+1)*Timedisc%pvar(Mesh%IMAX,j,:) &
               - i*Timedisc%pvar(Mesh%IMAX-1,j,:)
          END DO
       END DO
    END IF

    ! farfield boundary conditions at eastern boundary if requested
    IF (GetType(Timedisc%boundary(EAST)).EQ.FARFIELD) THEN
       ! set farfield data using initial conditions
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             ! linear extrapolation of initial data
             Timedisc%boundary(EAST)%data(i,j,:) = (i+1)*Timedisc%pvar(Mesh%IMAX,j,:) &
               - i*Timedisc%pvar(Mesh%IMAX-1,j,:)
          END DO
       END DO
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    CALL Info(Mesh," DATA-----> initial condition: " // "self-similar solution")
    WRITE (info_str, '(F5.2)') KAPPA
    CALL Info(Mesh,"            kappa:             " // TRIM(info_str))
    WRITE (info_str, '(ES8.2)') MBH0
    CALL Info(Mesh,"            initial BH mass:   " // TRIM(info_str))
    ! compute disk mass (within computational domain)
    mass = SUM(Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
         * Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
    WRITE (info_str, '(ES8.2)') mass
    CALL Info(Mesh,"            initial disk mass: " // TRIM(info_str))

  END SUBROUTINE InitData
 END PROGRAM Init
