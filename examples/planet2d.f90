!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: planet2d.f90                                                      #
!#                                                                           #
!# Copyright (C) 2006-2019                                                   #
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
!> \example planet2d.f90

!----------------------------------------------------------------------------!
PROGRAM planet2d
  USE fosite_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: GN      = 6.6742D-11        ! Newtons grav. constant !
  REAL, PARAMETER    :: CC      = 2.99792458D+8     ! speed of light         !
  REAL, PARAMETER    :: RG      = 8.31447           ! molar gas constant     !
  REAL, PARAMETER    :: AU      = 1.49597870691E+11 ! astronomical unit [m]  !
  REAL, PARAMETER    :: YEAR    = 3.15576E+7        ! Julian year [sec]      !
  REAL, PARAMETER    :: DAY     = 8.6400E+4         ! Day [sec]              !
  REAL, PARAMETER    :: ERAD    = 6.37E+6           ! radius of earth        !
  ! simulation parameter
  REAL, PARAMETER    :: TSIM    = 10*DAY            ! simulation stop time   !
  INTEGER, PARAMETER :: XRES    = 1                 ! theta-resolution       !
  INTEGER, PARAMETER :: YRES    = 32                ! phi-resolution         !
  INTEGER, PARAMETER :: ZRES    = 32                ! phi-resolution         !
  ! geometrical parameter
  REAL, PARAMETER    :: GPAR    = 1.0*ERAD          ! planet-radius          !
  REAL, PARAMETER    :: THETA0  = 0.00              ! axis-plane-angle       !
  REAL, PARAMETER    :: PHI0    = 0.0               !                        !
!  REAL, PARAMETER    :: OMEGA   = -2.99E-7         ! ang. rotation [rad/s]  !
  REAL, PARAMETER    :: OMEGA   = 0.0               ! ang. rotation [rad/s]  !
!  REAL, PARAMETER    :: FSUN    = -6.229E-7        ! freq.(!) day-night[/s] !
  REAL, PARAMETER    :: FSUN    = 0.0               ! freq.(!) day-night[/s] !
  REAL, PARAMETER    :: RD      = 1*AU              ! distance star-planet   !
                                                    ! planet in AU           !
  REAL, PARAMETER    :: PLANET_YEAR = YEAR          ! trop. yr of the planet !
  ! gaseous parameter
  REAL, PARAMETER    :: GAMMA   = 1.4               ! ratio of specific heats!
  REAL, PARAMETER    :: P0      = 1.014E+5          ! surf. press.[N/m**3]   !
  REAL, PARAMETER    :: MU      = 2.897E-2          ! molar mass of atmosp.  !
  REAL, PARAMETER    :: T_0     = 287.76            ! med. init. temp. [K]   !
  REAL, PARAMETER    :: RHO0    = P0*MU/(T_0*RG)    ! surf. dens.[kg/m**3]   !
  ! other atmospherical parameter
  REAL, PARAMETER    :: ALBEDO  = 0.3               ! albedo of the planet   !
  REAL, PARAMETER    :: INTENSITY= 1.4E+3           ! inten. at 1 AU         !
  REAL, PARAMETER    :: HEIGHT  = 2.0E+4            ! height of the atmo.    !
  ! other parameters                                ! opac.: dt/dp=-1/g*kappa!
  REAL, PARAMETER    :: GACC    = 9.81              ! grav. accel. [m/s**2]  !
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 10                ! num. output data sets  !
  CHARACTER(LEN=256), PARAMETER &                   ! output data dir        !
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &                   ! output data file name  !
                     :: OFNAME  = 'planet2d'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)

CALL Sim%InitFosite()
CALL MakeConfig(Sim, Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
CALL Sim%Run()
CALL Sim%Finalize()

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes, heating, sources,&
                               cooling, rotframe
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    !mesh settings
    mesh => Dict( &
         "meshtype"         / MIDPOINT, &
         "geometry"         / SPHERICAL, &
         "omega"            / OMEGA, &
         "inum"             / XRES, &
         "jnum"             / YRES, &
         "knum"             / ZRES, &
         "xmin"             / GPAR, &
         "xmax"             / GPAR, &
         "ymin"             / (0.01), &
         "ymax"             / (PI-0.01), &
         "zmin"             / (-PI), &
         "zmax"             / (PI), &
         "gparam"           / GPAR, &
         "dz"               / HEIGHT, &
         "output/rotation"  / 0, &
         "output/volume"    / 1, &
         "output/dAz"       / 1)

    ! boundary conditions
    boundary => Dict( &
         "western"          / REFLECTING, &
         "eastern"          / REFLECTING, &
         "southern"         / REFLECTING, &
         "northern"         / REFLECTING, &
         "bottomer"         / PERIODIC, &
         "topper"           / PERIODIC)

    ! physics settings
    physics => Dict( &
         "problem"          / EULER, &
         "gamma"            / GAMMA, &
         "dpmax"            / 1.0)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
         "fluxtype"         / KT, &
         "order"            / LINEAR, &
!         "variables"        / CONSERVATIVE, &
         "variables"        / PRIMITIVE, &
         "limiter"          / VANLEER)

    ! rotating frame for a sphere
    rotframe => Dict( &
         "stype"            / ROTATING_FRAME, &
         "gparam"           / GPAR, &
         "issphere"         / 1, &
         "x"                / 0.0, &
         "y"                / 0.0)

    ! cooling in infrared
    cooling => Dict( &
         "stype"            / PLANET_COOLING, &
         "output/Qcool"     / 1, &
         "output/T_s"       / 1, &
         "output/RHO_s"     / 1, &
         "distance"         / RD, &
         "output/P_s"       / 1, &
         "albedo"           / ALBEDO, &
         "mu"               / MU, &
         "intensity"        / INTENSITY, &
         "dz"               / HEIGHT, &
         "T_0"              / T_0, &
         "cvis"             / 0.1, &
         "gacc"             / GACC, &
         "gamma"            / GAMMA)

    ! heating by a star
    heating => Dict( &
         "stype"            / PLANET_HEATING, &
         "output/Qstar"     / 1, &
         "distance"         / RD, &
         "year"             / PLANET_YEAR, &
         "theta0"           / THETA0, &
         "phi0"             / PHI0, &
         "omegasun"         / FSUN, &
!          "R_planet"       / GPAR, &
         "albedo"           / ALBEDO, &
         "mu"               / MU, &
         "intensity"        / INTENSITY,&
         "dz"               / HEIGHT, &
!          "a_eff"          / 5.35e-7, &
         "cvis"             / 0.1,&
         "gacc"             / GACC, &
         "gamma"            / GAMMA)

    sources => Dict( &
         "rotframe"         / rotframe, &
         "cooling"          / cooling, &
         "heating"          / heating)

    ! time discretization settings
    timedisc => Dict( &
         "method"           / MODIFIED_EULER, &
         "order"            / 3, &
         "cfl"              / 0.4, &
         "stoptime"         / TSIM, &
         "dtlimit"          / 1.0E-15, &
         "maxiter"          / 1000000)

    datafile => Dict(&
         "fileformat"       / VTK, &
         "filename"         / (TRIM(ODIR) // TRIM(OFNAME)), &
         "count"            / ONUM)

    config => Dict( &
         "mesh"             / mesh, &
         "physics"          / physics, &
         "boundary"         / boundary, &
         "fluxes"           / fluxes, &
         "timedisc"         / timedisc, &
         "sources"          / sources, &
!         "logfile"          / logfile, &
         "datafile"         / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),             INTENT(IN)    :: Physics
    CLASS(mesh_base),                INTENT(IN)    :: Mesh
    CLASS(marray_compound), POINTER, INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL              :: theta1,phi1,theta2,phi2,vtheta,vphi,xlen,SIGMA
    REAL              :: R0,P1,PHISTART,THETASTART
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM) :: dv
    !------------------------------------------------------------------------!
    ! some starting constants
    PHISTART   = 0.0
    THETASTART = PI/3.
    P1 = P0-40.0E+2
    R0 = 0.4*GPAR
    SIGMA   = 1.0E+5

    ! different starting condition
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      p%density%data1d(:)   = P0/GACC
      p%pressure%data1d(:)  = Physics%Constants%RG* &
        p%density%data1d(:)*T_0/(MU*(1.+(GAMMA-1.)/GAMMA))
      !------------------------------------------------------------------------!
      !------------------- random velocity distribution -----------------------!
      CALL RANDOM_NUMBER (dv)
      p%data4d(:,:,:,Physics%XVELOCITY) = &
           p%data4d(:,:,:,Physics%XVELOCITY) + (dv(:,:,:,1)-0.5)*2.
      p%data4d(:,:,:,Physics%YVELOCITY) = &
           p%data4d(:,:,:,Physics%YVELOCITY) + (dv(:,:,:,2)-0.5)*2.
    END SELECT

    !------------------------------------------------------------------------!
    !--------------------- rotframe test ------------------------------------!
!     !    (velocity strips along latitude)
!     xlen = ABS(Mesh%xmax-Mesh%xmin)
!     WHERE ((Mesh%bcenter(:,:,1).GT.(Mesh%xmin+0.37*xlen)).AND. &
!          (Mesh%bcenter(:,:,1).LT.(Mesh%xmin+0.4*xlen)))
!        Timedisc%pvar(:,:,Physics%YVELOCITY) = -80.0
!     ELSEWHERE ((Mesh%bcenter(:,:,1).LT.(Mesh%xmin+(PI-0.37*xlen))).AND. &
!          (Mesh%bcenter(:,:,1).GT.(Mesh%xmin+(PI-0.4*xlen))))
!        Timedisc%pvar(:,:,Physics%YVELOCITY) = +80.0
!     ELSEWHERE ((Mesh%bcenter(:,:,1).LT.(Mesh%xmin+(PI-0.2*xlen))).AND. &
!          (Mesh%bcenter(:,:,1).GT.(Mesh%xmin+(PI-0.23*xlen))))
!        Timedisc%pvar(:,:,Physics%YVELOCITY) = -80.0
!     ELSEWHERE ((Mesh%bcenter(:,:,1).GT.(Mesh%xmin+0.2*xlen)).AND. &
!          (Mesh%bcenter(:,:,1).LT.(Mesh%xmin+0.23*xlen)))
!        Timedisc%pvar(:,:,Physics%YVELOCITY) = +80.0
!     END WHERE
!
    !------------------------------------------------------------------------!
    !---------------------- low density area --------------------------------!
!     WHERE ((Mesh%bcenter(:,:,1).LE.(THETASTART+R0/GPAR)).AND.&
!           (Mesh%bcenter(:,:,2).GE.(PHISTART-R0/GPAR/SIN(THETASTART))).AND.&
!           (Mesh%bcenter(:,:,1).GE.(THETASTART-R0/GPAR)).AND.&
!           (Mesh%bcenter(:,:,2).LE.(PHISTART+R0/GPAR/SIN(THETASTART))))
!        ! behind the shock front
!        Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
!     ELSEWHERE
!       Timedisc%pvar(:,:,Physics%PRESSURE)  = SIGMA*&
!               Timedisc%pvar(:,:,Physics%DENSITY)**GAMMA
!
!     END WHERE
    !------------------------------------------------------------------------!

    CALL Physics%Convert2Conservative(pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 2D planetary")

  END SUBROUTINE InitData

END PROGRAM planet2d
