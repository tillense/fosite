!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: bindisk.f90                                                       #
!#                                                                           #
!# Copyright (C) 2010-2017                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Anna Feiler      <afeiler@astrophysik.uni-kiel.de>                        #
!# Roman Avramenko  <ravramenko@astrophysik.uni-kiel.de>                     #
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
!> Program and data initialization for circumbinary accretion disk as MMSN
!!
!! \author Jannes Klee
!! \author Tobias Illenseer
!! \author Anna Feiler
!! \author Roman Avramenko
!!
!! \example mmsn.f90
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite_mod
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: GN      = 1.0             ! Newtons grav. constant    !
  REAL, PARAMETER    :: AU      = 1.0             ! astronomical unit         !
  REAL, PARAMETER    :: MSUN    = 1.0             ! solar mass [kg]              !
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 10              ! simulation time [binary orbits]
  REAL, PARAMETER    :: VALPHA  = 1e-3            ! alpha viscosity parameter !
  ! 1. binary system
!  REAL, PARAMETER    :: MBH1    = 0.690*MSUN
!  REAL, PARAMETER    :: MBH2    = 0.203*MSUN
  REAL, PARAMETER    :: MBH1    = 0.4465*MSUN
  REAL, PARAMETER    :: MBH2    = 0.4465*MSUN
!  REAL, PARAMETER    :: EXCENT  = 0.159          ! excentricity              !
  REAL, PARAMETER    :: EXCENT  = 0.0            ! excentricity              !
  REAL, PARAMETER    :: SEMMA   = 0.224*AU       ! semi mayor axis           !
  ! 2. disk
  REAL, PARAMETER    :: Rgap    = 2.5*SEMMA      ! estimated gap size        !
  REAL               :: XSCALE  = 1.0            ! scaling factor
  REAL               :: HRATIO  = 0.05           ! scaling factor
  ! gas paramter
  REAL, PARAMETER    :: MU      = 2.35e-3        ! mean molecular mass [kg/mol] !
  REAL, PARAMETER    :: RG      = 8.31447        ! molar gas constant        !
  ! mesh settings
  REAL, PARAMETER    :: GPAR = AU                ! geometry scaling paramete !
  REAL, PARAMETER    :: RMIN = 1.5*SEMMA         ! inner radius of the disk  !
  REAL, PARAMETER    :: RMAX = 5.0*AU            ! outer radius of the grid  !
  INTEGER, PARAMETER :: XRES = 128               ! x-resolution              !
  INTEGER, PARAMETER :: YRES = 128               ! y-resolution              !
  INTEGER, PARAMETER :: ZRES = 1                 ! y-resolution              !
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 100               ! number of output time st  !
  CHARACTER(LEN=256), PARAMETER &
                     :: ODIR  = "./"
  CHARACTER(LEN=256), PARAMETER &                ! output data file name     !
                     :: OFNAME = 'mmsn'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  REAL               :: OMEGA,PERIOD
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)
CALL Sim%InitFosite()
CALL MakeConfig(Sim, Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc, Sim%Fluxes, Sim%Sources)
CALL Sim%Run()
CALL Sim%Finalize()
DEALLOCATE(Sim)

CONTAINS

   SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: Sim
    TYPE(Dict_TYP),POINTER       :: config
    TYPE(Dict_TYP),POINTER       :: mesh, physics, fluxes, boundary,grav, &
                                    sources, binary, vis, timedisc, datafile, &
                                    rotframe
    !------------------------------------------------------------------------!
    ! some derived simulation parameters
    OMEGA  = SQRT(GN*(MBH1+MBH2)/SEMMA)/SEMMA
    PERIOD = 2.*PI / OMEGA
!     OMEGA  = 0.0

    ! mesh settings
    ! stellar orbits must be inside the central hole of the mesh
    mesh => Dict( &
              "meshtype"        / MIDPOINT, &
              "geometry"        / LOGCYLINDRICAL, &
!              "geometry"        / CYLINDRICAL, &
              "inum"            / XRES, &
              "jnum"            / YRES, &
              "knum"            / ZRES, &
              "xmin"            / LOG(RMIN/GPAR), &
              "xmax"            / LOG(RMAX/GPAR), &
!              "xmin"            / (RMIN/GPAR), &
!              "xmax"            / (RMAX/GPAR), &
              "ymin"            / 0.0, &
              "ymax"            / (2.0*PI), &
              "zmin"            / 0.0, &
              "zmax"            / 0.0, &
              "use_fargo"       / 1, &
              "fargo"           / 1, &
              "decomposition"   / (/ -1, 1, 1/), &
              "gparam"          / GPAR)
    IF (OMEGA.GT.0.0) CALL SetAttr(mesh,"omega",OMEGA)

    ! physics settings
    physics => Dict( &
              "problem"         / EULER_ISOTHERM, &
              "units"           / GEOMETRICAL,&
              "output/fcsound"  / 1, &
              "output/bccsound" / 1)


    ! boundary conditions
    boundary => Dict( &
              "western"         / CUSTOM,&
              "eastern"         / CUSTOM,&
              "southern"        / PERIODIC, &
              "northern"        / PERIODIC, &
              "bottomer"        / REFLECTING,&
              "topper"          / REFLECTING &
              )


    ! numerical fluxes and reconstruction method
    fluxes => Dict( &
              "order"           / LINEAR, &
              "fluxtype"        / KT, &
              "variables"       / PRIMITIVE, &
              "limiter"         / VANLEER, &
              "theta"           / 1.2)


    ! viscosity source term
    vis => Dict( &
              "stype"           / VISCOSITY, &
              "vismodel"        / ALPHA, &
              "cvis"            / 0.1,&
              "dynconst"        / VALPHA, &
              "output/dynvis"   / 1)

    rotframe => Dict( &
              "stype"           / ROTATING_FRAME)

    ! gravitational acceleration due to binary system
    binary => Dict( &
              "gtype"           / POINTMASS_BINARY, &
              "mass1"           / MBH1, &
              "mass2"           / MBH2, &
              "mesh_rot"        / 1, &
!              "excentricity"    / EXCENT, &
              "output/binpos"   / 1, &
              "output/omega"    / 1, &
              "semimayoraxis"   / SEMMA)

    ! source term due to all gravity terms
    grav => Dict( &
              "stype"           / GRAVITY, &
              "binary"          / binary,&
!              "self/gtype"      / SPECTRAL, &
!              "self/green"      / 1, &
              "output/height"   / 1, &
              "output/potential"/ 1, &
              "output/accel"    / 1)


    sources => Dict( &
              "vis"             / vis, &
              "grav"            / grav)
!     IF (OMEGA.GT.0.0) CALL SetAttr(sources,"rotframe",rotframe)

    ! time discretization settings
    timedisc => Dict( &
              "method"          / SSPRK,&
              "cfl"             / 0.3, &
              "stoptime"        / (PERIOD*TSIM), &
              "dtlimit"         / 1.0E-40,&
              "tol_rel"         / 1.0E-3, &
              "tol_abs"         / (/ 1.0E-08, 1., 1.0E-08 /), &
              "rhstype"         / 1, &
              "maxiter"         / 100000000)

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"      / VTK, &
              "filename"        / (TRIM(ODIR) // TRIM(OFNAME)), &
              "count"           / ONUM)

    config => Dict( &
              "physics"         / physics, &
              "fluxes"          / fluxes, &
              "mesh"            / mesh, &
              "boundary"        / boundary, &
              "sources"         / sources, &
              "timedisc"        / timedisc, &
              "datafile"        / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc,Fluxes,Sources)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE sources_gravity_mod, ONLY : sources_gravity
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN)     :: Mesh
    CLASS(physics_base),  INTENT(INOUT)  :: Physics
    CLASS(timedisc_base), INTENT(INOUT)  :: Timedisc
    CLASS(fluxes_base),   INTENT(INOUT)  :: Fluxes
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp => null()
    CLASS(sources_gravity), POINTER :: gp => null()
    INTEGER           :: i
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    REAL              :: Sigma0
    REAL, DIMENSION(:,:,:), POINTER :: r, Sigma
    REAL, DIMENSION(:,:,:,:), POINTER :: r_faces
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX)    :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,6)  :: fcsound
    CHARACTER(LEN=32) :: info_str
    !------------------------------------------------------------------------!
    ! get gravitational acceleration
    IF (ALLOCATED(Sources)) &
      sp => Sources%GetSourcesPointer(GRAVITY)
    IF (.NOT.ASSOCIATED(sp)) &
      CALL Physics%Error("mmsn::InitData","no gravity term initialized")

    ! set some pointers for convenience
    r => Mesh%RemapBounds(Mesh%radius%bcenter)
    r_faces => Mesh%RemapBounds(Mesh%radius%faces)

    ! Set sound speed
    bccsound = HRATIO*SQRT((MBH1+MBH2)*Physics%constants%GN/r(:,:,:))
    DO i =1,6
      fcsound(:,:,:,i) = HRATIO*SQRT((MBH1+MBH2)*Physics%constants%GN/r_faces(:,:,:,i))
    END DO

    ! set isothermal sound speeds
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_eulerisotherm)
      CALL phys%SetSoundSpeeds(Mesh,bccsound)
      CALL phys%SetSoundSpeeds(Mesh,fcsound)
    CLASS DEFAULT
      ! abort
      CALL phys%Error("InitData","physics not supported")
    END SELECT

    ! boundary conditions
    ! custom boundary conditions at western boundary if requested
    SELECT TYPE(bwest => Timedisc%Boundary%boundary(WEST)%p)
    CLASS IS (boundary_custom)
      CALL bwest%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER/))
    END SELECT
    SELECT TYPE(beast => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS (boundary_custom)
      CALL beast%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER/))
    END SELECT

    ! choose Sigma0 in a way that 2% of the mass of the binary is contained in 30 au
!    Sigma0     = 0.0008*(MBH1+MBH2)*SQRT(30.0/AU)/(4.*PI)
!    Sigma0     = 0.02*(MBH1+MBH2)/(4.*PI*(-(30.0/GPAR)**(-0.5)+(RMIN/GPAR)**(-0.5)))
!    Sigma0     = 0.02*(MBH1+MBH2)*SQRT(30.0/AU)/(4.*PI)
    Sigma0 = 0.00311
    SELECT TYPE (pvar => Timedisc%pvar)
    CLASS IS(statevector_eulerisotherm)
      ! 1. set surface density and initialize gravitational acceleration
      pvar%density%data3d(:,:,:) = Fgap(r(:,:,:),Rgap)*Sigma0*0.1*XSCALE*r(:,:,:)**(-1.5)
!       pvar%velocity%data4d(:,:,:,:) = 0.0
!       IF (Mesh%FARGO.GT.0) Timedisc%w(:,:) = 0.0
!       ! get conservative variables
!       CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
      ! 2. azimuthal velocity: balance initial radial gravitational acceleration
      !    with centrifugal acceleration
      CALL gp%UpdateGravity(Mesh,Physics,Fluxes,pvar,0.0,0.0)
      pvar%velocity%data4d(:,:,:,1:Physics%VDIM) = &
        Timedisc%GetCentrifugalVelocity(Mesh,Physics,Fluxes,Sources,(/0.,0.,1./),gp%accel%data4d)
      ! transform velocities to rotating frame
      ! ATTENTION: this works only if the second velocity is the azimuthal velocity
      pvar%velocity%data4d(:,:,:,2) = pvar%velocity%data4d(:,:,:,2) - Mesh%OMEGA*r(:,:,:)
    CLASS DEFAULT
      CALL Physics%Error("mmsn::InitData","unsupported state vector")
    END SELECT

    ! get conservative variables
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    IF (Mesh%fargo%GetType().EQ.2) &
       Timedisc%w(:,:) = SQRT(Physics%constants%GN*(MBH1+MBH2)/Mesh%radius%bcenter(:,Mesh%JMIN,:))

    ! print some information on stdout
    CALL Mesh%Info(" DATA-----> initial condition: " // "MMSN - Setup")
    WRITE (info_str, '(ES10.2)') XSCALE
    CALL Mesh%Info("                    disk mass: " // TRIM(info_str) // " MMSN")

  END SUBROUTINE InitData

  ELEMENTAL REAL FUNCTION Fgap(R, Rgap)
    IMPLICIT NONE
    REAL, INTENT(IN) :: R, Rgap

    Fgap = 1./(1. + exp(-(R-Rgap)/(0.1*Rgap)))

  END FUNCTION

END PROGRAM Init
