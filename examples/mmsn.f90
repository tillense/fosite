!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
  REAL, PARAMETER    :: GN      = 1.0             ! Newtons grav. constant    !
  REAL, PARAMETER    :: AU      = 1.0             ! astronomical unit         !
  REAL, PARAMETER    :: MSUN    = 1.0             ! solar mass [kg]              !
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 100            ! simulation time [binary orbits]
  REAL, PARAMETER    :: VALPHA  = 1e-3           ! alpha viscosity parameter !
  ! 1. binary system
  REAL, PARAMETER    :: MBH1    = 0.690*MSUN
  REAL, PARAMETER    :: MBH2    = 0.203*MSUN
!  REAL, PARAMETER    :: MBH1    = 0.4465*MSUN
!  REAL, PARAMETER    :: MBH2    = 0.4465*MSUN
  REAL, PARAMETER    :: EXCENT  = 0.159          ! excentricity              !
!  REAL, PARAMETER    :: EXCENT  = 0.0           ! excentricity              !
  REAL, PARAMETER    :: SEMMA   = 0.224*AU       ! semi mayor axis           !
  ! 2. disk
  REAL, PARAMETER    :: Rgap    = 2.5*SEMMA      ! estimated gap size        !
  REAL               :: XSCALE  = 1.0            ! scaling factor
  REAL               :: HRATIO  = 0.02           ! scaling factor
  ! gas paramter
  REAL, PARAMETER    :: MU      = 2.35e-3        ! mean molecular mass [kg/mol] !
  REAL, PARAMETER    :: RG      = 8.31447        ! molar gas constant        !
  ! mesh settings
  REAL, PARAMETER    :: GPAR = AU                ! geometry scaling paramete !
  REAL, PARAMETER    :: RMIN = 1.5*SEMMA         ! inner radius of the disk  !
  REAL, PARAMETER    :: RMAX = 5.0*AU            ! outer radius of the grid  !
  INTEGER, PARAMETER :: XRES = 64               ! x-resolution              !
  INTEGER, PARAMETER :: YRES = 64               ! y-resolution              !
  INTEGER, PARAMETER :: ZRES = 1                 ! y-resolution              !
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 100               ! number of output time st  !
  CHARACTER(LEN=256), PARAMETER &
                     :: ODIR  = "./"
  CHARACTER(LEN=256), PARAMETER &                ! output data file name     !
                     :: OFNAME = 'mmsn'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  REAL               :: CSISO,OMEGA,PERIOD
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)
CALL Sim%InitFosite()
CALL MakeConfig(Sim, Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc, Sim%Fluxes, Sim%Sources, &
              Sim%Timedisc%pvar, Sim%Timedisc%cvar)
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
                                    pmass
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=9)  :: geo_str,r1_str,r2_str,d_str
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    ! some derived simulation parameters
    OMEGA  = SQRT(GN*(MBH1+MBH2)/SEMMA)/SEMMA
    PERIOD = 2.*PI / OMEGA
    CSISO = HRATIO*SQRT((MBH1+MBH2)*GN/1.0)        ! r(:,:,:)
    OMEGA  = 0.0

    ! mesh settings
    ! stellar orbits must be inside the central hole of the mesh
    mesh => Dict( &
              "meshtype"        / MIDPOINT, &
              "geometry"        / LOGCYLINDRICAL, &
!              "geometry"        / CYLINDRICAL, &
              "inum"            / XRES, &
              "jnum"            / YRES, &
              "knum"            / ZRES, &
              "xmin"            / (RMIN/GPAR), &
              "xmax"            / (RMAX/GPAR), &
              "ymin"            / 0.0, &
              "ymax"            / (2.0*PI), &
              "zmin"            / 0.0, &
              "zmax"            / 0.0, &
              "omega"           / OMEGA, &
              "fargo"           / 0, &
              "decomposition"   / (/ -1, 1/), &
              "gparam"          / GPAR)


    ! physics settings
    physics => Dict( &
              "problem"         / EULER2D_ISOTHERM, &
!              "cs"              / CSISO, &
              "units"           / GEOMETRICAL,&
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


    ! gravitational acceleration due to binary system
    binary => Dict( &
              "gtype"           / POINTMASS_BINARY, &
              "mass1"           / MBH1, &
              "mass2"           / MBH2, &
              "mesh_rot"        / 1, &
              "excentricity"    / EXCENT, &
              "output/binpos"   / 1, &
              "output/omega"    / 1, &
              "semimayoraxis"   / SEMMA)

    pmass => Dict( &
              "gtype"           / POINTMASS, &
              "mass"            / MBH1)


    ! source term due to all gravity terms
    grav => Dict( &
              "stype"           / GRAVITY, &
!              "binary"         / binary,&
              "pointmass"       / pmass,&
!              "self/gtype"      / SPECTRAL, &
!              "self/green"      / 1, &
!              "energy"          / 0, &
              "output/height"   / 0, &
              "output/potential" / 1, &
              "output/accel"    / 1)


    sources => Dict( &
              "grav"            / grav,&
              "vis"             / vis )


    ! time discretization settings
    timedisc => Dict( &
              "method"          / SSPRK,&
              "cfl"             / 0.3, &
              "stoptime"        / (PERIOD*TSIM), &
              "dtlimit"         / 1.0E-40,&
              "tol_rel"         / 1.0E-3, &
              "tol_abs"         / (/ 1.0E-16, 1., 1.0E-16 /), &
              "rhstype"         / 1, &
!              "output/bflux"    / 1,&
!              "output/rhs"      / 1,&
!              "output/xmomentum"/ 1,&
!              "output/ymomentum"/ 1,&
              "maxiter"         / 100000000)

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"      / VTK, &
              "unit"            / 5555, &
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


  SUBROUTINE InitData(Mesh,Physics,Timedisc,Fluxes,Sources,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN)     :: Mesh
    CLASS(physics_base),  INTENT(INOUT)  :: Physics
    CLASS(timedisc_base), INTENT(INOUT)  :: Timedisc
    CLASS(fluxes_base),   INTENT(INOUT)  :: Fluxes
    CLASS(sources_base),  POINTER        :: Sources
    TYPE(Dict_TYP), POINTER              :: config
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                          INTENT(OUT)    :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp
    INTEGER           :: i,j,dir,ig,im
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    REAL              :: mdisk_temp
    REAL              :: s,q, BINMASS, Sigma0
    REAL, DIMENSION(:,:,:), POINTER :: r, Sigma
    REAL, DIMENSION(:,:,:,:), POINTER :: r_faces
    LOGICAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: sum_mask
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX)    :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,6)  :: fcsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX)    :: Temp
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,6)  :: Temp_faces
    CHARACTER(LEN=32) :: info_str
    !------------------------------------------------------------------------!
    ! set some pointers for convenience
    r => Mesh%RemapBounds(Mesh%radius%bcenter)
    r_faces => Mesh%RemapBounds(Mesh%radius%faces)
    ! pointer to density array
    Sigma => Mesh%RemapBounds(pvar(:,:,:,Physics%DENSITY))

    ! chose Sigma0 in a way that 2% of the mass of the binary is contained in 30 au
!    Sigma0     = 0.0008*(MBH1+MBH2)*SQRT(30.0/AU)/(4.*PI)
!    Sigma0     = 0.02*(MBH1+MBH2)/(4.*PI*(-(30.0/GPAR)**(-0.5)+(RMIN/GPAR)**(-0.5)))
!    Sigma0     = 0.02*(MBH1+MBH2)*SQRT(30.0/AU)/(4.*PI)
    Sigma0 = 0.00311
    Sigma = Fgap(r(:,:,:),Rgap)*Sigma0*0.1*XSCALE*r(:,:,:)**(-1.5)

    ! Set sound speed
    bccsound = HRATIO*SQRT((MBH1+MBH2)*Physics%constants%GN/r(:,:,:))
    DO i =1,6
      fcsound(:,:,:,i) = HRATIO*SQRT((MBH1+MBH2)*Physics%constants%GN/r_faces(:,:,:,i))
    END DO

    ! set isothermal sound speeds
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler2dit)
      CALL phys%SetSoundSpeeds(Mesh,bccsound)
      CALL phys%SetSoundSpeeds(Mesh,fcsound)
    CLASS DEFAULT
      ! abort
      CALL phys%Error("InitData","physics not supported")
    END SELECT

    ! 2. azimuthal velocity: balance initial radial acceleration with centrifugal acceleration
    ! get gravitational acceleration
    sp => Sources
    DO
      IF (ASSOCIATED(sp).EQV..FALSE.) RETURN
      SELECT TYPE(sp)
      CLASS IS(sources_gravity)
        gp => sp
        EXIT
      END SELECT
      sp => sp%next
    END DO

    IF (ASSOCIATED(sp)) THEN
       CALL gp%UpdateGravity(Mesh,Physics,Fluxes,pvar,0.0,0.0)
    ELSE
       CALL Sim%Error("InitData","no gravity term initialized")
    END IF

    ! ATTENTION: Don't use GetCentrifugalVelocity without the optional acceleration array!
    ! This would yield undefined data, because GetCentrifugalVelocity calls ComputeRHS
    ! which calls CenterBoundary. Since the FARFIELD boundary conditions are not
    ! initialized at this stage (see below), the result is undefined.
    pvar(:,:,:,Physics%XVELOCITY:Physics%XVELOCITY+Physics%DIM-1) = &
        Timedisc%GetCentrifugalVelocity(Mesh,Physics,Fluxes,Sources,(/0.,0.,1./),gp%accel)


    ! transform velocities to rotating frame
    Timedisc%pvar(:,:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,:,Physics%YVELOCITY) &
        - Mesh%OMEGA*r(:,:,:)

    ! get conservative variables
    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar,Timedisc%cvar)

    IF (Mesh%FARGO.EQ.2) &
       Timedisc%w(:,:) = SQRT(Physics%constants%GN*(MBH1+MBH2)/r(:,Mesh%JMIN,:))

    ! boundary conditions
    ! custom boundary conditions at western boundary if requested
    IF ((Timedisc%Boundary%boundary(WEST)%p%GetType()).EQ.CUSTOM) THEN
       Timedisc%Boundary%boundary(WEST)%p%cbtype(:,:,Physics%DENSITY) = CUSTOM_NOGRAD
       Timedisc%Boundary%boundary(WEST)%p%cbtype(:,:,Physics%XVELOCITY) = CUSTOM_OUTFLOW
       Timedisc%Boundary%boundary(WEST)%p%cbtype(:,:,Physics%YVELOCITY) = CUSTOM_KEPLER
    END IF
    IF ((Timedisc%Boundary%boundary(EAST)%p%GetType()).EQ.CUSTOM) THEN
       Timedisc%Boundary%boundary(EAST)%p%cbtype(:,:,Physics%DENSITY) = CUSTOM_NOGRAD
       Timedisc%Boundary%boundary(EAST)%p%cbtype(:,:,Physics%XVELOCITY) = CUSTOM_OUTFLOW
       Timedisc%Boundary%boundary(EAST)%p%cbtype(:,:,Physics%YVELOCITY) = CUSTOM_KEPLER
    END IF

    ! print some information on stdout
    CALL Mesh%Info(" DATA-----> initial condition:   " // "MMSN - Setup")
    WRITE (info_str, '(ES8.2)') XSCALE
    CALL Mesh%Info("  disk mass:   " // TRIM(info_str) // " MMSN")

  END SUBROUTINE InitData

  ELEMENTAL REAL FUNCTION Fgap(R, Rgap)
    IMPLICIT NONE
    REAL, INTENT(IN) :: R, Rgap

    Fgap = 1./(1. + exp(-(R-Rgap)/(0.1*Rgap)))

  END FUNCTION

END PROGRAM Init
