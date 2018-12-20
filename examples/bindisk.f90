!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: bindisk.f90                                                       #
!#                                                                           #
!# Copyright (C) 2010-2012                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Anna Feiler      <afeiler@astrophysik.uni-kiel.de>                        #
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
! Program and data initialization for circumbinary accretion disk
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
  USE functions
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
  REAL, PARAMETER :: GN      = 6.6742D-11     ! Newtons grav. constant       !
  REAL, PARAMETER :: CC      = 2.99792458D+8  ! speed of light               !
  REAL, PARAMETER :: RG      = 8.31447        ! molar gas constant           !
  REAL, PARAMETER :: MSUN    = 1.989D+30      ! solar mass [kg]              !
  REAL, PARAMETER :: RSUN    = 6.96342D8      ! solar radius [m]             !
  REAL, PARAMETER :: LSUN    = 3.85D+26       ! solar luminosity [W]         !
  REAL, PARAMETER :: AU      = 1.49597870691E+11 ! astronomical unit [m]     !
  REAL, PARAMETER :: PARSEC  = 0.5/PI*1.296E+6 * AU ! parsec [m]             !
  REAL, PARAMETER :: YEAR    = 3.15576E+7     ! Julian year [sec]            !
  REAL, PARAMETER :: SB      = 5.6704E-8      ! [W/m^2/K^4]

  ! simulation parameters
  REAL, PARAMETER :: TSIM    = 2D3*YEAR       ! simulation time              !
  INTEGER, PARAMETER :: ONUM = 100            ! number of output time steps  !
  CHARACTER(LEN=256), PARAMETER :: ODIR &     ! output directory             !
                                 = "./"        
  CHARACTER(LEN=256), PARAMETER &             ! output data file name        !
                     :: OFNAME = 'binary' 

  ! 1. binary system
  REAL, PARAMETER :: MS1     = 3.0 * MSUN     ! mass of primary component    !
  REAL, PARAMETER :: MS2     = 1.0 * MSUN     ! mass of secondary (<=MS1)    !
  REAL, PARAMETER :: EXCENT  = 0.05            ! excentricity                 !
  REAL, PARAMETER :: SEMMA   = 30.*AU         ! semi mayor axis              !
  ! 2. disk
  REAL, PARAMETER :: MDISK   = 0.01* MSUN     ! initial disk mass            !
  REAL, PARAMETER :: T_ISO   = 50.0           ! temperature (isothermal)     !  
  REAL, PARAMETER :: RMIN    = 1.5*SEMMA       ! inner radius of the disk     !
  ! initial condition: isothermal gaussian ring
  REAL, PARAMETER :: FWHM    = 100.0 * AU     ! FWHM of the initial ring     !
  REAL, PARAMETER :: R0      = 1.5*FWHM        !  location of the initial ring !
  REAL, PARAMETER :: RDISK   = 4.*FWHM        ! outer radius of the disk     !

  REAL, PARAMETER :: MU      = 2.35e-3         ! mean molecular mass [kg/mol] !
  REAL, PARAMETER :: CpCv    = 1.4             ! ratio of specific heats      !
  REAL, PARAMETER :: VALPHA  = 1e-2            ! alpha viscosity parameter    !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = LOGPOLAR        ! geometry of the mesh         !
  REAL, PARAMETER    :: GPAR = 0.1             ! geometry scaling parameter
  REAL, PARAMETER :: RMAX    = 3.0E+3*AU       ! outer radius of the grid     !    
  REAL, PARAMETER :: RHOMIN  = 1.e-5           ! minimum density on the mesh  !
  !!!!! ATTENTION: use even grid resolutions for VTK-output !!!!!!!!!!!!!!!!!!
  INTEGER, PARAMETER :: XRES = 80           ! x-resolution                 !
  INTEGER, PARAMETER :: YRES = 80           ! y-resolution                 !
  ! output file parameter
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)
 
  CALL MakeConfig(Sim, Sim%config)

  CALL SetupFosite(Sim)
 
  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,Sim%Fluxes, Sim%config)

  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

   SUBROUTINE MakeConfig(Sim, config)
    USE geometry_bipolar
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER &
                      :: config
    TYPE(Dict_TYP),POINTER &
                      :: mesh, physics, fluxes, boundary,grav,&
                         sources, binary, vis, timedisc, datafile
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=9)  :: geo_str,r1_str,r2_str,d_str
    INTEGER           :: i,j
    REAL              :: csiso
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!

 
    ! mesh settings
    ! stellar orbits have to be inside of the central hole of the mesh
     
     mesh => Dict("meshtype" / MIDPOINT, &
              "geometry" / LOGPOLAR, &
              "inum" / XRES, &
              "jnum" / YRES, &
              "xmin" / LOG(RMIN/GPAR), &
              "xmax" / LOG(RMAX/GPAR), &
              "ymin" / 0.0, &
              "ymax" / (2.0*PI), &
              "gparam" / GPAR)
             

    csiso = SQRT(RG/MU*T_ISO) ! isothermal sound speed  !
   
    ! physics settings
    physics => Dict(&
                     "problem" / EULER2D_ISOTHERM, &
                      "cs"      / csiso,& 
                      "units" / SI,&
                      "gamma" /  CpCv,&
                      "mu"    /MU, &
                      "output/bccsound" / 0)
                     

    ! boundary conditions
    boundary => Dict("western"  / NO_GRADIENTS,&
                     "eastern"  / NO_GRADIENTS,&
                     "southern" / PERIODIC, &
                     "northern" / PERIODIC)


    ! reconstruction method

    fluxes => Dict("order"     / LINEAR, &
                   "fluxtype"  / KT, &
                   "variables" / PRIMITIVE, & ! vars. to use for reconstruction!
                   "limiter"   / MINMOD, &    ! one of: minmod, monocent,...   !
                   "theta"     / 1.4)         ! optional parameter for limiter !


    ! viscosity source term
    vis => Dict("stype"      / VISCOSITY, &
                "vismodel"   / ALPHA_ALT, &
                "cvis"       / 0.1,&
                "dynconst"   / VALPHA)


    ! gravitational acceleration due to binary system
    binary => Dict("gtype" / POINTMASS_BINARY, &
                   "mass1" / MS1, &
                   "mass2" / MS2, &
                   "excentricity"  / EXCENT, &
                   "output/binpos" / 1, &
                   "output/omega"  / 1, &
                   "output/accel_bin" / 1, &
                   "output/height" / 1, &
                   "semimayoraxis" / SEMMA)


    ! source term due to all gravity terms
    grav => Dict( "stype"    / GRAVITY, &
                  "binary" / binary,&
                  "output/height"/ 1, &
                  "output/accel"/ 1)

   
    sources => Dict(&
              "grav" / grav,&     
              "vis" / vis )


    ! time discretization settings
    timedisc => Dict("method" /DORMAND_PRINCE,& 
               "cfl" / 0.3, &
               "stoptime" / TSIM, &
               "dtlimit" / 1.0E-6,&
               "tol_rel" / 1.0E-1, &
               "tol_abs" / (/ 1.0E-16, 1., 1.0E-16 /), &
               "output/bflux" / 1,&
               "output/rhs"/1,& 
               "maxiter" / 100000000)

    ! initialize data input/output
    datafile => Dict(&
!                 "fileformat" / HDF, &
!                 "fileformat" / GNUPLOT , &
              "fileformat" / XDMF, &
!               "unit" / 5555, &
               "filename" / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count" / ONUM)

    config => Dict("physics"  / physics, &
                   "fluxes"   / fluxes, &
                   "mesh"     / mesh, &
                   "boundary" / boundary, &
                   "sources"  / sources, &
                   "timedisc" / timedisc, &
                   "datafile" / datafile)
 
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc,Fluxes,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Timedisc_TYP) :: Timedisc
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Dict_TYP),POINTER &
                       :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,dir,ig,im
#ifdef PARALLEL
    INTEGER           :: ierror
    REAL              :: mdisk_local
#endif
    REAL              :: M_temp
    REAL              :: s, r, Sigma0
    CHARACTER(LEN=80) :: initial
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh
    INTENT(INOUT)      :: Timedisc,Physics
    !------------------------------------------------------------------------!

    ! set defaults
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHOMIN
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    ! initial surface density: Gaussian ring at R0 with mass MDISK
    ! between RMIN and RDISK
  
    s = 0.5 * FWHM / SQRT(2*LOG(2.0))  ! standard deviation of the initial Gaussian ring

    Sigma0  = 1.
    initial= 'gaussian ring'


   ! Calculate initial values for the primitive variables
!---------------------------------------------------------------------------------------------------!   
   ! Surface density:
   M_temp=0.0
          DO i = Mesh%IMIN,Mesh%IMAX
             DO j = Mesh%JMIN,Mesh%JMAX 

                r = Mesh%radius%bcenter(i,j)
                Timedisc%pvar(i,j,Physics%DENSITY)= Sigma0 * &
                                                    EXP(-0.5*((r-R0)/s)**2)+RHOMIN

                IF((r.LE.RDISK)) THEN
                   M_temp = M_temp + Timedisc%pvar(i,j,Physics%DENSITY)*Mesh%volume(i,j)
                END IF

             END DO
          END DO
     ! rescale density and pressure with correct disk mass
#ifdef PARALLEL
   mdisk_local = M_temp
   CALL MPI_AllReduce(mdisk_local,M_temp,1,DEFAULT_MPI_REAL,MPI_SUM, &
                     Mesh%comm_cart,ierror)          
#endif
   Timedisc%pvar(:,:,Physics%DENSITY)  = MDISK/M_temp*Timedisc%pvar(:,:,Physics%DENSITY)+RHOMIN
!---------------------------------------------------------------------------------------------------!
!velocity initial condition, resulting velocities are defined for curvlinear coordinates 
   CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

   Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = &
      GetCentrifugalVelocity(Timedisc,Mesh,Physics,Fluxes,(/0.,0.,1./)) 

!---------------------------------------------------------------------------------------------------!
  
 CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

!---------------------------------------------------------------------------------------------------!

    CALL Info(Mesh, " DATA-----> initial condition: " // &
          TRIM(initial))

  END SUBROUTINE InitData

END PROGRAM Init
