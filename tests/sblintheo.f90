!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sblintheo.f90                                                     #
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
!> \test Program and data initialization for the linear theory test in a SB.
!!
!! \author Jannes Klee
!!
!! References:
!!
!! - \cite gammie1996 Charles F. Gammie (1996). "Linear Theory of Magnetized,
!!     Viscous, Self-gravitating Gas Disks"
!!     The Astrophysical Journal 462: 725
!!
!! - \cite gammie2001 Charles F. Gammie (2001). "Nonlinear outcome of gravitational
!!     instability in cooling, gaseous disks"
!!     The Astrophysical Journal 553: 174-183
!!
!! - \cite paardekooper2012 Sijme-Jan Paardekooper (2012). "Numerical convergence
!!     in self-gravitating shearing sheet simulations and the stochastic nature
!!     of disc fragmentation"
!!     Monthly Notices of the Royal Astronomical Society, Volume 421, Issue 4,
!!     pp. 3286-3299.
!!
!! The test can be either done thermal or isothermal and with or without
!! gravity.
!!
!! Test for the following modules:
!!    1. shearing boundaries (\link boundary_shearing \endlink)
!!    2. fictious forces (\link sources_shearbox \endlink)
!!    3. selfgravity (\link gravity_sboxspectral \endlink)
!!
!!
!! Within the linear theory test the evolution of a shearing wave
!! \f[
!!    \Sigma = \Sigma_0 + \delta \Sigma\cos{\mathbf{k \cdot x}}
!! \f]
!! is examined. Additionally the background shear \f$ v_x = 0,
!! v_y = -q \Omega x\f$ is set up and either the speed of sound
!! (isothermal) or the pressure (thermal) is determined by the
!! Toomre Criterion:
!! \f[
!!    Q \sim 1 = \frac{\Omega c_s}{\pi G \Sigma_0} =
!!          \frac{\Omega \sqrt{\gamma \frac{P}{\Sigma_0}}}{\pi G \Sigma_0}.
!! \f]
!! There are analytical solutions for the thermal and isothermal case
!! \cite gammie2001 \cite paardekooper2012 \cite gammie1996 . A
!! python script solving cases can be found in folder tools.
!! \todo{upload python script!}
!----------------------------------------------------------------------------!
PROGRAM sblintheo
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: GN         = 1.0            ! grav. constant [GEOM]  !
  ! simulation parameter
  REAL, PARAMETER    :: OMEGA      = 1.0            ! rotation at fid. point !
  REAL, PARAMETER    :: SIGMA0     = 1.0            ! mean surf.dens.        !
  REAL, PARAMETER    :: DELSIGMA   = SIGMA0*5.0e-4  ! disturbance            !
  REAL, PARAMETER    :: TSIM       = 10.0/OMEGA     ! simulation time        !
  REAL, PARAMETER    :: GAMMA      = 2.0            ! dep. on vert. struct.  !
  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter     !
  REAL, PARAMETER    :: SOUNDSPEED = PI*GN*SIGMA0/OMEGA ! det. by. Toomre    !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN
  INTEGER, PARAMETER :: XRES       = 64            ! amount of cells in x-  !
  INTEGER, PARAMETER :: YRES       = 64            ! y-direction (rho/phi)  !
  INTEGER, PARAMETER :: ZRES       = 1
  REAL               :: DOMAINX    = 40.0           ! domain size [GEOM]     !
  REAL               :: DOMAINY    = 40.0           ! domain size [GEOM]     !
  ! fargo 0=off, 3=on (for SB)
  INTEGER, PARAMETER :: FARGO      = 3              ! 3 = Shearingbox        !
  ! number of output time steps
!  INTEGER, PARAMETER :: ONUM       = 1
  INTEGER, PARAMETER :: ONUM       = 100
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sblintheo"
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  REAL               :: maximum
  !--------------------------------------------------------------------------!

!TAP_PLAN(1)

ALLOCATE(Sim)

CALL Sim%InitFosite()
CALL MakeConfig(Sim, Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
CALL Sim%Run()

! search for the amplitude
!maximum = MAXVAL(Sim%Timedisc%pvar(:,:,:,Sim%Physics%DENSITY)) - SIGMA0
!TAP_CHECK_CLOSE(maximum, 0.0316463969505, 0.005, "Last max. < 0.005 deviation")

CALL Sim%Finalize()
DEALLOCATE(Sim)

!TAP_DONE

CONTAINS

  !> \public Configuration of the initial data.
  SUBROUTINE MakeConfig(Sim,config)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP), POINTER :: config
    !--------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: mesh,physics,fluxes,boundary,&
                               grav,shearingbox,sources,timedisc,&
                               datafile
    REAL                    :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
    !--------------------------------------------------------------------------!
    DOMAINX    = DOMAINX*GN*SIGMA0/(OMEGA*OMEGA)
    DOMAINY    = DOMAINY*GN*SIGMA0/(OMEGA*OMEGA)
    XMIN       = -0.5*DOMAINX
    XMAX       = +0.5*DOMAINX
    YMIN       = -0.5*DOMAINY
    YMAX       = +0.5*DOMAINY
    ZMIN       = 0.0
    ZMAX       = 0.0

    ! physics settings
    physics =>  Dict(&
                "problem"     / EULER2D_ISOTHERM, &
                "cs"          / SOUNDSPEED, &
                "units"       / GEOMETRICAL &
                )

    ! mesh settings
    mesh =>     Dict(&
                "meshtype"    / MIDPOINT, &
                "geometry"    / MGEO, &
                "inum"        / XRES, &
                "jnum"        / YRES, &
                "knum"        / ZRES, &
                "xmin"        / XMIN, &
                "xmax"        / XMAX, &
                "ymin"        / YMIN, &
                "ymax"        / YMAX, &
                "zmin"        / ZMIN, &
                "zmax"        / ZMAX, &
                "fargo"       / FARGO, &
                "shift_dir"   / 2, &
                "omega"       / OMEGA, &
                "decomposition"/ (/ 1, -1, 1/), &
                "output/rotation" / 0, &
                "output/volume"   / 0, &
                "output/bh"   / 0, &
                "output/dl"   / 0,  &
                "Q"           / Q &
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
                "western"     / SHEARING, &
                "eastern"     / SHEARING, &
                "southern"    / PERIODIC, &
                "northern"    / PERIODIC, &
 !               "western"     / PERIODIC, &
 !               "eastern"     / PERIODIC, &
 !               "southern"    / SHEARING, &
 !               "northern"    / SHEARING, &
                "bottomer"    / REFLECTING, &
                "topper"      / REFLECTING &
                )

    ! gravity settings (source term)
    grav =>     Dict(&
                "stype"               / GRAVITY, &
                "self/gtype"          / SBOXSPECTRAL, &
!                "output/accel"        / 1, &
                "self/output/phi"     / 0, &
                "self/output/accel_x" / 0, &
                "self/output/accel_y" / 0 &
                )

    ! shearing box fictious forces
    shearingbox => Dict(&
                "stype"           / SHEARBOX &
                )

    ! sources settings (contains source terms)
    sources =>  Dict(&
!                "grav"        / grav, &
                "shearing"    / shearingbox &
                )

    ! time discretization settings
    timedisc => Dict( &
              "method"        / DORMAND_PRINCE, &
              "cfl"           / 0.4, &
              "stoptime"      / TSIM, &
              "dtlimit"       / 1.0E-40, &
              "tol_rel"       / 0.1, &
              "maxiter"       / 100000000)

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"    / VTK, &
              "filename"      / (TRIM(ODIR) // TRIM(OFNAME)), &
              "count"         / ONUM)

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

  !> \public Initializes the data.
 SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                        INTENT(OUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! local variable declaration
    REAL              :: kx, ky
    REAL              :: SOUNDSPEED
    !------------------------------------------------------------------------!


    !---------------------- linear theory test ------------------------------!
    ! initial density
    kx = -2*(2*PI/DOMAINX)
    ky = 2*PI/DOMAINY
    IF(Mesh%WE_shear)THEN
      pvar(:,:,:,Physics%DENSITY)    = SIGMA0 + DELSIGMA*COS( &
                                              kx*Mesh%bcenter(:,:,:,1) + &
                                              ky*Mesh%bcenter(:,:,:,2))
      pvar(:,:,:,Physics%XVELOCITY)  = 0.0
      pvar(:,:,:,Physics%YVELOCITY)  = -Q*OMEGA*Mesh%bcenter(:,:,:,1)
    ELSE
      pvar(:,:,:,Physics%DENSITY)    = SIGMA0 + DELSIGMA*COS( &
                                              kx*Mesh%bcenter(:,:,:,2) - &
                                              ky*Mesh%bcenter(:,:,:,1))
      pvar(:,:,:,Physics%XVELOCITY)  = Q*OMEGA*Mesh%bcenter(:,:,:,2)
      pvar(:,:,:,Physics%YVELOCITY)  = 0.0
    END IF

    ! initial soundspeed (isothermal) or pressure (thermal)
    ! determined by Toomre-criterion
    SOUNDSPEED = PI*Physics%Constants%GN*SIGMA0/OMEGA
    IF (Physics%PRESSURE.GT.0) pvar(:,:,:,Physics%PRESSURE)   = &
              PI*PI*GN*GN*SIGMA0**3./(GAMMA*OMEGA**2)
    !------------------------------------------------------------------------!

    CALL Physics%Convert2Conservative(Mesh,pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Linear Theory Test - Shearingsheet")
  END SUBROUTINE InitData
END PROGRAM sblintheo
