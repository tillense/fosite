!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sblintheo.f90                                                     #
!#                                                                           #
!# Copyright (C) 2015-2018                                                   #
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
!> \test Program and data initialization for the linear theory test in a shearingsheet.
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
!! gravity. The standard test implementation here solves the isothermal case
!! with gravity until the first maximum in the amplitude is reached.
!!
!! Test for the following modules:
!!    1. shearing boundaries (\ref boundary_shearing_mod )
!!    2. fictious forces (\ref sources_shearbox_mod )
!!    3. selfgravity (\ref gravity_sboxspectral_mod )
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
!! python script solving these cases can be found in folder tools/scripts.
!----------------------------------------------------------------------------!
PROGRAM sblintheo
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: GN         = 1.0            ! grav. constant [GEOM]  !
  ! make a shearingsheet simulation and denote the direction of shearing
  ! (applies boundaries, ficitious forces and fargo automatically)
  INTEGER, PARAMETER :: SHEARSHEET_DIRECTION = 1
  ! simulation parameter
  REAL, PARAMETER    :: OMEGA      = 1.0            ! rotation at fid. point !
  REAL, PARAMETER    :: SIGMA0     = 1.0            ! mean surf.dens.        !
  REAL, PARAMETER    :: DELSIGMA   = SIGMA0*5.0e-4  ! disturbance            !
  REAL, PARAMETER    :: TSIM       = 4.54845485/OMEGA ! simulation time (first maximum) !
  REAL, PARAMETER    :: GAMMA      = 2.0            ! dep. on vert. struct.  !
  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter     !
  ! set (initial) speed of sound using the Toomre criterion and
  ! assuming marginal stabily; in non-isothermal simulations this is
  ! used to set the initial pressure (see InitData) whereas in isothermal
  ! simulations this is constant throughout the simulation
  REAL, PARAMETER    :: SOUNDSPEED = PI*GN*SIGMA0/OMEGA
  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN
  INTEGER, PARAMETER :: XRES       = 128            ! amount of cells in x-  !
  INTEGER, PARAMETER :: YRES       = 128            ! y-direction (rho/phi)  !
  INTEGER, PARAMETER :: ZRES       = 1
  REAL               :: DOMAINX    = 40.0           ! domain size [GEOM]     !
  REAL               :: DOMAINY    = 40.0           ! domain size [GEOM]     !
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM       = 10
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sblintheo"
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  REAL               :: maximum
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)

  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
  CALL Sim%Run()

  ! check amplitude
  maximum = MAXVAL(Sim%Timedisc%pvar%data4d(:,:,:,Sim%Physics%DENSITY)) - SIGMA0

  CALL Sim%Finalize()
  DEALLOCATE(Sim)

  TAP_CHECK_CLOSE(maximum, 0.02395931, 0.0003, "First max. < 3e-4 deviation")
  TAP_DONE

CONTAINS

  !> \public Configuration of the initial data.
  SUBROUTINE MakeConfig(Sim,config)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP), POINTER :: config
    !--------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: mesh,physics,fluxes,&
                               grav,sources,timedisc,datafile
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
                "problem"     / EULER_ISOTHERM, &
                "cs"          / SOUNDSPEED, &
                "units"       / GEOMETRICAL &
                )

    ! mesh settings
    mesh =>     Dict(&
                "meshtype"    / MIDPOINT, &
                "geometry"    / MGEO, &
                "shearingbox" / SHEARSHEET_DIRECTION, &
                "inum"        / XRES, &
                "jnum"        / YRES, &
                "knum"        / ZRES, &
                "xmin"        / XMIN, &
                "xmax"        / XMAX, &
                "ymin"        / YMIN, &
                "ymax"        / YMAX, &
                "zmin"        / ZMIN, &
                "zmax"        / ZMAX, &
                "omega"       / OMEGA, &
                "decomposition"/ (/ 1, -1, 1/), &
                "output/rotation" / 0, &
                "output/volume"   / 0, &
                "output/bh"   / 0, &
                "output/dl"   / 0,  &
                "Qshear"      / Q &
                )

    ! fluxes settings
    fluxes =>   Dict(&
                "order"       / LINEAR, &
                "fluxtype"    / KT, &
                "variables"   / PRIMITIVE, &
                "limiter"     / VANLEER &
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

    ! add gravitational sources
    sources =>  Dict(&
                "grav"        / grav &
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
                  "sources"     / sources, &
                  "timedisc"    / timedisc, &
                  "datafile"    / datafile &
                  )
  END SUBROUTINE MakeConfig

  !> \public Initializes the data.
  !> Set initial conditions
  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(marray_compound), POINTER, INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! local variable declaration
    REAL              :: kx, ky
    !------------------------------------------------------------------------!

    !---------------------- linear theory test ------------------------------!
    kx = -2*(2*PI/DOMAINX)
    ky = 2*PI/DOMAINY

    ! initial condition
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_eulerisotherm)
      IF(Mesh%shear_dir.EQ.2)THEN
        p%density%data3d(:,:,:) = SIGMA0 &
          + DELSIGMA*COS(kx*Mesh%bcenter(:,:,:,1) + ky*Mesh%bcenter(:,:,:,2))
        p%velocity%data2d(:,1) = 0.0
        p%velocity%data4d(:,:,:,2) = -Q*OMEGA*Mesh%bcenter(:,:,:,1)
      ELSE IF(Mesh%shear_dir.EQ.1)THEN
        p%density%data3d(:,:,:) = SIGMA0 &
          + DELSIGMA*COS(kx*Mesh%bcenter(:,:,:,2) - ky*Mesh%bcenter(:,:,:,1))
        p%velocity%data4d(:,:,:,1) = Q*OMEGA*Mesh%bcenter(:,:,:,2)
        p%velocity%data2d(:,2) = 0.0
      END IF
    TYPE IS(statevector_euler) ! non-isothermal HD
      IF(Mesh%shear_dir.EQ.2)THEN
        p%density%data3d(:,:,:) = SIGMA0 &
          + DELSIGMA*COS(kx*Mesh%bcenter(:,:,:,1) + ky*Mesh%bcenter(:,:,:,2))
        p%velocity%data2d(:,1) = 0.0
        p%velocity%data4d(:,:,:,2) = -Q*OMEGA*Mesh%bcenter(:,:,:,1)
      ELSE IF(Mesh%shear_dir.EQ.1)THEN
        p%density%data3d(:,:,:) = SIGMA0 &
          + DELSIGMA*COS(kx*Mesh%bcenter(:,:,:,2) - ky*Mesh%bcenter(:,:,:,1))
        p%velocity%data4d(:,:,:,1) = Q*OMEGA*Mesh%bcenter(:,:,:,2)
        p%velocity%data2d(:,2) = 0.0
      END IF
      p%pressure%data3d(:,:,:) = SOUNDSPEED**2 * SIGMA0 / GAMMA
    CLASS DEFAULT
      CALL Physics%Error("shear::InitData","only (non-)isothermal HD supported")
    END SELECT

    CALL Physics%Convert2Conservative(pvar,cvar)

    !------------------------------------------------------------------------!

    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Linear Theory Test - Shearingsheet")
  END SUBROUTINE InitData
END PROGRAM sblintheo
