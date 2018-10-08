!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: shear.f03                                                         #
!#                                                                           #
!# Copyright (C) 2008-2018                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> \test Test shear flow
!! \author Jannes Klee
!----------------------------------------------------------------------------!
PROGRAM RTI
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: GN         = 1.0            ! grav. constant [GEOM]  !
  INTEGER, PARAMETER :: UNITS      = GEOMETRICAL
  ! simulation parameter
  REAL, PARAMETER    :: OMEGA      = 1.0
  REAL, PARAMETER    :: SIGMA0     = 1.0
  REAL, PARAMETER    :: TSIM       = 10./OMEGA
  REAL, PARAMETER    :: GAMMA      = 2.0
  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter     !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN
  INTEGER, PARAMETER :: XRES       = 128
  INTEGER, PARAMETER :: YRES       = 128
  INTEGER, PARAMETER :: ZRES       = 1
  REAL               :: DOMAINX    = 320.0
  REAL               :: DOMAINY    = 320.0
  INTEGER, PARAMETER :: FARGO      = 3              ! 3 = Shearingbox        !
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM       = 100
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sbox"
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  !--------------------------------------------------------------------------!
!  TAP_PLAN(1)

  ALLOCATE(Sim)

  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
  CALL Sim%Run()

  CALL Sim%Finalize()
  DEALLOCATE(Sim)

!  TAP_CHECK(.TRUE.,"Simulation finished")
!  TAP_DONE

  CONTAINS

  !> Set all configurations and safe it in a dictionary
  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)            :: Sim
    TYPE(Dict_TYP), POINTER  :: config
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: mesh, physics, boundary, datafile, &
                                sources, timedisc, fluxes, shear
    !------------------------------------------------------------------------!
    INTENT(INOUT)            :: Sim
    REAL                     :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
    !------------------------------------------------------------------------!
    DOMAINX    = DOMAINX*GN*SIGMA0/(OMEGA*OMEGA)
    DOMAINY    = DOMAINY*GN*SIGMA0/(OMEGA*OMEGA)
    XMIN       = -0.5*DOMAINX
    XMAX       = +0.5*DOMAINX
    YMIN       = -0.5*DOMAINY
    YMAX       = +0.5*DOMAINY
    ZMIN       = 0.0
    ZMAX       = 0.0

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
                "omega"       / OMEGA, &
                "fargo"       / FARGO, &
                "output/rotation" / 0, &
                "output/volume"   / 0, &
                "output/bh"   / 0, &
                "output/dl"   / 0  &
                )

    ! physics settings
    physics =>  Dict(&
                "problem"     / EULER2D, &
                "gamma"       / GAMMA, &
                "units"       / UNITS &
                )

    ! flux calculation and reconstruction method
    fluxes => Dict( &
              "fluxtype"      / KT, &
              "order"         / LINEAR, &
              "variables"     / PRIMITIVE, &
              "limiter"       / VANLEER &
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
                "western"     / PERIODIC, &
                "eastern"     / PERIODIC, &
                "southern"    / PERIODIC, &
                "northern"    / PERIODIC, &
                "bottomer"    / REFLECTING, &
                "topper"      / REFLECTING &
                )

    ! shearing box fictious forces
    shear => Dict( &
                "stype"           / SHEARBOX, &
                "output/accel_x"  / 0, &
                "output/accel_y"  / 0 &
                )

    ! collect sources in dictionary
    sources => Dict( &
              "shearing"      / shear)

    ! time discretization settings
    timedisc => Dict( &
              "method"        / MODIFIED_EULER, &
              "order"         / 3, &
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

    ! collect all above dicts in the configuration dict
    config => Dict( &
              "mesh"          / mesh, &
              "physics"       / physics, &
              "boundary"      / boundary, &
              "fluxes"        / fluxes, &
              "sources"       / sources, &
              "timedisc"      / timedisc, &
              "datafile"      / datafile)
  END SUBROUTINE MakeConfig


  !> Set initial conditions
  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                        INTENT(OUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                      :: y0
    !------------------------------------------------------------------------!
    ! Test for shearing and boundary module

    pvar(:,:,:,Physics%XVELOCITY) = 0.0
    pvar(:,:,:,Physics%YVELOCITY) = -Mesh%Q*Mesh%bcenter(:,:,:,1)*Mesh%Omega
    pvar(:,:,:,Physics%PRESSURE)  = 1e-1

    WHERE ((ABS(Mesh%bcenter(:,:,:,1)).LT.0.15*DOMAINX.AND. &
            ABS(Mesh%bcenter(:,:,:,2)).LT.0.15*DOMAINY))
      pvar(:,:,:,Physics%DENSITY)   = SIGMA0
    ELSEWHERE
      pvar(:,:,:,Physics%DENSITY)  = SIGMA0*1e-2
    END WHERE

    CALL Physics%Convert2Conservative(Mesh,pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Shearing patch")
  END SUBROUTINE InitData
END PROGRAM RTI
