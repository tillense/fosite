!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gauss3d.f90                                                       #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> 3D Gaussian pressure or density pulse with and without rotation
!! \author Tobias Illenseer
!!
!! \cite illenseer2009
!----------------------------------------------------------------------------!
PROGRAM gauss3d
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameter
  REAL, PARAMETER     :: TSIM     = 0.3      ! simulation time
  REAL, PARAMETER     :: GAMMA    = 1.4      ! ratio of specific heats
  REAL, PARAMETER     :: CSISO    = &
!                                     0.0      ! non-isothermal simulation
                                    1.127    ! isothermal simulation
                                             !   with CSISO as sound speed
  ! initial condition (dimensionless units)
  REAL, PARAMETER     :: RHO0     = 1.0      ! ambient density
  REAL, PARAMETER     :: RHO1     = 1.0      ! peak density above RHO0
  REAL, PARAMETER     :: RWIDTH   = 0.06     ! half width of the Gaussian
  REAL, PARAMETER     :: P0       = 1.0      ! ambient pressure
  REAL, PARAMETER     :: P1       = 1.0      ! peak pressure above P0
  REAL, PARAMETER     :: PWIDTH   = 0.06     ! half width of the Gaussian
  REAL, PARAMETER     :: OMEGA0   = 0.0      ! angular velocity
  REAL, PARAMETER     :: ETA      = 0.0      ! dynamic viscosity (0.0 disables)
  ! location of the initial pulse (cartesian coordinates)
  REAL, PARAMETER     :: X0       = 0.0      ! x-position
  REAL, PARAMETER     :: Y0       = 0.0      ! y-position
  REAL, PARAMETER     :: Z0       = 0.0      ! z-position
  ! mesh settings
  INTEGER, PARAMETER  :: MGEO     = CARTESIAN! geometry
!   INTEGER, PARAMETER  :: MGEO     = CYLINDRICAL
!   INTEGER, PARAMETER  :: MGEO     = SPHERICAL
  INTEGER, PARAMETER  :: RES      = 30       ! x-/y-/z-resolution
  REAL, PARAMETER     :: RMAX     = 0.5      ! maximal radial extent
  REAL, PARAMETER     :: GPAR     = 0.8      ! geometry scaling parameter
  ! output parameters
  INTEGER, PARAMETER  :: ONUM     = 10       ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &            ! output data dir
                      :: ODIR     = './'
  CHARACTER(LEN=256), PARAMETER &            ! output data file name
                      :: OFNAME   = 'gauss3D'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  LOGICAL :: ok
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL Sim%Run()
  ok = .NOT.Sim%aborted
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

  TAP_CHECK(ok,"stoptime reached")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : Asinh, Acosh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP),POINTER  :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                 :: bc(6),xres,yres,zres
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, sources, &
                               timedisc, fluxes
    REAL                    :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1         = -RMAX
       x2         =  RMAX
       y1         = -RMAX
       y2         =  RMAX
       z1         = -RMAX
       z2         =  RMAX
       xres       = RES
       yres       = RES
       zres       = RES
       bc(WEST)   = NO_GRADIENTS
       bc(EAST)   = NO_GRADIENTS
       bc(SOUTH)  = NO_GRADIENTS
       bc(NORTH)  = NO_GRADIENTS
       bc(BOTTOM) = NO_GRADIENTS
       bc(TOP)    = NO_GRADIENTS
    CASE(CYLINDRICAL)
       x1 =  0.01
       x2 =  RMAX
       y1 =  0.0
       y2 = 2*PI
       z1 = -RMAX
       z2 =  RMAX
       xres       = RES / 2
       yres       = RES / 2
       zres       = RES
       bc(WEST)   = AXIS
       bc(EAST)   = AXIS
       bc(SOUTH)  = PERIODIC
       bc(NORTH)  = PERIODIC
       bc(BOTTOM) = NO_GRADIENTS
       bc(TOP)    = NO_GRADIENTS
    CASE(SPHERICAL)
       x1 = 0.01
       x2 = RMAX
       y1 = 0.0
       y2 = PI
       z1 = 0.0
       z2 = 2*PI
       xres       = RES / 2
       yres       = RES / 4
       zres       = RES / 2
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
       bc(BOTTOM)= PERIODIC
       bc(TOP)   = PERIODIC
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","geometry not supported for this test")
    END SELECT

    !mesh settings
    mesh => Dict( &
            "meshtype"  / MIDPOINT, &
            "geometry"  / MGEO, &
            "inum"      / xres, &
            "jnum"      / yres, &
            "knum"      / zres, &
            "xmin"      / x1,   &
            "xmax"      / x2,   &
            "ymin"      / y1,   &
            "ymax"      / y2,   &
            "zmin"      / z1,   &
            "zmax"      / z2,   &
!             "decomposition"/ (/ 1, 1, 2/), &
            "output/radius" / 1, &
            "gparam"    / GPAR)

    ! boundary conditions
    boundary => Dict( &
            "western"   / bc(WEST), &
            "eastern"   / bc(EAST), &
            "southern"  / bc(SOUTH), &
            "northern"  / bc(NORTH), &
            "bottomer"  / bc(BOTTOM), &
            "topper"    / bc(TOP))

    ! physics settings
    IF (CSISO.GT.TINY(CSISO)) THEN
      physics => Dict("problem" / EULER_ISOTHERM, &
                      "cs"      / CSISO)             ! isothermal speed of sound
    ELSE
      physics => Dict("problem"   / EULER, &
                      "gamma"     / GAMMA)           ! ratio of specific heats
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict( &
            "order"     / LINEAR, &
            "fluxtype"  / KT, &
            "variables" / PRIMITIVE, &
            "limiter"   / VANLEER, &
            "output/slopes" / 0)

    NULLIFY(sources)

    ! time discretization settings
    timedisc => Dict( &
            "method"    / MODIFIED_EULER, &
            "order"     / 3, &
            "tol_rel"   / 10.0, & ! > 1 disables adaptive step size
            "cfl"       / 0.4, &
            "stoptime"  / TSIM, &
            "dtlimit"   / 1.0E-8, &
            "maxiter"   / 10000000)

    datafile => Dict( &
            "fileformat" / VTK, &
            "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
            "count"      / ONUM)

    config => Dict( &
            "mesh"     / mesh, &
            "physics"  / physics, &
            "boundary" / boundary, &
            "fluxes"   / fluxes, &
            "datafile" / datafile, &
            "timedisc" / timedisc)

    IF (ASSOCIATED(sources)) &
        CALL SetAttr(config, "sources", sources)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    USE geometry_base_mod
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)         :: Physics
    CLASS(mesh_base)            :: Mesh
    CLASS(timedisc_base)        :: Timedisc
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX)   &
                                :: radius
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                                :: posvec
    !------------------------------------------------------------------------!
    INTENT(IN)                  :: Mesh,Physics
    INTENT(INOUT)               :: Timedisc
    !------------------------------------------------------------------------!
    IF (ANY((/ABS(X0),ABS(Y0),ABS(Z0)/).GT.TINY(Z0))) THEN
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,:,1) = X0
       posvec(:,:,:,2) = Y0
       posvec(:,:,:,3) = Z0
       CALL Mesh%geometry%Convert2Curvilinear(Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - posvec(:,:,:,:)
       ! compute its absolute value
       radius(:,:,:)   = SQRT(posvec(:,:,:,1)**2+posvec(:,:,:,2)**2+posvec(:,:,:,3)**2)
    ELSE
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:,:)   = Mesh%radius%bcenter(:,:,:)
       posvec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:)
    END IF

    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_eulerisotherm) ! isothermal HD
      ! Gaussian density pulse
      pvar%density%data3d(:,:,:)  = RHO0 + RHO1*EXP(-LOG(2.0) &
               * (radius(:,:,:)/RWIDTH)**2)
      ! vanishing velocities
      pvar%velocity%data1d(:) = 0.0
    TYPE IS(statevector_euler) ! non-isothermal HD
      ! constant density
      pvar%density%data1d(:)  = RHO0
      ! Gaussian pressure pulse
      pvar%pressure%data3d(:,:,:) = P0 + P1*EXP(-LOG(2.0) &
               * (radius(:,:,:)/RWIDTH)**2)
      ! vanishing velocities
      pvar%velocity%data1d(:) = 0.0
    CLASS DEFAULT
      CALL Physics%Error("gauss3d::InitData","unsupported physics")
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    CALL Mesh%Info(" DATA-----> initial condition: 3D Gaussian pulse")

  END SUBROUTINE InitData

END PROGRAM gauss3d
