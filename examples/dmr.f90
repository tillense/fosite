!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: KHI.f90                                                           #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
!> Program and data initialization for the Double Mach Reflection Problem
!! References:
!! [1] Woodward, P. and Colella, P. (1984). The numerical simulation of
!!     two-dimensional fluid flow with strong shocks"
!!     Journal of Computational Physics, vol. 54, April 1984, p. 115-173.
!! [2] Stone, James M. et al (2008). "Athena: A New Code for Astrophysical MHD"
!!     The Astrophysical Journal Supplement Series, Volume 178, Issue 1,
!!     pp. 137-177.
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE sources_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 0.25          ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4           ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: XRES = 260            ! resolution
  INTEGER, PARAMETER :: YRES = 80

  REAL, PARAMETER    :: ETA  = 0.

  INTEGER, PARAMETER :: ONUM = 100            ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME = 'dmr'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim%config)

!  CALL PrintDict(Sim%config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4),sgbc
    DOUBLE PRECISION, PARAMETER :: GN = 6.6742E-11     ! [m^3/kg/s^2]
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               sources, timedisc, fluxes, vis
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
           "inum" / XRES, &
           "jnum" / YRES, &
           "xmin" / 0.0, &
           "xmax" / 4.0, &
           "ymin" / 0.0, &
           "ymax" / 1.0, &
           "output/volume" / 1, &
           "output/bary"   / 1,&
           "output/bh"    / 1,&
           "output/dl"    / 1 &
    )

    ! physics settings
    physics => Dict( &
              "problem" / EULER2D, &
              "gamma"   / GAMMA &           ! ratio of specific heats        !
    )

    ! flux calculation and reconstruction method
    fluxes => Dict( &
             "fluxtype"  / HLLC, &
             "order"     / LINEAR, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / VANLEER, &       ! one of: minmod, monocent,...   !
             "theta"     / 1.2 &            ! optional parameter for limiter !
    )

    ! boundary conditions
    boundary => Dict( &
               "western"  / FIXED, &
               "eastern"  / NO_GRADIENTS, &
               "southern" / CUSTOM, &
               "northern" / DMR &
    )

    ! viscosity source term
    vis => Dict( &
          "stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / ETA, &
          "bulkconst" / (-2./3.*ETA), &
          "output/dynvis" / 0, &
          "output/stress" / 1, &
          "output/kinvis" / 0, &
          "output/bulkvis" / 0 &
    )

    sources => Dict("." / 0)
    IF (ETA.GT.TINY(ETA)) THEN
        sources => Dict("vis" / vis)
    END IF

    ! time discretization settings
    timedisc => Dict( &
               "method"   / SSPRK, &
               "cfl"      / 0.4, &
               "stoptime" / TSIM, &
               "dtlimit"  / 1.0E-10, &
               "maxiter"  / 10000000, &
               "tol_rel" / 1.E-2, &
               "tol_abs" / (/ 0., 1.E-4, 1.E-4, 0. /), &
               "output/pressure" / 1, &
               "output/density" / 1, &
               "output/xvelocity" / 1, &
               "output/yvelocity" / 1 &
    )

    ! initialize data input/output
    datafile => Dict(&
        "fileformat" / XDMF, &
        "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
        "count"      / ONUM &
    )

    config => Dict( &
             "mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
             "datafile" / datafile &
    )

  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    REAL              :: rho0, rho1, u0, u1, v0, v1, P0, P1, X0
    REAL              :: xlen
    INTEGER           :: n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! downstream
    rho0 = 8.0
    u0   =  8.25*COS(PI/6.)
    v0   = -8.25*SIN(PI/6.)
    P0   = 116.5

    ! upstream
    rho1 = 1.4
    u1   = 0.0
    v1   = 0.0
    P1   = 1.0

    X0 = 1./6.

    WHERE (Mesh%bcenter(:,:,1).LT.(X0+Mesh%bcenter(:,:,2)/SQRT(3.)))
       Timedisc%pvar(:,:,Physics%DENSITY) = rho0
       Timedisc%pvar(:,:,Physics%XVELOCITY) = u0
       Timedisc%pvar(:,:,Physics%YVELOCITY) = v0
       Timedisc%pvar(:,:,Physics%PRESSURE) = P0
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY) = rho1
       Timedisc%pvar(:,:,Physics%XVELOCITY) = u1
       Timedisc%pvar(:,:,Physics%YVELOCITY) = v1
       Timedisc%pvar(:,:,Physics%PRESSURE) = P1
    END WHERE

    ! Boundary conditions

    IF(GetType(Timedisc%boundary(SOUTH)).EQ.CUSTOM) THEN
    DO i=Mesh%IMIN,Mesh%IMAX
      IF(Mesh%bcenter(i,Mesh%JMIN-1,1).LT.X0) THEN
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%DENSITY) = CUSTOM_FIXED
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%XVELOCITY) = CUSTOM_FIXED
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%YVELOCITY) = CUSTOM_FIXED
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%PRESSURE) = CUSTOM_FIXED
        DO j=1,Mesh%GNUM
          Timedisc%boundary(SOUTH)%data(i,j,Physics%DENSITY) = rho0
          Timedisc%boundary(SOUTH)%data(i,j,Physics%XVELOCITY) = u0
          Timedisc%boundary(SOUTH)%data(i,j,Physics%YVELOCITY) = v0
          Timedisc%boundary(SOUTH)%data(i,j,Physics%PRESSURE) = P0
        END DO
      ELSE
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%DENSITY) = CUSTOM_REFLECT
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%XVELOCITY) = CUSTOM_REFLECT
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%YVELOCITY) = CUSTOM_REFLNEG
        Timedisc%boundary(SOUTH)%cbtype(i,Physics%PRESSURE) = CUSTOM_REFLECT
      END IF
    END DO
    END IF

    IF(GetType(Timedisc%boundary(WEST)).EQ.FIXED) THEN
    Timedisc%Boundary(WEST)%fixed(:,Physics%DENSITY)   = .TRUE.
    Timedisc%Boundary(WEST)%fixed(:,Physics%XVELOCITY) = .TRUE.
    Timedisc%Boundary(WEST)%fixed(:,Physics%YVELOCITY) = .TRUE.
    Timedisc%Boundary(WEST)%fixed(:,Physics%PRESSURE)  = .TRUE.
    DO i=1,Mesh%GNUM
      Timedisc%boundary(WEST)%data(i,:,Physics%DENSITY) = rho0
      Timedisc%boundary(WEST)%data(i,:,Physics%XVELOCITY) = u0
      Timedisc%boundary(WEST)%data(i,:,Physics%YVELOCITY) = v0
      Timedisc%boundary(WEST)%data(i,:,Physics%PRESSURE) = P0
    END DO
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
            "Double Mach Reflection")

  END SUBROUTINE InitData

END PROGRAM Init
