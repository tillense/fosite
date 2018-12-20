!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: implosion.f90                                                     #
!#                                                                           #
!# Copyright (C) 2014                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> Program and data initialization for the Implosion Problem
!! References:
!! [1] Liska, R. and Wendroff, B. (2003). "Comparison of several difference
!!     schemes on 1D and 2D test problems for the euler equations"
!!     SIAM J. SCI. COMPUT., Vol. 25, No. 3, pp. 995–1017
!! [2] W. Hui, P. Li, and Z. Li. "A unified coordinate system for solving the
!!     two-dimensional euler equations"
!!     J. Comp. Phys., 153 (1999), pp. 596–637
!! [3] S. Chang, X. Wang, and C. Chow. "The space-time conservation element
!!     and solution element method: A new high resolution and genuinely 
!!     multidimensional paradigm for solving conservation laws"
!!     J. Comp. Phys., 160 (1999)
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
  REAL, PARAMETER    :: TSIM  = 5.0           ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4           ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: XRES = 200            ! resolution
  INTEGER, PARAMETER :: YRES = 200

  INTEGER, PARAMETER :: ONUM = 100            ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME = 'implosion'
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
           "xmax" / 1.0, &
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
               "western"  / REFLECTING, &
               "eastern"  / REFLECTING, &
               "southern" / REFLECTING, &
               "northern" / REFLECTING &
    )

    sources => Dict("." / 0)

    ! time discretization settings
    timedisc => Dict( &
               "method"   / SSPRK, &
               "cfl"      / 0.4, &
               "stoptime" / TSIM, &
               "dtlimit"  / 1.0E-10, &
               "maxiter"  / 10000000, &
               "tol_rel" / 1.E-3, &
               "tol_abs" / (/ 0., 1.E-4, 1.E-4, 0. /), &
               "output/pressure" / 1, &
               "output/density" / 1, &
               "output/xvelocity" / 1, &
               "output/yvelocity" / 1 &
    )

    ! initialize data input/output
    datafile => Dict(&
        "fileformat" / HDF, &
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
    REAL              :: rho0, rho1, P0, P1, X0
    REAL              :: xlen
    INTEGER           :: n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! bottom left
    rho0 = 0.125
    P0   = 0.5

    ! upper right
    rho1 = 1.0
    P1   = 1.0

    X0 = 0.5

    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    WHERE(Mesh%bcenter(:,:,1)+Mesh%bcenter(:,:,2).LT.X0)
       Timedisc%pvar(:,:,Physics%DENSITY) = rho0
       Timedisc%pvar(:,:,Physics%PRESSURE) = P0
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY) = rho1
       Timedisc%pvar(:,:,Physics%PRESSURE) = P1
    END WHERE

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
            "Implosion test")

  END SUBROUTINE InitData

END PROGRAM Init
