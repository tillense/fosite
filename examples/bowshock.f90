!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: dmr.f90                                                           #
!#                                                                           #
!# Copyright (C) 2014                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                       #
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
!> Program and data initialization for a bow shock.
!! References:
!! [1] Developed by Jake Simon and David Nidever (University of Virginia) for
!!     Athena3D.
!!     http://www.astro.virginia.edu/VITA/ATHENA/bow-shock.html
!----------------------------------------------------------------------------!
PROGRAM bowshock
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
  REAL, PARAMETER    :: TSIM  = 10.0          ! simulation time
  REAL, PARAMETER    :: GAMMA = 5./3.         ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: XRES = 200*3          ! resolution
  INTEGER, PARAMETER :: YRES = 300*3

  INTEGER, PARAMETER :: ONUM = 100            ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME = 'bowshock'
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
           "xmin" / (-1.0), &
           "xmax" / 1.0, &
           "ymin" / (-1.5), &
           "ymax" / 1.5, &
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
               "southern" / NO_GRADIENTS, &
               "northern" / NO_GRADIENTS &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: R
    REAL              :: R0, X0, Y0
    REAL              :: xlen
    INTEGER           :: n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    ! Size and position of the high density sphere
    R0 = 0.0625
    X0 = -0.75
    Y0 = 0.

    R = SQRT((Mesh%bccart(:,:,1)-X0)**2+(Mesh%bccart(:,:,2)-Y0)**2)

    WHERE(R.LT.R0)
       Timedisc%pvar(:,:,Physics%DENSITY) = 100.0
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY) = 0.01
    END WHERE

    Timedisc%pvar(:,:,Physics%PRESSURE) = 0.0015
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    ! There is a uniform supersonic x velocity in all regions except for within
    ! the overdense region and to the right of this overdense region.
    ! Specifically, vx = 0 for the region defined by r ≤ 0.0625 and also defined
    ! by x > -0.75, |y| ≤ 0.0625. Everywhere else, vx = 1.0.
    WHERE((R.LT.R0).OR.((Mesh%bccart(:,:,1).GT.X0).AND.(ABS(Mesh%bccart(:,:,2)).LT.R0)))
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    ELSEWHERE
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 1.
    END WHERE

    ! Boundary conditions
    IF(GetType(Timedisc%boundary(WEST)).EQ.FIXED) THEN
      Timedisc%Boundary(WEST)%fixed(:,:)   = .TRUE.
      DO i=1,Mesh%GNUM
        Timedisc%boundary(WEST)%data(i,:,:) &
          = Timedisc%pvar(Mesh%IMIN,:,:)
      END DO
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
            "Bow Shock")

  END SUBROUTINE InitData

END PROGRAM bowshock
