!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: riemann1d.f90                                                     #
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
!> Acoustic Wave
!! \author Manuel Jung
!!
!! References:
!! [1] Bryson, S., and Levy, D.: On the Total Variation of High-Order
!!     Semi-Discrete Central Schemes for Conservation Laws
!!     Journal of Scientific Compiting, Vol 27, Nos. 1-3, June 2006
!----------------------------------------------------------------------------!
PROGRAM acousticwave
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
  USE solutions
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  REAL, PARAMETER :: TSIM    = 0.18
  REAL, PARAMETER :: GAMMA   = 1.4
  ! mesh settings
  INTEGER, PARAMETER :: XRES = 1600       ! x-resolution
  INTEGER, PARAMETER :: YRES = 1         ! y-resolution
  ! output parameters
  CHARACTER(LEN=256), PARAMETER &        ! output file name
                     :: OFNAME = 'acousticwave'
  INTEGER, PARAMETER :: ONUM = 100
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL Init(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               timedisc, fluxes
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
               "inum" / XRES, &             ! resolution in x and            !
               "jnum" / YRES, &             !   y direction                  !
               "xmin" / 0.0, &
               "xmax" / 1.0, &
               "ymin" / 0.0, &
               "ymax" / 1.0)

    physics => Dict(&
       "problem" / EULER2D, &
       "gamma"   / GAMMA)                   ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict(&
             "fluxtype"  / HLLC, &
             "order"     / LINEAR, &
             "variables" / PRIMITIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! time discretization settings
    timedisc => Dict( &
           "method"   / SSPRK, &
           "cfl"      / 0.4, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-10, &
           "maxiter"  / 100000000)

    ! boundary conditions
    boundary => Dict( &
           "western"  / NO_GRADIENTS, &
           "eastern"  / NO_GRADIENTS, &
           "southern" / NO_GRADIENTS, &
           "northern" / NO_GRADIENTS)

    ! initialize data input/output
    datafile => Dict( &
           "fileformat" / XDMF, &
           "filename"   / (TRIM(OFNAME)), &
           "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE Init(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    WHERE(Mesh%bccart(:,:,1).LE.0.1)
      Timedisc%pvar(:,:,Physics%DENSITY) = 3.857143
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 2.629369
      Timedisc%pvar(:,:,Physics%PRESSURE) = 10.3333
    ELSEWHERE
      Timedisc%pvar(:,:,Physics%DENSITY) = 1. + 0.2 *SIN(50.*Mesh%bccart(:,:,1))
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
      Timedisc%pvar(:,:,Physics%PRESSURE) = 1.
    END WHERE

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: Acoustic wave")

  END SUBROUTINE

END PROGRAM acousticwave
