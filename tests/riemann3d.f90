!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: riemann3d.f90                                                     #
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
!> 3D Riemann problem
!! \author Tobias Illenseer
!!
!! References:
!! [1] Langseth, J. O., LeVeque, R. J.: A wave propagation method for three-
!!     dimensional hyperbolic conservation laws, J. Comput. Phys., 165(1) (2000)
!!     DOI: 10.1006/jcph.2000.6606
!----------------------------------------------------------------------------!
PROGRAM riemann3d
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 0.7        ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL ! geometry
  INTEGER, PARAMETER :: XRES  = 100        ! x-resolution
  INTEGER, PARAMETER :: YRES  = 150        ! y-resolution
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 7           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'riemann3d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Simulation finished")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: z0,z1,r0,r1
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CYLINDRICAL, &
               "inum" / XRES, &
               "jnum" / YRES, &
               "xmin" / 0.0, &
               "xmax" / 1.0, &
               "ymin" / 0.0, &
               "ymax" / 1.5, &
             "gparam" / 1.0)

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / 1.4, &           ! ratio of specific heats        !
              "mu"      / 0.029)         ! mean molecular weight          !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! boundary conditions
    boundary => Dict("western"  / REFLECTING, &
               "eastern"  / REFLECTING, &
               "southern" / AXIS, &
               "northern" / REFLECTING)

   ! time discretization settings
   timedisc => Dict("method"   / MODIFIED_EULER, &
              "order"    / 3, &
              "cfl"      / 0.4, &
              "stoptime" / TSIM, &
              "dtlimit"  / 1.0E-9, &
              "maxiter"  / 1000000)

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: rmax,x0,y0
    REAL              :: P_in, P_out
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! radius and position of the bubble
    rmax = 0.2
    x0 = 0.0
    y0 = 0.4
    ! inside and ambient pressure
    P_in  = 5.0
    P_out = 1.0

    ! density
    Timedisc%pvar(:,:,Physics%DENSITY) = 1.0
    ! velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.
    ! pressure
    WHERE (((Mesh%bccart(:,:,1)-x0)**2+(Mesh%bccart(:,:,2)-y0)**2).LT.rmax**2)
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P_in
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%PRESSURE)  = P_out
    END WHERE

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: Spherical pressure discontinuity between walls")
  END SUBROUTINE InitData

END PROGRAM riemann3d
