!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: vortex3d.f90                                                      #
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
!> 3D isentropic vortex
!! \author Tobias Illenseer
!!
!! References:
!! [1] Yee, H. C. et al.: Low-dissipative high-order shock-capturing methods
!!     using characteristic-based filters, J. Comput. Phys. 150 (1999), 199-238
!!     DOI: 10.1006/jcph.1998.6177
!----------------------------------------------------------------------------!
PROGRAM vortex3d
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
  REAL, PARAMETER    :: TSIM    = 100.0    ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHOINF  = 1.       ! ambient density
  REAL, PARAMETER    :: PINF    = 1.       ! ambient pressure
  REAL, PARAMETER    :: VSTR    = 5.0      ! nondimensional vortex strength
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL ! geometry
  INTEGER, PARAMETER :: XRES = 1          ! x-resolution
  INTEGER, PARAMETER :: YRES = 400         ! y-resolution
  REAL, PARAMETER    :: RMAX = 5.0         ! outer radius
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter     !
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'vortex3d' 
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
    TYPE(Dict_TYP), POINTER :: physics, boundary, datafile, logfile, &
                               timedisc, fluxes, mesh
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
               "inum" / XRES, &
               "jnum" / YRES, &
               "xmin" / (-1.0), &
               "xmax" / 1.0, &
               "ymin" / 0.0, &
               "ymax" / RMAX)

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA)                  ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! boundary conditions
    boundary => Dict("western"  / NO_GRADIENTS, &
               "eastern"  / NO_GRADIENTS, &
               "southern" / AXIS, &
               "northern" / NO_GRADIENTS)

    ! time discretization settings
    timedisc => Dict(&
         "method"   / MODIFIED_EULER, &
         "order"    / 3, &
         "cfl"      / 0.4, &
         "stoptime" / TSIM, &
         "dtlimit"  / (1.0E-8*TSIM), &
         "maxiter"  / 100000)

    ! initialize data input/output
    datafile => Dict("fileformat" / GNUPLOT, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "filecycles" / 0, &
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
    INTEGER           :: i,j
    REAL              :: r,r2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! axis distance in 2.5D rotationally symmetric coordinates
          ! is always given by the scale factor hz
          r  = Mesh%hz%bcenter(i,j)
          r2 = r*r
          ! density
          ! ATTENTION: there's a factor of 1/PI missing in the density
          ! formula  eq. (3.3) in [1]
          Timedisc%pvar(i,j,Physics%DENSITY) = RHOINF * (1.0 - (GAMMA-1.0)*VSTR**2 / &
               (8*PI**2*GAMMA) * EXP(1.-r2) )**(1./(GAMMA-1.))
          ! pressure
          Timedisc%pvar(i,j,Physics%PRESSURE) = PINF &
               * (Timedisc%pvar(i,j,Physics%DENSITY)/RHOINF)**GAMMA
          ! rotational velocity
          Timedisc%pvar(i,j,Physics%ZVELOCITY) = 0.5*VSTR/PI*r*EXP(0.5*(1.-r2))
          ! other velocities
          Timedisc%pvar(i,j,Physics%XVELOCITY) = 0.
          Timedisc%pvar(i,j,Physics%YVELOCITY) = 0.
        END DO
    END DO
    
    ! set specific angular momentum for EULER3D_ROTAMT
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) THEN
       Timedisc%pvar(:,:,Physics%ZVELOCITY) = &
            Timedisc%pvar(:,:,Physics%ZVELOCITY)*Mesh%hz%bcenter(:,:)
    END IF
    
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 3D isentropic vortex")
  END SUBROUTINE InitData

END PROGRAM vortex3d
