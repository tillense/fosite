!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gresho.f90                                                        #
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
!> Gresho Vortex
!! \author Manuel Jung
!!
!! References:
!! [1] Liska, R. and Wendroff, B. (2003). "Comparison of several difference
!!     schemes on 1D and 2D test problems for the euler equations"
!!     SIAM J. SCI. COMPUT., Vol. 25, No. 3, pp. 995–1017
!! [2] P. Gresho (1990). "On the theory of semi-implicit projection methods
!!     for viscous incompressible flow and its implementation via
!!     finite-element method that also introduces a nearly consistent mass
!!     matrix: Part 2: Applications"
!!     Int. J. Numer. Methods Fluids, 11 (1990), pp. 621–659.
!! [3] M. S. Sahota, P. J. O’Rourke, and M. C. Cline (200). "Assesment of
!!     Vorticity and Angular Momentum Errors in CHAD for the Gresho Problem"
!!     Technical report LA-UR-00-2217, Los Alamos National Laboratory,
!!     Los Alamos, NM, 2000.
!----------------------------------------------------------------------------!
PROGRAM gresho
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 3.0      ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats
  ! initial condition (dimensionless units)
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES = 100         ! y-resolution
  ! physics settings
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'gresho'
  !--------------------------------------------------------------------------!
  REAL               :: sigma
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)
  CALL MakeConfig(Sim, Sim%config)
  CALL SetupFosite(Sim)
  CALL Init(Sim%Mesh,Sim%Physics,Sim%Timedisc)
  CALL RunFosite(Sim)
  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : ASINH
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes, sources, rotframe
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    !mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / (-1.), &
           "xmax"     / 1., &
           "ymin"     / (-1.), &
           "ymax"     / 1.)


    ! boundary conditions
    boundary => Dict(&
      "western" / PERIODIC, &
      "eastern" / PERIODIC, &
      "southern" / PERIODIC, &
      "northern" / PERIODIC)

    physics => Dict(&
      "problem"   / EULER2D, &
      "gamma"     / GAMMA)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
              "fluxtype"  / HLLC, &
              "order"     / LINEAR, &
              "variables" / PRIMITIVE, &
              "limiter"   / VANLEER, &               ! on cartesian mesh, but is
              "theta"     / 1.2)                     ! optional parameter for limiter !

    NULLIFY(sources)

    ! time discretization settings
    timedisc => Dict(&
         "method"   / MODIFIED_EULER, &
         "order"    / 3, &
         "cfl"      / 0.4, &
         "stoptime" / TSIM, &
         "dtlimit"  / 1.0E-10, &
         "maxiter"  / 10000000, &
!          "output/fluxes" / 1, &
         "output/iangularmomentum" / 1, &
         "output/rothalpy" / 1)

    ! initialize data input/output
     datafile => Dict(&
          "fileformat" / XDMF, &
          "filename"   / TRIM(OFNAME), &
          "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)

    ! add sources terms
    IF (ASSOCIATED(sources)) &
        CALL SetAttr(config, "sources", sources)


  END SUBROUTINE MakeConfig


  SUBROUTINE Init(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: vphi,p
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: ephi
    INTEGER           :: k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!

    WHERE(Mesh%radius%bcenter.LT.0.2)
      vphi = 5.*Mesh%radius%bcenter
      p = 5. + 12.5*Mesh%radius%bcenter**2
    ELSEWHERE(Mesh%radius%bcenter.LT.0.4)
      vphi = 2. - 5.*Mesh%radius%bcenter
      p = 9. - 4.*LOG(0.2) + 12.5*Mesh%radius%bcenter**2 - 20.*Mesh%radius%bcenter &
          + 4.*LOG(Mesh%radius%bcenter)
    ELSEWHERE
      vphi = 0.
      p = 3. + 4. * LOG(2.)
    END WHERE

    ! curvilinear components of azimuthal unit vector
    ! (maybe with respect to shifted origin)
    ! from ephi = ez x er = ez x posvec/radius = ez x (rxi*exi + reta*eeta)/r
    !             = rxi/r*(ez x exi) + reta/r*(ez x eeta) = rxi/r*eeta - reta/r*exi
    ! because (ez,exi,eeta) is right handed orthonormal set of basis vectors
    ephi(:,:,1) = -Mesh%posvec%bcenter(:,:,2)/Mesh%radius%bcenter(:,:)
    ephi(:,:,2) = Mesh%posvec%bcenter(:,:,1)/Mesh%radius%bcenter(:,:)

    Timedisc%pvar(:,:,Physics%DENSITY) = 1.
    Timedisc%pvar(:,:,Physics%XVELOCITY) = vphi * ephi(:,:,1)
    Timedisc%pvar(:,:,Physics%YVELOCITY) = vphi * ephi(:,:,2)
    Timedisc%pvar(:,:,Physics%PRESSURE) = p

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: Gresho Vortex")

  END SUBROUTINE Init
END PROGRAM gresho
