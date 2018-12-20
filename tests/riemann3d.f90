!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: riemann3d.f90                                                     #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
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
!> \test 3D Riemann problem
!! \author Tobias Illenseer
!!
!! References:
!! - \cite langseth2000 Langseth, J. O., LeVeque, R. J.: A wave propagation method for three-
!!     dimensional hyperbolic conservation laws, J. Comput. Phys., 165(1) (2000)
!!     DOI: 10.1006/jcph.2000.6606
!----------------------------------------------------------------------------!
PROGRAM riemann3d
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 0.7        ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL ! geometry
  INTEGER, PARAMETER :: XRES  = 150        ! x-resolution
  INTEGER, PARAMETER :: YRES  = 1         ! y-resolution
  INTEGER, PARAMETER :: ZRES  = 100        ! y-resolution
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 7           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'riemann3d'
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
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP),POINTER  :: config
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
               "knum" / ZRES, &
               "xmin" / 0.0, &
               "xmax" / 1.5, &
               "ymin" / (-PI), &
               "ymax" / PI, &
               "zmin" / 0.0, &
               "zmax" / 1.0, &
             "gparam" / 1.0)

    ! physics settings
    physics => Dict("problem" / EULER, &
              "gamma"   / 1.4, &           ! ratio of specific heats        !
              "mu"      / 0.029)         ! mean molecular weight          !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! boundary conditions
    boundary => Dict(&
               "western"  / AXIS, &
               "eastern"  / REFLECTING, &
               "southern" / PERIODIC, &
               "northern" / PERIODIC, &
               "bottomer" / REFLECTING, &
               "topper"   / REFLECTING)

   ! time discretization settings
   timedisc => Dict("method"   / MODIFIED_EULER, &
              "order"    / 3, &
              "cfl"      / 0.4, &
              "stoptime" / TSIM, &
              "dtlimit"  / 1.0E-9, &
              "maxiter"  / 1000000)

    ! initialize data input/output
    datafile => Dict( &
            "fileformat" / VTK, &
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
    USE geometry_base_mod
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)         :: Physics
    CLASS(mesh_base)            :: Mesh
    CLASS(timedisc_base)        :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: rmax,x0,y0,z0
    REAL              :: P_in, P_out
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! radius and cartesian position of the bubble
    rmax = 0.2
    x0 = 0.0
    y0 = 0.0
    z0 = 0.4
    
    ! inside and ambient pressure
    P_in  = 5.0
    P_out = 1.0

    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_euler) ! non-isothermal HD
      ! constant density
      pvar%density%data1d(:)  = 1.0
      ! vanishing velocities
      pvar%velocity%data1d(:) = 0.0
      ! pressure pulse
      WHERE (((Mesh%bccart(:,:,:,1)-x0)**2+(Mesh%bccart(:,:,:,2)-y0)**2 &
             +(Mesh%bccart(:,:,:,3)-z0)**2).LT.rmax**2)
        pvar%pressure%data3d(:,:,:) = P_in
      ELSEWHERE
        pvar%pressure%data3d(:,:,:) = P_out
      END WHERE
    CLASS DEFAULT
      CALL Physics%Error("riemann3d::InitData","physics not supported")
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    CALL Physics%Info(" DATA-----> initial condition: Spherical pressure discontinuity between walls")
  END SUBROUTINE InitData

END PROGRAM riemann3d
