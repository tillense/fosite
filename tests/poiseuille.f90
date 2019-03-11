!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: poiseuille.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006-2019                                                   #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
!> \test Poisseuille flow through a pipe
!! \author Bjoern Sperling
!! \author Tobias Illenseer
!!
!! \brief Program and data initialization for testing the law of Hagen-Poiseuille
!!
!----------------------------------------------------------------------------!
PROGRAM poiseuille
  USE fosite_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER :: TSIM  = 30.0                ! simulation time            !
  REAL, PARAMETER :: RE    = 4.375              ! Reynolds Number            !
  REAL, PARAMETER :: PIN   = 24.0               ! inflow pressure            !
  REAL, PARAMETER :: POUT  = 23.0               ! outflow pressure           !
  REAL, PARAMETER :: RHO0  = 1.0                ! initial density            !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!  INTEGER, PARAMETER :: MGEO = CARTESIAN
  INTEGER, PARAMETER :: XRES = 10               ! radial resolution          !
  INTEGER, PARAMETER :: YRES = 10               ! azimuthal resolution       !
  INTEGER, PARAMETER :: ZRES = 10               ! resolution along the tube  !
  REAL, PARAMETER    :: LTUBE = 10.0            ! length of the tube         !
  REAL, PARAMETER    :: RTUBE = 1.0             ! radius of the tube         !
  !--------------------------------------------------------------------------!
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'poiseuille'
  !--------------------------------------------------------------------------!
  CLASS(fosite),ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!
  ALLOCATE(SIM)
  CALL SIM%InitFosite()

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL Sim%Setup()

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL Sim%Run()
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, sources, &
                               vis, timedisc, fluxes
    REAL              :: x1,x2,y1,y2,z1,z2
    REAL              :: dvis, bvis
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! geometry dependent setttings
    SELECT CASE(MGEO)
    CASE(CYLINDRICAL)
      x1 = 0.0
      x2 = RTUBE
      y1 = 0.0
      y2 = 2*PI
      z1 = 0.0
      z2 = LTUBE
      bc(WEST)  = AXIS
      bc(EAST)  = NOSLIP
      bc(SOUTH) = PERIODIC
      bc(NORTH) = PERIODIC
      bc(BOTTOM)= FIXED
      bc(TOP)   = FIXED
    CASE(CARTESIAN)
      x1 = 0.0
      x2 = RTUBE
      y1 = 0.0
      y2 = RTUBE
      z1 = 0.0
      z2 = LTUBE
      bc(WEST)  = NOSLIP
      bc(EAST)  = NOSLIP
      bc(SOUTH) = NOSLIP
      bc(NORTH) = NOSLIP
      bc(BOTTOM)= FIXED
      bc(TOP)   = FIXED
    END SELECT

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "knum"     / ZRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "zmin"     / z1, &
           "zmax"     / z2)

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH), &
               "bottomer" / bc(BOTTOM), &
               "topper"   / bc(TOP))

    ! physics settings
    physics => Dict("problem" / EULER, &
              "gamma"   / 1.4)           ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! viscosity source term
    ! dynamic viscosity
    dvis = SQRT(0.25 * (PIN-POUT) * RTUBE**3 * RHO0 / LTUBE / RE)
    ! bulk viscosity
    bvis = -2./3. * dvis
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / dvis, &
          "bulkconst" / bvis)

    sources => Dict("vis" / vis)

    ! time discretization settings
    timedisc => Dict("method"   / MODIFIED_EULER, &
               "order"    / 3, &
               "cfl"      / 0.4, &
               "stoptime" / TSIM, &
               "dtlimit"  / 1.0E-8, &
               "maxiter"  / 1000000)

    ! initialize data input/output
    datafile => Dict( &
               "fileformat" / XDMF, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_base) :: Physics
    CLASS(Mesh_base)    :: Mesh
    CLASS(Timedisc_base):: Timedisc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
    Timedisc%pvar%data4d(:,:,:,Physics%DENSITY)   = RHO0
    Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = PIN + (POUT-PIN)/LTUBE &
         * Mesh%bccart(:,:,:,3)

    ! fixed boundary conditions
    ! inflow
    SELECT TYPE(bbottom => Timedisc%Boundary%boundary(BOTTOM)%p)
    CLASS IS (boundary_fixed)
       bbottom%data(:,:,:,Physics%DENSITY)    = RHO0
       bbottom%data(:,:,:,Physics%XVELOCITY)  = 0.0
       bbottom%data(:,:,:,Physics%YVELOCITY)  = 0.0
       bbottom%data(:,:,:,Physics%ZVELOCITY)  = 0.0
       bbottom%data(:,:,:,Physics%PRESSURE)   = PIN
       ! imposed density, pressure and tangential velocities;
       ! extrapolated normal velocity
       bbottom%fixed(:,:,Physics%DENSITY)   = .TRUE.
       bbottom%fixed(:,:,Physics%XVELOCITY) = .TRUE.
       bbottom%fixed(:,:,Physics%YVELOCITY) = .TRUE.
       bbottom%fixed(:,:,Physics%ZVELOCITY) = .FALSE.
       bbottom%fixed(:,:,Physics%PRESSURE)  = .TRUE.
    END SELECT
    ! outflow
    SELECT TYPE(btop => Timedisc%Boundary%boundary(TOP)%p)
    CLASS IS (boundary_fixed)
       btop%data(:,:,:,Physics%DENSITY)    = RHO0
       btop%data(:,:,:,Physics%XVELOCITY)  = 0.0
       btop%data(:,:,:,Physics%YVELOCITY)  = 0.0
       btop%data(:,:,:,Physics%ZVELOCITY)  = 0.0
       btop%data(:,:,:,Physics%PRESSURE)   = POUT
       ! imposed pressure;
       ! extrapolated density, tangential and normal velocities
       btop%fixed(:,:,Physics%DENSITY)   = .FALSE.
       btop%fixed(:,:,Physics%XVELOCITY) = .FALSE.
       btop%fixed(:,:,Physics%YVELOCITY) = .FALSE.
       btop%fixed(:,:,Physics%ZVELOCITY) = .FALSE.
       btop%fixed(:,:,Physics%PRESSURE)  = .TRUE.
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "tube with pressure gradient")

  END SUBROUTINE InitData

END PROGRAM poiseuille
