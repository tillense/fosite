!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: poiseuille.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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

!----------------------------------------------------------------------------!
!> Program and data initialization for a test of Hagen-Poiseuille equ. in a tube
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
  REAL, PARAMETER :: TSIM  = 30.0                ! simulation time            !
  REAL, PARAMETER :: RE    = 4.375              ! Reynolds Number            !
  REAL, PARAMETER :: PIN   = 24.0               ! inflow pressure            !
  REAL, PARAMETER :: POUT  = 23.0               ! outflow pressure           !
  REAL, PARAMETER :: RHO0  = 1.0                ! initial density            !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!  INTEGER, PARAMETER :: MGEO = CARTESIAN
  INTEGER, PARAMETER :: XRES = 40               ! radial resolution          !
  INTEGER, PARAMETER :: YRES = 80               ! resolution along the tube  !
  REAL, PARAMETER    :: LTUBE = 10.0            ! length of the tube         !
  REAL, PARAMETER    :: RTUBE = 1.0             ! radius of the tube         !
  !--------------------------------------------------------------------------!
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 30           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'poiseuille' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
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
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               vis, timedisc, fluxes
    REAL              :: x1,x2,y1,y2
    REAL              :: dvis, bvis
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! geometry dependent setttings
    SELECT CASE(MGEO)
     CASE(CYLINDRICAL)
       x1 = 0.0
       x2 = LTUBE
       y1 = 0.0
       y2 = RTUBE
    CASE(CARTESIAN)
       x1 = 0.0
       x2 = LTUBE
       y1 = 0.0
       y2 = RTUBE
    END SELECT

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2)
    
    ! compute viscosity constants
    ! dynamic viscosity
    dvis = SQRT(0.25 * (PIN-POUT) * RTUBE**3 * RHO0 / LTUBE / RE)
    ! bulk viscosity
    bvis = -2./3. * dvis

    ! geometry dependent setttings
    SELECT CASE(MGEO)
     CASE(CYLINDRICAL)
       bc(WEST)  = FIXED
       bc(EAST)  = FIXED
       bc(SOUTH) = AXIS
       bc(NORTH) = NOSLIP
    CASE(CARTESIAN)
       bc(WEST)  = FIXED
       bc(EAST)  = FIXED
       bc(SOUTH) = NOSLIP
       bc(NORTH) = NOSLIP
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / 1.4)           ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! viscosity source term
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

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY,, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log',, &
!!$         filecycles = 1)

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition    
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0. 
    Timedisc%pvar(:,:,Physics%PRESSURE) = PIN + (POUT-PIN)/LTUBE &
         * Mesh%bccart(:,:,2)

    ! fixed boundary conditions 
    ! inflow
    IF (GetType(Timedisc%Boundary(WEST)).EQ.FIXED) THEN
       Timedisc%Boundary(WEST)%data(:,:,Physics%DENSITY)    = RHO0
       Timedisc%Boundary(WEST)%data(:,:,Physics%XVELOCITY)  = 0.0
       Timedisc%Boundary(WEST)%data(:,:,Physics%YVELOCITY)  = 0.0
       Timedisc%Boundary(WEST)%data(:,:,Physics%ZVELOCITY)  = 0.0
       Timedisc%Boundary(WEST)%data(:,:,Physics%PRESSURE)   = PIN
       ! imposed density, pressure and tangential velocities;
       ! extrapolated normal velocity
       Timedisc%Boundary(WEST)%fixed(:,Physics%DENSITY)   = .TRUE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%XVELOCITY) = .FALSE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%YVELOCITY) = .TRUE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%ZVELOCITY) = .TRUE.
       Timedisc%Boundary(WEST)%fixed(:,Physics%PRESSURE)  = .TRUE.
    ENDIF
    ! outflow
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FIXED) THEN
       Timedisc%Boundary(EAST)%data(:,:,Physics%DENSITY)    = RHO0
       Timedisc%Boundary(EAST)%data(:,:,Physics%XVELOCITY)  = 0.0
       Timedisc%Boundary(EAST)%data(:,:,Physics%YVELOCITY)  = 0.0
       Timedisc%Boundary(EAST)%data(:,:,Physics%ZVELOCITY)  = 0.0
       Timedisc%Boundary(EAST)%data(:,:,Physics%PRESSURE)   = POUT
       ! imposed pressure; 
       ! extrapolated density, tangential and normal velocities
       Timedisc%Boundary(EAST)%fixed(:,Physics%DENSITY)   = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%XVELOCITY) = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%YVELOCITY) = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%ZVELOCITY) = .FALSE.
       Timedisc%Boundary(EAST)%fixed(:,Physics%PRESSURE)  = .TRUE.
    ENDIF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "tube with pressure gradient")

  END SUBROUTINE InitData

END PROGRAM Init
