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
!!
!! \example poiseuille.f90
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
    CASE(CARTESIAN)
       x1 = 0.0
       x2 = RTUBE
       y1 = 0.0
       y2 = RTUBE
       z1 = 0.0
       z2 = LTUBE
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

    ! compute viscosity constants
    ! dynamic viscosity
    dvis = SQRT(0.25 * (PIN-POUT) * RTUBE**3 * RHO0 / LTUBE / RE)
    ! bulk viscosity
    bvis = -2./3. * dvis

    ! geometry dependent setttings
    SELECT CASE(MGEO)
     CASE(CYLINDRICAL)
       bc(WEST)  = AXIS
       bc(EAST)  = NOSLIP
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
       bc(BOTTOM)= FIXED
       bc(TOP)   = FIXED
    CASE(CARTESIAN)
       bc(WEST)  = NOSLIP
       bc(EAST)  = NOSLIP
       bc(SOUTH) = NOSLIP
       bc(NORTH) = NOSLIP
       bc(BOTTOM)= FIXED
       bc(TOP)   = FIXED
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH), &
               "bottomer" / bc(BOTTOM), &
               "topper"   / bc(TOP))

    ! physics settings
    physics => Dict("problem" / EULER3D, &
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
    datafile => Dict( &
               "fileformat" / BINARY, &
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
    CLASS(Physics_base) :: Physics
    CLASS(Mesh_base)    :: Mesh
    CLASS(Timedisc_base):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,3) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
    Timedisc%pvar(:,:,:,Physics%DENSITY)   = RHO0
    Timedisc%pvar(:,:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,:,Physics%ZVELOCITY) = 0.
    Timedisc%pvar(:,:,:,Physics%PRESSURE) = PIN + (POUT-PIN)/LTUBE &
         * Mesh%bccart(:,:,:,3)

    ! fixed boundary conditions
    ! inflow
    IF ((Timedisc%Boundary%Boundary(BOTTOM)%p%GetType()).EQ.FIXED) THEN
       Timedisc%Boundary%Boundary(BOTTOM)%p%data(:,:,:,Physics%DENSITY)    = RHO0
       Timedisc%Boundary%Boundary(BOTTOM)%p%data(:,:,:,Physics%XVELOCITY)  = 0.0
       Timedisc%Boundary%Boundary(BOTTOM)%p%data(:,:,:,Physics%YVELOCITY)  = 0.0
       Timedisc%Boundary%Boundary(BOTTOM)%p%data(:,:,:,Physics%ZVELOCITY)  = 0.0
       Timedisc%Boundary%Boundary(BOTTOM)%p%data(:,:,:,Physics%PRESSURE)   = PIN
       ! imposed density, pressure and tangential velocities;
       ! extrapolated normal velocity
       Timedisc%Boundary%Boundary(BOTTOM)%p%fixed(:,:,Physics%DENSITY)   = .TRUE.
       Timedisc%Boundary%Boundary(BOTTOM)%p%fixed(:,:,Physics%XVELOCITY) = .TRUE.
       Timedisc%Boundary%Boundary(BOTTOM)%p%fixed(:,:,Physics%YVELOCITY) = .TRUE.
       Timedisc%Boundary%Boundary(BOTTOM)%p%fixed(:,:,Physics%ZVELOCITY) = .FALSE.
       Timedisc%Boundary%Boundary(BOTTOM)%p%fixed(:,:,Physics%PRESSURE)  = .TRUE.
    END IF
    ! outflow
    IF ((Timedisc%Boundary%Boundary(TOP)%p%GetType()).EQ.FIXED) THEN
       Timedisc%Boundary%Boundary(TOP)%p%data(:,:,:,Physics%DENSITY)    = RHO0
       Timedisc%Boundary%Boundary(TOP)%p%data(:,:,:,Physics%XVELOCITY)  = 0.0
       Timedisc%Boundary%Boundary(TOP)%p%data(:,:,:,Physics%YVELOCITY)  = 0.0
       Timedisc%Boundary%Boundary(TOP)%p%data(:,:,:,Physics%ZVELOCITY)  = 0.0
       Timedisc%Boundary%Boundary(TOP)%p%data(:,:,:,Physics%PRESSURE)   = POUT
       ! imposed pressure;
       ! extrapolated density, tangential and normal velocities
       Timedisc%Boundary%Boundary(TOP)%p%fixed(:,:,Physics%DENSITY)   = .FALSE.
       Timedisc%Boundary%Boundary(TOP)%p%fixed(:,:,Physics%XVELOCITY) = .FALSE.
       Timedisc%Boundary%Boundary(TOP)%p%fixed(:,:,Physics%YVELOCITY) = .FALSE.
       Timedisc%Boundary%Boundary(TOP)%p%fixed(:,:,Physics%ZVELOCITY) = .FALSE.
       Timedisc%Boundary%Boundary(TOP)%p%fixed(:,:,Physics%PRESSURE)  = .TRUE.
    END IF

    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "tube with pressure gradient")

  END SUBROUTINE InitData

END PROGRAM poiseuille
