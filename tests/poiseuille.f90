!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: poiseuille.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006-2023                                                   #
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
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER :: TSIM  = 2.0E+0             ! simulation time [TVIS]     !
  REAL, PARAMETER :: RE    = 1.0E+0             ! Reynolds number            !
  REAL, PARAMETER :: MA    = 1.0E-2             ! Mach number << 1           !
  REAL, PARAMETER :: RHO0  = 1.0E-0             ! initial density            !
  REAL, PARAMETER :: ETA   = 1.0E-0             ! dynamic viscosity (const.) !
  REAL, PARAMETER :: GAM   = 1.4                ! ratio of spec. heat coeff. !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!  INTEGER, PARAMETER :: MGEO = CARTESIAN
  INTEGER, PARAMETER :: XRES = 30               ! radial resolution          !
  INTEGER, PARAMETER :: YRES = 1                ! azimuthal resolution       !
  INTEGER, PARAMETER :: ZRES = 5                ! resolution along the tube  !
  REAL, PARAMETER    :: LTUBE = 10.0            ! length of the tube         !
  REAL, PARAMETER    :: RTUBE = 1.0             ! radius of the tube         !
  !--------------------------------------------------------------------------!
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10               ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &               ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &               ! output data file name
                     :: OFNAME = 'poiseuille'
  ! derived parameters
  REAL               :: TVIS                    ! viscous timescale          !
  REAL               :: UMAX                    ! max. velocity of lam. flow !
  REAL               :: PIN,POUT                ! inlet/outlet pressure      !
  !--------------------------------------------------------------------------!
  CLASS(fosite),ALLOCATABLE   :: Sim
  LOGICAL :: ok
  !--------------------------------------------------------------------------!

TAP_PLAN(1)

  ALLOCATE(SIM)
  CALL SIM%InitFosite()

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL Sim%Setup()

  ! set initial condition
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
    CLASS(Fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, sources, &
                               vis, timedisc, fluxes
    REAL              :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! compute derived parameters
    ! set velocity scale according to the given Reynolds number
    UMAX = ETA/RHO0 * RE/RTUBE
    ! viscous time scale
    TVIS = RHO0/ETA  * RTUBE**2
    ! inlet pressure determined by Mach number
    PIN  = RHO0/GAM * (UMAX/MA)**2
    ! outlet pressure (for given length and maximum velocity of laminar solution
    POUT = PIN * (1.0 - 4*GAM * (MA**2/RE) * (LTUBE/RTUBE) )
    IF (POUT .LT. 0.0) &
      CALL Sim%Error("poiseuille::MakeConfig","negative outlet pressure")

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
!      bc(EAST)  = NOSLIP
      bc(EAST)  = CUSTOM
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
!            "decomposition" / (/1,-1,1/), &
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
              "gamma"   / GAM)           ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / VANLEER, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! viscosity source term
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / ETA, &
          "bulkconst" / (-2./3.*ETA), &
          "cvis" / 0.4) ! viscous Courant number

    sources => Dict("vis" / vis)

    ! time discretization settings
    timedisc => Dict("method"   / SSPRK, &
               "order"    / 3, &
               "cfl"      / 0.4, &
               "stoptime" / (TSIM*TVIS), &
               "dtlimit"  / (1.0E-7*TVIS), &
               "maxiter"  / 100000000)

    ! initialize data input/output
    datafile => Dict( &
               "fileformat" / VTK, &
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
    TYPE(marray_base) :: dv
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER           :: i,n,clock
    CHARACTER(LEN=32) :: info_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initialize velocity perturbations
    dv = marray_base(3)
    CALL RANDOM_SEED(size=n)
    ALLOCATE(seed(n))
    ! seed the rng with a mix from current time and mpi rank
    CALL SYSTEM_CLOCK(count=clock)
    seed = clock + (Timedisc%GetRank()+1) * (/(i-1, i=1,n)/)
    CALL RANDOM_NUMBER(dv%data1d)
    DEALLOCATE(seed)

    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    CLASS IS(statevector_euler)
      pvar%density%data1d(:)      = RHO0
      pvar%pressure%data3d(:,:,:) = PIN + (POUT-PIN)/LTUBE * Mesh%bccart(:,:,:,3)
      pvar%velocity%data1d(:)     = (dv%data1d(:)-0.5) * 0.01
    CLASS DEFAULT
      CALL Timedisc%Error("poiseuille::InitData","physics not supported")
    END SELECT

    ! fixed boundary conditions
    ! inflow
    SELECT TYPE(bbottom => Timedisc%Boundary%boundary(BOTTOM)%p)
    CLASS IS (boundary_fixed)
      bbottom%data(:,:,:,:)                = 0.0
      bbottom%data(:,:,:,Physics%DENSITY)  = RHO0
      bbottom%data(:,:,:,Physics%PRESSURE) = PIN
      ! imposed density, pressure and tangential velocities;
      ! extrapolated normal velocity
      bbottom%fixed(:,:,:)                 = .TRUE.
      bbottom%fixed(:,:,Physics%ZVELOCITY) = .FALSE.
    END SELECT
    ! outflow
    SELECT TYPE(btop => Timedisc%Boundary%boundary(TOP)%p)
    CLASS IS (boundary_fixed)
      btop%data(:,:,:,:)                   = 0.0
      btop%data(:,:,:,Physics%DENSITY)     = RHO0
      btop%data(:,:,:,Physics%PRESSURE)    = POUT
      ! imposed pressure;
      ! extrapolated density, tangential and normal velocities
      btop%fixed(:,:,:)                    = .FALSE.
      btop%fixed(:,:,Physics%PRESSURE)     = .TRUE.
    END SELECT
    ! noslip boundary at outer radius
    SELECT TYPE(b => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS(boundary_noslip)
      b%data(:,:,:,Physics%ZVELOCITY)      = 0.0
    CLASS IS(boundary_custom)
      CALL b%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_REFLECT,CUSTOM_REFLNEG,CUSTOM_FIXED,CUSTOM_FIXED,CUSTOM_REFLECT/))
      b%data(:,:,:,:)                      = 0.0
      b%data(:,:,:,Physics%DENSITY)        = RHO0 ! actually not used
      b%data(:,:,:,Physics%PRESSURE)       = POUT ! actually not used
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "tube with pressure gradient")
    WRITE (info_str, '(ES10.2)') RE
    CALL Physics%Info("              Reynolds number: " // TRIM(info_str))
    WRITE (info_str, '(ES10.2)') MA
    CALL Physics%Info("                  Mach number: " // TRIM(info_str))
    WRITE (info_str, '(ES10.2)') UMAX
    CALL Physics%Info("        max. laminar velocity: " // TRIM(info_str) // " m/s")
    WRITE (info_str, '(ES10.2)') PIN
    CALL Physics%Info("               inlet pressure: " // TRIM(info_str) // " Pa")
    WRITE (info_str, '(ES10.2)') 1.-POUT/PIN
    CALL Physics%Info("       rel. pressure gradient: " // TRIM(info_str))
    WRITE (info_str, '(ES10.2)') TVIS
    CALL Physics%Info("           viscous time scale: " // TRIM(info_str) // " s")

  END SUBROUTINE InitData

END PROGRAM poiseuille
