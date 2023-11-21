!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: RTI.f90                                                           #
!#                                                                           #
!# Copyright (C) 2008-2023                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> \test Rayleight-Taylor instability
!! \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \todo implement 3D test in cylindrical geometry with
!!       and without axial symmetry
!!
!! References:
!! - \cite strutt1883 Rayleigh, Lord (1883) "Investigation of the character
!!     of the equilibrium of an incompressible heavy fluid of variable density"
!!     Proceedings of the London Mathematical Society 14: 170–177
!!     DOI: 10.1112/plms/s1-14.1.170
!! - \cite taylor1950  Taylor, Sir Geoffrey Ingram (1950). "The instability of liquid surfaces
!!     when accelerated in a direction perpendicular to their planes"
!!     Proceedings of the Royal Society of London. Series A, (1065): 192–196
!!     DOI: 10.1098/rspa.1950.0052
!! - \cite sharp1984 D.H. Sharp "An overview of Rayleigh-Taylor instability"
!!     Physica D: Nonlinear Phenomena, vol. 12, Issues 1-3, July 1984, Pages 3-10
!!     DOI: 10.1016/0167-2789(84)90510-4
!----------------------------------------------------------------------------!
PROGRAM RTI
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL,    PARAMETER  :: TSIM    = 10.0     ! simulation time                !
  REAL,    PARAMETER  :: DYNVIS  = 0.0      ! dynamic viscosity constant     !
  REAL,    PARAMETER  :: BULKVIS = 0.0      ! bulk viscosity constant        !
!  REAL,    PARAMETER  :: DYNVIS  = 1.0E-4
!  REAL,    PARAMETER  :: BULKVIS = -6.67E-5
  ! initial condition (SI units)
  REAL,    PARAMETER  :: RHO0    = 2.0      ! density: upper region          !
  REAL,    PARAMETER  :: RHO1    = 1.0      ! density: lower region          !
  REAL,    PARAMETER  :: YACC    = 0.2      ! grav. acceleration             !
  REAL,    PARAMETER  :: P0      = 1.2      ! pressure at the top            !
  REAL,    PARAMETER  :: A0      = 0.02     ! amplitude of init. disturbance !
  ! mesh settings
  INTEGER, PARAMETER  :: XRES    = 50       ! resolution in x                !
  INTEGER, PARAMETER  :: YRES    = 100      ! resolution in y                !
  INTEGER, PARAMETER  :: ZRES    = 1        ! z-resolution
  REAL, PARAMETER     :: WIDTH   = 1.0      ! width of comp. domain          !
  REAL, PARAMETER     :: HEIGHT  = 2.0      ! height of comp. domaina        !
  ! output parameters
  INTEGER, PARAMETER  :: ONUM    = 10       ! number of output data sets     !
  CHARACTER(LEN=256), PARAMETER &
                      :: ODIR    = './'     ! output data dir                !
  CHARACTER(LEN=256), PARAMETER &
                      :: OFNAME  = 'RTI'    ! output data file name          !
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  !--------------------------------------------------------------------------!

!  TAP_PLAN(1)

  ALLOCATE(Sim)

  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
  CALL Sim%Run()
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

!  TAP_CHECK(.TRUE.,"Simulation finished")
!  TAP_DONE

  CONTAINS

  !> Set all configurations and safe it in a dictionary
  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)            :: Sim
    TYPE(Dict_TYP), POINTER  :: config
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: mesh, physics, boundary, datafile, &
                                sources, timedisc, fluxes, caccel, vis
    REAL                     :: zmax
    !------------------------------------------------------------------------!
    INTENT(INOUT)            :: Sim
    !------------------------------------------------------------------------!
    IF (ZRES.GT.1) THEN
      zmax = 0.5*WIDTH
    ELSE
      zmax = 0.0
    END IF
    ! mesh settings
    mesh => Dict( &
              "meshtype"      / MIDPOINT,   &    ! use midpoint rule         !
              "geometry"      / CARTESIAN,  &    ! cartesian grid            !
              "inum"          / XRES,       &    ! resolution in x-direction !
              "jnum"          / YRES,       &    ! resolution in y-direction !
              "knum"          / ZRES,       &    ! resolution in z-direction !
              "xmin"          / (-0.5*WIDTH), &  ! minimum value in x-dir.   !
              "xmax"          / (0.5*WIDTH), &   ! maximum value in x-dir.   !
              "ymin"          / 0.,         &    ! minimum value in y-dir.   !
              "ymax"          / HEIGHT,     &    ! maximum value in y-dir.   !
              "zmin"          / (-zmax),    &    ! minimum value in z-dir.   !
              "zmax"          / zmax)            ! maximum value in z-dir.   !

    ! physics settings
    physics => Dict( &
              "problem"       / EULER, &         ! standard hydrodynamics    !
              "gamma"         / 1.4)             ! ratio of specific heats   !

    ! flux calculation and reconstruction method
    fluxes => Dict( &
              "fluxtype"      / KT, &            ! use Kurganov-Tadmor flux  !
              "order"         / LINEAR, &        ! linear reconstruction     !
              "variables"     / CONSERVATIVE, &  ! vars. for reconstruction  !
              "limiter"       / MONOCENT, &      ! type of the limiter       !
              "theta"         / 1.2)             ! optional param. for MC    !

    ! boundary conditions
    boundary => Dict( &
              "western"       / REFLECTING, &    ! reflecting boundary cond. !
              "eastern"       / REFLECTING, &    ! reflecting boundary cond. !
              "southern"      / REFLECTING, &    ! reflecting boundary cond. !
              "northern"      / REFLECTING, &    ! reflecting boundary cond. !
              "bottomer"      / REFLECTING, &    ! reflecting boundary cond. !
              "topper"        / REFLECTING)      ! reflecting boundary cond. !

    ! c-acceleration term
    caccel => Dict( &
              "stype"         / C_ACCEL, &       ! source 1: acceleration    !
              "yaccel"        / (-YACC))         ! constant acc.in y-dir.    !


    ! viscosity source term
    vis => Dict( &
              "stype"         / VISCOSITY, &     ! viscosity  source         !
              "vismodel"      / MOLECULAR, &     ! visc. model: molecular    !
              "dynconst"      / DYNVIS, &        ! const. dynamic viscosity  !
              "bulkconst"     / BULKVIS)         ! const. bulk viscosity     !

    ! collect sources in dictionary
    sources => Dict( &
!              "vis"           / vis, &           ! incl. visc. from above    !
              "caccel"        / caccel)          ! incl. accel. from above   !

    IF ((DYNVIS.GT.TINY(DYNVIS)).OR.(BULKVIS.GT.TINY(BULKVIS))) &
        CALL SetAttr(sources, "vis", vis)

    ! time discretization settings
    timedisc => Dict( &
              "method"        / MODIFIED_EULER, &! use modified euler        !
              "order"         / 3, &             ! third order accuracy      !
              "cfl"           / 0.4, &           ! courant-number            !
              "stoptime"      / TSIM, &          ! simulation stop-time      !
              "dtlimit"       / 1.0E-4, &        ! smallest allowed timestep !
              "maxiter"       / 100000)          ! max. iters before abort   !

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"    / VTK, &
              "filename"      / (TRIM(ODIR) // TRIM(OFNAME)), &
              "count"         / ONUM)

    ! collect all above dicts in the configuration dict
    config => Dict( &
              "mesh"          / mesh, &          ! all mesh-settings         !
              "physics"       / physics, &       ! all physics settings      !
              "boundary"      / boundary, &      ! all bounary settings      !
              "fluxes"        / fluxes, &        ! all fluxes settings       !
              "sources"       / sources, &       ! all sources               !
              "timedisc"      / timedisc, &      ! all timedisc settings     !
              "datafile"      / datafile)        ! all input/output settings !
  END SUBROUTINE MakeConfig


  !> \public set initial conditions
  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(marray_compound), POINTER, INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                      :: y0
    !------------------------------------------------------------------------!
    ! y-coordiante of the contact surface between the two fluids:
    ! half of the height with a little dip on the axis (x=0,y=0)
    y0(:,:,:) = 0.5*Mesh%ymax - A0*EXP(-0.5/(0.4*WIDTH)**2 &
         * (Mesh%bcenter(:,:,:,1)**2+Mesh%bcenter(:,:,:,3)**2))

    ! initial hydrostatic stratification
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler)
      WHERE (Mesh%bcenter(:,:,:,2).GT.y0(:,:,:))
        ! upper fluid
        p%density%data3d(:,:,:)  = RHO0
        p%pressure%data3d(:,:,:) = P0 + YACC * RHO0 * (Mesh%ymax-Mesh%bcenter(:,:,:,2))
      ELSEWHERE
        ! lower fluid
        p%density%data3d(:,:,:)  = RHO1
        p%pressure%data3d(:,:,:) = P0 + YACC * (RHO1 * (y0(:,:,:)-Mesh%bcenter(:,:,:,2)) &
              + RHO0 * (Mesh%ymax-y0(:,:,:)))
      END WHERE
      ! velocity vanishes everywhere
      p%velocity%data1d(:) = 0.
    CLASS DEFAULT
      CALL Physics%Error("RTI::InitData","only non-isothermal HD supported")
    END SELECT

    CALL Physics%Convert2Conservative(pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Rayleigh–Taylor instability")
  END SUBROUTINE InitData
END PROGRAM RTI
