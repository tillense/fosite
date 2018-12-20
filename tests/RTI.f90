!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: RTI.f90                                                           #
!#                                                                           #
!# Copyright (C) 2008-2012                                                   # 
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \test Rayleight-Taylor instability
!! \author Björn Sperling
!! \author Tobias Illenseer
!!
!! References:
!! [1] Rayleigh, Lord (1883) "Investigation of the character of the equilibrium
!!     of an incompressible heavy fluid of variable density"
!!     Proceedings of the London Mathematical Society 14: 170–177
!!     DOI: 10.1112/plms/s1-14.1.170
!! [2] Taylor, Sir Geoffrey Ingram (1950). "The instability of liquid surfaces
!!     when accelerated in a direction perpendicular to their planes"
!!     Proceedings of the Royal Society of London. Series A, (1065): 192–196
!!     DOI: 10.1098/rspa.1950.0052
!! [3] D.H. Sharp "An overview of Rayleigh-Taylor instability"
!!     Physica D: Nonlinear Phenomena, vol. 12, Issues 1-3, July 1984, Pages 3-10
!!     DOI: 10.1016/0167-2789(84)90510-4
!----------------------------------------------------------------------------!
PROGRAM RTI
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
  USE common_dict
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL,    PARAMETER  :: TSIM    = 10.0     ! simulation time
  REAL,    PARAMETER  :: DYNVIS  = 0.0      ! dynamic viscosity constant
  REAL,    PARAMETER  :: BULKVIS = 0.0      ! bulk viscosity constant
!!$  REAL,    PARAMETER  :: DYNVIS  = 1.0E-4   
!!$  REAL,    PARAMETER  :: BULKVIS = -6.67E-5
  ! initial condition (SI units)
  REAL,    PARAMETER  :: RHO0    = 2.0      ! density: upper region
  REAL,    PARAMETER  :: RHO1    = 1.0      ! density: lower region
  REAL,    PARAMETER  :: YACC    = 0.2      ! grav. acceleration
  REAL,    PARAMETER  :: P0      = 1.2      ! pressure at the top
  REAL,    PARAMETER  :: A0      = 0.02     ! amplitude of initial disturbance
  ! mesh settings
  INTEGER, PARAMETER  :: XRES    = 50       ! resolution in x
  INTEGER, PARAMETER  :: YRES    = 100      ! resolution in y
  REAL, PARAMETER     :: WIDTH   = 1.0      ! width of comp. domain
  REAL, PARAMETER     :: HEIGHT  = 2.0      ! height of comp. domain
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &           ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &           ! output data file name
                     :: OFNAME = 'RTI' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)    :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
  
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
    INTEGER           :: bc(4),sgbc
    DOUBLE PRECISION, PARAMETER :: GN = 6.6742E-11     ! [m^3/kg/s^2]
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               timedisc, fluxes, caccel, vis
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
               "inum" / XRES, &          ! resolution in x and            !
               "jnum" / YRES, &          !   y direction                  !             
               "xmin" / 0., &
               "xmax" / WIDTH, &
               "ymin" / 0., &
               "ymax" / HEIGHT)

    ! physics settings
    physics => Dict("problem" / EULER2D, &
              "gamma"   / 1.4)           ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! boundary conditions
    boundary => Dict("western"  / REFLECTING, &
               "eastern"  / REFLECTING, &
               "southern" / REFLECTING, &
               "northern" / REFLECTING)

    caccel => Dict("stype"  / C_ACCEL, & 
             "yaccel" / (-YACC))           ! acceleration in y-direction


    ! viscosity source term
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / DYNVIS, &
          "bulkconst" / BULKVIS)
    
    sources => Dict("caccel" / caccel)

    IF ((DYNVIS.GT.TINY(DYNVIS)).OR.(BULKVIS.GT.TINY(BULKVIS))) &
        CALL SetAttr(sources, "vis", vis)

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.4, &
           "stoptime" / TSIM, &
           "dtlimit"  / 1.0E-4, &
           "maxiter"  / 100000)

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


  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)&
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX):: y0
    REAL, DIMENSION(2) :: accel
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh, Physics
    INTENT(OUT)       :: pvar, cvar
    !------------------------------------------------------------------------!
    ! this marks the line between the two fluids
    y0(:,:) = 0.5*Mesh%ymax + A0*COS(2*PI*Mesh%bcenter(:,:,1)/Mesh%xmax)

    ! initial hydrostatic stratification
    WHERE (Mesh%bcenter(:,:,2).GT.y0(:,:))
       ! upper fluid
       pvar(:,:,Physics%DENSITY)  = RHO0
       pvar(:,:,Physics%PRESSURE) = P0 + YACC * RHO0 * (Mesh%ymax-Mesh%bcenter(:,:,2))
    ELSEWHERE
       ! lower fluid
       pvar(:,:,Physics%DENSITY)  = RHO1
       pvar(:,:,Physics%PRESSURE) = P0 + YACC * (RHO1 * (y0(:,:)-Mesh%bcenter(:,:,2)) &
            + RHO0 * (Mesh%ymax-y0(:,:)))
    END WHERE

    ! velocity vanishes everywhere
    pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = 0.

    CALL Convert2Conservative(Physics,Mesh,pvar,cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Rayleigh–Taylor instability")

  END SUBROUTINE InitData
END PROGRAM RTI
