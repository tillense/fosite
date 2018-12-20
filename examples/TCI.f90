!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: TCI.f90                                                           #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
!> Program and data initialization for test of Taylor-Couette experiment
!! References:
!! [1] Taylor, G.I. (1923). "Stability of a Viscous Liquid contained between
!!     Two Rotating Cylinders". Phil. Trans. Royal Society A223: 289–343.
!!     DOI: 10.1098/rsta.1923.0008
!! [2] Pfister, G. and Schmidt, H. and Cliffe, K.~A. and Mulllin, T. (1988)
!!     "Bifurcation phenomena in Taylor–Couette flow in a very short annulus"
!!     J. Fluid Mech., vol. 191, p. 1-18
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
  REAL, PARAMETER    :: TSIM       = 1.0E-2    ! simulation time [TVIS]
  REAL, PARAMETER    :: GAMMA      = 1.4       ! ratio of specific heats
  REAL, PARAMETER    :: RHO0       = 1.0E+0    ! density 
  REAL, PARAMETER    :: P0         = 1.0E+3    ! pressure (high for incompress)
  REAL, PARAMETER    :: TAYLOR     = 2.0E+3    ! Taylor number (at RMIN)
  REAL, PARAMETER    :: ETA        = 1.0E-2    ! dynamic viscosity (constant)
  REAL, PARAMETER    :: MUHAT      = 0.0       ! = Omega_OUT / Omega_IN
  ! mesh settings
  INTEGER, PARAMETER :: RRES       = 50        ! resolution in r-direction
  INTEGER, PARAMETER :: ZRES       = 100       ! resolution in z-direction
  REAL, PARAMETER    :: RMIN       = 0.5       ! inner radius
  REAL, PARAMETER    :: RMAX       = 1.5       ! outer radius
  REAL, PARAMETER    :: HEIGHT     = 2.0       ! height of the cylinder
  ! output parameter
  INTEGER, PARAMETER :: ONUM       = 100       ! number of output time steps
  CHARACTER(LEN=64), PARAMETER &               ! output data dir
                     :: ODIR = "./"
  CHARACTER(LEN=64), PARAMETER &               ! output file name
                     :: OFNAME = "TCI"
  ! derived parameters
  REAL               :: TVIS                   ! viscous timescale
  REAL               :: OMEGA_IN,OMEGA_OUT     ! inner/outer angular velocity
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Fluxes, Sim%Physics, Sim%Timedisc)
  
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
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               timedisc, fluxes, vis
    REAL              :: x1,x2,y1,y2
    INTEGER           :: testnum
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! compute some simulation parameters
    TVIS      = RHO0 / ETA * RMAX**2
    OMEGA_IN  = TAYLOR * ETA / RHO0 / SQRT(RMIN*(RMAX-RMIN)**3)
    OMEGA_OUT = MUHAT * OMEGA_IN

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CYLINDRICAL, &
               "inum" / ZRES, &                  ! resolution in x and            !
               "jnum" / RRES, &                  !   y direction                  !             
               "xmin" / (-0.5*HEIGHT), &
               "xmax" / (0.5*HEIGHT), &
               "ymin" / RMIN, &
               "ymax" / RMAX)
    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA, &                 ! ratio of specific heats        !
              "dpmax"   / 1.0E+3)                 ! for advanced time step control !

   ! flux calculation and reconstruction method
   fluxes => Dict("order"     / LINEAR, &
            "fluxtype"  / KT, &
            "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
            "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
            "theta"     / 1.2)          ! optional parameter for limiter !

    ! boundary conditions
    boundary => Dict("western"  / PERIODIC, &
               "eastern"  / PERIODIC, &
               "southern" / NOSLIP, &
               "northern" / NOSLIP)

    ! viscosity source term
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / ETA, &
          "bulkconst" / (-2.0/3.0*ETA), &
          "cvis"      / 0.5)

    NULLIFY(sources)
    IF (ETA.GT.TINY(ETA)) &
        CALL SetAttr(sources, "vis", vis)

    ! time discretization settings
    timedisc => Dict(&
         "method"    / MODIFIED_EULER, &
         "order"     / 3, &
         "cfl"       / 0.4, &
         "stoptime"  / (TSIM*TVIS), &
         "dtlimit"   / (1.0E-12*TVIS), &
         "maxiter"   / 10000000)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,, &
!!$         fileformat = BINARY,, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // "log",, &
!!$         dtwall     = 1800,, &
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

  SUBROUTINE InitData(Mesh,Fluxes,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=32) :: info_str   
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.0
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.0
    Timedisc%pvar(:,:,Physics%PRESSURE)  = P0

    ! boundary data
    IF (GetType(Timedisc%Boundary(SOUTH)).EQ.NOSLIP) THEN
       Timedisc%Boundary(SOUTH)%data(:,:,Physics%ZVELOCITY) = RMIN * OMEGA_IN
    ENDIF
    IF (GetType(Timedisc%Boundary(NORTH)).EQ.NOSLIP) THEN
       Timedisc%Boundary(NORTH)%data(:,1,Physics%ZVELOCITY) = RMAX * OMEGA_OUT
    ENDIF

    ! add velocity perturbations
    CALL RANDOM_SEED
    CALL RANDOM_NUMBER(dv)
    Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
         + (dv(:,:,1)-0.5)*0.01
    Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
         + (dv(:,:,2)-0.5)*0.01

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: constant density and pressure")
    WRITE (info_str, '(ES8.2)') TVIS
    CALL Info(Mesh,"            viscous time scale " // TRIM(info_str) // " s")

  END SUBROUTINE InitData

END PROGRAM Init
