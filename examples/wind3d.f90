!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: wind3d.f90                                                        #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
!> 3D (with rotational symmmetry) AGN wind simulation
!----------------------------------------------------------------------------!

!**************************************!
!* IMPORTANT:                         *!
!* - compile with autodouble          *!
!* - use primitive reconstruction     *!
!**************************************!

PROGRAM Init
  USE fosite
  USE constants_common, ONLY : C,GN,KB,NA
  USE constants_generic
  USE physics_generic
  USE fluxes_generic
  USE geometry_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE timedisc_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! some constants
  REAL, PARAMETER :: MSUN = 1.989E+30        ! solar mass [kg]
  REAL, PARAMETER :: YEAR = 3.155693E+07     ! one year [s]
  REAL, PARAMETER :: RG   = KB * NA          ! universal gas constant
  !--------------------------------------------------------------------------!
  ! basic simulation parameters
  REAL, PARAMETER :: TSIM   = 1.0E-1*YEAR    ! simulation time in [s]
  REAL, PARAMETER :: MBH    = 1.0E+6 * MSUN  ! mass of central point mass
  REAL, PARAMETER :: MDOT   = 1.0e-2 * MSUN/YEAR ! accretion rate
  REAL, PARAMETER :: RS     = 2*GN*MBH/C**2  ! Schwarzschild radius
  REAL, PARAMETER :: R0DISK = 3*RS           ! inner radius of disk
  REAL, PARAMETER :: GAMMA  = 5./3.          ! ratio of specific heats
  REAL, PARAMETER :: MU     = 0.602E-03      ! mean molecular weight
  REAL, PARAMETER :: RHOINF = 1.0E-16        ! density at infinity
  REAL, PARAMETER :: TINF   = 1.0E+8         ! temperature at infinity
  REAL, PARAMETER :: PINF   = RG/MU*TINF*RHOINF ! pressure at infinity
  ! mesh geometry
!   INTEGER, PARAMETER :: MGEO = CYLINDRICAL
  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
  ! computational domain
  REAL, PARAMETER :: RMIN = 3.0E+0*RS        ! inner radius [m]
  REAL, PARAMETER :: RMAX = 1.0E+3*RS        ! outer radius [m]
  ! resolution
  INTEGER, PARAMETER :: XRES = 50
  INTEGER, PARAMETER :: YRES = 50
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 100           ! number of output time steps  !
  CHARACTER(LEN=256), PARAMETER :: ODIR &! output directory
                                 = "./"
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'wind3d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP) :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Timedisc, Sim%Mesh, Sim%Physics, Sim%Fluxes)
  
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
                               grav, timedisc, fluxes, thomson, cool
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    SELECT CASE(MGEO)
    CASE(CYLINDRICAL)
       ! mesh settings
       mesh => Dict("meshtype" / MIDPOINT, &
              "geometry" / MGEO, &
                  "inum" / XRES, &       ! resolution in z and            !
                  "jnum" / YRES, &       !   r direction                  !             
                  "xmin" / 0.0, &
                  "xmax" / RMAX, &
                  "ymin" / RMIN, &
                  "ymax" / RMAX)
    CASE(SPHERICAL)
       ! mesh settings
       mesh => Dict("meshtype" / MIDPOINT, &
              "geometry" / MGEO, &
                  "inum" / XRES, &       ! resolution in r and            !
                  "jnum" / YRES, &        !   theta direction              !             
                  "xmin" / RMIN, &
                  "xmax" / RMAX, &
                  "ymin" / 0.0, &
                  "ymax" / (0.5 * PI))
    CASE(OBLATE_SPHEROIDAL)
       ! mesh settings
       mesh => Dict("meshtype" / MIDPOINT, &
              "geometry" / MGEO, &
                  "inum" / XRES, &        ! resolution in r and            !
                  "jnum" / YRES, &        !   theta direction              !             
                  "xmin" / 0.0, &
                  "xmax" / 1.4, &
                  "ymin" / 0.0, &
                  "ymax" / (0.5 * PI), &
                "gparam" / RMAX)  ! optional geometry parameter    !
    CASE(TANCYLINDRICAL)
       ! mesh settings
       mesh => Dict("meshtype" / MIDPOINT, &
              "geometry" / MGEO, &
                  "inum" / XRES, &        ! resolution in r and            !
                  "jnum" / YRES, &        !   theta direction              !             
                  "xmin" / 0.0, &
                  "xmax" / ATAN(5.0), &   ! ZMAX = 5*RMAX
                  "ymin" / RMIN, &
                  "ymax" / RMAX, &
                "gparam" / RMAX)          ! geometry parameter: vertical scale
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram", "geometry not supported")
    END SELECT

    SELECT CASE(MGEO)
    CASE(CYLINDRICAL,TANCYLINDRICAL)
       ! boundary conditions
       boundary => Dict("western"  / NOSLIP, &
                  "eastern"  / ABSORBING, &
                  "southern" / AXIS, &
                  "northern" / ABSORBING)
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram", "geometry not supported")
    END SELECT

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA, &         ! ratio of specific heats        !
              "mu"      / MU)              ! mean molecular weight          !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MINMOD, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! source term due to a point mass
    grav => Dict("stype" / GRAVITY, &
            "pmass/gtype" / POINTMASS, &      ! grav. accel. of a point mass   !
            "pmass/mass"  / MBH, &            ! black hole mass [kg]           !
            "output/accel" / 1)

    ! account for energy losses due to radiative cooling
    thomson => Dict("stype" / DISK_THOMSON, &
              "mass"  / MBH, &                   ! mass of central black hole
              "mdot"  / MDOT, &                  ! accretion rate
              "rin"   / R0DISK, &                ! inner radius of disk
              "rout"  / (10 * RMAX), &           ! outer radius of disk
              "output/accel" / 1)

    ! account for energy losses due to radiative cooling
    cool => Dict("stype" / COOLING, &
           "cvis"  / 0.5)             ! CFL number for cooling time    !

    sources => Dict("grav"   / grav, &
               "thomson" / thomson, &
               "cool"    / cool)

    ! time discretization settings
    timedisc => Dict( &
!            "method"    / MODIFIED_EULER, "order" / 3, &
           "method"    / RK_FEHLBERG, "order" / 5, &
           "cfl"       / 0.4, &
           "stoptime"  / TSIM, &
           "tol_rel"   / 1.0e-03, &
           "tol_abs"   / (/ 0.0, 1.0E-10, 1.0E-10, 1.0E-10, 0.0/), &
           "dtlimit"   / (1.0e-6*TSIM), &
           "maxiter"   / 10000000)

    ! initialize data input/output
!     datafile => Dict("fileformat" / VTK, "unit" / 5555, "filecycles" / (ONUM+1), &
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


  SUBROUTINE InitData(Timedisc,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,dir
    REAL              :: a,r
    REAL              :: Ldisk,Ledd,cs2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc,Physics,Fluxes
    !------------------------------------------------------------------------!
    ! disk luminosity
    Ldisk = 0.5*MDOT*Physics%constants%GN*MBH / R0DISK  ! disk luminosity
    Ledd  = 4*PI*Physics%constants%GN*MBH*Physics%Constants%C/Physics%Constants%KE

    ! initial condition, adiabatic sphere in balance
    cs2 = GAMMA*PINF/RHOINF
    a = (GAMMA-1.0)*Physics%constants%GN*MBH/cs2*(1.-Ldisk/Ledd)

    Timedisc%pvar(:,:,Physics%DENSITY) = RHOINF*(1.+a/Mesh%radius%bcenter(:,:))**(1./(GAMMA-1.0))
    Timedisc%pvar(:,:,Physics%PRESSURE) = PINF*(Timedisc%pvar(:,:,Physics%DENSITY)/RHOINF)**GAMMA

    ! set balanced velocity field
    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%XVELOCITY+Physics%DIM) = &
           GetCentrifugalVelocity(Timedisc,Mesh,Physics,Fluxes,(/1.,0.,0./))

    ! boundary conditions at the disks surface
    SELECT CASE(GetType(Mesh%geometry))
    CASE(CYLINDRICAL,TANCYLINDRICAL)
       ! set the velocity due to the centrifugal force
       ! z=0 boundary (disk surface)
       IF (GetType(Timedisc%boundary(WEST)).EQ.NOSLIP) THEN
          ! no radial velocity 
          Timedisc%boundary(WEST)%data(:,:,Physics%YVELOCITY) = 0.0
          ! Keplerian rotation
          DO j=Mesh%JMIN,Mesh%JMAX
             DO i=1,Mesh%GNUM
                Timedisc%boundary(WEST)%data(i,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
                     = SQRT(Physics%constants%GN*MBH/Mesh%hz%bcenter(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX))
             END DO
          END DO
       END IF
       ! set farfield data to initial condition
       IF (GetType(Timedisc%Boundary(WEST)).EQ.FARFIELD) THEN
          DO j=Mesh%JMIN,Mesh%JMAX
             DO i=1,Mesh%GNUM
                Timedisc%Boundary(WEST)%data(i,j,:) = Timedisc%pvar(Mesh%IMIN-i,j,:)
                Timedisc%boundary(WEST)%data(i,j,Physics%ZVELOCITY) &
                     = SQRT(Physics%constants%GN*MBH/Mesh%hz%bcenter(Mesh%IMIN-i,j))
            END DO
          END DO
       END IF
    END SELECT

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: AGN wind")
  END SUBROUTINE InitData

END PROGRAM Init
