!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gauss2d.f90                                                       #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> 2D Gaussian pressure pulse
!! \author Tobias Illenseer
!----------------------------------------------------------------------------!
PROGRAM gauss2d
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
  ! simulation parameter
  REAL, PARAMETER :: TSIM   = 0.3         ! simulation time
  REAL, PARAMETER :: GAMMA  = 1.4         ! ratio of specific heats
  REAL, PARAMETER :: CSISO  = 0.0         ! if .ne. 0.0 -> isothermal simulation
                                          !   with CSISO as sound speed
  ! initial condition (dimensionless units)
  REAL, PARAMETER :: RHO0   = 1.0         ! ambient density
  REAL, PARAMETER :: P0     = 1.0         ! ambient pressure
  REAL, PARAMETER :: AMP    = 1.0         ! amplitude of the pulse
  REAL, PARAMETER :: PWIDTH = 0.06        ! half width of the Gaussian
  REAL, PARAMETER :: ETA    = 0.0         ! dynamic viscosity (0.0 disables)
  REAL, PARAMETER :: OMEGA  = 0.0         ! angular velocity (0.0 disables)
  REAL, PARAMETER :: X0     = 0.0         ! center of pressure pulse given 
  REAL, PARAMETER :: Y0     = 0.0         ! in cartesian coordinates (x,y)
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN  ! geometry of the mesh
!!$  INTEGER, PARAMETER :: MGEO = POLAR
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
!!$  INTEGER, PARAMETER :: MGEO = BIPOLAR
!!$  INTEGER, PARAMETER :: MGEO = ELLIPTIC
  INTEGER, PARAMETER :: XRES = 50          ! resolution
  INTEGER, PARAMETER :: YRES = 50 
  REAL, PARAMETER    :: RMIN = 1.0E-2      ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX = 0.5*1.4142  ! outer radius 0.5*SQRT(2)
  REAL, PARAMETER    :: GPAR = 0.5         ! geometry scaling parameter
  ! physics settings
!!$  LOGICAL, PARAMETER :: WITH_IAR = .TRUE.  ! use EULER2D_IAMROT
  LOGICAL, PARAMETER :: WITH_IAR = .FALSE.
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10        ! number of output time steps
  CHARACTER(LEN=256), PARAMETER :: ODIR &! output directory
                                 = "./"
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'gauss2d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
  TAP_DONE

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
    TYPE(Dict_TYP), POINTER :: physics, boundary, datafile, logfile, sources, &
                               timedisc, fluxes, rotframe, vis, mesh
    REAL,DIMENSION(Sim%Mesh%IMAX-Sim%Mesh%IMIN+1,&
                   Sim%Mesh%JMAX-Sim%Mesh%JMIN+1) :: dummy
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 = -0.5
       x2 = 0.5
       y1 = -0.5
       y2 = 0.5
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NOSLIP
       bc(EAST)  = NOSLIP
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(0.5*SQRT(2.0)/GPAR)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(0.5*SQRT(2.0)/GPAR)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       x1 = ASINH(RMIN/GPAR)
       x2 = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(BIPOLAR)
       x1 = 0.0
       x2 = 2*PI
       y1 = -ASINH(GPAR/RMIN - 1.0)
       y2 = -y1
       bc(WEST)  = PERIODIC
       bc(EAST)  = PERIODIC      
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
    CASE(ELLIPTIC)
       x1 = RMIN
       x2 = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = 2*PI
       bc(WEST)  = FOLDED
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram", &
            " geometry should be one of cartesian,polar,logpolar,tanpolar,sinhpolar or bipolar")
    END SELECT
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "omega"    / OMEGA, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "gparam"   / GPAR)

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    IF (CSISO.GT.TINY(CSISO)) THEN
       physics => Dict("problem" / EULER2D_ISOTHERM, &
                 "cs"      / CSISO)                      ! isothermal sound speed  !
    ELSE
       IF (WITH_IAR) THEN
          physics => Dict("problem"   / EULER2D_IAMROT, &
                    "gamma"     / GAMMA)                 ! ratio of specific heats !
       ELSE
          physics => Dict("problem"   / EULER2D, &
                    "gamma"     / GAMMA)
       END IF
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
!!$             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "variables" / PRIMITIVE, &        ! vars. to use for reconstruction!
             "limiter"   / OSPRE, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    NULLIFY(sources)

    ! viscosity source term
    IF (ETA.GT.TINY(ETA)) THEN
        vis => Dict("stype"    / VISCOSITY, &
                    "vismodel" / MOLECULAR, &
                    "dynconst" / ETA)
        CALL SetAttr(sources, "vis", vis)
    END IF

    ! inertial forces
    IF (.NOT.WITH_IAR .AND. (OMEGA.GT.TINY(OMEGA))) THEN
        rotframe => Dict("stype" / ROTATING_FRAME)
        CALL SetAttr(sources, "rotframe", rotframe)
    END IF

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.4, &
           "stoptime" / TSIM, &
           "tol_rel"  / 0.01, &
           "dtlimit"  / 1.0E-5, &
           "output/iangularmomentum" / 1, &
           "output/rothalpy" / 1, &
           "maxiter"  / 1000000)

    ! initialize data input/output
    datafile => Dict( &
         "fileformat" / GNUPLOT, "filecycles" / 0, &
         "decimals" / 6, "cartcoords" / 1, &
         "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
         "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)

    IF (ASSOCIATED(sources)) &
        CALL SetAttr(config, "sources", sources)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP)  :: Physics
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Timedisc_TYP) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Sources_TYP), POINTER :: sp
    INTEGER           :: i,j,dir,ig
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: radius
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: ephi,posvec
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics
    INTENT(INOUT)      :: Timedisc
    !------------------------------------------------------------------------!
    IF (ABS(X0).LE.TINY(X0).AND.ABS(Y0).LE.TINY(Y0)) THEN
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:) = Mesh%radius%bcenter(:,:)
       posvec(:,:,:) = Mesh%posvec%bcenter(:,:,:)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,1) = X0
       posvec(:,:,2) = Y0
       CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing 
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:) = Mesh%posvec%bcenter(:,:,:) - posvec(:,:,:)
       ! compute its absolute value
       radius(:,:) = SQRT(posvec(:,:,1)**2+posvec(:,:,2)**2)
    END IF

    ! curvilinear components of azimuthal unit vector
    ! (maybe with respect to shifted origin)
    ! from ephi = ez x er = ez x posvec/radius = ez x (rxi*exi + reta*eeta)/r
    !             = rxi/r*(ez x exi) + reta/r*(ez x eeta) = rxi/r*eeta - reta/r*exi
    ! because (ez,exi,eeta) is right handed orthonormal set of basis vectors
    ephi(:,:,1) = -posvec(:,:,2)/radius(:,:)
    ephi(:,:,2) = posvec(:,:,1)/radius(:,:)

    ! velocities
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    SELECT CASE(GetType(Physics))
    CASE(EULER2D,EULER2D_IAMROT)
       ! non-isothermal setup with constant density and pressure pulse
       Timedisc%pvar(:,:,Physics%DENSITY)   = RHO0
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
          Timedisc%pvar(i,j,Physics%PRESSURE) = P0 + AMP*EXP(-LOG(2.0) * &
               (radius(i,j)/PWIDTH)**2)
       END FORALL
    CASE(EULER2D_ISOTHERM)
       ! in isothermal configurations the pressure is proportional to
       ! the density; thus the pulse is applied to the density
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
          Timedisc%pvar(i,j,Physics%DENSITY) = RHO0 + AMP*EXP(-LOG(2.0) * &
               (radius(i,j)/PWIDTH)**2)
       END FORALL
    END SELECT

    ! check for rotating reference frame
    sp => GetSourcesPointer(Physics%sources,ROTATING_FRAME)
    IF ((ASSOCIATED(sp).EQV..TRUE.).OR.(WITH_IAR.AND.(OMEGA.GT.TINY(OMEGA)))) THEN
       ! set curvilinear velocity components with respect to rotating frame
       FORALL (i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX)
          Timedisc%pvar(i,j,Physics%XVELOCITY:Physics%YVELOCITY) = &
            -OMEGA * radius(i,j) * ephi(i,j,1:2)
       END FORALL
    END IF

    ! boundary conditions
    DO dir=WEST,EAST
       DO j=Mesh%JMIN,Mesh%JMAX
          DO ig=1,Mesh%GNUM
             SELECT CASE(dir)
             CASE(WEST)
                i = Mesh%IMIN-ig
             CASE(EAST)
                i = Mesh%IMAX+ig
             END SELECT
             SELECT CASE(GetType(Timedisc%Boundary(dir)))
             CASE(NOSLIP)
                Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%YVELOCITY)
             CASE(CUSTOM)
                Timedisc%boundary(dir)%cbtype(j,Physics%DENSITY) = CUSTOM_REFLECT
                Timedisc%boundary(dir)%cbtype(j,Physics%XVELOCITY) = CUSTOM_REFLNEG
                Timedisc%boundary(dir)%cbtype(j,Physics%YVELOCITY) = CUSTOM_FIXED
                Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%YVELOCITY)
             END SELECT
          END DO
       END DO
    END DO
    DO dir=SOUTH,NORTH
       DO i=Mesh%IMIN,Mesh%IMAX
          DO ig=1,Mesh%GNUM
             SELECT CASE(dir)
             CASE(SOUTH)
                j = Mesh%JMIN-ig
             CASE(NORTH)
                j = Mesh%JMAX+ig
             END SELECT
             SELECT CASE(GetType(Timedisc%Boundary(dir)))
             CASE(NOSLIP)
                Timedisc%Boundary(dir)%data(i,ig,Physics%XVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%XVELOCITY)
             CASE(CUSTOM)
                Timedisc%boundary(dir)%cbtype(i,Physics%DENSITY) = CUSTOM_REFLECT
                Timedisc%boundary(dir)%cbtype(i,Physics%XVELOCITY) = CUSTOM_FIXED
                Timedisc%boundary(dir)%cbtype(i,Physics%YVELOCITY) = CUSTOM_REFLNEG
                Timedisc%Boundary(dir)%data(i,ig,Physics%XVELOCITY) = &
                    Timedisc%pvar(i,j,Physics%XVELOCITY)
             END SELECT
          END DO
       END DO
    END DO

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "2D gaussian pressure pulse")

  END SUBROUTINE InitData

END PROGRAM gauss2d
