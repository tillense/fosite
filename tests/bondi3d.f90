!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: bondi3d.f90                                                       #
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
!> 3D Bondi accretion
!! \author Tobias Illenseer
!!
!! References:
!! [1] Bondi, H.: On spherically symmetrical accretion,
!!     Mon. Not. Roy. Astron. Soc., 112 (1951)
!!     ADS link: http://adsabs.harvard.edu/abs/1952MNRAS.112..195B
!! [2] Padmanabhan, T.:Theoretical Astrophysics, Vol. I: Astrophysical
!!     Processes, Cambridge University Press (2000), Chapter 8.9.2
!!
!! \warning compile with autodouble
!----------------------------------------------------------------------------!
PROGRAM bondi3d
  USE fosite
  USE constants_common, ONLY : GN        ! gravitational constant in SI units
  USE constants_generic
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
  ! general constants
  REAL, PARAMETER    :: MSUN = 1.989E+30   ! solar mass [kg]                 !
  ! simulation parameters
  REAL, PARAMETER    :: TSIM = 10.0        ! simulation time [TAU] (free fall)
  REAL, PARAMETER    :: ACCMASS = 1.0*MSUN ! mass of the accreting object    !
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats         !
  ! boundary conditions
  REAL, PARAMETER    :: RHOINF = 3.351E-17 ! density at infinity [kg/m^3]    !
  REAL, PARAMETER    :: CSINF = 537.0      ! sound speed at infinity [m/s]   !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = SPHERICAL   ! geometry of the mesh            !
!!$  INTEGER, PARAMETER :: MGEO = OBLATE_SPHEROIDAL
!   INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = TANCYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = SINHSPHERICAL
  INTEGER, PARAMETER :: XRES = 100         ! x-resolution                    !
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution                    !
  REAL, PARAMETER    :: RIN  = 0.1         ! inner/outer radii in terms of   !
  REAL, PARAMETER    :: ROUT = 2.0         !   the Bondi radius RB, ROUT > 1 !
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'bondi3d'
  ! some derives quandities
  REAL               :: RB                 ! Bondi radius
  REAL               :: TAU                ! free fall time scale
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
  TAP_DONE

CONTAINS

 SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               grav, pmass, timedisc, fluxes
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! derived constants
    RB  = GN * ACCMASS / CSINF**2  ! bondi radius [m]      !
    TAU = RB / CSINF                        ! free fall time scale [s]       !

    ! geometry dependent setttings
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       x1 = RIN * RB
       x2 = ROUT * RB
       y1 = 0.0
       y2 = PI
    CASE(SINHSPHERICAL)
       x1 = LOG(RIN+SQRT(RIN**2+1.0))   ! = ASINH(RIN)
       x2 = LOG(ROUT+SQRT(ROUT**2+1.0)) ! = ASINH(ROUT)
       y1 = 0.0
       y2 = PI
    CASE(CYLINDRICAL)
       x1 = -ROUT * RB
       x2 = ROUT * RB
       y1 = RIN * RB
       y2 = ROUT * RB
    CASE(OBLATE_SPHEROIDAL)
       x1 = LOG(RIN+SQRT(RIN**2+1.0))   ! = ASINH(RIN)
       x2 = LOG(ROUT+SQRT(ROUT**2-1.0)) ! = ACOSH(ROUT)
       y1 = -0.5*PI
       y2 = 0.5*PI
    CASE(TANCYLINDRICAL)
       x1 = ATAN(-ROUT)
       x2 = ATAN(ROUT)
       y1 = RIN * RB
       y2 = ROUT * RB
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for Bondi accretion")
    END SELECT

    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "gparam"   / RB)

    ! geometry dependent setttings
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = FARFIELD
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(SINHSPHERICAL)
       bc(WEST)  = ABSORBING
       bc(EAST)  = FARFIELD
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(CYLINDRICAL)
       bc(WEST)  = FARFIELD
       bc(EAST)  = FARFIELD
       bc(SOUTH) = ABSORBING
       bc(NORTH) = FARFIELD
    CASE(OBLATE_SPHEROIDAL)
       bc(WEST)  = FARFIELD
       bc(EAST)  = FARFIELD
       bc(SOUTH) = AXIS
       bc(NORTH) = AXIS
    CASE(TANCYLINDRICAL)
       bc(WEST)  = FARFIELD
       bc(EAST)  = FARFIELD
       bc(SOUTH) = FARFIELD
       bc(NORTH) = FARFIELD
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for Bondi accretion")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / GAMMA)                 ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! gravity term due to a point mass
    pmass => Dict("gtype"    / POINTMASS, &   ! grav. accel. of a point mass     !
            "mass"     / ACCMASS, &           ! mass of the accreting object[kg] !
            "outbound" / 0)                   ! disable accretion

   ! source term due to all gravity terms
    grav => Dict("stype"    / GRAVITY, &
                 "pmass" / pmass)
              
    ! time discretization settings
    timedisc => Dict("method"    / MODIFIED_EULER, &
               "order"     / 3, &
               "cfl"       / 0.4, &
               "stoptime"  / (TSIM * TAU), &
               "dtlimit"   / (1.0E-6 * TAU), &
               "maxiter"   / 1000000)

    ! initialize log input/output
!     logfile => Dict("fileformat"  / BINARY, &
!              +  "filecycles" / 1, &
!                "filename"    / (TRIM(ODIR) // TRIM(OFNAME))

    ! initialize data input/output
    datafile => Dict("fileformat"  / GNUPLOT, &
               "filecycles" / 0, &
               "filename"    / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"       / ONUM)
    
    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources/grav"  / grav, &
             "timedisc" / timedisc, &
!              "logfile"  / logfile, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Fluxes,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: r,rho,vr,cs2
    INTEGER           :: i,j
    CHARACTER(LEN=64) :: info_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition: use data at infinity
    Timedisc%pvar(:,:,Physics%DENSITY)   = RHOINF
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%ZVELOCITY) = 0.
    Timedisc%pvar(:,:,Physics%PRESSURE)  = RHOINF * CSINF**2 / GAMMA
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    ! boundary conditions
    ! outflow condition, this is only for transition; in the end the flow
    ! becomes supersonic and all variables are extrapolated from the interior
    IF (GetType(Timedisc%Boundary(WEST)).EQ.FARFIELD) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             r = Mesh%radius%bcenter(Mesh%IMIN-i,j)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             Timedisc%Boundary(WEST)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(WEST)%data(i,j,Physics%XVELOCITY) = vr*Mesh%posvec%bcenter(Mesh%IMIN-i,j,1)/r
             Timedisc%Boundary(WEST)%data(i,j,Physics%YVELOCITY) = vr*Mesh%posvec%bcenter(Mesh%IMIN-i,j,2)/r
             Timedisc%Boundary(WEST)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(WEST)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
    END IF
    ! subsonic inflow according to Bondi's solution
    ! calculate Bondi solution for y=ymin..ymax at xmax
    IF (GetType(Timedisc%Boundary(EAST)).EQ.FARFIELD) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             r = Mesh%radius%bcenter(Mesh%IMAX+i,j)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(EAST)%data(i,j,Physics%XVELOCITY) = vr*Mesh%posvec%bcenter(Mesh%IMAX+i,j,1)/r
             Timedisc%Boundary(EAST)%data(i,j,Physics%YVELOCITY) = vr*Mesh%posvec%bcenter(Mesh%IMAX+i,j,2)/r
             Timedisc%Boundary(EAST)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(EAST)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
    END IF
    IF (GetType(Timedisc%Boundary(SOUTH)).EQ.FARFIELD) THEN
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             r = Mesh%radius%bcenter(i,Mesh%JMIN-j)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%XVELOCITY) = vr*Mesh%posvec%bcenter(i,Mesh%JMIN-j,1)/r
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%YVELOCITY) = vr*Mesh%posvec%bcenter(i,Mesh%JMIN-j,2)/r
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(SOUTH)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
    END IF
    IF (GetType(Timedisc%Boundary(NORTH)).EQ.FARFIELD) THEN
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             r = Mesh%radius%bcenter(i,Mesh%JMAX+j)
             CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
             cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
             Timedisc%Boundary(NORTH)%data(i,j,Physics%DENSITY)   = rho
             Timedisc%Boundary(NORTH)%data(i,j,Physics%XVELOCITY) = vr*Mesh%posvec%bcenter(i,Mesh%JMAX+j,1)/r
             Timedisc%Boundary(NORTH)%data(i,j,Physics%YVELOCITY) = vr*Mesh%posvec%bcenter(i,Mesh%JMAX+j,2)/r
             Timedisc%Boundary(NORTH)%data(i,j,Physics%ZVELOCITY) = 0.
             Timedisc%Boundary(NORTH)%data(i,j,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
       END DO
    END IF
    CALL Info(Mesh," DATA-----> initial condition: 3D Bondi accretion")
    WRITE(info_str,"(ES9.3)") RB
    CALL Info(Mesh, "                               " // "Bondi radius:       " &
         // TRIM(info_str) // " m")
    WRITE(info_str,"(ES9.3)") TAU
    CALL Info(Mesh, "                               " // "Free fall time:     " &
         // TRIM(info_str) // " s")
  END SUBROUTINE InitData


  SUBROUTINE bondi(r,gamma,rhoinf,csinf,rho,vr)
    USE roots
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! computes the Bondi solution for spherically symmetric accretion        !
    ! uses funcd, GetRoot                                                    !
    ! INPUT paramter:                                                        !
    !   r     : radius in units of the Bondi radius r_b = G*M/csinf**2       !
    !   gamma : ratio of specific heats (1 < gamma < 5/3)                    !
    !   rhoinf: density at infinity                                          !
    !   csinf : speed of sound at infinity                                   !
    ! OUTPUT paramter:                                                       !
    !   rho   : density @ r                                                  !
    !   vr    : radial velocity @ r                                          !
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,gamma,rhoinf,csinf
    REAL, INTENT(OUT) :: rho,vr
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE funcd(x,fx,dfx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx,dfx
       END SUBROUTINE funcd
    END INTERFACE
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: xacc = 1.0E-6     ! accuracy for root finding
    REAL :: gp1,gm1,g35,rc,chi,lambda,psi,gr
    INTEGER :: err
    COMMON /funcd_parameter/ gm1, gr
    !------------------------------------------------------------------------!
    ! for convenience
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    g35 = 0.5 * (5.0-3.0*gamma)

    ! critical radius
    rc = 0.5 * g35
    ! critical dimensionless accretion rate
    lambda = 0.25 * g35**(-g35/gm1)

    ! Newton-Raphson to solve Bondis equations for psi
    chi = r**2 / lambda
    gr  = chi**(2.*gm1/gp1) * (1./gm1 + 1./r)
    IF (r.LT.rc) THEN
       CALL GetRoot_newton(funcd,1.0,gr,psi,err)
    ELSE
       CALL GetRoot_newton(funcd,1.0E-6,1.0,psi,err)
    END IF
    
    ! return values
    rho = rhoinf * chi**(-2./gp1) / psi        ! density
    vr  = -csinf * chi**(-gm1/gp1) * psi       ! radial velocity
  END SUBROUTINE bondi


END PROGRAM bondi3d


! for exact Bondi solution at the outer boundary
PURE SUBROUTINE funcd(y,fy,dfy,plist)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
  REAL, INTENT(IN)  :: y
  REAL, INTENT(OUT) :: fy,dfy
  REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
  !------------------------------------------------------------------------!
  REAL :: gm1,gr
  COMMON /funcd_parameter/ gm1,gr
  !------------------------------------------------------------------------!
  fy  = 0.5*y*y + y**(-gm1) / gm1 - gr
  dfy = y - y**(-gm1-1.)
END SUBROUTINE funcd
