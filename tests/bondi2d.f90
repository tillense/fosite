!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: bondi2d.f03                                                       #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
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
!> \test Bondi 2D accretion
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief Program and data initialization for 2D Bondi accretion
!!
!! References:
!! - \cite bondi1952 Bondi, H.: On spherically symmetrical accretion,
!!    Mon. Not. Roy. Astron. Soc., 112 (1952)
!!    ADS link: http://adsabs.harvard.edu/abs/1952MNRAS.112..195B
!! - \cite padmanabhan2000 Padmanabhan, T.:Theoretical Astrophysics,
!!    Vol. I: Astrophysical Processes, Cambridge University Press (2000),
!!    Chapter 8.9.2
!!
!! \warning compile with autodouble
!----------------------------------------------------------------------------!
PROGRAM bondi2d
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  INTERFACE
     PURE SUBROUTINE funcd(x,fx,dfx,plist)
       IMPLICIT NONE
       REAL, INTENT(IN)  :: x
       REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
       REAL, INTENT(OUT) :: fx,dfx
     END SUBROUTINE funcd
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: MSUN    = 1.989E+30   ! solar mass [kg]
  REAL, PARAMETER    :: GN      = 6.67408E-11 ! gravitional constant (SI)
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 20.0        ! simulation time [TAU] (free fall)
  REAL, PARAMETER    :: ACCMASS = 1.0*MSUN    ! mass of the accreting object
  REAL, PARAMETER    :: GAMMA   = 1.4         ! ratio of specific heats
  ! boundary conditions
  REAL, PARAMETER    :: RHOINF  = 1.0E-20     ! density at infinity [kg/m^3]
  REAL, PARAMETER    :: CSINF   = 1.0E+04     ! sound speed at infinity [m/s]
  ! mesh settings
  INTEGER, PARAMETER :: MGEO    = POLAR       ! geometry
  INTEGER, PARAMETER :: XRES    = 50          ! x-resolution
  INTEGER, PARAMETER :: YRES    = 10          ! y-resolution
  INTEGER, PARAMETER :: ZRES    = 1           ! z-resolution
  REAL, PARAMETER    :: RIN     = 0.1         ! inner/outer radii in terms of
  REAL, PARAMETER    :: ROUT    = 2.0         !   the Bondi radius RB, ROUT > 1
  REAL, PARAMETER    :: GPAR    = 1.0         ! geometry scaling parameter in [RB]
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME  = 'bondi2d'
  ! some derives quandities
  REAL               :: RB                    ! Bondi radius
  REAL               :: TAU                   ! free fall time scale
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
!  CALL PrintDict(config)
  CALL Sim%Setup()
  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc)
  CALL Sim%Run()

  DEALLOCATE(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: Sim
    TYPE(Dict_TYP), POINTER      :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                      :: bc(6)
    TYPE(Dict_TYP), POINTER      :: mesh, physics, boundary, datafile, &
                                    grav, pmass, timedisc, fluxes
    REAL                         :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    ! derived constants
    RB  = GN * ACCMASS / CSINF**2           ! bondi radius [m]               !
    TAU = RB / CSINF                        ! free fall time scale [s]       !

    ! mesh settings
    SELECT CASE(MGEO)
    CASE(POLAR)
       x1 = RIN * RB
       x2 = ROUT * RB
       y1 = 0.0
       y2 = 2*PI
       z1 = 0.0
       z2 = 0.0
    CASE DEFAULT
       CALL Sim%Physics%Error("InitProgram","mesh geometry not supported for 2D Bondi accretion")
    END SELECT
    mesh => Dict( &
             "meshtype"  / MIDPOINT, &
             "geometry"  / MGEO, &
             "inum"      / XRES, &
             "jnum"      / YRES, &
             "knum"      / ZRES, &
             "xmin"      / x1, &
             "xmax"      / x2, &
             "ymin"      / y1, &
             "ymax"      / y2, &
             "zmin"      / z1, &
             "zmax"      / z2, &
             "gparam"    / (GPAR*RB))

    SELECT CASE(MGEO)
    CASE(POLAR)
       bc(WEST)   = ABSORBING
       bc(EAST)   = FIXED
       bc(SOUTH)  = PERIODIC
       bc(NORTH)  = PERIODIC
       bc(BOTTOM) = NO_GRADIENTS
       bc(TOP)    = NO_GRADIENTS
    CASE DEFAULT
       CALL Sim%Physics%Error("InitProgram","mesh geometry not supported for 2D Bondi accretion")
    END SELECT

    ! boundary conditions
    boundary => Dict( &
             "western"   / bc(WEST), &
             "eastern"   / bc(EAST), &
             "southern"  / bc(SOUTH), &
             "northern"  / bc(NORTH), &
             "bottomer"  / bc(BOTTOM), &
             "topper"    / bc(TOP))

    ! physics settings
    physics => Dict( &
             "problem"  / EULER2D, &
             "gamma"    / GAMMA)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
             "order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &
             "limiter"   / MONOCENT, &
             "theta"     / 1.2)

    ! gravity term due to a point mass
    pmass => Dict( &
             "gtype"     / POINTMASS, &
             "mass"      / ACCMASS, &
             "outbound"  / 0)

   ! source term due to all gravity terms
    grav => Dict( &
             "stype"     / GRAVITY, &
             "pmass"     / pmass)

    ! time discretization settings
    timedisc => Dict( &
             "method"    / MODIFIED_EULER, &
             "order"     / 3, &
             "cfl"       / 0.4, &
             "stoptime"  / (TSIM * TAU), &
             "dtlimit"   / (1.0E-6 * TAU), &
             "maxiter"   / 1000000)

    ! initialize data input/output
    datafile => Dict( &
            "fileformat" / VTK, &
            "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
            "count"      / ONUM)

    config => Dict( &
            "mesh"       / mesh, &
            "physics"    / physics, &
            "boundary"   / boundary, &
            "fluxes"     / fluxes, &
            "sources/grav"  / grav, &
            "timedisc"   / timedisc, &
            "datafile"   / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Fluxes,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(fluxes_base),   INTENT(IN)    :: Fluxes
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL                 :: r,rho,vr,cs2
    INTEGER              :: i,j,k
    CHARACTER(LEN=64)    :: info_str
    !------------------------------------------------------------------------!
    ! Bondi solution at the Bondi radius RB
!!$    CALL bondi(1.0,GAMMA,RHOINF,CSINF,rho,vr)
!!$    cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
    ! take values at infinity as initial condition
    rho = RHOINF
    cs2 = CSINF**2

    ! initial condition
    Timedisc%pvar%data4d(:,:,:,Physics%DENSITY)   = rho
    Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE)  = rho * cs2 / GAMMA
    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar%data4d,Timedisc%cvar%data4d)

    ! boundary condition: subsonic inflow according to Bondi's solution
    ! calculate Bondi solution for y=ymin..ymax at xmax
    IF ((Timedisc%Boundary%Boundary(EAST)%p%GetType()).EQ.FIXED) THEN
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
            ! get distance to the origin for each boundary cell
            r = Mesh%radius%bcenter(Mesh%IMAX+i,j,k)
            CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
            cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
            ! set boundary data to either primitive or conservative values
            ! depending on the reconstruction
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%DENSITY)   = rho
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%XVELOCITY) = vr
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%YVELOCITY) = 0.
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
          ! this tells the boundary routine which values to fix (.TRUE.)
          ! and which to extrapolate (.FALSE.)
          Timedisc%Boundary%Boundary(EAST)%p%fixed(j,k,:) = (/ .TRUE., .FALSE., .TRUE., .TRUE. /)
        END DO
      END DO
    END IF

    CALL Mesh%Info(" DATA-----> initial condition: 2D Bondi accretion")
    WRITE(info_str,"(ES9.3)") RB
    CALL Mesh%Info("                               " // "Bondi radius:       " &
         // TRIM(info_str) // " m")
    WRITE(info_str,"(ES9.3)") TAU
    CALL Mesh%Info( "                               " // "Free fall time:     " &
         // TRIM(info_str) // " s")
  END SUBROUTINE InitData

  SUBROUTINE bondi(r,gamma,rhoinf,csinf,rho,vr)
    USE roots
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    ! computes the Bondi solution for rotationally symmetric accretion in    !
    !   planar geometry                                                      !
    ! uses funcd, GetRoot_newton                                             !
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
    REAL, PARAMETER :: xacc = 1.0E-6     ! accuracy for root finding
    REAL :: gp1,gm1,rc,chi,lambda,psi,gr
    REAL :: params(2)
    INTEGER :: err
!     COMMON /funcd_parameter/ gm1, gr
    !------------------------------------------------------------------------!
    ! for convenience
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0

    ! critical radius
    rc = (3.0-gamma) / 2.0

    ! critical dimensionless accretion rate
    lambda = rc**(-rc/gm1)

    ! Newton-Raphson to solve Bondis equations for psi
    chi = r / lambda
    gr  = chi**(2.*gm1/gp1) * (1./gm1 + 1./r)
    params(1) = gm1
    params(2) = gr
    IF (r.LT.rc) THEN
       CALL GetRoot_newton(funcd,1.0,gr,psi,err,plist=params)
    ELSE
       CALL GetRoot_newton(funcd,1.0E-5,1.0,psi,err,plist=params)
    END IF

    ! return values
    rho = rhoinf * chi**(-2./gp1) / psi        ! density
    vr  = -csinf * chi**(-gm1/gp1) * psi       ! radial velocity
  END SUBROUTINE bondi

END PROGRAM bondi2d

! find the root of fy to compute the exact Bondi solution
PURE SUBROUTINE funcd(y,fy,dfy,plist)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
  REAL, INTENT(IN)  :: y
  REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
  REAL, INTENT(OUT) :: fy,dfy
  !------------------------------------------------------------------------!
  ! plist(1) = gamma - 1
  ! plist(2) = gr
  !------------------------------------------------------------------------!
  fy  = 0.5*y*y + y**(-plist(1)) / plist(1) - plist(2)
  dfy = y - y**(-plist(1)-1.)
END SUBROUTINE funcd




