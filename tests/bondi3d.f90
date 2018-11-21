!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: bondi3d.f03                                                       #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
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
!> \test Bondi 3D accretion
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief Program and data initialization for 3D Bondi accretion
!!
!! References:
!! - \cite bondi1952 Bondi, H.: On spherically symmetrical accretion,
!!     Mon. Not. Roy. Astron. Soc., 112 (1951)
!!     ADS link: http://adsabs.harvard.edu/abs/1952MNRAS.112..195B
!! - \cite padmanabhan2000 Padmanabhan, T.:Theoretical Astrophysics, Vol. I:
!!     Astrophysical Processes, Cambridge University Press (2000), Chapter 8.9.2
!!
!! \warning compile with autodouble
!----------------------------------------------------------------------------!
PROGRAM bondi3d
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: MSUN    = 1.989E+30   ! solar mass [kg]              !
  REAL, PARAMETER    :: GN      = 6.67408E-11 ! gravitional constant (SI)    !
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 2.0        ! simulation time [TAU]        !
  REAL, PARAMETER    :: ACCMASS = 1.0*MSUN    ! mass of the accreting object !
  REAL, PARAMETER    :: GAMMA   = 1.4         ! ratio of specific heats      !
  ! boundary conditions
  REAL, PARAMETER    :: RHOINF  = 3.351E-17   ! density at infinity [kg/m^3] !
  REAL, PARAMETER    :: CSINF   = 537.0       ! sound speed at infinity [m/s]!
  ! mesh settings
  INTEGER, PARAMETER :: MGEO    = SPHERICAL   ! geometry of the mesh         !
  INTEGER, PARAMETER :: XRES    = 50         ! x-resolution                 !
  INTEGER, PARAMETER :: YRES    = 6          ! y-resolution                 !
  INTEGER, PARAMETER :: ZRES    = 12          ! z-resolution                 !
  REAL, PARAMETER    :: RIN     = 0.1         ! inner/outer radii in terms of!
  REAL, PARAMETER    :: ROUT    = 2.0         ! the Bondi radius RB, ROUT > 1!
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME  = 'bondi3d'
  ! some derives quandities
  REAL               :: RB                    ! Bondi radius
  REAL               :: TAU                   ! free fall time scale
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  CLASS(marray_compound), POINTER :: pvar_exact
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
!  CALL PrintDict(Sim%config)
  CALL Sim%Setup()
  ! set initial condition
  CALL InitData(Sim%Mesh,Sim%Physics,Sim%Timedisc)
  ! compute numerical solution
  CALL Sim%Run()
  ! compute exact solution
  CALL new_statevector(Sim%Physics,pvar_exact,PRIMITIVE)
  CALL ExactSolution(Sim%Mesh,Sim%Physics,pvar_exact)
!> \todo implement custom or farfield boundary conditions and compare with analytic solution
  CALL pvar_exact%Destroy()
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

  TAP_CHECK(.TRUE.,"Finished simulation")
  TAP_DONE

CONTAINS

 SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP), POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                 :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               grav, pmass, timedisc, fluxes
    REAL                    :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: Sim
    !------------------------------------------------------------------------!
    ! derived constants
    RB  = GN * ACCMASS / CSINF**2            ! bondi radius [m]              !
    TAU = RB / CSINF                         ! free fall time scale [s]      !

    ! geometry dependent setttings
    SELECT CASE(MGEO)
    CASE(SPHERICAL)
      x1 = RIN * RB
      x2 = ROUT * RB
      y1 = 0.0
      y2 = PI
      z1 = 0.0
      z2 = 2*PI
      bc(WEST)   = NO_GRADIENTS
!       bc(EAST)   = FIXED
      bc(EAST)   = REFLECTING
      bc(SOUTH)  = AXIS
      bc(NORTH)  = AXIS
      bc(BOTTOM) = PERIODIC
      bc(TOP)    = PERIODIC
    CASE DEFAULT
       CALL Sim%Physics%Error("bondi3d::MakeConfig","geometry currently not supported")
    END SELECT

    mesh => Dict( &
          "meshtype"    / MIDPOINT, &
          "geometry"    / MGEO, &
          "inum"        / XRES, &
          "jnum"        / YRES, &
          "knum"        / ZRES, &
          "xmin"        / x1, &
          "xmax"        / x2, &
          "ymin"        / y1, &
          "ymax"        / y2, &
          "zmin"        / z1, &
          "zmax"        / z2, &
          "output/radius" / 1, &
          "gparam"      / RB)

    ! boundary conditions
    boundary => Dict( &
          "western"     / bc(WEST), &
          "eastern"     / bc(EAST), &
          "southern"    / bc(SOUTH), &
          "northern"    / bc(NORTH), &
          "bottomer"    / bc(BOTTOM), &
          "topper"      / bc(TOP))

    ! physics settings
    physics => Dict( &
          "problem"     / EULER, &
          "gamma"       / GAMMA)

    ! flux calculation and reconstruction method
    fluxes => Dict( &
          "order"       / LINEAR, &
          "fluxtype"    / KT, &
          "variables"   / PRIMITIVE, &
          "limiter"     / MONOCENT, &
          "theta"       / 1.2)

    ! gravity term due to a point mass
    pmass => Dict( &
          "gtype"       / POINTMASS, &
          "mass"        / ACCMASS, &
          "outbound"    / 0)

   ! source term due to all gravity terms
    grav => Dict( &
          "stype"       / GRAVITY, &
          "pmass"       / pmass)

    ! time discretization settings
    timedisc => Dict( &
          "method"      / MODIFIED_EULER, &
          "order"       / 3, &
          "cfl"         / 0.4, &
          "stoptime"    / (TSIM * TAU), &
          "dtlimit"     / (1.0E-6 * TAU), &
          "maxiter"     / 1000000)

    ! initialize data input/output
    datafile => Dict( &
          "fileformat"  / VTK, &
          "filename"    / (TRIM(ODIR) // TRIM(OFNAME)), &
          "count"       / ONUM)

    config => Dict( &
          "mesh"        / mesh, &
          "physics"     / physics, &
          "boundary"    / boundary, &
          "fluxes"      / fluxes, &
          "sources/grav"  / grav, &
          "timedisc"    / timedisc, &
          "datafile"    / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL              :: r,rho,vr,cs2
    INTEGER           :: i,j,k
    CHARACTER(LEN=64) :: info_str
    !------------------------------------------------------------------------!
    ! initial condition: use data at infinity
    Timedisc%pvar%data4d(:,:,:,Physics%DENSITY)   = RHOINF
    Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
    Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE)  = RHOINF * CSINF**2 / GAMMA

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    IF ((Timedisc%Boundary%Boundary(EAST)%p%GetType()).EQ.FIXED) THEN
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
            r = Mesh%radius%bcenter(Mesh%IMAX+i,j,k)
            CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
            cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%DENSITY)   = rho
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%XVELOCITY) = vr*Mesh%posvec%bcenter(Mesh%IMAX+i,j,k,1)/r
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%YVELOCITY) = vr*Mesh%posvec%bcenter(Mesh%IMAX+i,j,k,2)/r
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%ZVELOCITY) = 0.
            Timedisc%Boundary%Boundary(EAST)%p%data(i,j,k,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
        END DO
      END DO
    END IF


    ! boundary conditions
    ! outflow condition, this is only for transition; in the end the flow
    ! becomes supersonic and all variables are extrapolated from the interior
    CALL Mesh%Info(" DATA-----> initial condition: 3D Bondi accretion")
    WRITE(info_str,"(ES9.3)") RB
    CALL Mesh%Info("                               " // "Bondi radius:       " &
         // TRIM(info_str) // " m")
    WRITE(info_str,"(ES9.3)") TAU
    CALL Mesh%Info("                               " // "Free fall time:     " &
         // TRIM(info_str) // " s")
  END SUBROUTINE InitData


  SUBROUTINE ExactSolution(Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(marray_compound), POINTER     :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    REAL    :: r,rho,vr,cs2
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_euler) ! non-isothermal HD
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            r = Mesh%radius%bcenter(i,j,k)
            CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
            cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
            p%density%data3d(i,j,k)  = rho
            p%velocity%data4d(i,j,k,1:Physics%VDIM) = vr*Mesh%posvec%bcenter(i,j,k,1:Physics%VDIM)/r
            p%pressure%data3d(i,j,k) = rho * cs2 / GAMMA
          END DO
        END DO
      END DO
    CLASS DEFAULT
      CALL Physics%Error("bondi3d::ExactSolution","only non-isothermal HD supported")
    END SELECT

  END SUBROUTINE ExactSolution


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
