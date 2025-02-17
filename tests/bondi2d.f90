!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: bondi2d.f90                                                       #
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
  INTEGER, PARAMETER :: MGEO    = CYLINDRICAL ! geometry
  INTEGER, PARAMETER :: XRES    = 100         ! x-resolution
  INTEGER, PARAMETER :: YRES    = 6           ! y-resolution
  INTEGER, PARAMETER :: ZRES    = 1           ! z-resolution
  REAL, PARAMETER    :: RIN     = 0.1         ! inner/outer radii in terms of
  REAL, PARAMETER    :: ROUT    = 2.0         !   the Bondi radius RB, ROUT > 1
  REAL, PARAMETER    :: GPAR    = 1.0         ! geometry scaling parameter in [RB]
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &             ! output data dir
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &             ! output data file name
                     :: OFNAME  = 'bondi2d'
  ! some derives quandities
  REAL               :: RB                    ! Bondi radius
  REAL               :: TAU                   ! free fall time scale
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  REAL, DIMENSION(:), ALLOCATABLE :: sigma
  REAL :: sum_numer, sum_denom
  INTEGER :: n,DEN,VEL,PRE
  LOGICAL :: ok
  !--------------------------------------------------------------------------!

  ALLOCATE(Sim)
  CALL Sim%InitFosite()

#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
TAP_PLAN(4)
#ifdef PARALLEL
  END IF
#endif

  CALL MakeConfig(Sim, Sim%config)
!  CALL PrintDict(config)
  CALL Sim%Setup()
  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc)
  ! run the simulation
  CALL Sim%Run()
  ok = .NOT.Sim%aborted
  ! compare with exact solution if requested
  IF (ASSOCIATED(Sim%Timedisc%solution)) THEN
    ALLOCATE(sigma(Sim%Physics%VNUM))
    DO n=1,Sim%Physics%VNUM
      ! use L1 norm to estimate the deviation from the exact solution:
      !   Σ |pvar - pvar_exact| / Σ |pvar_exact|
      sum_numer = SUM(ABS(Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n) &
                        -Sim%Timedisc%solution%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n)))
      sum_denom = SUM(ABS(Sim%Timedisc%solution%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,&
                           Sim%Mesh%JMIN:Sim%Mesh%JMAX,Sim%Mesh%KMIN:Sim%Mesh%KMAX,n)))
#ifdef PARALLEL
      IF (Sim%GetRank().GT.0) THEN
        CALL MPI_Reduce(sum_numer,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
        CALL MPI_Reduce(sum_denom,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
      ELSE
        CALL MPI_Reduce(MPI_IN_PLACE,sum_numer,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
        CALL MPI_Reduce(MPI_IN_PLACE,sum_denom,1,DEFAULT_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,Sim%ierror)
#endif
      sigma(n) = sum_numer / sum_denom
#ifdef PARALLEL
      END IF
#endif
    END DO
  ELSE
    sigma(:) = 0.0
  END IF

  DEN = Sim%Physics%DENSITY
  VEL = Sim%Physics%XVELOCITY
  PRE = Sim%Physics%PRESSURE

#ifdef PARALLEL
  IF (Sim%GetRank().EQ.0) THEN
#endif
TAP_CHECK(ok,"stoptime reached")
! These lines are very long if expanded. So we can't indent it or it will be cropped.
TAP_CHECK_SMALL(sigma(DEN),1.0E-02,"density deviation < 1%")
TAP_CHECK_SMALL(sigma(VEL),1.0E-02,"radial velocity deviation < 1%")
! skip azimuthal velocity deviation, because exact value is 0
TAP_CHECK_SMALL(sigma(PRE),1.0E-02,"pressure deviation < 1%")
TAP_DONE
#ifdef PARALLEL
  END IF
#endif

  IF (ALLOCATED(sigma)) DEALLOCATE(sigma)
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

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
    CASE(CYLINDRICAL)
       x1 = RIN * RB
       x2 = ROUT * RB
       y1 = 0.0
       y2 = 2*PI
       z1 = 0.0
       z2 = 0.0
       bc(WEST)   = NO_GRADIENTS
       bc(EAST)   = CUSTOM
       bc(SOUTH)  = PERIODIC
       bc(NORTH)  = PERIODIC
       bc(BOTTOM) = NO_GRADIENTS
       bc(TOP)    = NO_GRADIENTS
    CASE DEFAULT
       CALL Sim%Physics%Error("bondi2d::MakeConfig","geometry currently not supported")
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
             "problem"  / EULER, &
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
             "output/solution" / 1, &
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
    INTEGER              :: i,j,k,l
    CHARACTER(LEN=64)    :: info_str
    !------------------------------------------------------------------------!
    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_euler)
      ! constant density and pressure, vanishing velocity
      pvar%density%data1d(:)  = RHOINF
      pvar%velocity%data1d(:) = 0.
      pvar%pressure%data1d(:) = RHOINF * CSINF**2 / GAMMA
    CLASS DEFAULT
      CALL Physics%Error("bondi2d::InitData","only non-isothermal HD supported")
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    ! store exact stationary solution if requested,
    ! i.e. output/solution == 1 in timedisc config
    IF (ASSOCIATED(Timedisc%solution)) THEN
      ! compute the stationary 2D planar bondi solution
      SELECT TYPE(pvar => Timedisc%solution)
      TYPE IS(statevector_euler)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              r = Mesh%radius%bcenter(i,j,k)
              CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
              cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
              pvar%density%data3d(i,j,k)  = rho
              ! set the minimum of the exact velocity to some value > 0
              ! for comparison with the numerical results
!NEC$ UNROLL(3)
              DO l=1,Physics%VDIM
                pvar%velocity%data4d(i,j,k,l) = SIGN(&
                  MAX(EPSILON(cs2),ABS(vr*Mesh%posvec%bcenter(i,j,k,l)/r)), &
                    vr*Mesh%posvec%bcenter(i,j,k,l))
              END DO
              pvar%pressure%data3d(i,j,k) = rho * cs2 / GAMMA
            END DO
          END DO
        END DO
      END SELECT
    END IF

    ! boundary condition: subsonic inflow according to Bondi's solution
    ! calculate Bondi solution for y=ymin..ymax at xmax
    SELECT TYPE(beast => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS (boundary_custom)
      ! this tells the boundary routine the actual boundary condition for each variable
      CALL beast%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_FIXED,CUSTOM_NOGRAD,CUSTOM_FIXED,CUSTOM_FIXED/))
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GINUM
            ! get distance to the origin for each boundary cell
            r = Mesh%radius%bcenter(Mesh%IMAX+i,j,k)
            CALL bondi(r/RB,GAMMA,RHOINF,CSINF,rho,vr)
            cs2 = CSINF**2 * (rho/RHOINF)**(GAMMA-1.0)
            beast%data(i,j,k,Physics%DENSITY)   = rho
            DO l=1,Physics%VDIM
              beast%data(i,j,k,Physics%XVELOCITY+l-1) &
                = vr*Mesh%posvec%bcenter(Mesh%IMAX+i,j,k,l)/r
            END DO
            beast%data(i,j,k,Physics%PRESSURE)  = rho * cs2 / GAMMA
          END DO
        END DO
      END DO
    END SELECT

    ! print some information
    CALL Mesh%Info(" DATA-----> initial condition: 2D Bondi accretion")
    WRITE(info_str,"(ES11.3)") RB
    CALL Mesh%Info("                               " // "Bondi radius:       " &
         // TRIM(info_str) // " m")
    WRITE(info_str,"(ES11.3)") TAU
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




