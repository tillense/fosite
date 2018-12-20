!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: orbitingcylinders.f90                                             #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> self-gravitating orbiting cylinder test from
!! Chan, Chi-kwan; Psaltis, Dimitrios; Özel, Feryal, 2006
!! Spectral Methods for Time-dependent Studies of Accretion Flows. II. 
!! Two-dimensional Hydrodynamic Disks with Self-Gravity
!! http://adsabs.harvard.edu/abs/2006ApJ...645..506C
!----------------------------------------------------------------------------!
!#define HAVE_FGSL
PROGRAM Init
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
#ifdef HAVE_FGSL
  USE fgsl
#endif
  USE functions, ONLY : Ei,Bessel_I0,Bessel_I1,Bessel_K0,Bessel_K1,Bessel_K0e
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! problem setup
  REAL, PARAMETER :: sigma       = 0.1
  REAL, DIMENSION(2), PARAMETER :: x1 = (/ 1.0, 0. /)
  REAL, DIMENSION(2), PARAMETER :: x2 = (/ 1.0, PI /)
  REAL, DIMENSION(2), PARAMETER :: x3 = (/ 0.9, 3./4.*PI /)
  REAL, DIMENSION(2), PARAMETER :: x4 = (/ 1.0, -PI/2. /)

  ! general constants
  REAL, PARAMETER :: GN          = 1.
  REAL, PARAMETER :: P0          = 1.
  REAL, PARAMETER :: P           = 10.
  REAL, PARAMETER :: t           = P*P0

  REAL, PARAMETER :: flaring     = 0.05
  REAL, PARAMETER :: RG          = 8.31
  REAL, PARAMETER :: MU          = 6.02E-04
  REAL, PARAMETER :: GAMMA       = 5./3.

  REAL, PARAMETER :: RMIN        = 0.4
  REAL, PARAMETER :: RMAX        = 2.0
  INTEGER, PARAMETER :: YRES     = 1024
  INTEGER, PARAMETER :: XRES     = YRES/4
  INTEGER, PARAMETER :: ONUM     = P * 1
  REAL, PARAMETER    :: omega    = 0.0
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP) :: Sim
  INTEGER          :: green
  LOGICAL          :: break
  REAL, DIMENSION(:,:), POINTER :: numpot, anapot
  REAL :: a,b
  !--------------------------------------------------------------------------!
 

  DO green=1,3

    CALL InitFosite(Sim)

    CALL MakeConfig(Sim%config)

    SELECT CASE(green)
    CASE(3)
      CALL SetAttr(Sim%config, "/mesh/rmin", 0.2)
      CALL SetAttr(Sim%config, "/mesh/rmax", 1.8)
      CALL SetAttr(Sim%config, "/mesh/inum", 64)
      CALL SetAttr(Sim%config, "/mesh/jnum", 192)
    END SELECT

    !CALL PrintDict(Sim%config)

    CALL SetupFosite(Sim)

    ! set initial condition
    CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics)

    SELECT CASE(green)
    CASE(1,2)
      CALL FirstStepFosite(Sim)
      print *,"max local relative error: ",MAXVAL(ABS((numpot-anapot)/anapot))
      print *,"global relative error: ",SUM(ABS(numpot-anapot))/SUM(ABS(anapot))
    CASE(3)
      CALL RunFosite(Sim)
    END SELECT

  END DO

  CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              fluxes,physics,grav
    !------------------------------------------------------------------------!
    CHARACTER(LEN=1)  :: str
    !------------------------------------------------------------------------!
    physics => Dict( &
        "problem" / EULER2D,&
        "mu" / MU, &
        "gamma" / GAMMA, &
        "units" / GEOMETRICAL)

    fluxes => Dict( &
        "fluxtype"  / KT, &!HLLC, &
        "order" / LINEAR, &
        "variables" / PRIMITIVE, &
        "limiter" / MONOCENT, &
        "theta" / 1.2)

    mesh => Dict( &
        "meshtype" / MIDPOINT, &
        "geometry" / POLAR, &
        "xmin" / RMIN, &
        "xmax" / RMAX, &
!        "geometry" / LOGPOLAR, &
!        "xmin" / LOG(RMIN/1.), &
!        "xmax" / LOG(RMAX/1.), &
        "gparam" / 1., &
        "omega" / omega, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "ymin" / (0.), &
        "ymax" / (2.*PI), &
        "decomposition" / (/ -1, 1/), &
        "output/dl" / 1, &
        "output/radius" / 1, &
        "output/volume" / 1 )

    boundary => Dict( &
        "western" / NOSLIP, &
        "eastern" / NOSLIP, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC)

    grav => Dict( "stype" / GRAVITY, &
        "output/accel" / 1, &
        "self/gtype" / SPECTRAL, &
        "self/green" / green, &
        "self/sigma" / sigma, &
        "self/output/potential" / 1)

    timedisc => Dict( &
        "method" / SSPRK, &
        "rhstype" / 1, &
!        "dtmax" / 8., &
!        "fargo" / 1, &
!        "method" / RK_FEHLBERG, &
!        "order" / 5, &
        "tol_rel" / 1.0E-3, &
        "cfl" / 0.3, &
        "stoptime" / t, &
        "dtlimit" / 1.0E-10, &
        "maxiter" / 2000000000, &
        "output/solution" / 1)

    WRITE(str,'(I1)')green

    datafile => Dict( &
        "fileformat" / XDMF, &
        "filename" / ("orbitingcylinders_" // str), &
        "unit" / 5555, &
        "count" / ONUM)
    
    config => Dict( &
        "physics" / physics, &
        "fluxes" / fluxes, &
        "mesh" / mesh, &
        "boundary" / boundary, &
        "sources/gravity" / grav, &
        "timedisc" / timedisc, &
        "datafile" / datafile)

  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Timedisc,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,k,dir,ig
    REAL              :: cs,s
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: r,r1,r2,r3,y
    TYPE(Sources_TYP), POINTER :: sp
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc,Physics
    !------------------------------------------------------------------------!
    r = Mesh%radius%bcenter
    sp => GetSourcesPointer(Physics%sources,GRAVITY)
    anapot => Timedisc%solution(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1)
    numpot => sp%glist%pot(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1)
    
    SELECT CASE(green)
    CASE(1)
      r1 = rsq(Mesh,Mesh%radius%bcenter,x1)
      r2 = rsq(Mesh,Mesh%radius%bcenter,x3)
      r3 = rsq(Mesh,Mesh%radius%bcenter,x4)
      Timedisc%pvar(:,:,Physics%DENSITY) &
        = 2.0 / (2.*PI*sigma**2) * EXP(-SQRT(r1)/sigma) &
        + 0.5 / (2.*PI*sigma**2) * EXP(-SQRT(r2)/sigma) &
        + 1.0 / (2.*PI*sigma**2) * EXP(-SQRT(r3)/sigma)

      r1 = SQRT(rsq(Mesh,Mesh%radius%faces(:,:,1),x1)) / (2.*sigma)
      r2 = SQRT(rsq(Mesh,Mesh%radius%faces(:,:,1),x3)) / (2.*sigma)
      r3 = SQRT(rsq(Mesh,Mesh%radius%faces(:,:,1),x4)) / (2.*sigma)
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          Timedisc%solution(i,j,1) &
#ifdef HAVE_FGSL
         = -2.0 * r1(i,j)/sigma*(fgsl_sf_bessel_ic0(r1(i,j))*fgsl_sf_bessel_kc1(r1(i,j))&
                  -fgsl_sf_bessel_ic1(r1(i,j))*fgsl_sf_bessel_kc0(r1(i,j))) &
           -0.5 * r2(i,j)/sigma*(fgsl_sf_bessel_ic0(r2(i,j))*fgsl_sf_bessel_kc1(r2(i,j))&
                  -fgsl_sf_bessel_ic1(r2(i,j))*fgsl_sf_bessel_kc0(r2(i,j))) &
           -1.0 * r3(i,j)/sigma*(fgsl_sf_bessel_ic0(r3(i,j))*fgsl_sf_bessel_kc1(r3(i,j))&
                  -fgsl_sf_bessel_ic1(r3(i,j))*fgsl_sf_bessel_kc0(r3(i,j)))
#else
         = -2.0 * r1(i,j)/sigma*(BESSEL_I0(r1(i,j))*BESSEL_K1(r1(i,j))&
                  -BESSEL_I1(r1(i,j))*BESSEL_K0(r1(i,j))) &
           -0.5 * r2(i,j)/sigma*(BESSEL_I0(r2(i,j))*BESSEL_K1(r2(i,j))&
                  -BESSEL_I1(r2(i,j))*BESSEL_K0(r2(i,j))) &
           -1.0 * r3(i,j)/sigma*(BESSEL_I0(r3(i,j))*BESSEL_K1(r3(i,j))&
                  -BESSEL_I1(r3(i,j))*BESSEL_K0(r3(i,j)))
#endif
        END DO
      END DO

      Timedisc%pvar(:,:,Physics%PRESSURE) = 1.0
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
      Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.0
    CASE(2)
      r1 = rsq(Mesh,Mesh%radius%bcenter,x1)
      r2 = rsq(Mesh,Mesh%radius%bcenter,x3)
      r3 = rsq(Mesh,Mesh%radius%bcenter,x4)
      Timedisc%pvar(:,:,Physics%DENSITY) &
        = 2.0 / (2.*PI*sigma**2) * EXP(-r1/(2.*sigma**2)) &
        + 0.5 / (2.*PI*sigma**2) * EXP(-r2/(2.*sigma**2)) &
        + 1.0 / (2.*PI*sigma**2) * EXP(-r3/(2.*sigma**2))

      r1 = rsq(Mesh,Mesh%radius%faces(:,:,1),x1)
      r2 = rsq(Mesh,Mesh%radius%faces(:,:,1),x3)
      r3 = rsq(Mesh,Mesh%radius%faces(:,:,1),x4)

      Timedisc%solution(:,:,1) &
        = - 2.0 / SQRT(r1) * ERF(SQRT(r1 / 2.) / sigma) &
          - 0.5 / SQRT(r2) * ERF(SQRT(r2 / 2.) / sigma) &
          - 1.0 / SQRT(r3) * ERF(SQRT(r3 / 2.) / sigma)
      Timedisc%pvar(:,:,Physics%PRESSURE) = 1.0
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
      Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.0
    CASE(3)
      r1 = rsq(Mesh,Mesh%curv%bcenter,x1)
      r2 = rsq(Mesh,Mesh%curv%bcenter,x2)
      Timedisc%pvar(:,:,Physics%DENSITY)  = 0.02/(3.2*PI) &
        + 0.99 * (  EXP(-0.5*(r1/sigma**2))/(2.*PI*sigma**2) &
                  + EXP(-0.5*(r2/sigma**2))/(2.*PI*sigma**2))

      Timedisc%pvar(:,:,Physics%PRESSURE) = 0.02/(3.2*PI) &
        + GN / (2.*PI*sigma**2) &
           * (  Ei(-(r1/sigma**2)) - Ei(-0.5*(r1/sigma**2)) &
              + Ei(-(r2/sigma**2)) - Ei(-0.5*(r2/sigma**2)))
         
      Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
      Timedisc%pvar(:,:,Physics%YVELOCITY) = r - r*omega
    END SELECT

    DO dir=WEST,EAST
        IF(GetType(Timedisc%Boundary(dir)).EQ.NOSLIP) THEN
            DO j=Mesh%JMIN,Mesh%JMAX
                DO ig=1,Mesh%GNUM
                    SELECT CASE (dir)
                    CASE(WEST)
                        i = Mesh%IMIN-ig
                    CASE(EAST)
                        i = Mesh%IMAX+ig
                    END SELECT
                    Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) &
                      = Timedisc%pvar(i,j,Physics%YVELOCITY)
                END DO
            END DO
        END IF
    END DO

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE InitData

  FUNCTION rsq(Mesh,r,x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP) :: Mesh
    REAL, DIMENSION(2) :: x
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: r
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: res
    !------------------------------------------------------------------------!
    INTENT(IN) :: x
    !------------------------------------------------------------------------!
    res = r**2 + x(1)**2 - 2.*r*x(1) * COS(Mesh%curv%bcenter(:,:,2)-x(2))
  END FUNCTION rsq

END PROGRAM Init
