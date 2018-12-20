!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: keplerianvortex.f90                                               #
!#                                                                           #
!# Copyright (C) 2014                                                        #
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
! Vortex transport in a keplerian disk
! [1] Mignone, A. et al.: "A conservative orbital advection scheme for
!     simulations of magnetized shear flows with the PLUTO code"
!     Astronomy & Astrophysics, Volume 545, id.A152, 16 pp.
! [2] Bodo, G. et al.: "Stability amd nonlinear adjustment of vortices in
!     Keplerian flows"
!     A & A 475, 51-61 (2007)
!----------------------------------------------------------------------------!
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
  USE physics_euler2Diamt, ONLY: SetPotential_euler2Dia
  USE common_dict
  USE timedisc_generic, ONLY: GetCentrifugalVelocity
  IMPLICIT NONE
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER :: GN          = 1.0
  REAL, PARAMETER :: MCENTRAL    = 1.0
  REAL, PARAMETER :: Sigma0      = 1.0
  REAL, PARAMETER :: M           = 10.
  REAL, PARAMETER :: omega       = 0.

  REAL, PARAMETER :: P0          = 2.0*PI/1.0
  REAL, PARAMETER :: P           = 3
  REAL, PARAMETER :: t           = P*P0

  REAL, PARAMETER :: GAMMA       = 5./3.

  REAL, PARAMETER :: RMIN        = 0.5
  REAL, PARAMETER :: RMAX        = 1.5

  INTEGER, PARAMETER :: GravEnergySource = 0
  INTEGER, PARAMETER :: XRES     = 128
  INTEGER, PARAMETER :: YRES     = 4*XRES
  INTEGER, PARAMETER :: ONUM     = P * 10
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)          :: Sim
  !--------------------------------------------------------------------------!

CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(Sim%config)

CALL SetupFosite(Sim)

! set initial condition
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics,Sim%Fluxes)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              sources,fluxes,pmass,physics,damping,vis,pot
    !------------------------------------------------------------------------!
    physics => Dict( &
        "problem" / EULER2D_IAMT, &
        "gamma" / GAMMA, &
        "units" / GEOMETRICAL)

    fluxes => Dict( &
        "fluxtype" / HLLC, &
        "order" / LINEAR, &
        "variables" / PRIMITIVE, &
        "limiter" / VANLEER)

    mesh => Dict( &
        "meshtype" / MIDPOINT, &
        "geometry" / POLAR, &
        "omega" / omega, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "xmin" / RMIN, &
        "xmax" / RMAX, &
        "ymin" / (-PI), &
        "ymax" / (PI), &
        "decomposition" / (/ 1, -1/), &
        "output/volume" / 1,&
        "output/bh" / 1,&
        "output/rotation" / 1,&
        "output/dl" / 1)

    boundary => Dict( &
        "western" / CUSTOM, &
        "eastern" / CUSTOM, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC)

    pmass => Dict( &
        "gtype" / POINTMASS, & 
        "potential" / NEWTON, &
        "mass" / MCENTRAL, &
        "cvis" / 0.9, &
        "outbound" / 0)!, &
!        "output/accel" / 1 )

    damping => Dict( &
        "stype" / WAVE_DAMPING, &
        "r0"    / (0.5), &
        "r1"    / (1.9), &
        "tau0"  / (0.1*2.0*PI / SQRT(GN*MCENTRAL/(0.4)**3)), &
        "tau1"  / (0.1*2.0*PI / SQRT(GN*MCENTRAL/(2.0)**3)))

    sources => Dict( &
!        "damping" / damping, &
        "grav/stype" / GRAVITY, &
        "grav/energy" / GravEnergySource, &
        "grav/output/accel" / 1, &
        "grav/pmass"   / pmass)

    timedisc => Dict( &
        "method" / SSPRK, &
!         "fargo" / 1, &
        "tol_rel" / 1.0E-3, &
        "tol_abs" / (/ 0., 1.E-3, 1.E-3, 0. /), &
        "cfl" / 0.4, &
        "stoptime" / t, &
        "dtlimit" / 1.0E-10, &
        "maxiter" / 2000000000, &
        "output/xmomentum" / 1, &
        "output/ymomentum" / 1, &
        "output/energy" / 1, &
        "output/geometrical_sources" / 1, &
        "output/external_sources" / 1, &
        !"output/bflux" / 1, &
        "output/fluxes" / 1, &
        "output/rhs" / 1)


    datafile => Dict( &
        "fileformat" / XDMF, &
        "filename" / "keplerianvortex", &
        "count" / ONUM)

    config => Dict( &
        "physics" / physics, &
        "fluxes" / fluxes, &
        "mesh" / mesh, &
        "boundary" / boundary, &
        "sources" / sources, &
        "timedisc" / timedisc, &
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: dvr, dvphi, x, y, phi, r
    REAL              :: h, kappa
    INTEGER           :: dir,j,ig,i,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc, Physics
    !------------------------------------------------------------------------!
    x = Mesh%bccart(:,:,1)
    y = Mesh%bccart(:,:,2)

    IF(GravEnergySource.EQ.0) THEN
      CALL SetPotential_euler2Dia(Physics,Mesh,-GN * MCENTRAL / Mesh%radius%bcenter)
      CALL SetPotential_euler2Dia(Physics,Mesh,-GN * MCENTRAL / Mesh%radius%faces)
    END IF

    Timedisc%cvar(:,:,Physics%DENSITY)   = Sigma0
    Timedisc%cvar(:,:,Physics%XMOMENTUM:Physics%YMOMENTUM) = 0.

    CALL Convert2Primitive(Physics,Mesh,Timedisc%cvar,Timedisc%pvar)

    IF(Physics%PRESSURE.GT.0) &
      Timedisc%pvar(:,:,Physics%PRESSURE)  = 1./(GAMMA*M**2)

    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = &
      Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) &
      + GetCentrifugalVelocity(Timedisc,&
          Mesh,Physics,Sim%Fluxes,(/0.,0.,1./))

    DO dir=WEST,EAST
      SELECT CASE(GetType(Timedisc%Boundary(dir)))
      CASE(CUSTOM)
        Timedisc%Boundary(dir)%cbtype(:,:) = CUSTOM_FIXED
        DO j=Mesh%JMIN,Mesh%JMAX
          DO ig=1,Mesh%GNUM
            SELECT CASE (dir)
            CASE(WEST)
              i = Mesh%IMIN-ig
            CASE(EAST)
              i = Mesh%IMAX+ig
            END SELECT
            Timedisc%Boundary(dir)%data(ig,j,:) &
              = Timedisc%pvar(i,j,:)
          END DO
        END DO
          Timedisc%boundary(dir)%cbtype(:,Physics%DENSITY) = CUSTOM_REFLECT
          Timedisc%boundary(dir)%cbtype(:,Physics%XVELOCITY) = CUSTOM_REFLNEG
          Timedisc%boundary(dir)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
            Timedisc%Boundary(dir)%cbtype(:,Physics%PRESSURE) = CUSTOM_REFLECT
      END SELECT
    END DO


    h = 1./(2.*M)
    kappa = -1.
    phi = Mesh%center(:,:,2)
    x = x - 1./SQRT(2.)
    y = y - 1./SQRT(2.)
    r = SQRT(x**2 + y**2)

    dvr   = kappa * EXP(-r**2/h**2) * (-y*COS(phi) + x*SIN(phi))
    dvphi = kappa * EXP(-r**2/h**2) * ( y*SIN(phi) + x*COS(phi))

    Timedisc%pvar(:,:,Physics%XVELOCITY) &
      = Timedisc%pvar(:,:,Physics%XVELOCITY) + dvr
    Timedisc%pvar(:,:,Physics%YVELOCITY) &
      = Timedisc%pvar(:,:,Physics%YVELOCITY) + dvphi

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

  END SUBROUTINE InitData
END PROGRAM Init
