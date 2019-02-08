!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: keplerianvortex.f90                                               #
!#                                                                           #
!# Copyright (C) 2014-2018                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!> Vortex transport in a keplerian disk \n
!!
!! \example keplerianvortex.f90
!!
!! References:
!! - \cite mignone2012 Mignone, A. et al.: "A conservative orbital advection scheme for
!!     simulations of magnetized shear flows with the PLUTO code"
!!     Astronomy & Astrophysics, Volume 545, id.A152, 16 pp.
!! - \cite bodo2007 Bodo, G. et al.: "Stability amd nonlinear adjustment of vortices in
!!     Keplerian flows"
!!     A & A 475, 51-61 (2007)
!----------------------------------------------------------------------------!
PROGRAM keplerianvortex
  USE fosite_mod
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

  INTEGER, PARAMETER :: GravEnergySource = 1
  INTEGER, PARAMETER :: XRES     = 64
  INTEGER, PARAMETER :: YRES     = 4*XRES
  INTEGER, PARAMETER :: ONUM     = P * 10
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE          :: Sim
  !--------------------------------------------------------------------------!

ALLOCATE(Sim)
CALL Sim%InitFosite()
CALL MakeConfig(Sim%config)
CALL Sim%Setup()
CALL InitData(Sim%Timedisc, Sim%Mesh, Sim%Physics, Sim%Fluxes)
CALL Sim%Run()
CALL Sim%Finalize()
DEALLOCATE(Sim)

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
        "problem" / EULER, &
        "gamma" / GAMMA, &
        "units" / GEOMETRICAL)

    fluxes => Dict( &
        "fluxtype" / KT, &
        "order" / LINEAR, &
        "variables" / PRIMITIVE, &
        "limiter" / VANLEER)

    mesh => Dict( &
        "meshtype" / MIDPOINT, &
        "geometry" / CYLINDRICAL, &
        "omega" / omega, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "knum" / 1, &
        "xmin" / RMIN, &
        "xmax" / RMAX, &
        "ymin" / (-PI), &
        "ymax" / (PI), &
        "zmin" / 0.0, &
        "zmax" / 0.0, &
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

!    damping => Dict( &
!        "stype" / WAVE_DAMPING, &
!        "zone0/rmin" / RMIN, &
!        "zone0/rmax" / (0.5), &
!        "zone0/tau"  / (0.1*2.0*PI / SQRT(GN*MCENTRAL/(0.4)**3)), &
!        "zone1/rmin" / (1.9), &
!        "zone1/rmax" / RMAX, &
!        "zone1/tau"  / (0.1*2.0*PI / SQRT(GN*MCENTRAL/(2.0)**3)))

    sources => Dict( &
!        "damping" / damping, &
        "grav/stype" / GRAVITY, &
!        "grav/energy" / GravEnergySource, &
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
!        "rhstype" / 1, &
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
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                      :: dvr, dvphi, x, y, phi, r
    REAL              :: h, kappa
    INTEGER           :: dir,j,ig,i,k
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp
    !------------------------------------------------------------------------!
    x = Mesh%bccart(:,:,:,1)
    y = Mesh%bccart(:,:,:,2)

    IF(GravEnergySource.EQ.0) THEN
      CALL SetPotential_euler2Dia(Physics,Mesh,-GN * MCENTRAL / Mesh%radius%bcenter)
      CALL SetPotential_euler2Dia(Physics,Mesh,-GN * MCENTRAL / Mesh%radius%faces)
    END IF

    ! get gravitational acceleration
    sp => Sim%Sources
    DO
      IF (ASSOCIATED(sp).EQV..FALSE.) RETURN
      SELECT TYPE(sp)
      CLASS IS(sources_gravity)
        gp => sp
        EXIT
      END SELECT
      sp => sp%next
    END DO

    SELECT TYPE (pvar => Timedisc%pvar)
    CLASS IS(statevector_euler)
      SELECT TYPE (cvar => Timedisc%cvar)
      CLASS IS(statevector_euler)
        cvar%density%data3d(:,:,:) = Sigma0
        cvar%momentum%data4d(:,:,:,1:Physics%VDIM) = 0.

        CALL Physics%Convert2Primitive(cvar,pvar)

        pvar%pressure%data3d(:,:,:) = 1./(GAMMA*M**2)

!        CALL gp%UpdateGravity(Mesh,Physics,Fluxes,pvar,0.0,0.0)

        pvar%velocity%data4d(:,:,:,1:Physics%VDIM) = &
          pvar%velocity%data4d(:,:,:,1:Physics%VDIM) &
          + Timedisc%GetCentrifugalVelocity(Mesh,Physics,Sim%Fluxes,Sim%Sources,(/0.,0.,1./))!,gp%accel%data4d)

      END SELECT
    END SELECT

    ! boundary conditions
    ! custom boundary conditions at western boundary if requested
    SELECT TYPE(bwest => Timedisc%Boundary%boundary(WEST)%p)
    CLASS IS (boundary_custom)
      CALL bwest%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_REFLECT,CUSTOM_REFLNEG,CUSTOM_KEPLER,CUSTOM_REFLECT/))
    END SELECT
    SELECT TYPE(beast => Timedisc%Boundary%boundary(EAST)%p)
    CLASS IS (boundary_custom)
      CALL beast%SetCustomBoundaries(Mesh,Physics, &
        (/CUSTOM_REFLECT,CUSTOM_REFLNEG,CUSTOM_KEPLER,CUSTOM_REFLECT/))
    END SELECT

    h = 1./(2.*M)
    kappa = -1.
    phi = Mesh%center(:,:,:,2)
    x = x - 1./SQRT(2.)
    y = y - 1./SQRT(2.)
    r = SQRT(x**2 + y**2)

    dvr   = kappa * EXP(-r**2/h**2) * (-y*COS(phi) + x*SIN(phi))
    dvphi = kappa * EXP(-r**2/h**2) * ( y*SIN(phi) + x*COS(phi))

    SELECT TYPE (pvar => Timedisc%pvar)
    CLASS IS(statevector_euler)
      SELECT TYPE (cvar => Timedisc%cvar)
      CLASS IS(statevector_euler)
        pvar%velocity%data4d(:,:,:,1) = pvar%velocity%data4d(:,:,:,1) + dvr
        pvar%velocity%data4d(:,:,:,2) = pvar%velocity%data4d(:,:,:,2) + dvphi
      END SELECT
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

  END SUBROUTINE InitData
END PROGRAM keplerianvortex
