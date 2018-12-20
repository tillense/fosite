!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: planeteu.f90                                                      #
!#                                                                           #
!# Copyright (C) 2011-2012                                                   #
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
!> For more information about the eu-planet project see:
!! http://ttt.astro.su.se/groups/planets/comparison//index.html
!! The initial conditions are described here:
!!http://ttt.astro.su.se/groups/planets/comparison//description.html
!! New URL: http://www.mps.mpg.de/data/outgoing/deval/comparison/
!!
!! -# GN*(MSUN + MPLANET) = RPLANET = 1
!! -# m = MPLANET/(MSUN + MPLANET)
!!    m = 10E-3 => Jupiter
!!    m = 10E-4 => Neptun
!! -# a) No viscosity
!!    b) Physical viscosity of 1.0e-5 in these units
!!       (=> alpha-vis = 6e-3..2.5e-3)
!! -# unsoftened gravity from the star
!! -# smoothed potential from the planet:
!!    Phi_planet = - MPLANET/SQRT(r**2 + epsilon**2)
!!    epsilon = 0.6*H (H: disc semi-thickness)
!! -# R = [0.4, 2.5] * RPLANET
!! -# RES = [128, 384] => planet at a cell center
!! -# Rotating or non-rotating reference frame
!! -# inner and outer boundaries should be reflecting
!! -# Wave killing zones..?
!! -# initial conditions:
!!    a) uniform density sigma: sigma * pi * RPLANET**2 / MSUN = 0.002
!!    b) Keplerian rotation
!!    c) introduce planet smoothly sin**2(pi * t / (10 * P_orbit))
!!       => planet reaches its full mass over five orbital rotations
!! -# required output after 2,5,10,20,50,100,200,500 orbital periods
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
  USE common_dict
  USE sources_rotframe, ONLY : Convert2RotatingFrame_rotframe
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  ! problem setup
!  REAL, PARAMETER    :: m = 1.E-10
  REAL, PARAMETER    :: m = 1.0E-3  ! Jupiter
!  REAL, PARAMETER    :: m = 1.0E-4  ! Neptun
  LOGICAL, PARAMETER :: enable_vis = .FALSE.

  ! general constants
  REAL, PARAMETER :: GN          = 1.0
  REAL, PARAMETER :: mcentral    = 1.0 - m
  REAL, PARAMETER :: a           = 1.0
  REAL, PARAMETER :: Sigma0      = 2.0E-3 * mcentral / (PI*a*a)
  REAL, PARAMETER :: omega       = 1.0 ! == SQRT(GN*mcentral/(a*a*a))

  REAL, PARAMETER :: P0          = 2.0*PI/omega
  REAL, PARAMETER :: P           = 10
  REAL, PARAMETER :: t           = P*P0

  REAL, PARAMETER :: flaring     = 0.05
!  REAL, PARAMETER :: csiso       = flaring * 1.0 ! 1.0 == SQRT(GN*mcentral/a)
  REAL, PARAMETER :: csiso       = 0.
  REAL, PARAMETER :: RG          = 8.31
  REAL, PARAMETER :: MU          = 6.02E-04
  REAL, PARAMETER :: GAMMA       = 5./3.
  REAL, PARAMETER :: eps         = 0.6 * (flaring * a)

  REAL, PARAMETER :: RMIN        = 0.4 * a
  REAL, PARAMETER :: RMAX        = 2.5 * a
  REAL, PARAMETER :: GPAR        = 1.0

  INTEGER, PARAMETER :: XRES     = 128
  INTEGER, PARAMETER :: YRES     = 3*XRES
  INTEGER, PARAMETER :: ONUM     = P * 1

  REAL            :: hillradius 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)          :: Sim
  !--------------------------------------------------------------------------!
  
hillradius = a*(m/(3.0*mcentral))**(1.0/3.0)

CALL InitFosite(Sim)

CALL MakeConfig(Sim%config)

!CALL PrintDict(Sim%config)

CALL SetupFosite(Sim)

! set initial condition
CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics)

CALL RunFosite(Sim)

CALL CloseFosite(Sim)

CONTAINS
  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,sources,&
                              fluxes,grav,pcentral,psat,vis,physics,rotframe,&
                              damping,binary
    REAL    :: a1,a2,a3
    !------------------------------------------------------------------------!
    physics => Dict( &
!        "problem" / EULER2D_ISOTHERM, &
        "problem" / EULER2D_ISOIAMT, &
        "cs" / csiso, &
!        "problem" / EULER2D_IAMT, &
        "mu" / MU, &
        "gamma" / GAMMA, &
        "dpmax" / 1.0E+1, &
        "pmin" / 1.0E-30, &
        "units" / GEOMETRICAL)

    fluxes => Dict( &
        "fluxtype"  / KT, &
        "order" / LINEAR, &
        "variables" / PRIMITIVE, &
!        "limiter" / MINMOD, &
!        "limiter" / SWEBY, &
!        "limiter" / MONOCENT, &
!        "limiter" / OSPRE, &
        "limiter" / VANLEER, &
        "theta" / 1.2)

    a1 = (1.3-a)/SINH(1.0)
    a2 = a
    a3 = 5.0

    mesh => Dict( &
        "meshtype" / MIDPOINT, &
        "geometry" / POLAR, &
        "omega" / omega, &
        "xmin" / RMIN, &
        "xmax" / RMAX, &
!        "geometry" / LOGPOLAR, &
!        "xmin" / LOG(RMIN/0.1), &
!        "xmax" / LOG(RMAX/0.1), &
!        "gparam" / 0.1, &
!        "geometry" / SINHPOLAR, &
!        "xmin" / (ASINH((RMIN-a2)/a1)+a3), &
!        "xmax" / (ASINH((RMAX-a2)/a1)+a3), &
!        "gparam" / a1, &
!        "gparam2" / a2, &
!        "gparam3" / a3, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "ymin" / (PI*(-1.0+1.0/XRES)), &
        "ymax" / (PI*( 1.0+1.0/XRES)), &
        "output/volume" / 1 )

    boundary => Dict( &
!        "western" / REFLECTING, &
!        "eastern" / REFLECTING, &
!        "western" / NO_GRADIENTS, &
!        "eastern" / NO_GRADIENTS, &
        "western" / NOSLIP, &
        "eastern" / NOSLIP, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC)

    pcentral => Dict( &
        "gtype" / POINTMASS, & 
        "potential" / NEWTON, &
        "mass" / mcentral, &
        "cvis" / 0.9, &
        "outbound" / 0 )

    psat => Dict( &
        "gtype" / POINTMASS, &
        "potential" / NEWTON, &
        "mass" / m, &
        "cvis" / 0.9, &
        "x" / a, &
        "y" / 0.0, &
        "outbound" / 0, &
!        "hillradius" / (1.0E-6 * hillradius), &
        "switchon" / (5.0 * P0), &
        "softening" / eps)

    binary => Dict( &
        "gtype" / POINTMASS_BINARY, &
        "mass1" / m, &
        "mass2" / mcentral, &
        "semimayoraxis" / a, &
        "excentricity" / 0.0, &
        "cvis"  / 0.5, &
        "softening1" / eps)

    grav => Dict("stype" / GRAVITY, &
            "cvis" / 0.5, &
            "output/accel" / 1, &
           "pcentral" / pcentral, &
           "psat" / psat)!, &
!           "binary" / binary)

    rotframe => Dict( &
        "stype" / ROTATING_FRAME)

    damping => Dict( &
        "stype" / WAVE_DAMPING, &
        "r0"    / (0.5*a), &
        "r1"    / (2.1*a), &
        "tau0"  / (2.0*PI / SQRT(GN*mcentral/(0.4*a)**3)), &
        "tau1"  / (2.0*PI / SQRT(GN*mcentral/(2.5*a)**3)))
    
    sources => Dict( &
        "grav" / grav, &
!        "rotframe" / rotframe &
        "damping" / damping &
        )

    IF(enable_vis) THEN
        vis => Dict( &
            "stype" / VISCOSITY, &
            "vismodel" / PRINGLE, &
            "dynconst" / 1.0E-5, &
            "bulkconst" / 0.0, &
            "cvis" / 0.3)
        CALL SetAttr(sources, "vis", vis)
    END IF

    timedisc => Dict( &
!        "method" / MODIFIED_EULER, &
        "method" / SSPRK, &
!        "method" / RK_FEHLBERG, &
!        "order" / 5, &
!        "tol_rel" / 1.0E-8, &
        "fargo" / 1, &
        "cfl" / 0.75, &
        "stoptime" / t, &
        "dtlimit" / 1.0E-10, &
        "maxiter" / 2000000000, &
        "output/xmomentum" / 1, &
        "output/ymomentum" / 1, &
!        "output/energy" / 1, &
        "output/geometrical_sources" / 1, &
        "output/external_sources" / 1, &
        "output/bflux" / 1)

    datafile => Dict( &
        "fileformat" / XDMF, &
        "filename" / "planeteu", &
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


  SUBROUTINE InitData(Timedisc,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,k,dir,ig
    REAL              :: r,p
    TYPE(Sources_TYP), POINTER :: sp
    REAL              :: cs
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc,Physics
    !------------------------------------------------------------------------!

    Timedisc%pvar(:,:,Physics%DENSITY) = Sigma0
    ! Since the initial density is constant everywhere, there is no pressure
    ! gradient force contributing to the initial YVELOCITY
    Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.0
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.0

    ! Account for pressure support by scaling the azimuthal velocity
    IF(csiso.LE.0.) THEN
        p = (1. - flaring**2)**0.5
    ELSE
        p = 1.
    END IF

    DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
            r = Mesh%radius%bcenter(i,j)
            Timedisc%pvar(i,j,Physics%YVELOCITY) = SQRT(GN*mcentral/r)*p &
                                                   - r*omega
        END DO
    END DO

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
                    r = Mesh%radius%bcenter(i,j)
                    Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) &
                                = SQRT(GN*mcentral/r)*p-r*omega
                END DO
            END DO
        END IF
    END DO


    IF(csiso.LE.0.) THEN
      CALL SetSoundSpeeds(Physics,Mesh,flaring * SQRT(GN*mcentral/Mesh%radius%bcenter(:,:)))
      CALL SetSoundSpeeds(Physics,Mesh,flaring * SQRT(GN*mcentral/Mesh%radius%faces(:,:,:)))
      IF(Physics%PRESSURE.GT.0) &
        Timedisc%pvar(:,:,Physics%PRESSURE) &
          = Timedisc%pvar(:,:,Physics%DENSITY) * flaring**2 &
            * GN*mcentral/Mesh%radius%bcenter(:,:)
    END IF

!    sp => GetSourcesPointer(Physics%sources,ROTATING_FRAME)
!    IF(ASSOCIATED(sp)) &
!        CALL Convert2RotatingFrame_rotframe(sp, Mesh, Physics, Timedisc%pvar)

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
  END SUBROUTINE InitData
END PROGRAM Init
