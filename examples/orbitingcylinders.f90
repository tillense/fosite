!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: orbitingcylinders.f03                                             #
!#                                                                           #
!# Copyright (C) 2013-2018                                                   #
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
!> self-gravitating orbiting cylinder test for spectral solver
!!
!! \author Manuel Jung
!!
!! \example orbitingcylinders.f03
!!
!! This part is adapted (and translated) from the PhD-thesis from Jung (2016,
!! \cite jung2016 ). It contains two different types of tests with cylinders.
!! Both originally come from Chan (2006, \cite chan2006 ):
!! 1. A static simulation of the gravitational potential in a
!!    - razor-thin disk.
!!    - gaussian distribution in vertical direction.
!! 2. A dynamic evolution of the cylinders.
!!
!! \section static_spectral_test Density-Potential Pairs
!!
!! This part compares the solution of the potential by the gravitational
!! \link gravity_spectral spectral solver \endlink with an existing analytical
!! result. Since the module only calculates the potential within the simulation
!! area, the potential produced by the boundaries need to be neglectable.
!! However, the potential from outside the potential can calculated with
!! additional gravitational modules.
!!
!! Assuming a density distributions around a central point
!! \f[
!!    \Sigma_{(r_k,\varphi_k)}(r,\varphi) = \frac{1}{2 \pi \sigma^2}
!!    \exp{\left( - \frac{R_k^2}{\sigma} \right)},
!! \f]
!! with \f$ \sigma \f$ a measurement for the width of and \f$ R_k \f$
!! the distance of a coordinate point form the center of the distribution.
!! For a razor thin disk, this yields the potential (see \cite chan2006):
!! \f[
!!    \Phi_{(r_k, \varphi_k)} = - \frac{G}{\sigma} \left(I_0(y_k)K_1(y_k) -
!!    I_1(y_k)K_0(y_k)\right),
!! \f]
!! where \f$ y_k = \frac{R_k}{2 \sigma} \f$. \f$ I_n \f$ and \f$ K_n \f$ are
!! the modified Bessel functions of first and second kind with integral order
!! (see \link functions functions \endlink).
!!
!! The density distribution in our example is
!! \f[
!!    \Sigma(r,\varphi) = 2 \Sigma_{(1,10^{-3})}(r,\varphi) + 0.5
!!      \Sigma_{(1,\pi+10^{-3})} (r,\varphi) + \Sigma_{(0.9,\frac{3}{4}\pi)}
!!      (r,\varphi),
!! \f]
!! with the according potential
!! \f[
!!    \Phi(r,\varphi) = 2 \Phi_{(1,10^{-3})}(r,\varphi) + 0.5
!!      \Phi_{(1,\pi+10^{-3})} (r,\varphi) + \Phi_{(0.9,\frac{3}{4}\pi)}
!!      (r,\varphi).
!! \f]
!!
!! \todo Make distinction between these two tests clear. What is with
!! \f$ \sigma \f$.
!!
!! Additionally, the implementation is tested for a vertical gaussian
!! distribution, with a scale height similar to \f$ \sigma =0.1 \f$. with the
!! same setup as above. This gives us a solution for the potential
!! \f[
!!    \Phi(r_k, \varphi_k) = - \frac{1}{R_k} \mathrm{erf}\left(
!!        \frac{R_k}{\sqrt{2} \sigma} \right).
!! \f]
!!
!! Below the relative error for the potential is shown,
!! solved for a razor-thin disk (left) and a vertical gaussian distribution
!! (right). For the first case we see a maximum error of \f$ 10^{-3} \f$
!! at steep gradients and \f$ 10^{-5} \f$ for the second case.
!!
!! \image html http://www.astrophysik.uni-kiel.de/fosite/orbitingcylinders_relativeerror.png "relative error of numerical to potential" \n
!!
!!
!!
!! \section rotating_cylinders Self-Gravitating Rotating Cylinders
!!
!! A time-dependent test from \cite chan2006, which generates two cylinders
!! with gaussian density distribution in polar grid placed opposite towards
!! each other. A background field adds a rigid rotation. The thereby arising
!! centrifugal force has to be balanced from the self-gravitation of the
!! mass-distribution. This is a very sophisticated setup for the accuracy of
!! the acceleration by self-gravitation and the general angular-momentum
!! conservation.
!! The setup is
!! \f[
!!    \varrho(r,\varphi) = \frac{10^{-2}}{\pi \left( r_{\mathrm{max}}^2 -
!!    r_{\mathrm{min}}^2 \right)} + 0.99 \left( \frac{1}{2 \pi \sigma^2}
!!    \exp{\left(- \frac{R_1}{2 \sigma^2} \right)} + \frac{1}{2 \pi \sigma}
!!    \exp{\left(- \frac{R_2}{2 \sigma^2} \right)}\right) \\
!!    p(r,\varphi) = \frac{10^{-2}}{\pi \left( r_{\mathrm{max}}^2 -
!!    r_{\mathrm{min}}^2 \right)} + \frac{G}{2 \pi \sigma^2}
!!    \left(\mathrm{Ei}\left(- \frac{R_1}{\sigma^2} \right) -
!!    \mathrm{Ei}\left(- \frac{R_1}{2 \sigma^2} \right) +
!!    \mathrm{Ei}\left(- \frac{R_2}{\sigma^2} \right) -
!!    \mathrm{Ei}\left(- \frac{R_2}{2 \sigma^2} \right)\right),
!! \f]
!! and \f$ \mathbf{v} = r \mathbf{e}_{\varphi} \f$.
!! The function \f$ \mathrm{Ei}(x) := \int_{-\infty}^x \frac{\exp{t}}{2} \mathrm{d}t \f$
!! is the exponential integral, \f$ R_1 \f$, \f$ R_2 \f$ are the radii from
!! the center of the cylinders and \f$ \sigma = 0.1 \f$ a measure for
!! spreading of the cylinders.
!!
!!  | Additional Settings:                                         |
!!  | ------------------                                           |
!!  | angular velocity: \f$ \Omega = 1.0 \f$                       |
!!  | heat capacity ratio: \f$ \gamma = \frac{5}{3} \f$            |
!!  | resolution: \f$ N_r \times N_{\varphi} = 256 \times 1024 \f$ |
!!  | simulation time: \f$ t_{\mathrm{sim}} = 100 \f$              |
!!  | gravitational constant: \f$ G = 1.0 \f$                      |
!!  | \f$ r_{\mathrm{min}}, r_{\mathrm{max}} = 0.2, 1.8 \f$        |
!!
!! The images below show the results at \f$ t = 0, 100 \f$. The cylinders retain
!! their shape very well. Within \f$ 100 \f$ orbits of the cylinders in the
!! rotating frame there is a loss of \f$ 0.06 \% \f$ of total angular
!! momentum. Compared with \cite chan2006, we have a \f$ 200 \times \f$ better
!! angular momentum conservation in one time step.
!!
!! \image html http://www.astrophysik.uni-kiel.de/fosite/orbitingcylinders.png "colormesh of density" \n
!! \image html http://www.astrophysik.uni-kiel.de/fosite/orbitingcylinders_cut.png "azimuthal cut through density field" \n
!! \image html http://www.astrophysik.uni-kiel.de/fosite/orbitingcylinders_pressure.png "colormesh of the pressure" \n
!!
!!
!! References:
!! - \cite chan2006 Chan, Chi-kwan; Psaltis, Dimitrios; Özel, Feryal, 2006
!! "Spectral Methods for Time-dependent Studies of Accretion Flows. II.
!! Two-dimensional Hydrodynamic Disks with Self-Gravity"
!! http://adsabs.harvard.edu/abs/2006ApJ...645..506C
!! - \cite jung2016 Jung, Manuel. “Multiskalensimulationen von
!!    Schwarzlochakkretion in Balkengalaxien.” PhD-thesis,
!!    Christian-Albrechts Universität Kiel, CAU Kiel, 2016.
!!    http://macau.uni-kiel.de/receive/dissertation_diss_00019153.
!----------------------------------------------------------------------------!
!#define HAVE_FGSL
PROGRAM orbitingcylinders
  USE fosite_mod
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
  INTEGER, PARAMETER :: YRES     = 128
  INTEGER, PARAMETER :: XRES     = YRES/4
  INTEGER, PARAMETER :: ONUM     = P * 1
  REAL, PARAMETER    :: omega    = 1.0
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  INTEGER :: green
  REAL, DIMENSION(:,:,:), POINTER :: numpot, anapot
  !--------------------------------------------------------------------------!

  ALLOCATE(Sim)
  DO green=3,3


    CALL Sim%InitFosite()
    CALL MakeConfig(Sim, Sim%config)

    SELECT CASE(green)
    CASE(3)
      CALL SetAttr(Sim%config, "/mesh/rmin", 0.2)
      CALL SetAttr(Sim%config, "/mesh/rmax", 1.8)
      CALL SetAttr(Sim%config, "/mesh/inum", 64)
      CALL SetAttr(Sim%config, "/mesh/jnum", 192)
    END SELECT

    CALL Sim%Setup()
    CALL InitData(Sim%Timedisc,Sim%Mesh,Sim%Physics,Sim%Sources)

    SELECT CASE(green)
    CASE(1,2)
      CALL Sim%FirstStep()
      print *,"max local relative error: ",MAXVAL(ABS((numpot-anapot)/anapot))
      print *,"global relative error: ",SUM(ABS(numpot-anapot))/SUM(ABS(anapot))
    CASE(3)
      CALL Sim%Run()
      CALL Sim%Finalize()
    END SELECT
  END DO
  DEALLOCATE(Sim)


CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP),POINTER :: mesh,boundary,timedisc,datafile,&
                              fluxes,physics,grav,rotframe,sources
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
        "geometry" / CYLINDRICAL, &
        "xmin" / RMIN, &
        "xmax" / RMAX, &
!        "geometry" / LOGPOLAR, &
!        "xmin" / LOG(RMIN/1.), &
!        "xmax" / LOG(RMAX/1.), &
        "gparam" / 1., &
        "omega" / OMEGA, &
        "inum" / XRES, &
        "jnum" / YRES, &
        "knum" / 1, &
        "ymin" / (0.), &
        "ymax" / (2.*PI), &
        "zmin" / (0.), &
        "zmax" / (0.), &
        "decomposition" / (/ -1, 1/))!, &
!        "output/dl" / 1, &
!        "output/radius" / 1, &
!        "output/volume" / 1 )

    boundary => Dict( &
        "western" / NOSLIP, &
        "eastern" / NOSLIP, &
        "southern" / PERIODIC, &
        "northern" / PERIODIC, &
        "bottomer" / REFLECTING, &
        "topper" / REFLECTING)

    grav => Dict( &
        "stype" / GRAVITY, &
        "output/accel" / 1, &
        "self/gtype" / SPECTRAL, &
        "self/green" / green, &
        "self/sigma" / sigma)
!        "self/output/potential" / 1)

    rotframe => Dict( &
               "stype" / ROTATING_FRAME, &
               "x"     / 0.0, &
               "y"     / 0.0, &
               "z"     / 0.0)

    sources => Dict( &
!              "rotframe"        / rotframe, &
              "grav"            / grav)

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
        "fileformat" / VTK, &
        "filename" / ("orbitingcylinders_" // str), &
        "unit" / 5555, &
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


  SUBROUTINE InitData(Timedisc,Mesh,Physics,Sources)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: Timedisc
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_base),  POINTER       :: Sources
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,k,dir,ig
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: r,r1,r2,r3
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp
    !------------------------------------------------------------------------!
    r = Mesh%radius%bcenter
    sp => Sources
    DO
      IF (ASSOCIATED(sp).EQV..FALSE.) RETURN
      SELECT TYPE(sp)
      CLASS IS(sources_gravity)
        gp => sp
        EXIT
      END SELECT
      sp => sp%next
    END DO


    anapot => Timedisc%solution%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1)
    numpot => gp%glist%pot(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1)

    SELECT CASE(green)
    CASE(1)
      r1 = rsq(Mesh,Mesh%radius%bcenter,x1)
      r2 = rsq(Mesh,Mesh%radius%bcenter,x3)
      r3 = rsq(Mesh,Mesh%radius%bcenter,x4)
      Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) &
        = 2.0 / (2.*PI*sigma**2) * EXP(-SQRT(r1)/sigma) &
        + 0.5 / (2.*PI*sigma**2) * EXP(-SQRT(r2)/sigma) &
        + 1.0 / (2.*PI*sigma**2) * EXP(-SQRT(r3)/sigma)

      r1 = SQRT(rsq(Mesh,Mesh%radius%faces(:,:,:,1),x1)) / (2.*sigma)
      r2 = SQRT(rsq(Mesh,Mesh%radius%faces(:,:,:,1),x3)) / (2.*sigma)
      r3 = SQRT(rsq(Mesh,Mesh%radius%faces(:,:,:,1),x4)) / (2.*sigma)
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            Timedisc%solution%data4d(i,j,k,1) &
#ifdef HAVE_FGSL
           = -2.0 * r1(i,j,k)/sigma*(fgsl_sf_bessel_ic0(r1(i,j,k))*fgsl_sf_bessel_kc1(r1(i,j,k))&
                    -fgsl_sf_bessel_ic1(r1(i,j,k))*fgsl_sf_bessel_kc0(r1(i,j,k))) &
             -0.5 * r2(i,j,k)/sigma*(fgsl_sf_bessel_ic0(r2(i,j,k))*fgsl_sf_bessel_kc1(r2(i,j,k))&
                    -fgsl_sf_bessel_ic1(r2(i,j,k))*fgsl_sf_bessel_kc0(r2(i,j,k))) &
             -1.0 * r3(i,j,k)/sigma*(fgsl_sf_bessel_ic0(r3(i,j,k))*fgsl_sf_bessel_kc1(r3(i,j,k))&
                    -fgsl_sf_bessel_ic1(r3(i,j,k))*fgsl_sf_bessel_kc0(r3(i,j,k)))
#else
           = -2.0 * r1(i,j,k)/sigma*(BESSEL_I0(r1(i,j,k))*BESSEL_K1(r1(i,j,k))&
                    -BESSEL_I1(r1(i,j,k))*BESSEL_K0(r1(i,j,k))) &
             -0.5 * r2(i,j,k)/sigma*(BESSEL_I0(r2(i,j,k))*BESSEL_K1(r2(i,j,k))&
                    -BESSEL_I1(r2(i,j,k))*BESSEL_K0(r2(i,j,k))) &
             -1.0 * r3(i,j,k)/sigma*(BESSEL_I0(r3(i,j,k))*BESSEL_K1(r3(i,j,k))&
                    -BESSEL_I1(r3(i,j,k))*BESSEL_K0(r3(i,j,k)))
#endif
          END DO
        END DO
      END DO

      Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.0
      Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.0
      Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.0
    CASE(2)
      r1 = rsq(Mesh,Mesh%radius%bcenter,x1)
      r2 = rsq(Mesh,Mesh%radius%bcenter,x3)
      r3 = rsq(Mesh,Mesh%radius%bcenter,x4)
      Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) &
        = 2.0 / (2.*PI*sigma**2) * EXP(-r1/(2.*sigma**2)) &
        + 0.5 / (2.*PI*sigma**2) * EXP(-r2/(2.*sigma**2)) &
        + 1.0 / (2.*PI*sigma**2) * EXP(-r3/(2.*sigma**2))

      r1 = rsq(Mesh,Mesh%radius%faces(:,:,:,1),x1)
      r2 = rsq(Mesh,Mesh%radius%faces(:,:,:,1),x3)
      r3 = rsq(Mesh,Mesh%radius%faces(:,:,:,1),x4)

      Timedisc%solution%data4d(:,:,:,1) &
        = - 2.0 / SQRT(r1) * ERF(SQRT(r1 / 2.) / sigma) &
          - 0.5 / SQRT(r2) * ERF(SQRT(r2 / 2.) / sigma) &
          - 1.0 / SQRT(r3) * ERF(SQRT(r3 / 2.) / sigma)
      Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.0
      Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.0
      Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.0
    CASE(3)
      r1 = rsq(Mesh,Mesh%curv%bcenter,x1)
      r2 = rsq(Mesh,Mesh%curv%bcenter,x2)
      Timedisc%pvar%data4d(:,:,:,Physics%DENSITY)  = 0.02/(3.2*PI) &
        + 0.99 * (  EXP(-0.5*(r1/sigma**2))/(2.*PI*sigma**2) &
                  + EXP(-0.5*(r2/sigma**2))/(2.*PI*sigma**2))

      Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.02/(3.2*PI) &
        + GN / (2.*PI*sigma**2) &
           * (  Ei(-(r1/sigma**2)) - Ei(-0.5*(r1/sigma**2)) &
              + Ei(-(r2/sigma**2)) - Ei(-0.5*(r2/sigma**2)))

      Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.0
      Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = r - r*omega
    END SELECT

    DO dir=WEST,EAST
      SELECT TYPE(bound => Timedisc%Boundary%boundary(dir)%p)
      CLASS IS (boundary_noslip)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            DO ig=1,Mesh%GNUM
              SELECT CASE (dir)
              CASE(WEST)
                  i = Mesh%IMIN-ig
              CASE(EAST)
                  i = Mesh%IMAX+ig
              END SELECT
              bound%data(ig,j,k,Physics%YVELOCITY) &
                = Timedisc%pvar%data4d(i,j,k,Physics%YVELOCITY)
            END DO
          END DO
        END DO
      END SELECT
    END DO

    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar%data4d,Timedisc%cvar%data4d)
  END SUBROUTINE InitData

  FUNCTION rsq(Mesh,r,x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base) :: Mesh
    REAL, DIMENSION(2) :: x
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: r
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: res
    !------------------------------------------------------------------------!
    INTENT(IN) :: x
    !------------------------------------------------------------------------!
    res = r**2 + x(1)**2 - 2.*r*x(1) * COS(Mesh%curv%bcenter(:,:,:,2)-x(2))
  END FUNCTION rsq

END PROGRAM orbitingcylinders
