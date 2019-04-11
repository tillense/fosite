!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_planetheating.f90                                         #
!#                                                                           #
!# Copyright (C) 2013-2019                                                   #
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
!> \addtogroup sources
!! - parameters of \link sources_planetheating \endlink as key-values
!! \key{cvis,REAL,safety factor for numerical stability, 0.1}
!! \key{intensity,REAL,intensity of the star at 1 AU}
!! \key{albedo,REAL,albedo of the planet}
!! \key{distance,REAL,semi-major axis of planet-star}
!! \key{mu,REAL,molar mass of planetary atmsophere}
!! \key{omegasun,REAL,angular velocity of planet,0.0}
!! \key{year,REAL,time of a year}
!! \key{gacc,REAL,gravitational accelaration on the surface of the planet}
!! \key{gamma,REAL,ratio of heats for planet}
!! \key{tetha0,REAL,beginning (and maximum) inclindation to ecliptic,0.0}
!! \key{phi0,REAL,beginning azimuthal angle,0.0}
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!!
!! \brief heating of a planet by a star
!!
!! \warning use SI units
!!
!! \extends sources_common
!! \ingroup sources
!!
!! Program and data initialization for the heating of a planet with a shallow
!! atmosphere by a star.
!!
!! Includes:
!! - different intensities (distance + luminosity of star)
!! - years and seasons (inclination of ecliptic)
!! - days
!!
!! Assumes:
!! - optically thin atmosphere for infalling radiation directly at the
!!    surface
!! - parallel infall of radiation (star not too close)
!!
!! \todo include eccentricity
!!
!! See subroutine \link updateplanetheating \endlink for detailed
!! prescription.
!!
!! References:
!! \cite pierrehumbert2010
!----------------------------------------------------------------------------!
MODULE sources_planetheating_mod
  USE sources_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE marray_compound_mod
  USE marray_base_mod
  USE logging_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(sources_base) :: sources_planetheating
    CHARACTER(LEN=32) :: source_name = "thin atmosphere heating"
    TYPE(marray_base) :: Qstar          !< energy term due to stellar cooling
    TYPE(marray_base) :: T_s            !< surface temperature (black body)
    TYPE(marray_base) :: cos1,sin1      !< helping arrays for precomputation
    REAL              :: intensity      !< intensity of the star
    REAL              :: albedo         !< albedo of the star
    REAL              :: distance       !< distance of the star
    REAL              :: omegasun       !< rotational omega of the star
    REAL              :: year           !< a year
    REAL              :: theta0, phi0   !< states were heating should start
  CONTAINS
    PROCEDURE :: InitSources_planetheating
    PROCEDURE :: InfoSources
    PROCEDURE :: SetOutput
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcTimestep_single
    PROCEDURE :: UpdatePlanetHeating
    PROCEDURE :: Finalize
  END TYPE
  PUBLIC :: &
       ! types
       sources_planetheating
  !--------------------------------------------------------------------------!

CONTAINS


  !> \public Constructor of the heating module for a planetary atmosphere.
  SUBROUTINE InitSources_planetheating(this,Mesh,Physics,Fluxes,config,IO)
    USE physics_euler_mod, ONLY : physics_euler
    USE geometry_spherical_mod, ONLY : geometry_spherical
    USE constants_si_mod, ONLY : constants_si
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(fluxes_base),  INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),      POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitLogging(stype,this%source_name)

    ! some sanity checks
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("InitSources_planetheating","physics not supported")
    END SELECT

    SELECT TYPE (const => Physics%constants)
    CLASS IS(constants_si)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("InitSources_planetheating","only SI units supported")
    END SELECT

    this%Qstar = marray_base()
    this%T_s = marray_base()
    this%cos1 = marray_base(2) !TODO check if 2 is really necessary
    this%sin1 = marray_base(2) !TODO check if 2 is really necessary

    ! Courant number, i.e. safety factor for numerical stability
    CALL GetAttr(config, "cvis", this%cvis, 0.1)
    ! intensity of the planets star at 1 AU
    CALL GetAttr(config, "intensity", this%intensity)
    ! albedo of the planet
    CALL GetAttr(config, "albedo", this%albedo)
    ! distance planet-star
    CALL GetAttr(config, "distance", this%distance)
    ! day-night omega
    CALL GetAttr(config, "omegasun", this%omegasun, 0.0)
    ! year
    CALL GetAttr(config, "year", this%year)
    ! beginning points
    CALL GetAttr(config, "theta0", this%theta0, 0.0)
    CALL GetAttr(config, "phi0", this%phi0, 0.0)

    SELECT TYPE (mgeo => Mesh%Geometry)
    CLASS IS(geometry_spherical)
      ! TODO: ATTENTION CHECK THAT THIS REGION IS CORRECT AFTER TRANSITION TO NEW FOSITE
      this%cos1%data4d(:,:,:,1) = COS(Mesh%bcenter(:,:,:,2))
!      this%cos1(:,:,2) = COS(Mesh%bcenter(:,:,2))
      this%sin1%data4d(:,:,:,1) = SIN(Mesh%bcenter(:,:,:,2))
!      this%sin1(:,:,2) = SIN(Mesh%bcenter(:,:,2))
    CLASS DEFAULT
      CALL this%Error("InitSources_planetheating","only spherical geometry allowed")
    END SELECT

    ! set initial time < 0
    ! TODO: Check were this time is set? Is it correct? Why not timedisc time?
    this%time = -1.0

    ! initialize arrays
    this%Qstar%data1d(:) = 0.0
    this%T_s%data1d(:)  = 0.0

    ! initialise output
    CALL this%SetOutput(Mesh,config,IO)

    ! call InitSources in base
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

  END SUBROUTINE InitSources_planetheating

  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating), INTENT(IN) :: this
    CLASS(mesh_base),           INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: param_str
    REAL :: tmp_out
    REAL, PARAMETER   :: AU      = 1.49597870691E+11    ! astr. unit [m]     !
    !------------------------------------------------------------------------!
    tmp_out = this%year/(3.6e3*24*365)
    WRITE (param_str,'(ES8.2)') tmp_out
    CALL this%Info("            sid. year [yr]:    " // TRIM(param_str))
    tmp_out = 1./(this%omegasun*3.6e3*24.0)
    WRITE (param_str,'(ES8.2)') ABS(tmp_out)
    CALL this%Info("            day [d]:           " // TRIM(param_str))
    tmp_out = this%distance/AU
    WRITE (param_str,'(ES8.2)') tmp_out
    CALL this%Info("            distance [au]:     " // TRIM(param_str))
  END SUBROUTINE InfoSources


  !> Add the calculated sources to the energy equation.
  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_base),          INTENT(INOUT) :: Physics
    CLASS(fluxes_base),           INTENT(IN)    :: Fluxes
    CLASS(sources_base),          INTENT(INOUT) :: Sources
    REAL,                         INTENT(IN)    :: time, dt
    CLASS(marray_compound),       INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    SELECT TYPE(s => sterm)
    TYPE IS (statevector_euler)
      s%density%data1d(:)   = 0.0
      s%momentum%data1d(:) = 0.0

      CALL this%UpdatePlanetHeating(Mesh,Sources,time)

      ! radiative heating by the central stars
      s%energy%data1d(:) = this%Qstar%data1d(:)
    END SELECT
  END SUBROUTINE ExternalSources_single


  !> Caculates the timestep corresponding to the heating.
  !!
  !! The timescale is calculated by \f$ t \sim Q_{\mathrm{star}}/P \f$, where
  !! \f$ Q_{\mathrm{star}} \f$ is the heating term and \f$ P \f$ the pressure.
  SUBROUTINE CalcTimestep_single(this,Mesh,Physics,Fluxes,pvar,cvar,time,dt)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)       :: Mesh
    CLASS(physics_base), INTENT(INOUT)    :: Physics
    CLASS(fluxes_base),  INTENT(IN)       :: Fluxes
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    REAL,                INTENT(IN)       :: time
    REAL,                INTENT(OUT)      :: dt
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    ! maximum of inverse cooling timescale t_cool ~ P/Q_cool
    dt = HUGE(invdt)
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_euler)
      invdt = MAXVAL(ABS(this%Qstar%data1d(:) / p%pressure%data1d(:)), &
                     MASK=Mesh%without_ghost_zones%mask1d(:))
      IF (invdt.GT.TINY(invdt)) dt = this%cvis / invdt
    END SELECT
  END SUBROUTINE CalcTimestep_single


  !> Updates the heating for a given time
  !!
  !! \image html sphere_trafo.jpg
  !! \image latex sphere_trafo.pdf
  !!
  !! Heating is done via
  !! \f[
  !!    Q_{\mathrm{star}} = Q_{\mathrm{star,AU}} \left(
  !!        \frac{a_{\mathrm{AU}}}{a}\right)^2
  !!        \sin{\theta'}\cos{\varphi'}\left(1- \alpha \right),
  !! \f]
  !! where \f$ Q_{\mathrm{star,AU}} \f$ is the intensity at \f$ 1\,\mathrm{AU}
  !! \f$ distance, \f$ a_{\mathrm{AU}} \f$ and \f$ a \f$ are the semi-major
  !! axes and \f$ \alpha \f$ is the albedo. The angles \f$ \theta' \f$ and
  !! \f$ \varphi' \f$ are transformed angles that correspond to the view
  !! factor.
  !!
  !! In order to derive the new angles you need to use spherical trigonometry.
  !! It yields:
  !! \f[
  !!    \theta' = \arccos{(\cos{\theta}\cos{\theta_0(t)})+\sin{\theta}
  !!              \sin{\theta_0(t)}\cos{(\varphi-\varphi_0(t))}}   \\
  !!    \varphi' = \arctan{\left( \frac{\sin \theta \sin{(\varphi -
  !!                \varphi_0(t))}}{\cos{\theta}\sin{\theta_0(t)}-
  !!                \sin{\theta}\cos{\theta_0(t)}\cos{(\varphi
  !!                - \varphi_0(t))}} \right)}.
  !! \f]
  !! where \f$ \theta_0 \f$ is the angle where the sun beam falls
  !! vertically onto the planet. The value vor \f$ \varphi_0 \f$ just changes
  !! were the day begins and can be safely set to zero.
  !!
  !! In order to simulate seasons \f$ \theta_0(t) = \theta_0(0)
  !! \cos{2\pi\frac{t}{t_a}} \f$, with \f$ t_a \f$ the time of a year is
  !! applied. For days \f$ \varphi(t) = \varphi_0(0) + 2 \pi f t \f$ is
  !! calculated. Here f is the frequency [\f$ s^{-1} \f$].
  SUBROUTINE UpdatePlanetHeating(this,Mesh,Sources,time)
    USE functions
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(sources_base),          INTENT(INOUT) :: Sources
    REAL,                         INTENT(IN)    :: time
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    REAL              :: theta2,phi1,phi2
    REAL              :: distance_change
    !------------------------------------------------------------------------!
    ! TODO: This should be definitely done more generically
    REAL, PARAMETER   :: AU      = 1.49597870691E+11    ! astr. unit [m]     !
    REAL, PARAMETER   :: DAY     = 8.6400E+4            ! Day [sec]          !
    !------------------------------------------------------------------------!
    ! heating by the star
    IF (time.NE.this%time) THEN
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            !--------------------------------------------------------------!
            phi1=Mesh%bcenter(i,j,k,2) + 2.*PI*this%omegasun*this%time
            ! transformation
            theta2 = ACOS(this%cos1%data4d(i,j,k,1)*COS(this%theta0*&
                     COS(2*PI*this%time/this%year)) + &
                     this%sin1%data4d(i,j,k,1)*SIN(this%theta0*&
                     COS(2*PI*this%time/this%year))*COS(phi1-this%phi0))

            ! modulo function to prevent phi2>2*PI,phi2<0
            phi2   = MODULO(ATAN2(this%sin1%data4d(i,j,k,1)*SIN(phi1-this%phi0),&
                     (-this%cos1%data4d(i,j,k,1)*SIN(this%theta0*COS(2.*PI*this%time/&
                     this%year)) + this%sin1%data4d(i,j,k,1)*COS(phi1-this%phi0)*&
                     COS(this%theta0*COS(2*PI*this%time/this%year)))),2*PI)

            ! calculate heating source                                     !
            !--------------------------------------------------------------!
            IF (phi2.LE.PI/2..OR.phi2.GE.3./2.*PI) THEN
              ! for sin,cos have a look at spherical trigonometry
              ! the albedo can contain things like scattering etc. (not implemented)

              this%Qstar%data3d(i,j,k) = this%intensity*&
                   (AU/this%distance)**(2.)*(1.-(this%albedo))*SIN(theta2)*COS(phi2)
            ELSE
              this%Qstar%data3d(i,j,k) = 0.0
            END IF
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE UpdatePlanetHeating

  !> Sets the output parameters.
  !!
  !! Output:
  !! 1. heating term: \f$ Q_{\mathrm{star}} \f$
  SUBROUTINE SetOutput(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating) :: this
    CLASS(mesh_base), INTENT(IN) :: Mesh
    TYPE(Dict_TYP),   POINTER    :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!
    ! heating source term
    CALL GetAttr(config, "output/Qstar", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "Qstar", &
              this%Qstar%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
  END SUBROUTINE SetOutput


  !> Deallocation
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetheating), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Qstar%Destroy()
    CALL this%T_s%Destroy()
    CALL this%cos1%Destroy()
    CALL this%sin1%Destroy()
  END SUBROUTINE Finalize
END MODULE sources_planetheating_mod
