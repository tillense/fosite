!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_diskcooling.f03                                           #
!#                                                                           #
!# Copyright (C) 2011-2018                                                   #
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
!> \addtogroup sources
!! - parameters of \link sources_diskcooling \endlink as key-values
!! \key{method,INTEGER,cooling model (GRAY or GAMMIE)}
!! \key{cvis,REAL,safety factor for numerical stability, 0.1}
!! \key{Tmin,REAL,set a minimum temperature, 1.0E-30}
!! \key{rhomin,REAL,set a minimum density, 1.0E-30}
!! \key{b_cool,REAL,cooling parameter, 1.0}
!! \key{b_start,REAL,starting cooling parameter, 1.0}
!! \key{b_final,REAL,final cooling parameter, 1.0}
!! \key{t_start,REAL,time to start decreasing value of beta from b_start, 0.0}
!! \key{dt_bdec,REAL,time over which b_cool should reach its final value
!!      beginning at t_start, 0.0}
!! \key{output/Qcool,INTEGER,enable(=1) output of cooling function,0}
!----------------------------------------------------------------------------!
!> \author Anna Feiler
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief source terms module for cooling of geometrically thin
!! accretion disks
!!
!! Supported methods:
!! - \link lambda_gray Gray cooling \endlink according to Hubeny \cite hubeny1990
!!   using opacities from Bell & Lin \cite bell1994 . The Rosseland mean opacities
!!   are then computed using the interpolation formula of Gail 2003 (private communication).
!! - \link lambda_gammie Simple cooling \endlink model according to Gammie
!!   \cite gammie2001 with a constant coupling between dynamical and cooling time scale.
!!
!! \warning use SI units for gray cooling
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_diskcooling_mod
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
  INTEGER, PARAMETER :: GRAY   = 1
  INTEGER, PARAMETER :: GAMMIE = 2
  INTEGER, PARAMETER :: GAMMIE_SB = 3
  !--------------------------------------------------------------------------!
  REAL, PARAMETER :: SQRT_THREE = 1.73205080757
  REAL, PARAMETER :: SQRT_TWOPI = 2.50662827463

! Rosseland mean opacity constants in SI units;
! taken from Bell & Lin, ApJ, 427, 1994
! kappa_i= kappa_0i * rho**rexp(i) * T**Texp(i)

!!$  DOUBLE PRECISION, PARAMETER :: kappa0(8) = (/ &
!!$       2.00D-05, & ! ice grains                     [m^2/kg/K^2]
!!$       2.00D+15, & ! evaporation of ice grains      [m^2/kg*K^7]
!!$       1.00D-02, & ! metal grains                 [m^2/kg/K^0.5]
!!$       2.00D+77, & ! evaporation of metal grains [m^5/kg^2*K^24]
!!$       1.00D-11, & ! molecules                [m^4/kg^(5/3)/K^3]
!!$       1.00D-38, & ! H-scattering            [m^3/kg^(4/3)/K^10]
!!$       1.50D+16, & ! bound-free and free-free [m^5/kg^2*K^(5/2)]
!!$       KE /)       ! electron scattering                [m^2/kg]
  REAL, PARAMETER :: logkappa0(8) = (/ &
       -10.8197782844, & ! ice grains                     [m^2/kg/K^2]
       35.2319235755, &  ! evaporation of ice grains      [m^2/kg*K^7]
       -4.60517018599, & ! metal grains                 [m^2/kg/K^0.5]
       177.992199341, &  ! evaporation of metal grains [m^5/kg^2*K^24]
       -25.3284360229, & ! molecules                [m^4/kg^(5/3)/K^3]
       -87.4982335338, & ! H-scattering            [m^3/kg^(4/3)/K^10]
       37.246826596, &   ! bound-free and free-free [m^5/kg^2*K^(5/2)]
       -3.3581378922 /)  ! electron scattering                [m^2/kg]
  REAL, PARAMETER :: Texp(8) = (/ 2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0 /)
  REAL, PARAMETER :: rexp(8) = (/ 0.0, 0.0, 0.0, 1.0, 2./3., 1./3., 1.0, 0.0 /)
  REAL, PARAMETER :: T0 = 3000      ! temperature constant (opacity interpolation)

  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(sources_base) :: sources_diskcooling
    CHARACTER(LEN=32) :: source_name = "thin accretion disk cooling"
    CLASS(logging_base), ALLOCATABLE :: cooling
    REAL                            :: b_cool       !< cooling parameter (Gammie)
    REAL                            :: b_start      !< starting cooling
    REAL                            :: b_final      !< final cooling parameter
    REAL                            :: t_start      !< cooling starting time
    REAL                            :: dt_bdec      !< cooling decrease time
    TYPE(marray_base) :: Qcool        !< energy sink due to cooling
    TYPE(marray_base) :: ephir        !< azimuthal unit vector / radius
    REAL                            :: T_0          !< equilibrium temp
    REAL                            :: rho_0        !< minimum density
  CONTAINS
    PROCEDURE :: InitSources_diskcooling
    PROCEDURE :: InfoSources
    PROCEDURE :: SetOutput
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcTimestep_single
    PROCEDURE :: UpdateCooling
    PROCEDURE :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       sources_diskcooling, &
       ! constants
       GRAY, GAMMIE, GAMMIE_SB
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of disk cooling module
  SUBROUTINE InitSources_diskcooling(this,Mesh,Physics,Fluxes,config,IO)
    USE physics_euler_mod, ONLY : physics_euler
    USE geometry_cylindrical_mod, ONLY : geometry_cylindrical
    USE geometry_cartesian_mod, ONLY : geometry_cartesian
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling) :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),       POINTER    :: config,IO
    !------------------------------------------------------------------------!
    INTEGER :: stype
    INTEGER :: cooling_func
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitLogging(stype,this%source_name)

    ! some sanity checks
    ! isothermal modules are excluded
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler)
      ! do nothing
    CLASS DEFAULT
      ! abort
      CALL this%Error("InitSources_diskcooling","physics not supported")
    END SELECT

    SELECT CASE(cooling_func)
    CASE(GAMMIE)
      ! check for geometry
      SELECT TYPE (mgeo => Mesh%Geometry)
      CLASS IS(geometry_cylindrical)
        ! do nothing
      CLASS DEFAULT
        CALL this%Error("InitSources_diskcooling","geometry not supported")
      END SELECT
    CASE(GAMMIE_SB)
      SELECT TYPE (mgeo => Mesh%Geometry)
      TYPE IS(geometry_cartesian)
        ! do nothing
      CLASS DEFAULT
        CALL this%Error("InitSources_diskcooling","geometry not supported")
      END SELECT
    END SELECT

    ! check for dimensions
    IF (Physics%VDIM.GT.2) THEN
      CALL this%Error("InitSources_diskcooling","only 2D simulations supported")
    END IF

    ! get cooling method
    CALL GetAttr(config,"method",cooling_func)
    ALLOCATE(logging_base::this%cooling)
    SELECT CASE(cooling_func)
    CASE(GRAY)
      IF (Physics%constants%GetType().NE.SI) &
         CALL this%Error("InitSources_diskcooling","only SI units supported for gray cooling")
      CALL this%cooling%InitLogging(GRAY,"gray cooling")

      this%Qcool = marray_base()

    CASE(GAMMIE)
      CALL this%cooling%InitLogging(GAMMIE,"Gammie cooling")

      this%Qcool = marray_base()
      this%ephir = marray_base(2)

    CASE(GAMMIE_SB)
      CALL this%cooling%InitLogging(GAMMIE_SB,"Gammie cooling (SB)")

      this%Qcool = marray_base()

    CASE DEFAULT
       CALL this%Error("UpdateCooling","Cooling function not supported!")
    END SELECT

    ! Courant number, i.e. safety factor for numerical stability
    CALL GetAttr(config, "cvis", this%cvis, 0.1)

    ! minimum temperature
    CALL GetAttr(config, "Tmin", this%T_0, 1.0E-30)

    ! minimum density
    CALL GetAttr(config, "Rhomin", this%rho_0, 1.0E-30)

    ! cooling time
    CALL GetAttr(config, "b_cool", this%b_cool, 1.0E+00)

    ! initial and final cooling time
    CALL GetAttr(config, "b_start", this%b_start, 0.0)
    CALL GetAttr(config, "b_final", this%b_final, this%b_cool)

    ! starting point for beta
    CALL GetAttr(config, "t_start", this%t_start, 0.0)

    ! timescale for changing beta
    CALL GetAttr(config, "dt_bdec", this%dt_bdec, -1.0)

    ! set initial time < 0
    this%time = -1.0

    ! initialize arrays
    this%Qcool%data1d(:)  = 0.0
    SELECT CASE(cooling_func)
    CASE(GRAY)
      ! do nothing
    CASE(GAMMIE)
      this%ephir%data4d(:,:,:,1) = -Mesh%posvec%bcenter(:,:,:,2)/Mesh%radius%bcenter(:,:,:)**2
      this%ephir%data4d(:,:,:,2) = Mesh%posvec%bcenter(:,:,:,1)/Mesh%radius%bcenter(:,:,:)**2
    CASE(GAMMIE_SB)
      ! do nothing
    END SELECT

    !initialise output
    CALL this%SetOutput(Mesh,config,IO)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

  END SUBROUTINE InitSources_diskcooling


  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling), INTENT(IN) :: this
    CLASS(mesh_base),           INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: param_str
    !------------------------------------------------------------------------!
    CALL this%Info("            cooling function:  " // TRIM(this%cooling%GetName()))
    SELECT CASE(this%cooling%GetType())
    CASE(GRAY)
       WRITE (param_str,'(ES8.2)') this%T_0
       CALL this%Info("            min. temperature:  " // TRIM(param_str))
       WRITE (param_str,'(ES8.2)') this%rho_0
       CALL this%Info("            minimum density:   " // TRIM(param_str))
    CASE(GAMMIE)
       WRITE (param_str,'(ES8.2)') this%b_cool
       CALL this%Info("            cooling parameter: " // TRIM(param_str))
    CASE(GAMMIE_SB)
       WRITE (param_str,'(ES8.2)') this%b_cool
       CALL this%Info("            cooling parameter: " // TRIM(param_str))
       IF(this%dt_bdec.GE.0.0) THEN
        WRITE (param_str,'(ES8.2)') this%b_start
        CALL this%Info("            initial cooling parameter: " // TRIM(param_str))
        WRITE (param_str,'(ES8.2)') this%b_final
        CALL this%Info("            final cooling parameter: " // TRIM(param_str))
        WRITE (param_str,'(ES8.2)') this%t_start
        CALL this%Info("            starting b_dec time: " // TRIM(param_str))
        WRITE (param_str,'(ES8.2)') this%dt_bdec
        CALL this%Info("            operating time: " // TRIM(param_str))
       END IF
    END SELECT
  END SUBROUTINE InfoSources


  SUBROUTINE SetOutput(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling)  :: this
    CLASS(mesh_base), INTENT(IN) :: Mesh
    TYPE(Dict_TYP), POINTER     :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!

    !cooling energy source term
    CALL GetAttr(config, "output/Qcool", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "Qcool", &
         this%Qcool%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

  END SUBROUTINE SetOutput


  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)         :: Mesh
    CLASS(physics_base),INTENT(INOUT)   :: Physics
    CLASS(sources_base),INTENT(INOUT)   :: Sources
    CLASS(fluxes_base),INTENT(IN)       :: Fluxes
    REAL,INTENT(IN)                     :: time, dt
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    sterm%data2d(:,Physics%DENSITY) = 0.0
    sterm%data2d(:,Physics%XMOMENTUM) = 0.0
    sterm%data2d(:,Physics%YMOMENTUM) = 0.0

    SELECT TYPE(phys => Physics)
    TYPE IS (physics_euler)
       CALL this%UpdateCooling(Mesh,phys,Sources,time,pvar)
    END SELECT
    ! energy loss due to radiation processes
    sterm%data2d(:,Physics%ENERGY) = -this%Qcool%data1d(:)
  END SUBROUTINE ExternalSources_single


  SUBROUTINE CalcTimestep_single(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    REAL,                INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar,cvar
    REAL,                INTENT(OUT)   :: dt
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    ! maximum of inverse cooling timescale t_cool ~ P/Q_cool
    invdt = MAXVAL(ABS(this%Qcool%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) &
         / pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%PRESSURE)))
    IF (invdt.GT.TINY(invdt)) THEN
       dt = this%cvis / invdt
    ELSE
       dt = HUGE(invdt)
    END IF
  END SUBROUTINE CalcTimestep_single

  !> \private Updates the cooling function at each time step.
  SUBROUTINE UpdateCooling(this,Mesh,Physics,Sources,time,pvar)
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling)         :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_euler),INTENT(INOUT) :: Physics
    CLASS(sources_base), INTENT(INOUT) :: Sources
    CLASS(marray_compound),INTENT(INOUT):: pvar
    REAL,                INTENT(IN)    :: time
    !------------------------------------------------------------------------!
    REAL              :: muRgamma,Qfactor
    !------------------------------------------------------------------------!
    ! calculate value for beta
    IF (this%dt_bdec.GE.0.0) THEN
      IF (time .LT. this%t_start) THEN
        this%b_cool = this%b_start
      ELSE IF ((time .GE. this%t_start) .AND. ((time-this%t_start) .LT. &
      this%dt_bdec)) THEN
        this%b_cool = this%b_start +  &
            (this%b_final-this%b_start)/this%dt_bdec*(time - this%t_start)
      ELSE IF ((time-this%t_start) .GE. this%dt_bdec) THEN
        this%b_cool = this%b_final
      END IF
    END IF

    ! energy loss due to radiation processes
    SELECT CASE(this%cooling%GetType())
       CASE(GRAY)
          ! some constants
          muRgamma = Physics%mu/(Physics%Constants%RG*Physics%gamma)
          Qfactor  = 8./3.*Physics%Constants%SB
          ! compute gray cooling term
          this%Qcool%data1d(:) = Lambda_gray(pvar%data2d(:,Physics%DENSITY),Sources%height%data1d(:), &
               muRgamma*Physics%bccsound%data1d(:)*Physics%bccsound%data1d(:), &
               this%rho_0,this%T_0,Qfactor)
       CASE(GAMMIE)
          ! compute Gammie cooling term with
          ! t_cool = b_cool / Omega
          ! and Omega = ephi*v / r
          this%Qcool%data3d(:,:,:) = Lambda_gammie(pvar%data4d(:,:,:,Physics%PRESSURE) / (Physics%gamma-1.), &
              ABS((this%ephir%data4d(:,:,:,1)*pvar%data4d(:,:,:,Physics%XVELOCITY) &
                  +this%ephir%data4d(:,:,:,2)*(pvar%data4d(:,:,:,Physics%YVELOCITY)&
                    +Mesh%OMEGA*Mesh%radius%bcenter(:,:,:)))) / this%b_cool)
       CASE(GAMMIE_SB)
          ! in sb t_cool = b_cool
          this%Qcool%data1d(:) = Lambda_gammie(pvar%data2d(:,Physics%PRESSURE) / (Physics%gamma-1.), &
              Mesh%OMEGA/this%b_cool)
    END SELECT

  END SUBROUTINE UpdateCooling

  ELEMENTAL FUNCTION RosselandMeanOpacity_new(logrho,logT) RESULT(kappa_R)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: logrho,logT
    REAL             :: kappa_R
    !------------------------------------------------------------------------!
    REAL :: kappa4(8)
    REAL :: Tfactor
    !------------------------------------------------------------------------!
    ! compute kappa_i^4 for all cooling mechanisms
!NEC$ unroll(8)
    kappa4(:) = EXP(MAX(-40.0,MIN(40.0, &
         4.*(logkappa0(:)+rexp(:)*logrho+Texp(:)*logT))))
    ! compute (T/T0)**10
    Tfactor = EXP(MAX(-40.0,MIN(40.0,10.*(logT-LOG(T0)))))
    ! compute Rosseland mean using Gails interpolation formula
    kappa_R = 1. /(TINY(kappa_R) + &
           (1./kappa4(1) + 1./(1.+Tfactor)/(kappa4(2)+kappa4(3)))**0.25 &
         + (1./(kappa4(4)+kappa4(5)+kappa4(6)) + 1./(kappa4(7)+kappa4(8)))**0.25)
  END FUNCTION RosselandMeanOpacity_new

  !> \private Gray cooling
  !!
  !! The cooling function is given by
  !! \f[
  !!    \Lambda= 2\sigma T_{eff}^4
  !! \f]
  !! where \f$ \sigma \f$ is the Stefan-Boltzmann constant (see e. g. Pringle
  !! \cite pringle1981 ). If the disk is optically thick for its own radiation,
  !! one can use the radiation diffusion approximation and relate the effective
  !! temperature to the midplane temperature according to
  !! \f[
  !!    T_{eff}^4 = \frac{8}{3} \frac{T_c^4}{\tau_{eff}}
  !! \f]
  !! where \f$ \tau_{eff} \f$ is an effective optical depth (see e. g. Hubeny
  !! \cite hubeny1990 ).
  !!
  ELEMENTAL FUNCTION Lambda_gray(Sigma,h,Tc,rho0,T0,Qf) RESULT(Qcool)
    IMPLICIT NONE
    !---------------------------------------------------------------------!
    REAL, INTENT(IN) :: Sigma,h,Tc,rho0,T0,Qf
    REAL :: Qcool
    REAL :: logrho,logT,kappa,tau,tau_eff,T0Tc
    !---------------------------------------------------------------------!
    ! logarithm of midplane density
    ! log(SQRT(2*Pi)^-1 * Sigma / H ) = -log(SQRT(2*Pi) * H / Sigma)
    logrho = -LOG(MAX(rho0,SQRT_TWOPI * h / Sigma))
    ! logarithm of midplane temperature
    logT = LOG(MAX(T0,Tc))
    ! compute Rosseland mean absorption coefficient using Gails formula
!CDIR IEXPAND
    kappa = RosselandMeanOpacity_new(logrho,logT)
    ! optical depth
    tau = 0.5*kappa*Sigma
    ! effective optical depth
    tau_eff = 0.5*tau + 1./(3.*tau) + 1./SQRT_THREE
    ! temperature ratio
    T0Tc = T0/Tc
    ! cooling term
    Qcool = Qf/tau_eff * EXP(4.0*logT) * (1.0-T0Tc*T0Tc*T0Tc*T0Tc) ! = Tc**4-T0**4
  END FUNCTION Lambda_gray

  !> \private Gammie cooling
  !!
  !! The cooling function is given by
  !! \f[
  !!    \Lambda= -E_{int}/t_{cool}
  !! \f]
  !! with the cooling time scale \f$ t_{cool} = b \Omega^{-1} \f$
  ELEMENTAL FUNCTION Lambda_gammie(Eint,t_cool_inv) RESULT(Qcool)
    IMPLICIT NONE
    !---------------------------------------------------------------------!
    REAL, INTENT(IN) :: Eint,t_cool_inv
    REAL :: Qcool
    !---------------------------------------------------------------------!
    Qcool = Eint * t_cool_inv
  END FUNCTION Lambda_gammie

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_diskcooling), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    SELECT CASE(this%GetType())
    CASE(GRAY)
      CALL this%Qcool%Destroy()
    CASE(GAMMIE)
      CALL this%Qcool%Destroy()
      CALL this%ephir%Destroy()
    CASE(GAMMIE_SB)
      CALL this%Qcool%Destroy()
    END SELECT
  END SUBROUTINE Finalize

END MODULE sources_diskcooling_mod
