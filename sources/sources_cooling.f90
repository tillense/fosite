!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_cooling.f90                                               #
!#                                                                           #
!# Copyright (C) 2009-2024                                                   #
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
!> \addtogroup sources
!! - parameters of \link sources_cooling_mod sources_cooling\endlink as key-values
!! \key{cvis,REAL,safety factor for numerical stability, 0.1}
!! \key{switchon,REAL,soft switch on,-1.0}
!! \key{Tmin,REAL,set a minimum temperature, 1.0E-30}
!! \key{output/Qcool,INTEGER,enable(=1) output of cooling function,0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief source terms module for simple optically thin cooling
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_cooling_mod
  USE sources_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE marray_compound_mod
  USE marray_base_mod
  USE sources_c_accel_mod
  USE logging_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "optically thin cooling"
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(sources_base) :: sources_cooling
    TYPE(marray_base), ALLOCATABLE  :: Q        !< energy sink due to cooling
    REAL      :: switchon
    REAL      :: Tmin !< temperature minimum
  CONTAINS
    PROCEDURE :: InitSources
    PROCEDURE :: SetOutput
    PROCEDURE :: ExternalSources
    PROCEDURE :: CalcTimestep
    PROCEDURE :: UpdateCooling
    FINAL :: Finalize
  END TYPE
 
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    ! types
    sources_cooling
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources(this,Mesh,Physics,Fluxes,config,IO)
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_cooling), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),       POINTER    :: config,IO
    !------------------------------------------------------------------------!
    INTEGER :: stype,err
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    ! call basic initialization subroutine
    CALL this%InitSources_base(stype,source_name)

    ! Courant number, i.e. safety factor for numerical stability
    CALL GetAttr(config, "cvis", this%cvis, 0.1)

    ! minimum temperature
    CALL GetAttr(config, "Tmin", this%Tmin, 1.0E-30)

    ! soft switch on time scale for activating cooling (negative disables)
    CALL GetAttr(config, "switchon", this%switchon, -1.0)

    ! some sanity checks
    ! isothermal modules are excluded
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler)
      ! do nothing
    CLASS DEFAULT
      ! abort
      CALL this%Error("sources_cooling::InitSources","physics not supported")
    END SELECT

    ALLOCATE(this%Q,STAT=err)
    IF (err.NE.0) CALL this%Error("sources_cooling::InitSources","memory allocation failed")
    this%Q = marray_base()

    ! set initial time < 0
    this%time = -1.0

    ! initialize cooling function
    this%Q%data1d(:)  = 0.0
    
    ! register ouput arrays
    CALL this%SetOutput(Mesh,config,IO)
  END SUBROUTINE InitSources


  SUBROUTINE SetOutput(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_cooling), INTENT(INOUT) :: this
    CLASS(mesh_base), INTENT(IN) :: Mesh
    TYPE(Dict_TYP), POINTER     :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!
    ! register cooling term for output
    CALL GetAttr(config, "output/Qcool", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "Qcool", &
         this%Q%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
  END SUBROUTINE SetOutput


  SUBROUTINE ExternalSources(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_cooling), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)         :: Mesh
    CLASS(physics_base),INTENT(INOUT)   :: Physics
    CLASS(sources_base),INTENT(INOUT)   :: Sources
    CLASS(fluxes_base),INTENT(IN)       :: Fluxes
    REAL,INTENT(IN)                     :: time, dt
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    SELECT TYPE(s => sterm)
    TYPE IS (statevector_euler)
      s%density%data1d(:) = 0.0
      s%momentum%data1d(:) = 0.0

      SELECT TYPE(phys => Physics)
      TYPE IS (physics_euler)
        SELECT TYPE(p => pvar)
        TYPE IS (statevector_euler)
          CALL this%UpdateCooling(Mesh,phys,time,p)
        END SELECT
      END SELECT

      ! energy loss due to radiation processes
      s%energy%data1d(:) = -this%Q%data1d(:)

    END SELECT
  END SUBROUTINE ExternalSources

 
  !> \public caculates the limiting time step due to cooling
  !!
  !! The timescale is calculated by \f$ t \sim Q_{\mathrm{cool}}/P \f$, where
  !! \f$ Q_{\mathrm{cool}} \f$ is the heating term and \f$ P \f$ the pressure.
  SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes,pvar,cvar,time,dt,dtcause)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_cooling), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    REAL,                INTENT(IN)    :: time
    REAL,                INTENT(INOUT) :: dt
    INTEGER,             INTENT(OUT)   :: dtcause
    !------------------------------------------------------------------------!
    REAL              :: invdt
    !------------------------------------------------------------------------!
    ! maximum of inverse cooling timescale t_cool ~ P/Q_cool
    dt = HUGE(invdt)
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_euler)
      SELECT TYPE(phys => Physics)
      TYPE IS (physics_euler)
        SELECT TYPE(p => pvar)
        TYPE IS (statevector_euler)
          CALL this%UpdateCooling(Mesh,phys,time,p)
          invdt = MAXVAL(ABS(this%Q%data1d(:) / p%pressure%data1d(:)), &
                        MASK=Mesh%without_ghost_zones%mask1d(:))
          IF (invdt.GT.TINY(invdt)) dt = this%cvis / invdt
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcTimestep


  !> \private Updates the cooling function at each time step.
  SUBROUTINE UpdateCooling(this,Mesh,Physics,time,pvar)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_cooling),       INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_euler),         INTENT(IN)    :: Physics
    REAL,                         INTENT(IN)    :: time
    CLASS(statevector_euler),     INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    REAL              :: muRgamma,Namu,scaling
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_euler)
      ! modify Qcool during switchon phase
      IF (time.GE.0.0.AND.time.LE.this%switchon) THEN
        ! compute new scaling factor < 1 during switchon phase
        scaling = SIN(0.5*PI*time/this%switchon)**2
      ELSE
        scaling = 1.0
      END IF
      ! conversion factor cs2 -> T / Kelvin
      muRgamma = Physics%mu /(Physics%Constants%RG*Physics%gamma) &
        / Physics%Constants%cf_temperature
      ! conversion factor mass density -> particle density / m^3
      Namu = Physics%Constants%NA/Physics%mu &
        / Physics%Constants%cf_density * Physics%Constants%cf_mass
      ! compute cooling source term ~ n^2 * Λ(T) with
      ! particle density n / m^3 and temperature T / K
      ! return value of lambda is given SI units i.e. W/m^3
      WHERE (Mesh%without_ghost_zones%mask1d(:))
        this%Q%data1d(:) = scaling * Physics%Constants%cf_energy &
          * (Namu * p%density%data1d(:))**2 &
          * lambda(muRgamma * Physics%bccsound%data1d(:)**2)
      ELSEWHERE
        this%Q%data1d(:) = 0.0
      END WHERE
    END SELECT
  END SUBROUTINE UpdateCooling


  !> \private simple optically thin cooling function
  !! ATTENTION: input data should be in SI-units, output is provided in W/m^3
  ELEMENTAL FUNCTION lambda(T) RESULT(L)
    IMPLICIT NONE
    REAL, INTENT(IN) :: T
    REAL :: L
    IF (T.LT.1.26E+4) THEN
       ! disable cooling for T < 1.26e4 K, i.e. set L to 0
       ! this is a very very rough approximation
       ! dont't trust this function for T < 1.26e4 !!!
       L = 0.0
    ELSE IF ((T.GE.1.26E+4).AND.(T.LT.1.5E+5)) THEN
       L = 3.7D-38 * EXP(0.58*LOG(T))
    ELSE IF ((T.GE.1.5e5).AND.(T.LT.2.2E+7)) THEN
       L = 1.2D-31 * EXP(-0.68*LOG(T))
    ELSE
       L = 1.2D-39 * EXP(0.41*LOG(T))
    END IF
  END FUNCTION lambda


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_cooling), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%Q)
  END SUBROUTINE Finalize

END MODULE sources_cooling_mod
