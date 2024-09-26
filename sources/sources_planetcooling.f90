!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_planetcooling.f90                                         #
!#                                                                           #
!# Copyright (C) 2011-2024                                                   #
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
!! - parameters of \link sources_planetcooling \endlink as key-values
!! \warning use SI units for initialization parameters
!! \key{cvis,REAL,safety factor for numerical stability, 0.9}
!! \key{intensity,REAL,intensity of the star at 1 AU}
!! \key{albedo,REAL,albedo of planetary atmosphere}
!! \key{distance,REAL,semi-major axis of planetary orbit}
!! \key{T_0,REAL,long term equilibrium surface temperature of the planet}
!! \key{gacc,REAL,gravitational accelaration on the planets surface}
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!! \author Tobias Illenseer
!!
!! \brief gray cooling of planetary atmospheres
!!
!! Program and data initialization for gray cooling of geometrically thin
!! planetary atmospheres. The cooling term is implemented as a sink in
!! the energy equation and accounts for outgoing longwave radiation (OLR).
!! For more details see \link updateplanetcooling \endlink .
!!
!! References:
!! - \cite pierrehumbert2010 R. T. Pierrehumbert, Principles of Planetary Climate,
!!     Cambridge University Press (2010)
!!
!! \extends sources_cooling
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_planetcooling_mod
  USE sources_base_mod
  USE sources_cooling_mod
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
  CHARACTER(LEN=32), PARAMETER :: source_name = "cooling of planetary atmosphere"
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(sources_cooling) :: sources_planetcooling
    TYPE(marray_base), ALLOCATABLE &
                      :: T_s, &         !< temperature at the surface
                         rho_s, &       !< density at the surface
                         P_s            !< pressure at the surface
    REAL              :: T_0            !< equilibrium temperature
    REAL              :: distance       !< distance of the star
    REAL              :: intensity      !< intensity of the star at 1 au
    REAL              :: albedo         !< albedo of the planet
    REAL              :: gacc           !< gravitational acceleration at surface
    REAL              :: tau_inf        !< optical depth in infinity
    REAL              :: const1, const2 !< two constants
  CONTAINS
    PROCEDURE :: InitSources
    PROCEDURE :: InfoSources
    PROCEDURE :: SetOutput
    PROCEDURE :: UpdateCooling
    FINAL :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       sources_planetcooling
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of the cooling module for a planetary atmosphere
  SUBROUTINE InitSources(this,Mesh,Physics,Fluxes,config,IO)
    USE physics_euler_mod, ONLY : physics_euler
    USE constants_si_mod, ONLY : constants_si
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetcooling), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(fluxes_base),  INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),      POINTER    :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    REAL              :: intensity
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitSources_base(stype,source_name)

    ! some sanity checks
    ! some sanity checks
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_euler)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("sources_planetheating::InitSources","physics not supported")
    END SELECT

    SELECT TYPE (const => Physics%constants)
    CLASS IS(constants_si)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("sources_planetheating::InitSources","only SI units supported")
    END SELECT

    ALLOCATE(this%Q,this%T_s)
    this%Q = marray_base()
    this%T_s = marray_base()

    ! Courant number, i.e. safety factor for numerical stability
    CALL GetAttr(config, "cvis", this%cvis, 0.9)
    ! intensity of the planets star at 1 AU
    CALL GetAttr(config, "intensity", this%intensity)
    ! distance planet-star
    CALL GetAttr(config, "distance", this%distance)
    ! albedo of the planet
    CALL GetAttr(config, "albedo", this%albedo)
    ! equilibrium temperature
    CALL GetAttr(config,"T_0", this%T_0)
    ! gravitational acceleration
    CALL GetAttr(config, "gacc", this%gacc)

    ! set initial time < 0
    this%time = -1.0

    this%Q%data1d(:) = 0.0
    this%T_s%data1d(:)   = this%T_0

    ! initial calculation of optical thickness
    SELECT TYPE(phys => Physics)
    CLASS IS(physics_euler)
      intensity = (this%intensity/4.)*(phys%Constants%AU/this%distance)**(2.)
      this%tau_inf= (this%T_0 * (Gamma(1.0+4.*(phys%gamma-1.)/phys%gamma) &
        / ((1.-this%albedo)*intensity/phys%Constants%SB))**0.25)**(phys%gamma/(phys%gamma-1.))
      this%const1 = phys%mu/(phys%gamma*phys%Constants%RG)
      this%const2 = phys%Constants%SB*this%tau_inf**(-4.*(phys%gamma-1.)/phys%gamma) &
                  *Gamma(1.0+4.*(phys%gamma-1.)/phys%gamma)
    END SELECT
    CALL this%SetOutput(Mesh,config,IO)

    ! call InitSources in base
    CALL this%InfoSources(Mesh)
  END SUBROUTINE InitSources


  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetcooling), INTENT(IN) :: this
    CLASS(mesh_base),             INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: param_str
    REAL, PARAMETER   :: AU      = 1.49597870691E+11    ! astr. unit [m]     !
    !------------------------------------------------------------------------!
    WRITE (param_str,'(ES10.2)') (this%intensity/4.)*(AU/this%distance)**(2.)
    CALL this%Info("            intensity:         " // TRIM(param_str) // " W/m^2")
    WRITE (param_str,'(ES10.2)') this%T_0
    CALL this%Info("            mean equil. temp.: " // TRIM(param_str) // " K")
    WRITE (param_str,'(ES10.2)') this%tau_inf
    CALL this%Info("            opt. depth:        " // TRIM(param_str))
    WRITE (param_str,'(ES10.2)') this%albedo
    CALL this%Info("            albedo:            " // TRIM(param_str))
  END SUBROUTINE InfoSources


  !> Updates the cooling for a given time
  !!
  !! The cooling term is in its essence a modified Stefan-Boltzman law, but
  !! it includes many assumptions about the considered atmosphere
  !! (see \cite pierrehumbert2010 199 pp.):
  !!  - all heating/cooling processes are within the troposphere
  !!  - the atmosphere is plane-parallel
  !!  - the vertical stratification is dry-adiabatic
  !!  - all the greenhouse gases are integrated and assumed to be in infrared
  !!  - a not to small optical depth \f$ \tau_{\infty} \f$ (\f$ \geq 5 \f$)
  !!    (in infrared)
  !!
  !! This two-stream approximation leads to the outgoing-longwave radiation
  !! (OLR)
  !! \f[
  !!    I_{+,\infty} = \sigma T_{\mathrm{s}}^4
  !!              \tau_{\infty}^{4\frac{\gamma - 1}{\gamma}}
  !!              \Gamma \left( 1 + \frac{4(\gamma -1)}{\gamma} \right),
  !! \f]
  !! with \f$ T_{\mathrm{s}} \f$ the surface temperature and  \f$ \Gamma \f$
  !! the Gamma function
  !!
  !! \image html olr.jpg
  SUBROUTINE UpdateCooling(this,Mesh,Physics,time,pvar)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    USE functions, ONLY : LnGamma
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetcooling), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_euler),         INTENT(IN)    :: Physics
    REAL,                         INTENT(IN)    :: time
    CLASS(statevector_euler),     INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    REAL              :: Qcool
    REAL              :: const1,const2       ! needed due to performance issues
    !------------------------------------------------------------------------!
    ! calculate the cooling function if time has changed
    IF (time.NE.this%time) THEN
      ! calculate planet-surface temperature using integrated
      ! pressure and density
      this%T_s%data1d(:) = this%const1 * pvar%pressure%data1d(:)/pvar%density%data1d(:)

      ! calculate surface density and pressure if output is requested
      ! (irrelevant for cooling function)
      IF (ALLOCATED(this%P_s)) this%P_s%data1d(:) = pvar%density%data1d(:)*this%gacc
      IF (ALLOCATED(this%rho_s)) this%rho_s%data1d(:) = Physics%mu/(Physics%Constants%RG) &
        *this%P_s%data1d(:)/this%T_s%data1d(:)

      ! cooling source
      WHERE (Mesh%without_ghost_zones%mask1d(:))
        this%Q%data1d(:) = this%const2*this%T_s%data1d(:)**4
      ELSEWHERE
        this%Q%data1d(:) = 0.0
      END WHERE

      ! update time
      this%time=time
    END IF
  END SUBROUTINE UpdateCooling


  !> Sets the output parameters.
  !!
  !! Output:
  !! 1. cooling term: \f$ Q_{\mathrm{cool}} \f$
  !! 2. surface temperature: \f$ T_{\mathrm{s}} \f$
  !! 3. surface pressure: \f$ P_{\mathrm{s}} \f$
  !! 4. surface density: \f$ \varrho_{\mathrm{s}} \f$
  !!
  !! \todo allocate output arrays P_s and RHO_s only if requested
  SUBROUTINE SetOutput(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_planetcooling) :: this
    CLASS(mesh_base),             INTENT(IN) :: Mesh
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!

    !cooling source term
    CALL GetAttr(config, "output/Qcool", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "Qcool", &
              this%Q%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    ! temperature
    CALL GetAttr(config, "output/T_s", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "T_s", &
              this%T_s%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    ! surface pressure
    CALL GetAttr(config, "output/P_s", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
      ALLOCATE(this%P_s)
      this%P_s = marray_base()
      this%P_s%data1d(:)   = 0.0
      CALL SetAttr(IO, "P_s", &
              this%P_s%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    END IF

    ! surface density
    CALL GetAttr(config, "output/rho_s", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
      ALLOCATE(this%rho_s)
      this%rho_s = marray_base()
      this%rho_s%data1d(:) = 0.0
      CALL SetAttr(IO, "rho_s", &
              this%rho_s%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    END IF
  END SUBROUTINE SetOutput


  !> Destructor
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_planetcooling), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%T_s)
    IF (ALLOCATED(this%P_s)) DEALLOCATE(this%P_s)
    IF (ALLOCATED(this%rho_s)) DEALLOCATE(this%rho_s)
    ! this%Q is deallocated in parent class destructor
  END SUBROUTINE Finalize

END MODULE sources_planetcooling_mod
