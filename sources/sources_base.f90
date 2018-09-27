!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2013                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup sources
!! - general parameters of sources group as key-values
!! \key{stype,INTEGER,Type of source}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!!
!! \brief generic source terms module providing functionaly common
!! to all source terms
!!
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE physics_base_mod!, GeometricalSources_Physics => GeometricalSources, &
     !  ExternalSources_Physics => ExternalSources
  USE fluxes_base_mod
!  USE mesh_common, ONLY : Mesh_TYP
!  USE timedisc_common, ONLY : Timedisc_TYP
!  USE sources_c_accel, InitSources_common => InitSources, &
!       CloseSources_common => CloseSources
!  USE sources_generic_mod
!  USE sources_diskthomson
!  USE sources_viscosity_mod
!  USE sources_wave_damping
!  USE sources_cooling
!  USE sources_stellarheating
!  USE sources_rotframe
!  USE sources_sgs
!  USE sources_diskcooling
!  USE sources_planetheating
!  USE sources_planetcooling
!  USE sources_forcing
!  USE sources_shearbox
!  USE gravity_generic
!  USE physics_generic_mod 
!  USE fluxes_generic
!  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE,ABSTRACT, EXTENDS(logging_base) :: sources_base
     !> \name Variables
     CLASS(sources_base), POINTER    :: next => null() !< next source in list
     !TYPE(Gravity_TYP), POINTER      :: glist => null()!< gravity list
     CLASS(sources_base), POINTER     :: viscosity    !< molecular,alpha,beta
     !TYPE(Common_TYP)                :: cooling      !< gray,gammie cooling func.
     REAL                            :: time         !< simulation time
     INTEGER                         :: timeid       !<    update of this id?
     !> 0: no src term in energy equation
     !! 1: src term in energy equation
     LOGICAL                         :: addtoenergy
!FIXME
     LOGICAL                         :: update_disk_height !< enable/disable computation of disk scale height
     REAL                            :: mass         !< mass for diskthomson
     REAL                            :: mdot         !< disk accretion rate
     REAL                            :: dynconst,bulkconst ! viscosity const.
     REAL                            :: cvis         !< viscous Courant no.
     REAL                            :: eps1,eps2    !< softening parameter
     REAL                            :: gparam       !< geometry parameter
     !> effective surface of the dust grains
     REAL                            :: a_eff
     !> optical properties of the disk surface
     REAL                            :: kappa
     REAL                            :: T_star       !< temperature star
     REAL                            :: R_star       !< radius star
     REAL                            :: Qabs         !< dust absorption efficiency
     REAL                            :: T_sublim_min !< dust starts to sublimate
     !> dust density = 0.0 due to sublimation
     REAL                            :: T_sublim_max
     !> \name
     !!#### wave_damping

     !> inner and outer wave damping boundaries
     REAL, DIMENSION(2)              :: r
     !> time of a orbital period at the inner and outer boundaries
     REAL, DIMENSION(2)              :: tau
     !> \name
     !!#### diskcooling
     REAL                            :: b_cool       !< cooling parameter (Gammie)
     REAL, DIMENSION(:,:), POINTER   :: Qcool        !< energy sink due to cooling
     REAL, DIMENSION(:,:,:), POINTER :: ephir        !< azimuthal unit vector / radius
     !> \name
     !!####
     INTEGER                         :: dust_type    !< select mc3d dust catalogue
     !> 1: binary primary component or single star heating
     !! 2: secondary
     INTEGER                         :: star
     REAL, DIMENSION(:,:,:,:), POINTER :: accel,accart !< acceleration
     REAL, DIMENSION(:,:,:,:), POINTER :: bcposvec,bccart !< position vector
     REAL, DIMENSION(:,:,:), POINTER   :: radius       !< distance to origin
     REAL, DIMENSION(:,:,:), POINTER   :: invr         !< 1./radius
     REAL, DIMENSION(:,:,:), POINTER   :: cs           !< speed of sound
     REAL, DIMENSION(:,:,:), POINTER   :: omega        !< angular velocity
     REAL, DIMENSION(:,:,:,:), POINTER :: omega2       !< Omega Kepler squared
     REAL, DIMENSION(:,:,:), POINTER :: gxr3         !< = GN*x/radius**3
     REAL, DIMENSION(:,:,:), POINTER   :: cellmass     !< rho*dV
     ! --------------------------------------------------------------------- !
     !> \name
     !!#### planet_heating/cooling
     REAL, DIMENSION(:,:), POINTER   :: Energy       !< energy
     REAL, DIMENSION(:,:), POINTER   :: RHO_s        !< surface density
     REAL, DIMENSION(:,:), POINTER   :: P_s          !< surface pressure
     REAL, DIMENSION(:,:), POINTER   :: T_s          !< temperatur
     REAL                            :: tau_inf      !< optical depth
     REAL                            :: T_0          !< equilibrium temp
     REAL                            :: rho_0        !< minimum density
     REAL, DIMENSION(:,:), POINTER   :: T_init       !< initial temperatur
     REAL                            :: intensity
     REAL                            :: albedo
     REAL                            :: distance     !< distance planet<->star
     REAL                            :: R_planet     !< radius planet
     REAL                            :: mu           !< molar mass
     REAL                            :: theta0       !< shifting theta-angle
     REAL                            :: phi0         !< shifting phi-angle
     REAL                            :: omegasun     !< day-night omega
     REAL                            :: year         !< trop. year of a planet
     REAL, DIMENSION(:,:,:), POINTER :: centproj     !< rot.frame centr.3d->2d
     REAL, DIMENSION(:,:,:), POINTER :: cos1
     REAL, DIMENSION(:,:,:), POINTER :: sin1
     REAL                            :: c_p          !< spec. heat cap.
     REAL                            :: gacc         !< grav. acceleration
     REAL                            :: gamma        !< rat. spec. heats
     INTEGER                         :: use_envelope !< enable vicosity envelope
     ! --------------------------------------------------------------------- !
     REAL                            :: Q            !> shearing parameter
     REAL, DIMENSION(:,:,:), POINTER   :: height,h_ext !< disk scale height
     REAL, DIMENSION(:,:), POINTER   :: Sigma_dust   !< dust surface desnity
     REAL, DIMENSION(:,:), POINTER   :: invheight2   !< 1/h**2
     REAL, DIMENSION(:,:), POINTER   :: Qstar        !< stellar heating source
     REAL, DIMENSION(:,:), POINTER   :: H_tau1       !< z(tau=1)
     REAL                            :: tau_c        !< cooling time
     REAL, DIMENSION(:,:), POINTER   :: rescale      !< convert tau_z to tau_s
     !> angle between disk surface and line of sight to the star
     REAL, DIMENSION(:,:), POINTER   :: flaring_angle
     !> inverse 3d distance to heating central object
     REAL, DIMENSION(:,:), POINTER   :: invdr_3d
     !> vector from heating central object to mesh points, 3dim
     REAL, DIMENSION(:,:,:), POINTER :: Distance_3d
     REAL, DIMENSION(:,:,:), POINTER :: D_3d_norm    !< nornmalized Distance_3d
     !> normal vector of the disk surface in curvilinear coordinates
     REAL, DIMENSION(:,:,:), POINTER :: n
     REAL, DIMENSION(:,:,:), POINTER :: n_cart       !< cartesian n
     REAL, DIMENSION(:,:,:), POINTER   :: dynvis, &    !< dynamic viscosity
                                        kinvis, &    !< kinematic viscosity
                                        bulkvis, &   !< bulk viscosity
                                        envelope
     REAL, DIMENSION(:,:,:), POINTER :: pot          !< gravitational potential
     !> components of the stress tensor
     REAL, DIMENSION(:,:,:), POINTER   :: btxx,btyy,&
          btzz,btxy,btxz,btyz,tmp,tmp2,tmp3
     REAL, DIMENSION(:,:,:), POINTER :: tmp4         !<    temp array
     REAL, DIMENSION(:,:), POINTER   :: Sxx,Syy,&
          Szz,Sxy,Sxz,Syz
     !> source terms of sgs module
     REAL, DIMENSION(:,:), POINTER   :: diff,rhoeps,&
                                        sigma
     REAL, DIMENSION(:,:,:), POINTER :: ftxx,ftyy,&
          ftzz,ftxy,ftxz,ftyz
     REAL, DIMENSION(:,:), POINTER   :: delta        !< half-width of the filter
     REAL, DIMENSION(:,:,:), POINTER :: cent         !< rot. frame centrifugal
     REAL, DIMENSION(:,:,:), POINTER :: init_pvar    !< pvar state from init
     REAL, DIMENSION(:,:,:), POINTER  :: ptr3        !< 3d pointer
#ifdef HAVE_FFTW
     REAL, DIMENSION(:,:,:), POINTER :: fk, rand     !< forcing (fourier)
     !> \name
     !!#### forcing
     REAL                            :: L,K0,F0,CHI,&
                                        T,invsqrtN,&
                                        stoptime
     TYPE(C_PTR)                     :: plan_r2r     !< fftw plan (real to real)
     REAL(C_DOUBLE), POINTER         :: temp_c(:,:)  !< fftw temp storage
     REAL(C_DOUBLE), POINTER         :: Ftemp_c(:,:)
#endif
  CONTAINS
    PROCEDURE :: InitSources_all
    PROCEDURE (InitiateSources), DEFERRED :: InitiateSources
!    PROCEDURE :: InitSources_common
    GENERIC :: InitSources => InitiateSources!, InitSources_common
    PROCEDURE :: InfoSources_all
    PROCEDURE (InfoSources), DEFERRED :: InfoSources
    PROCEDURE :: MallocSources
    PROCEDURE :: CloseSources_all
    PROCEDURE (Close_Sources), DEFERRED :: Close_Sources
    GENERIC :: CloseSources => Close_Sources!, CloseSources_common
    PROCEDURE :: GeometricalSources
    PROCEDURE :: ExternalSources_all
    PROCEDURE (ExternalSources), DEFERRED :: ExternalSources
    PROCEDURE :: CalcTimestep_all
    PROCEDURE (CalcTimestep), DEFERRED :: CalcTimestep
!    PROCEDURE :: GetSourcesPointer
  END TYPE sources_base
  ABSTRACT INTERFACE
    SUBROUTINE InitiateSources(this,Mesh,Physics,Fluxes, config,IO) !Timedisc,config,IO)
      IMPORT Mesh_base, Physics_base,Fluxes_base, sources_base, Dict_TYP
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(Sources_base) :: this
      CLASS(Mesh_base)    :: Mesh
      CLASS(Physics_base) :: Physics
      CLASS(Fluxes_base)  :: Fluxes
    !  CLASS(Timedisc_base) :: Timedisc
      TYPE(Dict_TYP),POINTER :: config,IO
    !------------------------------------------------------------------------!
      INTENT(IN)        :: Mesh,Physics,Fluxes
    END SUBROUTINE
    SUBROUTINE InfoSources(this,Mesh)
      IMPORT sources_base, Mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(Sources_base),INTENT(IN) :: this
      CLASS(Mesh_base),INTENT(IN)    :: Mesh
    END SUBROUTINE
    SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt)
      IMPORT Sources_base, Mesh_base, Physics_base, Fluxes_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(Sources_base),INTENT(INOUT) :: this
      CLASS(Mesh_base),INTENT(IN)         :: Mesh
      CLASS(Physics_base),INTENT(INOUT)   :: Physics
      CLASS(Fluxes_base),INTENT(IN)       :: Fluxes
      REAL,INTENT(IN)                     :: time
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum) &
                      :: pvar,cvar
      REAL              :: dt
      !------------------------------------------------------------------------!
      INTENT(IN)        :: pvar,cvar
      INTENT(OUT)       :: dt
    END SUBROUTINE
    SUBROUTINE Close_Sources(this)
      IMPORT Sources_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(Sources_base) :: this
      !------------------------------------------------------------------------!
      INTENT(INOUT)     :: this
    END SUBROUTINE
    SUBROUTINE ExternalSources(this,Mesh,Physics,Fluxes,time,pvar,cvar,sterm)
      IMPORT Sources_base, Mesh_base, Physics_base, Fluxes_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(Sources_base),INTENT(INOUT) :: this
      CLASS(Mesh_base),INTENT(IN)         :: Mesh
      CLASS(Physics_base),INTENT(INOUT)   :: Physics
      CLASS(Fluxes_base),INTENT(IN)       :: Fluxes
      REAL,INTENT(IN)                     :: time
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
      !------------------------------------------------------------------------!
      INTENT(IN)        :: cvar,pvar
      INTENT(OUT)       :: sterm
    END SUBROUTINE
  END INTERFACE
  ! tempory storage for source terms
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: temp_sterm
  ! flags for source terms
  INTEGER, PARAMETER :: GRAVITY          = 1
  INTEGER, PARAMETER :: DISK_THOMSON     = 2
  INTEGER, PARAMETER :: VISCOSITY        = 3
  INTEGER, PARAMETER :: C_ACCEL          = 4
  INTEGER, PARAMETER :: COOLING          = 5
  INTEGER, PARAMETER :: ROTATING_FRAME   = 20
  INTEGER, PARAMETER :: SGS              = 23
  INTEGER, PARAMETER :: DISK_COOLING     = 24
  INTEGER, PARAMETER :: WAVE_DAMPING     = 25
  INTEGER, PARAMETER :: FORCING          = 26
  INTEGER, PARAMETER :: PLANET_HEATING   = 27
  INTEGER, PARAMETER :: PLANET_COOLING   = 28
  INTEGER, PARAMETER :: STELLAR_HEATING  = 29
  INTEGER, PARAMETER :: SHEARBOX         = 30
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       sources_base, &
       ! constants
       VISCOSITY!, & !GRAVITY, DISK_THOMSON, C_ACCEL, COOLING, &
       !ROTATING_FRAME, SGS, DISK_COOLING, WAVE_DAMPING, FORCING, &
       !POINTMASS, POINTMASS_BINARY, MONOPOL, &
       !NEWTON, WIITA, STELLAR_HEATING, &
      ! MOLECULAR, ALPHA, BETA, PRINGLE, ALPHA_ALT!, &
       !MULTIGRID, SPECTRAL, POTENTIAL, &
       !RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL, &
       !SPHERMULTEXPAN, CYLINMULTEXPAN, &
       !PLANET_HEATING, PLANET_COOLING, &
       !MRN_DUST, DUST_100, DUST_1000,&
       !GRAY,GAMMIE,GAMMIE_SB, &
       !SHEARBOX
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_all(list,Mesh,Fluxes,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),INTENT(IN) :: list
    CLASS(Mesh_base),INTENT(IN)    :: Mesh
    CLASS(Fluxes_base),INTENT(IN)  :: Fluxes
    CLASS(Physics_base),INTENT(INOUT) :: Physics
!    CLASS(Timedisc_base) :: Timedisc
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: sp
    TYPE(Dict_TYP),POINTER :: dir,src,IOsrc,gsrc => null(),gdir => null()
    INTEGER           :: update_disk_height = 0
    !------------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL list%Error("InitSources","physics and/or mesh module uninitialized")
    ! allocate common memory for all sources
    IF (.NOT.ALLOCATED(temp_sterm)) THEN
       CALL list%MallocSources(Mesh,Physics)
    END IF
    dir => config
    DO WHILE(ASSOCIATED(dir))
      NULLIFY(IOsrc)
      IF(HasChild(dir)) THEN
        src => GetChild(dir)
        CALL GetAttr(src, "stype", stype)

        SELECT CASE(stype)
!        CASE(GRAVITY)
!           ! skip initialization of gravity modules here and initialize them
!           ! at the end to make sure gravity is the first source term in the list
!           gdir => dir
!           gsrc => src
!        CASE(DISK_THOMSON)
!           ! radiational acceleration due to Thomson scattering
!           ! of accretion disk radiation
!           CALL InitSources_diskthomson(list,Mesh,Physics,src,IOsrc)
        CASE(VISCOSITY)
           ! viscous diffusion and heating
           CALL list%InitSources(Mesh,Physics,Fluxes,src,IOsrc)
!        CASE(C_ACCEL)
!           ! constant acceleration in x- and y-direction
!           CALL InitSources_c_accel(list,Mesh,Physics,src)
!        CASE(WAVE_DAMPING)
!           ! wave damping for planet eu experiments
!           CALL InitSources_wave_damping(list,Mesh,Physics,src)
!        CASE(COOLING)
!           ! simple cooling function
!           CALL InitSources_cooling(list,Mesh,Physics,src)
!        CASE(ROTATING_FRAME)
!           ! inertial forces due to rotating reference frame
!           CALL InitSources_rotframe(list,Mesh,Physics,src)
!        CASE(SGS)
!           CALL InitSources_sgs(list,Mesh,Physics,Fluxes,src,IOsrc)
!        CASE(DISK_COOLING)
!           CALL InitSources_diskcooling(list,Mesh,Physics,Timedisc,src,IOsrc)
!           IF (GetType(list%cooling).EQ.GRAY) update_disk_height = 1
!        CASE(FORCING)
!           CALL InitSources_forcing(list,Mesh,Physics,Fluxes,src,IOsrc)
!        CASE(PLANET_COOLING)
!          CALL InitSources_planetcooling(list,Mesh,Physics,src,IOsrc)
!        CASE(PLANET_HEATING)
!          CALL InitSources_planetheating(list,Mesh,Physics,src,IOsrc)
!        CASE(STELLAR_HEATING)
!           ! heating by the central star
!           CALL InitSources_stellarheating(list,Mesh,Physics,Timedisc,src,IOsrc)
!           update_disk_height = 1
!        CASE(SHEARBOX)
!           ! skip initialization and run it after the gravity but before the
!           ! other parts
!          CALL InitSources_shearbox(list,Mesh,Physics,src,IOsrc)
        CASE DEFAULT
           CALL list%Error("InitSources", "unknown source term")
        END SELECT
        IF(ASSOCIATED(IOsrc)) &
          CALL SetAttr(IO, GetKey(dir), IOsrc)
      END IF
      dir => GetNext(dir)
    END DO
!
!    ! finally initialize gravity
!    IF(ASSOCIATED(gsrc)) THEN
!       NULLIFY(IOsrc)
!       CALL SetAttr(gsrc,"update_disk_height", update_disk_height)
!       CALL InitGravity(list,Mesh,Fluxes,Physics,Timedisc%Boundary,GRAVITY,gsrc,IOsrc)
!       IF(ASSOCIATED(IOsrc)) &
!         CALL SetAttr(IO, GetKey(gdir), IOsrc)
!    END IF
!
!!    IF(ASSOCIATED(ssrc)) THEN
!!       NULLIFY(IOsrc)
!!       CALL InitSources_shearbox(list,Mesh,Physics,src,IOsrc)
!!       IF(ASSOCIATED(IOsrc)) &
!!         CALL SetAttr(IO, GetKey(sdir), IOsrc)
!!    END IF

    ! print some information
    CALL list%InfoSources(Mesh)
  END SUBROUTINE InitSources_all


  SUBROUTINE MallocSources(list,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),INTENT(IN) :: list !, POINTER :: list
    CLASS(Mesh_base),INTENT(IN)    :: Mesh
    CLASS(Physics_base),INTENT(IN) :: Physics
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    ! temporay storage
    ALLOCATE(temp_sterm(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum), &
         STAT=err)
    IF (err.NE.0) CALL list%Error("MallocSources_generic", "Unable allocate memory!")
  END SUBROUTINE MallocSources


  SUBROUTINE InfoSources_all(list, Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),Target,INTENT(IN) :: list !, POINTER :: list,srcptr
    CLASS(Sources_base),POINTER    :: srcptr
    CLASS(Mesh_base),INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh
    !------------------------------------------------------------------------!
    ! go through all source terms in the list
    srcptr => list
    ! skip gravity
    IF (srcptr%GetType().EQ.GRAVITY) srcptr => srcptr%next
    DO WHILE (ASSOCIATED(srcptr))
       CALL srcptr%Info(" SOURCES--> source term:       " // srcptr%GetName())
       SELECT CASE(srcptr%GetType())
!       CASE(DISK_THOMSON)
!          CALL InfoSources_diskthomson(srcptr)
       CASE(VISCOSITY)
          CALL srcptr%InfoSources(Mesh)
!       CASE(ROTATING_FRAME)
!          CALL InfoSources_rotframe(srcptr,Mesh)
!       CASE(DISK_COOLING)
!          CALL InfoSources_diskcooling(srcptr)
!       CASE(FORCING)
!          CALL InfoSources_forcing(srcptr)
!       CASE(SHEARBOX)
!          CALL InfoSources_shearbox(srcptr)
!       CASE DEFAULT
!          ! do nothing
       END SELECT
       srcptr => srcptr%next
    END DO
  END SUBROUTINE InfoSources_all


  SUBROUTINE GeometricalSources(this,Physics,Mesh,Fluxes,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_base), INTENT(IN)   :: this
    CLASS(Physics_base),INTENT(INOUT)  :: Physics
    CLASS(Mesh_base),INTENT(IN)     :: Mesh
    CLASS(Fluxes_base),INTENT(IN)   :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum) &
         :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! calculate geometrical sources depending on the integration rule
    SELECT CASE(Mesh%GetType())
    CASE(MIDPOINT)
       ! use center values for midpoint rule
       CALL Physics%GeometricalSources(Mesh,pvar,cvar,sterm)
  !  CASE(TRAPEZOIDAL)
  !     ! use reconstructed corner values for trapezoidal rule
  !     CALL GeometricalSources_physics(Physics,Mesh,Fluxes%prim,Fluxes%cons,sterm)
    END SELECT
  END SUBROUTINE GeometricalSources


  SUBROUTINE ExternalSources_all(this,Mesh,Fluxes,Physics,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),Target,INTENT(IN) :: this !, POINTER :: this
    CLASS(Mesh_base),INTENT(IN)    :: Mesh
    CLASS(Fluxes_base),INTENT(IN)  :: Fluxes
    CLASS(Physics_base),INTENT(INOUT) :: Physics
    REAL              :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: srcptr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: time,dt,pvar,cvar
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! reset sterm
    sterm(:,:,:,:) = 0.
    ! go through all source terms in the list
    srcptr => this
    DO WHILE (ASSOCIATED(srcptr))
       ! call specific subroutine

       SELECT CASE(srcptr%GetType())
!       CASE(GRAVITY)
!          CALL GravitySources(srcptr,Mesh,Physics,Fluxes,time,dt,pvar,cvar,temp_sterm)
!       CASE(DISK_THOMSON)
!          CALL ExternalSources_diskthomson(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
       CASE(VISCOSITY)
          CALL srcptr%ExternalSources(Mesh,Physics,Fluxes,time,pvar,cvar,temp_sterm)
!       CASE(C_ACCEL)
!          CALL ExternalSources_c_accel(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
!       CASE(WAVE_DAMPING)
!          CALL ExternalSources_wave_damping(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
!       CASE(COOLING)
!          CALL ExternalSources_cooling(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
!       CASE(ROTATING_FRAME)
!          CALL ExternalSources_rotframe(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
!       CASE(SGS)
!          CALL ExternalSources_sgs(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
!       CASE(DISK_COOLING)
!          CALL ExternalSources_diskcooling(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
!       CASE(STELLAR_HEATING)
!          CALL ExternalSources_stellarheating(srcptr,Mesh,Physics,Fluxes,time,pvar,cvar,temp_sterm)
!       CASE(FORCING)
!          CALL ExternalSources_forcing(srcptr,Mesh,Physics,time,dt,pvar,cvar,temp_sterm)
!       CASE(PLANET_COOLING)
!          CALL ExternalSources_planetcooling(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
!       CASE(PLANET_HEATING)
!          CALL ExternalSources_planetheating(srcptr,Mesh,Physics,time,pvar,cvar,temp_sterm)
!       CASE(SHEARBOX)
!          CALL ExternalSources_shearbox(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
!       CASE DEFAULT
!          CALL Error(srcptr,"ExternalSources", "unknown source term")
       END SELECT

       ! add to the sources
       sterm(:,:,:,:) = sterm(:,:,:,:) + temp_sterm(:,:,:,:)
       ! next source term
       srcptr => srcptr%next
    END DO
    ! reset ghost cell data
    sterm(Mesh%IGMIN:Mesh%IMIN-Mesh%IP1,:,:,:) = 0.0
    sterm(Mesh%IMAX+Mesh%IP1:Mesh%IGMAX,:,:,:) = 0.0
    sterm(:,Mesh%JGMIN:Mesh%JMIN-Mesh%JP1,:,:) = 0.0
    sterm(:,Mesh%JMAX+Mesh%JP1:Mesh%JGMAX,:,:) = 0.0
    sterm(:,:,Mesh%KGMIN:Mesh%KMIN-Mesh%KP1,:) = 0.0
    sterm(:,:,Mesh%KMAX+Mesh%KP1:Mesh%KGMAX,:) = 0.0
  END SUBROUTINE ExternalSources_all


  SUBROUTINE CalcTimestep_all(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt,dtcause)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),Target,INTENT(IN) :: this !, POINTER :: this
    CLASS(Mesh_base),INTENT(IN)    :: Mesh
    CLASS(Physics_base),INTENT(INOUT) :: Physics
    CLASS(Fluxes_base),INTENT(IN)  :: Fluxes
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL              :: dt
    INTEGER           :: dtcause
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: srcptr
    REAL              :: dt_new
    INTEGER           :: hc
    !------------------------------------------------------------------------!
    INTENT(IN)        :: time,pvar,cvar
    INTENT(INOUT)     :: dt,dtcause
    !------------------------------------------------------------------------!

    ! go through all source terms in the list
    srcptr => this
    DO WHILE(ASSOCIATED(srcptr))
       ! call specific subroutine
       SELECT CASE(srcptr%GetType())
       CASE(DISK_THOMSON,C_ACCEL,ROTATING_FRAME,WAVE_DAMPING,GRAVITY,&
         PLANET_HEATING,SHEARBOX)
          ! do nothing
          dt_new = dt
       CASE(VISCOSITY)
          CALL srcptr%CalcTimestep(Mesh,Physics,Fluxes,time,pvar,cvar,dt_new)
!       CASE(COOLING)
!          CALL CalcTimestep_cooling(srcptr,Mesh,Physics,time,pvar,dt_new)
!       CASE(DISK_COOLING)
!          CALL CalcTimestep_diskcooling(srcptr,Mesh,Physics,Fluxes,time,pvar,dt_new)
!       CASE(STELLAR_HEATING)
!           CALL CalcTimestep_StellarHeating(srcptr,Mesh,Physics,Fluxes,time,pvar,dt_new)
!       CASE(SGS)
!          CALL CalcTimestep_sgs(srcptr,Mesh,Physics,time,pvar,cvar,dt_new)
!       CASE(FORCING)
!          CALL CalcTimestep_forcing(srcptr,Mesh,Physics,pvar,cvar,dt_new)
!       CASE(PLANET_COOLING)
!          CALL CalcTimestep_planetcooling(srcptr,Mesh,Physics,time,pvar,dt_new)
!       CASE DEFAULT
!          CALL Error(srcptr,"CalcTimestep", "unknown source term")
       END SELECT
       ! who was it?
       IF (dt_new .LT. dt) dtcause=srcptr%GetType()
       dt = MIN(dt,dt_new)
       ! next source term
       srcptr => srcptr%next
    END DO
  END SUBROUTINE CalcTimestep_all


  SUBROUTINE CloseSources_all(this,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_base),target,INTENT(IN) :: this !, POINTER :: this
    CLASS(Fluxes_base),INTENT(IN)  :: Fluxes
    !------------------------------------------------------------------------!
    CLASS(Sources_base), POINTER :: srcptr
    !------------------------------------------------------------------------!
    ! call deallocation procedures for all source terms
    DO
       srcptr => this
       IF (.NOT.ASSOCIATED(srcptr)) EXIT
       srcptr => srcptr%next
       IF (.NOT.srcptr%Initialized()) &
            CALL srcptr%Error("CloseSources","not initialized")
!        IF (.NOT.ASSOCIATED(this)) EXIT
!        srcptr => this%next
!        IF (.NOT.Initialized(srcptr)) &
!             CALL this%Error("CloseSources","not initialized")
       ! call specific deconstructor
       SELECT CASE(srcptr%GetType())
!       CASE(GRAVITY)
!          CALL CloseGravity(srcptr)
!       CASE(DISK_THOMSON)
!          CALL CloseSources_diskthomson(srcptr)
       CASE(VISCOSITY)
          CALL srcptr%CloseSources()
!       CASE(C_ACCEL)
!          CALL CloseSources_c_accel(srcptr,Fluxes)
!       CASE(WAVE_DAMPING)
!          CALL CloseSources_wave_damping(srcptr,Fluxes)
!       CASE(COOLING)
!          CALL CloseSources_cooling(srcptr)
!       CASE(ROTATING_FRAME)
!          CALL CloseSources_rotframe(srcptr)
!       CASE(SGS)
!          CALL CloseSources_sgs(srcptr)
!       CASE(DISK_COOLING)
!          CALL CloseSources_diskcooling(srcptr)
!       CASE(STELLAR_HEATING)
!          CALL CloseSources_stellarheating(srcptr)
!       CASE(FORCING)
!          CALL CloseSources_forcing(srcptr)
!       CASE(PLANET_COOLING)
!          CALL CloseSources_planetcooling(srcptr)
!       CASE(PLANET_HEATING)
!          CALL CloseSources_planetheating(srcptr)
!       CASE(SHEARBOX)
!          CALL CloseSources_shearbox(srcptr)
       END SELECT
       ! deallocate source term structure
       DEALLOCATE(srcptr)
    END DO
    ! release temporary storage
    DEALLOCATE(temp_sterm)
  END SUBROUTINE CloseSources_all


END MODULE sources_base_mod
