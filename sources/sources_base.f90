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
  USE physics_base_mod
  USE fluxes_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(logging_base) :: sources_base
     !> \name Variables
     CLASS(sources_base), POINTER    :: next => null() !< next source in list
     !TYPE(Gravity_TYP), POINTER      :: glist => null()!< gravity list
     !TYPE(Common_TYP)                :: viscosity    !< molecular,alpha,beta
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
     REAL, DIMENSION(:,:,:), POINTER :: accel,accart !< acceleration
     REAL, DIMENSION(:,:,:), POINTER :: bcposvec,bccart !< position vector
     REAL, DIMENSION(:,:), POINTER   :: radius       !< distance to origin
     REAL, DIMENSION(:,:), POINTER   :: invr         !< 1./radius
     REAL, DIMENSION(:,:), POINTER   :: cs           !< speed of sound
     REAL, DIMENSION(:,:), POINTER   :: omega        !< angular velocity
     REAL, DIMENSION(:,:,:), POINTER :: omega2       !< Omega Kepler squared
     REAL, DIMENSION(:,:,:), POINTER :: gxr3         !< = GN*x/radius**3
     REAL, DIMENSION(:,:), POINTER   :: cellmass     !< rho*dV
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
     REAL, DIMENSION(:,:), POINTER   :: height,h_ext !< disk scale height
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
     REAL, DIMENSION(:,:), POINTER   :: dynvis, &    !< dynamic viscosity
                                        kinvis, &    !< kinematic viscosity
                                        bulkvis, &   !< bulk viscosity
                                        envelope
     REAL, DIMENSION(:,:,:), POINTER :: pot          !< gravitational potential
     !> components of the stress tensor
     REAL, DIMENSION(:,:), POINTER   :: btxx,btyy,&
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
!  CONTAINS
!    PROCEDURE :: InitSources
!    PROCEDURE :: MallocSources
!    PROCEDURE :: CloseSources
!    PROCEDURE :: GeometricalSources
!    PROCEDURE :: ExternalSources
!    PROCEDURE :: CalcTimestep
!    PROCEDURE :: GetSourcesPointer
  END TYPE sources_base
  ! tempory storage for source terms
!  REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: temp_sterm
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
       sources_base!, &
       ! constants
       !GRAVITY, DISK_THOMSON, VISCOSITY, C_ACCEL, COOLING, &
       !ROTATING_FRAME, SGS, DISK_COOLING, WAVE_DAMPING, FORCING, &
       !POINTMASS, POINTMASS_BINARY, MONOPOL, &
       !NEWTON, WIITA, STELLAR_HEATING, &
       !MOLECULAR, ALPHA, BETA, PRINGLE, ALPHA_ALT, &
       !MULTIGRID, SPECTRAL, POTENTIAL, &
       !RED_BLACK_GAUSS_SEIDEL,BLOCK_GAUSS_SEIDEL,GAUSS_SEIDEL, &
       !SPHERMULTEXPAN, CYLINMULTEXPAN, &
       !PLANET_HEATING, PLANET_COOLING, &
       !MRN_DUST, DUST_100, DUST_1000,&
       !GRAY,GAMMIE,GAMMIE_SB, &
       !SHEARBOX
  !--------------------------------------------------------------------------!

!CONTAINS

!  SUBROUTINE InitSources(list,Mesh,Fluxes,Physics,Timedisc,config,IO)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: list
!    TYPE(Mesh_TYP)    :: Mesh
!    TYPE(Fluxes_TYP)  :: Fluxes
!    TYPE(Physics_TYP) :: Physics
!    TYPE(Timedisc_TYP) :: Timedisc
!    TYPE(Dict_TYP),POINTER :: config,IO
!    INTEGER           :: stype
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: sp
!    TYPE(Dict_TYP),POINTER :: dir,src,IOsrc,gsrc => null(),gdir => null()
!    INTEGER           :: update_disk_height = 0
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Fluxes
!    INTENT(INOUT)     :: Timedisc,Physics
!    !------------------------------------------------------------------------!
!    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
!         CALL Error(list,"InitSources","physics and/or mesh module uninitialized")
!    ! allocate common memory for all sources
!    IF (.NOT.ALLOCATED(temp_sterm)) THEN
!       CALL MallocSources(list,Mesh,Physics)
!    END IF
!    dir => config
!    DO WHILE(ASSOCIATED(dir))
!      NULLIFY(IOsrc)
!
!      IF(HasChild(dir)) THEN
!        src => GetChild(dir)
!        CALL GetAttr(src, "stype", stype)
!
!        SELECT CASE(stype)
!        CASE(GRAVITY)
!           ! skip initialization of gravity modules here and initialize them
!           ! at the end to make sure gravity is the first source term in the list
!           gdir => dir
!           gsrc => src
!        CASE(DISK_THOMSON)
!           ! radiational acceleration due to Thomson scattering
!           ! of accretion disk radiation
!           CALL InitSources_diskthomson(list,Mesh,Physics,src,IOsrc)
!        CASE(VISCOSITY)
!           ! viscous diffusion and heating
!           CALL InitSources_viscosity(list,Mesh,Physics,Fluxes,Timedisc,src,IOsrc)
!           IF (GetType(list%viscosity).EQ.ALPHA_ALT) update_disk_height = 1
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
!        CASE DEFAULT
!           CALL Error(list,"InitSources", "unknown source term")
!        END SELECT
!        IF(ASSOCIATED(IOsrc)) &
!          CALL SetAttr(IO, GetKey(dir), IOsrc)
!      END IF
!      dir => GetNext(dir)
!    END DO
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
!
!    ! print some information
!    CALL InfoSources(list, Mesh)
!  END SUBROUTINE InitSources
!
!
!  SUBROUTINE MallocSources(list,Mesh,Physics)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: list
!    TYPE(Mesh_TYP)    :: Mesh
!    TYPE(Physics_TYP) :: Physics
!    !------------------------------------------------------------------------!
!    INTEGER           :: err
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Physics
!    !------------------------------------------------------------------------!
!    ! temporay storage
!    ALLOCATE(temp_sterm(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
!         STAT=err)
!    IF (err.NE.0) CALL Error(list, "MallocSources_generic", "Unable allocate memory!")
!  END SUBROUTINE MallocSources
!
!
!  SUBROUTINE InfoSources(list, Mesh)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: list,srcptr
!    TYPE(Mesh_TYP)    :: Mesh
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh
!    !------------------------------------------------------------------------!
!    ! go through all source terms in the list
!    srcptr => list
!    ! skip gravity
!    IF (GetType(srcptr).EQ.GRAVITY) srcptr => srcptr%next
!    DO WHILE (ASSOCIATED(srcptr))
!       CALL Info(srcptr, " SOURCES--> source term:       " // GetName(srcptr))
!!CDIR IEXPAND
!       SELECT CASE(GetType(srcptr))
!       CASE(DISK_THOMSON)
!          CALL InfoSources_diskthomson(srcptr)
!       CASE(VISCOSITY)
!          CALL InfoSources_viscosity(srcptr)
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
!       END SELECT
!       srcptr => srcptr%next
!    END DO
!  END SUBROUTINE InfoSources
!
!
!  SUBROUTINE GeometricalSources(Physics,Mesh,Fluxes,pvar,cvar,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP)  :: Physics
!    TYPE(Mesh_TYP)     :: Mesh
!    TYPE(Fluxes_TYP)   :: Fluxes
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
!         :: pvar,cvar,sterm
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Fluxes,pvar,cvar
!    INTENT(INOUT)     :: Physics
!    INTENT(OUT)       :: sterm
!    !------------------------------------------------------------------------!
!    ! calculate geometrical sources depending on the integration rule
!!CDIR IEXPAND
!    SELECT CASE(GetType(Mesh))
!    CASE(MIDPOINT)
!       ! use center values for midpoint rule
!       CALL GeometricalSources_physics(Physics,Mesh,pvar,cvar,sterm)
!    CASE(TRAPEZOIDAL)
!       ! use reconstructed corner values for trapezoidal rule
!       CALL GeometricalSources_physics(Physics,Mesh,Fluxes%prim,Fluxes%cons,sterm)
!    END SELECT
!  END SUBROUTINE GeometricalSources
!
!
!  SUBROUTINE ExternalSources(this,Mesh,Fluxes,Physics,time,dt,pvar,cvar,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    TYPE(Fluxes_TYP)  :: Fluxes
!    TYPE(Physics_TYP) :: Physics
!    REAL              :: time,dt
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
!                      :: cvar,pvar,sterm
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: srcptr
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Fluxes,time,dt,pvar,cvar
!    INTENT(INOUT)     :: Physics
!    INTENT(OUT)       :: sterm
!    !------------------------------------------------------------------------!
!    ! reset sterm
!    sterm(:,:,:) = 0.
!    ! go through all source terms in the list
!    srcptr => this
!    DO WHILE (ASSOCIATED(srcptr))
!       ! call specific subroutine
!
!!CDIR IEXPAND
!       SELECT CASE(GetType(srcptr))
!       CASE(GRAVITY)
!          CALL GravitySources(srcptr,Mesh,Physics,Fluxes,time,dt,pvar,cvar,temp_sterm)
!       CASE(DISK_THOMSON)
!          CALL ExternalSources_diskthomson(srcptr,Mesh,Physics,pvar,cvar,temp_sterm)
!       CASE(VISCOSITY)
!          CALL ExternalSources_viscosity(srcptr,Mesh,Physics,Fluxes,time,pvar,cvar,temp_sterm)
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
!       END SELECT
!
!       ! add to the sources
!       sterm(:,:,:) = sterm(:,:,:) + temp_sterm(:,:,:)
!       ! next source term
!       srcptr => srcptr%next
!    END DO
!    ! reset ghost cell data
!    sterm(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
!    sterm(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
!    sterm(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
!    sterm(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
!  END SUBROUTINE ExternalSources
!
!
!  SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt,dtcause)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    TYPE(Physics_TYP) :: Physics
!    TYPE(Fluxes_TYP)  :: Fluxes
!    REAL              :: time
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
!                      :: pvar,cvar
!    REAL              :: dt
!    INTEGER           :: dtcause
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: srcptr
!    REAL              :: dt_new
!    INTEGER           :: hc
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Fluxes,time,pvar,cvar
!    INTENT(INOUT)     :: dt,dtcause,Physics
!    !------------------------------------------------------------------------!
!
!    ! go through all source terms in the list
!    srcptr => this
!    DO WHILE(ASSOCIATED(srcptr))
!       ! call specific subroutine
!!CDIR IEXPAND
!       SELECT CASE(GetType(srcptr))
!       CASE(DISK_THOMSON,C_ACCEL,ROTATING_FRAME,WAVE_DAMPING,GRAVITY,&
!         PLANET_HEATING,SHEARBOX)
!          ! do nothing
!          dt_new = dt
!       CASE(VISCOSITY)
!          CALL CalcTimestep_viscosity(srcptr,Mesh,Physics,Fluxes,time,pvar,cvar,dt_new)
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
!       END SELECT
!       ! who was it?
!       IF (dt_new .LT. dt) dtcause=GetType(srcptr)
!       dt = MIN(dt,dt_new)
!       ! next source term
!       srcptr => srcptr%next
!    END DO
!  END SUBROUTINE CalcTimestep
!
!
!  SUBROUTINE CloseSources(this,Fluxes)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: this
!    TYPE(Fluxes_TYP)  :: Fluxes
!    !------------------------------------------------------------------------!
!    TYPE(Sources_TYP), POINTER :: srcptr
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Fluxes
!    !------------------------------------------------------------------------!
!    ! call deallocation procedures for all source terms
!    DO
!       srcptr => this
!       IF (.NOT.ASSOCIATED(srcptr)) EXIT
!       this => srcptr%next
!       IF (.NOT.Initialized(srcptr)) &
!            CALL Error(this,"CloseSources","not initialized")
!       ! call specific deconstructor
!!CDIR IEXPAND
!       SELECT CASE(GetType(srcptr))
!       CASE(GRAVITY)
!          CALL CloseGravity(srcptr)
!       CASE(DISK_THOMSON)
!          CALL CloseSources_diskthomson(srcptr)
!       CASE(VISCOSITY)
!          CALL CloseSources_viscosity(srcptr)
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
!       END SELECT
!       ! deallocate source term structure
!       DEALLOCATE(srcptr)
!    END DO
!    ! release temporary storage
!    DEALLOCATE(temp_sterm)
!  END SUBROUTINE CloseSources

END MODULE sources_base_mod
