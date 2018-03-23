!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007 - 2016                                                 #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!> \addtogroup physics
!! - general physics settings
!! \key{problem,INTEGER,advection problem
!!      (see \link physics_generic \endlink for a list of currently supported
!!       advection problems)}
!! \key{units,INTEGER,unit system
!!      (see \link constants_generic \endlink for a list of currently supported
!!       unit systems)}
!! \key{mu,REAL,mean molecular weight (default is for air
!!      at normal conditions),0.029}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief generic module for the advection problem
!!
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_base_mod
  USE constants_generic_mod
  USE logging_base_mod
  USE mesh_base_mod
!  USE sources_base
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, ABSTRACT, EXTENDS(logging_base) :: physics_base
     !> \name Variables
!     CLASS(logging_base)    :: advproblem            !< advection problem
     CLASS(constants_base), ALLOCATABLE :: constants             !< physical constants
!     TYPE(sources_base, POINTER &
!                         :: sources => null()     !< list of source terms
     REAL                :: gamma,&               !< ratio of spec. heats
                            time,&                !< simulation time
                            mu, &                 !< mean molecular weight
                            csiso, &              !< isothermal sound speed
                            eps                   !< softening length
     INTEGER             :: VNUM, &               !< number of variables
                            DIM, &                !< Dimension (1, 2 or 3)
                            DENSITY,PRESSURE, &
                            ENERGY,SGSPRESSURE, &
                            SGSENERGY, &
                            XVELOCITY,XMOMENTUM, &
                            YVELOCITY,YMOMENTUM,&
                            ZVELOCITY,ZMOMENTUM   !< array indicies for primitive and conservative variables
     LOGICAL             :: transformed_yvelocity !< .TRUE. if SubtractBackgroundVelocity was called before
                                                  !! .FALSE. otherwise
     LOGICAL             :: supports_absorbing    !< absorbing boundary conditions supported
                            !! \details .TRUE. if absorbing boundary conditions are supported by the physics module
     LOGICAL             :: supports_farfield     !< farfield boundary conditions supported
                            !! \details .TRUE. if farfield boundary conditions are supported by the physics module
     LOGICAL             :: advanced_wave_speeds  !< use Roe averages for min/max wave speed estimates
     CHARACTER(LEN=16), DIMENSION(:), POINTER &
                         :: pvarname,cvarname     !< names of variables
     REAL, DIMENSION(:,:,:), POINTER &
                         :: bccsound, &           !< bary centered speed of sound
                            bcradius, &           !< distance to the origin bary center values
                            divposvec, &          !< divergence of the position vector
                            bphi, &               !< bary centered constant gravitational potential
                            tmp,tmp1,tmp2,tmp3, &
                            tmp4,tmp5             !< temporary storage
     REAL, DIMENSION(:,:,:,:), POINTER &
                         :: fcsound, &            !< speed of sound faces
                            fradius, &            !< distance to the origin face values
                            bcposvec, &           !< curvilinear components of the position vector bary center values
                            w => NULL(), &        !< fargo bulk velocity
                            fphi, &               !< face centered constant gravitational potential
                            hy                    !< chy or fhy depending on reconstruction
     REAL, DIMENSION(:,:,:,:,:), POINTER &
                         :: fcent, &              !< centrifugal force
                            fposvec               !< curvilinear components of the position vector face values
!------------------------------------------------------------!
  CONTAINS
    PROCEDURE :: InitPhysics
    PROCEDURE (ExternalSources),              DEFERRED :: ExternalSources
    !------Convert2Primitve--------!
    PROCEDURE (Convert2Primitive_center),     DEFERRED :: Convert2Primitive_center
    PROCEDURE (Convert2Primitive_centsub),    DEFERRED :: Convert2Primitive_centsub
    PROCEDURE (Convert2Primitive_faces),      DEFERRED :: Convert2Primitive_faces
    PROCEDURE (Convert2Primitive_facesub),    DEFERRED :: Convert2Primitive_facesub
    GENERIC   :: Convert2Primitive => &
                   Convert2Primitive_center, &
                   Convert2Primitive_centsub, &
                   Convert2Primitive_faces, &
                   Convert2Primitive_facesub
    !------Convert2Conservative----!
    PROCEDURE (Convert2Conservative_center),  DEFERRED :: Convert2Conservative_center
    PROCEDURE (Convert2Conservative_centsub), DEFERRED :: Convert2Conservative_centsub
    PROCEDURE (Convert2Conservative_faces),   DEFERRED :: Convert2Conservative_faces
    PROCEDURE (Convert2Conservative_facesub), DEFERRED :: Convert2Conservative_facesub
    GENERIC   :: Convert2Conservative => &
                   Convert2Conservative_center, &
                   Convert2Conservative_centsub, &
                   Convert2Conservative_faces, &
                   Convert2Conservative_facesub
    !------Soundspeed Routines-----!
    PROCEDURE :: SetSoundSpeeds_center
    PROCEDURE :: SetSoundSpeeds_faces
    GENERIC   :: SetSoundSpeeds => &
                   SetSoundSpeeds_center, &
                   SetSoundSpeeds_faces
    PROCEDURE (UpdateSoundSpeed_center),     DEFERRED :: UpdateSoundSpeed_center
    PROCEDURE (UpdateSoundSpeed_faces),      DEFERRED :: UpdateSoundSpeed_faces
    GENERIC   :: UpdateSoundSpeed => &
                   UpdateSoundSpeed_center, &
                   UpdateSoundSpeed_faces
    !------Wavespeed Routines-----!
    PROCEDURE (CalcWaveSpeeds_center),        DEFERRED :: CalcWaveSpeeds_center
    PROCEDURE (CalcWaveSpeeds_faces),         DEFERRED :: CalcWaveSpeeds_faces
    GENERIC   :: CalculateWaveSpeeds => &
                   CalcWaveSpeeds_center, &
                   CalcWaveSpeeds_faces
    !------Flux Routines-----------!
    PROCEDURE (CalcFluxesX),                  DEFERRED :: CalcFluxesX
    GENERIC   :: CalculateFluxesX => CalcFluxesX
    PROCEDURE (CalcFluxesY),                  DEFERRED :: CalcFluxesY
    GENERIC   :: CalculateFluxesY => CalcFluxesY
    PROCEDURE (CalcFluxesZ),                  DEFERRED :: CalcFluxesZ
    GENERIC   :: CalculateFluxesZ => CalcFluxesZ
    !------Geometry Routines-------!
    PROCEDURE (GeometricalSources_center),    DEFERRED :: GeometricalSources_center!, GeometricalSources_faces
    GENERIC   :: GeometricalSources => &
                   GeometricalSources_center!, &
                   !GeometricalSources_faces
    PROCEDURE (ReflectionMasks),              DEFERRED :: ReflectionMasks

! these routines are only necessary for special boundaries, fluxes
!   PROCEDURE ::  CalculateCharSystemX
!   PROCEDURE ::  CalculateCharSystemY
!   PROCEDURE ::  CalculateBoundaryDataX
!   PROCEDURE ::  CalculateBoundaryDataY
!   PROCEDURE ::  CalculatePrim2RiemannX
!   PROCEDURE ::  CalculatePrim2RiemannY
!   PROCEDURE ::  CalculateRiemann2PrimX
!   PROCEDURE ::  CalculateRiemann2PrimY
!   PROCEDURE ::  CalculateStresses
!   PROCEDURE ::  ViscositySources
!   PROCEDURE ::  SGSSources
!   PROCEDURE ::  CalculateSGSTensor
   PROCEDURE (MASKS), DEFERRED :: AxisMasks
!   PROCEDURE :: GetSoundSpeed_adiabatic

    PROCEDURE :: FinalizePhysics
  END TYPE physics_base

  ABSTRACT INTERFACE
    PURE SUBROUTINE Convert2Primitive_center(this,Mesh,cvar,pvar)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(IN)  :: cvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(OUT) :: pvar
    END SUBROUTINE
    PURE SUBROUTINE Convert2Primitive_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,cvar,pvar)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: i1,i2,j1,j2,k1,k2
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(IN)  :: cvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(OUT) :: pvar
    END SUBROUTINE
    PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,cons,prim)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: cons
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: prim
    END SUBROUTINE
    PURE SUBROUTINE Convert2Primitive_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,cons,prim)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: i1,i2,j1,j2,k1,k2
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: cons
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: prim
    END SUBROUTINE
    PURE SUBROUTINE Convert2Conservative_center(this,Mesh,pvar,cvar)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(IN)  :: pvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(OUT) :: cvar
    END SUBROUTINE
    PURE SUBROUTINE Convert2Conservative_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,pvar,cvar)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: i1,i2,j1,j2,k1,k2
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(IN)  :: pvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum), &
                           INTENT(OUT) :: cvar
    END SUBROUTINE
    PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,prim,cons)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: prim
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: cons
    END SUBROUTINE
    PURE SUBROUTINE Convert2Conservative_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,prim,cons)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: i1,i2,j1,j2,k1,k2
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: prim
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: cons
    END SUBROUTINE
    PURE SUBROUTINE ExternalSources(this,Mesh,accel,pvar,cvar,sterm)
      IMPORT physics_base, mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS), &
                           INTENT(IN)  :: accel
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                           INTENT(IN)  :: pvar,cvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                           INTENT(OUT)  :: sterm
    END SUBROUTINE
    PURE SUBROUTINE Masks(this,reflX,reflY,reflZ)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN)            :: this
      LOGICAL, DIMENSION(this%VNUM), INTENT(OUT) :: reflX,reflY,reflZ
    END SUBROUTINE
    PURE SUBROUTINE UpdateSoundSpeed_center(this,Mesh,pvar)
      IMPORT physics_base,mesh_base
      IMPLICIT NONE
      CLASS(physics_base),INTENT(INOUT) :: this
      CLASS(mesh_base),   INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
        INTENT(IN)                      :: pvar
    END SUBROUTINE
    PURE SUBROUTINE UpdateSoundSpeed_faces(this,Mesh,prim)
      IMPORT physics_base,mesh_base
      IMPLICIT NONE
      CLASS(physics_base),INTENT(INOUT) :: this
      CLASS(mesh_base),   INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                          INTENT(IN)    :: prim
    END SUBROUTINE
    PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar,minwav,maxwav)
      IMPORT physics_base,mesh_base
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                           INTENT(IN)    :: pvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS), &
                           INTENT(OUT)   :: minwav,maxwav
    END SUBROUTINE
    PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim,cons,minwav,maxwav)
      IMPORT physics_base,mesh_base
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                           INTENT(IN)    :: prim,cons
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS), &
                           INTENT(OUT)   :: minwav,maxwav
    END SUBROUTINE
    PURE SUBROUTINE CalcFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
      IMPORT physics_base,mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: nmin,nmax
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: prim,cons
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: xfluxes
    END SUBROUTINE
    PURE SUBROUTINE CalcFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
      IMPORT physics_base,mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: nmin,nmax
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: prim,cons
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: yfluxes
    END SUBROUTINE
    PURE SUBROUTINE CalcFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
      IMPORT physics_base,mesh_base
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: nmin,nmax
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(IN)  :: prim,cons
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum), &
                           INTENT(OUT) :: zfluxes
    END SUBROUTINE
    PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
      IMPORT physics_base,mesh_base
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                           INTENT(IN)    :: pvar,cvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                           INTENT(OUT)   :: sterm
    END SUBROUTINE
    PURE SUBROUTINE ReflectionMasks(this,reflX,reflY,reflZ)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN)  :: this
      LOGICAL, DIMENSION(this%VNUM), &
                           INTENT(OUT) :: reflX,reflY,reflZ
    END SUBROUTINE
  END INTERFACE
  !--------------------------------------------------------------------------!
  ! flags for advection problems
  INTEGER, PARAMETER :: EULER2D             = 1
  INTEGER, PARAMETER :: EULER2D_ISOTHERM    = 2
  INTEGER, PARAMETER :: EULER3D_ROTSYM      = 3
  INTEGER, PARAMETER :: EULER3D_ROTAMT      = 4
  INTEGER, PARAMETER :: EULER3D_ROTSYMSGS   = 5
  INTEGER, PARAMETER :: EULER2D_SGS         = 7
  INTEGER, PARAMETER :: EULER3D_ROTAMTSGS   = 8
  INTEGER, PARAMETER :: EULER2D_ISOIAMT     = 9
  INTEGER, PARAMETER :: EULER2D_IAMT        = 11
  INTEGER, PARAMETER :: EULER2D_IAMROT      = 12
  INTEGER, PARAMETER :: EULER2D_ISOIAMROT   = 13
  INTEGER, PARAMETER :: EULER3D_ISOTHERM    = 14
  INTEGER, PARAMETER :: EULER3D             = 15
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       physics_base, &
       ! constants - flags for  identification in dictionary by an integer
       EULER2D, EULER2D_ISOTHERM, EULER3D_ROTSYM, EULER3D_ROTAMT, &
       EULER2D_SGS, EULER3D_ROTSYMSGS, EULER3D_ROTAMTSGS, &
       SI, CGS, GEOMETRICAL, EULER2D_ISOIAMT, &
       EULER2D_IAMT, EULER2D_IAMROT, EULER2D_ISOIAMROT, &
       EULER3D_ISOTHERM, EULER3D, &
       ! methods - only elemental functions
       SetWaveSpeeds, &
       CalcWaveSpeeds, &
       MomentumSourcesX, &
       MomentumSourcesY, &
       MomentumSourcesZ
  !--------------------------------------------------------------------------!

CONTAINS

  !> Initialization for the base physical object
  !!
  !! - initialize constants and unit system
  !! - read out attributes from config dictrionary (getAttr(..))
  !! - allocation of arrays
  !! - specific tweaks (update soundspeed, etc.)
  !! - print infostring to terminal
  SUBROUTINE InitPhysics(this,Mesh,config,IO,problem,pname,vnum)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)        :: Mesh
    TYPE(Dict_TYP),POINTER &
                                       :: config, IO
    INTEGER                            :: problem,vnum
    CHARACTER(LEN=32)                  :: pname
    !------------------------------------------------------------------------!
    INTEGER                            :: units
    INTEGER                            :: err, valwrite
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: problem,vnum
    !------------------------------------------------------------------------!
    CALL this%InitLogging(problem,pname)

    ! check initialization of Mesh
    IF (.NOT.Mesh%Initialized()) &
         CALL this%Error("InitPhysics","mesh module uninitialized")

    ! units
    CALL GetAttr(config, "units", units, SI)
    CALL new_constants(this%constants, units)

    !CALL GetAttr(config, "problem", problem)

    ! ratio of specific heats
    CALL GetAttr(config, "gamma", this%gamma, 1.4)

    ! mean molecular weight
    CALL GetAttr(config, "mu", this%mu, 0.029)

    ! isothermal sound speed
    CALL GetAttr(config, "cs", this%csiso, 0.0)

    ! enable advanced wave speed estimates (computationally more expensive)
    ! uses Roe averages between cell boundaries
    CALL GetAttr(config, "advanced_wave_speeds", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       this%advanced_wave_speeds = .TRUE.
    ELSE
       this%advanced_wave_speeds = .FALSE.
    END IF

    ! softening parameter to smooth out singularity near center of rotation
    ! (only necessary, if it's inside the computational domain)
    ! set to 0.0 to disable
    ! the softening length is the product of this parameter and the
    ! size of the grid cell next to the center of rotation; thus a value larger
    ! than 1.0 leads to larger softening whereas smaller values will
    ! probably cause odd behaviour due to the 1/r singularity;
    ! if the minimal r on the computational domain is larger than
    ! the size of the associated grid cell, softening is disabled, because
    ! the center of rotation lies outside of the computational domain
    CALL GetAttr(config, "softening", this%eps, 1.0)

    this%vnum = vnum
    ! allocate memory for arrays common to all physics modules
    ALLOCATE(this%pvarname(this%vnum),this%cvarname(this%vnum),                                           &
             this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                 &
             this%tmp1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp3(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp4(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp5(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%bccsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),            &
             this%fcsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces), &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitPhysics", "Unable to allocate memory.")

    this%tmp(:,:,:)  = 0.
    this%tmp1(:,:,:) = 0.
    this%tmp2(:,:,:) = 0.
    this%tmp3(:,:,:) = 0.
    this%tmp4(:,:,:) = 0.
    this%tmp5(:,:,:) = 0.

    IF(this%csiso.GT.0.) THEN
      this%bccsound(:,:,:)  = this%csiso
      this%fcsound(:,:,:,:) = this%csiso
    ELSE
      this%bccsound(:,:,:)  = 0.
      this%fcsound(:,:,:,:) = 0.
    END IF

   ! enable/disable absorbing and farfield boundary conditions
   ! TODO Not yet tested for 3D
   SELECT CASE(problem)
    CASE(EULER2D,EULER2D_ISOTHERM,&
         EULER2D_ISOIAMT,EULER2D_IAMT,EULER2D_SGS, &
         EULER3D_ROTSYM,EULER3D_ROTSYMSGS,EULER3D_ROTAMT,EULER3D_ROTAMTSGS)
       this%supports_absorbing = .TRUE.
    CASE DEFAULT
       this%supports_absorbing = .FALSE.
    END SELECT
    SELECT CASE(problem)
    CASE(EULER2D,EULER2D_ISOTHERM,EULER2D_SGS, &
         EULER3D_ROTSYM,EULER3D_ROTSYMSGS,EULER3D_ROTAMT,EULER3D_ROTAMTSGS)
       this%supports_farfield  = .TRUE.
    CASE DEFAULT
       this%supports_farfield  = .FALSE.
    END SELECT

    ! reset source term pointer
    !NULLIFY(this%sources)

    ! check if output of sound speeds is requested
    CALL GetAttr(config, "output/bccsound", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "bccsound",&
                    this%bccsound(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    CALL GetAttr(config, "output/fcsound", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL Setattr(IO, "fcsound",&
                    this%fcsound(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:))

    this%time = -1.
    ! print some information
    CALL this%Info(" PHYSICS--> advection problem: " // TRIM(this%GetName()))

  END SUBROUTINE InitPhysics

  !> Sets soundspeeds at cell-centers
  PURE SUBROUTINE SetSoundSpeeds_center(this,Mesh,bccsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)    :: bccsound
    !------------------------------------------------------------------------!
    this%bccsound(:,:,:) = bccsound(:,:,:)
  END SUBROUTINE SetSoundSpeeds_center

  !> Sets soundspeeds at cell-faces
  PURE SUBROUTINE SetSoundSpeeds_faces(this,Mesh,fcsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES), &
                         INTENT(IN)    :: fcsound
    !------------------------------------------------------------------------!
    this%fcsound(:,:,:,:) = fcsound(:,:,:,:)
  END SUBROUTINE SetSoundSpeeds_faces

  !> \todo NOT VERIFIED - most probably wrong: taken from CalcGeometricalSources in euler3Disotherm
  !! at least order is not consistent
  !!
  !! global elemental routine
  ELEMENTAL FUNCTION MomentumSourcesX(my,mz,vx,vy,vz,P,cxyx,cyxy,czxz,cxzx) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: my,mz,vx,vy,vz,P,cxyx,cyxy,czxz,cxzx
    REAL :: st
    !------------------------------------------------------------------------!
    st  = -my * (cxyx * vx - cyxy * vy) + mz * (czxz * vz - cxzx * vx) + (cyxy + czxz) * P
  END FUNCTION MomentumSourcesX

  !> \todo NOT VERIFIED - most probably wrong: taken from CalcGeometricalSources in euler3Disotherm
  !! at least order is not consistent
  !!
  !! global elemental routine
  ELEMENTAL FUNCTION MomentumSourcesY(mz,mx,vx,vy,vz,P,cxyx,cyxy,czyz,cyzy) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mx,mz,vx,vy,vz,P,cxyx,cyxy,czyz,cyzy
    REAL :: st
    !------------------------------------------------------------------------!
    st = mx * (cxyx * vx - cyxy * vy) + mz * (czyz * vz - cyzy * vy) + (cxyx + czyz) * P
  END FUNCTION MomentumSourcesY

  !> \todo NOT VERIFIED - most probably wrong: taken from CalcGeometricalSources in euler3Disotherm
  !! at least order is not consistent
  !!
  !! global elemental routine
  ELEMENTAL FUNCTION MomentumSourcesZ(mx,my,vx,vy,vz,P,cxzx,czxz,czyz,cyzy) RESULT(st)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mx,my,vx,vy,vz,P,cxzx,czxz,czyz,cyzy
    REAL :: st
    !------------------------------------------------------------------------!
    st  =  mx * (cxzx * vx - czxz * vz) + my * (cyzy * vy - czyz * vz) + (cxzx + cyzy) * P
  END FUNCTION MomentumSourcesZ

  !> \todo NOT VERIFIED
  !!
  !! global elemental routine
  ELEMENTAL SUBROUTINE CalcWaveSpeeds(cs,v,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: minwav,maxwav
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    minwav = MIN(0.,v-cs)
    maxwav = MAX(0.,v+cs)
  END SUBROUTINE CalcWaveSpeeds

  !> \todo NOT VERIFIED
  !!
  !! global elemental routine
  ELEMENTAL SUBROUTINE SetWaveSpeeds(cs,v,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: minwav,maxwav
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    minwav = v-cs
    maxwav = v+cs
  END SUBROUTINE SetWaveSpeeds

  !> Destructor
  SUBROUTINE FinalizePhysics(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("ClosePhysics","not initialized")
    ! deallocate pointer variables used in all physics modules
    DEALLOCATE(this%tmp,this%tmp1,this%tmp2,this%tmp3,this%tmp4,this%tmp5, &
               this%bccsound,this%fcsound,this%pvarname,this%cvarname)
  END SUBROUTINE FinalizePhysics

END MODULE physics_base_mod
