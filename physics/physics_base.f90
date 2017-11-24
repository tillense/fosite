!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_generic.f90                                               #
!#                                                                           #
!# Copyright (C) 2007 - 2016                                                 #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!!
!! \brief generic module for the advection problem
!!
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_base_mod
  USE constants_generic_mod
  USE logging_base_mod
  USE mesh_base_mod
!  USE sources_common, ONLY : Sources_TYP
!  USE mesh_generic, ONLY : Mesh_TYP, Initialized, CARTESIAN
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, ABSTRACT, EXTENDS(logging_base) :: physics_base
     !> \name Variables
     CLASS(constants_base), ALLOCATABLE :: consts   !< physical constants
     REAL                :: gamma,&               !< ratio of spec. heats
                            time,&                !< simulation time       
                            mu, &                 !< mean molecular weight
                            csiso, &              !< isothermal sound speed
                            eps                   !< softening length
     INTEGER             :: VNUM, &               !< number of variables
                            DIM, &                !< Dimension (2 or 3)
                            DENSITY,PRESSURE,ENERGY,SGSPRESSURE,SGSENERGY, &
                            XVELOCITY,XMOMENTUM,&
                            YVELOCITY,YMOMENTUM,&
                            ZVELOCITY,ZMOMENTUM   !< array indicies for primitive and conservative variables
     LOGICAL             :: supports_absorbing    !< absorbing boundary conditions supported
                            !! \details .TRUE. if absorbing boundary conditions are supported by the physics module
     LOGICAL             :: supports_farfield     !< farfield boundary conditions supported
                            !! \details .TRUE. if farfield boundary conditions are supported by the physics module
     CHARACTER(LEN=16), DIMENSION(:), POINTER &
                         :: pvarname,cvarname     !< names of variables
     REAL, DIMENSION(:,:,:), POINTER &
                         :: bccsound, &           !< bary centered speed of sound
                            amin, amax, &
                            bmin, bmax, &        
                            cmin, cmax,&          !< wave speeds
                            bcradius, &           !< distance to the origin bary center values
                            divposvec, &          !< divergence of the position vector
                            bphi, &               !< bary centered constant gravitational potential
                            tmp                   !< temporary storage
     REAL, DIMENSION(:,:,:,:), POINTER &
                         :: fcsound, &            !< speed of sound faces
                            fradius, &            !< distance to the origin face values
                            tmin, tmax, &         !< temporary storage
                            bcposvec, &           !< curvilinear components of the position vector bary center values
                            w => NULL(), &        !< fargo bulk velocity
                            fphi, &               !< face centered constant gravitational potential
                            hy                    !< chy or fhy depending on reconstruction
     REAL, DIMENSION(:,:,:,:,:), POINTER &
                         :: fcent, &              !< centrifugal force
                            fposvec               !< curvilinear components of the position vector face values
  CONTAINS
    PROCEDURE :: InitPhysics
    PROCEDURE :: MaxWaveSpeeds
    PROCEDURE :: CalculateFluxesX
    PROCEDURE :: CalculateFluxesY
    PROCEDURE :: CalculateFluxesZ
!   PROCEDURE ::  CalculateCharSystemX
!   PROCEDURE ::  CalculateCharSystemY
!   PROCEDURE ::  CalculateBoundaryDataX
!   PROCEDURE ::  CalculateBoundaryDataY
!   PROCEDURE ::  CalculatePrim2RiemannX
!   PROCEDURE ::  CalculatePrim2RiemannY
!   PROCEDURE ::  CalculateRiemann2PrimX
!   PROCEDURE ::  CalculateRiemann2PrimY
!   PROCEDURE ::  CalculateStresses
    PROCEDURE (CalcSoundSpeeds_center), DEFERRED :: CalcSoundSpeeds_center
    PROCEDURE (CalcSoundSpeeds_faces) , DEFERRED :: CalcSoundSpeeds_faces
    GENERIC   :: CalcSoundSpeeds => CalcSoundSpeeds_center, &
                   CalcSoundSpeeds_faces
    PROCEDURE :: ExternalSources
!   PROCEDURE ::  ViscositySources&
!   PROCEDURE ::  SGSSources
!   PROCEDURE ::  CalculateSGSTensor
!    PROCEDURE (MASKS), DEFERRED :: ReflectionMasks
!    PROCEDURE (MASKS), DEFERRED :: AxisMasks
!    PROCEDURE :: GetSoundSpeed_adiabatic
    PROCEDURE :: FinalizePhysics
    PROCEDURE :: CalcWaveSpeeds_center
    PROCEDURE :: CalcWaveSpeeds_faces
    PROCEDURE (CalcWaveSpeeds), DEFERRED :: CalcWaveSpeeds
    GENERIC   :: CalculateWaveSpeeds => &
                   CalcWaveSpeeds, &
                   CalcWaveSpeeds_center, &
                   CalcWaveSpeeds_faces
    PROCEDURE :: GeometricalSources_center!, GeometricalSources_faces
    PROCEDURE (CalcGeometricalSources), DEFERRED :: CalcGeometricalSources
    GENERIC   :: GeometricalSources => &
                   CalcGeometricalSources, &
                   GeometricalSources_center!, &
                   !GeometricalSources_faces
    PROCEDURE :: Convert2Primitive_center
    PROCEDURE :: Convert2Primitive_centsub
    PROCEDURE :: Convert2Primitive_faces
    PROCEDURE :: Convert2Primitive_facesub
    PROCEDURE (Cons2Prim), DEFERRED :: Cons2Prim
    GENERIC   :: Convert2Primitive => &
                   Cons2Prim, &
                   Convert2Primitive_center, &
                   Convert2Primitive_centsub, &
                   Convert2Primitive_faces, &
                   Convert2Primitive_facesub
    PROCEDURE :: Convert2Conservative_center
    PROCEDURE :: Convert2Conservative_centsub
    PROCEDURE :: Convert2Conservative_faces
    PROCEDURE :: Convert2Conservative_facesub
    PROCEDURE (Prim2Cons), DEFERRED :: Prim2Cons
    GENERIC   :: Convert2Conservative => &
                   Prim2Cons, &
                   Convert2Conservative_center, &
                   Convert2Conservative_centsub, &
                   Convert2Conservative_faces, &
                   Convert2Conservative_facesub
    PROCEDURE :: SetSoundSpeeds_center
    PROCEDURE :: SetSoundSpeeds_faces
    GENERIC   :: SetSoundSpeeds => &
                   SetSoundSpeeds_center, &
                   SetSoundSpeeds_faces
    PROCEDURE (CalcFlux), DEFERRED :: CalcFlux
  END TYPE physics_base

  ABSTRACT INTERFACE
    PURE SUBROUTINE Masks(this,reflX,reflY,reflZ)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN)            :: this
      LOGICAL, DIMENSION(this%VNUM), INTENT(OUT) :: reflX,reflY,reflZ
    END SUBROUTINE
  END INTERFACE
  !ABSTRACT INTERFACE
  !  PURE SUBROUTINE CalcSoundSpeeds_center(this,Mesh,pvar)
  !  IMPORT :: physics_base,mesh_base
  !  IMPLICIT NONE
  !  CLASS(physics_base),INTENT(INOUT) :: this
  !  CLASS(mesh_base), INTENT(IN) :: Mesh
  !  REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), &
  !       INTENT(IN) :: pvar
  !  END SUBROUTINE
  !END INTERFACE
  !ABSTRACT INTERFACE
  !  PURE SUBROUTINE CalcSoundSpeeds_faces(this,Mesh,prim)
  !  IMPORT physics_base, mesh_base
  !  CLASS(physics_base),INTENT(INOUT) :: this
  !  CLASS(mesh_base), INTENT(IN) :: Mesh
  !  REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM), &
  !       INTENT(IN) :: prim
  !  END SUBROUTINE
  !END INTERFACE
  ABSTRACT INTERFACE
    ELEMENTAL SUBROUTINE CalcFlux(this,cs,rho,v,m1,m2,m3,f1,f2,f3,f4)
      IMPORT physics_base
      IMPLICIT NONE
      CLASS(physics_base), INTENT(IN) :: this
      REAL, INTENT(IN)                :: cs,rho,v,m1,m2,m3
      REAL, INTENT(OUT)               :: f1, f2, f3, f4
    END SUBROUTINE
    ELEMENTAL SUBROUTINE CalcWaveSpeeds(this,cs,v,amin,amax)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN) :: this
      REAL, INTENT(IN)                :: cs,v
      REAL, INTENT(OUT)               :: amin,amax
    END SUBROUTINE
    ELEMENTAL SUBROUTINE CalcGeometricalSources(this,mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz,srho,smx,smy,smz)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN) :: this
      REAL, INTENT(IN)                :: mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz
      REAL, INTENT(OUT)               :: srho, smx, smy, smz
    END SUBROUTINE
    ELEMENTAL SUBROUTINE Prim2Cons(this,rho_in,u,v,w,rho_out,mu,mv,mw)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN) :: this
      REAL, INTENT(IN)                :: rho_in,u,v,w
      REAL, INTENT(OUT)               :: rho_out,mu,mv,mw
    END SUBROUTINE
    ELEMENTAL SUBROUTINE Cons2Prim(this,rho_in,mu,mv,mw,rho_out,u,v,w)
      IMPORT physics_base
      CLASS(physics_base), INTENT(IN) :: this
      REAL, INTENT(IN)                :: rho_in,mu,mv,mw
      REAL, INTENT(OUT)               :: rho_out,u,v,w
    END SUBROUTINE
    PURE SUBROUTINE CalcSoundSpeeds_center(this,Mesh,pvar)
      IMPORT physics_base,mesh_base
      IMPLICIT NONE
      CLASS(physics_base),INTENT(INOUT) :: this
      CLASS(mesh_base),INTENT(IN)       :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
        INTENT(IN)                      :: pvar
    END SUBROUTINE
    PURE SUBROUTINE CalcSoundSpeeds_faces(this,Mesh,prim)
      IMPORT physics_base,mesh_base
      IMPLICIT NONE
      CLASS(physics_base),INTENT(INOUT) :: this
      CLASS(mesh_base),INTENT(IN)       :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%VNUM), &
        INTENT(IN)                      :: prim
    END SUBROUTINE
  END INTERFACE
  !INTERFACE GetSoundSpeed_adiabatic
  !   MODULE PROCEDURE GetSoundSpeed_euler2D
  !END INTERFACE
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
  INTEGER, PARAMETER :: EULER3D_ISOTH       = 14
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       physics_base, &
       ! constants
       EULER2D, EULER2D_ISOTHERM, EULER3D_ROTSYM, EULER3D_ROTAMT, & 
       EULER2D_SGS, EULER3D_ROTSYMSGS, EULER3D_ROTAMTSGS, &
       SI, CGS, GEOMETRICAL, EULER2D_ISOIAMT, &
       EULER2D_IAMT, EULER2D_IAMROT, EULER2D_ISOIAMROT, &
       EULER3D_ISOTH
     ! methods
  !--------------------------------------------------------------------------!

CONTAINS

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
    INTEGER                            :: err, i, valwrite
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: problem,vnum
    !------------------------------------------------------------------------!
    CALL this%InitLogging(problem,pname)

    ! check initialization of Mesh
    IF (.NOT.Mesh%Initialized()) &
         CALL this%Error("InitPhysics","mesh module uninitialized")

    ! units
    CALL GetAttr(config, "units", units, SI)
    CALL new_constants(this%consts, units)

    !CALL GetAttr(config, "problem", problem)

    ! ratio of specific heats
    CALL GetAttr(config, "gamma", this%gamma, 1.4)

    ! mean molecular weight
    CALL GetAttr(config, "mu", this%mu, 0.029)

    ! isothermal sound speed
    CALL GetAttr(config, "cs", this%csiso, 0.0)

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
    ALLOCATE(this%pvarname(this%vnum),this%cvarname(this%vnum),                                        &
             this%amin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
             this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
             this%bmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
             this%bmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
             this%cmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
             this%cmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
              this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
             this%tmin(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4),           &
             this%tmax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4),           &
         this%bccsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),             &
          this%fcsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces), &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitPhysics", "Unable to allocate memory.")
    
    this%amax(:,:,:)   = 0.
    this%bmin(:,:,:)   = 0.
    this%bmax(:,:,:)   = 0.
    this%tmp(:,:,:)    = 0.
    this%tmin(:,:,:,:) = 0.
    this%tmax(:,:,:,:) = 0.

    IF(this%csiso.GT.0.) THEN
      this%bccsound(:,:,:)  = this%csiso
      this%fcsound(:,:,:,:) = this%csiso
    ELSE
      this%bccsound(:,:,:)  = 0.
      this%fcsound(:,:,:,:) = 0.
    END IF

!    SELECT CASE(problem)
!    CASE(EULER2D)
!       CALL InitPhysics_euler2D(this,Mesh,problem)
!    CASE(EULER2D_ISOTHERM)
!       CALL InitPhysics_euler2Dit(this,Mesh,problem)
!    CASE(EULER3D_ROTSYM)
!       CALL InitPhysics_euler3Drs(this,Mesh,problem)
!    CASE(EULER3D_ROTAMT)
!       CALL InitPhysics_euler3Dra(this,Mesh,problem)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL InitPhysics_euler3DrsSGS(this,Mesh,problem)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL InitPhysics_euler3DraSGS(this,Mesh,problem)
!    CASE(EULER2D_SGS)
!       CALL InitPhysics_euler2Dsgs(this,Mesh,problem)
!    CASE(EULER2D_IAMT)
!       CALL InitPhysics_euler2Dia(this,Mesh,problem)
!    CASE(EULER2D_IAMROT)
!       CALL InitPhysics_euler2Diar(this,Mesh,problem)
!    CASE(EULER2D_ISOIAMT)
!       CALL InitPhysics_euler2Ditia(this,Mesh,problem)
!    CASE(EULER2D_ISOIAMROT)
!       CALL InitPhysics_euler2Ditiar(this,Mesh,problem)
!    CASE DEFAULT
!       CALL Error(this, "InitPhysics", "Unknown advection problem.")
!    END SELECT

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

    CALL GetAttr(config, "output/wave_speeds", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      CALL SetAttr(IO, "amin", &
                   this%amin(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CALL SetAttr(IO, "amax", &
                   this%amax(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CALL SetAttr(IO, "bmin", &
                   this%bmin(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CALL SetAttr(IO, "bmax", &
                   this%bmax(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CALL SetAttr(IO, "cmin", &
                   this%cmin(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      CALL SetAttr(IO, "cmax", &
                   this%cmax(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
      END IF

    this%time = -1.
    ! print some information
    CALL this%Info(" PHYSICS--> advection problem: " // TRIM(this%GetName()))

  END SUBROUTINE InitPhysics


  PURE SUBROUTINE SetSoundSpeeds_center(this,Mesh,bccsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this 
    CLASS(mesh_base),INTENT(IN)        :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                                       :: bccsound
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: bccsound
    !------------------------------------------------------------------------!
    this%bccsound(:,:,:) = bccsound(:,:,:)
  END SUBROUTINE SetSoundSpeeds_center


  PURE SUBROUTINE SetSoundSpeeds_faces(this,Mesh,fcsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this 
    CLASS(mesh_base),INTENT(IN)        :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces) &
                                       :: fcsound
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: fcsound
    !------------------------------------------------------------------------!
    this%fcsound(:,:,:,:) = fcsound(:,:,:,:)
  END SUBROUTINE SetSoundSpeeds_faces
  

!  PURE SUBROUTINE CalcSoundSpeeds_center(this,Mesh,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_base),INTENT(INOUT) :: this
!    CLASS(mesh_base),INTENT(IN) :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
!         :: pvar
!    !------------------------------------------------------------------------!
!    INTEGER           :: i,j
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: pvar
!    !------------------------------------------------------------------------!
!!CDIR COLLAPSE
!    DO j=Mesh%JGMIN,Mesh%JGMAX
!      DO i=Mesh%IGMIN,Mesh%IGMAX
!        this%bccsound(i,j) = this%GetSoundSpeed(&
!          this%gamma,&
!          pvar(i,j,this%DENSITY),&
!          pvar(i,j,this%PRESSURE))
!      END DO
!    END DO
!  END SUBROUTINE CalcSoundSpeeds_center
!
!
!  PURE SUBROUTINE CalcSoundSpeeds_faces(this,Mesh,prim)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(physics_base),INTENT(INOUT) :: this
!    CLASS(mesh_base),INTENT(IN) :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
!         :: prim
!    !------------------------------------------------------------------------!
!    INTEGER           :: i,j,k
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: prim
!    !------------------------------------------------------------------------!
!!CDIR COLLAPSE
!    DO k=1,4
!      DO j=Mesh%JGMIN,Mesh%JGMAX
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          this%fcsound(i,j,k) = this%GetSoundSpeed(&
!            this%gamma,&
!            prim(i,j,k,this%DENSITY),&
!            prim(i,j,k,this%PRESSURE))
!        END DO
!      END DO
!    END DO
!  END SUBROUTINE CalcSoundSpeeds_faces



  PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)       :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER                           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)                        :: pvar
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
         DO i=Mesh%IGMIN,Mesh%IGMAX
          ! x-direction
!CDIR IEXPAND
            CALL this%CalcWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%XVELOCITY),&
                 this%amin(i,j,k),this%amax(i,j,k))
          ! y-direction
!CDIR IEXPAND
          CALL this%CalcWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%YVELOCITY),&
               this%bmin(i,j,k),this%bmax(i,j,k))
          ! z-direction
!CDIR IEXPAND
          CALL this%CalcWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%ZVELOCITY),&
               this%cmin(i,j,k),this%cmax(i,j,k))
         END DO
      END DO
    END DO  
  END SUBROUTINE CalcWaveSpeeds_center


  PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this 
    CLASS(mesh_base),INTENT(IN)        :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%VNUM) &
                                       :: prim
    !------------------------------------------------------------------------!
    INTEGER                            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: prim
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell interfaces
!CDIR COLLAPSE
   DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
         DO i=Mesh%IGMIN,Mesh%IGMAX
            ! western
!CDIR IEXPAND
            CALL this%CalcWaveSpeeds(this%fcsound(i,j,k,1),         &
                                     prim(i,j,k,1,this%XVELOCITY),  &
                                     this%tmin(i,j,k,1),            &
                                     this%tmax(i,j,k,1)             &
                                    )
          ! eastern
!CDIR IEXPAND
          CALL this%CalcWaveSpeeds(this%fcsound(i,j,k,2),        &
                                   prim(i,j,k,2,this%XVELOCITY), &
                                   this%amin(i,j,k),             &
                                   this%amax(i,j,k)              &
                                  )
          ! southern
!CDIR IEXPAND
          CALL this%CalcWaveSpeeds(this%fcsound(i,j,k,3),        &
                                   prim(i,j,k,3,this%YVELOCITY), &
                                   this%tmin(i,j,k,2),           &
                                   this%tmax(i,j,k,2)            &
                                  )
          ! northern
!CDIR IEXPAND
          CALL this%CalcWaveSpeeds(this%fcsound(i,j,k,4),        &
                                   prim(i,j,k,4,this%YVELOCITY), &
                                   this%bmin(i,j,k),             &
                                   this%bmax(i,j,k)              &
                                  )
          ! bottom
!CDIR IEXPAND
          CALL this%CalcWaveSpeeds(this%fcsound(i,j,k,5),        &
                                   prim(i,j,k,5,this%ZVELOCITY), &
                                   this%tmin(i,j,k,3),           &
                                   this%tmax(i,j,k,3)            &
                                  )
          ! top
!CDIR IEXPAND
          CALL this%CalcWavespeeds(this%fcsound(i,j,k,6),        &
                                   prim(i,j,k,6,this%ZVELOCITY), &
                                   this%cmin(i,j,k),             &
                                   this%cmax(i,j,k)              &
                                  )
          END DO
       END DO
    END DO
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
    DO k=Mesh%KGMIN,Mesh%KGMAX
       DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
          DO i=Mesh%IMIN-1,Mesh%IMAX
             ! western interfaces
             this%amin(i,j,k) = MIN(this%tmin(i+1,j,k,1),this%amin(i,j,k))
             ! eastern interfaces
             this%amax(i,j,k) = MAX(this%tmax(i+1,j,k,1),this%amax(i,j,k))
          END DO
       END DO
    END DO
!CDIR COLLAPSE
   DO k=Mesh%KGMIN,Mesh%KGMAX
       DO j=Mesh%JMIN-1,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! southern interfaces
             this%bmin(i,j,k) = MIN(this%tmin(i,j+1,k,2),this%bmin(i,j,k))
             ! northern interfaces
             this%bmax(i,j,k) = MAX(this%tmax(i,j+1,k,2),this%bmax(i,j,k))
         END DO
      END DO
  END DO
!CDIR COLLAPSE
    DO k=Mesh%KMIN-1,Mesh%KMAX
       DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
             !bottom interfaces
             this%cmin(i,j,k) = MIN(this%tmin(i,j,k+1,3),this%cmin(i,j,k))
             !top interfaces
             this%cmax(i,j,k) = MAX(this%tmax(i,j,k+1,3),this%cmax(i,j,k))
          END DO
       END DO
    END DO
  END SUBROUTINE CalcWaveSpeeds_faces


  PURE SUBROUTINE MaxWaveSpeeds(this,Mesh,time,pvar,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this 
    CLASS(mesh_base), INTENT(IN)       :: Mesh
    REAL                               :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum) &
                                       :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                                       :: amax
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: pvar,time
    INTENT(OUT)                        :: amax
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%CalculateWaveSpeeds(Mesh,pvar)

    amax(:,:,:,1) = MAX(this%amax,-this%amin)
    amax(:,:,:,2) = MAX(this%bmax,-this%bmin)
    amax(:,:,:,3) = MAX(this%cmax,-this%cmin)
  END SUBROUTINE MaxWaveSpeeds


  PURE SUBROUTINE CalculateFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this 
    CLASS(mesh_base), INTENT(IN)    :: Mesh
    INTEGER                         :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: prim,cons,xfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: nmin,nmax,prim,cons
    INTENT(OUT)                     :: xfluxes
    !------------------------------------------------------------------------!
    CALL this%CalcFlux(this%fcsound(:,:,:,nmin:nmax),        &
                       prim(:,:,:,nmin:nmax,this%DENSITY),   &
                       prim(:,:,:,nmin:nmax,this%XVELOCITY), &
                       cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
                       cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
                       cons(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                    xfluxes(:,:,:,nmin:nmax,this%DENSITY),   &
                    xfluxes(:,:,:,nmin:nmax,this%XMOMENTUM), &
                    xfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
                    xfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM)  &
                      )
  END SUBROUTINE CalculateFluxesX


  PURE SUBROUTINE CalculateFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this 
    CLASS(mesh_base), INTENT(IN)    :: Mesh
    INTEGER                         :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: prim,cons,yfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: nmin,nmax,prim,cons
    INTENT(OUT)                     :: yfluxes
    !------------------------------------------------------------------------!
    CALL this%CalcFlux(this%fcsound(:,:,:,nmin:nmax),        &
                       prim(:,:,:,nmin:nmax,this%DENSITY),   &
                       prim(:,:,:,nmin:nmax,this%YVELOCITY), &
                       cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
                       cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
                       cons(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                    yfluxes(:,:,:,nmin:nmax,this%DENSITY),   &
                    yfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
                    yfluxes(:,:,:,nmin:nmax,this%XMOMENTUM), &
                    yfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM)  &
                      )
  END SUBROUTINE CalculateFluxesY

  PURE SUBROUTINE CalculateFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this 
    CLASS(mesh_base), INTENT(IN)    :: Mesh
    INTEGER                         :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: prim,cons,zfluxes
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: nmin,nmax,prim,cons
    INTENT(OUT)                     :: zfluxes
    !------------------------------------------------------------------------!
    CALL this%CalcFlux(this%fcsound(:,:,:,nmin:nmax),                &
                               prim(:,:,:,nmin:nmax,this%DENSITY),   &
                               prim(:,:,:,nmin:nmax,this%ZVELOCITY), &
                               cons(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                               cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
                               cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
                            zfluxes(:,:,:,nmin:nmax,this%DENSITY),   &
                            zfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                            zfluxes(:,:,:,nmin:nmax,this%XMOMENTUM), &
                            zfluxes(:,:,:,nmin:nmax,this%YMOMENTUM)  &
                      )
  END SUBROUTINE CalculateFluxesZ



!  PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i,dir,pvar,lambda,xvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i,dir
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,xvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i,dir,pvar
!    INTENT(INOUT)     :: lambda
!    INTENT(OUT)       :: xvar
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcCharSystemX_euler2D(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcCharSystemX_euler2Dit(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcCharSystemX_euler3Drs(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcCharSystemX_euler3Dra(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcCharSystemX_euler3DrsSGS(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcCharSystemX_euler3DraSGS(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER2D_SGS)
!       CALL CalcCharSystemX_euler2Dsgs(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER2D_IAMT)
!       CALL CalcCharSystemX_euler2Dia(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcCharSystemX_euler2Diar(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER2D_ISOIAMT)
!       CALL CalcCharSystemX_euler2Ditia(this,Mesh,i,dir,pvar,lambda,xvar)
!    CASE(EULER2D_ISOIAMROT)
!!!$ FIXME, not implemented
!!       CALL CalcCharSystemX_euler2Ditiar(this,Mesh,i,dir,pvar,lambda,xvar)
!    END SELECT
!  END SUBROUTINE CalculateCharSystemX
!
!
!  PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j,dir,pvar,lambda,xvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j,dir
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,xvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j,dir,pvar
!    INTENT(INOUT)     :: lambda
!    INTENT(OUT)       :: xvar
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcCharSystemY_euler2D(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcCharSystemY_euler2Dit(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcCharSystemY_euler3Drs(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcCharSystemY_euler3Dra(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcCharSystemY_euler3DrsSGS(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcCharSystemY_euler3DraSGS(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER2D_SGS)
!       CALL CalcCharSystemY_euler2Dsgs(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER2D_IAMT)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcCharSystemY_euler2Diar(this,Mesh,j,dir,pvar,lambda,xvar)
!    CASE(EULER2D_ISOIAMT)
!    CASE(EULER2D_ISOIAMROT)
!    END SELECT
!  END SUBROUTINE CalculateCharSystemY
!
!
!  PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i,dir,xvar,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i,dir
!    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i,dir,xvar
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcBoundaryDataX_euler2D(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcBoundaryDataX_euler2Dit(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcBoundaryDataX_euler3Drs(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcBoundaryDataX_euler3Dra(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcBoundaryDataX_euler3DrsSGS(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcBoundaryDataX_euler3DraSGS(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER2D_SGS)
!       CALL CalcBoundaryDataX_euler2Dsgs(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER2D_IAMT)
!       CALL CalcBoundaryDataX_euler2Dia(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcBoundaryDataX_euler2Diar(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER2D_ISOIAMT)
!       CALL CalcBoundaryDataX_euler2Ditia(this,Mesh,i,dir,xvar,pvar)
!    CASE(EULER2D_ISOIAMROT)
!!!$ FIXME, not implemented
!!       CALL CalcBoundaryDataX_euler2Ditiar(this,Mesh,i,dir,xvar,pvar)
!    END SELECT
!  END SUBROUTINE CalculateBoundaryDataX
!
!
!  PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j,dir,xvar,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j,dir
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: xvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j,dir,xvar
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcBoundaryDataY_euler2D(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcBoundaryDataY_euler2Dit(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcBoundaryDataY_euler3Drs(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcBoundaryDataY_euler3Dra(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcBoundaryDataY_euler3Drssgs(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcBoundaryDataY_euler3Drasgs(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER2D_SGS)
!       CALL CalcBoundaryDataY_euler2Dsgs(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER2D_IAMT)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcBoundaryDataY_euler2Diar(this,Mesh,j,dir,xvar,pvar)
!    CASE(EULER2D_ISOIAMT)
!    CASE(EULER2D_ISOIAMROT)
!    END SELECT
!  END SUBROUTINE CalculateBoundaryDataY
!
!  PURE SUBROUTINE CalculatePrim2RiemannX(this,Mesh,i,pvar,lambda,Rinv)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,Rinv
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i,pvar
!    INTENT(INOUT)     :: lambda
!    INTENT(OUT)       :: Rinv
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcPrim2RiemannX_euler2D(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcPrim2RiemannX_euler2Dit(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcPrim2RiemannX_euler2Diar(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER2D_SGS)
!       CALL CalcPrim2RiemannX_euler2Dsgs(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcPrim2RiemannX_euler3Drs(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcPrim2RiemannX_euler3Drssgs(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcPrim2RiemannX_euler3Dra(this,Mesh,i,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcPrim2RiemannX_euler3Drasgs(this,Mesh,i,pvar,lambda,Rinv)
!    END SELECT
!  END SUBROUTINE CalculatePrim2RiemannX
!
!  PURE SUBROUTINE CalculatePrim2RiemannY(this,Mesh,j,pvar,lambda,Rinv)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,Rinv
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j,pvar
!    INTENT(INOUT)     :: lambda
!    INTENT(OUT)       :: Rinv
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcPrim2RiemannY_euler2D(this,Mesh,j,pvar,lambda,Rinv)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcPrim2RiemannY_euler2Dit(this,Mesh,j,pvar,lambda,Rinv)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcPrim2RiemannY_euler2Diar(this,Mesh,j,pvar,lambda,Rinv)
!    CASE(EULER2D_SGS)
!       CALL CalcPrim2RiemannY_euler2Dsgs(this,Mesh,j,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcPrim2RiemannY_euler3Drs(this,Mesh,j,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcPrim2RiemannY_euler3Drssgs(this,Mesh,j,pvar,lambda,Rinv)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcPrim2RiemannY_euler3Dra(this,Mesh,j,pvar,lambda,Rinv) 
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcPrim2RiemannY_euler3Drasgs(this,Mesh,j,pvar,lambda,Rinv) 
!    END SELECT
!  END SUBROUTINE CalculatePrim2RiemannY
!
! PURE SUBROUTINE CalculateRiemann2PrimX(this,Mesh,i,Rinv,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX,this%VNUM) :: Rinv
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i,Rinv
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcRiemann2PrimX_euler2D(this,Mesh,i,Rinv,pvar)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcRiemann2PrimX_euler2Dit(this,Mesh,i,Rinv,pvar)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcRiemann2PrimX_euler2Diar(this,Mesh,i,Rinv,pvar)
!    CASE(EULER2D_SGS)
!       CALL CalcRiemann2PrimX_euler2Dsgs(this,Mesh,i,Rinv,pvar)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcRiemann2PrimX_euler3Drs(this,Mesh,i,Rinv,pvar)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcRiemann2PrimX_euler3Drssgs(this,Mesh,i,Rinv,pvar)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcRiemann2PrimX_euler3Dra(this,Mesh,i,Rinv,pvar)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcRiemann2PrimX_euler3Drasgs(this,Mesh,i,Rinv,pvar)
!    END SELECT
!  END SUBROUTINE CalculateRiemann2PrimX
!
!  PURE SUBROUTINE CalculateRiemann2PrimY(this,Mesh,j,Rinv,pvar)
!    IMPLICIT NONE
!   !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,this%VNUM) :: Rinv
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j,Rinv
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcRiemann2PrimY_euler2D(this,Mesh,j,Rinv,pvar)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcRiemann2PrimY_euler2Dit(this,Mesh,j,Rinv,pvar)
!    CASE(EULER2D_IAMROT)
!!!$ FIXME, not implemented
!!!$       CALL CalcRiemann2PrimY_euler2Diar(this,Mesh,j,Rinv,pvar)
!    CASE(EULER2D_SGS)
!       CALL CalcRiemann2PrimY_euler2Dsgs(this,Mesh,j,Rinv,pvar)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcRiemann2PrimY_euler3Drs(this,Mesh,j,Rinv,pvar)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcRiemann2PrimY_euler3Drssgs(this,Mesh,j,Rinv,pvar)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcRiemann2PrimY_euler3Dra(this,Mesh,j,Rinv,pvar)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcRiemann2PrimY_euler3Drasgs(this,Mesh,j,Rinv,pvar)
!    END SELECT
!  END SUBROUTINE CalculateRiemann2PrimY
!
!  PURE SUBROUTINE CalculateStresses(this,Mesh,pvar,dynvis,bulkvis, &
!       btxx,btxy,btxz,btyy,btyz,btzz)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
!         dynvis,bulkvis,btxx,btxy,btxz,btyy,btyz,btzz
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,pvar,dynvis,bulkvis
!    INTENT(INOUT)     :: this
!    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
!    !------------------------------------------------------------------------!
!!CDIR IEXPAND
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!       CALL CalcStresses_euler2D(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    CASE(EULER2D_ISOTHERM)
!       CALL CalcStresses_euler2Dit(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    CASE(EULER3D_ROTSYM)
!       CALL CalcStresses_euler3Drs(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btxz,btyy,btyz,btzz)
!    CASE(EULER3D_ROTAMT)
!       CALL CalcStresses_euler3Dra(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btxz,btyy,btyz,btzz)
!    CASE(EULER3D_ROTSYMSGS)
!       CALL CalcStresses_euler3Drssgs(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btxz,btyy,btyz,btzz)
!    CASE(EULER3D_ROTAMTSGS)
!       CALL CalcStresses_euler3Drasgs(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btxz,btyy,btyz,btzz)
!    CASE(EULER2D_SGS)
!       CALL CalcStresses_euler2Dsgs(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    CASE(EULER2D_IAMT)
!       CALL CalcStresses_euler2Dia(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    CASE(EULER2D_IAMROT)
!       CALL CalcStresses_euler2Diar(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    CASE(EULER2D_ISOIAMT)
!       CALL CalcStresses_euler2Ditia(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    CASE(EULER2D_ISOIAMROT)
!       CALL CalcStresses_euler2Ditiar(this,Mesh,pvar,dynvis,bulkvis, &
!            btxx,btxy,btyy)
!    END SELECT
!  END SUBROUTINE CalculateStresses


  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)        :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum) &
                                       :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER                            :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: pvar,cvar
    INTENT(OUT)                        :: sterm
    !------------------------------------------------------------------------!
    ! compute geometrical source only for non-cartesian mesh except for the
    ! EULER2D_IAMROT case for which geometrical sources are always necessary.
    IF ((Mesh%Geometry%GetType().NE.CARTESIAN).OR. &
        (this%GetType().EQ.EULER2D_IAMROT).OR. &
        (this%GetType().EQ.EULER2D_ISOIAMROT)) THEN
 
 
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
         DO i=Mesh%IGMIN,Mesh%IGMAX
            CALL this%CalcGeometricalSources(cvar(i,j,k,this%XMOMENTUM),                       &
                                             cvar(i,j,k,this%YMOMENTUM),                       &
                                             cvar(i,j,k,this%ZMOMENTUM),                       &
                                             pvar(i,j,k,this%XVELOCITY),                       &
                                             pvar(i,j,k,this%YVELOCITY),                       &
                                             pvar(i,j,k,this%ZVELOCITY),                       &
                                             pvar(i,j,k,this%DENSITY)*this%bccsound(i,j,k)**2, &
                                Mesh%cxyx%bcenter(i,j,k),                                      &
                                Mesh%cxzx%bcenter(i,j,k),                                      &
                                Mesh%cyxy%bcenter(i,j,k),                                      &
                                Mesh%cyzy%bcenter(i,j,k),                                      &
                                Mesh%czxz%bcenter(i,j,k),                                      &
                                Mesh%czyz%bcenter(i,j,k),                                      &
                                            sterm(i,j,k,this%DENSITY),                         &
                                            sterm(i,j,k,this%XMOMENTUM),                       &
                                            sterm(i,j,k,this%YMOMENTUM),                       &
                                            sterm(i,j,k,this%ZMOMENTUM)                        &
                                           )
         END DO
      END DO
   END DO
    ! reset ghost cell data
    sterm(Mesh%IGMIN:Mesh%IMIN-1,:,:,:) = 0.0
    sterm(Mesh%IMAX+1:Mesh%IGMAX,:,:,:) = 0.0
    sterm(:,Mesh%JGMIN:Mesh%JMIN-1,:,:) = 0.0
    sterm(:,Mesh%JMAX+1:Mesh%JGMAX,:,:) = 0.0
    sterm(:,:,Mesh%KGMIN:Mesh%KMIN-1,:) = 0.0
    sterm(:,:,Mesh%KMAX+1:Mesh%KGMAX,:) = 0.0
    END IF
  END SUBROUTINE GeometricalSources_center

  
!  PURE SUBROUTINE GeometricalSources_faces(this,Mesh,prim,cons,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%vnum) &
!         :: prim,cons
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%vnum) &
!         :: sterm
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,prim,cons
!    INTENT(INOUT)     :: this
!    INTENT(OUT)       :: sterm
!    !------------------------------------------------------------------------!
!    ! compute geometrical source only for non-cartesian mesh except for the
!    ! EULER2D_IAMROT case for which geometrical sources are always necessary.
!    IF ((GetType(Mesh%geometry).NE.CARTESIAN).OR. &
!        (GetType(this).EQ.EULER2D_IAMROT)) THEN
!    ! calculate geometrical sources depending on the advection problem
!!CDIR COLLAPSE
!    DO j=Mesh%JGMIN,Mesh%JGMAX
!       DO i=Mesh%IGMIN,Mesh%IGMAX
!          ! no geometrical density sources
!          sterm(i,j,this%DENSITY)   = 0.
!          ! momentum sources (sum up corner values, don't use SUM function,
!          ! because it prevents COLLAPSING and causes poor vectorization
!          sterm(i,j,this%XMOMENTUM) = MomentumSourcesX_euler2Dit(&
!              cons(i,j,1,this%YMOMENTUM),prim(i,j,1,this%XVELOCITY),&
!              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%DENSITY)*this%fcsound(i,j,1)**2, &
!              Mesh%cxyx%corners(i,j,1),Mesh%cyxy%corners(i,j,1),Mesh%czxz%corners(i,j,1)) &
!            + MomentumSourcesX_euler2Dit(&
!              cons(i,j,2,this%YMOMENTUM),prim(i,j,2,this%XVELOCITY),&
!              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%DENSITY)*this%fcsound(i,j,2)**2, &
!              Mesh%cxyx%corners(i,j,2),Mesh%cyxy%corners(i,j,2),Mesh%czxz%corners(i,j,2)) &
!            + MomentumSourcesX_euler2Dit(&
!              cons(i,j,3,this%YMOMENTUM),prim(i,j,3,this%XVELOCITY),&
!              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%DENSITY)*this%fcsound(i,j,3)**2, &
!              Mesh%cxyx%corners(i,j,3),Mesh%cyxy%corners(i,j,3),Mesh%czxz%corners(i,j,3)) &
!            + MomentumSourcesX_euler2Dit(&
!              cons(i,j,4,this%YMOMENTUM),prim(i,j,4,this%XVELOCITY),&
!              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%DENSITY)*this%fcsound(i,j,4)**2, &
!              Mesh%cxyx%corners(i,j,4),Mesh%cyxy%corners(i,j,4),Mesh%czxz%corners(i,j,4))
!
!          sterm(i,j,this%YMOMENTUM) = MomentumSourcesY_euler2Dit(&
!              cons(i,j,1,this%XMOMENTUM),prim(i,j,1,this%XVELOCITY),&
!              prim(i,j,1,this%YVELOCITY),prim(i,j,1,this%DENSITY)*this%fcsound(i,j,1)**2, &
!              Mesh%cxyx%corners(i,j,1),Mesh%cyxy%corners(i,j,1),Mesh%czyz%corners(i,j,1)) &
!            + MomentumSourcesY_euler2Dit(&
!              cons(i,j,2,this%XMOMENTUM),prim(i,j,2,this%XVELOCITY),&
!              prim(i,j,2,this%YVELOCITY),prim(i,j,2,this%DENSITY)*this%fcsound(i,j,2)**2, &
!              Mesh%cxyx%corners(i,j,2),Mesh%cyxy%corners(i,j,2),Mesh%czyz%corners(i,j,2)) &
!            + MomentumSourcesY_euler2Dit(&
!              cons(i,j,3,this%XMOMENTUM),prim(i,j,3,this%XVELOCITY),&
!              prim(i,j,3,this%YVELOCITY),prim(i,j,3,this%DENSITY)*this%fcsound(i,j,3)**2, &
!              Mesh%cxyx%corners(i,j,3),Mesh%cyxy%corners(i,j,3),Mesh%czyz%corners(i,j,3)) &
!            + MomentumSourcesY_euler2Dit(&
!              cons(i,j,4,this%XMOMENTUM),prim(i,j,4,this%XVELOCITY),&
!              prim(i,j,4,this%YVELOCITY),prim(i,j,4,this%DENSITY)*this%fcsound(i,j,4)**2, &
!              Mesh%cxyx%corners(i,j,4),Mesh%cyxy%corners(i,j,4),Mesh%czyz%corners(i,j,4))
!       END DO
!    END DO
!    ! reset ghost cell data
!    sterm(Mesh%IGMIN:Mesh%IMIN-1,:,:) = 0.0
!    sterm(Mesh%IMAX+1:Mesh%IGMAX,:,:) = 0.0
!    sterm(:,Mesh%JGMIN:Mesh%JMIN-1,:) = 0.0
!    sterm(:,Mesh%JMAX+1:Mesh%JGMAX,:) = 0.0
!    END IF
!  END SUBROUTINE GeometricalSources_faces


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                                    :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                    :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER                         :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: accel,pvar,cvar
    INTENT(OUT)                     :: sterm
    !------------------------------------------------------------------------!
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             sterm(i,j,k,this%DENSITY)   = 0.
             sterm(i,j,k,this%XMOMENTUM) = pvar(i,j,k,this%DENSITY) * accel(i,j,k,1)
             sterm(i,j,k,this%YMOMENTUM) = pvar(i,j,k,this%DENSITY) * accel(i,j,k,2)
             sterm(i,j,k,this%ZMOMENTUM) = pvar(i,j,k,this%DENSITY) * accel(i,j,k,3)
          END DO
       END DO
    END DO
  END SUBROUTINE ExternalSources


!  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
!         pvar,sterm
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
!         btxx,btxy,btxz,btyy,btyz,btzz
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz
!    INTENT(INOUT)     :: this
!    INTENT(OUT)       :: sterm
!    !------------------------------------------------------------------------!
!!CDIR IEXPAND
!    SELECT CASE(GetType(this))
!    CASE(EULER2D)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2D(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    CASE(EULER2D_ISOTHERM)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    CASE(EULER3D_ROTSYM)
!!CDIR IEXPAND
!       CALL ViscositySources_euler3Drs(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
!            btyz,btzz,sterm)
!    CASE(EULER3D_ROTAMT)
!!CDIR IEXPAND
!       CALL ViscositySources_euler3Dra(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
!            btyz,btzz,sterm)
!    CASE(EULER3D_ROTSYMSGS)
!!CDIR IEXPAND
!       CALL ViscositySources_euler3DrsSGS(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
!            btyz,btzz,sterm)
!    CASE(EULER3D_ROTAMTSGS)
!!CDIR IEXPAND
!       CALL ViscositySources_euler3DraSGS(this,Mesh,pvar,btxx,btxy,btxz,btyy, &
!            btyz,btzz,sterm)
!    CASE(EULER2D_SGS)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2Dsgs(this,Mesh,pvar,btxx,btxy,btyy, &
!            sterm)
!    CASE(EULER2D_IAMT)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2Dia(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    CASE(EULER2D_IAMROT)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2Diar(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    CASE(EULER2D_ISOIAMT)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2Ditia(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    CASE(EULER2D_ISOIAMROT)
!!CDIR IEXPAND
!       CALL ViscositySources_euler2Ditiar(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    END SELECT
!  END SUBROUTINE ViscositySources
!
! PURE SUBROUTINE SGSSources(this,Mesh,Sources,pvar,cvar,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    TYPE(Sources_TYP) :: Sources
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
!         pvar,cvar,sterm
!!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
!!         btxx,btxy,btxz,btyy,btyz,btzz
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,pvar,cvar
!    INTENT(INOUT)     :: this,Sources
!    INTENT(OUT)       :: sterm
!    !------------------------------------------------------------------------!
!!CDIR IEXPAND
!    SELECT CASE(GetType(this))
!    CASE(EULER3D_ROTSYMSGS)
!!CDIR IEXPAND
!       CALL SGSSources_euler3DrsSGS(this,Mesh,Sources,pvar,cvar,sterm)
!    CASE(EULER3D_ROTAMTSGS)
!!CDIR IEXPAND
!       CALL SGSSources_euler3DraSGS(this,Mesh,Sources,pvar,cvar,sterm)
!    CASE(EULER2D_SGS)
!!CDIR IEXPAND
!       CALL SGSSources_euler2Dsgs(this,Mesh,Sources,pvar,cvar,sterm)
!    END SELECT
!  END SUBROUTINE SGSSources
!
!  PURE SUBROUTINE CalculateSGSTensor(this,Mesh,Sources,C,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    TYPE(Sources_TYP) :: Sources
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
!         C,pvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,C,pvar
!    INTENT(INOUT)     :: this,Sources
!    !------------------------------------------------------------------------!
!!CDIR IEXPAND
!    SELECT CASE(GetType(this))
!    CASE(EULER3D_ROTSYMSGS)
!!CDIR IEXPAND
!       CALL CalcSGSTensor_euler3DrsSGS(this,Mesh,Sources,C,pvar)
!    CASE(EULER3D_ROTAMTSGS)
!!CDIR IEXPAND
!       CALL CalcSGSTensor_euler3DraSGS(this,Mesh,Sources,C,pvar)
!    CASE(EULER2D_SGS)
!!CDIR IEXPAND
!       CALL CalcSGSTensor_euler2Dsgs(this,Mesh,Sources,C,pvar)
!    END SELECT
!  END SUBROUTINE CalculateSGSTensor

  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum) &
                                    :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: cvar
    INTENT(OUT)                     :: pvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Primitive_centsub(Mesh,Mesh%IGMIN,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,cvar,pvar)
  END SUBROUTINE Convert2Primitive_center

  
  PURE SUBROUTINE Convert2Primitive_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    INTEGER                         :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum) &
                                    :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: i1,i2,j1,j2,k1,k2,cvar
    INTENT(OUT)                     :: pvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Primitive(cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                                cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                                cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                                cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM), &
                                pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                                pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                                pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                                pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY)  &
                               )
  END SUBROUTINE Convert2Primitive_centsub


  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: cons
    INTENT(OUT)                     :: prim
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Primitive_facesub(Mesh,Mesh%IGMIN,Mesh%IGMAX, &
                                             Mesh%JGMIN,Mesh%JGMAX, &
                                             Mesh%KGMIN,Mesh%KGMAX, &
                                        cons,prim)
  END SUBROUTINE Convert2Primitive_faces


  PURE SUBROUTINE Convert2Primitive_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    INTEGER                         :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: i1,i2,j1,j2,k1,k2,cons
    INTENT(OUT)                     :: prim
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Primitive(cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                                cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                                cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                                cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM), &
                                prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                                prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                                prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                                prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY)  &
                               )
  END SUBROUTINE Convert2Primitive_facesub


  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum) &
                                    :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: pvar
    INTENT(OUT)                     :: cvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Conservative_centsub(Mesh,Mesh%IGMIN,Mesh%IGMAX, &
                                                Mesh%JGMIN,Mesh%JGMAX, &
                                                Mesh%KGMIN,Mesh%KGMAX, &
                                           pvar,cvar)
  END SUBROUTINE Convert2Conservative_center


  PURE SUBROUTINE Convert2Conservative_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    INTEGER                         :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%vnum) &
                                    :: cvar,pvar
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: i1,i2,j1,j2,k1,k2,pvar
    INTENT(OUT)                     :: cvar
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Conservative(pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                                   pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                                   pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                                   pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY), &
                                   cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , & 
                                   cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                                   cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                                   cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM)  &
                                  )
  END SUBROUTINE Convert2Conservative_centsub


  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: prim
    INTENT(OUT)                     :: cons
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Conservative_facesub(Mesh,Mesh%IGMIN,Mesh%IGMAX, &
                                                Mesh%JGMIN,Mesh%JGMAX, &
                                                Mesh%KGMIN,Mesh%KGMAX, &
                                           prim,cons)
  END SUBROUTINE Convert2Conservative_faces
  

  PURE SUBROUTINE Convert2Conservative_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(IN) :: this
    CLASS(mesh_base),INTENT(IN)     :: Mesh
    INTEGER                         :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%vnum) &
                                    :: cons,prim
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: i1,i2,j1,j2,k1,k2,prim
    INTENT(OUT)                     :: cons
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Conservative(prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                                   prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                                   prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                                   prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY), &
                                   cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                                   cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                                   cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                                   cons(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY)  &
                                  )
  END SUBROUTINE Convert2Conservative_facesub
  

  SUBROUTINE FinalizePhysics(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("ClosePhysics","not initialized")
    ! deallocate pointer variables used in all physics modules
    DEALLOCATE(this%amin,this%amax,this%bmin,this%bmax,this%cmin,this%cmax,this%tmp, &
         this%tmin,this%tmax,this%bccsound,this%fcsound, &
         this%pvarname,this%cvarname)
  END SUBROUTINE FinalizePhysics

END MODULE physics_base_mod
