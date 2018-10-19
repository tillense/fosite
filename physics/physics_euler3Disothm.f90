!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_euler3Disothm.f90                                        #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
!# Lars Bösch       <lboesch@astrophysik.uni-kiel.de>                        #
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
!! - isothermal gas dynamics
!!   \key{cs,REAL,isothermal sound speed,0.0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!! \author Lars Boesch
!! \author Jannes Klee
!!
!! \brief basic module for 3D isothermal Euler equations
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler3Dit_mod
  USE physics_euler2dit_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER           :: num_var = 4          ! number of variables !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D isotherm"
  !--------------------------------------------------------------------------!
  TYPE,  EXTENDS(physics_euler2Dit) :: physics_euler3Dit
  CONTAINS
    PROCEDURE :: InitPhysics_euler3Dit
    !------Convert2Primitive-------!
    PROCEDURE :: Convert2Primitive_center
    PROCEDURE :: Convert2Primitive_centsub
    PROCEDURE :: Convert2Primitive_faces
    PROCEDURE :: Convert2Primitive_facesub
    !------Convert2Conservative----!
    PROCEDURE :: Convert2Conservative_center
    PROCEDURE :: Convert2Conservative_centsub
    PROCEDURE :: Convert2Conservative_faces
    PROCEDURE :: Convert2Conservative_facesub
    !------Soundspeed Routines-----!
    PROCEDURE :: UpdateSoundSpeed_center
    PROCEDURE :: UpdateSoundSpeed_faces
    !------Wavespeed Routines------!
    PROCEDURE :: CalcWaveSpeeds_center
    PROCEDURE :: CalcWaveSpeeds_faces
    !------Flux Routines-----------!
    PROCEDURE :: CalcFluxesX
    PROCEDURE :: CalcFluxesY
    PROCEDURE :: CalcFluxesZ
    !------fargo routines----------!
    PROCEDURE :: AddBackgroundVelocityX
    PROCEDURE :: SubtractBackgroundVelocityX
    PROCEDURE :: AddBackgroundVelocityY
    PROCEDURE :: SubtractBackgroundVelocityY
    PROCEDURE :: FargoSources

    PROCEDURE :: GeometricalSources_center
    PROCEDURE :: ExternalSources

    PROCEDURE :: ReflectionMasks                      ! for reflecting boundaries
!    PROCEDURE :: CalcIntermediateStateX_euler3Dit    ! for HLLC
!    PROCEDURE :: CalcIntermediateStateY_euler3Dit    ! for HLLC
    PROCEDURE :: CalculateCharSystemX           ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemY           ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemZ           ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataX         ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataY         ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataZ         ! for absorbing boundaries
!    PROCEDURE :: CalcPrim2RiemannX_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcPrim2RiemannY_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimX_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimY_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcRoeAverages_euler2Dit           ! for advanced wavespeeds
  !  GENERIC :: CalcStresses_euler3D => CalcStresses_euler3Dit
    PROCEDURE :: ViscositySources
    PROCEDURE :: ViscositySources_euler3Dit
    PROCEDURE :: CalcStresses_euler
!    PROCEDURE :: ExternalSources_euler2Dit
!    PROCEDURE :: FargoSources_euler2Dit
!    PROCEDURE :: GeometricalSources_faces
    PROCEDURE :: AxisMasks
!    PROCEDURE :: CalcFlux_euler2Dit, &
    PROCEDURE     :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
   PUBLIC :: &
       ! types
       physics_euler3Dit, &
       SetEigenValues, &
       SetBoundaryData, &
       SetCharVars

  !--------------------------------------------------------------------------!

CONTAINS

  !> Intialization of isothermal physics
  !!
  !! - calls intialization of base routines of physics
  !! - set array indices, names and number of dimensions
  SUBROUTINE InitPhysics_euler3Dit(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    TYPE(Dict_TYP), POINTER,  INTENT(IN)    :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: err
    !------------------------------------------------------------------------!
    CALL this%InitPhysics(Mesh,config,IO,EULER3D_ISOTHERM,problem_name,num_var)

    ! isothermal sound speed
    CALL GetAttr(config, "cs", this%csiso, 0.0)
   
    ! set array indices
    this%DENSITY   = 1                                 ! mass density        !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    this%ZVELOCITY = 4                                 ! no z-velocity       !
    this%ZMOMENTUM = 4                                 ! no z-momentum       !
    this%PRESSURE  = 0                                 ! no pressure         !
    this%ENERGY    = 0                                 ! no total energy     !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%pvarname(this%ZVELOCITY) = "zvelocity"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%cvarname(this%ZMOMENTUM) = "zmomentum"
    this%DIM = 3
        
    ! allocate memory for arrays common to all physics modules
    ALLOCATE(this%bccsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),            &
             this%fcsound(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces), &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitPhysics_euler3dit", "Unable to allocate memory.")

    IF(this%csiso.GT.0.) THEN
      this%bccsound(:,:,:)  = this%csiso
      this%fcsound(:,:,:,:) = this%csiso
    ELSE
      this%bccsound(:,:,:)  = 0.
      this%fcsound(:,:,:,:) = 0.
    END IF

  END SUBROUTINE InitPhysics_euler3Dit

  !> Converts to primitives at cell centers
  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL this%Convert2Primitive_centsub(Mesh,Mesh%IGMIN,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,cvar,pvar)
  END SUBROUTINE Convert2Primitive_center

  !> Converts to primitive variables at cell centers
  PURE SUBROUTINE Convert2Primitive_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    CALL Cons2Prim(cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                                cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                                cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                                cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM), &
                                pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                                pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                                pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                                pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY)  &
                               )
  END SUBROUTINE Convert2Primitive_centsub

  !> Converts to conservative variables at faces
  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL this%Convert2Primitive_facesub(Mesh,Mesh%IGMIN,Mesh%IGMAX, &
                                             Mesh%JGMIN,Mesh%JGMAX, &
                                             Mesh%KGMIN,Mesh%KGMAX, &
                                        cons,prim)
  END SUBROUTINE Convert2Primitive_faces

  !> Converts to conservative variables at faces
  PURE SUBROUTINE Convert2Primitive_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    CALL Cons2Prim(cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                   cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                   prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY)  &
                   )
  END SUBROUTINE Convert2Primitive_facesub

  !> Convert from primtive to conservative variables at cell-centers
  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL this%Convert2Conservative_centsub(Mesh,Mesh%IGMIN,Mesh%IGMAX, &
                                                Mesh%JGMIN,Mesh%JGMAX, &
                                                Mesh%KGMIN,Mesh%KGMAX, &
                                           pvar,cvar)
  END SUBROUTINE Convert2Conservative_center


  !> Convert from primtive to conservative variables at cell-centers
  PURE SUBROUTINE Convert2Conservative_centsub(this,Mesh,i1,i2,j1,j2,k1,k2,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Prim2Cons(pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                   pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                   pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                   pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY), &
                   cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                   cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                   cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                   cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM)  &
                   )
  END SUBROUTINE Convert2Conservative_centsub


  !> Convert to from primitve to conservative variables at cell-faces
  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: prim
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL this%Convert2Conservative_facesub(Mesh,Mesh%IGMIN,Mesh%IGMAX, &
                                                Mesh%JGMIN,Mesh%JGMAX, &
                                                Mesh%KGMIN,Mesh%KGMAX, &
                                           prim,cons)
  END SUBROUTINE Convert2Conservative_faces


  !> Convert to from primitve to conservative variables at cell-faces
  PURE SUBROUTINE Convert2Conservative_facesub(this,Mesh,i1,i2,j1,j2,k1,k2,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: prim
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    CALL Prim2Cons(prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                   prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                   prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                   cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                   cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM)  &
                   )
  END SUBROUTINE Convert2Conservative_facesub

  !> Empty routine for isothermal simulation
  !!
  !! Will be overwritten in physics with energy equation
  !! \todo Have a look for a nicer solution of this issue
  PURE SUBROUTINE UpdateSoundSpeed_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE UpdateSoundSpeed_center

  !> Empty routine for isothermal simulation
  !!
  !! Will be overwritten in physics with energy equation
  !! \todo Have a look for a nicer solution of this issue
  PURE SUBROUTINE UpdateSoundSpeed_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)    :: prim
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE UpdateSoundSpeed_faces

  !> Calculates wave speeds at cell-centers
  PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)    :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS), &
                              INTENT(OUT)   :: minwav,maxwav
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
         DO i=Mesh%IGMIN,Mesh%IGMAX
          ! x-direction
          CALL SetWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%XVELOCITY),&
               minwav(i,j,k,1),maxwav(i,j,k,1))
          ! y-direction
          CALL SetWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%YVELOCITY),&
               minwav(i,j,k,2),maxwav(i,j,k,2))
          ! z-direction
          CALL SetWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%ZVELOCITY),&
               minwav(i,j,k,3),maxwav(i,j,k,3))
         END DO
      END DO
    END DO
  END SUBROUTINE CalcWaveSpeeds_center

  !> Calculates wave speeds at cell-faces
  PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim,cons,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)    :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS), &
                              INTENT(OUT)   :: minwav,maxwav
    !------------------------------------------------------------------------!
    REAL                                    :: uRoe, aminRoe, amaxRoe
    INTEGER                                 :: i,j,k
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell interfaces
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! western
          CALL SetWaveSpeeds(this%fcsound(i,j,k,1), &
                             prim(i,j,k,1,this%XVELOCITY), &
                             this%tmp(i,j,k),this%tmp1(i,j,k))
          ! eastern
          CALL SetWaveSpeeds(this%fcsound(i,j,k,2), &
                             prim(i,j,k,2,this%XVELOCITY), &
                             minwav(i,j,k,1),maxwav(i,j,k,1))
          ! southern
          CALL SetWaveSpeeds(this%fcsound(i,j,k,3), &
                             prim(i,j,k,3,this%YVELOCITY), &
                             this%tmp2(i,j,k),this%tmp3(i,j,k))
          ! northern
          CALL SetWaveSpeeds(this%fcsound(i,j,k,4), &
                             prim(i,j,k,4,this%YVELOCITY), &
                             minwav(i,j,k,2),maxwav(i,j,k,2))
          ! bottom
          CALL SetWaveSpeeds(this%fcsound(i,j,k,5), &
                             prim(i,j,k,5,this%ZVELOCITY), &
                             this%tmp4(i,j,k),this%tmp5(i,j,k))
          ! top
          CALL SetWavespeeds(this%fcsound(i,j,k,6), &
                             prim(i,j,k,6,this%ZVELOCITY), &
                             minwav(i,j,k,3),maxwav(i,j,k,3))
        END DO
      END DO
    END DO
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
    IF (this%advanced_wave_speeds) THEN
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IMIN-1,Mesh%IMAX
          ! western & eastern interfaces
          ! get Roe averaged x-velocity
          CALL SetRoeAverages(prim(i,j,k,2,this%DENSITY),prim(i+1,j,k,1,this%DENSITY), &
                              prim(i,j,k,2,this%XVELOCITY),prim(i+1,j,k,1,this%XVELOCITY), &
                              uRoe)
          ! compute Roe averaged wave speeds
          CALL SetWaveSpeeds(this%fcsound(i,j,k,2),uRoe,aminRoe,amaxRoe)
          minwav(i,j,k,1) = MIN(aminRoe,this%tmp(i+1,j,k),minwav(i,j,k,1))
          maxwav(i,j,k,1) = MAX(amaxRoe,this%tmp1(i+1,j,k),maxwav(i,j,k,1))
        END DO
      END DO
    END DO
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JMIN-1,Mesh%JMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! southern & northern interfaces
          ! get Roe averaged y-velocity
          CALL SetRoeAverages(prim(i,j,k,4,this%DENSITY),prim(i,j+1,k,3,this%DENSITY), &
                              prim(i,j,k,4,this%YVELOCITY),prim(i,j+1,k,3,this%YVELOCITY), &
                              uRoe)
          ! compute Roe averaged wave speeds
          CALL SetWaveSpeeds(this%fcsound(i,j,k,4),uRoe,aminRoe,amaxRoe)
          minwav(i,j,k,2) = MIN(aminRoe,this%tmp2(i,j+1,k),minwav(i,j,k,2))
          maxwav(i,j,k,2) = MAX(amaxRoe,this%tmp3(i,j+1,k),maxwav(i,j,k,2))
        END DO
      END DO
    END DO
    DO k=Mesh%KMIN-1,Mesh%KMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! topper & bottomer interfaces
          ! get Roe averaged z-velocity
          CALL SetRoeAverages(prim(i,j,k,6,this%DENSITY),prim(i,j,k+1,5,this%DENSITY), &
                               prim(i,j,k,6,this%ZVELOCITY),prim(i,j,k+1,5,this%ZVELOCITY), &
                               uRoe)
          ! compute Roe averaged wave speeds
          CALL SetWaveSpeeds(this%fcsound(i,j,k,6),uRoe,aminRoe,amaxRoe)
          minwav(i,j,k,3) = MIN(aminRoe,this%tmp4(i,j,k+1),minwav(i,j,k,3))
          maxwav(i,j,k,3) = MAX(amaxRoe,this%tmp5(i,j,k+1),maxwav(i,j,k,3))
        END DO
      END DO
    END DO
    ELSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IMIN-1,Mesh%IMAX
          ! western & eastern interfaces
            minwav(i,j,k,1) = MIN(0.0,this%tmp(i+1,j,k) ,minwav(i,j,k,1))
            maxwav(i,j,k,1) = MAX(0.0,this%tmp1(i+1,j,k),maxwav(i,j,k,1))
          END DO
        END DO
      END DO
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JMIN-1,Mesh%JMAX
  !NEC$ IVDEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! southern & northern interfaces
            minwav(i,j,k,2) = MIN(0.0,this%tmp2(i,j+1,k),minwav(i,j,k,2))
            maxwav(i,j,k,2) = MAX(0.0,this%tmp3(i,j+1,k),maxwav(i,j,k,2))
          END DO
        END DO
      END DO
      DO k=Mesh%KMIN-1,Mesh%KMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
  !NEC$ IVDEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! bottom & top interfaces
          minwav(i,j,k,3) = MIN(0.0,this%tmp4(i,j,k+1),minwav(i,j,k,3))
          maxwav(i,j,k,3) = MAX(0.0,this%tmp5(i,j,k+1),maxwav(i,j,k,3))
        END DO
      END DO
    END DO
    END IF
  END SUBROUTINE CalcWaveSpeeds_faces

  !> Calculates geometrical sources at cell-center
  PURE SUBROUTINE GeometricalSources_center(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)    :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT)   :: sterm
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k
    !------------------------------------------------------------------------!
    ! compute geometrical source only for non-cartesian mesh except for the
    ! EULER2D_IAMROT case for which geometrical sources are always necessary.
    IF ((Mesh%Geometry%GetType().NE.CARTESIAN).OR. &
        (this%GetType().EQ.EULER2D_IAMROT).OR.     &
        (this%GetType().EQ.EULER2D_ISOIAMROT)) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            CALL CalcGeometricalSources(cvar(i,j,k,this%XMOMENTUM),                       &
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

  !> Calculates geometrical sources at cell-faces
!  PURE SUBROUTINE GeometricalSources_faces(this,Mesh,prim,cons,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
!         :: prim,cons
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
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


  !> Calculate Fluxes in x-direction
  PURE SUBROUTINE CalcFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: xfluxes
    !------------------------------------------------------------------------!
    CALL SetFlux(this%fcsound(:,:,:,nmin:nmax),           &
                 prim(:,:,:,nmin:nmax,this%DENSITY),      &
                 prim(:,:,:,nmin:nmax,this%XVELOCITY),    &
                 cons(:,:,:,nmin:nmax,this%XMOMENTUM),    &
                 cons(:,:,:,nmin:nmax,this%YMOMENTUM),    &
                 cons(:,:,:,nmin:nmax,this%ZMOMENTUM),    &
                 xfluxes(:,:,:,nmin:nmax,this%DENSITY),   &
                 xfluxes(:,:,:,nmin:nmax,this%XMOMENTUM), &
                 xfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 xfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM))
  END SUBROUTINE CalcFluxesX

  !> Calculate Fluxes in y-direction
  PURE SUBROUTINE CalcFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: yfluxes
    !------------------------------------------------------------------------!
    CALL SetFlux(this%fcsound(:,:,:,nmin:nmax),           &
                 prim(:,:,:,nmin:nmax,this%DENSITY),      &
                 prim(:,:,:,nmin:nmax,this%YVELOCITY),    &
                 cons(:,:,:,nmin:nmax,this%YMOMENTUM),    &
                 cons(:,:,:,nmin:nmax,this%XMOMENTUM),    &
                 cons(:,:,:,nmin:nmax,this%ZMOMENTUM),    &
                 yfluxes(:,:,:,nmin:nmax,this%DENSITY),   &
                 yfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 yfluxes(:,:,:,nmin:nmax,this%XMOMENTUM), &
                 yfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM))
  END SUBROUTINE CalcFluxesY

  !> Calculate Fluxes in z-direction
  PURE SUBROUTINE CalcFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: nmin,nmax
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: zfluxes
    !------------------------------------------------------------------------!
    CALL SetFlux(this%fcsound(:,:,:,nmin:nmax),           &
                 prim(:,:,:,nmin:nmax,this%DENSITY),      &
                 prim(:,:,:,nmin:nmax,this%ZVELOCITY),    &
                 cons(:,:,:,nmin:nmax,this%ZMOMENTUM),    &
                 cons(:,:,:,nmin:nmax,this%YMOMENTUM),    &
                 cons(:,:,:,nmin:nmax,this%XMOMENTUM),    &
                 zfluxes(:,:,:,nmin:nmax,this%DENSITY),   &
                 zfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                 zfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 zfluxes(:,:,:,nmin:nmax,this%XMOMENTUM))
  END SUBROUTINE CalcFluxesZ

  PURE SUBROUTINE AddBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                              INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_yvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             pvar(i,j,k,this%YVELOCITY) = pvar(i,j,k,this%YVELOCITY) + w(i,k)
             cvar(i,j,k,this%YMOMENTUM) = cvar(i,j,k,this%YMOMENTUM) &
                                      + cvar(i,j,k,this%DENSITY)*w(i,k)
          END DO
        END DO
      END DO
      this%transformed_yvelocity = .FALSE.
    END IF
  END SUBROUTINE AddBackgroundVelocityY


  PURE SUBROUTINE AddBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                              INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_xvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             pvar(i,j,k,this%XVELOCITY) = pvar(i,j,k,this%XVELOCITY) + w(j,k)
             cvar(i,j,k,this%XMOMENTUM) = cvar(i,j,k,this%XMOMENTUM) &
                                      + cvar(i,j,k,this%DENSITY)*w(j,k)
          END DO
        END DO
      END DO
      this%transformed_xvelocity = .FALSE.
    END IF
  END SUBROUTINE AddBackgroundVelocityX


  PURE SUBROUTINE SubtractBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                              INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_yvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            pvar(i,j,k,this%YVELOCITY) = pvar(i,j,k,this%YVELOCITY) - w(i,k)
            cvar(i,j,k,this%YMOMENTUM) = cvar(i,j,k,this%YMOMENTUM) &
                                     - cvar(i,j,k,this%DENSITY)*w(i,k)
          END DO
        END DO
      END DO
      this%transformed_yvelocity = .TRUE.
    END IF
  END SUBROUTINE SubtractBackgroundVelocityY

  PURE SUBROUTINE SubtractBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                              INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_xvelocity) THEN
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            pvar(i,j,k,this%XVELOCITY) = pvar(i,j,k,this%XVELOCITY) - w(j,k)
            cvar(i,j,k,this%XMOMENTUM) = cvar(i,j,k,this%XMOMENTUM) &
                                     - cvar(i,j,k,this%DENSITY)*w(j,k)
          END DO
        END DO
      END DO
      this%transformed_xvelocity = .TRUE.
    END IF
  END SUBROUTINE SubtractBackgroundVelocityX


  !> \public sources terms for fargo advection
  !!
  !! If the background velocity \f$\vec{w}=w\,\hat{e}_\eta\f$ with
  !! \f$w\f$ independent of \f$\eta\f$ and \f$t\f$ is subtracted from
  !! the overall velocity of the flow, an additional source term occurs
  !! in the \f$\eta\f$-momentum equation:
  !! \f[
  !!     S_\mathrm{Fargo} = -\varrho u_\xi \frac{1}{h_\xi} \partial_\xi \left(h_\xi w\right)
  !!                      = -\varrho u_\xi w \,\partial_\xi \left(\ln{|h_\xi w|}\right)
  !! \f]
  PURE SUBROUTINE FargoSources(this,Mesh,w,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)    :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                              INTENT(IN)    :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM+this%PNUM), &
                              INTENT(INOUT) :: sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    DO k=Mesh%KMIN,Mesh%KMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
!            sterm(i,j,this%XMOMENTUM) = cvar(i,j,this%YMOMENTUM) * w(i) &
!                   * Mesh%cyxy%bcenter(i,j)
           sterm(i,j,k,this%YMOMENTUM) = -cvar(i,j,k,this%XMOMENTUM) &
                  * 0.5 * (w(i+1,k)-w(i-1,k)) / Mesh%dlx(i,j,k)
        END DO
      END DO
    END DO
  END SUBROUTINE FargoSources


  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NDIMS), &
                              INTENT(IN)  :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: sterm
    !------------------------------------------------------------------------!
    INTEGER                         :: i,j,k
    !------------------------------------------------------------------------!
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

  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)      :: Mesh
    REAL,                   INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                           :: pvar
    REAL,                   INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                                           :: btxx,btxy,btxz,btyy,btyz,btzz
    REAL,                   INTENT(OUT), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: sterm
   !------------------------------------------------------------------------!
   CALL this%ViscositySources_euler3Dit(Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
 
   !compute scalar product of v and tau (x-component)
   this%tmp(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxx(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btxy(:,:,:) &
                    + pvar(:,:,:,this%ZVELOCITY)*btxz(:,:,:)

   !compute scalar product of v and tau (y-component)
   this%tmp1(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxy(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btyy(:,:,:) &
                    + pvar(:,:,:,this%ZVELOCITY)*btyz(:,:,:)

   !compute scalar product of v and tau (z-component)
   this%tmp2(:,:,:) = pvar(:,:,:,this%XVELOCITY)*btxz(:,:,:) &
                    + pvar(:,:,:,this%YVELOCITY)*btyz(:,:,:) &
                    + pvar(:,:,:,this%ZVELOCITY)*btzz(:,:,:)
   ! compute vector divergence of scalar product v and tau
   CALL Mesh%Divergence(this%tmp(:,:,:),this%tmp1(:,:,:), &
      this%tmp2(:,:,:),sterm(:,:,:,this%ENERGY))
 END SUBROUTINE ViscositySources

  !> Characteristic variables for absorbing boundary conditions
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit),INTENT(IN)   :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: i,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i1,i2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at i
    CALL SetEigenValues( &
          this%bccsound(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
    ! compute characteristic variables
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
    CALL SetCharVars( &
          this%fcsound(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,WEST), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
  END SUBROUTINE CalculateCharSystemX


  !> Characteristic variables for absorbing boundary conditions
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: j,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, DIMENSION(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
    INTEGER           :: j1,j2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at j
    CALL SetEigenValues( &
          this%bccsound(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
          lambda(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4))
    ! compute characteristic variables
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
    CALL SetCharVars( &
          this%fcsound(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,SOUTH), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4))
  END SUBROUTINE CalculateCharSystemY


  !> Characteristic variables for absorbing boundary conditions
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateCharSystemZ(this,Mesh,k,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)    :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    INTEGER,                INTENT(IN)    :: k,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                          :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), &
                            INTENT(OUT)   :: xvar
    !------------------------------------------------------------------------!
    INTEGER           :: k1,k2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at k
    CALL SetEigenValues( &
          this%bccsound(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%ZVELOCITY), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4))
    ! compute characteristic variables
    k1 = k + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    k2 = MAX(k,k1)
    k1 = MIN(k,k1)
    CALL SetCharVars( &
          this%fcsound(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,BOTTOM), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%YVELOCITY), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4))
  END SUBROUTINE CalculateCharSystemZ


  !> Calculate boundary data for absorbing boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3dit), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: i1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER                            :: i2,fidx
    !------------------------------------------------------------------------!
    i2 = i1 + SIGN(1,dir)  ! i +/- 1 depending on the sign of dir
    IF (i2.LT.i1) THEN
      fidx = 1  !WEST
    ELSE
      fidx = 2  !EAST
    END IF
    CALL SetBoundaryData( &
          this%fcsound(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,fidx), &
          1.0*SIGN(1,dir), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY))
  END SUBROUTINE CalculateBoundaryDataX


  !> Calculate boundary data for absorbing boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3dit), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: j1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%KGMIN:Mesh%KMAX,Mesh%IGMIN:Mesh%IGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER                            :: j2,fidx
    !------------------------------------------------------------------------!
    j2 = j1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
    IF (j2.LT.j1) THEN
      fidx = 3  !SOUTH
    ELSE
      fidx = 4  !NORTH
    END IF
    CALL SetBoundaryData( &
          this%fcsound(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,fidx), &
          1.0*SIGN(1,dir), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,3), &
          xvar(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,4), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY))
  END SUBROUTINE CalculateBoundaryDataY


  !> Calculate boundary data for absorbing boundaries
  !\todo NOT VERIFIED
  PURE SUBROUTINE CalculateBoundaryDataZ(this,Mesh,k1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3dit), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    INTEGER,                INTENT(IN) :: k1,dir
    REAL,                   INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
                                       :: xvar
    REAL,                   INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER                            :: k2,fidx
    !------------------------------------------------------------------------!
    k2 = k1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
    IF (k2.LT.k1) THEN
      fidx = 5  !BOTTOM
    ELSE
      fidx = 6  !TOP
    END IF
    CALL SetBoundaryData( &
          this%fcsound(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,fidx), &
          1.0*SIGN(1,dir), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%YVELOCITY), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,4), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%ZVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%YVELOCITY))
  END SUBROUTINE CalculateBoundaryDataZ

  ELEMENTAL SUBROUTINE SetCharVars(cs,rho1,rho2,u1,u2,v1,v2,w1,w2,&
       xvar1,xvar2,xvar3,xvar4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho1,rho2,u1,u2,v1,v2,w1,w2
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3,xvar4
    !------------------------------------------------------------------------!
    REAL :: dlnrho,du
    !------------------------------------------------------------------------!
    dlnrho = LOG(rho2/rho1)
    du = u2-u1
    ! characteristic variables
    xvar1 = cs*dlnrho - du
    xvar2 = v2-v1
    xvar3 = cs*dlnrho + du
    xvar4 = w2-w1
  END SUBROUTINE SetCharVars


  ELEMENTAL SUBROUTINE SetBoundaryData(cs,dir,rho1,u1,v1,w1,xvar1,xvar2, &
       xvar3,xvar4,rho2,u2,v2,w2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,dir,rho1,u1,v1,w1,xvar1,xvar2,xvar3,xvar4
    REAL, INTENT(OUT) :: rho2,u2,v2,w2
    !------------------------------------------------------------------------!
    ! extrapolate boundary values using characteristic variables
    rho2 = rho1 * EXP(dir*0.5*(xvar3+xvar1)/cs)
    u2   = u1 + dir*0.5*(xvar3-xvar1)
    v2   = v1 + dir*xvar2
    w2   = w1 + dir*xvar4
  END SUBROUTINE SetBoundaryData




  ! identical to isothermal case
  PURE SUBROUTINE ViscositySources_euler3Dit(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
   !------------------------------------------------------------------------!
    CLASS(Physics_euler3Dit),INTENT(IN)  :: this
    CLASS(Mesh_base),INTENT(IN)        :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: &
          pvar,sterm
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: &
                    btxx,btxy,btxz,btyy,btyz,btzz
   !------------------------------------------------------------------------!
   !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,btxx,btxy,btxz,btyy,btyz,btzz
    INTENT(OUT)       :: sterm
   !------------------------------------------------------------------------!
   ! mean values of stress tensor components across the cell interfaces

   ! viscosity source terms
    sterm(:,:,:,this%DENSITY) = 0.0 

   ! compute viscous momentum sources
   ! divergence of stress tensor with symmetry btyx=btxy
    CALL Mesh%Divergence(btxx,btxy,btxz,btxy,btyy,btyz,btxz,btyz,btzz,sterm(:,:,:,this%XMOMENTUM), &
                         sterm(:,:,:,this%YMOMENTUM),sterm(:,:,:,this%ZMOMENTUM))
  END SUBROUTINE ViscositySources_euler3Dit

  ! identical to isothermal case. 
  PURE SUBROUTINE CalcStresses_euler(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(Mesh_base), INTENT(IN)          :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: &
         dynvis,bulkvis,btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,dynvis,bulkvis
    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    ! compute components of the stress tensor at cell bary centers
    ! inside the computational domain including one slice of ghost cells

    ! compute bulk viscosity first and store the result in this%tmp
    CALL Mesh%Divergence(pvar(:,:,:,this%XVELOCITY),pvar(:,:,:,this%YVELOCITY),pvar(:,:,:,this%ZVELOCITY),this%tmp(:,:,:))
    this%tmp(:,:,:) = bulkvis(:,:,:)*this%tmp(:,:,:)

!NEC$ OUTERLOOP_UNROLL(8)
  DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
    DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!NEC$ IVDEP
       DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
          ! compute the diagonal elements of the stress tensor
          btxx(i,j,k) = dynvis(i,j,k) * &
                ((pvar(i+1,j,k,this%XVELOCITY) - pvar(i-1,j,k,this%XVELOCITY)) / Mesh%dlx(i,j,k) &
               + 2.0 * Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) &
               + 2.0 * Mesh%cxzx%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) ) &
               + this%tmp(i,j,k)

          btyy(i,j,k) = dynvis(i,j,k) * &
               ( (pvar(i,j+1,k,this%YVELOCITY) - pvar(i,j-1,k,this%YVELOCITY)) / Mesh%dly(i,j,k) &
               + 2.0 * Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               + 2.0 * Mesh%cyzy%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) ) &
               + this%tmp(i,j,k)

          btzz(i,j,k) = dynvis(i,j,k) * &
               ( (pvar(i,j,k+1,this%ZVELOCITY) - pvar(i,j,k-1,this%ZVELOCITY)) / Mesh%dlz(i,j,k) &
               + 2.0 * Mesh%czxz%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               + 2.0 * Mesh%czyz%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) ) &
               + this%tmp(i,j,k)

          ! compute the off-diagonal elements (no bulk viscosity)
          btxy(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i+1,j,k,this%YVELOCITY) - pvar(i-1,j,k,this%YVELOCITY)) / Mesh%dlx(i,j,k) &
               + (pvar(i,j+1,k,this%XVELOCITY) - pvar(i,j-1,k,this%XVELOCITY)) / Mesh%dly(i,j,k) ) &
               - Mesh%cxyx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) &
               - Mesh%cyxy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

          btxz(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i+1,j,k,this%ZVELOCITY) - pvar(i-1,j,k,this%ZVELOCITY)) / Mesh%dlx(i,j,k) &
               + (pvar(i,j,k+1,this%XVELOCITY) - pvar(i,j,k-1,this%XVELOCITY)) / Mesh%dlz(i,j,k) ) &
               - Mesh%czxz%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) &
               - Mesh%cxzx%bcenter(i,j,k) * pvar(i,j,k,this%XVELOCITY) )

          btyz(i,j,k) = dynvis(i,j,k) * ( 0.5 * &
               ( (pvar(i,j,k+1,this%YVELOCITY) - pvar(i,j,k-1,this%YVELOCITY)) / Mesh%dlz(i,j,k) &
               + (pvar(i,j+1,k,this%ZVELOCITY) - pvar(i,j-1,k,this%ZVELOCITY)) / Mesh%dly(i,j,k) ) &
               - Mesh%czyz%bcenter(i,j,k) * pvar(i,j,k,this%ZVELOCITY) &
               - Mesh%cyzy%bcenter(i,j,k) * pvar(i,j,k,this%YVELOCITY) )

       END DO
    END DO
  END DO
  END SUBROUTINE CalcStresses_euler

  !> Mask for reflecting boundaries
  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit),      INTENT(IN)  :: this
    LOGICAL, DIMENSION(this%VNUM), INTENT(OUT) :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .FALSE.
    reflX(this%ZVELOCITY) = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .FALSE.
    ! top / bottom boundary
    reflZ(this%DENSITY)   = .FALSE.
    reflZ(this%XVELOCITY) = .FALSE.
    reflZ(this%YVELOCITY) = .FALSE.
    reflZ(this%ZVELOCITY) = .TRUE.
  END SUBROUTINE ReflectionMasks

  ! TODO: \warning not clear since 3D version if this is correct. Most probably
  ! axis boundaries can be applied always in two dimensions. Now only x-y plane
  PURE SUBROUTINE AxisMasks(this,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
    LOGICAL, DIMENSION(this%VNUM), &
                            INTENT(OUT) :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .TRUE.
    reflX(this%ZVELOCITY) = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
    reflY(this%ZVELOCITY) = .FALSE.
    ! bottomer / topper boundary
    reflZ(this%DENSITY)   = .FALSE.
    reflZ(this%XVELOCITY) = .FALSE.
    reflZ(this%YVELOCITY) = .FALSE.
    reflZ(this%ZVELOCITY) = .FALSE.
  END SUBROUTINE AxisMasks

  ELEMENTAL SUBROUTINE SetEigenValues(cs,v,l1,l2,l3,l4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: l1,l2,l3,l4
    !------------------------------------------------------------------------!
    ! all eigenvalues of the isothermal euler problem
    l1 = v - cs
    l2 = v
    l3 = v
    l4 = v + cs
  END SUBROUTINE SetEigenValues

  !> set minimal and maximal wave speeds
  ELEMENTAL SUBROUTINE SetWaveSpeeds(cs,v,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: minwav,maxwav
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    minwav = MIN(0.,v-cs)
    maxwav = MAX(0.,v+cs)
  END SUBROUTINE SetWaveSpeeds

  !> Set fluxes
  !!
  !! non-global elemtal subroutine
  ELEMENTAL SUBROUTINE SetFlux(cs,rho,v,m1,m2,m3,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,v,m1,m2,m3
    REAL, INTENT(OUT) :: f1, f2, f3, f4
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + rho*cs*cs
    f3 = m2*v
    f4 = m3*v
  END SUBROUTINE SetFlux

  !> Convert to from conservative to primitive variables
  !!
  !! non-global elemental routine
  ELEMENTAL SUBROUTINE Cons2Prim(rho_in,mu,mv,mw,rho_out,u,v,w)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,mu,mv,mw
    REAL, INTENT(OUT) :: rho_out,u,v,w
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u       = mu * inv_rho
    v       = mv * inv_rho
    w       = mw * inv_rho
  END SUBROUTINE Cons2Prim

  !> Convert to from primitive to conservative variables at cell-faces
  !!
  !! non-global elemental routine
  ELEMENTAL SUBROUTINE Prim2Cons(rho_in,u,v,w,rho_out,mu,mv,mw)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u,v,w
    REAL, INTENT(OUT) :: rho_out,mu,mv,mw
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    mw = rho_in * w
  END SUBROUTINE Prim2Cons

  ! TODO: Not verified
  !!
  !! non-global elemental routine
  ELEMENTAL SUBROUTINE SetRoeAverages(rhoL,rhoR,vL,vR,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rhoL,rhoR,vL,vR
    REAL, INTENT(OUT) :: v
    !------------------------------------------------------------------------!
    REAL :: sqrtrhoL,sqrtrhoR
    !------------------------------------------------------------------------!
    sqrtrhoL = SQRT(rhoL)
    sqrtrhoR = SQRT(rhoR)
    v = 0.5*(sqrtrhoL*vL + sqrtrhoR*vR) / (sqrtrhoL + sqrtrhoR)
  END SUBROUTINE SetRoeAverages


  !> momentum source terms due to inertial forces
  !! P is the isothermal pressure rho*cs*cs
  ELEMENTAL SUBROUTINE CalcGeometricalSources(mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz,srho,smx,smy,smz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz
    REAL, INTENT(OUT) :: srho, smx, smy, smz
    !------------------------------------------------------------------------!
    srho =  0.
    smx  = -my * (cxyx * vx - cyxy * vy) + mz * (czxz * vz - cxzx * vx) + (cyxy + czxz) * P
    smy  =  mx * (cxyx * vx - cyxy * vy) + mz * (czyz * vz - cyzy * vy) + (cxyx + czyz) * P
    smz  =  mx * (cxzx * vx - czxz * vz) + my * (cyzy * vy - czyz * vz) + (cxzx + cyzy) * P
  END SUBROUTINE CalcGeometricalSources


  !> Destructor
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE physics_euler3Dit_mod
