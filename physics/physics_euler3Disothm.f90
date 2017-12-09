!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_euler3D_isothm.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2017                                                   #
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
  USE physics_base_mod
  USE mesh_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 4              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 3D isotherm"
  !--------------------------------------------------------------------------!
  TYPE,  EXTENDS(physics_base) :: physics_euler3Dit
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

    PROCEDURE :: GeometricalSources_center
    PROCEDURE :: ExternalSources

!    PROCEDURE :: CalcIntermediateStateX_euler3Dit    ! for HLLC
!    PROCEDURE :: CalcIntermediateStateY_euler3Dit    ! for HLLC
!    PROCEDURE :: CalcCharSystemX_euler3Dit           ! for absorbing boundaries
!    PROCEDURE :: CalcCharSystemY_euler3Dit           ! for absorbing boundaries
!    PROCEDURE :: CalcBoundaryDataX_euler3Dit         ! for absorbing boundaries
!    PROCEDURE :: CalcBoundaryDataY_euler3Dit         ! for absorbing boundaries
!    PROCEDURE :: CalcPrim2RiemannX_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcPrim2RiemannY_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimX_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimY_euler3Dit         ! for farfield boundaries
!    PROCEDURE :: CalcRoeAverages_euler2Dit           ! for advanced wavespeeds
!    PROCEDURE :: CalcStresses_euler3Dit
!    PROCEDURE :: ViscositySources_euler3Dit
!    PROCEDURE :: ExternalSources_euler2Dit
!    PROCEDURE :: FargoSources_euler2Dit
!    PROCEDURE :: GeometricalSources_faces
!    PROCEDURE :: ReflectionMasks
!    PROCEDURE :: AxisMasks
!    PROCEDURE :: SetEigenValues_euler2Dit
!    PROCEDURE :: SetBoundaryData_euler3Dit
!    PROCEDURE :: SetCharVars_euler3Dit
!    PROCEDURE :: CalcFlux_euler2Dit, &
!    MomentumSourcesX_euler2Dit, &
!    MomentumSourcesY_euler2Dit, &
!    AddBackgroundVelocity_euler2Dit, &
!    SubtractBackgroundVelocity_euler2Dit, &

    FINAL     :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
   PUBLIC :: &
       ! types
       physics_euler3Dit
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
    INTEGER           :: err
    !------------------------------------------------------------------------!
    CALL this%InitPhysics(Mesh,config,IO,EULER3D_ISOTH,problem_name,num_var)
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
!CDIR IEXPAND
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
    INTEGER                                :: i,j,k
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
    INTEGER           :: i,j,k,l
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE UpdateSoundSpeed_faces

  !> Calculates wave speeds at cell-centers
  PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar,amin,amax,bmin,bmax,cmin,cmax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)    :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(OUT)   :: amin,amax,bmin,bmax,cmin,cmax
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
         DO i=Mesh%IGMIN,Mesh%IGMAX
          ! x-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%XVELOCITY),&
               amin(i,j,k),amax(i,j,k))
          ! y-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%YVELOCITY),&
               bmin(i,j,k),bmax(i,j,k))
          ! z-direction
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%bccsound(i,j,k),pvar(i,j,k,this%ZVELOCITY),&
               cmin(i,j,k),cmax(i,j,k))
         END DO
      END DO
    END DO
  END SUBROUTINE CalcWaveSpeeds_center

  !> Calculates wave speeds at cell-faces
  PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim,cons,amin,amax,bmin,bmax,cmin,cmax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)    :: prim,cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(OUT)   :: amin,amax,bmin,bmax,cmin,cmax
    !------------------------------------------------------------------------!
    REAL                                    :: uRoe, aminRoe, amaxRoe
    INTEGER                                 :: i,j,k
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell interfaces
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! western
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,1), &
                             prim(i,j,k,1,this%XVELOCITY), &
                             this%tmp(i,j,k),this%tmp1(i,j,k))
          ! eastern
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,2), &
                             prim(i,j,k,2,this%XVELOCITY), &
                             amin(i,j,k),amax(i,j,k))
          ! southern
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,3), &
                             prim(i,j,k,3,this%YVELOCITY), &
                             this%tmp2(i,j,k),this%tmp3(i,j,k))
          ! northern
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,4), &
                             prim(i,j,k,4,this%YVELOCITY), &
                             bmin(i,j,k),bmax(i,j,k))
          ! bottom
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,5), &
                             prim(i,j,k,5,this%ZVELOCITY), &
                             this%tmp4(i,j,k),this%tmp5(i,j,k))
          ! top
!CDIR IEXPAND
          CALL SetWavespeeds(this%fcsound(i,j,k,6), &
                             prim(i,j,k,6,this%ZVELOCITY), &
                             cmin(i,j,k),cmax(i,j,k))
        END DO
      END DO
    END DO
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
    IF (this%advanced_wave_speeds) THEN
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
        DO i=Mesh%IMIN-1,Mesh%IMAX
          ! western & eastern interfaces
          ! get Roe averaged x-velocity
!CDIR IEXPAND
          CALL SetRoeAverages(prim(i,j,k,2,this%DENSITY),prim(i+1,j,k,1,this%DENSITY), &
                              prim(i,j,k,2,this%XVELOCITY),prim(i+1,j,k,1,this%XVELOCITY), &
                              uRoe)
          ! compute Roe averaged wave speeds
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,2),uRoe,aminRoe,amaxRoe)
          amin(i,j,k) = MIN(aminRoe,this%tmp(i+1,j,k),amin(i,j,k))
          amax(i,j,k) = MAX(amaxRoe,this%tmp1(i+1,j,k),amax(i,j,k))
        END DO
      END DO
    END DO
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JMIN-1,Mesh%JMAX
!CDIR NODEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! southern & northern interfaces
          ! get Roe averaged y-velocity
!CDIR IEXPAND
          CALL SetRoeAverages(prim(i,j,k,4,this%DENSITY),prim(i,j+1,k,3,this%DENSITY), &
                              prim(i,j,k,4,this%YVELOCITY),prim(i,j+1,k,3,this%YVELOCITY), &
                              uRoe)
          ! compute Roe averaged wave speeds
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,4),uRoe,aminRoe,amaxRoe)
          bmin(i,j,k) = MIN(aminRoe,this%tmp2(i,j+1,k),bmin(i,j,k))
          bmax(i,j,k) = MAX(amaxRoe,this%tmp3(i,j+1,k),bmax(i,j,k))
        END DO
      END DO
    END DO
!CDIR COLLAPSE
    DO k=Mesh%KGMIN-1,Mesh%KGMAX
      DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! topper & bottomer interfaces
          ! get Roe averaged z-velocity
!CDIR IEXPAND
          CALL SetRoeAverages(prim(i,j,k,6,this%DENSITY),prim(i,j,k+1,5,this%DENSITY), &
                               prim(i,j,k,6,this%ZVELOCITY),prim(i,j,k+1,5,this%ZVELOCITY), &
                               uRoe)
          ! compute Roe averaged wave speeds
!CDIR IEXPAND
          CALL SetWaveSpeeds(this%fcsound(i,j,k,6),uRoe,aminRoe,amaxRoe)
          cmin(i,j,k) = MIN(aminRoe,this%tmp3(i,j,k+1),bmin(i,j,k))
          cmax(i,j,k) = MAX(amaxRoe,this%tmp4(i,j,k+1),bmax(i,j,k))
        END DO
      END DO
    END DO
    ELSE
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
        DO i=Mesh%IMIN-1,Mesh%IMAX
          ! western & eastern interfaces
          amin(i,j,k) = MIN(0.0,this%tmp(i+1,j,k) ,amin(i,j,k))
          amax(i,j,k) = MAX(0.0,this%tmp1(i+1,j,k),amax(i,j,k))
        END DO
      END DO
    END DO
!CDIR COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JMIN-1,Mesh%JMAX
!CDIR NODEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! southern & northern interfaces
          bmin(i,j,k) = MIN(0.0,this%tmp2(i,j+1,k),bmin(i,j,k))
          bmax(i,j,k) = MAX(0.0,this%tmp3(i,j+1,k),bmax(i,j,k))
        END DO
      END DO
    END DO
!CDIR COLLAPSE
    DO k=Mesh%KMIN-1,Mesh%KMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! bottom & top interfaces
          cmin(i,j,k) = MIN(0.0,this%tmp4(i,j,k+1),cmin(i,j,k))
          cmax(i,j,k) = MAX(0.0,this%tmp5(i,j,k+1),cmax(i,j,k))
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
    INTEGER                            :: i,j,k
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
    CALL SetFlux(this%fcsound(:,:,:,nmin:nmax), &
                 prim(:,:,:,nmin:nmax,this%DENSITY), &
                 prim(:,:,:,nmin:nmax,this%XVELOCITY), &
                 cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
                 cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 cons(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                 xfluxes(:,:,:,nmin:nmax,this%DENSITY), &
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
    CALL SetFlux(this%fcsound(:,:,:,nmin:nmax), &
                 prim(:,:,:,nmin:nmax,this%DENSITY),&
                 prim(:,:,:,nmin:nmax,this%YVELOCITY), &
                 cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
                 cons(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                 yfluxes(:,:,:,nmin:nmax,this%DENSITY), &
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
    CALL SetFlux(this%fcsound(:,:,:,nmin:nmax), &
                 prim(:,:,:,nmin:nmax,this%DENSITY), &
                 prim(:,:,:,nmin:nmax,this%ZVELOCITY), &
                 cons(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                 cons(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 cons(:,:,:,nmin:nmax,this%XMOMENTUM), &
                 zfluxes(:,:,:,nmin:nmax,this%DENSITY), &
                 zfluxes(:,:,:,nmin:nmax,this%ZMOMENTUM), &
                 zfluxes(:,:,:,nmin:nmax,this%YMOMENTUM), &
                 zfluxes(:,:,:,nmin:nmax,this%XMOMENTUM))
  END SUBROUTINE CalcFluxesZ

  ! momentum and energy sources due to external force
  PURE SUBROUTINE ExternalSources(this,Mesh,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
                              INTENT(IN)  :: accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: sterm
    !------------------------------------------------------------------------!
    INTEGER                         :: i,j,k
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

  !> Maks for reflecting boundaries
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

  !> momentum source terms due to inertial forces
  !! P is the isothermal pressure rho*cs*cs
  !!
  !! \todo the syntax of this part was changed during transition to 3D.
  !!       Have a look at physics_euler3D and revert if the solution here is not
  !!       so nice.
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

  !> Destructor
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(physics_euler3Dit), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizePhysics()
  END SUBROUTINE Finalize

END MODULE physics_euler3Dit_mod
