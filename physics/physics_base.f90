!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_base.f90                                                  #
!#                                                                           #
!# Copyright (C) 2007 - 2018                                                 #
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
!> \defgroup physics Physics
!! \{
!! \brief Family of physics modules
!! \}
!----------------------------------------------------------------------------!
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
!> \brief Basic physics module
!!
!! \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! This module provides the abstract class for all physics classes.
!!
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_base_mod
  USE constants_generic_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> named integer constants for flavour of state vectors
  ENUM, BIND(C)
    ENUMERATOR :: UNDEFINED=0, PRIMITIVE=1, CONSERVATIVE=2
  END ENUM
  TYPE, ABSTRACT, EXTENDS(logging_base) :: physics_base
     !> \name Classes
     CLASS(constants_base), ALLOCATABLE :: constants  !< physical constants
     !> \name
     !! #### Variables
     REAL                :: time,&                !< simulation time
                            mu, &                 !< mean molecular weight
                            eps                   !< softening length
     INTEGER             :: VNUM, &               !< number of variables
                            PNUM, &               !< number of passive variables
                            VDIM, &               !< vector dimensions (1, 2 or 3)
                            DENSITY,PRESSURE, &
                            ENERGY,SGSPRESSURE, &
                            SGSENERGY, &
                            XVELOCITY,XMOMENTUM, &
                            YVELOCITY,YMOMENTUM,&
                            ZVELOCITY,ZMOMENTUM   !< array indicies for primitive and conservative variables
     LOGICAL             :: transformed_xvelocity !< .TRUE. if SubtractBackgroundVelocity was called before
     LOGICAL             :: transformed_yvelocity !< .TRUE. if SubtractBackgroundVelocity was called before
     LOGICAL             :: transformed_zvelocity !< .TRUE. if SubtractBackgroundVelocity was called before
                                                  !! .FALSE. otherwise
     LOGICAL             :: supports_absorbing    !< absorbing boundary conditions supported
                            !! \details .TRUE. if absorbing boundary conditions are supported by the physics module
     LOGICAL             :: supports_farfield     !< farfield boundary conditions supported
                            !! \details .TRUE. if farfield boundary conditions are supported by the physics module
     LOGICAL             :: advanced_wave_speeds  !< use Roe averages for min/max wave speed estimates
     CHARACTER(LEN=16), DIMENSION(:), POINTER &
                         :: pvarname,cvarname     !< names of variables
     !> \name
     !! #### Arrays
     REAL, DIMENSION(:,:,:), POINTER :: &
                            bcradius, &           !< distance to the origin bary center values
                            divposvec, &          !< divergence of the position vector
                            bphi, &               !< bary centered constant gravitational potential
                            tmp,tmp1,tmp2,tmp3, &
                            tmp4,tmp5             !< temporary storage
     REAL, DIMENSION(:,:,:,:), POINTER :: &
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
    PROCEDURE :: PrintConfiguration
    PROCEDURE (new_statevector),              DEFERRED :: new_statevector
    PROCEDURE (ExternalSources),              DEFERRED :: ExternalSources
    PROCEDURE (EnableOutput),                 DEFERRED :: EnableOutput
    !------Convert2Primitve--------!
    PROCEDURE (Convert2Primitive_new),        DEFERRED :: Convert2Primitive_new
    PROCEDURE (Convert2Primitive_center),     DEFERRED :: Convert2Primitive_center
    PROCEDURE (Convert2Primitive_centsub),    DEFERRED :: Convert2Primitive_centsub
    PROCEDURE (Convert2Primitive_faces),      DEFERRED :: Convert2Primitive_faces
    PROCEDURE (Convert2Primitive_facesub),    DEFERRED :: Convert2Primitive_facesub
    GENERIC   :: Convert2Primitive => &
                   Convert2Primitive_new, &
                   Convert2Primitive_center, &
                   Convert2Primitive_centsub, &
                   Convert2Primitive_faces, &
                   Convert2Primitive_facesub
    !------Convert2Conservative----!
    PROCEDURE (Convert2Conservative_new),     DEFERRED :: Convert2Conservative_new
    PROCEDURE (Convert2Conservative_center),  DEFERRED :: Convert2Conservative_center
    PROCEDURE (Convert2Conservative_centsub), DEFERRED :: Convert2Conservative_centsub
    PROCEDURE (Convert2Conservative_faces),   DEFERRED :: Convert2Conservative_faces
    PROCEDURE (Convert2Conservative_facesub), DEFERRED :: Convert2Conservative_facesub
    GENERIC   :: Convert2Conservative => &
                   Convert2Conservative_new, &
                   Convert2Conservative_center, &
                   Convert2Conservative_centsub, &
                   Convert2Conservative_faces, &
                   Convert2Conservative_facesub
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
    !------Fargo Routines-----------!
    PROCEDURE (AddBackgroundVelocityX),       DEFERRED :: AddBackgroundVelocityX
    PROCEDURE (SubtractBackgroundVelocityX),  DEFERRED :: SubtractBackgroundVelocityX
    PROCEDURE (AddBackgroundVelocityY),       DEFERRED :: AddBackgroundVelocityY
    PROCEDURE (SubtractBackgroundVelocityY),  DEFERRED :: SubtractBackgroundVelocityY
    PROCEDURE (AddBackgroundVelocityZ),       DEFERRED :: AddBackgroundVelocityZ
    PROCEDURE (SubtractBackgroundVelocityZ),  DEFERRED :: SubtractBackgroundVelocityZ
    PROCEDURE (AddFargoSources),              DEFERRED :: AddFargoSources
    !------Geometry Routines-------!
    PROCEDURE (GeometricalSources),           DEFERRED :: GeometricalSources
    PROCEDURE (Masks),                        DEFERRED :: ReflectionMasks
    PROCEDURE (Masks),                        DEFERRED :: AxisMasks

! these routines are only necessary for special boundaries, fluxes
   PROCEDURE   (CalculateCharSystemX), DEFERRED :: CalculateCharSystemX
   PROCEDURE   (CalculateCharSystemY), DEFERRED :: CalculateCharSystemY
   PROCEDURE   (CalculateCharSystemZ), DEFERRED :: CalculateCharSystemZ
   PROCEDURE   (CalculateBoundaryDataX), DEFERRED :: CalculateBoundaryDataX
   PROCEDURE   (CalculateBoundaryDataY), DEFERRED :: CalculateBoundaryDataY
   PROCEDURE   (CalculateBoundaryDataZ), DEFERRED :: CalculateBoundaryDataZ
   GENERIC   :: CalcCharSystemX => CalculateCharSystemX
   GENERIC   :: CalcCharSystemY => CalculateCharSystemY
   GENERIC   :: CalcCharSystemZ => CalculateCharSystemZ
   GENERIC   :: CalcBoundaryDataX => CalculateBoundaryDataX
   GENERIC   :: CalcBoundaryDataY => CalculateBoundaryDataY
   GENERIC   :: CalcBoundaryDataZ => CalculateBoundaryDataZ
!   PROCEDURE ::  CalculatePrim2RiemannX
!   PROCEDURE ::  CalculatePrim2RiemannY
!   PROCEDURE ::  CalculateRiemann2PrimX
!   PROCEDURE ::  CalculateRiemann2PrimY
!   PROCEDURE ::  SGSSources
!   PROCEDURE ::  CalculateSGSTensor
!   PROCEDURE :: GetSoundSpeed_adiabatic

    PROCEDURE (Finalize), DEFERRED :: Finalize
    PROCEDURE :: Finalize_base
  END TYPE physics_base

  ABSTRACT INTERFACE
    SUBROUTINE EnableOutput(this,Mesh,config,IO)
      USE common_dict
      IMPORT physics_base, mesh_base
      IMPLICIT NONE
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),        INTENT(IN)   :: Mesh
      TYPE(Dict_TYP), POINTER, INTENT(IN)   :: config, IO
    END SUBROUTINE
    SUBROUTINE new_statevector(this,new_sv,flavour,num)
      IMPORT physics_base, marray_compound
      IMPLICIT NONE
      CLASS(physics_base), INTENT(IN) :: this
      CLASS(marray_compound), POINTER :: new_sv
      INTEGER, OPTIONAL, INTENT(IN)   :: flavour,num
    END SUBROUTINE
    PURE SUBROUTINE Convert2Primitive_new(this,cvar,pvar)
      IMPORT physics_base, marray_compound
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(marray_compound), INTENT(INOUT) :: cvar,pvar
    END SUBROUTINE
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
    PURE SUBROUTINE Convert2Conservative_new(this,pvar,cvar)
      IMPORT physics_base, marray_compound
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
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
    PURE SUBROUTINE ExternalSources(this,accel,pvar,cvar,sterm)
      IMPORT physics_base, mesh_base, marray_base, marray_compound
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(marray_base),   INTENT(IN)      :: accel
      CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    END SUBROUTINE
    PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar,minwav,maxwav)
      IMPORT physics_base,mesh_base,marray_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      CLASS(marray_compound), INTENT(INOUT) :: pvar
      TYPE(marray_base), INTENT(INOUT)     :: minwav,maxwav
    END SUBROUTINE
    PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim,cons,minwav,maxwav)
      IMPORT physics_base,mesh_base,marray_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                           INTENT(IN)    :: prim,cons
      TYPE(marray_base), INTENT(INOUT)   :: minwav,maxwav
    END SUBROUTINE
    PURE SUBROUTINE CalcFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: nmin,nmax
      CLASS(marray_compound), INTENT(INOUT) :: prim,cons,xfluxes
    END SUBROUTINE
    PURE SUBROUTINE CalcFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: nmin,nmax
      CLASS(marray_compound), INTENT(INOUT) :: prim,cons,yfluxes
    END SUBROUTINE
    PURE SUBROUTINE CalcFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(IN)  :: this
      CLASS(mesh_base),    INTENT(IN)  :: Mesh
      INTEGER,             INTENT(IN)  :: nmin,nmax
      CLASS(marray_compound), INTENT(INOUT) :: prim,cons,zfluxes
    END SUBROUTINE
    PURE SUBROUTINE AddBackgroundVelocityX(this,Mesh,w,pvar,cvar)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    END SUBROUTINE AddBackgroundVelocityX
    PURE SUBROUTINE AddBackgroundVelocityY(this,Mesh,w,pvar,cvar)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    END SUBROUTINE AddBackgroundVelocityY
    PURE SUBROUTINE AddBackgroundVelocityZ(this,Mesh,w,pvar,cvar)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    END SUBROUTINE AddBackgroundVelocityZ
    PURE SUBROUTINE SubtractBackgroundVelocityX(this,Mesh,w,pvar,cvar)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    END SUBROUTINE SubtractBackgroundVelocityX
    PURE SUBROUTINE SubtractBackgroundVelocityY(this,Mesh,w,pvar,cvar)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    END SUBROUTINE SubtractBackgroundVelocityY
    PURE SUBROUTINE SubtractBackgroundVelocityZ(this,Mesh,w,pvar,cvar)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    END SUBROUTINE SubtractBackgroundVelocityZ
    PURE SUBROUTINE GeometricalSources(this,Mesh,pvar,cvar,sterm)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    END SUBROUTINE
    PURE SUBROUTINE AddFargoSources(this,Mesh,w,pvar,cvar,sterm)
      IMPORT physics_base,mesh_base,marray_compound
      CLASS(physics_base), INTENT(IN)    :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(IN)    :: w
      CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    END SUBROUTINE AddFargoSources
    PURE SUBROUTINE Masks(this,Mesh,reflX,reflY,reflZ)
      IMPORT physics_base, mesh_base
      CLASS(physics_base),           INTENT(IN)  :: this
      CLASS(mesh_base),              INTENT(IN)  :: Mesh       
      LOGICAL, DIMENSION(this%VNUM), INTENT(OUT) :: reflX,reflY,reflZ
    END SUBROUTINE
    PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i,dir,pvar,lambda,xvar)
      IMPORT physics_base, mesh_base
      !----------------------------------------------------------------------!
      CLASS(physics_base), INTENT(IN) :: this
      CLASS(mesh_base),       INTENT(IN)    :: Mesh
      INTEGER,                INTENT(IN)    :: i,dir
      REAL,                   INTENT(IN), &
        DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                            :: pvar
      REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), INTENT(OUT) :: lambda
      REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), INTENT(OUT)   :: xvar
   END SUBROUTINE
   PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j,dir,pvar,lambda,xvar)
      IMPORT physics_base, mesh_base
      !----------------------------------------------------------------------!
      CLASS(physics_base), INTENT(IN) :: this
      CLASS(mesh_base),       INTENT(IN)    :: Mesh
      INTEGER,                INTENT(IN)    :: j,dir
      REAL,                   INTENT(IN), &
        DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                            :: pvar
      REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), INTENT(OUT) :: lambda
      REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), INTENT(OUT)   :: xvar
   END SUBROUTINE
   PURE SUBROUTINE CalculateCharSystemZ(this,Mesh,k,dir,pvar,lambda,xvar)
      IMPORT physics_base, mesh_base
      !----------------------------------------------------------------------!
      CLASS(physics_base), INTENT(IN) :: this
      CLASS(mesh_base),       INTENT(IN)    :: Mesh
      INTEGER,                INTENT(IN) :: k,dir
      REAL,                   INTENT(IN), &
        DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                            :: pvar
      REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), INTENT(OUT) :: lambda
      REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), INTENT(OUT)   :: xvar
   END SUBROUTINE
   PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i1,dir,xvar,pvar)
      IMPORT physics_base, mesh_base
     CLASS(physics_base), INTENT(IN)    :: this
     CLASS(mesh_base),       INTENT(IN)    :: Mesh
     INTEGER,                INTENT(IN)    :: i1,dir
     REAL,                   INTENT(INOUT), &
        DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                            :: pvar
     REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), INTENT(IN)   :: xvar
   END SUBROUTINE
   PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j1,dir,xvar,pvar)
      IMPORT physics_base, mesh_base
     CLASS(physics_base), INTENT(IN)    :: this
     CLASS(mesh_base),       INTENT(IN)    :: Mesh
     INTEGER,                 INTENT(IN)   :: j1,dir
     REAL,                   INTENT(INOUT), &
        DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                            :: pvar
     REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), INTENT(IN)   :: xvar
   END SUBROUTINE
   PURE SUBROUTINE CalculateBoundaryDataZ(this,Mesh,k1,dir,xvar,pvar)
      IMPORT physics_base, mesh_base
     CLASS(physics_base), INTENT(IN)    :: this
     CLASS(mesh_base),       INTENT(IN)    :: Mesh
     INTEGER,                INTENT(IN)    :: k1,dir
     REAL,                   INTENT(INOUT), &
        DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                            :: pvar
     REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM), INTENT(IN)   :: xvar
   END SUBROUTINE
   SUBROUTINE Finalize(this)
     IMPORT physics_base
     IMPLICIT NONE
     CLASS(physics_base),INTENT(INOUT)  :: this
   END SUBROUTINE

 END INTERFACE
  !--------------------------------------------------------------------------!
  ! flags for advection problems
!  INTEGER, PARAMETER :: EULER2D             = 1
  INTEGER, PARAMETER :: EULER_ISOTHERM      = 16 !> \todo should become 1 in the future,
                                                 !! if all isothermal modules are merged
  INTEGER, PARAMETER :: EULER               = 17 !> \todo should become 2 in the future,
                                                 !! if euler2D/euler3D modules are merged
!  INTEGER, PARAMETER :: EULER2D_ISOTHERM    = 2
!  INTEGER, PARAMETER :: EULER3D_ROTSYM      = 3
!  INTEGER, PARAMETER :: EULER3D_ROTAMT      = 4
!  INTEGER, PARAMETER :: EULER3D_ROTSYMSGS   = 5
!  INTEGER, PARAMETER :: EULER2D_SGS         = 7
!  INTEGER, PARAMETER :: EULER3D_ROTAMTSGS   = 8
!  INTEGER, PARAMETER :: EULER2D_ISOIAMT     = 9
!  INTEGER, PARAMETER :: EULER2D_IAMT        = 11
!  INTEGER, PARAMETER :: EULER2D_IAMROT      = 12
!  INTEGER, PARAMETER :: EULER2D_ISOIAMROT   = 13
!  INTEGER, PARAMETER :: EULER3D_ISOTHERM    = 14
!  INTEGER, PARAMETER :: EULER3D             = 15
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       physics_base, &
       ! constants - flags for  identification in dictionary by an integer
       EULER_ISOTHERM, EULER, &
       SI, CGS, GEOMETRICAL,  &
       UNDEFINED, PRIMITIVE, CONSERVATIVE
  !--------------------------------------------------------------------------!

CONTAINS

  !> Initialization for the base physical object
  !!
  !! - initialize constants and unit system
  !! - read out attributes from config dictrionary (getAttr(..))
  !! - allocation of arrays
  !! - specific tweaks (update soundspeed, etc.)
  !! - print infostring to terminal
  SUBROUTINE InitPhysics(this,Mesh,config,IO,problem,pname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)        :: Mesh
    TYPE(Dict_TYP),POINTER &
                                       :: config, IO
    INTEGER                            :: problem
    CHARACTER(LEN=32)                  :: pname
    !------------------------------------------------------------------------!
    INTEGER                            :: units
    INTEGER                            :: err, valwrite
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: problem
    !------------------------------------------------------------------------!
    CALL this%InitLogging(problem,pname)

    ! check initialization of Mesh
    IF (.NOT.Mesh%Initialized()) &
         CALL this%Error("InitPhysics","mesh module uninitialized")

    ! units
    CALL GetAttr(config, "units", units, SI)
    CALL new_constants(this%constants, units)

    !CALL GetAttr(config, "problem", problem)

    ! mean molecular weight
    CALL GetAttr(config, "mu", this%mu, 0.029)

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

    ! determine physical vector dimensions based on dimimensionality of the grid
    ! whether rotationa symmetry is assumed
    this%VDIM = Mesh%NDIMS
    IF (Mesh%ROTSYM.GT.0) this%VDIM = this%VDIM + 1

    ! set this to appropriate values in derived classes
    this%VNUM = 0  ! number of hydrodynamical variables in state vector
    this%PNUM = 0  ! number of passive scalars in state vector

    ! allocate memory for arrays common to all physics modules
    ALLOCATE(this%tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                 &
             this%tmp1(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp3(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp4(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             this%tmp5(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),                &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitPhysics", "Unable to allocate memory.")

    this%tmp(:,:,:)  = 0.
    this%tmp1(:,:,:) = 0.
    this%tmp2(:,:,:) = 0.
    this%tmp3(:,:,:) = 0.
    this%tmp4(:,:,:) = 0.
    this%tmp5(:,:,:) = 0.

   ! enable/disable absorbing and farfield boundary conditions
   ! TODO Not yet tested for 3D
   SELECT CASE(problem)
!    CASE()
!       this%supports_absorbing = .TRUE.
    CASE DEFAULT
       this%supports_absorbing = .FALSE.
    END SELECT
    SELECT CASE(problem)
!    CASE()
!       this%supports_farfield  = .TRUE.
    CASE DEFAULT
       this%supports_farfield  = .FALSE.
    END SELECT

    ! no background velocity subtracted (important for fargo advection)
    this%transformed_xvelocity = .FALSE.
    this%transformed_yvelocity = .FALSE.
    this%transformed_zvelocity = .FALSE.

    ! reset source term pointer
    !NULLIFY(this%sources)

    this%time = -1.
  END SUBROUTINE InitPhysics

  SUBROUTINE PrintConfiguration(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Info(" PHYSICS--> advection problem: " // TRIM(this%GetName()))
  END SUBROUTINE PrintConfiguration

  !> Destructor
  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("ClosePhysics","not initialized")
    ! deallocate pointer variables used in all physics modules
    DEALLOCATE(this%tmp,this%tmp1,this%tmp2,this%tmp3,this%tmp4,this%tmp5, &
               this%pvarname,this%cvarname)
  END SUBROUTINE Finalize_base

END MODULE physics_base_mod
