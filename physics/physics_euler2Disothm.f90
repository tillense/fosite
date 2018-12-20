!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: physics_euler2Disothm.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2016                                                   #
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
!!
!! \brief basic module for 2D isothermal Euler equations
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler2Dit_mod
  USE physics_base_mod
  USE field_base_mod
!  USE sources_common, ONLY : Sources_TYP
  USE mesh_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: num_var = 3              ! number of variables       !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler 2D isotherm"
  !--------------------------------------------------------------------------!

  TYPE, EXTENDS(FieldSet) :: FieldSet_isoth_pvar
    CLASS(FieldS), ALLOCATABLE :: density
    CLASS(FieldS), ALLOCATABLE :: xvelocity
    CLASS(FieldS), ALLOCATABLE :: yvelocity
  END TYPE

  TYPE, EXTENDS(FieldSet) :: FieldSet_isoth_cvar
    CLASS(FieldS), ALLOCATABLE :: density
    CLASS(FieldS), ALLOCATABLE :: xmomentum
    CLASS(FieldS), ALLOCATABLE :: ymomentum
  END TYPE

  TYPE, EXTENDS(physics_base) :: physics_euler2Dit
    CLASS(FieldSet_isoth_pvar), ALLOCATABLE :: pvar
    CLASS(FieldSet_isoth_cvar), ALLOCATABLE :: cvar
  CONTAINS
    PROCEDURE :: InitPhysics_euler2Dit
    PROCEDURE :: CalcWaveSpeeds
    PROCEDURE :: CalcGeometricalSources
    PROCEDURE :: ReflectionMasks
    PROCEDURE :: AxisMasks
    PROCEDURE :: CalcFlux
    PROCEDURE :: CalcSoundSpeeds_center
    PROCEDURE :: CalcSoundSpeeds_faces
    PROCEDURE :: Cons2Prim
    PROCEDURE :: Prim2Cons
    !PROCEDURE :: ViscositySources_euler2Dit
    !PROCEDURE :: SetEigenValues_euler2Dit
    !PROCEDURE :: SetBoundaryData_euler2Dit
    !PROCEDURE :: SetCharVars_euler2Dit
    !PROCEDURE :: CalcCharSystemX_euler2Dit
    !PROCEDURE :: CalcCharSystemY_euler2Dit
    !PROCEDURE :: CalcBoundaryDataX_euler2Dit
    !PROCEDURE :: CalcBoundaryDataY_euler2Dit
    !PROCEDURE :: CalcPrim2RiemannX_euler2Dit
    !PROCEDURE :: CalcPrim2RiemannY_euler2Dit
    !PROCEDURE :: CalcRiemann2PrimX_euler2Dit
    !PROCEDURE :: CalcRiemann2PrimY_euler2Dit
    !PROCEDURE :: CalcStresses_euler2Dit
    FINAL     :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       physics_euler2Dit
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler2Dit(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN) :: Mesh
    INTEGER           :: problem
    TYPE(Dict_TYP),POINTER &
                      :: config, IO
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    CALL this%InitPhysics(Mesh,config,IO,EULER2D_ISOTHERM,problem_name,num_var)
    !IF (PRESENT(pname).AND.PRESENT(nvar)) THEN
    !   CALL InitPhysics(this,problem,pname,nvar)
    !ELSE IF (PRESENT(pname).OR.PRESENT(nvar)) THEN
    !   CALL Error(this, "InitPhysics_euler2Dit", "Both or no optional " &
    !    // "arguments at all have to be defined.")
    !ELSE
    !   CALL InitPhysics(this,problem,problem_name,num_var)
    !END IF

    ALLOCATE( &
    this%pvar%density%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
    this%pvar%xvelocity%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
    this%pvar%yvelocity%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
    this%cvar%density%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
    this%cvar%xmomentum%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
    this%cvar%ymomentum%data(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
    STAT = err)

    this%pvar = 0.0
    this%cvar = 0.0

    ! set array indices
    this%DENSITY   = 1                                 ! mass density        !
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    this%YVELOCITY = 3                                 ! y-velocity          !
    this%YMOMENTUM = 3                                 ! y-momentum          !
    this%ZVELOCITY = 0                                 ! no z-velocity       !
    this%ZMOMENTUM = 0                                 ! no z-momentum       !
    this%PRESSURE  = 0                                 ! no pressure         !
    this%ENERGY    = 0                                 ! no total energy     !
    ! set names for primitive and conservative variables
    this%pvarname(this%DENSITY)   = "density"
    this%pvarname(this%XVELOCITY) = "xvelocity"
    this%pvarname(this%YVELOCITY) = "yvelocity"
    this%cvarname(this%DENSITY)   = "density"
    this%cvarname(this%XMOMENTUM) = "xmomentum"
    this%cvarname(this%YMOMENTUM) = "ymomentum"
    this%DIM = 2

  END SUBROUTINE InitPhysics_euler2Dit


!  PURE SUBROUTINE CalcCharSystemX_euler2Dit(this,Mesh,i,dir,pvar,lambda,xvar)
!    USE boundary_common, ONLY : WEST
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i,dir
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: lambda,xvar
!    !------------------------------------------------------------------------!
!    INTEGER           :: i1,i2
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i,dir,pvar
!    INTENT(INOUT)     :: lambda
!    INTENT(OUT)       :: xvar
!    !------------------------------------------------------------------------!
!    ! compute eigenvalues at i
!    CALL SetEigenValues_euler2Dit(this%bccsound(i,:),pvar(i,:,this%XVELOCITY), &
!         lambda(:,1),lambda(:,2),lambda(:,3))
!    ! compute characteristic variables using cell mean values of adjacent
!    ! cells to calculate derivatives and the isothermal speed of sound
!    ! at the intermediate cell face
!    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
!    i2 = MAX(i,i1)
!    i1 = MIN(i,i1)
!    CALL SetCharVars_euler2Dit(this%fcsound(i2,:,WEST),pvar(i1,:,this%DENSITY), &
!         pvar(i2,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
!         pvar(i2,:,this%XVELOCITY),pvar(i1,:,this%YVELOCITY), &
!         pvar(i2,:,this%YVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3))
!  END SUBROUTINE CalcCharSystemX_euler2Dit
!
!
!  PURE SUBROUTINE CalcCharSystemY_euler2Dit(this,Mesh,j,dir,pvar,lambda,xvar)
!    USE boundary_common, ONLY : SOUTH
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j,dir
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: lambda,xvar
!    !------------------------------------------------------------------------!
!    INTEGER           :: j1,j2
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j,dir,pvar
!    INTENT(INOUT)     :: lambda
!    INTENT(OUT)       :: xvar
!    !------------------------------------------------------------------------!
!    ! compute eigenvalues at j
!    CALL SetEigenValues_euler2Dit(this%bccsound(:,j),pvar(:,j,this%YVELOCITY), &
!         lambda(:,1),lambda(:,2),lambda(:,3))
!    ! compute characteristic variables using cell mean values of adjacent
!    ! cells to calculate derivatives and the isothermal speed of sound
!    ! at the intermediate cell face
!    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
!    j2 = MAX(j,j1)
!    j1 = MIN(j,j1)
!    CALL SetCharVars_euler2Dit(this%fcsound(:,j2,SOUTH),pvar(:,j1,this%DENSITY), &
!         pvar(:,j2,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
!         pvar(:,j2,this%YVELOCITY),pvar(:,j1,this%XVELOCITY), &
!         pvar(:,j2,this%XVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3))
!  END SUBROUTINE CalcCharSystemY_euler2Dit
!
!
!  PURE SUBROUTINE CalcBoundaryDataX_euler2Dit(this,Mesh,i1,dir,xvar,pvar)
!    USE boundary_common, ONLY : WEST,EAST
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i1,dir
!    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    !------------------------------------------------------------------------!
!    INTEGER           :: i2,fidx
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i1,dir,xvar
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!    i2 = i1 + SIGN(1,dir)  ! i +/- 1 depending on the sign of dir
!    IF (i2.LT.i1) THEN
!       fidx = WEST
!    ELSE
!       fidx = EAST
!    END IF
!    CALL SetBoundaryData_euler2Dit(this%fcsound(i1,:,fidx),1.0*SIGN(1,dir), &
!         pvar(i1,:,this%DENSITY),pvar(i1,:,this%XVELOCITY), &
!         pvar(i1,:,this%YVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3), &
!         pvar(i2,:,this%DENSITY),pvar(i2,:,this%XVELOCITY),pvar(i2,:,this%YVELOCITY))
!  END SUBROUTINE CalcBoundaryDataX_euler2Dit
!
!
!  PURE SUBROUTINE CalcBoundaryDataY_euler2Dit(this,Mesh,j1,dir,xvar,pvar)
!    USE boundary_common, ONLY : SOUTH,NORTH
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j1,dir
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: xvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    !------------------------------------------------------------------------!
!    INTEGER           :: j2,fidx
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j1,dir,xvar
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!    j2 = j1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
!    IF (j2.LT.j1) THEN
!       fidx = SOUTH
!    ELSE
!       fidx = NORTH
!    END IF
!    CALL SetBoundaryData_euler2Dit(this%fcsound(:,j1,fidx),1.0*SIGN(1,dir), &
!         pvar(:,j1,this%DENSITY),pvar(:,j1,this%YVELOCITY), &
!         pvar(:,j1,this%XVELOCITY),xvar(:,1),xvar(:,2),xvar(:,3), &
!         pvar(:,j2,this%DENSITY),pvar(:,j2,this%YVELOCITY),pvar(:,j2,this%XVELOCITY))
!  END SUBROUTINE CalcBoundaryDataY_euler2Dit
!
!
!  PURE SUBROUTINE CalcStresses_euler2Dit(this,Mesh,pvar,dynvis,bulkvis, &
!       btxx,btxy,btyy)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
!         dynvis,bulkvis,btxx,btxy,btyy
!    !------------------------------------------------------------------------!
!    INTEGER           :: i,j
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,pvar,dynvis,bulkvis
!    INTENT(INOUT)     :: this
!    INTENT(OUT)       :: btxx,btxy,btyy
!    !------------------------------------------------------------------------!
!    ! compute components of the stress tensor at cell bary centers
!    ! inside the computational domain including one slice of ghost cells
!
!    ! compute bulk viscosity first and store the result in this%tmp
!!CDIR IEXPAND
!    CALL Divergence(Mesh,pvar(:,:,this%XVELOCITY),pvar(:,:,this%YVELOCITY),this%tmp(:,:))
!    this%tmp(:,:) = bulkvis(:,:)*this%tmp(:,:)
!
!!CDIR OUTERUNROLL=8
!    DO j=Mesh%JMIN-1,Mesh%JMAX+1
!!CDIR NODEP
!       DO i=Mesh%IMIN-1,Mesh%IMAX+1
!          ! compute the diagonal elements of the stress tensor
!          btxx(i,j) = dynvis(i,j) * &
!               ( (pvar(i+1,j,this%XVELOCITY) - pvar(i-1,j,this%XVELOCITY)) / Mesh%dlx(i,j) &
!               + 2.0 * Mesh%cxyx%bcenter(i,j) * pvar(i,j,this%YVELOCITY) ) &
!               + this%tmp(i,j)
!               
!          btyy(i,j) = dynvis(i,j) * &
!               ( (pvar(i,j+1,this%YVELOCITY) - pvar(i,j-1,this%YVELOCITY)) / Mesh%dly(i,j) &
!               + 2.0 * Mesh%cyxy%bcenter(i,j) * pvar(i,j,this%XVELOCITY) ) &
!               + this%tmp(i,j)
!
!          ! compute the off-diagonal elements (no bulk viscosity)
!          btxy(i,j) = dynvis(i,j) * ( 0.5 * &
!               ( (pvar(i+1,j,this%YVELOCITY) - pvar(i-1,j,this%YVELOCITY)) / Mesh%dlx(i,j) &
!               + (pvar(i,j+1,this%XVELOCITY) - pvar(i,j-1,this%XVELOCITY)) / Mesh%dly(i,j) ) &
!               - Mesh%cxyx%bcenter(i,j) * pvar(i,j,this%XVELOCITY) &
!               - Mesh%cyxy%bcenter(i,j) * pvar(i,j,this%YVELOCITY) )
!       END DO
!    END DO
!  END SUBROUTINE CalcStresses_euler2Dit
!
!
!  PURE SUBROUTINE CalcPrim2RiemannX_euler2Dit(this,Mesh,i,pvar,lambda,Rinv)
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
!    ! compute eigenvalues at i
!!CDIR IEXPAND
!    CALL SetEigenValues_euler2Dit(this%bccsound(i,:),pvar(i,:,this%XVELOCITY), &
!         lambda(:,1),lambda(:,2),lambda(:,3))
!    ! compute Riemann invariants
!!CDIR IEXPAND
!    CALL Prim2Riemann_euler2Dit(this%bccsound(i,:),pvar(i,:,this%DENSITY), &
!         pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY), &
!         Rinv(:,1),Rinv(:,2),Rinv(:,3))
!  END SUBROUTINE CalcPrim2RiemannX_euler2Dit
!
!
!  PURE SUBROUTINE CalcPrim2RiemannY_euler2Dit(this,Mesh,j,pvar,lambda,Rinv)
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
!    ! compute eigenvalues at j
!!CDIR IEXPAND
!    CALL SetEigenValues_euler2Dit(this%bccsound(:,j),pvar(:,j,this%YVELOCITY),&
!         lambda(:,1),lambda(:,2),lambda(:,3))
!    ! compute Riemann invariants
!!CDIR IEXPAND
!    CALL Prim2Riemann_euler2Dit(this%bccsound(:,j),pvar(:,j,this%DENSITY), &
!         pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY), &
!         Rinv(:,1),Rinv(:,2),Rinv(:,3))
!  END SUBROUTINE CalcPrim2RiemannY_euler2Dit
!
!  PURE SUBROUTINE CalcRiemann2PrimX_euler2Dit(this,Mesh,i,Rinv,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: i
!    REAL, DIMENSION(Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: Rinv
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    !------------------------------------------------------------------------!
!
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,i,Rinv
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!!CDIR IEXPAND
!    CALL Riemann2Prim_euler2Dit(this%bccsound(i,:),Rinv(:,1),Rinv(:,2),Rinv(:,3), &
!         pvar(i,:,this%DENSITY),pvar(i,:,this%XVELOCITY),pvar(i,:,this%YVELOCITY))
!  END SUBROUTINE CalcRiemann2PrimX_euler2Dit
!
!  PURE SUBROUTINE CalcRiemann2PrimY_euler2Dit(this,Mesh,j,Rinv,pvar)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    INTEGER           :: j
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,this%VNUM) :: Rinv
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: pvar
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,j,Rinv
!    INTENT(INOUT)     :: pvar
!    !------------------------------------------------------------------------!
!!CDIR IEXPAND
!    CALL Riemann2Prim_euler2Dit(this%bccsound(:,j),Rinv(:,1),Rinv(:,2),Rinv(:,3), &
!         pvar(:,j,this%DENSITY),pvar(:,j,this%YVELOCITY),pvar(:,j,this%XVELOCITY))
!  END SUBROUTINE CalcRiemann2PrimY_euler2Dit
!
!
!  PURE SUBROUTINE ViscositySources_euler2Dit(this,Mesh,pvar,btxx,btxy,btyy,sterm)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Physics_TYP) :: this
!    TYPE(Mesh_TYP)    :: Mesh
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: &
!         pvar,sterm
!    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: &
!         btxx,btxy,btyy
!    !------------------------------------------------------------------------!
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: this,Mesh,pvar,btxx,btxy,btyy
!    INTENT(OUT)       :: sterm
!    !------------------------------------------------------------------------!
!    ! mean values of stress tensor components across the cell interfaces
!
!    ! viscosity source terms
!    sterm(:,:,this%DENSITY) = 0.0 
!
!    ! compute viscous momentum sources
!    ! divergence of stress tensor with symmetry btyx=btxy
!!CDIR IEXPAND
!    CALL Divergence(Mesh,btxx,btxy,btxy,btyy,sterm(:,:,this%XMOMENTUM), &
!                    sterm(:,:,this%YMOMENTUM))
!  END SUBROUTINE ViscositySources_euler2Dit



  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .FALSE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .FALSE.
    reflY(this%YVELOCITY) = .TRUE.
  END SUBROUTINE ReflectionMasks


  PURE SUBROUTINE AxisMasks(this,reflX,reflY)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    LOGICAL, DIMENSION(this%VNUM) :: reflX,reflY
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: reflX,reflY
    !------------------------------------------------------------------------!
    ! western / eastern boundary
    reflX(this%DENSITY)   = .FALSE.
    reflX(this%XVELOCITY) = .TRUE.
    reflX(this%YVELOCITY) = .TRUE.
    ! southern / northern boundary
    reflY(this%DENSITY)   = .FALSE.
    reflY(this%XVELOCITY) = .TRUE.
    reflY(this%YVELOCITY) = .TRUE.
  END SUBROUTINE AxisMasks

  PURE SUBROUTINE CalcSoundSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN) :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) &
         :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE CalcSoundSpeeds_center


  PURE SUBROUTINE CalcSoundSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN) :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4,this%VNUM) &
         :: prim
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: prim
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE CalcSoundSpeeds_faces


  ELEMENTAL SUBROUTINE CalcWaveSpeeds(this,cs,v,amin,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: amin,amax
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    amin = MIN(0.,v-cs)
    amax = MAX(0.,v+cs)
  END SUBROUTINE CalcWaveSpeeds


!  ELEMENTAL SUBROUTINE SetEigenValues_euler2Dit(cs,v,l1,l2,l3)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: cs,v
!    REAL, INTENT(OUT) :: l1,l2,l3
!    !------------------------------------------------------------------------!
!    ! all eigenvalues of the isothermal euler problem
!    l1 = v - cs
!    l2 = v
!    l3 = v + cs
!  END SUBROUTINE SetEigenValues_euler2Dit
!
!
!  ELEMENTAL SUBROUTINE SetCharVars_euler2Dit(cs,rho1,rho2,u1,u2,v1,v2,&
!       xvar1,xvar2,xvar3)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: cs,rho1,rho2,u1,u2,v1,v2
!    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3
!    !------------------------------------------------------------------------!
!    REAL :: dlnrho,du
!    !------------------------------------------------------------------------!
!    dlnrho = LOG(rho2/rho1)
!    du = u2-u1
!    ! characteristic variables
!    xvar1 = cs*dlnrho - du
!    xvar2 = v2-v1
!    xvar3 = cs*dlnrho + du
!  END SUBROUTINE SetCharVars_euler2Dit
!
!
!  ELEMENTAL SUBROUTINE SetBoundaryData_euler2Dit(cs,dir,rho1,u1,v1,xvar1,xvar2, &
!       xvar3,rho2,u2,v2)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: cs,dir,rho1,u1,v1,xvar1,xvar2,xvar3
!    REAL, INTENT(OUT) :: rho2,u2,v2
!    !------------------------------------------------------------------------!
!    ! extrapolate boundary values using characteristic variables
!    rho2 = rho1 * EXP(dir*0.5*(xvar3+xvar1)/cs)
!    u2   = u1 + dir*0.5*(xvar3-xvar1)
!    v2   = v1 + dir*xvar2
!  END SUBROUTINE SetBoundaryData_euler2Dit
!
!
!  ELEMENTAL SUBROUTINE Prim2Riemann_euler2Dit(cs,rho,vx,vy,Rminus,Rvt,Rplus)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: cs,rho,vx,vy
!    REAL, INTENT(OUT) :: Rminus,Rvt,Rplus
!    !------------------------------------------------------------------------!
!    REAL :: cslnrho
!    !------------------------------------------------------------------------!
!    cslnrho = cs*LOG(rho)
!    ! compute 1st Riemann invariant (R+)
!    Rplus = vx + cslnrho     
!    ! compute 2st Riemann invariant (R-) 
!    Rminus = vx - cslnrho
!    ! tangential velocity
!    Rvt = vy   
!  END SUBROUTINE Prim2Riemann_euler2Dit
!
!
!  ELEMENTAL SUBROUTINE Riemann2Prim_euler2Dit(cs,Rminus,Rvt,Rplus,rho,vx,vy)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: cs,Rminus,Rvt,Rplus
!    REAL, INTENT(OUT) :: rho,vx,vy
!    !------------------------------------------------------------------------!
!    ! tangential velocity    
!    vy = Rvt  
!    ! normal velocity
!    vx = 0.5*(Rplus+Rminus)
!    ! density
!    rho = EXP(0.5*(Rplus-Rminus)/cs)
!  END SUBROUTINE Riemann2Prim_euler2Dit


  ELEMENTAL SUBROUTINE CalcFlux(this,cs,rho,v,m1,m2,f1,f2,f3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: cs,rho,v,m1,m2
    REAL, INTENT(OUT) :: f1, f2, f3
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + rho*cs*cs
    f3 = m2*v
  END SUBROUTINE CalcFlux


  ! momentum source terms due to inertial forces
  ! P is the isothermal pressure rho*cs*cs
  ELEMENTAL SUBROUTINE CalcGeometricalSources(this,mx,my,vx,vy,P,cxyx,cyxy,czxz,czyz,srho,smx,smy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: mx,my,vx,vy,P,cxyx,cyxy,czxz,czyz
    REAL, INTENT(OUT) :: srho, smx, smy
    !------------------------------------------------------------------------!
    srho = 0.
    smx = -my * (cxyx * vx - cyxy * vy) + (cyxy + czxz) * P
    smy = mx * (cxyx * vx - cyxy * vy) + (cxyx + czyz) * P
  END SUBROUTINE CalcGeometricalSources


  ELEMENTAL SUBROUTINE Cons2Prim(this,rho_in,mu,mv,rho_out,u,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: rho_in,mu,mv
    REAL, INTENT(OUT) :: rho_out,u,v
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u = mu * inv_rho
    v = mv * inv_rho
  END SUBROUTINE Cons2Prim

  
  ELEMENTAL SUBROUTINE Prim2Cons(this,rho_in,u,v,rho_out,mu,mv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler2Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: rho_in,u,v
    REAL, INTENT(OUT) :: rho_out,mu,mv
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
  END SUBROUTINE Prim2Cons
  

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(physics_euler2Dit), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizePhysics()
  END SUBROUTINE Finalize

END MODULE physics_euler2Dit_mod
