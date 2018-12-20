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
!!
!! \brief basic module for 2D isothermal Euler equations
!!
!! \extends physics_common
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_euler3Dit_mod
  USE physics_base_mod
!  USE sources_common, ONLY : Sources_TYP
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
    FINAL     :: Finalize
    PROCEDURE :: CalcWaveSpeeds
    !PROCEDURE :: CalcCharSystemX_euler3Dit
    !PROCEDURE :: CalcCharSystemY_euler3Dit
    !PROCEDURE :: CalcBoundaryDataX_euler3Dit
    !PROCEDURE :: CalcBoundaryDataY_euler3Dit
    !PROCEDURE :: CalcPrim2RiemannX_euler3Dit
    !PROCEDURE :: CalcPrim2RiemannY_euler3Dit
    !PROCEDURE :: CalcRiemann2PrimX_euler3Dit
    !PROCEDURE :: CalcRiemann2PrimY_euler3Dit
    !PROCEDURE :: CalcStresses_euler3Dit
    PROCEDURE :: CalcGeometricalSources
    !PROCEDURE :: ViscositySources_euler3Dit
    PROCEDURE :: ReflectionMasks
    !PROCEDURE :: AxisMasks
    PROCEDURE :: CalcFlux
    PROCEDURE :: CalcSoundSpeeds_center
    PROCEDURE :: CalcSoundSpeeds_faces
    !PROCEDURE :: SetEigenValues_euler3Dit
    !PROCEDURE :: SetBoundaryData_euler3Dit
    !PROCEDURE :: SetCharVars_euler3Dit
    PROCEDURE :: Cons2Prim
    PROCEDURE :: Prim2Cons
  END TYPE
  !--------------------------------------------------------------------------!
   PUBLIC :: &
       ! types
       physics_euler3Dit
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitPhysics_euler3Dit(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN) :: Mesh
    INTEGER           :: problem
    TYPE(Dict_TYP),POINTER &
                      :: config, IO
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    CALL this%InitPhysics(Mesh,config,IO,EULER3D_ISOTH,problem_name,num_var)
    !IF (PRESENT(pname).AND.PRESENT(nvar)) THEN
    !   CALL InitPhysics(this,problem,pname,nvar)
    !ELSE IF (PRESENT(pname).OR.PRESENT(nvar)) THEN
    !   CALL Error(this, "InitPhysics_euler2Dit", "Both or no optional " &
    !    // "arguments at all have to be defined.")
    !ELSE
    !   CALL InitPhysics(this,problem,problem_name,num_var)
    !END IF
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

  PURE SUBROUTINE ReflectionMasks(this,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
    LOGICAL, DIMENSION(this%VNUM)        :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    INTENT(OUT)                          :: reflX,reflY,reflZ
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



  PURE SUBROUTINE CalcSoundSpeeds_center(this,Mesh,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)            :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) &
                                           :: pvar
    !------------------------------------------------------------------------!
    INTEGER                                :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)                             :: pvar
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE CalcSoundSpeeds_center


  PURE SUBROUTINE CalcSoundSpeeds_faces(this,Mesh,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN) :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%nfaces,this%VNUM) &
         :: prim
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,l
    !------------------------------------------------------------------------!
    INTENT(IN)        :: prim
    !------------------------------------------------------------------------!
    ! Sound speed is constant - nothing to do.
  END SUBROUTINE CalcSoundSpeeds_faces

  ELEMENTAL SUBROUTINE CalcWaveSpeeds(this,cs,v,amin,amax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: amin,amax
    !------------------------------------------------------------------------!
    ! minimal and maximal wave speeds
    amin = MIN(0.,v-cs)
    amax = MAX(0.,v+cs)
  END SUBROUTINE CalcWaveSpeeds

  ELEMENTAL SUBROUTINE CalcFlux(this,cs,rho,v,m1,m2,m3,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: cs,rho,v,m1,m2,m3
    REAL, INTENT(OUT) :: f1, f2, f3, f4
    !------------------------------------------------------------------------!
    f1 = rho*v
    f2 = m1*v + rho*cs*cs
    f3 = m2*v
    f4 = m3*v
  END SUBROUTINE CalcFlux


  ! momentum source terms due to inertial forces
  ! P is the isothermal pressure rho*cs*cs
  ELEMENTAL SUBROUTINE CalcGeometricalSources(this,mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz,srho,smx,smy,smz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: mx,my,mz,vx,vy,vz,P,cxyx,cxzx,cyxy,cyzy,czxz,czyz
    REAL, INTENT(OUT) :: srho, smx, smy, smz
    !------------------------------------------------------------------------!
    srho =  0.
    smx  = -my * (cxyx * vx - cyxy * vy) + mz * (czxz * vz - cxzx * vx) + (cyxy + czxz) * P
    smy  =  mx * (cxyx * vx - cyxy * vy) + mz * (czyz * vz - cyzy * vy) + (cxyx + czyz) * P
    smz  =  mx * (cxzx * vx - czxz * vz) + my * (cyzy * vy - czyz * vz) + (cxzx + cyzy) * P
  END SUBROUTINE CalcGeometricalSources


  ELEMENTAL SUBROUTINE Cons2Prim(this,rho_in,mu,mv,mw,rho_out,u,v,w)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
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

  ELEMENTAL SUBROUTINE Prim2Cons(this,rho_in,u,v,w,rho_out,mu,mv,mw)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_euler3Dit), INTENT(IN) :: this
    REAL, INTENT(IN)  :: rho_in,u,v,w
    REAL, INTENT(OUT) :: rho_out,mu,mv,mw
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
    mv = rho_in * v
    mw = rho_in * w
  END SUBROUTINE Prim2Cons
  

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(physics_euler3Dit), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizePhysics()
  END SUBROUTINE Finalize

END MODULE physics_euler3Dit_mod
