!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: physics_eulerisotherm.f90                                         #
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
!! \brief physics module for 1D,2D and 3D isothermal Euler equations
!!
!! \extends physics_base
!! \ingroup physics
!----------------------------------------------------------------------------!
MODULE physics_eulerisotherm_mod
  USE logging_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER           :: num_var = 3          ! number of variables !
  CHARACTER(LEN=32), PARAMETER :: problem_name = "Euler isotherm"
  !--------------------------------------------------------------------------!
  TYPE,  EXTENDS(physics_base) :: physics_eulerisotherm
    !> \name Classes
    !! speed of sound
    CLASS(marray_base), ALLOCATABLE &
                         :: bccsound, &      !< at cell bary centers
                            fcsound          !< at cell faces
    !> \name
    !! #### Variables
    REAL                 :: csiso            !< isothermal sound speed
    LOGICAL :: vector_component_enabled(3)   !< to test which vector components are available
  CONTAINS
    PROCEDURE :: InitPhysics_eulerisotherm
    PROCEDURE :: PrintConfiguration_eulerisotherm
    PROCEDURE :: EnableOutput
    PROCEDURE :: new_statevector
    !------Convert2Primitive-------!
    PROCEDURE :: Convert2Primitive_new
    PROCEDURE :: Convert2Primitive_center
    PROCEDURE :: Convert2Primitive_centsub
    PROCEDURE :: Convert2Primitive_faces
    PROCEDURE :: Convert2Primitive_facesub
    !------Convert2Conservative----!
    PROCEDURE :: Convert2Conservative_new
    PROCEDURE :: Convert2Conservative_center
    PROCEDURE :: Convert2Conservative_centsub
    PROCEDURE :: Convert2Conservative_faces
    PROCEDURE :: Convert2Conservative_facesub
    !------Soundspeed Routines-----!
    PROCEDURE :: SetSoundSpeeds_center
    PROCEDURE :: SetSoundSpeeds_faces
    GENERIC   :: SetSoundSpeeds => SetSoundSpeeds_center, SetSoundSpeeds_faces
    !------Wavespeed Routines------!
    PROCEDURE :: CalcWaveSpeeds_center
    PROCEDURE :: CalcWaveSpeeds_faces
    !------Flux Routines-----------!
    PROCEDURE :: CalcFluxesX
    PROCEDURE :: CalcFluxesY
    PROCEDURE :: CalcFluxesZ
    !------Fargo Routines----------!
    PROCEDURE :: AddBackgroundVelocityX
    PROCEDURE :: AddBackgroundVelocityY
    PROCEDURE :: AddBackgroundVelocityZ
    PROCEDURE :: SubtractBackgroundVelocityX
    PROCEDURE :: SubtractBackgroundVelocityY
    PROCEDURE :: SubtractBackgroundVelocityZ
    PROCEDURE :: AddFargoSources

    PROCEDURE :: GeometricalSources
    PROCEDURE :: ExternalSources
    PROCEDURE :: ViscositySources

    PROCEDURE :: CalcStresses
!    PROCEDURE :: CalcIntermediateStateX_eulerisotherm    ! for HLLC
!    PROCEDURE :: CalcIntermediateStateY_eulerisotherm    ! for HLLC
    PROCEDURE :: CalculateCharSystemX           ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemY           ! for absorbing boundaries
    PROCEDURE :: CalculateCharSystemZ           ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataX         ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataY         ! for absorbing boundaries
    PROCEDURE :: CalculateBoundaryDataZ         ! for absorbing boundaries
!    PROCEDURE :: CalcPrim2RiemannX_eulerisotherm         ! for farfield boundaries
!    PROCEDURE :: CalcPrim2RiemannY_eulerisotherm         ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimX_eulerisotherm         ! for farfield boundaries
!    PROCEDURE :: CalcRiemann2PrimY_eulerisotherm         ! for farfield boundaries
!    PROCEDURE :: CalcRoeAverages_eulerisotherm           ! for advanced wavespeeds
!    PROCEDURE :: ExternalSources_eulerisotherm
    PROCEDURE :: ReflectionMasks                      ! for reflecting boundaries
    PROCEDURE :: AxisMasks

    PROCEDURE     :: Finalize
  END TYPE
  TYPE, EXTENDS(marray_compound) :: statevector_eulerisotherm
    INTEGER :: flavour = UNDEFINED
    TYPE(marray_base), POINTER &
                            :: density => null(), &
                               velocity => null(), &
                               momentum => null()
    CONTAINS
    PROCEDURE :: AssignMArray_0
  END TYPE
!  INTERFACE statevector_eulerisotherm
!    MODULE PROCEDURE CreateStateVector
!  END INTERFACE
  INTERFACE SetFlux
    MODULE PROCEDURE SetFlux1d, SetFlux2d, SetFlux3d
  END INTERFACE
  INTERFACE Cons2Prim
    MODULE PROCEDURE Cons2Prim1d, Cons2Prim2d, Cons2Prim3d
  END INTERFACE
  INTERFACE Prim2Cons
    MODULE PROCEDURE Prim2Cons1d, Prim2Cons2d, Prim2Cons3d
  END INTERFACE
  !--------------------------------------------------------------------------!
   PUBLIC :: &
       ! types
       physics_eulerisotherm, &
       statevector_eulerisotherm, &
       CreateStateVector_eulerisotherm
  !--------------------------------------------------------------------------!

CONTAINS

!----------------------------------------------------------------------------!
!> \par methods of class physics_eulerisotherm

  !> Intialization of isothermal physics
  !!
  !! - calls intialization of base routines of physics
  !! - set array indices, names and number of dimensions
  SUBROUTINE InitPhysics_eulerisotherm(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    TYPE(Dict_TYP), POINTER,  INTENT(IN)    :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: next_idx,err
    !------------------------------------------------------------------------!
    CALL this%InitPhysics(Mesh,config,IO,EULER_ISOTHERM,problem_name)

    ! set the total number of variables in a state vector
    this%VNUM = this%VDIM + 1

    ! get isothermal sound speed from configuration
    CALL GetAttr(config, "cs", this%csiso, 0.0)

    ! allocate memory for arrays used in eulerisotherm
    ALLOCATE(this%pvarname(this%VNUM),this%cvarname(this%VNUM),this%bccsound, &
             this%fcsound, &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitPhysics_eulerisotherm", "Unable to allocate memory.")

    !> \todo remove / improve in future version
    !! set array indices for 1st,2nd,3rd non-vanishing velocities
    !! this may actually not coincide with the x,y and z-velocities
    this%DENSITY   = 1                                 ! mass density        !
    this%pvarname(this%DENSITY)   = "density"
    this%cvarname(this%DENSITY)   = "density"
    this%XVELOCITY = 2                                 ! x-velocity          !
    this%XMOMENTUM = 2                                 ! x-momentum          !
    IF (this%VDIM.GE.2) THEN
      this%YVELOCITY = 3                               ! y-velocity          !
      this%YMOMENTUM = 3                               ! y-momentum          !
    ELSE
      this%YVELOCITY = 0                               ! no y-velocity       !
      this%YMOMENTUM = 0                               ! no y-momentum       !
    END IF
    IF (this%VDIM.EQ.3) THEN
      this%ZVELOCITY = 4                               ! z-velocity          !
      this%ZMOMENTUM = 4                               ! z-momentum          !
    ELSE
      this%ZVELOCITY = 0                               ! no z-velocity       !
      this%ZMOMENTUM = 0                               ! no z-momentum       !
    END IF
    this%PRESSURE  = 0                                 ! no pressure         !
    this%ENERGY    = 0                                 ! no total energy     !

    ! set names shown in the data file and 
    ! whether a particular vector component is available
    this%vector_component_enabled(1:3) = .FALSE.
    next_idx = 2
    IF (Mesh%INUM.GT.1.OR.Mesh%ROTSYM.EQ.1) THEN
      this%pvarname(next_idx) = "xvelocity"
      this%cvarname(next_idx) = "xmomentum"
      this%vector_component_enabled(1) = .TRUE.
      next_idx = next_idx + 1
    END IF
    IF (Mesh%JNUM.GT.1.OR.Mesh%ROTSYM.EQ.2) THEN
      this%pvarname(next_idx) = "yvelocity"
      this%cvarname(next_idx) = "ymomentum"
      this%vector_component_enabled(2) = .TRUE.
      next_idx = next_idx + 1
    END IF
    IF (Mesh%KNUM.GT.1.OR.Mesh%ROTSYM.EQ.3) THEN
      this%pvarname(next_idx) = "zvelocity"
      this%cvarname(next_idx) = "zmomentum"
      this%vector_component_enabled(3) = .TRUE.
    END IF

    ! create new mesh array bccsound
    this%bccsound = marray_base()
    this%fcsound = marray_base(Mesh%NFACES)

    IF(this%csiso.GT.0.) THEN
      this%bccsound%data1d(:) = this%csiso
      this%fcsound%data1d(:)  = this%csiso
    ELSE
      this%bccsound%data1d(:) = 0.
      this%fcsound%data1d(:)  = 0.
    END IF

    CALL this%EnableOutput(Mesh,config,IO)
  END SUBROUTINE InitPhysics_eulerisotherm

  SUBROUTINE PrintConfiguration_eulerisotherm(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%PrintConfiguration()
  END SUBROUTINE PrintConfiguration_eulerisotherm

  !> Enables output of certain arrays defined in this class
  SUBROUTINE EnableOutput(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN)   :: Mesh
    TYPE(Dict_TYP), POINTER, INTENT(IN)   :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: valwrite
    !------------------------------------------------------------------------!
    ! check if output of sound speeds is requested
    CALL GetAttr(config, "output/bccsound", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "bccsound",&
                    this%bccsound%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    CALL GetAttr(config, "output/fcsound", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL Setattr(IO, "fcsound",&
                    this%fcsound%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:))
  END SUBROUTINE EnableOutput

  !> \public allocate an initialize new isothermal state vector
  SUBROUTINE new_statevector(this,new_sv,flavour,num)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN) :: this
    CLASS(marray_compound), POINTER :: new_sv
    INTEGER, OPTIONAL, INTENT(IN)   :: flavour,num
    !------------------------------------------------------------------------!
    ALLOCATE(statevector_eulerisotherm::new_sv)
    SELECT TYPE(sv => new_sv)
    TYPE IS (statevector_eulerisotherm)
!       sv = statevector_eulerisotherm(this,flavour,num)
      CALL CreateStateVector_eulerisotherm(this,sv,flavour,num)
    END SELECT
  END SUBROUTINE new_statevector

  !> Sets soundspeeds at cell-centers
  PURE SUBROUTINE SetSoundSpeeds_center(this,Mesh,bccsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)    :: bccsound
    !------------------------------------------------------------------------!
    this%bccsound = bccsound
  END SUBROUTINE SetSoundSpeeds_center

  !> Sets soundspeeds at cell-faces
  PURE SUBROUTINE SetSoundSpeeds_faces(this,Mesh,fcsound)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES), &
                         INTENT(IN)    :: fcsound
    !------------------------------------------------------------------------!
    this%fcsound = fcsound
  END SUBROUTINE SetSoundSpeeds_faces
  
  !> Converts to primitives at cell centers using state vectors
  PURE SUBROUTINE Convert2Primitive_new(this,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN) :: this
    CLASS(marray_compound), INTENT(INOUT) :: cvar,pvar
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    TYPE IS (statevector_eulerisotherm)
      SELECT TYPE(p => pvar)
      TYPE IS (statevector_eulerisotherm)
        IF (c%flavour.EQ.CONSERVATIVE.AND.p%flavour.EQ.PRIMITIVE) THEN
          ! perform the transformation depending on dimensionality
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D velocity / momentum
            CALL Cons2Prim(c%density%data1d(:),c%momentum%data1d(:), &
                          p%density%data1d(:),p%velocity%data1d(:))
          CASE(2) ! 2D velocity / momentum
            CALL Cons2Prim(c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%momentum%data2d(:,2),p%density%data1d(:), &
                          p%velocity%data2d(:,1),p%velocity%data2d(:,2))
          CASE(3) ! 3D velocity / momentum
            CALL Cons2Prim(c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%momentum%data2d(:,2),c%momentum%data2d(:,3), &
                          p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%velocity%data2d(:,2),p%velocity%data2d(:,3))
          END SELECT
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Primitive_new

  !> Converts to primitives at cell centers
  PURE SUBROUTINE Convert2Primitive_center(this,Mesh,cvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
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
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Cons2Prim(cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                      cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                      pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                      pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY) &
                    )
    CASE(2) ! 2D velocity / momentum
      CALL Cons2Prim(cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                      cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                      cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                      pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                      pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                      pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY) &
                    )
    CASE(3) ! 3D velocity / momentum
      CALL Cons2Prim(cvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                      cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                      cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                      cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM), &
                      pvar(i1:i2,j1:j2,k1:k2,this%DENSITY),   &
                      pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                      pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                      pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY) &
                    )
    END SELECT
  END SUBROUTINE Convert2Primitive_centsub

  !> Converts to conservative variables at faces
  PURE SUBROUTINE Convert2Primitive_faces(this,Mesh,cons,prim)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
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
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: cons
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: prim
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Cons2Prim(cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY) &
                    )
    CASE(2) ! 2D velocity / momentum
      CALL Cons2Prim(cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY) &
                    )
    CASE(3) ! 3D velocity / momentum
      CALL Cons2Prim(cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY) &
                    )
    END SELECT
  END SUBROUTINE Convert2Primitive_facesub

  !> Converts to conservative at cell centers using state vectors
  PURE SUBROUTINE Convert2Conservative_new(this,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN) :: this
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    TYPE IS (statevector_eulerisotherm)
      SELECT TYPE(c => cvar)
      TYPE IS (statevector_eulerisotherm)
        IF (p%flavour.EQ.PRIMITIVE.AND.c%flavour.EQ.CONSERVATIVE) THEN
          ! perform the transformation depending on dimensionality
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D velocity / momentum
            CALL Prim2Cons(p%density%data1d(:),p%velocity%data1d(:), &
                          c%density%data1d(:),c%momentum%data1d(:))
          CASE(2) ! 2D velocity / momentum
            CALL Prim2Cons(p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%velocity%data2d(:,2),c%density%data1d(:), &
                          c%momentum%data2d(:,1),c%momentum%data2d(:,2))
          CASE(3) ! 3D velocity / momentum
            CALL Prim2Cons(p%density%data1d(:),p%velocity%data2d(:,1), &
                          p%velocity%data2d(:,2),p%velocity%data2d(:,3), &
                          c%density%data1d(:),c%momentum%data2d(:,1), &
                          c%momentum%data2d(:,2),c%momentum%data2d(:,3))
          END SELECT
        ELSE
          ! do nothing
        END IF
      END SELECT
    END SELECT
  END SUBROUTINE Convert2Conservative_new

  !> Convert from primtive to conservative variables at cell-centers
  PURE SUBROUTINE Convert2Conservative_center(this,Mesh,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
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
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(IN)  :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM), &
                              INTENT(OUT) :: cvar
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Prim2Cons(pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM) &
                    )
    CASE(2) ! 2D velocity / momentum
      CALL Prim2Cons(pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM) &
                    )
    CASE(3) ! 3D velocity / momentum
      CALL Prim2Cons(pvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     pvar(i1:i2,j1:j2,k1:k2,this%XVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%YVELOCITY), &
                     pvar(i1:i2,j1:j2,k1:k2,this%ZVELOCITY), &
                     cvar(i1:i2,j1:j2,k1:k2,this%DENSITY)  , &
                     cvar(i1:i2,j1:j2,k1:k2,this%XMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%YMOMENTUM), &
                     cvar(i1:i2,j1:j2,k1:k2,this%ZMOMENTUM) &
                    )
    END SELECT
  END SUBROUTINE Convert2Conservative_centsub


  !> Convert to from primitve to conservative variables at cell-faces
  PURE SUBROUTINE Convert2Conservative_faces(this,Mesh,prim,cons)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
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
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: i1,i2,j1,j2,k1,k2
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)  :: prim
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(OUT) :: cons
    !------------------------------------------------------------------------!
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity / momentum
      CALL Prim2Cons(prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM) &
                    )
    CASE(2) ! 2D velocity / momentum
      CALL Prim2Cons(prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM) &
                    )
    CASE(3) ! 3D velocity / momentum
      CALL Prim2Cons(prim(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     prim(i1:i2,j1:j2,k1:k2,:,this%XVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%YVELOCITY), &
                     prim(i1:i2,j1:j2,k1:k2,:,this%ZVELOCITY), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%DENSITY)  , &
                     cons(i1:i2,j1:j2,k1:k2,:,this%XMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%YMOMENTUM), &
                     cons(i1:i2,j1:j2,k1:k2,:,this%ZMOMENTUM) &
                    )
    END SELECT
  END SUBROUTINE Convert2Conservative_facesub

  !> Calculates wave speeds at cell-centers
  PURE SUBROUTINE CalcWaveSpeeds_center(this,Mesh,pvar,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(marray_compound), INTENT(INOUT)   :: pvar
    TYPE(marray_base), INTENT(INOUT)          :: minwav,maxwav
    !------------------------------------------------------------------------!
    INTEGER :: m,n
    !------------------------------------------------------------------------!
    ! compute minimal and maximal wave speeds at cell centers
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_eulerisotherm)
!NEC$ SHORTLOOP
      m = 1
      DO n=1,Mesh%NDIMS
        IF (Mesh%ROTSYM.EQ.n) m = m + 1 ! skip this velocity
        CALL SetWaveSpeeds(this%bccsound%data1d(:),p%velocity%data2d(:,m),&
                minwav%data2d(:,n),maxwav%data2d(:,n))
        m = m + 1
      END DO
    END SELECT
  END SUBROUTINE CalcWaveSpeeds_center

  !> Calculates wave speeds at cell-faces
  PURE SUBROUTINE CalcWaveSpeeds_faces(this,Mesh,prim,cons,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,this%VNUM), &
                              INTENT(IN)    :: prim,cons
    TYPE(marray_base), INTENT(INOUT)          :: minwav,maxwav
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k,m,n
    !------------------------------------------------------------------------!
    m = this%XVELOCITY
    n = 1
    IF (Mesh%INUM.GT.1) THEN
      ! compute minimal and maximal western/eastern wave speeds 
      CALL SetWaveSpeeds(this%fcsound%data4d(:,:,:,2*n-1),prim(:,:,:,2*n-1,m), &
                         this%tmp(:,:,:),this%tmp1(:,:,:))
      CALL SetWaveSpeeds(this%fcsound%data4d(:,:,:,2*n),prim(:,:,:,2*n,m), &
                         minwav%data4d(:,:,:,n),maxwav%data4d(:,:,:,n))
      ! determine the minimum and maximum at cell interfaces
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN+Mesh%IM1,Mesh%IMAX
            ! western & eastern interfaces
            minwav%data4d(i,j,k,n) = MIN(0.0,this%tmp(i+Mesh%IP1,j,k) ,minwav%data4d(i,j,k,n))
            maxwav%data4d(i,j,k,n) = MAX(0.0,this%tmp1(i+Mesh%IP1,j,k),maxwav%data4d(i,j,k,n))
          END DO
        END DO
      END DO
      n = n + 1
      m = m + 1
    ELSE
      IF (Mesh%ROTSYM.EQ.1) m = m + 1  ! increases the velocity index
    END IF

    IF (Mesh%JNUM.GT.1) THEN
      ! compute minimal and maximal southern/northern wave speeds 
      CALL SetWaveSpeeds(this%fcsound%data4d(:,:,:,2*n-1),prim(:,:,:,2*n-1,m), &
                         this%tmp(:,:,:),this%tmp1(:,:,:))
      CALL SetWaveSpeeds(this%fcsound%data4d(:,:,:,2*n),prim(:,:,:,2*n,m), &
                         minwav%data4d(:,:,:,n),maxwav%data4d(:,:,:,n))
      ! determine the minimum and maximum at cell interfaces
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JMIN+Mesh%JM1,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! southern & northern interfaces
            minwav%data4d(i,j,k,n) = MIN(0.0,this%tmp(i,j+Mesh%JP1,k),minwav%data4d(i,j,k,n))
            maxwav%data4d(i,j,k,n) = MAX(0.0,this%tmp1(i,j+Mesh%JP1,k),maxwav%data4d(i,j,k,n))
          END DO
        END DO
      END DO
      n = n + 1
      m = m + 1
    ELSE
      IF (Mesh%ROTSYM.EQ.2) m = m + 1  ! increases the velocity index
    END IF

    IF (Mesh%KNUM.GT.1) THEN
      ! compute minimal and maximal lower/upper wave speeds
      CALL SetWaveSpeeds(this%fcsound%data4d(:,:,:,2*n-1),prim(:,:,:,2*n-1,m), &
                         this%tmp(:,:,:),this%tmp1(:,:,:))
      CALL SetWaveSpeeds(this%fcsound%data4d(:,:,:,2*n),prim(:,:,:,2*n,m), &
                         minwav%data4d(:,:,:,n),maxwav%data4d(:,:,:,n))
      ! determine the minimum and maximum at cell interfaces
      DO k=Mesh%KMIN+Mesh%KM1,Mesh%KMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! bottom & top interfaces
            minwav%data4d(i,j,k,n) = MIN(0.0,this%tmp(i,j,k+Mesh%KP1),minwav%data4d(i,j,k,n))
            maxwav%data4d(i,j,k,n) = MAX(0.0,this%tmp1(i,j,k+Mesh%KP1),maxwav%data4d(i,j,k,n))
          END DO
        END DO
      END DO
    END IF

    !!! THIS IS MOST PROBABLY BROKEN !!!
    ! set minimal and maximal wave speeds at cell interfaces of neighboring cells
!    IF (this%advanced_wave_speeds) THEN
!    DO k=Mesh%KGMIN,Mesh%KGMAX
!      DO j=Mesh%JGMIN,Mesh%JGMAX
!!NEC$ IVDEP
!        DO i=Mesh%IMIN-1,Mesh%IMAX
!          ! western & eastern interfaces
!          ! get Roe averaged x-velocity
!          CALL SetRoeAverages(prim(i,j,k,2,this%DENSITY),prim(i+1,j,k,1,this%DENSITY), &
!                              prim(i,j,k,2,this%XVELOCITY),prim(i+1,j,k,1,this%XVELOCITY), &
!                              uRoe)
!          ! compute Roe averaged wave speeds
!          CALL SetWaveSpeeds(this%fcsound(i,j,k,2),uRoe,aminRoe,amaxRoe)
!          minwav(i,j,k,1) = MIN(aminRoe,this%tmp(i+1,j,k),minwav(i,j,k,1))
!          maxwav(i,j,k,1) = MAX(amaxRoe,this%tmp1(i+1,j,k),maxwav(i,j,k,1))
!        END DO
!      END DO
!    END DO
!    DO k=Mesh%KGMIN,Mesh%KGMAX
!      DO j=Mesh%JMIN-1,Mesh%JMAX
!!NEC$ IVDEP
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          ! southern & northern interfaces
!          ! get Roe averaged y-velocity
!          CALL SetRoeAverages(prim(i,j,k,4,this%DENSITY),prim(i,j+1,k,3,this%DENSITY), &
!                              prim(i,j,k,4,this%YVELOCITY),prim(i,j+1,k,3,this%YVELOCITY), &
!                              uRoe)
!          ! compute Roe averaged wave speeds
!          CALL SetWaveSpeeds(this%fcsound(i,j,k,4),uRoe,aminRoe,amaxRoe)
!          minwav(i,j,k,2) = MIN(aminRoe,this%tmp2(i,j+1,k),minwav(i,j,k,2))
!          maxwav(i,j,k,2) = MAX(amaxRoe,this%tmp3(i,j+1,k),maxwav(i,j,k,2))
!        END DO
!      END DO
!    END DO
!    DO k=Mesh%KGMIN-1,Mesh%KGMAX
!      DO j=Mesh%JMIN,Mesh%JMAX
!!NEC$ IVDEP
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          ! topper & bottomer interfaces
!          ! get Roe averaged z-velocity
!          CALL SetRoeAverages(prim(i,j,k,6,this%DENSITY),prim(i,j,k+1,5,this%DENSITY), &
!                               prim(i,j,k,6,this%ZVELOCITY),prim(i,j,k+1,5,this%ZVELOCITY), &
!                               uRoe)
!          ! compute Roe averaged wave speeds
!          CALL SetWaveSpeeds(this%fcsound(i,j,k,6),uRoe,aminRoe,amaxRoe)
!          minwav(i,j,k,3) = MIN(aminRoe,this%tmp3(i,j,k+Mesh%kp1),minwav(i,j,k,2))
!          maxwav(i,j,k,3) = MAX(amaxRoe,this%tmp4(i,j,k+Mesh%kp1),maxwav(i,j,k,2))
!        END DO
!      END DO
!    END DO
!    ELSE
!    END IF
  END SUBROUTINE CalcWaveSpeeds_faces

  !> Calculates geometrical sources
  PURE SUBROUTINE GeometricalSources(this,Mesh,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    ! compute geometrical source only for non-cartesian mesh
    IF (Mesh%Geometry%GetType().NE.CARTESIAN) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          SELECT TYPE(s => sterm)
          TYPE IS(statevector_eulerisotherm)
            ! no source terms
            s%density%data1d(:) = 0.0
            SELECT CASE(Mesh%VECTOR_COMPONENTS)
            CASE(VECTOR_X) ! 1D momentum in x-direction
              ! vy = vz = my = mz = 0
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),0.0,0.0, &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  0.0,0.0)
            CASE(VECTOR_Y) ! 1D momentum in y-direction
              ! vx = vz = mx = mz = 0
              s%momentum%data2d(:,1) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,p%velocity%data2d(:,1),0.0, &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  0.0,0.0)
            CASE(VECTOR_Z) ! 1D momentum in z-direction
              ! vx = vy = mx = my = 0
              s%momentum%data2d(:,1) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,0.0,p%velocity%data2d(:,1), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  0.0,0.0)
            CASE(IOR(VECTOR_X,VECTOR_Y)) ! 2D momentum in x-y-plane
              ! vz = mz = 0
              ! x-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2),0.0, &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  c%momentum%data2d(:,2),0.0)
              ! y-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2),0.0, &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  c%momentum%data2d(:,1),0.0)
            CASE(IOR(VECTOR_X,VECTOR_Z)) ! 2D momentum in x-z-plane
              ! vy = my = 0
              ! x-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),0.0,p%velocity%data2d(:,2), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  0.0,c%momentum%data2d(:,2))
              ! z-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),0.0,p%velocity%data2d(:,2), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  c%momentum%data2d(:,1),0.0)
            CASE(IOR(VECTOR_Y,VECTOR_Z)) ! 2D momentum in y-z-plane
              ! vx = mx = 0
              ! y-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  0.0,c%momentum%data2d(:,2))
              ! z-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  0.0,p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  0.0,c%momentum%data2d(:,1))
            CASE(IOR(IOR(VECTOR_X,VECTOR_Y),VECTOR_Z)) ! 3D momentum
              ! x-momentum
              s%momentum%data2d(:,1) = GetGeometricalSourceX( &
                  Mesh%cxyx%data2d(:,2),Mesh%cxzx%data2d(:,2), &
                  Mesh%cyxy%data2d(:,2),Mesh%czxz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%velocity%data2d(:,3), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  c%momentum%data2d(:,2),c%momentum%data2d(:,3))
              ! y-momentum
              s%momentum%data2d(:,2) = GetGeometricalSourceY( &
                  Mesh%cxyx%data2d(:,2),Mesh%cyxy%data2d(:,2), &
                  Mesh%cyzy%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%velocity%data2d(:,3), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  c%momentum%data2d(:,1),c%momentum%data2d(:,3))
              ! z-momentum
              s%momentum%data2d(:,3) = GetGeometricalSourceZ( &
                  Mesh%cxzx%data2d(:,2),Mesh%cyzy%data2d(:,2), &
                  Mesh%czxz%data2d(:,2),Mesh%czyz%data2d(:,2), &
                  p%velocity%data2d(:,1),p%velocity%data2d(:,2), &
                  p%velocity%data2d(:,3), &
                  p%density%data1d(:)*this%bccsound%data1d(:)**2, &
                  c%momentum%data2d(:,1),c%momentum%data2d(:,2))
            CASE DEFAULT
              ! return NaN
              s%momentum%data1d(:) = NAN_DEFAULT_REAL
            END SELECT
          END SELECT
        END SELECT
      END SELECT
      ! reset ghost cell data
      ! reset ghost cell data
      IF (Mesh%INUM.GT.1) THEN
        sterm%data4d(Mesh%IGMIN:Mesh%IMIN+Mesh%IM1,:,:,:) = 0.0
        sterm%data4d(Mesh%IMAX+Mesh%IP1:Mesh%IGMAX,:,:,:) = 0.0
      END IF
      IF (Mesh%JNUM.GT.1) THEN
        sterm%data4d(:,Mesh%JGMIN:Mesh%JMIN+Mesh%JM1,:,:) = 0.0
        sterm%data4d(:,Mesh%JMAX+Mesh%JP1:Mesh%JGMAX,:,:) = 0.0
      END IF
      IF (Mesh%KNUM.GT.1) THEN
        ! collapse the first 2 dimensions
        sterm%data3d(:,Mesh%KGMIN:Mesh%KMIN+Mesh%KM1,:) = 0.0
        sterm%data3d(:,Mesh%KMAX+Mesh%KP1:Mesh%KGMAX,:) = 0.0
      END IF
    END IF
  END SUBROUTINE GeometricalSources

  !> Calculate Fluxes in x-direction
  PURE SUBROUTINE CalcFluxesX(this,Mesh,nmin,nmax,prim,cons,xfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: nmin,nmax
    CLASS(marray_compound), INTENT(INOUT) :: prim,cons,xfluxes
    !------------------------------------------------------------------------!
    SELECT TYPE(p => prim)
    TYPE IS(statevector_eulerisotherm)
      SELECT TYPE(c => cons)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(f => xfluxes)
        TYPE IS(statevector_eulerisotherm)
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1))
          CASE(2) ! 2D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,2))
          CASE(3) ! 3D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,3),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,3))
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcFluxesX

  !> Calculate Fluxes in y-direction
  PURE SUBROUTINE CalcFluxesY(this,Mesh,nmin,nmax,prim,cons,yfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: nmin,nmax
    CLASS(marray_compound), INTENT(INOUT) :: prim,cons,yfluxes
    !------------------------------------------------------------------------!
    SELECT TYPE(p => prim)
    TYPE IS(statevector_eulerisotherm)
      SELECT TYPE(c => cons)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(f => yfluxes)
        TYPE IS(statevector_eulerisotherm)
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1))
          CASE(2) ! 2D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,1))
          CASE(3) ! 3D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,3),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,3))
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcFluxesY

  !> Calculate Fluxes in z-direction
  PURE SUBROUTINE CalcFluxesZ(this,Mesh,nmin,nmax,prim,cons,zfluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(mesh_base),         INTENT(IN)  :: Mesh
    INTEGER,                  INTENT(IN)  :: nmin,nmax
    CLASS(marray_compound), INTENT(INOUT) :: prim,cons,zfluxes
    !------------------------------------------------------------------------!
    SELECT TYPE(p => prim)
    TYPE IS(statevector_eulerisotherm)
      SELECT TYPE(c => cons)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(f => zfluxes)
        TYPE IS(statevector_eulerisotherm)
          SELECT CASE(this%VDIM)
          CASE(1) ! 1D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,1))
          CASE(2) ! 2D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,2),    &
                      f%momentum%data3d(:,nmin:nmax,1))
          CASE(3) ! 3D flux
            CALL SetFlux(this%fcsound%data2d(:,nmin:nmax), &
                      p%density%data2d(:,nmin:nmax),       &
                      p%velocity%data3d(:,nmin:nmax,3),    &
                      c%momentum%data3d(:,nmin:nmax,3),    &
                      c%momentum%data3d(:,nmin:nmax,1),    &
                      c%momentum%data3d(:,nmin:nmax,2),    &
                      f%density%data2d(:,nmin:nmax),       &
                      f%momentum%data3d(:,nmin:nmax,3),    &
                      f%momentum%data3d(:,nmin:nmax,1),    &
                      f%momentum%data3d(:,nmin:nmax,2))
          END SELECT
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE CalcFluxesZ

  !> compute viscous source terms
  PURE SUBROUTINE ViscositySources(this,Mesh,pvar,btxx,btxy,btxz,btyy,btyz,btzz,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(marray_compound), INTENT(INOUT) :: pvar,sterm
    REAL,                   INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                                           :: btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_eulerisotherm)
      SELECT TYPE(s => sterm)
      CLASS IS(statevector_eulerisotherm)
        ! no viscous sources in continuity equation
        s%density%data1d(:) = 0.0
        ! viscous momentum sources
        SELECT CASE(this%VDIM)
!         CASE(1) ! 1D velocities are currently not supported
!           CALL Mesh%Divergence(btxx,sterm(:,:,:,this%XMOMENTUM))
        CASE(2)
          ! divergence of stress tensor with symmetry btyx=btxy
          CALL Mesh%Divergence(btxx,btxy,btxy,btyy,s%momentum%data4d(:,:,:,1), &
                               s%momentum%data4d(:,:,:,2))
        CASE(3)
          ! divergence of stress tensor with symmetry btyx=btxy, btxz=btzx, btyz=btzy
          CALL Mesh%Divergence(btxx,btxy,btxz,btxy,btyy,btyz,btxz,btyz,btzz, &
               s%momentum%data4d(:,:,:,1),s%momentum%data4d(:,:,:,2), &
               s%momentum%data4d(:,:,:,3))
        CASE DEFAULT
          ! return NaN
          s%data1d(:) = NAN_DEFAULT_REAL
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE ViscositySources

  !> calculate components of the stress tensor
  !!
  !! The components are computed at cell bary centers inside the computational
  !! domain including one slice of ghost cells.
  PURE SUBROUTINE CalcStresses(this,Mesh,pvar,dynvis,bulkvis, &
       btxx,btxy,btxz,btyy,btyz,btzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base), INTENT(IN)          :: Mesh
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    CLASS(marray_base), INTENT(INOUT)     :: dynvis,bulkvis
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: &
         btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: btxx,btxy,btxz,btyy,btyz,btzz
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_eulerisotherm)
      SELECT CASE(Mesh%VECTOR_COMPONENTS)
!!!! 1D velocities are currently not supported !!!!
!         CASE(VECTOR_X) ! 1D velocity in x-direction
!         CASE(VECTOR_Y) ! 1D velocity in y-direction
!         CASE(VECTOR_Z) ! 1D velocity in z-direction
      CASE(IOR(VECTOR_X,VECTOR_Y)) ! 2D velocities in x-y-plane
        ! compute bulk viscosity first and store the result in this%tmp
        CALL Mesh%Divergence(p%velocity%data4d(:,:,:,1),p%velocity%data4d(:,:,:,2),&
                             this%tmp(:,:,:))
        this%tmp(:,:,:) = bulkvis%data3d(:,:,:)*this%tmp(:,:,:)
        DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
          DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!NEC$ IVDEP
            DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
              ! compute the diagonal elements of the stress tensor
              btxx(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,1) - p%velocity%data4d(i-Mesh%IP1,j,k,1)) &
                      / Mesh%dlx%data3d(i,j,k) &
                  + 2.0 * Mesh%cxyx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2)) + this%tmp(i,j,k)
              btyy(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i,j+Mesh%JP1,k,2) - p%velocity%data4d(i,j-Mesh%JP1,k,2)) &
                      / Mesh%dly%data3d(i,j,k) &
                  + 2.0 * Mesh%cyxy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) ) + this%tmp(i,j,k)
              ! compute the off-diagonal elements (no bulk viscosity)
              btxy(i,j,k) = dynvis%data3d(i,j,k) * ( 0.5 * &
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,2) - p%velocity%data4d(i-Mesh%IP1,j,k,2)) &
                      / Mesh%dlx%data3d(i,j,k) &
                  +  (p%velocity%data4d(i,j+Mesh%JP1,k,1) - p%velocity%data4d(i,j-Mesh%JP1,k,1)) &
                      / Mesh%dly%data3d(i,j,k) ) &
                  - Mesh%cxyx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) &
                  - Mesh%cyxy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) )
            END DO
          END DO
        END DO
      CASE(IOR(VECTOR_X,VECTOR_Z)) ! 2D velocities in x-z-plane
        ! compute bulk viscosity first and store the result in this%tmp
        CALL Mesh%Divergence(p%velocity%data4d(:,:,:,1),p%velocity%data4d(:,:,:,2),&
                             this%tmp(:,:,:))
        this%tmp(:,:,:) = bulkvis%data3d(:,:,:)*this%tmp(:,:,:)
        DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
          DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!NEC$ IVDEP
            DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
              ! compute the diagonal elements btxx, btzz => btyy of the stress tensor
              btxx(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,1) - p%velocity%data4d(i-Mesh%IP1,j,k,1)) &
                      / Mesh%dlx%data3d(i,j,k) &
                  + 2.0 * Mesh%cxzx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) ) & ! vz => vy
                  + this%tmp(i,j,k)
              btyy(i,j,k) = dynvis%data3d(i,j,k) * & ! btzz => btyy
                    ((p%velocity%data4d(i,j,k+Mesh%KP1,2) - p%velocity%data4d(i,j,k-Mesh%KP1,2)) & ! vz => vy
                      / Mesh%dlz%data3d(i,j,k) &
                  + 2.0 * Mesh%czxz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) ) &
                  + this%tmp(i,j,k)
              ! compute the off-diagonal elements btxz => btxy (no bulk viscosity)
              btxy(i,j,k) = dynvis%data3d(i,j,k) * ( 0.5 * & ! btxz => btxy
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,2) - p%velocity%data4d(i-Mesh%IP1,j,k,2)) & ! vz => vy
                      / Mesh%dlx%data3d(i,j,k) &
                  +  (p%velocity%data4d(i,j,k+Mesh%KP1,1) - p%velocity%data4d(i,j,k-Mesh%KP1,1)) &
                      / Mesh%dlz%data3d(i,j,k) ) &
                  - Mesh%czxz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) & ! vz => vy
                  - Mesh%cxzx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) )
            END DO
          END DO
        END DO
      CASE(IOR(VECTOR_Y,VECTOR_Z)) ! 2D velocities in y-z-plane
        CALL Mesh%Divergence(p%velocity%data4d(:,:,:,1),p%velocity%data4d(:,:,:,2),&
                             this%tmp(:,:,:))
        this%tmp(:,:,:) = bulkvis%data3d(:,:,:)*this%tmp(:,:,:)
        DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
          DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!NEC$ IVDEP
            DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
              ! compute the diagonal elements btyy => btxx, btzz => btyy of the stress tensor
              btxx(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i,j+Mesh%JP1,k,1) - p%velocity%data4d(i,j-Mesh%JP1,k,1)) &
                      / Mesh%dly%data3d(i,j,k) &
                  + 2.0 * Mesh%cyzy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) ) & ! vz => vy
                  + this%tmp(i,j,k)
              btyy(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i,j,k+Mesh%KP1,2) - p%velocity%data4d(i,j,k-Mesh%KP1,2)) & ! vz => vy
                      / Mesh%dlz%data3d(i,j,k) &
                  + 2.0 * Mesh%czyz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) ) &
                  + this%tmp(i,j,k)
              ! compute the off-diagonal elements btyz => btxy (no bulk viscosity)
              btxy(i,j,k) = dynvis%data3d(i,j,k) * ( 0.5 * &
                    ((p%velocity%data4d(i,j,k+Mesh%KP1,1) - p%velocity%data4d(i,j,k-Mesh%KP1,1)) & ! vy => vx
                      / Mesh%dlz%data3d(i,j,k) &
                  +  (p%velocity%data4d(i,j+Mesh%JP1,k,2) - p%velocity%data4d(i,j-Mesh%JP1,k,2)) & ! vz => vy
                      / Mesh%dly%data3d(i,j,k) ) &
                  - Mesh%czyz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) & ! vz => vy
                  - Mesh%cyzy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) ) ! vy => vx
            END DO
          END DO
        END DO
      CASE(IOR(IOR(VECTOR_X,VECTOR_Y),VECTOR_Z)) ! 3D velocities
        ! compute bulk viscosity first and store the result in this%tmp
        CALL Mesh%Divergence(p%velocity%data4d(:,:,:,1),p%velocity%data4d(:,:,:,2),&
                            p%velocity%data4d(:,:,:,3),this%tmp(:,:,:))
        this%tmp(:,:,:) = bulkvis%data3d(:,:,:)*this%tmp(:,:,:)
        DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
          DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX+Mesh%JP1
!NEC$ IVDEP
            DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX+Mesh%IP1
              ! compute the diagonal elements of the stress tensor
              btxx(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,1) - p%velocity%data4d(i-Mesh%IP1,j,k,1)) &
                      / Mesh%dlx%data3d(i,j,k) &
                  + 2.0 * Mesh%cxyx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) &
                  + 2.0 * Mesh%cxzx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,3) ) &
                  + this%tmp(i,j,k)
              btyy(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i,j+Mesh%JP1,k,2) - p%velocity%data4d(i,j-Mesh%JP1,k,2)) &
                      / Mesh%dly%data3d(i,j,k) &
                  + 2.0 * Mesh%cyxy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1)  &
                  + 2.0 * Mesh%cyzy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,3) ) &
                  + this%tmp(i,j,k)
              btzz(i,j,k) = dynvis%data3d(i,j,k) * &
                    ((p%velocity%data4d(i,j,k+Mesh%KP1,3) - p%velocity%data4d(i,j,k-Mesh%KP1,3)) &
                      / Mesh%dlz%data3d(i,j,k) &
                  + 2.0 * Mesh%czxz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) &
                  + 2.0 * Mesh%czyz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) ) &
                  + this%tmp(i,j,k)
              ! compute the off-diagonal elements (no bulk viscosity)
              btxy(i,j,k) = dynvis%data3d(i,j,k) * ( 0.5 * &
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,2) - p%velocity%data4d(i-Mesh%IP1,j,k,2)) &
                      / Mesh%dlx%data3d(i,j,k) &
                  +  (p%velocity%data4d(i,j+Mesh%JP1,k,1) - p%velocity%data4d(i,j-Mesh%JP1,k,1)) &
                      / Mesh%dly%data3d(i,j,k) ) &
                  - Mesh%cxyx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) &
                  - Mesh%cyxy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) )
              btxz(i,j,k) = dynvis%data3d(i,j,k) * ( 0.5 * &
                    ((p%velocity%data4d(i+Mesh%IP1,j,k,3) - p%velocity%data4d(i-Mesh%IP1,j,k,3)) &
                      / Mesh%dlx%data3d(i,j,k) &
                  +  (p%velocity%data4d(i,j,k+Mesh%KP1,1) - p%velocity%data4d(i,j,k-Mesh%KP1,1)) &
                      / Mesh%dlz%data3d(i,j,k) ) &
                  - Mesh%czxz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,3) &
                  - Mesh%cxzx%bcenter(i,j,k) * p%velocity%data4d(i,j,k,1) )
              btyz(i,j,k) = dynvis%data3d(i,j,k) * ( 0.5 * &
                    ((p%velocity%data4d(i,j,k+Mesh%KP1,2) - p%velocity%data4d(i,j,k-Mesh%KP1,2)) &
                      / Mesh%dlz%data3d(i,j,k) &
                  +  (p%velocity%data4d(i,j+Mesh%JP1,k,3) - p%velocity%data4d(i,j-Mesh%JP1,k,3)) &
                      / Mesh%dly%data3d(i,j,k) ) &
                  - Mesh%czyz%bcenter(i,j,k) * p%velocity%data4d(i,j,k,3) &
                  - Mesh%cyzy%bcenter(i,j,k) * p%velocity%data4d(i,j,k,2) )
            END DO
          END DO
        END DO
      CASE DEFAULT
        ! return NaN
        btxx(:,:,:) = NAN_DEFAULT_REAL
        btxy(:,:,:) = NAN_DEFAULT_REAL
        btxz(:,:,:) = NAN_DEFAULT_REAL
        btyy(:,:,:) = NAN_DEFAULT_REAL
        btyz(:,:,:) = NAN_DEFAULT_REAL
        btzz(:,:,:) = NAN_DEFAULT_REAL
      END SELECT
    END SELECT
  END SUBROUTINE CalcStresses

  !> compute momentum sources given an external force
  PURE SUBROUTINE ExternalSources(this,accel,pvar,cvar,sterm)
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN)  :: this
    CLASS(marray_base),       INTENT(IN)  :: accel
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER :: n
    !------------------------------------------------------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_eulerisotherm)
      SELECT TYPE(c => cvar)
      CLASS IS(statevector_eulerisotherm)
        SELECT TYPE(s => sterm)
        CLASS IS(statevector_eulerisotherm)
          s%density%data1d(:) = 0.0
!NEC$ UNROLL(3)
          DO n=1,this%VDIM
            s%momentum%data2d(:,n) = c%density%data1d(:) * accel%data2d(:,n)
          END DO
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE ExternalSources


  PURE SUBROUTINE AddBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_xvelocity) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                p%velocity%data4d(i,j,k,1) = p%velocity%data4d(i,j,k,1) + w(j,k)
                c%momentum%data4d(i,j,k,1) = c%momentum%data4d(i,j,k,1) &
                    + c%density%data3d(i,j,k)*w(j,k)
              END DO
            END DO
          END DO
          this%transformed_xvelocity = .FALSE.
        END SELECT
      END SELECT
      this%transformed_xvelocity = .FALSE.
    END IF
  END SUBROUTINE AddBackgroundVelocityX

  PURE SUBROUTINE AddBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_yvelocity) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                p%velocity%data4d(i,j,k,2) = p%velocity%data4d(i,j,k,2) + w(i,k)
                c%momentum%data4d(i,j,k,2) = c%momentum%data4d(i,j,k,2) &
                    + c%density%data3d(i,j,k)*w(i,k)
              END DO
            END DO
          END DO
          this%transformed_yvelocity = .FALSE.
        END SELECT
      END SELECT
    END IF
  END SUBROUTINE AddBackgroundVelocityY

  PURE SUBROUTINE AddBackgroundVelocityZ(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                              INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (this%transformed_zvelocity) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                p%velocity%data4d(i,j,k,3) = p%velocity%data4d(i,j,k,3) + w(i,j)
                c%momentum%data4d(i,j,k,3) = c%momentum%data4d(i,j,k,3) &
                    + c%density%data3d(i,j,k)*w(i,j)
              END DO
            END DO
          END DO
          this%transformed_zvelocity = .FALSE.
        END SELECT
      END SELECT
    END IF
  END SUBROUTINE AddBackgroundVelocityZ

  PURE SUBROUTINE SubtractBackgroundVelocityX(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_xvelocity) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                p%velocity%data4d(i,j,k,1) = p%velocity%data4d(i,j,k,1) - w(j,k)
                c%momentum%data4d(i,j,k,1) = c%momentum%data4d(i,j,k,1) &
                    - c%density%data3d(i,j,k)*w(j,k)
              END DO
            END DO
          END DO
          this%transformed_xvelocity = .TRUE.
        END SELECT
      END SELECT
    END IF
  END SUBROUTINE SubtractBackgroundVelocityX

  PURE SUBROUTINE SubtractBackgroundVelocityY(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_yvelocity) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                p%velocity%data4d(i,j,k,2) = p%velocity%data4d(i,j,k,2) - w(i,k)
                c%momentum%data4d(i,j,k,2) = c%momentum%data4d(i,j,k,2) &
                    - c%density%data3d(i,j,k)*w(i,k)
              END DO
            END DO
          END DO
          this%transformed_yvelocity = .TRUE.
        END SELECT
      END SELECT
    END IF
  END SUBROUTINE SubtractBackgroundVelocityY

  PURE SUBROUTINE SubtractBackgroundVelocityZ(this,Mesh,w,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    !------------------------------------------------------------------------!
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                              INTENT(IN)    :: w
    CLASS(marray_compound), INTENT(INOUT) ::  pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    IF (.NOT.this%transformed_zvelocity) THEN
      SELECT TYPE(p => pvar)
      TYPE IS(statevector_eulerisotherm)
        SELECT TYPE(c => cvar)
        TYPE IS(statevector_eulerisotherm)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                p%velocity%data4d(i,j,k,3) = p%velocity%data4d(i,j,k,3) - w(i,j)
                c%momentum%data4d(i,j,k,3) = c%momentum%data4d(i,j,k,3) &
                    - c%density%data3d(i,j,k)*w(i,j)
              END DO
            END DO
          END DO
          this%transformed_zvelocity = .TRUE.
        END SELECT
      END SELECT
    END IF
  END SUBROUTINE SubtractBackgroundVelocityZ



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
  PURE SUBROUTINE AddFargoSources(this,Mesh,w,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN) :: this
    CLASS(mesh_base), INTENT(IN)             :: Mesh
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), INTENT(IN) &
                                             :: w
    CLASS(marray_compound), INTENT(INOUT)    :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(c => cvar)
    CLASS IS(statevector_eulerisotherm)
      SELECT TYPE(s => sterm)
      CLASS IS(statevector_eulerisotherm)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              ! ATTENTION: fargo sources are added to the given sterm
              s%momentum%data4d(i,j,k,2) = s%momentum%data4d(i,j,k,2) &
                - c%momentum%data4d(i,j,k,1) * 0.5 * (w(i+1,k)-w(i-1,k)) &
                / Mesh%dlx%data3d(i,j,k)
            END DO
          END DO
        END DO
      END SELECT
    END SELECT
  END SUBROUTINE AddFargoSources

  !> return masks for reflecting boundaries
  !!
  !! At axis boundaries we change the sign of normal velocities at each boundary.
  PURE SUBROUTINE ReflectionMasks(this,Mesh,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm),  INTENT(IN)  :: this
    CLASS(mesh_base),              INTENT(IN)  :: Mesh       
    LOGICAL, DIMENSION(this%VNUM), INTENT(OUT) :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    reflX(:) = .FALSE.
    reflY(:) = .FALSE.
    reflZ(:) = .FALSE.
    SELECT CASE(this%VDIM)
    CASE(1) ! 1D velocity
      IF (Mesh%INUM.GT.1) reflX(2) = .TRUE. ! reflect vx at east/west-boundaries
      IF (Mesh%JNUM.GT.1) reflY(2) = .TRUE. ! reflect vy at south/north-boundaries
      IF (Mesh%KNUM.GT.1) reflZ(2) = .TRUE. ! reflect vz at bottom/top-boundaries
    CASE(2) ! 2D velocity
      IF (Mesh%KNUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.3) THEN
        ! transport in x-y-plane
        reflX(2) = .TRUE. ! reflect vx at east/west-boundaries
        reflY(3) = .TRUE. ! reflect vy at south/north-boundaries
      ELSE IF (Mesh%JNUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.2) THEN
        ! transport in x-z-plane
        reflX(2) = .TRUE. ! reflect vx at east/west-boundaries
        reflZ(3) = .TRUE. ! reflect vz at bottom/top-boundaries
      ELSE IF (Mesh%INUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.1) THEN
        ! transport in y-z-plane
        reflY(2) = .TRUE. ! reflect vy at south/north-boundaries
        reflZ(3) = .TRUE. ! reflect vz at bottom/top-boundaries
      END IF
    CASE(3) ! 3D velocity
      reflX(2) = .TRUE. ! reflect vx at east/west-boundaries
      reflY(3) = .TRUE. ! reflect vy at south/north-boundaries
      reflZ(4) = .TRUE. ! reflect vz at bottom/top-boundaries
    END SELECT
  END SUBROUTINE ReflectionMasks

  !> return masks for axis boundaries
  !! \warning Not rigorously tested!
  !!
  !! At axis boundaries we change the sign of normal velocities as in reflecting
  !! boundary conditions and in addition the sign of the tangential velocity in
  !! the plane perpendicular to the axis is changed.
  PURE SUBROUTINE AxisMasks(this,Mesh,reflX,reflY,reflZ)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm),  INTENT(IN)  :: this
    CLASS(mesh_base),              INTENT(IN)  :: Mesh       
    LOGICAL, DIMENSION(this%VNUM), INTENT(OUT) :: reflX,reflY,reflZ
    !------------------------------------------------------------------------!
    INTEGER :: aidx
    !------------------------------------------------------------------------!
    ! set masks according to reflecting boundaries, i. e. change sign
    ! of normal velocities
    CALL this%ReflectionMasks(Mesh,reflX,reflY,reflZ)
    ! get coordinate index of azimuthal angle (=0 if there is no azimuthal angle)
    aidx = Mesh%geometry%GetAzimuthIndex()
    IF (aidx.GT.0) THEN
      ! geometry has azimuthal angle -> determine which velocity changes sign
      SELECT CASE(this%VDIM)
      CASE(1) ! 1D velocity -> do nothing
      CASE(2) ! 2D velocity
      IF (Mesh%KNUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.3) THEN
        ! transport in x-y-plane
        SELECT CASE(aidx)
        CASE(1) ! 1st coordinate is the azimuthal angle
          reflY(2) = .TRUE. ! vx is tangential at south/north-boundaries
        CASE(2) ! 2nd coordinate is the azimuthal angle
          reflX(3) = .TRUE. ! vy is tangential at east/west-boundaries
        CASE(3) ! 3rd coordinate is the azimuthal angle
          ! do nothing
        END SELECT
      ELSE IF (Mesh%JNUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.2) THEN
        ! transport in x-z-plane
        SELECT CASE(aidx)
        CASE(1) ! 1st coordinate is the azimuthal angle
          reflZ(2) = .TRUE. ! vx is tangential at bottom/top-boundaries
        CASE(2) ! 2nd coordinate is the azimuthal angle
          ! do nothing
        CASE(3) ! 3rd coordinate is the azimuthal angle
          reflX(3) = .TRUE. ! vz is tangential at east/west-boundaries
        END SELECT
      ELSE IF (Mesh%INUM.EQ.1.AND..NOT.Mesh%ROTSYM.EQ.1) THEN
        ! transport in y-z-plane
        SELECT CASE(aidx)
        CASE(1) ! 1st coordinate is the azimuthal angle
          ! do nothing
        CASE(2) ! 2nd coordinate is the azimuthal angle
          reflZ(2) = .TRUE. ! vy is tangential at bottom/top-boundaries
        CASE(3) ! 3rd coordinate is the azimuthal angle
          reflY(3) = .TRUE. ! vz is tangential at south/north-boundaries
        END SELECT
      END IF
      CASE(3) ! 3D velocity
        SELECT CASE(aidx)
        CASE(1) ! 1st coordinate is the azimuthal angle
          reflZ(2) = .TRUE. ! vx is tangential at bottom/top-boundaries
        CASE(2) ! 2nd coordinate is the azimuthal angle
          reflX(3) = .TRUE. ! vy is tangential at east/west-boundaries
        CASE(3) ! 3rd coordinate is the azimuthal angle
          reflY(4) = .TRUE. ! vz is tangential at south/north-boundaries
        END SELECT
      END SELECT
    END IF
  END SUBROUTINE AxisMasks

  PURE SUBROUTINE CalculateCharSystemX(this,Mesh,i,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_eulerisotherm), INTENT(IN):: this
    CLASS(Mesh_base), INTENT(IN)   :: Mesh
    INTEGER, INTENT(IN)          :: i,dir
    REAL, INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    REAL, &
      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, INTENT(OUT), &
      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: xvar
    !------------------------------------------------------------------------!
    INTEGER           :: i1,i2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at i
    CALL SetEigenValues( &
          this%bccsound%data3d(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
          pvar(i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
    ! compute characteristic variables using cell mean values of adjacent
    ! cells to calculate derivatives and the isothermal speed of sound
    ! at the intermediate cell face
    i1 = i + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    i2 = MAX(i,i1)
    i1 = MIN(i,i1)
    CALL SetCharVars( &
          this%fcsound%data4d(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,WEST), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3))
  END SUBROUTINE CalculateCharSystemX


  PURE SUBROUTINE CalculateCharSystemY(this,Mesh,j,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_eulerisotherm),INTENT(IN) :: this
    CLASS(Mesh_base),INTENT(IN)    :: Mesh
    INTEGER,INTENT(IN)             :: j,dir
    REAL, INTENT(IN), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    REAL, &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, INTENT(OUT), &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: xvar
    !------------------------------------------------------------------------!
    INTEGER           :: j1,j2
    !------------------------------------------------------------------------!
    ! compute eigenvalues at j
    CALL SetEigenValues(this%bccsound%data3d(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX), &
          pvar(Mesh%IMIN:Mesh%IMAX,j,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
    ! compute characteristic variables using cell mean values of adjacent
    ! cells to calculate derivatives and the isothermal speed of sound
    ! at the intermediate cell face
    j1 = j + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
    j2 = MAX(j,j1)
    j1 = MIN(j,j1)
    CALL SetCharVars(this%fcsound%data4d(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,SOUTH), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
  END SUBROUTINE CalculateCharSystemY


  PURE SUBROUTINE CalculateCharSystemZ(this,Mesh,k,dir,pvar,lambda,xvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_eulerisotherm), INTENT(IN):: this
    CLASS(Mesh_base), INTENT(IN)   :: Mesh
    INTEGER, INTENT(IN)          :: k,dir
    REAL, INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    REAL, &
       DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM),INTENT(OUT) :: lambda
    REAL, INTENT(OUT), & 
     DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
    !------------------------------------------------------------------------!
!     INTEGER           :: k1,k2
    !------------------------------------------------------------------------!
    !TODO: Should not exist in 2D !
    ! compute eigenvalues at k
!    CALL SetEigenValues(this%bccsound%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k,this%ZVELOCITY), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          lambda(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3))
!    ! compute characteristic variables using cell mean values of adjacent
!    ! cells to calculate derivatives and the isothermal speed of sound
!    ! at the intermediate cell face
!    k1 = k + SIGN(1,dir) ! left handed if dir<0 and right handed otherwise
!    k2 = MAX(k,k1)
!    k1 = MIN(k,k1)
!    CALL SetCharVars(this%fcsound%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,BOTTOM), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%XVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%XVELOCITY), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3))
  END SUBROUTINE CalculateCharSystemZ

  PURE SUBROUTINE CalculateBoundaryDataX(this,Mesh,i1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_eulerisotherm), INTENT(IN):: this
    CLASS(Mesh_base), INTENT(IN)   :: Mesh
    INTEGER, INTENT(IN)          :: i1,dir
    REAL, INTENT(IN), &
      DIMENSION(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: xvar
    REAL, INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i2,fidx
    !------------------------------------------------------------------------!
    i2 = i1 + SIGN(1,dir)  ! i +/- 1 depending on the sign of dir
    IF (i2.LT.i1) THEN
       fidx = WEST
    ELSE
       fidx = EAST
    END IF
    CALL SetBoundaryData( &
          this%fcsound%data4d(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KGMIN:Mesh%KGMAX,fidx), &
          1.0*SIGN(1,dir), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,3), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          pvar(i2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY))
  END SUBROUTINE CalculateBoundaryDataX


  PURE SUBROUTINE CalculateBoundaryDataY(this,Mesh,j1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_eulerisotherm), INTENT(IN):: this
    CLASS(Mesh_base), INTENT(IN)   :: Mesh
    INTEGER, INTENT(IN)          :: j1,dir
    REAL, INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: xvar
    REAL, INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: j2,fidx
    !------------------------------------------------------------------------!
    j2 = j1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
    IF (j2.LT.j1) THEN
       fidx = SOUTH
    ELSE
       fidx = NORTH
    END IF
    CALL SetBoundaryData( &
          this%fcsound%data4d(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,fidx), &
          1.0*SIGN(1,dir), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j1,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,1), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,2), &
          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%KMIN:Mesh%KMAX,3), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%DENSITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%YVELOCITY), &
          pvar(Mesh%IMIN:Mesh%IMAX,j2,Mesh%KMIN:Mesh%KMAX,this%XVELOCITY))
  END SUBROUTINE CalculateBoundaryDataY

  PURE SUBROUTINE CalculateBoundaryDataZ(this,Mesh,k1,dir,xvar,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_eulerisotherm), INTENT(IN):: this
    CLASS(Mesh_base), INTENT(IN)   :: Mesh
    INTEGER, INTENT(IN)          :: k1,dir
    REAL, INTENT(IN), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,this%VNUM) :: xvar
    REAL, INTENT(INOUT), &
      DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,this%VNUM) :: pvar
    !------------------------------------------------------------------------!
!     INTEGER           :: k2,fidx
    !------------------------------------------------------------------------!
    !TODO: Should not exist in 2D !
    !    k2 = k1 + SIGN(1,dir)  ! j +/- 1 depending on the sign of dir
!    IF (k2.LT.k1) THEN
!       fidx = BOTTOM
!    ELSE
!       fidx = TOP
!    END IF
!    CALL SetBoundaryData( &
!          this%fcsound%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,fidx), &
!          1.0*SIGN(1,dir), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k1,this%XVELOCITY), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2), &
!          xvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,3), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%DENSITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%YVELOCITY), &
!          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k2,this%XVELOCITY))
  END SUBROUTINE CalculateBoundaryDataZ

  !> \public Destructor of the physics_eulerisotherm class
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%bccsound%Destroy() ! call destructor of the mesh array
    CALL this%fcsound%Destroy()
    DEALLOCATE(this%bccsound,this%fcsound)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

!----------------------------------------------------------------------------!
!> \par methods of class statevector_eulerisotherm

  !> \public Constructor of statevector_eulerisotherm
  !!
  !! \attention This is not a class member itself, instead its an ordinary
  !!            module procedure. The function name is overloaded with
  !!            the class name.
  FUNCTION CreateStateVector(Physics,flavour,num) RESULT(new_sv)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN) :: Physics
    INTEGER, OPTIONAL, INTENT(IN) :: flavour,num
    TYPE(statevector_eulerisotherm) :: new_sv
    !-------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized()) &
      CALL Physics%Error("physics_eulerisotherm::CreateStatevector", "Physics not initialized.")

    ! create a new empty compound of marrays
    new_sv = marray_compound(num)
    IF (PRESENT(flavour)) THEN
      SELECT CASE(flavour)
      CASE(PRIMITIVE,CONSERVATIVE)
        new_sv%flavour = flavour
      CASE DEFAULT
        new_sv%flavour = UNDEFINED
      END SELECT
    END IF
    SELECT CASE(new_sv%flavour)
    CASE(PRIMITIVE)
      ! allocate memory for density and velocity mesh arrays
      ALLOCATE(new_sv%density,new_sv%velocity)
      IF (PRESENT(num).AND.num.GT.0) THEN
        ! create a bunch of scalars and vectors
        new_sv%density  = marray_base(num)              ! num scalars
        new_sv%velocity = marray_base(num,Physics%VDIM) ! num vectors
      ELSE
        new_sv%density  = marray_base()                 ! one scalar
        new_sv%velocity = marray_base(Physics%VDIM)     ! one vector
      END IF
      ! append to compound
      CALL new_sv%AppendMArray(new_sv%density)
      CALL new_sv%AppendMArray(new_sv%velocity)
    CASE(CONSERVATIVE)
      ! allocate memory for density and momentum mesh arrays
      ALLOCATE(new_sv%density,new_sv%momentum)
      IF (PRESENT(num).AND.num.GT.0) THEN
        ! create a bunch of scalars and vectors
        new_sv%density  = marray_base(num)              ! num scalars
        new_sv%momentum = marray_base(num,Physics%VDIM) ! num vectors
      ELSE
        new_sv%density  = marray_base()                 ! one scalar
        new_sv%momentum = marray_base(Physics%VDIM)     ! one vector
      END IF
      ! append to compound
      CALL new_sv%AppendMArray(new_sv%density)
      CALL new_sv%AppendMArray(new_sv%momentum)
    CASE DEFAULT
      CALL Physics%Warning("physics_eulerisotherm::CreateStateVector", "Empty state vector created.")
    END SELECT
  END FUNCTION CreateStateVector

  SUBROUTINE CreateStateVector_eulerisotherm(Physics,new_sv,flavour,num)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(physics_eulerisotherm), INTENT(IN) :: Physics
    TYPE(statevector_eulerisotherm), INTENT(INOUT) :: new_sv
    INTEGER, OPTIONAL, INTENT(IN) :: flavour,num
    !-------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized()) &
      CALL Physics%Error("physics_eulerisotherm::CreateStatevector", "Physics not initialized.")

    ! create a new empty compound of marrays
    new_sv = marray_compound(num)
    IF (PRESENT(flavour)) THEN
      SELECT CASE(flavour)
      CASE(PRIMITIVE,CONSERVATIVE)
        new_sv%flavour = flavour
      CASE DEFAULT
        new_sv%flavour = UNDEFINED
      END SELECT
    END IF
    SELECT CASE(new_sv%flavour)
    CASE(PRIMITIVE)
      ! allocate memory for density and velocity mesh arrays
      ALLOCATE(new_sv%density,new_sv%velocity)
      ! create a bunch of scalars and vectors
      IF (PRESENT(num)) THEN
        new_sv%density  = marray_base(num)              ! num scalars
        new_sv%velocity = marray_base(num,Physics%VDIM) ! num vectors
      ELSE
        new_sv%density  = marray_base()                 ! one scalar
        new_sv%velocity = marray_base(Physics%VDIM)     ! one vector
      END IF      
      ! append to compound
      CALL new_sv%AppendMArray(new_sv%density)
      CALL new_sv%AppendMArray(new_sv%velocity)
    CASE(CONSERVATIVE)
      ! allocate memory for density and momentum mesh arrays
      ALLOCATE(new_sv%density,new_sv%momentum)
      ! create a bunch of scalars and vectors
      IF (PRESENT(num)) THEN
        new_sv%density  = marray_base(num)              ! num scalars
        new_sv%momentum = marray_base(num,Physics%VDIM) ! num vectors
      ELSE
        new_sv%density  = marray_base()                 ! one scalar
        new_sv%momentum = marray_base(Physics%VDIM)     ! one vector
      END IF
      ! append to compound
      CALL new_sv%AppendMArray(new_sv%density)
      CALL new_sv%AppendMArray(new_sv%momentum)
    CASE DEFAULT
      CALL Physics%Warning("physics_eulerisotherm::CreateStateVector_eulerisotherm", &
                           "Empty state vector created.")
    END SELECT
  END SUBROUTINE CreateStateVector_eulerisotherm

  !> assigns one state vector to another state vector
  SUBROUTINE AssignMArray_0(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(statevector_eulerisotherm),INTENT(INOUT) :: this
    CLASS(marray_base),INTENT(IN)    :: ma
    !------------------------------------------------------------------------!
    CALL this%marray_compound%AssignMArray_0(ma)
    IF (SIZE(this%data1d).GT.0) THEN
      SELECT TYPE(src => ma)
      CLASS IS(statevector_eulerisotherm)
        ! copy flavour
        this%flavour = src%flavour
        ! set pointer to the data structures in the compound
        !> \todo make this more generic, i.e. this should not depend
        !! on the position in the list of compound items
        this%density => this%GetItem(this%FirstItem()) ! density is the first item
        SELECT CASE(this%flavour)
        CASE(PRIMITIVE)
          ! velocity is the second item
          this%velocity => this%GetItem(this%NextItem(this%FirstItem()))
        CASE(CONSERVATIVE)
          ! momentum is the second item
          this%momentum => this%GetItem(this%NextItem(this%FirstItem()))
        CASE DEFAULT
          ! error, this should not happen
        END SELECT
      CLASS DEFAULT
        ! error, this should not happen
      END SELECT
    END IF
  END SUBROUTINE AssignMArray_0


!----------------------------------------------------------------------------!
!> \par elemental non-class subroutines / functions

  !> \private set minimal and maximal wave speeds
  ELEMENTAL SUBROUTINE SetWaveSpeeds(cs,v,minwav,maxwav)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: minwav,maxwav
    !------------------------------------------------------------------------!
    minwav = MIN(0.,v-cs)  ! minimal wave speed
    maxwav = MAX(0.,v+cs)  ! maximal wave speed
  END SUBROUTINE SetWaveSpeeds

  !> \private compute all eigenvalues
  ELEMENTAL SUBROUTINE SetEigenValues(cs,v,l1,l2,l3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,v
    REAL, INTENT(OUT) :: l1,l2,l3
    !------------------------------------------------------------------------!
    l1 = v - cs
    l2 = v
    l3 = v + cs
  END SUBROUTINE SetEigenValues

  !> \private set mass and 1D momentum flux for transport along the 1st dimension
  ELEMENTAL SUBROUTINE SetFlux1d(cs,rho,u,mu,f1,f2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,u,mu
    REAL, INTENT(OUT) :: f1,f2
    !------------------------------------------------------------------------!
    f1 = rho*u
    f2 = mu*u + rho*cs*cs
  END SUBROUTINE SetFlux1d

  !> \private set mass and 2D momentum flux for transport along the 1st dimension
  ELEMENTAL SUBROUTINE SetFlux2d(cs,rho,u,mu,mv,f1,f2,f3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,u,mu,mv
    REAL, INTENT(OUT) :: f1,f2,f3
    !------------------------------------------------------------------------!
    CALL SetFlux1d(cs,rho,u,mu,f1,f2)
    f3 = mv*u
  END SUBROUTINE SetFlux2d

  !> \private set mass and 3D momentum flux for transport along the 1st dimension
  ELEMENTAL SUBROUTINE SetFlux3d(cs,rho,u,mu,mv,mw,f1,f2,f3,f4)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho,u,mu,mv,mw
    REAL, INTENT(OUT) :: f1,f2,f3,f4
    !------------------------------------------------------------------------!
    CALL SetFlux2d(cs,rho,u,mu,mv,f1,f2,f3)
    f4 = mw*u
  END SUBROUTINE SetFlux3d

  !> \private Convert from 1D conservative to primitive variables
  ELEMENTAL SUBROUTINE Cons2Prim1d(rho_in,mu,rho_out,u)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,mu
    REAL, INTENT(OUT) :: rho_out,u
    !------------------------------------------------------------------------!
    rho_out = rho_in
    u       = mu / rho_in
  END SUBROUTINE Cons2Prim1d

  !> \private Convert from 2D conservative to primitive variables
  ELEMENTAL SUBROUTINE Cons2Prim2d(rho_in,mu,mv,rho_out,u,v)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,mu,mv
    REAL, INTENT(OUT) :: rho_out,u,v
    !------------------------------------------------------------------------!
    REAL :: inv_rho
    !------------------------------------------------------------------------!
    inv_rho = 1./rho_in
    rho_out = rho_in
    u       = mu * inv_rho
    v       = mv * inv_rho
  END SUBROUTINE Cons2Prim2d

  !> \private Convert from 3D conservative to primitive variables
  ELEMENTAL SUBROUTINE Cons2Prim3d(rho_in,mu,mv,mw,rho_out,u,v,w)
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
  END SUBROUTINE Cons2Prim3d

  !> \private Convert from 1D primitive to conservative variables
  ELEMENTAL SUBROUTINE Prim2Cons1d(rho_in,u,rho_out,mu)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u
    REAL, INTENT(OUT) :: rho_out,mu
    !------------------------------------------------------------------------!
    rho_out = rho_in
    mu = rho_in * u
  END SUBROUTINE Prim2Cons1d

  !> \private Convert from 2D primitive to conservative variables
  ELEMENTAL SUBROUTINE Prim2Cons2d(rho_in,u,v,rho_out,mu,mv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u,v
    REAL, INTENT(OUT) :: rho_out,mu,mv
    !------------------------------------------------------------------------!
    CALL Prim2Cons1d(rho_in,u,rho_out,mu)
    mv = rho_in * v
  END SUBROUTINE Prim2Cons2d

  !> \private Convert from 3D primitive to conservative variables
  ELEMENTAL SUBROUTINE Prim2Cons3d(rho_in,u,v,w,rho_out,mu,mv,mw)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: rho_in,u,v,w
    REAL, INTENT(OUT) :: rho_out,mu,mv,mw
    !------------------------------------------------------------------------!
    CALL Prim2Cons2d(rho_in,u,v,rho_out,mu,mv)
    mw = rho_in * w
  END SUBROUTINE Prim2Cons3d

  !> \private compute characteristic variables using adjacent primitve states
  ELEMENTAL SUBROUTINE SetCharVars(cs,rho1,rho2,u1,u2,v1,v2,&
       xvar1,xvar2,xvar3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,rho1,rho2,u1,u2,v1,v2
    REAL, INTENT(OUT) :: xvar1,xvar2,xvar3
    !------------------------------------------------------------------------!
    REAL :: dlnrho,du
    !------------------------------------------------------------------------!
    dlnrho = LOG(rho2/rho1)
    du = u2-u1
    ! characteristic variables
    xvar1 = cs*dlnrho - du
    xvar2 = v2-v1
    xvar3 = cs*dlnrho + du
  END SUBROUTINE SetCharVars

  !> \private extrapolate boundary values using primitve and characteristic variables
  ELEMENTAL SUBROUTINE SetBoundaryData(cs,dir,rho1,u1,v1,xvar1,xvar2, &
       xvar3,rho2,u2,v2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cs,dir,rho1,u1,v1,xvar1,xvar2,xvar3
    REAL, INTENT(OUT) :: rho2,u2,v2
    !------------------------------------------------------------------------!
    rho2 = rho1 * EXP(dir*0.5*(xvar3+xvar1)/cs)
    u2   = u1 + dir*0.5*(xvar3-xvar1)
    v2   = v1 + dir*xvar2
  END SUBROUTINE SetBoundaryData

  !> geometrical momentum source terms
  !! P is the either isothermal pressure rho*cs**2 or the real pressure.
  !!
  !! \attention These elemental functions exist multiple times for performance
  !!  reasons (inlining). Please keep this in mind for changes.
  !!  Other modules with this function:
  !!      - physics_euler_mod
  !!
  !! x-momentum geometrical source term
  ELEMENTAL FUNCTION GetGeometricalSourceX(cxyx,cxzx,cyxy,czxz,vx,vy,vz,P,my,mz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cxyx,cxzx,cyxy,czxz,vx,vy,vz,P,my,mz
    REAL :: GetGeometricalSourceX
    !------------------------------------------------------------------------!
    GetGeometricalSourceX = my * (cyxy*vy - cxyx*vx) + mz * (czxz*vz - cxzx*vx) &
                            + (cyxy + czxz) * P
  END FUNCTION GetGeometricalSourceX

  !> y-momentum geometrical source term
  ELEMENTAL FUNCTION GetGeometricalSourceY(cxyx,cyxy,cyzy,czyz,vx,vy,vz,P,mx,mz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cxyx,cyxy,cyzy,czyz,vx,vy,vz,P,mx,mz
    REAL :: GetGeometricalSourceY
    !------------------------------------------------------------------------!
    GetGeometricalSourceY = mz * (czyz*vz - cyzy*vy) + mx * (cxyx*vx - cyxy*vy) &
                            + (cxyx + czyz) * P
  END FUNCTION GetGeometricalSourceY

  !> z-momentum geometrical source term
  ELEMENTAL FUNCTION GetGeometricalSourceZ(cxzx,cyzy,czxz,czyz,vx,vy,vz,P,mx,my)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: cxzx,cyzy,czxz,czyz,vx,vy,vz,P,mx,my
    REAL :: GetGeometricalSourceZ
    !------------------------------------------------------------------------!
    GetGeometricalSourceZ = mx * (cxzx*vx - czxz*vz) + my * (cyzy*vy - czyz*vz) &
                            + (cxzx + cyzy) * P
  END FUNCTION GetGeometricalSourceZ

!  ! TODO: Not verified
!  !!
!  !! non-global elemental routine
!  ELEMENTAL SUBROUTINE SetRoeAverages(rhoL,rhoR,vL,vR,v)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: rhoL,rhoR,vL,vR
!    REAL, INTENT(OUT) :: v
!    !------------------------------------------------------------------------!
!    REAL :: sqrtrhoL,sqrtrhoR
!    !------------------------------------------------------------------------!
!    sqrtrhoL = SQRT(rhoL)
!    sqrtrhoR = SQRT(rhoR)
!    v = 0.5*(sqrtrhoL*vL + sqrtrhoR*vR) / (sqrtrhoL + sqrtrhoR)
!  END SUBROUTINE SetRoeAverages

END MODULE physics_eulerisotherm_mod
