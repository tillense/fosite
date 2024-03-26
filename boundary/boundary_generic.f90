!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2006-2016                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!> \addtogroup boundary
!! \key{western,INTEGER,boundary condition at western boundary}
!! \key{eastern,INTEGER,boundary condition at eastern boundary}
!! \key{southern,INTEGER,boundary condition at southern boundary}
!! \key{northern,INTEGER,boundary condition at northern boundary}
!! \key{bottomer,INTEGER,boundary condition at bottomer boundary}
!! \key{topper,INTEGER,boundary condition at topper boundary}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief Generic boundary module
!!
!! This module and its object holds the six boundaries, where every boundary
!! is its own object.
!----------------------------------------------------------------------------!
MODULE boundary_generic_mod
  USE logging_base_mod
  USE marray_compound_mod
  USE mesh_base_mod
  USE boundary_base_mod
  USE boundary_custom_mod
  USE boundary_reflecting_mod
  USE boundary_nogradients_mod
  USE boundary_periodic_mod
  USE boundary_inner_mod
  USE boundary_axis_mod
  USE boundary_absorbing_mod
  USE boundary_fixed_mod
  USE boundary_noslip_mod
  USE boundary_shearing_mod
  USE boundary_farfield_mod
  USE physics_base_mod
  USE common_dict
  !--------------------------------------------------------------------------!
  TYPE, PRIVATE                       :: boundary_p
    CLASS(boundary_base), ALLOCATABLE :: p
  END TYPE

  TYPE, EXTENDS(logging_base)         :: boundary_generic
    !> \name variables
    TYPE(boundary_p)                  :: Boundary(6)
  CONTAINS
    PROCEDURE :: InitBoundary
    PROCEDURE :: CenterBoundary
    PROCEDURE :: SetCornerEdges
#ifdef PARALLEL
    PROCEDURE :: MPIBoundaryCommunication
    PROCEDURE :: InitBoundary_MPI
    PROCEDURE :: MPIbuffer2pvar
    PROCEDURE :: MPIpvar2buffer
#endif
    PROCEDURE :: Finalize
  END TYPE boundary_generic
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       boundary_base, boundary_generic, new_boundary
       ! constants
  PRIVATE :: InitBoundary, CenterBoundary, Finalize
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE new_boundary(Boundary,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_generic), ALLOCATABLE    :: Boundary
    CLASS(mesh_base),        INTENT(INOUT)  :: Mesh
    CLASS(physics_base),     INTENT(IN)     :: Physics
    TYPE(Dict_TYP), POINTER                 :: config,IO
    !------------------------------------------------------------------------!
    ALLOCATE(Boundary)
    CALL Boundary%InitBoundary(Mesh,Physics,config,IO)
  END SUBROUTINE

  SUBROUTINE InitBoundary(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_generic), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(INOUT) :: Mesh
    CLASS(physics_base),     INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER                :: config,IO
    INTEGER         :: western, eastern, southern, northern, bottomer, topper
    !------------------------------------------------------------------------!
    INTEGER               :: new(6)
    INTEGER               :: western_std, eastern_std, northern_std, southern_std, &
                             topper_std, bottomer_std
    LOGICAL, DIMENSION(3) :: periods = .FALSE.
    INTEGER               :: dir
    !------------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("InitBoundary","physics and/or mesh module uninitialized")

    ! check for shearingsheet standard boundaries
    western_std = 0
    eastern_std = 0
    northern_std = 0
    southern_std = 0
    bottomer_std = 0
    topper_std = 0
    IF (Mesh%shear_dir.EQ.1) THEN
      western_std = PERIODIC
      eastern_std = PERIODIC
      southern_std = SHEARING
      northern_std = SHEARING
      bottomer_std = REFLECTING
      topper_std = REFLECTING
    ELSE IF (Mesh%shear_dir.EQ.2) THEN
      western_std = SHEARING
      eastern_std = SHEARING
      southern_std = PERIODIC
      northern_std = PERIODIC
      bottomer_std = REFLECTING
      topper_std = REFLECTING
    END IF

    CALL GetAttr(config, "western",   western, western_std)
    CALL GetAttr(config, "eastern",   eastern, eastern_std)
    CALL GetAttr(config, "southern", southern, northern_std)
    CALL GetAttr(config, "northern", northern, southern_std)
    CALL GetAttr(config, "bottomer", bottomer, bottomer_std)
    CALL GetAttr(config, "topper",     topper, topper_std)

    new(WEST)   = western
    new(EAST)   = eastern
    new(SOUTH)  = southern
    new(NORTH)  = northern
    new(BOTTOM) = bottomer
    new(TOP)    = topper

#ifdef PARALLEL
    ! define inner connections, where boundaries are no true ones
    IF (Mesh%mycoords(1).NE.0)              new(WEST)   = NONE
    IF (Mesh%mycoords(1).NE.Mesh%dims(1)-1) new(EAST)   = NONE
    IF (Mesh%mycoords(2).NE.0)              new(SOUTH)  = NONE
    IF (Mesh%mycoords(2).NE.Mesh%dims(2)-1) new(NORTH)  = NONE
    IF (Mesh%mycoords(3).NE.0)              new(BOTTOM) = NONE
    IF (Mesh%mycoords(3).NE.Mesh%dims(3)-1) new(TOP)    = NONE
#endif

    ! Check for correct shifting and boundaries
    IF (western.EQ.SHEARING.AND.eastern.EQ.SHEARING) THEN
      IF (.NOT.Mesh%shear_dir.EQ.2) &
        CALL this%Error("InitBoundary", &
          "Please apply shifting in second dimension, when applying shearing boundaries at western/eastern.")
#ifdef PARALLEL
      CALL this%Error("InitBoundary", &
        "Parallel mode is not allowed with shearing in West-East direction.")
#endif
    ELSE IF (southern.EQ.SHEARING.AND.northern.EQ.SHEARING) THEN
      IF (.NOT.Mesh%shear_dir.EQ.1) &
        CALL this%Error("InitBoundary", &
          "Please apply shifting in first dimension, when applying shearing boundaries at northern/southern.")
    ELSE IF (bottomer.EQ.SHEARING.AND.topper.EQ.SHEARING) THEN
        CALL this%Error("InitBoundary", &
          "shifting in topper/bottomer direction not allowed.")
    END IF

    ! initialize every boundary
    ! IMPORTANT: do this before anything else
    DO dir=WEST,TOP
      SELECT CASE(new(dir))
      CASE(ABSORBING)
        ALLOCATE(boundary_absorbing::this%Boundary(dir)%p)
      CASE(AXIS)
        ALLOCATE(boundary_axis::this%Boundary(dir)%p)
      CASE(CUSTOM)
        ALLOCATE(boundary_custom::this%Boundary(dir)%p)
      CASE(FARFIELD)
        ALLOCATE(boundary_farfield::this%Boundary(dir)%p)
      CASE(FIXED)
        ALLOCATE(boundary_fixed::this%Boundary(dir)%p)
      CASE(NO_GRADIENTS)
        ALLOCATE(boundary_nogradients::this%Boundary(dir)%p)
      CASE(NOSLIP)
        ALLOCATE(boundary_noslip::this%boundary(dir)%p)
      CASE(PERIODIC)
        ALLOCATE(boundary_periodic::this%Boundary(dir)%p)
      CASE(REFLECTING)
        ALLOCATE(boundary_reflecting::this%Boundary(dir)%p)
      CASE(SHEARING)
        ALLOCATE(boundary_shearing::this%Boundary(dir)%p)
#ifdef PARALLEL
      CASE(NONE)
        ALLOCATE(boundary_inner::this%Boundary(dir)%p)
#endif
      CASE DEFAULT
        CALL this%Error("boundary_generic::InitBoundary","unknown boundary condition")
      END SELECT

      SELECT TYPE(obj => this%Boundary(dir)%p)
      TYPE IS (boundary_absorbing)
        CALL obj%InitBoundary_absorbing(Mesh,Physics,dir,config)
      TYPE IS (boundary_axis)
        CALL obj%InitBoundary_axis(Mesh,Physics,dir,config)
      TYPE IS (boundary_custom)
        CALL obj%InitBoundary_custom(Mesh,Physics,dir,config)
      TYPE IS (boundary_farfield)
        CALL obj%InitBoundary_farfield(Mesh,Physics,dir,config)
      TYPE IS (boundary_fixed)
        CALL obj%InitBoundary_fixed(Mesh,Physics,dir,config)
      TYPE IS (boundary_nogradients)
        CALL obj%InitBoundary_nogradients(Mesh,Physics,dir,config)
      TYPE IS (boundary_noslip)
        CALL obj%InitBoundary_noslip(Mesh,Physics,dir,config)
      TYPE IS (boundary_periodic)
        CALL obj%InitBoundary_periodic(Mesh,Physics,dir,config)
      TYPE IS (boundary_reflecting)
        CALL obj%InitBoundary_reflecting(Mesh,Physics,dir,config)
      TYPE IS (boundary_shearing)
        CALL obj%InitBoundary_shearing(Mesh,Physics,dir,config)
#ifdef PARALLEL
      TYPE IS (boundary_inner)
        CALL obj%InitBoundary_inner(Mesh,Physics,dir,config)
#endif
      END SELECT
    END DO

    ! check periodicity
    IF ((western.EQ.PERIODIC.AND.eastern.EQ.PERIODIC) .OR. &
        (western.EQ.SHEARING.AND.eastern.EQ.SHEARING)) THEN
        periods(1) = .TRUE.
    ELSE IF (western.EQ.PERIODIC.NEQV.eastern.EQ.PERIODIC) THEN
       CALL this%boundary(WEST)%p%Error("InitBoundary", &
            "Opposite boundary should be periodic.")
    ELSE IF (western.EQ.SHEARING.NEQV.eastern.EQ.SHEARING) THEN
       CALL this%boundary(WEST)%p%Error("InitBoundary", &
            "Opposite boundary should be shearing.")
    END IF
    IF ((southern.EQ.PERIODIC.AND.northern.EQ.PERIODIC) .OR. &
        (southern.EQ.SHEARING.AND.northern.EQ.SHEARING)) THEN
       periods(2) = .TRUE.
    ELSE IF (southern.EQ.PERIODIC.NEQV.northern.EQ.PERIODIC) THEN
       CALL this%boundary(SOUTH)%p%Error("InitBoundary", &
            "Opposite boundary should be periodic.")
    ELSE IF (southern.EQ.SHEARING.NEQV.northern.EQ.SHEARING) THEN
       CALL this%boundary(SOUTH)%p%Error("InitBoundary", &
            "Opposite boundary should be shearing.")
    END IF
    IF ((bottomer.EQ.PERIODIC.AND.topper.EQ.PERIODIC) .OR. &
        (bottomer.EQ.SHEARING.AND.topper.EQ.SHEARING)) THEN
       periods(3) = .TRUE.
    ELSE IF (bottomer.EQ.PERIODIC.NEQV.topper.EQ.PERIODIC) THEN
       CALL this%boundary(BOTTOM)%p%Error("InitBoundary", &
            "Opposite boundary should be periodic.")
    ELSE IF (bottomer.EQ.SHEARING.NEQV.topper.EQ.SHEARING) THEN
       CALL this%boundary(BOTTOM)%p%Error("InitBoundary", &
            "Opposite boundary should be shearing.")
    END IF

#ifdef PARALLEL
    CALL this%InitBoundary_MPI(Mesh,Physics,periods)
#endif

  END SUBROUTINE InitBoundary


  !> Sets boundaries in all directions
  SUBROUTINE CenterBoundary(this,Mesh,Physics,time,pvar,cvar)
    USE boundary_shearing_mod, ONLY : boundary_shearing
    USE physics_eulerisotherm_mod, ONLY : statevector_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    REAL,                   INTENT(IN)    :: time
    CLASS(marray_compound), INTENT(INOUT) :: pvar, cvar
    !------------------------------------------------------------------------!
    CALL Physics%Convert2Primitive(Mesh%IMIN,Mesh%IMAX,Mesh%JMIN, &
             Mesh%JMAX,Mesh%KMIN,Mesh%KMAX,cvar,pvar)
    SELECT TYPE(c => cvar)
    CLASS IS(statevector_eulerisotherm)
      SELECT TYPE(p => pvar)
      CLASS IS(statevector_eulerisotherm)
        p%fargo_transformation_applied = c%fargo_transformation_applied
      END SELECT
    END SELECT

    this%err = 0
    ! set physical boundary conditions at western and eastern boundaries
    IF (Mesh%INUM.GT.1) THEN
      CALL this%Boundary(WEST)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                    "western boundary condition failed")
      CALL this%Boundary(EAST)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                    "eastern boundary condition failed")
    END IF

    ! set physical boundary conditions at southern and northern boundaries
    !
    ! Here an extra case for parallel execution is handled, when shearing
    ! boundaries are applied in northern/southern direction. The application
    ! of the boundary conditions is done after the MPI communication, because
    ! the shearing boundaries need to first apply periodic boundary conditions.
    ! This is done MPI below across the computational domains.
#ifdef PARALLEL
    SELECT TYPE(bound1 => this%Boundary(SOUTH)%p)
    TYPE IS (boundary_shearing)
      ! Apply shearing boundaries after MPI communication. Compare with
      ! code passage further below.
    CLASS DEFAULT
      IF (Mesh%JNUM.GT.1) THEN
        CALL this%Boundary(SOUTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
        IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                      "southern boundary condition failed")
      END IF
    END SELECT

    SELECT TYPE(bound2 => this%Boundary(NORTH)%p)
    TYPE IS (boundary_shearing)
      ! Apply shearing boundaries after MPI communication. Compare with
      ! code passage further below.
    CLASS DEFAULT
      IF (Mesh%JNUM.GT.1) THEN
        CALL this%Boundary(NORTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
        IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                      "northern boundary condition failed")
      END IF
    END SELECT
#else
    IF (Mesh%JNUM.GT.1) THEN
      CALL this%Boundary(SOUTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
        IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                      "southern boundary condition failed")
      CALL this%Boundary(NORTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
        IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                      "northern boundary condition failed")
    END IF
#endif

    ! set physical boundary conditions at top and bottom boundaries
    IF (Mesh%KNUM.GT.1) THEN
      CALL this%Boundary(BOTTOM)%p%SetBoundaryData(Mesh,Physics,time,pvar)
        IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                      "bottom boundary condition failed")
      CALL this%Boundary(TOP)%p%SetBoundaryData(Mesh,Physics,time,pvar)
        IF (this%err.GT.0) CALL this%Error("boundary_generic::CenterBoundary", &
                                      "top boundary condition failed")
    END IF

#ifdef PARALLEL
    CALL MPIBoundaryCommunication(this,Mesh,Physics,pvar%data4d)
#endif

#ifdef PARALLEL
    ! set physical boundary conditions at southern and northern boundaries
    ! Here an extra case for parallel execution is handled. When shearing
    ! boundaries are applied, periodic boundaries are assumed.
    SELECT TYPE(bound1 => this%Boundary(SOUTH)%p)
    TYPE IS (boundary_shearing)
      IF (Mesh%JNUM.GT.1) THEN
        CALL this%Boundary(SOUTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      END IF
    CLASS DEFAULT
      ! do nothing
    END SELECT
    SELECT TYPE(bound2 => this%Boundary(NORTH)%p)
    TYPE IS (boundary_shearing)
      IF (Mesh%JNUM.GT.1) THEN
        CALL this%Boundary(NORTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      END IF
    CLASS DEFAULT
      ! do nothing
    END SELECT
#endif

    CALL SetCornerEdges(this,Mesh,Physics,pvar%data4d)

    ! convert primitive variables in ghost cells
    IF (Mesh%INUM.GT.1) THEN
      CALL Physics%Convert2Conservative(Mesh%IGMIN,Mesh%IMIN-Mesh%IP1, &
          Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
      CALL Physics%Convert2Conservative(Mesh%IMAX+Mesh%IP1,Mesh%IGMAX, &
          Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    END IF
    IF (Mesh%JNUM.GT.1) THEN
      CALL Physics%Convert2Conservative(Mesh%IMIN,Mesh%IMAX,    &
            Mesh%JGMIN,Mesh%JMIN-Mesh%JP1,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
      CALL Physics%Convert2Conservative(Mesh%IMIN,Mesh%IMAX,    &
            Mesh%JMAX+Mesh%JP1,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    END IF
    IF (Mesh%KNUM.GT.1) THEN
      CALL Physics%Convert2Conservative(Mesh%IMIN,Mesh%IMAX,    &
          Mesh%JMIN,Mesh%JMAX,Mesh%KGMIN,Mesh%KMIN-Mesh%KP1,pvar,cvar)
      CALL Physics%Convert2Conservative(Mesh%IMIN,Mesh%IMAX,    &
          Mesh%JMIN,Mesh%JMAX,Mesh%KMAX+Mesh%KP1,Mesh%KGMAX,pvar,cvar)
    END IF
  END SUBROUTINE CenterBoundary


  !> Calculates the corner in 2D and the corners and edges in 3D
  !!
  !! This is a interpolation of corners & edges outside
  !! the computational domain, if they are undefined (e.g. there are
  !! no periodic or inner boundaries involved in the corner)
  !! this is also necessary, because we need some of these values in the
  !! viscosity module
  !! This part calculates the corner in 2D and the edges in 3D. Further
  !! below the corners in 3D are approximated.
  !! \attention Only the diagonal corner cells in 3D corners are approximated
  !!            because the others are not necessary by any module (they are
  !!            also set reasonably by the setting the edges, but could be
  !!            approximated better).
  SUBROUTINE SetCornerEdges(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    !------------------------------------------------------------------------!
    IF ((Mesh%INUM.NE.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.EQ.1).OR. &
        (Mesh%INUM.NE.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.NE.1)) THEN

      ! south-west
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(2).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
          DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
            DO i=j+1,Mesh%GNUM
              pvar(Mesh%IMIN-i,Mesh%JMIN-j,:,:) = pvar(Mesh%IMIN-i,Mesh%JMIN-j+1,:,:)
              pvar(Mesh%IMIN-j,Mesh%JMIN-i,:,:) = pvar(Mesh%IMIN-j+1,Mesh%JMIN-i,:,:)
            END DO
            pvar(Mesh%IMIN-j,Mesh%JMIN-j,:,:) = 0.5 * (pvar(Mesh%IMIN-j,Mesh%JMIN,:,:) &
                 + pvar(Mesh%IMIN,Mesh%JMIN-j,:,:))
          END DO
#ifdef PARALLEL
      END IF
#endif
      ! south-east
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(2).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
          DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
            DO i=j+1,Mesh%GNUM
              pvar(Mesh%IMAX+i,Mesh%JMIN-j,:,:) = pvar(Mesh%IMAX+i,Mesh%JMIN-j+1,:,:)
              pvar(Mesh%IMAX+j,Mesh%JMIN-i,:,:) = pvar(Mesh%IMAX+j-1,Mesh%JMIN-i,:,:)
            END DO
            pvar(Mesh%IMAX+j,Mesh%JMIN-j,:,:) = 0.5 * (pvar(Mesh%IMAX+j,Mesh%JMIN,:,:) &
                 + pvar(Mesh%IMAX,Mesh%JMIN-j,:,:))
          END DO
#ifdef PARALLEL
      END IF
#endif
      ! north-west
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(2).EQ.Mesh%dims(2)-1) THEN
#endif
!NEC$ SHORTLOOP
          DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
            DO i=j+1,Mesh%GNUM
              pvar(Mesh%IMIN-i,Mesh%JMAX+j,:,:) = pvar(Mesh%IMIN-i,Mesh%JMAX+j-1,:,:)
              pvar(Mesh%IMIN-j,Mesh%JMAX+i,:,:) = pvar(Mesh%IMIN-j+1,Mesh%JMAX+i,:,:)
            END DO
            pvar(Mesh%IMIN-j,Mesh%JMAX+j,:,:) = 0.5 * (pvar(Mesh%IMIN-j,Mesh%JMAX,:,:) &
                 + pvar(Mesh%IMIN,Mesh%JMAX+j,:,:))
          END DO
#ifdef PARALLEL
      END IF
#endif

      ! north-east
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(2).EQ.Mesh%dims(2)-1) THEN
#endif
!NEC$ SHORTLOOP
          DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
            DO i=j+1,Mesh%GNUM
              pvar(Mesh%IMAX+i,Mesh%JMAX+j,:,:) = pvar(Mesh%IMAX+i,Mesh%JMAX+j-1,:,:)
              pvar(Mesh%IMAX+j,Mesh%JMAX+i,:,:) = pvar(Mesh%IMAX+j-1,Mesh%JMAX+i,:,:)
            END DO
            pvar(Mesh%IMAX+j,Mesh%JMAX+j,:,:) = 0.5 * (pvar(Mesh%IMAX+j,Mesh%JMAX,:,:) &
                 + pvar(Mesh%IMAX,Mesh%JMAX+j,:,:))
          END DO
#ifdef PARALLEL
      END IF
#endif
    END IF

    IF ((Mesh%INUM.NE.1.AND.Mesh%JNUM.EQ.1.AND.Mesh%KNUM.NE.1).OR. &
        (Mesh%INUM.NE.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.NE.1)) THEN

      ! bottom-west
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
          DO i=1,Mesh%GNUM
!NEC$ SHORTLOOP
            DO k=i+1,Mesh%GNUM
              pvar(Mesh%IMIN-i,:,Mesh%KMIN-k,:) = pvar(Mesh%IMIN-i+1,:,Mesh%KMIN-k,:)
              pvar(Mesh%IMIN-k,:,Mesh%KMIN-i,:) = pvar(Mesh%IMIN-k,:,Mesh%KMIN-i+1,:)
            END DO
            pvar(Mesh%IMIN-i,:,Mesh%KMIN-i,:) = 0.5 * (pvar(Mesh%IMIN-i,:,Mesh%KMIN,:) &
                 + pvar(Mesh%IMIN,:,Mesh%KMIN-i,:))
          END DO
#ifdef PARALLEL
      END IF
#endif
      ! bottom-east
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
        DO i=1,Mesh%GNUM
!NEC$ SHORTLOOP
          DO k=i+1,Mesh%GNUM
            pvar(Mesh%IMAX+i,:,Mesh%KMIN-k,:) = pvar(Mesh%IMAX+i-1,:,Mesh%KMIN-k,:)
            pvar(Mesh%IMAX+k,:,Mesh%KMIN-i,:) = pvar(Mesh%IMAX+k,:,Mesh%KMIN-i+1,:)
          END DO
          pvar(Mesh%IMAX+i,:,Mesh%KMIN-i,:) = 0.5 * (pvar(Mesh%IMAX,:,Mesh%KMIN-i,:) &
               + pvar(Mesh%IMAX+i,:,Mesh%KMIN,:))
        END DO
#ifdef PARALLEL
      END IF
#endif

      ! top-west
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
        DO i=1,Mesh%GNUM
!NEC$ SHORTLOOP
          DO k=i+1,Mesh%GNUM
            pvar(Mesh%IMIN-i,:,Mesh%KMAX+k,:) = pvar(Mesh%IMIN-i+1,:,Mesh%KMAX+k,:)
            pvar(Mesh%IMIN-k,:,Mesh%KMAX+i,:) = pvar(Mesh%IMIN-k,:,Mesh%KMAX+i-1,:)
          END DO
          pvar(Mesh%IMIN-i,:,Mesh%KMAX+i,:) = 0.5 * (pvar(Mesh%IMIN,:,Mesh%KMAX+i,:) &
               + pvar(Mesh%IMIN-i,:,Mesh%KMAX,:))
        END DO
#ifdef PARALLEL
      END IF
#endif

      ! top-east
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
        DO i=1,Mesh%GNUM
!NEC$ SHORTLOOP
          DO k=i+1,Mesh%GNUM
            pvar(Mesh%IMAX+i,:,Mesh%KMAX+k,:) = pvar(Mesh%IMAX+i-1,:,Mesh%KMAX+k,:)
            pvar(Mesh%IMAX+k,:,Mesh%KMAX+i,:) = pvar(Mesh%IMAX+k,:,Mesh%KMAX+i-1,:)
          END DO
          pvar(Mesh%IMAX+i,:,Mesh%KMAX+i,:) = 0.5 * (pvar(Mesh%IMAX,:,Mesh%KMAX+i,:) &
               + pvar(Mesh%IMAX+i,:,Mesh%KMAX,:))
        END DO
#ifdef PARALLEL
      END IF
#endif
    END IF

    IF ((Mesh%INUM.EQ.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.NE.1).OR. &
        (Mesh%INUM.NE.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.NE.1)) THEN

      ! bottom-south
#ifdef PARALLEL
      IF(Mesh%mycoords(2).EQ.0.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
        DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
          DO k=j+1,Mesh%GNUM
            pvar(:,Mesh%JMIN-j,Mesh%KMIN-k,:) = pvar(:,Mesh%JMIN-j+1,Mesh%KMIN-k,:)
            pvar(:,Mesh%JMIN-k,Mesh%KMIN-j,:) = pvar(:,Mesh%JMIN-k,Mesh%KMIN-j+1,:)
          END DO
          pvar(:,Mesh%JMIN-j,Mesh%KMIN-j,:) = 0.5 * (pvar(:,Mesh%JMIN,Mesh%KMIN-j,:) &
               + pvar(:,Mesh%JMIN-j,Mesh%KMIN,:))
        END DO
#ifdef PARALLEL
      END IF
#endif

      ! bottom-north
#ifdef PARALLEL
      IF(Mesh%mycoords(2).EQ.Mesh%dims(2)-1.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
        DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
          DO k=j+1,Mesh%GNUM
            pvar(:,Mesh%JMAX+j,Mesh%KMIN-k,:) = pvar(:,Mesh%JMAX+j-1,Mesh%KMIN-k,:)
            pvar(:,Mesh%JMAX+k,Mesh%KMIN-j,:) = pvar(:,Mesh%JMAX+k,Mesh%KMIN-j+1,:)
          END DO
!NEC$ IVDEP
          pvar(:,Mesh%JMAX+j,Mesh%KMIN-j,:) = 0.5 * (pvar(:,Mesh%JMAX,Mesh%KMIN-j,:) &
               + pvar(:,Mesh%JMAX+j,Mesh%KMIN,:))
        END DO
#ifdef PARALLEL
      END IF
#endif

      ! top-south
#ifdef PARALLEL
      IF(Mesh%mycoords(2).EQ.0.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
        DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
          DO k=j+1,Mesh%GNUM
            pvar(:,Mesh%JMIN-j,Mesh%KMAX+k,:) = pvar(:,Mesh%JMIN-j+1,Mesh%KMAX+k,:)
            pvar(:,Mesh%JMIN-k,Mesh%KMAX+j,:) = pvar(:,Mesh%JMIN-k,Mesh%KMAX+j-1,:)
          END DO
          pvar(:,Mesh%JMIN-j,Mesh%KMAX+j,:) = 0.5 * (pvar(:,Mesh%JMIN,Mesh%KMAX+j,:) &
               + pvar(:,Mesh%JMIN-j,Mesh%KMAX,:))
        END DO
#ifdef PARALLEL
      END IF
#endif

      ! top-north
#ifdef PARALLEL
      IF(Mesh%mycoords(2).EQ.Mesh%dims(2)-1.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
      DO j=1,Mesh%GNUM
!NEC$ SHORTLOOP
        DO k=j+1,Mesh%GNUM
            pvar(:,Mesh%JMAX+j,Mesh%KMAX+k,:) = pvar(:,Mesh%JMAX+j-1,Mesh%KMAX+k,:)
            pvar(:,Mesh%JMAX+k,Mesh%KMAX+j,:) = pvar(:,Mesh%JMAX+k,Mesh%KMAX+j-1,:)
        END DO
          pvar(:,Mesh%JMAX+j,Mesh%KMAX+j,:) = 0.5 * (pvar(:,Mesh%JMAX,Mesh%KMAX+j,:) &
               + pvar(:,Mesh%JMAX+j,Mesh%KMAX,:))
      END DO
#ifdef PARALLEL
      END IF
#endif
    END IF

    ! Set corner cells (only 3D)
    IF (Mesh%INUM.NE.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.NE.1) THEN
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(2).EQ.0.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMIN-i,Mesh%JMIN-i,Mesh%KMIN-i,:) = (pvar(Mesh%IMIN-i,Mesh%JMIN,Mesh%KMIN,:) &
          + pvar(Mesh%IMIN,Mesh%JMIN-i,Mesh%KMIN,:) + pvar(Mesh%IMIN,Mesh%JMIN,Mesh%KMIN-i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(2).EQ.0.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMAX+i,Mesh%JMIN-i,Mesh%KMIN-i,:) = (pvar(Mesh%IMAX+i,Mesh%JMIN,Mesh%KMIN,:) &
          + pvar(Mesh%IMAX,Mesh%JMIN-i,Mesh%KMIN,:) + pvar(Mesh%IMAX,Mesh%JMIN,Mesh%KMIN-i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(2).EQ.Mesh%dims(2)-1.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMIN-i,Mesh%JMAX+i,Mesh%KMIN-i,:) = (pvar(Mesh%IMIN-i,Mesh%JMAX,Mesh%KMIN,:) &
          + pvar(Mesh%IMIN,Mesh%JMAX+i,Mesh%KMIN,:) + pvar(Mesh%IMIN,Mesh%JMAX,Mesh%KMIN-i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(2).EQ.Mesh%dims(2)-1.AND.Mesh%mycoords(3).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMAX+i,Mesh%JMAX+i,Mesh%KMIN-i,:) = (pvar(Mesh%IMAX+i,Mesh%JMAX,Mesh%KMIN,:) &
          + pvar(Mesh%IMAX,Mesh%JMAX+i,Mesh%KMIN,:) + pvar(Mesh%IMAX,Mesh%JMAX,Mesh%KMIN-i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(2).EQ.0.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMIN-i,Mesh%JMIN-i,Mesh%KMAX+i,:) = (pvar(Mesh%IMIN-i,Mesh%JMIN,Mesh%KMAX,:) &
          + pvar(Mesh%IMIN,Mesh%JMIN-i,Mesh%KMAX,:) + pvar(Mesh%IMIN,Mesh%JMIN,Mesh%KMAX+i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(2).EQ.0.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMAX+i,Mesh%JMIN-i,Mesh%KMAX+i,:) = (pvar(Mesh%IMAX+i,Mesh%JMIN,Mesh%KMAX,:) &
          + pvar(Mesh%IMAX,Mesh%JMIN-i,Mesh%KMAX,:) + pvar(Mesh%IMAX,Mesh%JMIN,Mesh%KMAX+i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0.AND.Mesh%mycoords(2).EQ.Mesh%dims(2)-1.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMIN-i,Mesh%JMAX+i,Mesh%KMAX+i,:) = (pvar(Mesh%IMIN-i,Mesh%JMAX,Mesh%KMAX,:) &
          + pvar(Mesh%IMIN,Mesh%JMAX+i,Mesh%KMAX,:) + pvar(Mesh%IMIN,Mesh%JMAX,Mesh%KMAX+i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1.AND.Mesh%mycoords(2).EQ.Mesh%dims(2)-1.AND.Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
#endif
!NEC$ SHORTLOOP
      DO i=1,Mesh%GNUM
        pvar(Mesh%IMAX+i,Mesh%JMAX+i,Mesh%KMAX+i,:) = (pvar(Mesh%IMAX+i,Mesh%JMAX,Mesh%KMAX,:) &
          + pvar(Mesh%IMAX,Mesh%JMAX+i,Mesh%KMAX,:) + pvar(Mesh%IMAX,Mesh%JMAX,Mesh%KMAX+i,:))/3.0
      END DO
#ifdef PARALLEL
      END IF
#endif
    END IF

  END SUBROUTINE


#ifdef PARALLEL
  !> \public initializes the MPI communication
  SUBROUTINE InitBoundary_MPI(this,Mesh,Physics,periods)
#ifdef HAVE_MPI_MOD
    USE mpi
#endif
    IMPLICIT NONE
#ifdef HAVE_MPIF_H
    include 'mpif.h'
#endif
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(INOUT) :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    LOGICAL, DIMENSION(3),  INTENT(IN)    :: periods
    !------------------------------------------------------------------------!
    INTEGER               :: comm_old
    INTEGER               :: ierr, dir
    LOGICAL, DIMENSION(SIZE(Mesh%dims)) :: remain_dims = .FALSE.
    !------------------------------------------------------------------------!
    ! create new cartesian communicator using Mesh%comm_cart
    ! and account for the periodicity
    ! IMPORTANT: disable reordering of nodes
    comm_old = Mesh%comm_cart
    CALL MPI_Cart_create(comm_old,SIZE(Mesh%dims),Mesh%dims,periods,.FALSE.,Mesh%comm_cart,ierr)

    ! save ranks of neighbor processes
    CALL MPI_Cart_shift(Mesh%comm_cart,0,1,Mesh%neighbor(WEST),Mesh%neighbor(EAST),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,1,1,Mesh%neighbor(SOUTH),Mesh%neighbor(NORTH),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,2,1,Mesh%neighbor(BOTTOM),Mesh%neighbor(TOP),ierr)

    ! create communicators for every column and row of the cartesian
    !	topology (used eg. for fargo shifts)
    remain_dims = (/ .FALSE., .TRUE., .TRUE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Icomm,ierr)
    remain_dims = (/ .TRUE., .FALSE., .TRUE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Jcomm,ierr)
    remain_dims = (/ .TRUE., .TRUE., .FALSE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Kcomm,ierr)

    ! allocate memory for boundary data buffers
    ALLOCATE(&
      this%boundary(WEST)%p%sendbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),   &
      this%boundary(WEST)%p%recvbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),   &
      this%boundary(EAST)%p%sendbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),   &
      this%boundary(EAST)%p%recvbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),   &
      this%boundary(SOUTH)%p%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),  &
      this%boundary(SOUTH)%p%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),  &
      this%boundary(NORTH)%p%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),  &
      this%boundary(NORTH)%p%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),  &
      this%boundary(BOTTOM)%p%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
      this%boundary(BOTTOM)%p%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
      this%boundary(TOP)%p%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM),    &
      this%boundary(TOP)%p%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM),    &
      STAT=ierr)
    IF (ierr.NE.0) THEN
       CALL this%boundary(WEST)%p%Error("boundary_generic::InitBoundary_MPI", &
            "Unable to allocate memory for data buffers.")
    END IF

    ! initialize all buffers with 0
    DO dir=WEST,TOP
      this%boundary(dir)%p%recvbuf = 0.
      this%boundary(dir)%p%sendbuf = 0.
    END DO
  END SUBROUTINE InitBoundary_MPI

  !> Handles the MPI communication of the inner (and physical periodic) boundaries
  SUBROUTINE MPIBoundaryCommunication(this,Mesh,Physics,pvar)
#ifdef HAVE_MPI_MOD
    USE mpi
#endif
    IMPLICIT NONE
#ifdef HAVE_MPIF_H
    include 'mpif.h'
#endif
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: ierr
    INTEGER                               :: req(4)
#ifdef MPI_USE_SENDRECV
    INTEGER                               :: status(MPI_STATUS_SIZE)
#else
    INTEGER                               :: status(MPI_STATUS_SIZE,4)
#endif
    !------------------------------------------------------------------------!
    ! NOTE: There are two different implementations of the MPI non-blocking
    ! communication, one uses MPI_Sendrecv (default) the other MPI_Irecv and
    ! MPI_Issend; to switch to the latter one has to remove -DMPI_USE_SENDRECV
    ! from the compile command which can be done by setting -DMPI_USE_SENDRECV=OFF
    ! when invoking cmake

    ! western <-> eastern MPI communication
#ifndef MPI_USE_SENDRECV
    ! receive boundary data from eastern neighbor (non-blocking)
    CALL MPI_Irecv(this%Boundary(EAST)%p%recvbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+WEST,Mesh%comm_cart,req(1),ierr)
#endif
    ! fill send buffer if western neighbor exists
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
      CALL this%MPIpvar2buffer(Mesh,Physics,WEST,pvar,this%Boundary(WEST)%p%sendbuf)
#ifdef MPI_USE_SENDRECV
    ! send boundary data to western and receive from eastern neighbor
    CALL MPI_Sendrecv(this%Boundary(WEST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,this%Boundary(EAST)%p%recvbuf,  &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
#else
    ! send boundary data to western neighbor
    CALL MPI_Issend(this%Boundary(WEST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,Mesh%comm_cart,req(2),ierr)

    ! receive boundary data from western neighbor
    CALL MPI_Irecv(this%Boundary(WEST)%p%recvbuf, &
      Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+EAST,Mesh%comm_cart,req(3),ierr)
#endif
    ! fill send buffer if eastern neighbor exists
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
      CALL this%MPIpvar2buffer(Mesh,Physics,EAST,pvar,this%Boundary(EAST)%p%sendbuf)
#ifdef MPI_USE_SENDRECV
    ! send boundary data to eastern and receive from western neighbor
    CALL MPI_Sendrecv(this%Boundary(EAST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,this%Boundary(WEST)%p%recvbuf,  &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
#else
    ! send boundary data to eastern neighbor
    CALL MPI_Issend(this%Boundary(EAST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,Mesh%comm_cart,req(4),ierr)
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
#endif
    ! copy data from receive buffers into ghosts cells
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
      CALL this%MPIbuffer2pvar(Mesh,Physics,WEST,this%Boundary(WEST)%p%recvbuf,pvar)
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
      CALL this%MPIbuffer2pvar(Mesh,Physics,EAST,this%Boundary(EAST)%p%recvbuf,pvar)

    ! southern <-> northern MPI communication
#ifndef MPI_USE_SENDRECV
    ! receive boundary data from northern neighbor
    CALL MPI_Irecv(this%Boundary(NORTH)%p%recvbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+SOUTH,Mesh%comm_cart,req(1),ierr)
#endif
    ! fill send buffer if southern neighbor exists
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
      CALL this%MPIpvar2buffer(Mesh,Physics,SOUTH,pvar,this%Boundary(SOUTH)%p%sendbuf)
#ifdef MPI_USE_SENDRECV
    ! send boundary data to southern and receive from northern neighbor
    CALL MPI_Sendrecv(this%Boundary(SOUTH)%p%sendbuf, &
        Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
        DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,this%Boundary(NORTH)%p%recvbuf,          &
        Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
        DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
#else
    ! send boundary data to southern neighbor
    CALL MPI_Issend(this%Boundary(SOUTH)%p%sendbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,Mesh%comm_cart,req(2),ierr)

    ! receive boundary data from southern neighbor
    CALL MPI_Irecv(this%Boundary(SOUTH)%p%recvbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+NORTH,Mesh%comm_cart,req(3),ierr)
#endif
    ! fill send buffer if northern neighbor exists
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
      CALL this%MPIpvar2buffer(Mesh,Physics,NORTH,pvar,this%Boundary(NORTH)%p%sendbuf)
#ifdef MPI_USE_SENDRECV
    ! send boundary data to northern and receive from southern neighbor
    CALL MPI_Sendrecv(this%Boundary(NORTH)%p%sendbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,this%Boundary(SOUTH)%p%recvbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
#else
    ! send boundary data to northern neighbor
    CALL MPI_Issend(this%Boundary(NORTH)%p%sendbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,Mesh%comm_cart,req(4),ierr)
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
#endif
    ! copy data from receive buffers into ghosts cells
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
      CALL this%MPIbuffer2pvar(Mesh,Physics,SOUTH,this%Boundary(SOUTH)%p%recvbuf,pvar)
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
      CALL this%MPIbuffer2pvar(Mesh,Physics,NORTH,this%Boundary(NORTH)%p%recvbuf,pvar)

    ! bottom <-> top MPI communication
#ifndef MPI_USE_SENDRECV
    ! receive boundary data from top neighbor
    CALL MPI_Irecv(this%Boundary(TOP)%p%recvbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),10+BOTTOM,Mesh%comm_cart,req(1),ierr)
#endif
    ! fill send buffer if bottomer neighbor exists
    IF (Mesh%neighbor(BOTTOM).NE.MPI_PROC_NULL) &
      CALL this%MPIpvar2buffer(Mesh,Physics,BOTTOM,pvar,this%Boundary(BOTTOM)%p%sendbuf)
#ifdef MPI_USE_SENDRECV
    ! send boundary data to bottom and receive from top neighbor
    CALL MPI_Sendrecv(this%Boundary(BOTTOM)%p%sendbuf,&
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),10+BOTTOM,this%Boundary(TOP)%p%recvbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
#else
    ! send boundary data to bottom neighbor
    CALL MPI_Issend(this%Boundary(BOTTOM)%p%sendbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),10+BOTTOM,Mesh%comm_cart,req(2),ierr)

    ! receive boundary data from bottom neighbor
    CALL MPI_Irecv(this%Boundary(BOTTOM)%p%recvbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),10+TOP,Mesh%comm_cart,req(3),ierr)
#endif
    ! send boundary data to top and receive from bottom neighbor
    IF (Mesh%neighbor(TOP).NE.MPI_PROC_NULL) &
      CALL this%MPIpvar2buffer(Mesh,Physics,TOP,pvar,this%Boundary(TOP)%p%sendbuf)
#ifdef MPI_USE_SENDRECV
    ! send boundary data to top and receive from bottom neighbor
    CALL MPI_Sendrecv(this%Boundary(TOP)%p%sendbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),10+TOP,this%Boundary(BOTTOM)%p%recvbuf,  &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,  &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
#else
    ! send boundary data to top neighbor
    CALL MPI_Issend(this%Boundary(TOP)%p%sendbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),10+TOP,Mesh%comm_cart,req(4),ierr)
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
#endif
    ! copy data from receive buffers into ghosts cells
    IF (Mesh%neighbor(BOTTOM).NE.MPI_PROC_NULL) &
      CALL this%MPIbuffer2pvar(Mesh,Physics,BOTTOM,this%Boundary(BOTTOM)%p%recvbuf,pvar)
    IF (Mesh%neighbor(TOP).NE.MPI_PROC_NULL) &
      CALL this%MPIbuffer2pvar(Mesh,Physics,TOP,this%Boundary(TOP)%p%recvbuf,pvar)
  END SUBROUTINE MPIBoundaryCommunication

  !> Copys buffer data to ghost cells
  SUBROUTINE MPIbuffer2pvar(this,Mesh,Physics,output_dir,buffer,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    INTEGER,                INTENT(IN)    :: output_dir
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: buffer
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                              INTENT(OUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,n
    !------------------------------------------------------------------------!
    SELECT CASE(output_dir)
    CASE(WEST)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
           DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ shortloop
             DO i=1,Mesh%GINUM
                pvar(Mesh%IGMIN+i-1,j,k,n) = buffer(i,j-Mesh%JGMIN+1,k-Mesh%KGMIN+1,n)
             END DO
           END DO
         END DO
       END DO
    CASE(EAST)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
           DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ shortloop
             DO i=1,Mesh%GINUM
                pvar(Mesh%IMAX+i,j,k,n) = buffer(i,j-Mesh%JGMIN+1,k-Mesh%KGMIN+1,n)
             END DO
           END DO
         END DO
       END DO
    CASE(SOUTH)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ shortloop
           DO j=1,Mesh%GJNUM
             DO i=Mesh%IGMIN,Mesh%IGMAX
                pvar(i,Mesh%JGMIN+j-1,k,n) = buffer(i-Mesh%IGMIN+1,j,k-Mesh%KGMIN+1,n)
             END DO
           END DO
         END DO
       END DO
    CASE(NORTH)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ shortloop
           DO j=1,Mesh%GJNUM
             DO i=Mesh%IGMIN,Mesh%IGMAX
                pvar(i,Mesh%JMAX+j,k,n) = buffer(i-Mesh%IGMIN+1,j,k-Mesh%KGMIN+1,n)
             END DO
           END DO
         END DO
       END DO
    CASE(BOTTOM)
!NEC$ shortloop
       DO n=1,Physics%VNUM
!NEC$ shortloop
         DO k=1,Mesh%GKNUM
           DO j=Mesh%JGMIN,Mesh%JGMAX
             DO i=Mesh%IGMIN,Mesh%IGMAX
                pvar(i,j,Mesh%KGMIN+k-1,n) = buffer(i-Mesh%IGMIN+1,j-Mesh%JGMIN+1,k,n)
             END DO
           END DO
         END DO
       END DO
    CASE(TOP)
!NEC$ shortloop
       DO n=1,Physics%VNUM
!NEC$ shortloop
         DO k=1,Mesh%GKNUM
           DO j=Mesh%JGMIN,Mesh%JGMAX
             DO i=Mesh%IGMIN,Mesh%IGMAX
                pvar(i,j,Mesh%KMAX+k,n) = buffer(i-Mesh%IGMIN+1,j-Mesh%JGMIN+1,k,n)
             END DO
           END DO
         END DO
       END DO
    END SELECT
  END SUBROUTINE MPIbuffer2pvar

  !> Copys ghost cell data to buffer
  SUBROUTINE MPIpvar2buffer(this,Mesh,Physics,input_dir,pvar,buffer)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    INTEGER,                INTENT(IN)    :: input_dir
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                              INTENT(IN)  :: pvar
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: buffer
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,n
    !------------------------------------------------------------------------!
    SELECT CASE(input_dir)
    CASE(WEST)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
           DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ shortloop
             DO i=1,Mesh%GINUM
                buffer(i,j-Mesh%JGMIN+1,k-Mesh%KGMIN+1,n) = pvar(Mesh%IMIN+i-1,j,k,n)
             END DO
           END DO
         END DO
       END DO
    CASE(EAST)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
           DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ shortloop
             DO i=1,Mesh%GINUM
                buffer(i,j-Mesh%JGMIN+1,k-Mesh%KGMIN+1,n) = pvar(Mesh%IMAX-Mesh%GINUM+i,j,k,n)
             END DO
           END DO
         END DO
       END DO
    CASE(SOUTH)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ shortloop
           DO j=1,Mesh%GJNUM
             DO i=Mesh%IGMIN,Mesh%IGMAX
               buffer(i-Mesh%IGMIN+1,j,k-Mesh%KGMIN+1,n) = pvar(i,Mesh%JMIN+j-1,k,n)
             END DO
           END DO
         END DO
       END DO
    CASE(NORTH)
!NEC$ shortloop
       DO n=1,Physics%VNUM
         DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ shortloop
           DO j=1,Mesh%GJNUM
             DO i=Mesh%IGMIN,Mesh%IGMAX
               buffer(i-Mesh%IGMIN+1,j,k-Mesh%KGMIN+1,n) = pvar(i,Mesh%JMAX-Mesh%GJNUM+j,k,n)
             END DO
           END DO
         END DO
       END DO
    CASE(BOTTOM)
!NEC$ shortloop
       DO n=1,Physics%VNUM
!NEC$ shortloop
         DO k=1,Mesh%GKNUM
           DO j=Mesh%JGMIN,Mesh%JGMAX
             DO i=Mesh%IGMIN,Mesh%IGMAX
               buffer(i-Mesh%IGMIN+1,j-Mesh%JGMIN+1,k,n) = pvar(i,j,Mesh%KMIN+k-1,n)
             END DO
           END DO
         END DO
       END DO
    CASE(TOP)
!NEC$ shortloop
       DO n=1,Physics%VNUM
!NEC$ shortloop
         DO k=1,Mesh%GKNUM
           DO j=Mesh%JGMIN,Mesh%JGMAX
             DO i=Mesh%IGMIN,Mesh%IGMAX
               buffer(i-Mesh%IGMIN+1,j-Mesh%JGMIN+1,k,n) = pvar(i,j,Mesh%KMAX-Mesh%GKNUM+k,n)
             END DO
           END DO
         END DO
       END DO
    END SELECT
  END SUBROUTINE MPIpvar2buffer
#endif

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    INTEGER                               :: dir
    !------------------------------------------------------------------------!
    ! loop over all boundaries
    DO dir=1,6
      CALL this%boundary(dir)%p%Finalize()
#ifdef PARALLEL
      ! deallocate MPI send/recv buffers
      DEALLOCATE(this%boundary(dir)%p%sendbuf,this%boundary(dir)%p%recvbuf)
#endif
    END DO
  END SUBROUTINE Finalize

END MODULE boundary_generic_mod
