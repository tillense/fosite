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
  USE physics_base_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  TYPE, PRIVATE                       :: boundary_p
    CLASS(boundary_base), ALLOCATABLE :: p
#ifdef PARALLEL
    !> \name variables in parallel mode
    REAL, DIMENSION(:,:,:,:), POINTER :: sendbuf, &     !< send buffer for boundary data
                                         recvbuf        !< receive buffer for boundary data
#endif
  END TYPE

  TYPE, EXTENDS(logging_base)         :: boundary_generic
    !> \name variables
    TYPE(boundary_p)                  :: Boundary(6)
    LOGICAL                           :: PhysicalCorner  !< Is the left corner physical?
#ifdef PARALLEL
    !> \name variables in parallel mode
    REAL, DIMENSION(:,:,:,:), POINTER :: sendbuf, &     !< send buffer for boundary data
                                         recvbuf        !< receive buffer for boundary data
#endif
  CONTAINS
    PROCEDURE :: InitBoundary
    PROCEDURE :: CenterBoundary
    PROCEDURE :: Finalize
  END TYPE boundary_generic
  !--------------------------------------------------------------------------!
  !> \name Public Attributes
  CHARACTER(LEN=32), DIMENSION(6), PARAMETER :: &
  direction_name = (/'  west', '  east', ' south', ' north', 'bottom', '   top' /)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       boundary_generic, new_boundary
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
    LOGICAL, DIMENSION(3) :: periods = .FALSE.
    INTEGER               :: dir
#ifdef PARALLEL
    INTEGER               :: comm_old
    INTEGER               :: ierr
    LOGICAL, DIMENSION(SIZE(Mesh%dims)) :: remain_dims = .FALSE.
#endif
    !------------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("InitBoundary","physics and/or mesh module uninitialized")

    CALL GetAttr(config, "western",   western)
    CALL GetAttr(config, "eastern",   eastern)
    CALL GetAttr(config, "southern", southern)
    CALL GetAttr(config, "northern", northern)
    CALL GetAttr(config, "bottomer", bottomer)
    CALL GetAttr(config, "topper",     topper)

    new(WEST)   = western
    new(EAST)   = eastern
    new(SOUTH)  = southern
    new(NORTH)  = northern
    new(BOTTOM) = bottomer
    new(TOP)    = topper

#ifdef PARALLEL
    ! define inner connections, where boundaries are no true ones
    IF (Mesh%mycoords(1).NE.0)               new(WEST) = NONE
    IF (Mesh%mycoords(1).NE.Mesh%dims(1)-1)  new(EAST) = NONE
    IF (Mesh%mycoords(2).NE.0)              new(SOUTH) = NONE
    IF (Mesh%mycoords(2).NE.Mesh%dims(2)-1) new(NORTH) = NONE
    IF (Mesh%mycoords(3).NE.0)             new(BOTTOM) = NONE
    IF (Mesh%mycoords(3).NE.Mesh%dims(3)-1)   new(TOP) = NONE
#endif

    ! Check for correct shifting and boundaries
    IF (Mesh%FARGO.EQ.3.AND.western.EQ.SHEARING.AND.eastern.EQ.SHEARING) THEN
      IF (.NOT.Mesh%WE_shear) &
        CALL this%Error("InitBoundary", &
          "Please apply shifting in second dimension, when applying shearing boundaries at western/eastern.")
#ifdef PARALLEL
      CALL this%Error("InitBoundary", &
        "Parallel mode is not allowed with shearing in West-East direction.")
#endif
    ELSE IF (Mesh%FARGO.EQ.3.AND.southern.EQ.SHEARING.AND.northern.EQ.SHEARING) THEN
      IF (.NOT.Mesh%SN_shear) &
        CALL this%Error("InitBoundary", &
          "Please apply shifting in first dimension, when applying shearing boundaries at northern/southern.")
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
      END SELECT

      SELECT TYPE(obj => this%Boundary(dir)%p)
      TYPE IS (boundary_absorbing)
        CALL obj%InitBoundary_absorbing(Mesh,Physics,dir,config)
      TYPE IS (boundary_axis)
        CALL obj%InitBoundary_axis(Mesh,Physics,dir,config)
      TYPE IS (boundary_custom)
        CALL obj%InitBoundary_custom(Mesh,Physics,dir,config)
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

    ! This sets this%Boundary(dir)%PhysicalCorner to .True., if
    ! there are no periodic boundary conditions and if
    ! is a corner without inner boundaries.
    ! At physical corners the corner value have to be
    ! interpolated in CenterBoundary.
    ! e.g. this(1)%PhysicalCorner describes the Southwest corner.
    !!
    ! \todo BUG: THIS IS WRONG, WE HAVE NOW 8 CORNERS IN 3D
!    DO dir=1,6
!      this%Boundary(dir)%p%PhysicalCorner = .FALSE.
!    END DO
!    ! TODO Help needed for parallelisation
!    ! TODO TODO TODO TODO TODO TODO TODO TODO TODO
!    IF(.NOT.ANY(periods)) THEN
!#ifdef PARALLEL
!      IF(Mesh%mycoords(1).EQ.0) THEN
!        IF(Mesh%mycoords(2).EQ.0) THEN
!          IF(Mesh%mycoords(3).EQ.0) THEN
!            this%Boundary(1)%p%PhysicalCorner = .TRUE.
!          END IF
!          IF(Mesh%mycoords(3).EQ.Mesh%dims(3)-1)
!            this%Boundary(4)%PhysicalCorner = .TRUE.
!      END IF
!        IF(Mesh%mycoords(2).EQ.Mesh%dims(2)-1)
!          this%Boundary(4)%PhysicalCorner = .TRUE.
!      END IF
!      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1) THEN
!        IF(Mesh%mycoords(2).EQ.0) &
!          this%Boundary(3)%PhysicalCorner = .TRUE.
!        IF(Mesh%mycoords(2).EQ.Mesh%dims(2)-1) &
!          this%Boundary(2)%PhysicalCorner = .TRUE.
!      END IF
!#else
!      DO dir=1,6
!        this%Boundary(dir)%p%PhysicalCorner = .TRUE.
!      END DO
!#endif
!    END IF


#ifdef PARALLEL
    ! create new cartesian communicator using Mesh%comm_cart
    ! and account for the periodicity
    ! IMPORTANT: disable reordering of nodes
    comm_old = Mesh%comm_cart
    CALL MPI_Cart_create(comm_old,SIZE(Mesh%dims),Mesh%dims,periods,.FALSE.,Mesh%comm_cart,ierr)

    ! save ranks of neighbor processes
    CALL MPI_Cart_shift(Mesh%comm_cart,0,1,Mesh%neighbor(WEST),Mesh%neighbor(EAST),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,1,1,Mesh%neighbor(SOUTH),Mesh%neighbor(NORTH),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,2,1,Mesh%neighbor(TOP),Mesh%neighbor(BOTTOM),ierr)

    ! create communicators for every column and row of the cartesian
    !	topology (used eg. for fargo shifts)
    remain_dims = (/ .FALSE., .TRUE., .TRUE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Icomm,ierr)
    remain_dims = (/ .TRUE., .FALSE., .TRUE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Jcomm,ierr)
    remain_dims = (/ .TRUE., .TRUE., .FALSE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Kcomm,ierr)

    ! allocate memory for boundary data buffers
    ALLOCATE(                                                                                                  &
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
       CALL this%boundary(WEST)%p%Error("InitBoundary", &
            "Unable to allocate memory for data buffers.")
    END IF
    DO dir=WEST,TOP
      this%boundary(dir)%p%recvbuf = 0.
      this%boundary(dir)%p%sendbuf = 0.
    END DO
#endif

  END SUBROUTINE InitBoundary


  SUBROUTINE CenterBoundary(this,Mesh,Physics,time,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    REAL,                   INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(INOUT) :: pvar, cvar
    !------------------------------------------------------------------------!
    INTEGER                               :: i,j,k
#ifdef PARALLEL
    INTEGER                               :: req(4)
    INTEGER                               :: ierr
!    REAL    :: mpi_buf2(Mesh%IGMIN:Mesh%IGMAX,1:Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!    REAL    :: mpi_buf2(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) !TODO: ONLY 2D
!    INTEGER :: status2(MPI_STATUS_SIZE)
#ifdef MPI_USE_SENDRECV
    INTEGER                               :: status(MPI_STATUS_SIZE)
#else
    INTEGER                               :: status(MPI_STATUS_SIZE,4)
#endif
#endif
    !------------------------------------------------------------------------!
    CALL Physics%Convert2Primitive(Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN, &
             Mesh%JMAX,Mesh%KMIN,Mesh%KMAX,cvar,pvar)

    ! set physical boundary conditions at western and eastern boundaries
    IF (Mesh%INUM.GT.1) THEN
      CALL this%Boundary(WEST)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      CALL this%Boundary(EAST)%p%SetBoundaryData(Mesh,Physics,time,pvar)
    END IF

!    IF (Mesh%SN_shear) THEN
!      IF(Mesh%dims(2).GT.1) THEN
!        mpi_buf2(Mesh%IGMIN:Mesh%IGMAX,1:Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
!          pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX-Mesh%GJNUM+1:Mesh%JMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!        CALL MPI_Sendrecv_replace(&
!          mpi_buf2,&
!          2*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM, &
!          DEFAULT_MPI_REAL, &
!          Mesh%neighbor(NORTH), 53+NORTH, &
!          Mesh%neighbor(SOUTH), MPI_ANY_TAG, &
!          Mesh%comm_cart, status2, ierr)
!        pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
!          mpi_buf2(Mesh%IGMIN:Mesh%IGMAX,1:Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!        mpi_buf2(Mesh%IGMIN:Mesh%IGMAX,1:Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
!          pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN:Mesh%JMIN+Mesh%GJNUM-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!        CALL MPI_Sendrecv_replace(&
!          mpi_buf2,&
!          2*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM, &
!          DEFAULT_MPI_REAL, &
!          Mesh%neighbor(SOUTH), 53+SOUTH, &
!          Mesh%neighbor(NORTH), MPI_ANY_TAG, &
!          Mesh%comm_cart, status2, ierr)
!        pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
!          mpi_buf2(Mesh%IGMIN:Mesh%IGMAX,1:Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!      ELSE
!        DO j = 1,Mesh%GJNUM
!          ! southern northern (periodic in first step - further shift-treatment below)
!          pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN-j,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
!            pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX-j+1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!          pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+j,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
!            pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN+j-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!        END DO
!      END IF
!    END IF


    ! set physical boundary conditions at top and bottom boundaries
    IF (Mesh%KNUM.GT.1) THEN
      CALL this%Boundary(BOTTOM)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      CALL this%Boundary(TOP)%p%SetBoundaryData(Mesh,Physics,time,pvar)
    END IF


!> \todo not verified
#ifdef PARALLEL
    ! NOTE: if you want to use MPI_Sendrecv instead of nonblocking
    ! MPI_Irecv and  MPI_Issend for exchange of ghost cell data,
    ! you must add -DMPI_USE_SENDRECV to the compile command

    ! initiate western/eastern MPI communication
#ifdef MPI_USE_SENDRECV
    ! send boundary data to western and receive from eastern neighbor
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) THEN
      this%Boundary(WEST)%p%sendbuf(:,:,:,:) = pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GINUM-1, &
        Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    CALL MPI_Sendrecv(this%Boundary(WEST)%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,this%Boundary(EAST)%p%recvbuf,  &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) THEN
      pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
        this%Boundary(EAST)%p%recvbuf(:,:,:,:)
    END IF
#else
    ! receive boundary data from eastern neighbor
    CALL MPI_Irecv(this%Boundary(EAST)%p%recvbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+WEST,Mesh%comm_cart,req(1),ierr)
    ! fill send buffer if western neighbor exists
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) THEN
      this%Boundary(WEST)%p%sendbuf(:,:,:,:) = &
        pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GINUM-1,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    ! send boundary data to western neighbor
    CALL MPI_Issend(this%Boundary(WEST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,Mesh%comm_cart,req(2),ierr)
#endif
#endif
#ifdef PARALLEL
#ifdef MPI_USE_SENDRECV
    ! send boundary data to eastern and receive from western neighbor
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) THEN
          this%Boundary(EAST)%p%sendbuf(:,:,:,:) = pvar(Mesh%IMAX-Mesh%GINUM+1:Mesh%IMAX, &
             Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    CALL MPI_Sendrecv(this%Boundary(EAST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,this%Boundary(WEST)%p%recvbuf,  &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) THEN
      pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
              this%Boundary(WEST)%p%recvbuf(:,:,:,:)
    END IF
#else
    ! receive boundary data from western neighbor
    CALL MPI_Irecv(this%Boundary(WEST)%p%recvbuf, &
      Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+EAST,Mesh%comm_cart,req(3),ierr)
    ! fill send buffer if eastern neighbor exists
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) THEN
      this%Boundary(EAST)%p%sendbuf(:,:,:,:) = &
         pvar(Mesh%IMAX-Mesh%GINUM+1:Mesh%IMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    ! send boundary data to eastern neighbor
    CALL MPI_Issend(this%Boundary(EAST)%p%sendbuf, &
         Mesh%GINUM*(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,Mesh%comm_cart,req(4),ierr)
#endif
#endif
#ifdef PARALLEL
#ifndef MPI_USE_SENDRECV
   ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
   ! copy data from recieve buffers into ghosts cells
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) THEN
        pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
          this%Boundary(WEST)%p%recvbuf(:,:,:,:)
    END IF
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) THEN
        pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
          this%Boundary(EAST)%p%recvbuf(:,:,:,:)
    END IF
#endif
#endif

#ifdef PARALLEL
   ! initiate southern/northern MPI communication
#ifdef MPI_USE_SENDRECV
   ! send boundary data to southern and receive from northern neighbor
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) THEN
        this%Boundary(SOUTH)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JMIN:Mesh%JMIN+Mesh%GJNUM-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    CALL MPI_Sendrecv(this%Boundary(SOUTH)%p%sendbuf, &
        Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
        DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,this%Boundary(NORTH)%recvbuf,          &
        Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
        DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) THEN
        pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
          this%Boundary(NORTH)%p%recvbuf(:,:,:,:)
    END IF
#else
   ! receive boundary data from northern neighbor
    CALL MPI_Irecv(this%Boundary(NORTH)%p%recvbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+SOUTH,Mesh%comm_cart,req(1),ierr)
    ! fill send buffer if southern neighbor exists
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) THEN
         this%Boundary(SOUTH)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JMIN:Mesh%JMIN+Mesh%GJNUM-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    ! send boundary data to southern neighbor
    CALL MPI_Issend(this%Boundary(SOUTH)%p%sendbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,Mesh%comm_cart,req(2),ierr)
#endif
#endif
#ifdef PARALLEL
#ifdef MPI_USE_SENDRECV
   ! send boundary data to northern and receive from southern neighbor
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) THEN
         this%Boundary(NORTH)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JMAX-Mesh%GJNUM+1:Mesh%JMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    CALL MPI_Sendrecv(this%Boundary(NORTH)%p%sendbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,this%Boundary(SOUTH)%p%recvbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) THEN
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
           this%Boundary(SOUTH)%p%recvbuf(:,:,:,:)
    END IF
#else
    ! receive boundary data from southern neighbor
    CALL MPI_Irecv(this%Boundary(SOUTH)%p%recvbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+NORTH,Mesh%comm_cart,req(3),ierr)
    ! fill send buffer if northern neighbor exists
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) THEN
         this%Boundary(NORTH)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JMAX-Mesh%GJNUM+1:Mesh%JMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    END IF
    ! send boundary data to northern neighbor
    CALL MPI_Issend(this%Boundary(NORTH)%p%sendbuf, &
         Mesh%GJNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,Mesh%comm_cart,req(4),ierr)
#endif
#endif
#ifdef PARALLEL
#ifndef MPI_USE_SENDRECV
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) THEN
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
         this%Boundary(SOUTH)%p%recvbuf(:,:,:,:)
    END IF
    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) THEN
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = &
           this%Boundary(NORTH)%p%recvbuf(:,:,:,:)
    END IF
#endif
#endif

    ! set physical boundary conditions at southern and northern boundaries
    IF (Mesh%JNUM.GT.1) THEN
      CALL this%Boundary(SOUTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
      CALL this%Boundary(NORTH)%p%SetBoundaryData(Mesh,Physics,time,pvar)
    END IF

#ifdef PARALLEL
   ! initiate bottomer/topper MPI communication
#ifdef MPI_USE_SENDRECV
   ! send boundary data to bottomer and receive from topper neighbor
    IF (Mesh%neighbor(BOTTOM).NE.MPI_PROC_NULL) THEN
        this%Boundary(BOTTOM)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMIN+Mesh%GKNUM-1,1:Physics%VNUM)
    END IF
    CALL MPI_Sendrecv(this%Boundary(BOTTOM)%p%sendbuf,&
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),10+BOTTOM,this%Boundary(TOP)%p%recvbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,   &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(TOP).NE.MPI_PROC_NULL) THEN
        pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMAX+1:Mesh%KGMAX,1:Physics%VNUM) = &
          this%Boundary(TOP)%p%recvbuf(:,:,:,:)
    END IF
#else
   ! receive boundary data from topper neighbor
    CALL MPI_Irecv(this%Boundary(TOP)%p%recvbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),10+BOTTOM,Mesh%comm_cart,req(1),ierr)
    ! fill send buffer if bottomer neighbor exists
    IF (Mesh%neighbor(BOTTOM).NE.MPI_PROC_NULL) THEN
         this%Boundary(BOTTOM)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMIN+Mesh%GKNUM-1,1:Physics%VNUM)
    END IF
    ! send boundary data to bottomer neighbor
    CALL MPI_Issend(this%Boundary(BOTTOM)%p%sendbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),10+BOTTOM,Mesh%comm_cart,req(2),ierr)
#endif
#endif
#ifdef PARALLEL
#ifdef MPI_USE_SENDRECV
   ! send boundary data to northern and receive from southern neighbor
    IF (Mesh%neighbor(TOP).NE.MPI_PROC_NULL) THEN
         this%Boundary(TOP)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JGMIN:Mesh%JGMAX,Mesh%KMAX-Mesh%GKNUM+1:Mesh%KMAX,1:Physics%VNUM)
    END IF
    CALL MPI_Sendrecv(this%Boundary(TOP)%p%sendbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),10+TOP,this%Boundary(BOTTOM)%p%recvbuf,  &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,  &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
    IF (Mesh%neighbor(BOTTOM).NE.MPI_PROC_NULL) THEN
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KMIN-1,1:Physics%VNUM) = &
           this%Boundary(BOTTOM)%p%recvbuf(:,:,:,:)
    END IF
#else
    ! receive boundary data from southern neighbor
    CALL MPI_Irecv(this%Boundary(BOTTOM)%p%recvbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(BOTTOM),10+TOP,Mesh%comm_cart,req(3),ierr)
    ! fill send buffer if northern neighbor exists
    IF (Mesh%neighbor(TOP).NE.MPI_PROC_NULL) THEN
         this%Boundary(TOP)%p%sendbuf(:,:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
           Mesh%JGMIN:Mesh%JGMAX,Mesh%KMAX-Mesh%GKNUM+1:Mesh%KMAX,1:Physics%VNUM)
    END IF
    ! send boundary data to northern neighbor
    CALL MPI_Issend(this%Boundary(TOP)%p%sendbuf, &
         Mesh%GKNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1)*Physics%VNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(TOP),10+TOP,Mesh%comm_cart,req(4),ierr)
#endif
#endif
#ifdef PARALLEL
#ifndef MPI_USE_SENDRECV
    ! wait for unfinished MPI communication
    CALL MPI_Waitall(4,req,status,ierr)
    IF (Mesh%neighbor(BOTTOM).NE.MPI_PROC_NULL) THEN
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KMIN-1,1:Physics%VNUM) = &
         this%Boundary(BOTTOM)%p%recvbuf(:,:,:,:)
    END IF
    IF (Mesh%neighbor(TOP).NE.MPI_PROC_NULL) THEN
         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMAX+1:Mesh%KGMAX,1:Physics%VNUM) = &
           this%Boundary(TOP)%p%recvbuf(:,:,:,:)
    END IF
#endif
#endif

    ! \todo This module should be not necessary anymore, since the
    !       remaining region vanishes.
    ! This is a interpolation of corners outside
    ! the computational domain, if they are undefined (e.g. there are
    ! no periodic or inner boundaries involved in the corner)
    ! this is also necessary, because we need some of these values in the
    ! viscosity module
    ! For 1D we only copy the values:
    IF (Mesh%JNUM.EQ.1 .AND. Mesh%KNUM.EQ.1) THEN
      ! for 1D simulations along the x-direction copy internal data
      ! into the ghost cells at the yz-faces
      DO k=1,Mesh%GKNUM
        DO j=1,Mesh%GJNUM
          pvar(:,Mesh%JMIN-j,Mesh%KMIN-k,:) = pvar(:,Mesh%JMIN,Mesh%KMIN,:)
          pvar(:,Mesh%JMAX+j,Mesh%KMAX+k,:) = pvar(:,Mesh%JMAX,Mesh%KMAX,:)
        END DO
      END DO
    ELSE IF (Mesh%INUM.EQ.1 .AND. Mesh%KNUM.EQ.1) THEN
      ! for 1D simulations along the y-direction copy internal data
      ! into the ghost cells at the xz-faces
      DO k=1,Mesh%GKNUM
        DO i=1,Mesh%GINUM
          pvar(Mesh%IMIN-i,:,Mesh%KMIN-k,:) = pvar(Mesh%IMIN,:,Mesh%KMIN,:)
          pvar(Mesh%IMAX+i,:,Mesh%KMAX+k,:) = pvar(Mesh%IMAX,:,Mesh%KMAX,:)
        END DO
      END DO
    ELSE IF (Mesh%INUM.EQ.1 .AND. Mesh%JNUM.EQ.1) THEN
      ! for 1D simulations along the z-direction copy internal data
      ! into the ghost cells at the xy-faces
      DO j=1,Mesh%GJNUM
        DO i=1,Mesh%GINUM
          pvar(Mesh%IMIN-i,Mesh%JMIN-j,:,:) = pvar(Mesh%IMIN,Mesh%JMIN,:,:)
          pvar(Mesh%IMAX+i,Mesh%JMAX+j,:,:) = pvar(Mesh%IMAX,Mesh%JMAX,:,:)
        END DO
      END DO
    ELSE IF (Mesh%KNUM.EQ.1.AND.Mesh%INUM.NE.1.AND.Mesh%JNUM.NE.1) THEN ! 2D case
        ! south west
        IF(this%Boundary(1)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO i=j+1,Mesh%GINUM
              ! copy data from adjacent ghost cells for off diagonal corner ghost cells
              pvar(Mesh%IMIN-i,Mesh%JMIN-j,:,:) = pvar(Mesh%IMIN-i,Mesh%JMIN-j+1,:,:)
              pvar(Mesh%IMIN-j,Mesh%JMIN-i,:,:) = pvar(Mesh%IMIN-j+1,Mesh%JMIN-i,:,:)
            END DO
            ! interpolate diagonal corner ghost cell values from adjacent cells in
          ! both directions
            pvar(Mesh%IMIN-j,Mesh%JMIN-j,:,:) = 0.5 * (pvar(Mesh%IMIN-j,Mesh%JMIN,:,:) &
                 + pvar(Mesh%IMIN,Mesh%JMIN-j,:,:))
          END DO
        END IF

        ! north west
        IF(this%Boundary(4)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO i=j+1,Mesh%GINUM
              pvar(Mesh%IMIN-i,Mesh%JMAX+j,:,:) = pvar(Mesh%IMIN-i,Mesh%JMAX+j-1,:,:)
              pvar(Mesh%IMIN-j,Mesh%JMAX+i,:,:) = pvar(Mesh%IMIN-j+1,Mesh%JMAX+i,:,:)
            END DO
            pvar(Mesh%IMIN-j,Mesh%JMAX+j,:,:) = 0.5 * (pvar(Mesh%IMIN-j,Mesh%JMAX,:,:) &
                 + pvar(Mesh%IMIN,Mesh%JMAX+j,:,:))
          END DO
        END IF

        ! north east
        IF(this%Boundary(2)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO i=j+1,Mesh%GINUM
              pvar(Mesh%IMAX+i,Mesh%JMAX+j,:,:) = pvar(Mesh%IMAX+i,Mesh%JMAX+j-1,:,:)
              pvar(Mesh%IMAX+j,Mesh%JMAX+i,:,:) = pvar(Mesh%IMAX+j-1,Mesh%JMAX+i,:,:)
            END DO
            pvar(Mesh%IMAX+j,Mesh%JMAX+j,:,:) = 0.5 * (pvar(Mesh%IMAX+j,Mesh%JMAX,:,:) &
                 + pvar(Mesh%IMAX,Mesh%JMAX+j,:,:))
          END DO
        END IF

        ! south east
        IF(this%Boundary(3)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO i=j+1,Mesh%GINUM
              pvar(Mesh%IMAX+i,Mesh%JMIN-j,:,:) = pvar(Mesh%IMAX+i,Mesh%JMIN-j+1,:,:)
              pvar(Mesh%IMAX+j,Mesh%JMIN-i,:,:) = pvar(Mesh%IMAX+j-1,Mesh%JMIN-i,:,:)
            END DO
            pvar(Mesh%IMAX+j,Mesh%JMIN-j,:,:) = 0.5 * (pvar(Mesh%IMAX+j,Mesh%JMIN,:,:) &
                 + pvar(Mesh%IMAX,Mesh%JMIN-j,:,:))
          END DO
        END IF
      ELSE IF (Mesh%INUM.EQ.1.AND.Mesh%JNUM.NE.1.AND.Mesh%KNUM.NE.1) THEN
        !TODO: Hier hatte Jannes JNUM durch GNUM ersetzt...warum? evtl GJNUM
        !gemeint?
         ! bottom south
        IF(this%Boundary(1)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO k=j+1,Mesh%GKNUM
              ! copy data from adjacent ghost cells for off diagonal corner ghost cells
              pvar(:,Mesh%JMIN-j,Mesh%KMIN-k,:) = pvar(:,Mesh%JMIN-j+1,Mesh%KMIN-k,:)
              pvar(:,Mesh%JMIN-k,Mesh%KMIN-j,:) = pvar(:,Mesh%JMIN-k,Mesh%KMIN-j+1,:)
            END DO
            ! interpolate diagonal corner ghost cell values from adjacent cells in
            ! both directions
            pvar(:,Mesh%JMIN-j,Mesh%KMIN-k,:) = 0.5 * (pvar(:,Mesh%JMIN,Mesh%KMIN-j,:) &
                 + pvar(:,Mesh%JMIN-j,Mesh%KMIN,:))  !hier habe ich was korrigiert! Ist das richtig?
          END DO
        END IF

        ! bottom north
        IF(this%Boundary(4)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO k=j+1,Mesh%GKNUM
              pvar(:,Mesh%JMAX+j,Mesh%KMIN-k,:) = pvar(:,Mesh%JMAX+j-1,Mesh%KMIN-k,:)
              pvar(:,Mesh%JMAX+k,Mesh%KMIN-j,:) = pvar(:,Mesh%JMAX+k,Mesh%KMIN-j+1,:)
            END DO
            pvar(:,Mesh%JMAX+j,Mesh%KMIN-j,:) = 0.5 * (pvar(:,Mesh%JMAX,Mesh%KMIN-j,:) &
                 + pvar(:,Mesh%JMAX+j,Mesh%KMIN,:))
          END DO
        END IF

        ! top north
        IF(this%Boundary(2)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO k=j+1,Mesh%GKNUM
              pvar(:,Mesh%JMAX+j,Mesh%KMAX+k,:) = pvar(:,Mesh%JMAX+j-1,Mesh%KMAX+k,:)
              pvar(:,Mesh%JMAX+k,Mesh%KMAX+j,:) = pvar(:,Mesh%JMAX+k,Mesh%KMAX+j-1,:)
            END DO
            pvar(:,Mesh%JMAX+j,Mesh%KMAX+j,:) = 0.5 * (pvar(:,Mesh%JMAX,Mesh%KMAX+j,:) &
                 + pvar(:,Mesh%JMAX+j,Mesh%KMAX,:))
          END DO
        END IF

        ! top south
        IF(this%Boundary(3)%p%PhysicalCorner) THEN
          DO j=1,Mesh%GJNUM
            DO k=j+1,Mesh%GKNUM
              pvar(:,Mesh%JMIN-j,Mesh%KMAX+k,:) = pvar(:,Mesh%JMIN-j+1,Mesh%KMAX+k,:)
              pvar(:,Mesh%JMIN-k,Mesh%KMAX+j,:) = pvar(:,Mesh%JMIN-k,Mesh%KMAX+j-1,:)
            END DO
            pvar(:,Mesh%JMIN-j,Mesh%KMAX+j,:) = 0.5 * (pvar(:,Mesh%JMIN,Mesh%KMAX+j,:) &
                 + pvar(:,Mesh%JMIN-j,Mesh%KMAX,:))
          END DO
        END IF
      ELSE IF (Mesh%JNUM.EQ.1.AND.Mesh%INUM.NE.1.AND.Mesh%KNUM.NE.1) THEN
        !TODO: Hier hatte Jannes JNUM durch GNUM ersetzt...warum? evtl GJNUM
        !gemeint?
         ! bottom west
        IF(this%Boundary(1)%p%PhysicalCorner) THEN
          DO i=1,Mesh%GINUM
            DO k=i+1,Mesh%GKNUM
              ! copy data from adjacent ghost cells for off diagonal corner ghost cells
              pvar(Mesh%IMIN-i,:,Mesh%KMIN-k,:) = pvar(Mesh%IMIN-i+1,:,Mesh%KMIN-k,:)
              pvar(Mesh%IMIN-k,:,Mesh%KMIN-i,:) = pvar(Mesh%IMIN-k,:,Mesh%KMIN-i+1,:)
            END DO
            ! interpolate diagonal corner ghost cell values from adjacent cells in
          ! both directions
            pvar(Mesh%IMIN-i,:,Mesh%KMIN-i,:) = 0.5 * (pvar(Mesh%IMIN-i,:,Mesh%KMIN,:) &
                 + pvar(Mesh%IMIN,:,Mesh%KMIN-i,:))
          END DO
        END IF

        ! bottom east
        IF(this%Boundary(4)%p%PhysicalCorner) THEN
          DO i=1,Mesh%GINUM
            DO k=i+1,Mesh%GKNUM
              pvar(Mesh%IMAX+i,:,Mesh%KMIN-k,:) = pvar(Mesh%IMAX+i-1,:,Mesh%KMIN-k,:)
              pvar(Mesh%IMAX+k,:,Mesh%KMIN-i,:) = pvar(Mesh%IMAX+k,:,Mesh%KMIN-i+1,:)
            END DO
            pvar(Mesh%IMAX+i,:,Mesh%KMIN-i,:) = 0.5 * (pvar(Mesh%IMAX,:,Mesh%KMIN-i,:) &
                 + pvar(Mesh%IMAX+i,:,Mesh%KMIN,:))
          END DO
        END IF

        ! top east
        IF(this%Boundary(2)%p%PhysicalCorner) THEN
          DO i=1,Mesh%GJNUM
            DO k=i+1,Mesh%GKNUM
              pvar(Mesh%IMAX+i,:,Mesh%KMAX+k,:) = pvar(Mesh%IMAX+i-1,:,Mesh%KMAX+k,:)
              pvar(Mesh%IMAX+k,:,Mesh%KMAX+i,:) = pvar(Mesh%IMAX+k,:,Mesh%KMAX+i-1,:)
            END DO
            pvar(Mesh%IMAX+i,:,Mesh%KMAX+i,:) = 0.5 * (pvar(Mesh%IMAX,:,Mesh%KMAX+i,:) &
                 + pvar(Mesh%IMAX+i,:,Mesh%KMAX,:))
          END DO
        END IF

        ! top west
        IF(this%Boundary(3)%p%PhysicalCorner) THEN
          DO i=1,Mesh%GJNUM
            DO k=i+1,Mesh%GKNUM
              pvar(Mesh%IMIN-i,:,Mesh%KMAX+k,:) = pvar(Mesh%IMIN-i+1,:,Mesh%KMAX+k,:)
              pvar(Mesh%IMIN-k,:,Mesh%KMAX+i,:) = pvar(Mesh%IMIN-k,:,Mesh%KMAX+i-1,:)
            END DO
            pvar(Mesh%IMIN-i,:,Mesh%KMAX+i,:) = 0.5 * (pvar(Mesh%IMIN,:,Mesh%KMAX+i,:) &
                 + pvar(Mesh%IMIN-i,:,Mesh%KMAX,:))
          END DO
        END IF
      !ELSE
      !  CALL this%warning('CenterBoundary', 'No Setting for corner-GC found! Only for 2D implemented!')
      END IF
      ! TODO Ghost cells have to to set for 3D

    ! convert primitive variables in ghost cells
    CALL Physics%Convert2Conservative(Mesh,Mesh%IGMIN,Mesh%IMIN-1, &
         Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMAX+1,Mesh%IGMAX, &
         Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,    &
         Mesh%JGMIN,Mesh%JMIN-1,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,    &
         Mesh%JMAX+1,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,    &
         Mesh%JMIN,Mesh%JMAX,Mesh%KGMIN,Mesh%KMIN-1,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,    &
         Mesh%JMIN,Mesh%JMAX,Mesh%KMAX+1,Mesh%KGMAX,pvar,cvar)
  END SUBROUTINE CenterBoundary


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    INTEGER                               :: dir
    !------------------------------------------------------------------------!
    ! loop over all boundaries
    DO dir=1,6
       CALL this%Boundary(dir)%p%Finalize()
    END DO
  END SUBROUTINE Finalize

END MODULE boundary_generic_mod
