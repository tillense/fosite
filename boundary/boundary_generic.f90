!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2006-2016                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief Generic boundary module
!!
!! This module provides the generic interface routines to all boundary
!! modules.
!!
!----------------------------------------------------------------------------!
MODULE boundary_generic_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE boundary_base_mod
!  USE boundary_reflecting_mod
  USE boundary_nogradients_mod
  USE boundary_periodic_mod
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
  TYPE, PRIVATE :: boundary_p
    CLASS(boundary_base), ALLOCATABLE :: p
  END TYPE

  TYPE, EXTENDS(logging_base) :: boundary_generic
    PRIVATE
    !> \name Variables
    !CLASS(boundary_base) :: Boundary(4)
    TYPE(boundary_p) :: Boundary(6)

    LOGICAL           :: PhysicalCorner  !< Is the left corner physical?
#ifdef PARALLEL
    !> \name Variables in Parallel Mode
    REAL,DIMENSION(:,:,:,:),POINTER :: &
                         sendbuf, &      !< send buffer for boundary data
                         recvbuf         !< receive buffer for boundary data
#endif
  CONTAINS
    PROCEDURE :: InitBoundary
    PROCEDURE :: CenterBoundary
    FINAL :: Finalize
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
    CLASS(Boundary_generic),ALLOCATABLE :: Boundary
    CLASS(mesh_base),INTENT(IN)         :: Mesh
    CLASS(physics_base),INTENT(IN)      :: Physics
    TYPE(Dict_TYP),POINTER              :: config,IO
    !------------------------------------------------------------------------!
    ALLOCATE(Boundary)
    CALL Boundary%InitBoundary(Mesh,Physics,config,IO)
  END SUBROUTINE

  SUBROUTINE InitBoundary(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)           :: Mesh
    CLASS(physics_base),INTENT(IN)        :: Physics
    TYPE(Dict_TYP),POINTER                :: config,IO
    INTEGER         :: western, eastern, southern, northern, bottomer, topper
    !------------------------------------------------------------------------!
    INTEGER               :: new(6)
    LOGICAL, DIMENSION(3) :: periods = .FALSE.
    INTEGER               :: dir
#ifdef PARALLEL
    INTEGER               :: comm_old
    INTEGER               :: sizeofreal, ignum, jgnum, kgnum, twoslices
    INTEGER               :: ierr
    LOGICAL, DIMENSION(SIZE(Mesh%dims)) :: remain_dims = .FALSE.
#endif
    !------------------------------------------------------------------------!
    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("InitBoundary","physics and/or mesh module uninitialized")

    CALL GetAttr(config, "western", western)
    CALL GetAttr(config, "eastern", eastern)
    CALL GetAttr(config, "southern", southern)
    CALL GetAttr(config, "northern", northern)
    CALL GetAttr(config, "bottomer", bottomer)
    CALL GetAttr(config, "topper", topper)

    new(WEST)   = western
    new(EAST)   = eastern
    new(SOUTH)  = southern
    new(NORTH)  = northern
    new(BOTTOM) = bottomer
    new(TOP)    = topper

#ifdef PARALLEL
    ! define connections
    IF (Mesh%mycoords(1).NE.0)  new(WEST) = NONE
    IF (Mesh%mycoords(1).NE.Mesh%dims(1)-1)  new(EAST) = NONE
    IF (Mesh%mycoords(2).NE.0)  new(SOUTH) = NONE
    IF (Mesh%mycoords(2).NE.Mesh%dims(2)-1)  new(NORTH) = NONE
    IF (Mesh%mycoords(3).NE.0) new(BOTTOM) = NONE
    IF (Mesh%mycoords(3).NE.Mesh%dims(3)-1) new(TOP) = NONE
#endif

    ! initialize every boundary
    ! IMPORTANT: do this before anything else
    DO dir=WEST,TOP
      SELECT CASE(new(dir))
!      CASE(REFLECTING)
!        ALLOCATE(boundary_reflecting::this%Boundary(dir)%p)
      CASE(NO_GRADIENTS)
        ALLOCATE(boundary_nogradients::this%Boundary(dir)%p)
      CASE(PERIODIC)
        ALLOCATE(boundary_periodic::this%Boundary(dir)%p)
      CASE DEFAULT
        CALL this%Error("new_boundary","Unkown boundary type.")
      END SELECT
      SELECT TYPE(obj => this%Boundary(dir)%p)
!      TYPE IS (boundary_reflecting)
!        CALL obj%InitBoundary_reflecting(Mesh,Physics,dir,config)
      TYPE IS (boundary_nogradients)
        CALL obj%InitBoundary_nogradients(Mesh,Physics,dir,config)
      TYPE IS (boundary_periodic)
        CALL obj%InitBoundary_periodic(Mesh,Physics,dir,config)
      END SELECT
    END DO

    ! check periodicity
    IF (western.EQ.PERIODIC.AND.eastern.EQ.PERIODIC) THEN
       periods(1) = .TRUE.
    ELSE IF (western.EQ.PERIODIC.NEQV.eastern.EQ.PERIODIC) THEN
       CALL this%Error("InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF
    IF (southern.EQ.PERIODIC.AND.northern.EQ.PERIODIC) THEN
       periods(2) = .TRUE.
    ELSE IF (southern.EQ.PERIODIC.NEQV.northern.EQ.PERIODIC) THEN
       CALL this%Error("InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF
    IF (top.EQ.PERIODIC.AND.bottom.EQ.PERIODIC) THEN
       periods(3) = .TRUE.
    ELSE if (top.EQ.PERIODIC.NEQV.bottom.EQ.PERIODIC) THEN
       CALL this%Error("InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF

    ! This sets this(dir)%PhysicalCorner to .True., if
    ! there are no periodic boundary conditions and if
    ! is a corner without inner boundaries.
    ! At physical corners the corner value have to be
    ! interpolated in CenterBoundary.
    ! e.g. this(1)%PhysicalCorner describes the Southwest corner.
    DO dir=1,6
      this%Boundary(dir)%p%PhysicalCorner = .FALSE.
    END DO
    ! TODO Help needed for parallelisation
    IF(.NOT.ANY(periods)) THEN
#ifdef PARALLEL
      IF(Mesh%mycoords(1).EQ.0) THEN
        IF(Mesh%mycoords(2).EQ.0) &
          this(1)%PhysicalCorner = .TRUE.
        IF(Mesh%mycoords(2).EQ.Mesh%dims(2)-1) &
          this(4)%PhysicalCorner = .TRUE.
      END IF
      IF(Mesh%mycoords(1).EQ.Mesh%dims(1)-1) THEN
        IF(Mesh%mycoords(2).EQ.0) &
          this(3)%PhysicalCorner = .TRUE.
        IF(Mesh%mycoords(2).EQ.Mesh%dims(2)-1) &
          this(2)%PhysicalCorner = .TRUE.
      END IF
#else
      DO dir=1,6
        this%Boundary(dir)%p%PhysicalCorner = .TRUE.
      END DO
#endif
    END IF


#ifdef PARALLEL
    ! create new cartesian communicator using Mesh%comm_cart
    ! and account for the periodicity
    ! IMPORTANT: disable reordering of nodes
    comm_old = Mesh%comm_cart
    CALL MPI_Cart_create(comm_old,SIZE(Mesh%dims),Mesh%dims,periods,.FALSE., &
         Mesh%comm_cart,ierr)

    ! save ranks of neighbor processes
    CALL MPI_Cart_shift(Mesh%comm_cart,0,1,Mesh%neighbor(WEST),Mesh%neighbor(EAST),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,1,1,Mesh%neighbor(SOUTH),Mesh%neighbor(NORTH),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,2,1,Mesh%neighbor(TOP), Mesh%neighbor(BOTTOM),ierr)

    ! create communicators for every column and row of the cartesian
    !	topology (used eg. for fargo shifts)
    remain_dims = (/ .FALSE., .TRUE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Icomm,ierr)
    remain_dims = (/ .TRUE., .FALSE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Jcomm,ierr)
    remain_dims = (/ .TRUE., .FALSE. /)
    CALL MPI_Cart_Sub(Mesh%comm_cart,remain_dims,Mesh%Kcomm,ierr)

    ! allocate memory for boundary data buffers
    ALLOCATE(this(WEST)%sendbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(WEST)%recvbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(EAST)%sendbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(EAST)%recvbuf(Mesh%GINUM,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(SOUTH)%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(SOUTH)%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(NORTH)%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(NORTH)%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%GJNUM,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
         this(BOTTOM)%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
         this(BOTTOM)%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
         this(TOP)%sendbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
         this(TOP)%recvbuf(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%GKNUM,Physics%VNUM), &
         STAT=ierr)
    IF (ierr.NE.0) THEN
       CALL Error(this(WEST),"InitBoundary", &
            "Unable to allocate memory for data buffers.")
    END IF
    DO dir=WEST,BOTTOM
      this(dir)%recvbuf = 0.
      this(dir)%sendbuf = 0.
    END DO
#endif

  END SUBROUTINE InitBoundary


  SUBROUTINE CenterBoundary(this,Mesh,Physics,time,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_generic),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)           :: Mesh
    CLASS(physics_base),INTENT(IN)        :: Physics
    REAL                                  :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM) &
                                          :: pvar, cvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
#ifdef PARALLEL
    INTEGER       :: req(4)
    INTEGER       :: ierr
#ifdef MPI_USE_SENDRECV
    INTEGER       :: status(MPI_STATUS_SIZE)
#else
    INTEGER       :: status(MPI_STATUS_SIZE,4)
#endif
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: time
    INTENT(INOUT) :: cvar
    !------------------------------------------------------------------------!
    CALL Physics%Convert2Primitive(Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN, &
             Mesh%JMAX,Mesh%KMIN,Mesh%KMAX,cvar,pvar)

    ! set physical boundary conditions at western and eastern boundaries
    IF (Mesh%INUM.GT.1) THEN
      CALL this%Boundary(WEST)%p%SetBoundaryData(Mesh,Physics,pvar)
      CALL this%Boundary(EAST)%p%SetBoundaryData(Mesh,Physics,pvar)
    END IF

    ! set physical boundary conditions at southern and northern boundaries
    IF (Mesh%JNUM.GT.1) THEN
      CALL this%Boundary(SOUTH)%p%SetBoundaryData(Mesh,Physics,pvar)
      CALL this%Boundary(NORTH)%p%SetBoundaryData(Mesh,Physics,pvar)
    END IF

    ! set physical boundary conditions at top and bottom boundaries
    IF (Mesh%KNUM.GT.1) THEN
      CALL this%Boundary(BOTTOM)%p%SetBoundaryData(Mesh,Physics,pvar)
      CALL this%Boundary(TOP)%p%SetBoundaryData(Mesh,Physics,pvar)
    END IF


! TODO Tried to extend the communication for 3D but no guarantee for correctness
!#ifdef PARALLEL
!    ! NOTE: if you want to use MPI_Sendrecv instead of nonblocking
!    ! MPI_Irecv and  MPI_Issend for exchange of ghost cell data,
!    ! you must add -DMPI_USE_SENDRECV to the compile command
!
!    ! initiate western/eastern MPI communication
!#ifdef MPI_USE_SENDRECV
!    ! send boundary data to western and receive from eastern neighbor
!    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
!        this(WEST)%sendbuf(:,:,:,:) = pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GINUM-1, &
!                                        Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
!    CALL MPI_Sendrecv(this(WEST)%sendbuf,Mesh%GJNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,this(EAST)%recvbuf, &
!         Mesh%GJNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(EAST), &
!         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
!    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM) = this(EAST)%recvbuf(:,:,:,:)
!#else
    ! receive boundary data from eastern neighbor
!    CALL MPI_Irecv(this(EAST)%recvbuf,Mesh%GJNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+WEST,Mesh%comm_cart,req(1),ierr)
    ! fill send buffer if western neighbor exists
!    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
!         this(WEST)%sendbuf(:,:,:,:) = pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GNUM-1, &
!                                          Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1:Physics%VNUM)
    ! send boundary data to western neighbor
!    CALL MPI_Issend(this(WEST)%sendbuf,Mesh%GJNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,Mesh%comm_cart,req(2),ierr)
!#endif
!#endif
!#ifdef PARALLEL
!#ifdef MPI_USE_SENDRECV
    ! send boundary data to eastern and receive from western neighbor
!   IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
!         this(EAST)%sendbuf(:,:,:) = pvar(Mesh%IMAX-Mesh%GNUM+1:Mesh%IMAX, &
!                                          Mesh%JGMIN:Mesh%JGMAX,1:Physics%VNUM)
!    CALL MPI_Sendrecv(this%Boundary(EAST)%p%sendbuf,Mesh%GNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,this(WEST)%recvbuf, &
!         Mesh%GNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(WEST), &
!         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
!    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JGMIN:Mesh%JGMAX,1:Physics%VNUM) = this%Boundary%p%(WEST)%recvbuf(:,:,:)
!#else
    ! receive boundary data from western neighbor
!    CALL MPI_Irecv(this%Boundary(WEST)%p%recvbuf,Mesh%GNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+EAST,Mesh%comm_cart,req(3),ierr)
    ! fill send buffer if eastern neighbor exists
!    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
!         this%Boundary(EAST)%p%%sendbuf(:,:,:) = pvar(Mesh%IMAX-Mesh%GNUM+1:Mesh%IMAX, &
!                                          Mesh%JGMIN:Mesh%JGMAX,1:Physics%VNUM)
    ! send boundary data to eastern neighbor
!    CALL MPI_Issend(this%Boundary(EAST)%p%sendbuf,Mesh%GNUM*(Mesh%JGMAX-Mesh%JGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,Mesh%comm_cart,req(4),ierr)
!#endif
!#endif
!#ifdef PARALLEL
!#ifndef MPI_USE_SENDRECV
    ! wait for unfinished MPI communication
!    CALL MPI_Waitall(4,req,status,ierr)
    ! copy data from recieve buffers into ghosts cells
!    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JGMIN:Mesh%JGMAX,1:Physics%VNUM) = this(WEST)%recvbuf(:,:,:)
!    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,1:Physics%VNUM) = this(EAST)%recvbuf(:,:,:)
!#endif
!#endif

!#ifdef PARALLEL
    ! initiate southern/northern MPI communication
!#ifdef MPI_USE_SENDRECV
    ! send boundary data to southern and receive from northern neighbor
!    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
!         this(SOUTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
!                                           Mesh%JMIN:Mesh%JMIN+Mesh%GNUM-1,1:Physics%VNUM)
!    CALL MPI_Sendrecv(this(SOUTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,this(NORTH)%recvbuf, &
!         Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(NORTH), &
!         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
!    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,1:Physics%VNUM) = this(NORTH)%recvbuf(:,:,:)
!#else
    ! receive boundary data from northern neighbor
!    CALL MPI_Irecv(this(NORTH)%recvbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+SOUTH,Mesh%comm_cart,req(1),ierr)
!    ! fill send buffer if southern neighbor exists
!    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
!         this(SOUTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
!                                           Mesh%JMIN:Mesh%JMIN+Mesh%GNUM-1,1:Physics%VNUM)
!    ! send boundary data to southern neighbor
!    CALL MPI_Issend(this(SOUTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,Mesh%comm_cart,req(2),ierr)
!#endif
!#endif
!#ifdef PARALLEL
!#ifdef MPI_USE_SENDRECV
    ! send boundary data to northern and receive from southern neighbor
!    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
!         this(NORTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
!                                           Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX,1:Physics%VNUM)

!    CALL MPI_Sendrecv(this(NORTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,this(SOUTH)%recvbuf, &
!         Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%VNUM,DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH), &
!         MPI_ANY_TAG,Mesh%comm_cart,status,ierr)
!    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,1:Physics%VNUM) = this(SOUTH)%recvbuf(:,:,:)
!#else
!    ! receive boundary data from southern neighbor
!    CALL MPI_Irecv(this(SOUTH)%recvbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+NORTH,Mesh%comm_cart,req(3),ierr)
!    ! fill send buffer if northern neighbor exists
!    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
!         this(NORTH)%sendbuf(:,:,:) = pvar(Mesh%IGMIN:Mesh%IGMAX, &
!                                           Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX,1:Physics%VNUM)
!    ! send boundary data to northern neighbor
!    CALL MPI_Issend(this(NORTH)%sendbuf,Mesh%GNUM*(Mesh%IGMAX-Mesh%IGMIN+1)*Physics%vnum, &
!         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,Mesh%comm_cart,req(4),ierr)
!#endif
!#endif
!#ifdef PARALLEL
!#ifndef MPI_USE_SENDRECV
!    ! wait for unfinished MPI communication
!    CALL MPI_Waitall(4,req,status,ierr)
!    IF (Mesh%neighbor(SOUTH).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1,1:Physics%VNUM) = this(SOUTH)%recvbuf(:,:,:)
!    IF (Mesh%neighbor(NORTH).NE.MPI_PROC_NULL) &
!         pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX,1:Physics%VNUM) = this(NORTH)%recvbuf(:,:,:)
!#endif
!#endif
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
      ELSE IF (Mesh%INUM.EQ.1.AND.Mesh%GNUM.NE.1.AND.Mesh%KNUM.NE.1) THEN
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
            pvar(:,Mesh%JMIN-j,Mesh%KMIN-k,:) = 0.5 * (pvar(Mesh%IMIN-j,Mesh%JMIN,:,:) &
                 + pvar(:,Mesh%JMIN-j,Mesh%KMIN,:))
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
              pvar(:,Mesh%JMAX+j,Mesh%KMAX+k,:) = pvar(:,Mesh%JMAX+j-1,Mesh%KMAX+i,:)
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
      ELSE IF (Mesh%GNUM.EQ.1.AND.Mesh%INUM.NE.1.AND.Mesh%KNUM.NE.1) THEN
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
              pvar(Mesh%IMAX+i,:,Mesh%KMIN-i,:) = pvar(Mesh%IMAX+k,:,Mesh%KMIN-i+1,:)
            END DO
            pvar(Mesh%IMAX+i,:,Mesh%KMIN-i,:) = 0.5 * (pvar(Mesh%IMAX,:,Mesh%KMIN-i,:) &
                 + pvar(Mesh%IMAX+i,:,Mesh%KMIN,:))
          END DO
        END IF

        ! top east
        IF(this%Boundary(2)%p%PhysicalCorner) THEN
          DO i=1,Mesh%GJNUM
            DO k=i+1,Mesh%GKNUM
              pvar(Mesh%IMAX+i,:,Mesh%KMAX+k,:) = pvar(Mesh%IMAX+i-1,:,Mesh%KMAX+i,:)
              pvar(Mesh%IMAX+i,:,Mesh%KMAX+i,:) = pvar(Mesh%IMAX+i,:,Mesh%KMAX+i-1,:)
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
      END IF
      ! TODO Ghost cells have to to set for 3D

    ! convert primitive variables in ghost cells
    CALL Physics%Convert2Conservative(Mesh,Mesh%IGMIN,Mesh%IMIN-1,&
         Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMAX+1,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,&
         Mesh%JGMIN,Mesh%JMIN-1,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,&
         Mesh%JMAX+1,Mesh%JGMAX,Mesh%KGMIN,Mesh%KGMAX,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,&
         Mesh%JMIN,Mesh%JMAX,Mesh%KGMIN,Mesh%KMIN-1,pvar,cvar)
    CALL Physics%Convert2Conservative(Mesh,Mesh%IMIN,Mesh%IMAX,&
         Mesh%JMIN,Mesh%JMAX,Mesh%KMAX+1,Mesh%KGMAX,pvar,cvar)


!  CONTAINS
!
!    ! set boundary data for (real) physical boundaries in direction "dir"
!    SUBROUTINE SetBoundaryData(dir)
!      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: dir
!!CDIR IEXPAND
!      SELECT CASE(GetType(this(dir)))
!      CASE(NO_GRADIENTS)
!         CALL CenterBoundary_nogradients(this(dir),Mesh,Physics,pvar)
!      CASE(PERIODIC)
!         ! do nothing in parallel version, because periodicity is
!         ! handled via MPI communication
!#ifndef PARALLEL
!         CALL CenterBoundary_periodic(this(dir),Mesh,Physics,pvar)
!#endif
!      CASE(REFLECTING)
!         CALL CenterBoundary_reflecting(this(dir),Mesh,Physics,pvar)
!      CASE(AXIS)
!         CALL CenterBoundary_axis(this(dir),Mesh,Physics,pvar)
!      CASE(FOLDED)
!         CALL CenterBoundary_folded(this(dir),Mesh,Physics,pvar)
!      CASE(FIXED)
!         CALL CenterBoundary_fixed(this(dir),Mesh,Physics,pvar)
!      CASE(EXTRAPOLATION)
!         CALL CenterBoundary_extrapolation(this(dir),Mesh,Physics,pvar)
!      CASE(NOH2D,NOH3D)
!         CALL CenterBoundary_noh(this(dir),Mesh,Physics,time,pvar)
!      CASE(NOSLIP)
!         CALL CenterBoundary_noslip(this(dir),Mesh,Physics,pvar)
!      CASE(CUSTOM)
!         CALL CenterBoundary_custom(this(dir),Mesh,Physics,pvar)
!      CASE(FARFIELD)
!         CALL CenterBoundary_farfield(this(dir),Mesh,Physics,pvar)
!      CASE(ABSORBING)
!         CALL CenterBoundary_absorbing(this(dir),Mesh,Physics,cvar,pvar)
!      CASE(DMR)
!         CALL CenterBoundary_dmr(this(dir),Mesh,Physics,time,pvar)
!      CASE(SHEARING)
!         CALL CenterBoundary_shearing(this(dir),Mesh,Physics,time,pvar)
!      END SELECT
!    END SUBROUTINE SetBoundaryData
!
  END SUBROUTINE CenterBoundary


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(boundary_generic), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    INTEGER                               :: dir
    !------------------------------------------------------------------------!
    ! loop over all boundaries
    DO dir=1,6
       DEALLOCATE(this%Boundary(dir)%p)
    END DO
  END SUBROUTINE Finalize

END MODULE boundary_generic_mod
