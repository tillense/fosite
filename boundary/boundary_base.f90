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
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief Base boundary module
!!
!----------------------------------------------------------------------------!
MODULE boundary_base_mod
  USE logging_base_mod
  USE mesh_base_mod
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
  PRIVATE
  TYPE, ABSTRACT, EXTENDS(logging_base) :: boundary_base
    !> boundary orientation: west, east, south, north
    CLASS(logging_base), ALLOCATABLE   :: direction
    INTEGER                            :: IMID,        & !< i index of cell in the middle
                                          JMID,        & !< j index of cell in the middle
                                          KMID,        & !< j index of cell in the middle
                                          nohdim         !< dimension of Noh problem
    LOGICAL                            :: first_call     !< used in far-field bc
    LOGICAL, DIMENSION(:), POINTER     ::              &
                                          reflX,       & !< mask array for reflecting bc
                                          reflY,       & !< mask array for reflecting bc
                                          reflZ          !< mask array for reflecting bc
    LOGICAL, DIMENSION(:,:,:), POINTER :: fixed          !< mask array for fixed bc
    INTEGER, DIMENSION(:,:,:), POINTER :: cbtype         !< custom boundary condition type
    REAL, DIMENSION(:,:,:), POINTER    :: invr,        & !< inverse distance to center
                                          Rscale,      & !< radial scaling constants
                                          invRscale,   & !< inverse radial scaling constants
                                          xvar,        & !< characteristic variables for absorbing bc
                                          lambda,      & !< eigenvalues for absorbing bc
                                          Rinv,        & !< Riemann invariants at the boundary
                                          RinvInf        !< far field Riemann invariants
    REAL                               :: shearval       !< shear box parameter
    LOGICAL                            :: PhysicalCorner !< Is the left corner physical?
    !> boundary data array, e.g. for fixed,
    !! custom and noslip and boundary conditions
    REAL, DIMENSION(:,:,:,:), POINTER  :: &
                         data
    !> \todo Probably obsolte variable. Remove it!
    REAL, DIMENSION(:,:,:,:), POINTER  :: &
                         accel => NULL() !< pointer to gravitational accel
#ifdef PARALLEL
    !> \name Variables in Parallel Mode
    REAL, DIMENSION(:,:,:,:), POINTER  :: sendbuf,     & !< send buffer for boundary data
                                          recvbuf        !< receive buffer for boundary data
#endif
  CONTAINS

    PROCEDURE :: InitBoundary
    PROCEDURE (SetBoundaryData), DEFERRED :: SetBoundaryData
    PROCEDURE :: FinalizeBoundary
    PROCEDURE :: GetDirection
  END TYPE boundary_base

  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  ABSTRACT INTERFACE
    PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,pvar)
      IMPORT boundary_base,mesh_base,physics_base
      IMPLICIT NONE
      CLASS(boundary_base), INTENT(INOUT)    :: this
      CLASS(mesh_base),     INTENT(IN)    :: Mesh
      CLASS(physics_base),  INTENT(IN)    :: Physics
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(INOUT) :: pvar
    END SUBROUTINE
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> \name Public Attributes
  CHARACTER(LEN=32), DIMENSION(6), PARAMETER :: &
        direction_name = (/'  west', '  east', ' south', ' north', 'bottom', '   top' /)
  !< string literal for each orientation
#ifdef PARALLEL
  INTEGER, PARAMETER :: NONE            = 0
#endif
  INTEGER, PARAMETER :: NO_GRADIENTS    = 1
  INTEGER, PARAMETER :: PERIODIC        = 2
  INTEGER, PARAMETER :: REFLECTING      = 3
  INTEGER, PARAMETER :: AXIS            = 4
  INTEGER, PARAMETER :: FOLDED          = 5
  INTEGER, PARAMETER :: FIXED           = 6
  INTEGER, PARAMETER :: EXTRAPOLATION   = 7
  INTEGER, PARAMETER :: NOH2D           = 8
  INTEGER, PARAMETER :: NOH3D           = 9
  INTEGER, PARAMETER :: NOSLIP          = 10
  INTEGER, PARAMETER :: CUSTOM          = 11
  INTEGER, PARAMETER :: FARFIELD        = 12
  INTEGER, PARAMETER :: ABSORBING       = 13
  INTEGER, PARAMETER :: DMR             = 14
  INTEGER, PARAMETER :: SHEARING        = 15
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       boundary_base, &
       ! constants
#ifdef PARALLEL
       NONE,&
#endif
       NO_GRADIENTS, PERIODIC, REFLECTING, AXIS, FOLDED, FIXED, EXTRAPOLATION, &
       NOH2D, NOH3D, NOSLIP, CUSTOM, FARFIELD, ABSORBING, DMR, SHEARING!, &
!       CUSTOM_NOGRAD, CUSTOM_PERIOD, CUSTOM_REFLECT, CUSTOM_REFLNEG, &
!       CUSTOM_EXTRAPOL, CUSTOM_FIXED, CUSTOM_LOGEXPOL, &
!       CUSTOM_OUTFLOW, CUSTOM_KEPLER, CUSTOM_ANGKEPLER, CUSTOM_POISSON
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary(this,Mesh,Physics,bctype,bcname,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_base), INTENT(INOUT)       :: this
    CLASS(mesh_base),     INTENT(IN)          :: Mesh
    CLASS(physics_base),  INTENT(IN)          :: Physics
    TYPE(Dict_TYP),       INTENT(IN), POINTER :: config
    INTEGER,              INTENT(IN)          :: bctype,dir
    CHARACTER(LEN=32),    INTENT(IN)          :: bcname
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, PARAMETER    :: strlen = 32
    CHARACTER(LEN=strlen) :: sendbuf
    CHARACTER(LEN=strlen) :: recvbuf
    INTEGER               :: status(MPI_STATUS_SIZE)
    INTEGER               :: ierror
#endif
    !------------------------------------------------------------------------!
    ! set boundary condition
    CALL this%InitLogging(bctype,bcname)
    ALLOCATE(this%direction)
    CALL this%direction%InitLogging(dir,direction_name(dir))
    ! check for wrong direction
    SELECT CASE(dir)
    CASE(WEST,EAST,NORTH,SOUTH,TOP,BOTTOM)
       ! ok
    CASE DEFAULT
       CALL this%Error("InitBoundary_common", "Unknown direction")
    END SELECT
    !! set direction
    !todo: Was passiert hiermit???
    !CALL InitCommon(this%direction,direction,direction_name(direction))

    IF(((Physics%GetType().EQ.EULER2D_ISOIAMT).OR.&
        (Physics%GetType().EQ.EULER2D_IAMT)).AND. &
       ((dir.EQ.NORTH).OR.(dir.EQ.SOUTH)).AND.    &
       (.NOT.((bctype.EQ.PERIODIC)                &
#ifdef PARALLEL
              .OR.(bctype.EQ.NONE)                &
#endif
              )))                                 &
      CALL this%Error("InitBoundary_one", "All IAMT Physics need periodic" &
        // " boundary conditions in NORTH/SOUTH direction")

    ! print some information
#ifdef PARALLEL
    ! send boundary information to the rank 0 process;
    ! we only need this to synchronize the output
    IF (this%GetRank() .EQ. 0 .AND. this%GetRank().EQ.Mesh%rank0_boundaries(dir)) THEN
       ! print output without communication
#endif
       CALL this%Info(" BOUNDARY-> condition:         " //  TRIM(this%direction%GetName()) &
            // " " // TRIM(this%GetName()), this%GetRank())
#ifdef PARALLEL
    ELSE IF (this%GetRank().EQ.Mesh%rank0_boundaries(dir)) THEN
       ! send info to root
       sendbuf = TRIM(this%direction%GetName())//" "//TRIM(this%GetName())
       CALL MPI_SEND(sendbuf,strlen,MPI_CHARACTER,0, &
                     0,MPI_COMM_WORLD,ierror)
    ELSE IF (this%GetRank().EQ.0) THEN
       ! receive input from rank0_boundaries(dir)
       CALL MPI_RECV(recvbuf,strlen,MPI_CHARACTER,Mesh%rank0_boundaries(dir),&
                     MPI_ANY_TAG,MPI_COMM_WORLD,status,ierror)
       CALL this%Info(" BOUNDARY-> condition:         " // TRIM(recvbuf),this%GetRank())
    END IF
#endif
  END SUBROUTINE InitBoundary

  !> \public Get the direction number.
  !! \return direction number
  PURE FUNCTION GetDirection(this) RESULT(dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_base), INTENT(IN) :: this !< \param [in] this boundary type
    INTEGER                          :: dir
    !------------------------------------------------------------------------!
    dir = this%direction%GetType()
  END FUNCTION GetDirection


  SUBROUTINE FinalizeBoundary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_base),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("CloseBoundary_one","not initialized")
#ifdef PARALLEL
    ! deallocate MPI send/recv buffers
    DEALLOCATE(this%sendbuf,this%recvbuf)
#endif
    DEALLOCATE(this%direction)
  END SUBROUTINE FinalizeBoundary

END MODULE boundary_base_mod
