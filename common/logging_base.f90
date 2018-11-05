!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: common_types.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief Basic fosite module
!!
!! This module defines the basic data types and methods inherited by all
!! fosite modules.
!!
!----------------------------------------------------------------------------!
MODULE logging_base_mod
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
  !> common data structure
  !!
  !! This type is imported by all basic data types in fosite. Its main purpose
  !! is the specification of distinct objects which could by identified by
  !! its module type number. In parallel execution it stores the rank of the
  !! current MPI process and the total number of processes.
  TYPE logging_base
    PRIVATE
    INTEGER           :: mod_type                !< module type number
    CHARACTER(LEN=32) :: mod_name                !< module type name
    INTEGER           :: err                     !< error code
    LOGICAL           :: isinitialized = .FALSE. !< init status
    INTEGER, POINTER  :: myrank                  !< rank of parallel process
    INTEGER, POINTER  :: ppnum                   !< number of parallel processes
    LOGICAL, POINTER  :: parinit                 !< init status of parallel process
  CONTAINS
    PROCEDURE :: InitLogging
    FINAL     :: Finalize
    PROCEDURE :: GetType
    PROCEDURE :: GetName
    PROCEDURE :: GetRank
    PROCEDURE :: GetNumProcs
    PROCEDURE :: Initialized
    PROCEDURE :: Info
    PROCEDURE :: Warning
    PROCEDURE :: Error
  END TYPE logging_base
  ! these variables should be the same for all objects
  ! of the current process
#ifdef PARALLEL
  INTEGER, SAVE          :: DEFAULT_MPI_REAL = MPI_REAL              !< default real type for MPI
  INTEGER, SAVE          :: DEFAULT_MPI_2REAL = MPI_2REAL            !< default 2real type for MPI
  INTEGER, SAVE          :: DEFAULT_MPI_COMPLEX = MPI_DOUBLE_COMPLEX !< default real type for MPI
  REAL, PARAMETER        :: dummy   = 1.0                            !< check default real type
#endif
  INTEGER, PARAMETER     :: STDERR  = 0                              !< fortran stderr unit
  INTEGER, PARAMETER     :: STDOUT  = 6                              !< fortran stdout unit
  INTEGER, SAVE, TARGET  :: myrank  = 0                              !< MPI rank
  INTEGER, SAVE, TARGET  :: ppnum   = 1                              !< MPI number of processes
  LOGICAL, SAVE, TARGET  :: parinit = .FALSE.                        !< MPI initialization status
  CHARACTER(LEN=1), SAVE :: prefix  = ' '                            !< preceds info output
  !--------------------------------------------------------------------------!
  PUBLIC ::                 &
       ! types
       logging_base,        &
#ifdef PARALLEL
       DEFAULT_MPI_REAL,    &
       DEFAULT_MPI_2REAL,   &
       DEFAULT_MPI_COMPLEX, &
#endif
       SetPrefix
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of common class.
  !!
  !! Sets the module type (number) and name. In addition it initializes
  !! other internal variables including those related to parallel execution
  !! mode such as the MPI rank and the total number of MPI processes.
  !! The default MPI real data types are determined.
  SUBROUTINE InitLogging(this,t,n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(INOUT) :: this !< \param [out] this module type and name
    INTEGER                            :: t    !< \param [in]  t module type identification number
    CHARACTER(LEN=*)                   :: n    !< \param [in]  n module type name
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: t,n
    !------------------------------------------------------------------------!
    this%mod_type = t
    this%mod_name = n
    this%myrank => myrank
    this%ppnum => ppnum
    this%parinit => parinit
#ifdef PARALLEL
    IF (.NOT.parinit) THEN
       CALL MPI_Comm_rank(mpi_comm_world,this%myrank,this%err)
       CALL MPI_Comm_size(mpi_comm_world,this%ppnum,this%err)
       this%parinit = .TRUE.
       ! determine the default MPI data type for real numbers
       SELECT CASE (SELECTED_REAL_KIND(PRECISION(dummy)))
       CASE(4)
          DEFAULT_MPI_REAL = MPI_REAL4
          DEFAULT_MPI_2REAL = MPI_2REAL
       CASE(8)
          DEFAULT_MPI_REAL = MPI_REAL8
          DEFAULT_MPI_2REAL = MPI_2DOUBLE_PRECISION
       CASE DEFAULT
          CALL Warning(this,"InitCommon","Cannot determine default MPI real types.")
       END SELECT
    END IF
#endif
    this%err = 0
    this%isinitialized  = .TRUE.
  END SUBROUTINE InitLogging


  !> \public Destructor of common class.
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(logging_base), INTENT(INOUT) :: this !< \param [in,out] this module type and name
    !------------------------------------------------------------------------!
    this%isinitialized  = .FALSE.
  END SUBROUTINE Finalize


  !> \public Get the module type number.
  !! \return module type number
  PURE FUNCTION GetType(this) RESULT(t)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN) :: this !< \param [in] this module type and name
    INTEGER                         :: t
    !------------------------------------------------------------------------!
    t = this%mod_type
  END FUNCTION GetType


  !> \public Set character preceding the info output.
  !!
  SUBROUTINE SetPrefix(val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=1), INTENT(IN) :: val !< \param [in] val preceding character
    !------------------------------------------------------------------------!
    prefix = val
  END SUBROUTINE SetPrefix


  !> \public Get the module name.
  !! \return module name
  PURE FUNCTION GetName(this) RESULT(n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN) :: this !< \param [in] this module type and name
    CHARACTER(LEN=32)               :: n
    !------------------------------------------------------------------------!
    n = this%mod_name
  END FUNCTION GetName


  !> \public Get the MPI rank.
  !! \return MPI rank
  !!
  !! This function returns the MPI rank, i.e. process number, of the current
  !! process in parallel mode. If fosite is compiled for serial execution
  !! the return value is always 0.
  PURE FUNCTION GetRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN) :: this !< \param [in] this module type and name
    INTEGER                         :: r
    !------------------------------------------------------------------------!
    r = this%myrank
  END FUNCTION GetRank


  !> \public Get the total number of MPI processes.
  !! \return number of MPI processes
  !!
  !! If the code is compiled for serial execution it returns always 1.
  PURE FUNCTION GetNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN) :: this !< \param [in] this module type and name
    INTEGER                         :: p
    !------------------------------------------------------------------------!
    p = this%ppnum
  END FUNCTION GetNumProcs


  !> \public Query initialization status
  !! \return initialization status
  !!
  !! It returns .TRUE. for initialized modules otherwise .FALSE.
  PURE FUNCTION Initialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN) :: this !< \param [in] this module type and name
    LOGICAL                         :: i
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    i = this%isinitialized .AND. this%parinit
#else
    i = this%isinitialized
#endif
  END FUNCTION Initialized


  !> \public Print information on standard output.
  !!
  !! Print the string given in `msg` preceded by the
  !! character previously set by \link SetPrefix \endlink . If the optional
  !! parameter `tostderr` is .TRUE. the information is printed to STDERR instead
  !! of STDOUT.
  !!
  !! In parallel mode it prints information only on the MPI process given by `rank`
  !! (defaults to 0 if omitted). In addition it prints the MPI rank if `node_info`
  !! is .TRUE.
  SUBROUTINE Info(this,msg,rank,node_info,tostderr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN)  :: this      !< \param [in] this module type and name
    CHARACTER(LEN=*), INTENT(IN)     :: msg       !< \param [in] msg info message
    INTEGER, OPTIONAL, INTENT(IN)    :: rank      !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL, INTENT(IN)    :: node_info !< \param [in] node_info enable rank output
    LOGICAL, OPTIONAL, INTENT(IN)    :: tostderr  !< \param [in] tostderr enable STDERR output
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER :: ierr
#endif
    INTEGER :: print_rank,output_unit
    LOGICAL :: print_node_info
    !------------------------------------------------------------------------!
    IF (PRESENT(rank)) THEN
       print_rank = rank
    ELSE
       print_rank = 0
    END IF
    IF (PRESENT(node_info)) THEN
       print_node_info = node_info
    ELSE
       print_node_info = .FALSE.
    END IF
    ! output unit for printing, defaults to STDOUT
    output_unit = STDOUT
    IF (PRESENT(tostderr).AND.tostderr) THEN
       output_unit = STDERR ! print on STDERR if requested
    ENDIF
#ifdef PARALLEL
    IF (.NOT.parinit) CALL MPI_Comm_rank(mpi_comm_world,myrank,ierr)
#endif
    ! use "myrank" here instead of "this%myrank"
    ! because "this" might be uninitialized
    IF (myrank.EQ.print_rank) THEN
      WRITE (output_unit,'(A)',ADVANCE='NO') prefix
#ifdef PARALLEL
      IF (print_node_info) &
        WRITE (output_unit,'(A,I4.4,A)',ADVANCE='NO') "NODE [", myrank, "] "
#endif
      WRITE (output_unit,'(A)') TRIM(msg)
      ! Flush the output buffer, to make sure it is written even if the
      ! program is invoked from a batch system.
      ! This is Fortran 2003 standard! In Fortran 90/95 flush is implemented
      ! as a subroutine, hence one needs a CALL statement, i.e.,
      ! CALL FLUSH(output_unit)
      FLUSH(output_unit)
    END IF
  END SUBROUTINE Info


  !> \public Print warning message on standard error.
  !!
  !! Print the warning message given by `msg` preceded by the name of the
  !! module procedure which triggers the warning.
  SUBROUTINE Warning(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN) :: this    !< \param [in] this module type and name
    CHARACTER(LEN=*), INTENT(IN)    :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*), INTENT(IN)    :: msg     !< \param [in] msg warning message
    INTEGER, OPTIONAL, INTENT(IN)   :: rank    !< \param [in] rank MPI rank
    !------------------------------------------------------------------------!
     CALL Info(this,"WARNING in " // TRIM(modproc) // ": " // TRIM(msg),&
               rank,node_info=.TRUE.,tostderr=.TRUE.)
  END SUBROUTINE Warning


  !> \public Print error message on standard error and terminate the program.
  !!
  !! Print the error message given by `msg` preceded by the name of the
  !! module procedure which triggers the error. After that abort
  !! program execution.
  SUBROUTINE Error(this,modproc,msg,rank,node_info)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(logging_base), INTENT(IN)  :: this      !< \param [in,out] this module type and name
    CHARACTER(LEN=*), INTENT(IN)     :: modproc   !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*), INTENT(IN)     :: msg       !< \param [in] msg warning message
    INTEGER, OPTIONAL, INTENT(IN)    :: rank      !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL, INTENT(IN)    :: node_info !< \param [in] node_info enable rank output
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER :: ierr
#endif
    !------------------------------------------------------------------------!
    CALL Info(this,"ERROR in " // TRIM(modproc) // ": " // TRIM(msg),&
         rank=rank,node_info=node_info,tostderr=.TRUE.)
    ! abort execution
#ifdef PARALLEL
    CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    STOP 99
#endif
  END SUBROUTINE Error

END MODULE logging_base_mod
