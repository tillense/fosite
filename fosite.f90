!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# program file: fosite.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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

!----------------------------------------------------------------------------!
!> Module fosite
!----------------------------------------------------------------------------!
MODULE fosite_mod
  USE logging_base_mod
  USE sources_generic_mod
  USE gravity_generic_mod
  USE physics_generic_mod
  USE fluxes_generic_mod
  USE reconstruction_generic_mod
  USE mesh_generic_mod
  USE boundary_generic_mod
  USE timedisc_generic_mod
  USE fileio_generic_mod
!  USE gravity_generic_mod
  USE integration
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
  ! fosite struc class (only part of fosite class)
  TYPE, PRIVATE :: FositeStruc_TYP
     CHARACTER(LEN=16)      :: name
     INTEGER                :: pos
     INTEGER                :: dim
     INTEGER                :: rank
  END TYPE FositeStruc_TYP
  INTEGER, PRIVATE, PARAMETER           :: MAXLEN = 500
  INTEGER, PRIVATE, PARAMETER           :: simtype = 1
  CHARACTER(LEN=32), PRIVATE, PARAMETER :: simname = "fosite"

  ! fosite main class
  TYPE, EXTENDS(logging_base) :: fosite
  !> Fosite data structure
  !!
  !! This is the basic data structure, which stores all other information
  !! inside of it.
  !> \name Variables
  TYPE(Dict_TYP),POINTER :: config => null()    ! global config           !
  TYPE(Dict_TYP),POINTER :: IO => null()        ! global In-/Output Dict  !
  CLASS(mesh_base),ALLOCATABLE       :: Mesh
  CLASS(fluxes_base),ALLOCATABLE     :: Fluxes
  CLASS(physics_base),ALLOCATABLE    :: Physics
  CLASS(fileio_base),ALLOCATABLE     :: Datafile
  CLASS(timedisc_base),ALLOCATABLE   :: Timedisc
  CLASS(fileio_base),ALLOCATABLE     :: Logfile
  CLASS(sources_base), POINTER &
                         :: Sources => null()   !< list of source terms   !
  INTEGER                :: iter
  LOGICAL                :: break
  DOUBLE PRECISION       :: wall_time           ! wall clock elapsed time !
  DOUBLE PRECISION       :: log_time            ! time for next log output!
  DOUBLE PRECISION       :: start_time          ! system clock start time !
  DOUBLE PRECISION       :: end_time            ! system clock end time   !
  DOUBLE PRECISION       :: run_time            ! = end_time - start_time !
  INTEGER                :: start_count         ! system clock count
  TYPE(Fositestruc_TYP),DIMENSION(:), POINTER &
                         :: structure           ! structure of variables  !
#ifdef PARALLEL
  INTEGER                :: ierror
  REAL                   :: dt_all              ! min timestep of all     !
                                                ! processes               !
#endif
  CHARACTER(MAXLEN)      :: buffer

  CONTAINS
    PROCEDURE :: InitFosite
    PROCEDURE :: Setup
    PROCEDURE :: FirstStep
    PROCEDURE :: Step
    PROCEDURE :: Run
    PROCEDURE :: PrintInfo
    PROCEDURE :: PrintBoundaryFluxes
    PROCEDURE :: PrintSummary
    PROCEDURE :: ComputeRunTime
    PROCEDURE :: Finalize
  END TYPE fosite
  !--------------------------------------------------------------------------!
  PUBLIC  :: fosite
  PRIVATE :: InitFosite, Setup, Run, Step, PrintInfo, PrintBoundaryFluxes, &
             PrintSummary, FirstStep, Finalize, ComputeRunTime
  !--------------------------------------------------------------------------!


CONTAINS

  SUBROUTINE InitFosite(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    !--------------------------------------------------------------------------!
#ifdef PARALLEL
    LOGICAL            :: already_initialized
#endif
    !--------------------------------------------------------------------------!
#ifdef PARALLEL
! initialize MPI library for parallel execution, if Fosite is not initialized
    IF(.NOT.this%Initialized()) &
      CALL MPI_Initialized(already_initialized,this%ierror)
      IF (.NOT.already_initialized) &
      CALL MPI_Init(this%ierror)
#endif
    !> if Fosite is already initialized, close it, but do not finalize MPI
    IF(this%Initialized()) &
      CALL this%Finalize(.FALSE.)

    CALL this%InitLogging(simtype,simname)

    IF (this%GetRank().EQ.0) THEN
        ! print some information
        WRITE(this%buffer, "(A)")&
            "+---------------------------------------------------------+"
        CALL this%Info(this%buffer)
        WRITE(this%buffer, "(A1,A29,A28,A1)")&
            "|",TRIM(simname),"","|"
        CALL this%Info(this%buffer)
        WRITE(this%buffer, "(A1,A35,A22,A1)")&
            "|",TRIM(VERSION),"","|"
        CALL this%Info(this%buffer)
        WRITE(this%buffer, "(A)")&
            "|          Solution of 3D advection problems              |"
        CALL this%Info(this%buffer)
        WRITE(this%buffer, "(A)")&
            "+---------------------------------------------------------+"
        CALL this%Info(this%buffer)
        CALL this%Info("Initializing simulation:")
    END IF
    CALL InitDict()
    this%iter = 0
  END SUBROUTINE InitFosite

  SUBROUTINE Setup(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite),  INTENT(INOUT) :: this
    TYPE(Dict_TYP), POINTER       :: dir, IOdir, config_copy
    !--------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
      CALL this%Error("Setup","Sim is uninitialized")

    CALL DeleteDict(this%IO)

    CALL GetAttr(this%config, "mesh", dir)
    IF(ASSOCIATED(dir)) THEN
      NULLIFY(IOdir)
      CALL new_mesh(this%Mesh, dir, IOdir)
      IF(ASSOCIATED(IOdir)) CALL SetAttr(this%IO, "mesh", IOdir)
    END IF

    CALL GetAttr(this%config, "physics", dir)
    IF(ASSOCIATED(dir)) THEN
      NULLIFY(IOdir)
      CALL new_physics(this%Physics, this%Mesh, dir, IOdir)
      IF(ASSOCIATED(IOdir)) CALL SetAttr(this%IO, "physics", IOdir)
    END IF

    CALL GetAttr(this%config, "fluxes", dir)
    IF(ASSOCIATED(dir)) THEN
      NULLIFY(IOdir)
      CALL new_fluxes(this%Fluxes, this%Mesh, this%Physics, dir, IOdir)
      IF(ASSOCIATED(IOdir)) CALL SetAttr(this%IO, "fluxes", IOdir)
    END IF

    CALL GetAttr(this%config, "timedisc", dir)
    IF(ASSOCIATED(dir)) THEN
      NULLIFY(IOdir)
      CALL new_timedisc(this%Timedisc, this%Mesh,this%Physics,dir,IOdir)
      IF(ASSOCIATED(IOdir)) CALL SetAttr(this%IO, "timedisc", IOdir)
    END IF

    CALL GetAttr(this%config, "boundary", dir)
    IF(ASSOCIATED(dir)) THEN
      NULLIFY(IOdir)
      CALL new_boundary(this%Timedisc%Boundary, this%Mesh, this%Physics, dir,IOdir)
      IF(ASSOCIATED(IOdir)) CALL SetAttr(this%IO, "boundary", IOdir)
    END IF

    CALL GetAttr(this%config, "sources", dir)
    IF(ASSOCIATED(dir)) THEN
      IF(HasChild(dir)) THEN
        NULLIFY(IOdir)
        CALL new_sources(this%Sources, &
          this%Mesh,this%Fluxes,this%Physics,dir,IOdir)
        IF(ASSOCIATED(IOdir)) CALL SetAttr(this%IO, "sources", IOdir)
      END IF
    END IF

    CALL GetAttr(this%config, "datafile", dir)
    IF(ASSOCIATED(dir)) THEN
      CALL new_fileio(this%Datafile, this%Mesh, this%Physics, this%Timedisc,&
                      this%Sources, dir,this%IO)!,this%config)
    END IF

    CALL GetAttr(this%config, "logfile", dir)
    IF(ASSOCIATED(dir)) THEN
      CALL new_fileio(this%Logfile, this%Mesh, this%Physics, this%Timedisc,&
                      this%Sources,dir,this%IO)
    END IF

    CALL SetAttr(this%config, "version", TRIM(VERSION))
    CALL CopyDict(this%config, config_copy)

    CALL SetAttr(this%IO, "config", config_copy)

  END SUBROUTINE Setup


  SUBROUTINE FirstStep(this)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    INTEGER                      :: datetime(8)
    !--------------------------------------------------------------------------!

#ifdef PARALLEL
    this%start_time = MPI_Wtime()
#else
    CALL CPU_TIME(this%start_time)
#endif

    CALL SYSTEM_CLOCK(this%start_count)

    ! make sure that the initial data is written to the log file
    this%wall_time = this%start_time
    this%log_time  = this%wall_time

    IF (this%GetRank().EQ.0) THEN
        CALL this%Info( &
            "===================================================================")
        CALL this%Info("Starting calculation...")
        CALL date_and_time(values = datetime)
        WRITE(this%buffer,"(A,I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2)")&
              "Time: ", &
              datetime(3), ".", datetime(2), ".", datetime(1), " ",&
              datetime(5), ":", datetime(6), ":", datetime(7)
       CALL this%Info(this%buffer)
       CALL this%Info( &
            "step     time        n           t     min(dt) due to    adj.")
       CALL this%Info( &
            "-------------------------------------------------------------------")
    END IF

    ! determine the background velocity if fargo advection type 1 is enabled
    IF (this%Mesh%FARGO.EQ.1) THEN
       ! make sure there is valid data at least in the i-ghost cells
       CALL this%Timedisc%Boundary%CenterBoundary(this%Mesh,this%Physics,&
                             0.0,this%Timedisc%pvar%data4d,this%Timedisc%cvar%data4d)
       CALL this%Timedisc%CalcBackgroundVelocity(this%Mesh,this%Physics, &
                             this%Timedisc%pvar%data4d,this%Timedisc%cvar%data4d,this%Timedisc%w)
    END IF

    ! do a complete update of all data
    CALL this%Timedisc%ComputeRHS(this%Mesh,this%Physics,this%Sources,this%Fluxes, &
         this%Timedisc%time,0.0,this%Timedisc%pvar%data4d,this%Timedisc%cvar%data4d, &
         CHECK_ALL,this%Timedisc%rhs)

    ! calculate timestep
    this%Timedisc%dt = this%Timedisc%CalcTimestep(this%Mesh,this%Physics,this%Sources,&
                            this%Fluxes,this%Timedisc%time,this%Timedisc%dtcause)

    SELECT TYPE(phys => this%Physics)
    CLASS IS(physics_eulerisotherm)
      IF(phys%csiso.GT.0.) THEN
        IF(ANY(phys%bccsound.NE.phys%csiso)) THEN
            CALL this%Error("FirstStep","isothermal sound speed set, but "&
            // "arrays bccsound and/or fcsound have been overwritten.")
        END IF
      END IF
    CLASS DEFAULT
      ! do nothing
    END SELECT


    ! store old values
    this%Timedisc%cold(:,:,:,:) = this%Timedisc%cvar%data4d(:,:,:,:)

    ! store initial data
    IF (this%Timedisc%time.EQ.0.0) THEN
        CALL this%Datafile%WriteDataset(this%Mesh,this%Physics,this%Fluxes,&
                          this%Timedisc,this%config,this%IO)
        IF (this%GetRank().EQ.0) CALL this%PrintInfo(0,0,0.0,0.0,0,this%Timedisc%n_adj)
    END IF

    IF(this%Timedisc%break) &
      CALL this%Error("FirstStep","Initial data invalid!")
  END SUBROUTINE FirstStep


  SUBROUTINE Run(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    !--------------------------------------------------------------------------!

    ! main loop
    DO WHILE((this%Timedisc%maxiter.LE.0).OR.(this%iter.LE.this%Timedisc%maxiter))
      IF(this%Step()) EXIT
    END DO

  END SUBROUTINE Run


  FUNCTION Step(this) RESULT(break)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    LOGICAL                      :: break
#ifdef PARALLEL
    REAL,DIMENSION(2)            :: dt_buf
#endif
!    TYPE(Gravity_TYP), POINTER   :: pmass
    !--------------------------------------------------------------------------!
    break = .FALSE.

    IF (this%iter.EQ.0) &
        CALL this%FirstStep()

    ! finish simulation if stop time is reached
    IF (ABS(this%Timedisc%stoptime-this%Timedisc%time)&
            .LE.1.0E-05*this%Timedisc%stoptime) THEN
      break = .TRUE.
      RETURN
    END IF

#ifdef PARALLEL
    ! In Fortran MPI_MINLOC is only able to have two values of the same kind!
    dt_buf(1) = this%Timedisc%dt
    dt_buf(2) = this%Timedisc%dtcause
    CALL MPI_Allreduce(MPI_IN_PLACE,dt_buf,1,DEFAULT_MPI_2REAL,MPI_MINLOC,&
         this%Mesh%comm_cart,this%ierror)
    this%Timedisc%dt = dt_buf(1)
    this%Timedisc%dtcause = dt_buf(2)

    this%run_time = MPI_Wtime() - this%log_time
    CALL MPI_Allreduce(this%run_time,this%wall_time,1,MPI_DOUBLE_PRECISION,&
                       MPI_MIN,this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%wall_time)
    this%wall_time = this%wall_time - this%log_time
#endif

    ! adjust timestep for output and calculate the wall clock time
    CALL this%Datafile%AdjustTimestep(this%Timedisc%time,&
                        this%Timedisc%dt,this%Timedisc%dtcause)

#ifdef PARALLEL
    ! In Fortran MPI_MINLOC is only able to have two values of the same kind!
    dt_buf(1) = this%Timedisc%dt
    dt_buf(2) = this%Timedisc%dtcause
    CALL MPI_Allreduce(MPI_IN_PLACE,dt_buf,1,DEFAULT_MPI_2REAL,MPI_MINLOC,&
         this%Mesh%comm_cart,this%ierror)
    this%Timedisc%dt = dt_buf(1)
    this%Timedisc%dtcause = dt_buf(2)
#endif

    ! advance the solution in time
    CALL this%Timedisc%IntegrationStep(this%Mesh,this%Physics,this%Sources, &
                                       this%Fluxes,this%iter,this%config,this%IO)

    ! write output to data file
    IF ((ABS(this%Datafile%time-this%Timedisc%time)&
            .LE.1.0E-5*this%Datafile%time).OR.&
        this%Timedisc%break) THEN
       CALL this%Datafile%WriteDataset(this%Mesh,this%Physics,this%Fluxes,&
                         this%Timedisc,this%config,this%IO)
       CALL this%PrintInfo(this%Datafile%step-1, this%iter,this%Timedisc%time,&
               this%Timedisc%dtmin,this%Timedisc%dtmincause,this%Timedisc%n_adj)

       ! Stop program if a break is requested from SolveODE
       IF(this%Timedisc%break.AND.this%GetRank().EQ.0) THEN
         CALL this%Warning("SolveODE", "Time step too small, aborting.",0)
         break = .True.
         RETURN
       END IF

       ! reset dt_min,dtmincause and n_adj
       this%Timedisc%dtmin = this%Timedisc%stoptime
       this%Timedisc%dtmincause = -99
       this%Timedisc%n_adj = 0
       !IF(GetRank(this).EQ.0) &
       !  WRITE(*,"(A,ES10.4,A,ES10.4)") "dtmean: ", this%Timedisc%dtmean, " +- ",&
       !    SQRT(this%Timedisc%dtstddev/(this%Timedisc%dtaccept-1))/this%Timedisc%dtmean
       this%Timedisc%dtmean = 0.
       this%Timedisc%dtstddev = 0.
       this%Timedisc%dtaccept = 0
       IF(this%Timedisc%write_error) &
         this%Timedisc%errorval = 0.
    END IF

    ! write output to log file
!    IF (ALLOCATED(this%Logfile)) THEN
!      IF(this%wall_time.GE.this%Logfile%dtwall) THEN
!        CALL this%Logfile%WriteDataset(this%Mesh,this%Physics,this%Fluxes,&
!                          this%Timedisc,this%config,this%IO)
!        this%log_time = this%wall_time + this%log_time
!      END IF
!    END IF

    ! calculate next timestep
    this%Timedisc%dt = this%Timedisc%CalcTimestep(this%Mesh,this%Physics,this%Sources, &
                            this%Fluxes,this%Timedisc%time,this%Timedisc%dtcause)
  END FUNCTION Step

  SUBROUTINE Finalize(this,mpifinalize_)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    LOGICAL,OPTIONAL, INTENT(IN) :: mpifinalize_
    !--------------------------------------------------------------------------!
    LOGICAL                      :: mpifinalize = .TRUE.
    !--------------------------------------------------------------------------!
    mpifinalize = .TRUE.
    CALL this%ComputeRunTime()
    CALL this%PrintBoundaryFluxes()
    CALL this%PrintSummary()

    CALL this%Datafile%Finalize()
    DEALLOCATE(this%Datafile)
    IF (ALLOCATED(this%Logfile)) THEN
      CALL this%Logfile%Finalize()
      DEALLOCATE(this%Logfile)
    END IF
    CALL this%Timedisc%Finalize()
    DEALLOCATE(this%Timedisc)

    IF (ASSOCIATED(this%Sources)) THEN
      CALL this%Sources%Finalize()
      DEALLOCATE(this%Sources)
    END IF
    CALL this%Physics%Finalize()
    DEALLOCATE(this%Physics)
    CALL this%Fluxes%Finalize()
    DEALLOCATE(this%Fluxes)
    CALL this%Mesh%Finalize()
    DEALLOCATE(this%Mesh)

    CALL DeleteDict(this%IO)
    CALL DeleteDict(this%config)

#ifdef PARALLEL
    IF(PRESENT(mpifinalize_)) THEN
      mpifinalize = mpifinalize_
    END IF
    IF(mpifinalize) THEN
      CALL MPI_Finalize(this%ierror)
    END IF
#endif
  END SUBROUTINE Finalize

  SUBROUTINE ComputeRunTime(this)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT)   :: this
    !--------------------------------------------------------------------------!

#ifdef PARALLEL
    this%end_time = MPI_Wtime() - this%start_time
    CALL MPI_Allreduce(this%end_time,this%run_time,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
        this%Mesh%comm_cart,this%ierror)
#else
    CALL CPU_TIME(this%end_time)
    this%run_time = this%end_time - this%start_time
#endif

  END SUBROUTINE ComputeRunTime


  SUBROUTINE PrintInfo(this,step,i,t,d,dc,na)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    INTEGER            :: step, i, dc, na
    CHARACTER(LEN=9)   :: dtcause
    REAL               :: t, d
    !------------------------------------------------------------------------!
    INTEGER            :: drt, c, c_rate, c_max
    !------------------------------------------------------------------------!
    INTENT(IN)         :: i,t,d
    !--------------------------------------------------------------------------!
    IF (this%GetRank().EQ.0) THEN

       SELECT CASE (dc)
        ! positive values represent source terms
        CASE(1:)
          WRITE(dtcause, "(A,I2.2,A)") " S",dc," "
        CASE(DTCAUSE_CFL)
          WRITE(dtcause, "(A)") " cfl "
        CASE(DTCAUSE_ERRADJ)
           WRITE(dtcause, "(A)") " err_adj "
        CASE(DTCAUSE_SMALLERR)
           WRITE(dtcause, "(A)") " err "
        CASE(DTCAUSE_FILEIO)
           ! output by this reason is suppressed by default
           WRITE(dtcause, "(A)") " fileio "
        CASE DEFAULT
           WRITE(dtcause, "(A,I3.2,A)") " ?", dc, "? "
       END SELECT

       CALL SYSTEM_CLOCK(c,c_rate, c_max)
       ! overflow every 24days => assumption: max one per simulation
       drt = (c - this%start_count)/c_rate
       IF (drt .LT. 0) THEN
         drt = c_max/c_rate + drt
       END IF
       WRITE(this%buffer,"(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I8,A,ES11.3,A,ES11.3,A,I5)")&
             step, " ", drt/3600, ":", mod(drt,3600)/60, ":", mod(drt,60),&
             " ", i, " ", t, " ", d, dtcause, na
       CALL this%Info(this%buffer)
    END IF
  END SUBROUTINE PrintInfo

  SUBROUTINE PrintBoundaryFluxes(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT)         :: this
    INTEGER                              :: k
    REAL, DIMENSION(this%Physics%VNUM,6) :: bflux
    !--------------------------------------------------------------------------!
    DO k=1,6
       bflux(:,k) = this%Fluxes%GetBoundaryFlux(this%Mesh,this%Physics,k)
    END DO
    IF (this%GetRank().EQ.0) THEN
       CALL this%Info("-------------------------------------------------------------------")
       CALL this%Info("total boundary fluxes:")
       CALL this%Info("                      west        east        south")
       DO k=1,this%Physics%VNUM
          WRITE(this%buffer,"(T2,A,T21,3(ES12.3))")TRIM(this%Physics%cvarname(k)), &
               bflux(k,WEST), bflux(k,EAST), bflux(k,SOUTH)
          CALL this%Info(this%buffer)
       END DO
       CALL this%Info("                      north       bottom      top")
       DO k=1,this%Physics%VNUM
          WRITE(this%buffer,"(T2,A,T21,3(ES12.3))")TRIM(this%Physics%cvarname(k)), &
               bflux(k,NORTH), bflux(k,BOTTOM), bflux(k,TOP)
          CALL this%Info(this%buffer)
       END DO
    END IF
  END SUBROUTINE PrintBoundaryFluxes

  SUBROUTINE PrintSummary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite), INTENT(INOUT) :: this
    !--------------------------------------------------------------------------!
    IF (this%GetRank().EQ.0) THEN
       CALL this%Info("===================================================================")
       IF ((this%Timedisc%maxiter.LE.0).OR.(this%iter.LT.this%Timedisc%maxiter)) THEN
          CALL this%Info("calculation finished correctly.")
       ELSE
          CALL this%Info("too many iterations, aborting!")
       END IF
       WRITE(this%buffer,"(A,F10.2,A)")" main loop runtime: ", this%run_time, " sec."
       CALL this%Info(this%buffer)
    END IF
  END SUBROUTINE PrintSummary

END MODULE fosite_mod
