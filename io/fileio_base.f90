!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_base.f90                                                   #
!#                                                                           #
!# Copyright (C) 2008-2023                                                   #
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
!> \addtogroup fileio
!! \brief Handles all Input/Output
!! - general parameters of fileio group as key-values
!! \key{count,INTEGER,number of output steps,1}
!! \key{filecycles,INTEGER,number of data files (=0 => one file and append data),count+1}
!! \key{stoptime,REAL,stop time for output,Timedisc%stoptime}
!! \key{fileformat,INTEGER,type of fileio}
!! \key{filepath,CHARACTER,file path,""}
!! \key{filename,CHARACTER,file name}
!! \key{dtwall,INTEGER,wall clock time between successive outputs,3600}
!! \key{unit,INTEGER,fortran unit number for I/O,lastunit+1}
#ifdef HAVE_NETCDF
!! \key{ncfmt,INTEGER,netcdf format type}
#endif
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief Generic file I/O module
!!
!! This module provides the generic interface routines to all file I/O
!! modules.
!!
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_base_mod
  USE logging_base_mod
  USE common_dict
  USE fluxes_base_mod
  USE physics_base_mod
  USE sources_base_mod
  USE mesh_base_mod
  USE timedisc_base_mod
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
  !--------------------------------------------------------------------------!
  ! Private Attributes section starts here:
  !> \name Private Attributes
  !>#### file name and extension lengths
  INTEGER, PARAMETER :: FEXTLEN = 4                  !< file name extension length
  INTEGER, PARAMETER :: FMLTLEN = 5                  !< length of multi process string (parallel mode only)
  INTEGER, PARAMETER :: FCYCLEN = 5                  !< length of timestep string
  INTEGER, PARAMETER :: FNAMLEN = 256                !< file name length (without any extension)
  INTEGER, PARAMETER :: FPATLEN = 1024               !< file path length (without file name)
  !> \name
  !!#### handling multiple files in parallel mode
#ifdef PARALLEL
  CHARACTER(LEN=FMLTLEN), SAVE :: fmextstr = ""      !< multi process string, overwritten below
#endif
  INTEGER, PARAMETER           :: MAXMLTFILES = 1000 !< max. number files per time step (parallel
                                                     !! mode with one file per node)
  !> \name
  !!#### handling multiple data files with time step in their names
  INTEGER, PARAMETER           :: MAXCYCLES = 10000  !< max. number of data files (not counting
                                                     !! multiple files per time step in parallel mode)
  CHARACTER(LEN=32), SAVE      :: cycfmt             !< format string for cycles
  !--------------------------------------------------------------------------!
  !> output-pointer for array data (binary,gnuplot,vtk)
  TYPE ValPtr_TYP
    REAL, DIMENSION(:,:,:), POINTER :: val
  END TYPE

  TYPE Output_TYP
    TYPE(ValPtr_TYP), DIMENSION(:), POINTER :: p
    CHARACTER(LEN=128)                      :: key
    CHARACTER(LEN=1024)                     :: path
    INTEGER(KIND=4)                         :: numbytes
  END TYPE Output_TYP

  !> output-pointer for time step scalar data (gnuplot)
  TYPE TSOutput_TYP
    REAL, POINTER           :: val
    CHARACTER(LEN=128)      :: key
  END TYPE TSOutput_TYP

  !> class for Fortran file handle
  TYPE, EXTENDS(logging_base) :: filehandle_fortran
    !> \name Variables
    INTEGER   :: fid                       !< unique ID for file access (Fortran i/o unit)
    LOGICAL   :: textfile                  !< true for text, i.e. ascii stream
    CHARACTER(LEN=FNAMLEN) :: filename     !< file name without extension
    CHARACTER(LEN=FPATLEN) :: path         !< file path without filename
    CHARACTER(LEN=FEXTLEN) :: extension    !< file name extension
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFilehandle
    PROCEDURE :: OpenFile
    PROCEDURE :: CloseFile
    PROCEDURE :: GetFilename
    PROCEDURE :: GetBasename
    PROCEDURE :: GetUnitNumber
    PROCEDURE :: GetFormat
    FINAL :: Finalize_fortran
  END TYPE filehandle_fortran

  !> class for MPI file handle
  TYPE, EXTENDS(filehandle_fortran) :: filehandle_mpi
  END TYPE filehandle_mpi

  !> FileIO base class
  TYPE, ABSTRACT,EXTENDS(logging_base) :: fileio_base
     !> \name Variables
     CLASS(filehandle_fortran), ALLOCATABLE &
                            :: datafile    !< file handle for data file
     CHARACTER(LEN=512)     :: linebuf     !< char buffer fo field data
     CHARACTER(LEN=512)     :: tslinebuf   !< char buffer for time step data
     LOGICAL                :: cartcoords  !< output cartesian coordinates
     INTEGER                :: unit        !< i/o unit
     INTEGER                :: step        !< counter for output steps
     INTEGER                :: count       !< number of output steps
     INTEGER                :: cycles      !< number of output files
     INTEGER                :: dtwall      !< wall clock time difference
     INTEGER                :: inum,jnum,& !< mesh extent
                               knum
#ifndef PARALLEL
     INTEGER(KIND=8)        :: offset      !< header offset (MPI see below), must be 8 byte for VTK
#endif
     REAL                   :: stoptime    !< final simulation time for data output
     REAL                   :: starttime   !< initial simulation time for data output
     REAL                   :: time        !< output time
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                               binout      !< binary data output buffer
     TYPE(Output_TYP),DIMENSION(:), POINTER :: &
                               output      !< list of output fields
     TYPE(TSOutput_TYP),DIMENSION(:), POINTER :: &
                               tsoutput    !< list of scalar time step output
     REAL, DIMENSION(:,:) , POINTER  :: &
                               bflux       !< boundary flux output buffer
     REAL                   :: walltime    !< adds output before walltime is
#ifdef PARALLEL
     !> \name Variables in Parallel Mode
     LOGICAL                :: multfiles   !< spread files across nodes
     INTEGER                :: handle      !< MPI file handle
     INTEGER                :: bufsize     !< output data buffer size
     INTEGER                :: cbufsize    !< corner output data buffer size
     INTEGER                :: basictype   !< data type for points
     INTEGER                :: filetype    !< data type for data i/o
     INTEGER, DIMENSION(MPI_STATUS_SIZE) :: &
                               status      !< MPI i/o status record
     INTEGER(KIND=MPI_OFFSET_KIND) :: offset !< skip header bytes
#endif
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFileio
    PROCEDURE (WriteHeader), DEFERRED :: WriteHeader
    PROCEDURE (WriteDataset), DEFERRED:: WriteDataset
    PROCEDURE :: AdjustTimestep
    PROCEDURE (Finalize), DEFERRED          :: Finalize
    PROCEDURE :: Finalize_base
!     PROCEDURE :: GetFilename
!     PROCEDURE :: MakeMultstr
    PROCEDURE :: IncTime
  END TYPE fileio_base

  ! Interfaces
  ABSTRACT INTERFACE
    SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
      IMPORT fileio_base,mesh_base,physics_base,fluxes_base,timedisc_base,Dict_TYP
      IMPLICIT NONE
      CLASS(fileio_base),   INTENT(INOUT) :: this
      CLASS(mesh_base),     INTENT(IN)    :: Mesh
      CLASS(physics_base),  INTENT(INOUT) :: Physics
      CLASS(fluxes_base),   INTENT(IN)    :: Fluxes
      CLASS(timedisc_base), INTENT(IN)    :: Timedisc
      TYPE(Dict_TYP), POINTER             :: Header,IO
    END SUBROUTINE
    SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
      IMPORT fileio_base,mesh_base,physics_base,Dict_TYP
      IMPLICIT NONE
      CLASS(fileio_base),  INTENT(INOUT)  :: this      !< \param [in,out] this fileio type
      CLASS(mesh_base),    INTENT(IN)     :: Mesh      !< \param [in] Mesh mesh type
      CLASS(physics_base), INTENT(IN)     :: Physics   !< \param [in] Physics physics type
      TYPE(Dict_TYP), POINTER             :: Header,IO
    END SUBROUTINE
    SUBROUTINE Finalize(this)
      IMPORT fileio_base
      IMPLICIT NONE
      CLASS(fileio_base),  INTENT(INOUT)  :: this
    END SUBROUTINE
  END INTERFACE
  !--------------------------------------------------------------------------!

  !> \name Public Attributes
  !! #### file status and access modes
!  INTEGER, PARAMETER :: FILE_EXISTS = B'00000001' !< file status for existing files
  INTEGER, PARAMETER :: READONLY = 1       !< readonly access
  INTEGER, PARAMETER :: READEND  = 2       !< readonly access at end
  INTEGER, PARAMETER :: REPLACE  = 3       !< read/write access replacing file
  INTEGER, PARAMETER :: APPEND   = 4       !< read/write access at end
  !> \}

  !> basic file formats
  INTEGER, PARAMETER :: FORTRANFILE = 1
  INTEGER, PARAMETER :: MPIFILE = 2
  !> file I/O types
  INTEGER, PARAMETER :: BINARY  = 1
  INTEGER, PARAMETER :: GNUPLOT = 2
  INTEGER, PARAMETER :: VTK     = 4
!  INTEGER, PARAMETER :: NPY     = 5
!  INTEGER, PARAMETER :: HDF     = 6
  INTEGER, PARAMETER :: XDMF    = 7
  !--------------------------------------------------------------------------!
  INTEGER, SAVE :: lastunit = 10

  INTERFACE filehandle_fortran
    MODULE PROCEDURE CreateFilehandle
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fileio_base, Output_TYP, TSOutput_TYP, ValPtr_TYP, &
       filehandle_fortran, &
       ! constants
       BINARY, GNUPLOT, VTK, XDMF, &
#ifdef PARALLEL
       filehandle_mpi, &
       DEFAULT_MPI_REAL,&
#endif
       READONLY, READEND, REPLACE, APPEND
!       FILE_EXISTS, &


CONTAINS


  !> constructor for Fortran file handle
  FUNCTION CreateFilehandle(filename,path,extension,textfile) RESULT(new_fh)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(filehandle_fortran) :: new_fh
    CHARACTER(LEN=*), INTENT(IN) :: filename !< \param [in] filename file name without extension
    CHARACTER(LEN=*), INTENT(IN) :: path     !< \param [in] path file path without filename
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: extension !< \param [in] extension file name extension
    LOGICAL, OPTIONAL, INTENT(IN)      :: textfile !< \param [in] textfile true for text data
    !-------------------------------------------------------------------!
    CALL new_fh%InitFilehandle(filename,path,extension,textfile)
  END FUNCTION CreateFilehandle

  !> basic initialization of Fortran file handle
  SUBROUTINE InitFilehandle(this,filename,path,extension,textfile)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(INOUT) :: this !< \param [inout] this file handle class
    CHARACTER(LEN=*), INTENT(IN) :: filename !< \param [in] filename file name without extension
    CHARACTER(LEN=*), INTENT(IN) :: path     !< \param [in] path file path without filename
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: extension !< \param [in] extension file name extension
    LOGICAL, OPTIONAL, INTENT(IN)      :: textfile !< \param [in] textfile true for text data
    !-------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) CALL this%InitLogging(FORTRANFILE,"fortran")
    ! check file name length
    IF (LEN_TRIM(filename).GT.FNAMLEN) &
       CALL this%Error("fileio_base::InitFilehandle","file name too long")
    this%filename = filename
    ! check file path length
    IF (LEN_TRIM(path).GT.FPATLEN) &
       CALL this%Error("fileio_base::InitFilehandle","file path too long")
    this%path = path
    ! check file name extension
    IF (PRESENT(extension)) THEN
      IF (LEN_TRIM(extension).GT.FEXTLEN-1) &
        CALL this%Error("fileio_base::InitFilehandle","file name extension too long")
      this%extension = extension
    ELSE
      this%extension = "dat"
    END IF
    IF (PRESENT(textfile)) THEN
      this%textfile = textfile
    ELSE
      this%textfile = .TRUE. ! default is textfile
    END IF
    ! format string for writing file names with explicit time step
    WRITE (cycfmt, "('(A,I',I1,'.',I1,',A)')") FCYCLEN-1,FCYCLEN-1
  END SUBROUTINE InitFilehandle
 
  !> \public Generic constructor for file I/O
  !!
  SUBROUTINE InitFileIO(this,Mesh,Physics,Timedisc,Sources,config,IO,fmtname,fext)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base),  INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),    INTENT(IN)    :: Mesh     !< \param [in] Mesh mesh type
    CLASS(physics_base), INTENT(IN)    :: Physics  !< \param [in] Physics physics type
    CLASS(timedisc_base),INTENT(IN)    :: Timedisc !< \param [in] Timedisc timedisc type
    CLASS(sources_base), INTENT(IN), POINTER :: Sources !< \param [in] Sources sources type
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: config  !< \param [in] config dict with I/O configuration
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: IO      !< \param [in] IO dict with pointers to I/O arrays
    CHARACTER(LEN=*),    INTENT(IN)    :: fmtname  !< format name
    CHARACTER(LEN=*),    INTENT(IN)    :: fext     !< fext file name extension
    !------------------------------------------------------------------------!
    INTEGER                        :: fileformat !< fmt fileio type number
    CHARACTER(LEN=MAX_CHAR_LEN)    :: fname      !< fname file name
    CHARACTER(LEN=MAX_CHAR_LEN)    :: fpath      !< fpath file path
    !INTEGER                        :: fcycles    !< fcycles number of file cycles
    INTEGER                        :: unit       !<  unit fortran i/o unit number
!    LOGICAL                        :: success
    CHARACTER(LEN=32)              :: timestamp
    INTEGER                        :: count_def, fcycles_def, dtwall_def
    INTEGER                        :: cartcoords
    REAL                           :: stoptime_def
    REAL                           :: time
    TYPE(Dict_TYP),POINTER         :: oldconfig => null()
    !------------------------------------------------------------------------!
    ! wall clock time between successive outputs
    ! this is mainly intended for log file outputs
    CALL GetAttr(config, "dtwall" , dtwall_def, 3600) !default is one hour

    ! number of output steps
    CALL GetAttr(config, "count" , count_def, 1)

    ! number of data files
    ! fcycles = 0     : one data file, append data
    ! fcycles = 1     : one data file, overwrite data
    ! fcycles = X > 1 : X data files
    CALL GetAttr(config, "filecycles", fcycles_def, count_def+1)

    ! stop time for output defaults to simulation stop time
    CALL GetAttr(config, "stoptime"  , stoptime_def, Timedisc%stoptime)
    CALL GetAttr(config, "fileformat", fileformat)
    CALL GetAttr(config, "filename"  , fname)
    fpath = ""
    CALL GetAttr(config, "filepath"  , fpath, fpath)
    IF (.NOT.this%datafile%Initialized()) THEN
      CALL this%datafile%InitFilehandle(fname,fpath,fext)
    END IF
    CALL GetAttr(config, "unit"      , unit , lastunit+1)
    lastunit = unit
    CALL GetAttr(config, "walltime"  , this%walltime, HUGE(1.0))

    ! mesh coordinates are cartesian by default (cartcoords == 1)
    ! set to 0 for curvilinear coordinates (currently supported in gnuplot and vtk)
    CALL GetAttr(config, "cartcoords", cartcoords, 1)
    IF (cartcoords.EQ.0) THEN
      this%cartcoords = .FALSE.
    ELSE
      this%cartcoords = .TRUE.
    END IF

    ! InitLogging after reading fileformat from config
    CALL this%InitLogging(fileformat,fmtname)

    ! check cycles
    IF (fcycles_def.GT.MAXCYCLES) &
       CALL this%Error("InitFileIO","file cycles exceed limits")

#ifdef PARALLEL
    ! turn on multiple file output if requested
    IF (this%multfiles) THEN
       ! check number of parallel processes
       IF (this%GetNumProcs().GT.MAXMLTFILES) &
          CALL this%Error("InitFileIO","number of processes for multiple file output exceeds limits")
       fmextstr = this%MakeMultstr()
    END IF
#endif

    this%cycles    = fcycles_def
    this%unit      = unit
    this%stoptime  = stoptime_def
    this%starttime = Timedisc%time  ! from beginning of simulation
    this%dtwall    = dtwall_def
    this%count     = count_def
    this%offset    = 0
    this%err       = 0

    CALL Getattr(config, "step", this%step, 0)

    ! compute the (actual) output time
    this%time = Timedisc%time

    ! print some information
    CALL this%Info(" FILEIO---> file type:         " // TRIM(this%GetName()))
 !   IF (success) THEN
       WRITE (timestamp,'(ES10.4)') Timedisc%time
       CALL this%Info("            time stamp:        " // TRIM(timestamp))
 !   END IF

    ! time for next output
    IF (Timedisc%time.GT.0.0) CALL this%IncTime()
    IF (ASSOCIATED(oldconfig)) CALL DeleteDict(oldconfig)
!!!!! old fosite code, may be obsolete
! #ifdef PARALLEL
!     ! check data type extents in files
!     ! first try to create a new dummy file
!     CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(IOR(MPI_MODE_RDWR,&
!          MPI_MODE_CREATE),IOR(MPI_MODE_EXCL,MPI_MODE_DELETE_ON_CLOSE)),&
!          MPI_INFO_NULL,this%handle,this%error)
!     ! maybe file exists
!     IF (this%error.NE.0) CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),&
!          MPI_MODE_RDONLY,MPI_INFO_NULL,this%handle,this%error)
!     ! then check the data type sizes
!     ! extent of integer in file
!     IF (this%error.EQ.0) CALL MPI_File_get_type_extent(this%handle,MPI_INTEGER,&
!          this%intext,this%error)
!     ! extent of real in file
!     IF (this%error.EQ.0) CALL MPI_File_get_type_extent(this%handle,DEFAULT_MPI_REAL,&
!          this%realext,this%error)
!     IF (this%error.NE.0) CALL Error(this,"InitFileIO","unable to check file properties")
!     CALL MPI_File_close(this%handle,this%error)
! #endif
  END SUBROUTINE InitFileIO

  !> \public Increments the counter for timesteps and sets the time for next output
  !!
  PURE SUBROUTINE IncTime(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(INOUT) :: this!< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    this%time = this%time + ABS(this%stoptime - this%starttime) / this%count
    this%step = this%step + 1
  END SUBROUTINE IncTime


  !> \public Adjust the current timestep
  !!
  !! Last timestep before output must fit to desired time for output.
 PURE SUBROUTINE AdjustTimestep(this,time,dt,dtcause)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(IN) :: this    !< \param [in] this fileio type
    REAL                           :: time    !< \param [in,out] time
    REAL                           :: dt      !< \param [in,out] dt timestep
    INTEGER                        :: dtcause !< \param [in,out] dtcause cause of smallest dt
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: time,dt,dtcause
    !------------------------------------------------------------------------!
    IF ((time+dt)/this%time.GT.1.0) THEN
       dt = this%time - time
       dtcause = DTCAUSE_FILEIO
    ELSE IF((time+1.5*dt)/this%time.GT.1.0) THEN
       dt = 0.5*(this%time - time)
       dtcause = DTCAUSE_FILEIO
    END IF
  END SUBROUTINE AdjustTimestep


!   !> \public Get a file label (multiples files in parallel mode)
!   !! without filenumber => use GetRank; with filenumber < 0 => empty label
!   !! \result file label
!   FUNCTION MakeMultstr(this,fn) RESULT (multstr)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(filehandle_fortran), INTENT(IN) :: this !< \param [in,out] this fileio type
!     INTEGER, OPTIONAL,  INTENT(IN) :: fn       !< \param [in] fn number of file
!     CHARACTER(LEN=FMLTLEN)         :: multstr
!     !------------------------------------------------------------------------!
! #ifdef PARALLEL
!     INTEGER                        :: fn_l
!     CHARACTER(LEN=32)              :: mextfmt
! #endif
!     !------------------------------------------------------------------------!
! #ifdef PARALLEL
!     IF (PRESENT(fn)) THEN
!       fn_l = fn
!     ELSE
!       fn_l = this%GetRank()
!     END IF
!     ! with fn < 0 you can suppress a label
!     IF (this%multfiles .AND. fn_l .GE. 0) THEN
!         WRITE (mextfmt, "('(A2,I',I1,'.',I1,')')") FMLTLEN-2,FMLTLEN-2
! 
!         WRITE (multstr, mextfmt) "-r", fn_l
!     ELSE
!       multstr = ""
!     END IF
! #else
!     multstr = ""
! #endif
!   END FUNCTION MakeMultstr


  !> \public Get the current file name without path
  !! e.g. important for vtk (pvts files)
  !! \result current file name
  FUNCTION GetBasename(this,step,cycles) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(IN) :: this   !< \param [in] this file handle
    INTEGER, INTENT(IN)                   :: step   !< \param [in] step time step
    INTEGER, OPTIONAL, INTENT(IN)         :: cycles !< \param [in] step number of file cycles
    CHARACTER(LEN=256)                    :: fname
    !------------------------------------------------------------------------!
    INTEGER :: fnum
    !------------------------------------------------------------------------!
    fnum = -1
    IF (PRESENT(cycles)) THEN
      IF (cycles.GT.0) &
          ! determine file number based on current step and number of files,
          ! i.e. cycles, at all
        fnum = MODULO(step,cycles)
    ELSE
      ! if cycles are not given then fnum is given by the current step
      fnum = step
    END IF
    IF (fnum.GT.-1) THEN
      ! generate a file name with fnum appended
      WRITE (fname, FMT=TRIM(cycfmt)) TRIM(this%filename) // "_", &
         fnum, "." // TRIM(this%extension)
    ELSE
      ! do not append any number to the basic file name, just the extension
      WRITE (fname,"(A)") TRIM(this%filename) // "." // TRIM(this%extension)
    END IF
  END FUNCTION GetBasename

!   !> \public Get the current file name
!   !! \result current file name
!   FUNCTION GetFilename(this,fn) RESULT (fname)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(fileio_base), INTENT(IN) :: this !< \param [in] this fileio type
!     INTEGER, OPTIONAL, INTENT(IN)  :: fn   !< \param [in] fn number of file
!     CHARACTER(LEN=256)             :: fname
!     !------------------------------------------------------------------------!
!     fname = TRIM(this%path) // TRIM(GetBasename(this,fn))
!   END FUNCTION GetFilename

  !> \public Generic deconstructor of the file I/O
  !!
  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base),INTENT(INOUT) :: this !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("fileio_base::Finalize_base","FileIO not initialized")
  END SUBROUTINE Finalize_base

  !> \public Open file to access Fortran stream
  !!
  SUBROUTINE OpenFile(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(INOUT):: this   !< \param [in,out] this fileio type
    INTEGER,           INTENT(IN)           :: action !< \param [in] action mode of file access
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)  :: sta,act,pos,frm
    !------------------------------------------------------------------------!
    SELECT CASE(action)
    CASE(READONLY)
       sta = 'OLD'
       act = 'READ'
       pos = 'REWIND'
    CASE(READEND)
       sta = 'OLD'
       act = 'READ'
       pos = 'APPEND'
    CASE(REPLACE)
       sta = 'REPLACE'
       act = 'WRITE'
       pos = 'REWIND'
    CASE(APPEND)
       sta = 'OLD'
       act = 'READWRITE'
       pos = 'APPEND'
    CASE DEFAULT
       CALL this%Error("fileio_base::OpenFile_serial","Unknown access mode.")
    END SELECT
    OPEN(UNIT=this%GetUnitNumber(),FILE=this%GetFilename(),STATUS=TRIM(sta), &
         ACCESS='STREAM',ACTION=TRIM(act),POSITION=TRIM(pos),FORM=this%GetFormat(), &
         IOSTAT=this%err)
  END SUBROUTINE OpenFile
 
  !> \public Open file for input/ouput in parallel mode
  !!
!   SUBROUTINE OpenFile_parallel(this,action)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(filehandle_mpi), INTENT(INOUT)    :: this   !< \param [in,out] this fileio type
!     INTEGER,           INTENT(IN)           :: action !< \param [in] action mode of file access
!     !------------------------------------------------------------------------!
! #ifdef PARALLEL
!     INTEGER(KIND=MPI_OFFSET_KIND) :: offset!< \param [in] offset offset for MPI
! #endif
!     !------------------------------------------------------------------------!
!     SELECT CASE(action)
! #ifdef PARALLEL
!     CASE(READONLY)
!        CALL MPI_File_open(MPI_COMM_WORLD,this%GetFilename(),MPI_MODE_RDONLY, &
!             MPI_INFO_NULL,this%fid,this%err)
!        this%offset = 0
!        CALL MPI_File_seek(this%fid,this%offset,MPI_SEEK_SET,this%err)
!     CASE(READEND)
!        CALL MPI_File_open(MPI_COMM_WORLD,this%GetFilename(),IOR(MPI_MODE_RDONLY,&
!             MPI_MODE_APPEND),MPI_INFO_NULL,this%fid,this%error)
!        ! opening in append mode doesn't seem to work for pvfs2, hence ...
!        offset = 0
!        CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
!        CALL MPI_File_sync(this%handle,this%error)
!     CASE(REPLACE)
!        CALL MPI_File_delete(GetFilename(this),MPI_INFO_NULL,this%error)
!        CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_WRONLY,&
!             MPI_MODE_CREATE),MPI_INFO_NULL,this%handle,this%error)
!     CASE(APPEND)
!        CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_RDWR,&
!             MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)
!        ! opening in append mode doesn't seem to work for pvfs2, hence ...
!        offset = 0
!        CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
!        CALL MPI_File_sync(this%handle,this%error)
! #endif
!     CASE DEFAULT
!        ! abort with error if either the access mode is wrong or someone tries parallel i/o in serial mode
!        CALL this%Error("fileio_base::OpenFile_parallel","Unknown access mode.")
!     END SELECT
!   END SUBROUTINE OpenFile_parallel

  !> \public get Fortran i/o unit number
  FUNCTION GetUnitNumber(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTEGER :: GetUnitNumber
    !------------------------------------------------------------------------!
    GetUnitNumber = this%fid
  END FUNCTION GetUnitNumber

  !> \public get Fortran i/o format string
  FUNCTION GetFormat(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    CHARACTER(LEN=16) :: GetFormat
    !------------------------------------------------------------------------!
    IF (this%textfile) THEN
      GetFormat = "FORMATTED"
    ELSE
      GetFormat = "UNFORMATTED"
    END IF
  END FUNCTION GetFormat

  !> \public get file name of Fortran stream
  FUNCTION GetFilename(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    CHARACTER(LEN=256) :: GetFilename
    !------------------------------------------------------------------------!
!     GetFilename = TRIM(this%path) !// TRIM(GetBasename(this,fn))
  END FUNCTION GetFilename

  !> \public close Fortran stream
  !!
  SUBROUTINE CloseFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_fortran), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
! #ifdef PARALLEL
!     CALL MPI_File_close(this%handle,this%err)
! #else
    CLOSE(UNIT=this%GetUnitNumber(),IOSTAT=this%err)
! #endif
  END SUBROUTINE CloseFile

  !> \public destructor of Fortran stream handle
  SUBROUTINE Finalize_fortran(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(filehandle_fortran), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
  END SUBROUTINE Finalize_fortran
  
END MODULE fileio_base_mod
