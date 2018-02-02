!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_generic.f90                                                #
!#                                                                           #
!# Copyright (C) 2008-2014                                                   #
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
  CHARACTER(LEN=FMLTLEN), SAVE :: fmextstr = ""      !< multi process string, overwritten below
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
    REAL, DIMENSION(:,:), POINTER           :: val
    TYPE(ValPtr_TYP), DIMENSION(:), POINTER :: p
    CHARACTER(LEN=128)                      :: key
    INTEGER(KIND=4)                         :: numbytes
  END TYPE Output_TYP

  !> output-pointer for time step scalar data (gnuplot)
  TYPE TSOutput_TYP
    REAL, POINTER           :: val
    CHARACTER(LEN=128)      :: key
  END TYPE TSOutput_TYP

  !> FileIO class
  TYPE, ABSTRACT,EXTENDS(logging_base) :: fileio_base
     !> \name Variables
     TYPE(logging_base)     :: format      !< i/o file format
     CHARACTER(LEN=512)     :: linebuf     !< char buffer fo field data
     CHARACTER(LEN=512)     :: tslinebuf   !< char buffer for time step data
     CHARACTER(LEN=FNAMLEN) :: filename    !< file name without extension
     CHARACTER(LEN=FPATLEN) :: path        !< file path without filename
     CHARACTER(LEN=FEXTLEN) :: extension   !< file name extension
     LOGICAL                :: cartcoords  !< output cartesian coordinates
     INTEGER                :: unit        !< i/o unit
     INTEGER                :: error_io    !< i/o error code
     INTEGER                :: step        !< counter for output steps
     INTEGER                :: count       !< number of output steps
     INTEGER                :: cycles      !< number of output files
     INTEGER                :: dtwall      !< wall clock time difference
     INTEGER                :: realsize    !< byte size of real numbers
     INTEGER                :: intsize     !< byte size of integer numbers
     INTEGER                :: inum,jnum,& !< mesh extent
                               knum
#ifndef PARALLEL
     INTEGER(KIND=8)        :: offset      !< header offset (MPI see below), must be 8 byte for VTK
#endif
     REAL                   :: stoptime    !< final simulation time for data output
     REAL                   :: time        !< output time
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                               binout      !< binary data output buffer
     !> \name
     !!#### VTK specific variables
     CHARACTER(LEN=32)      :: buf         !< buffer for character i/o
     CHARACTER(LEN=12)      :: realfmt     !< real format string
     CHARACTER(LEN=14)      :: endianness  !< endianness string
     CHARACTER(LEN=14)      :: endian      !< endianness string
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                                vtkcoords   !< VTK coordinate data
     !> \name
     !!#### GNUPLOT specific variables
     CHARACTER(LEN=512)     :: heading     !< char buffer for heading (field data)
     CHARACTER(LEN=512)     :: tsheading   !< char buffer for heading (time step data)
     CHARACTER(LEN=64)      :: fmtstr      !< format string
     CHARACTER(LEN=64)      :: linefmt     !< output line format string
     INTEGER                :: COLS        !< number of output columns
     INTEGER                :: TSCOLS      !< number of output columns for time step data
     INTEGER                :: MAXCOLS     !< upper limit for output cols
                                           !! (MAXCOLS < LEN(linebuf)/FLEN)
     INTEGER                :: DECS        !< decimal places for real number output
     INTEGER                :: FLEN        !< output field length
     INTEGER                :: linelen     !< length of a line
     INTEGER                :: tslinelen   !< length of a line for time step data
#ifdef HAVE_HDF5_MOD
     !> \name
     !!#### HDF specific variables
     INTEGER(HID_T)         :: fid         !< file id
     INTEGER(HID_T)         :: xferid      !< xfer id
#endif
     TYPE(Output_TYP),DIMENSION(:), POINTER :: &
                               output      !< list of output fields
     TYPE(TSOutput_TYP),DIMENSION(:), POINTER :: &
                               tsoutput    !< list of scalar time step output
     REAL, DIMENSION(:,:) , POINTER  :: &
                               bflux       !< boundary flux output buffer
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
     !> \name
     !!#### GNUPLOT specific variables
     INTEGER                :: blocknum    !< number of output blocks
     INTEGER(KIND=MPI_ADDRESS_KIND) :: &
                               realext, &  !< real data type extent
                               intext      !< integer data type extent
     CHARACTER, DIMENSION(:,:), POINTER :: &
                               outbuf      !< output buffer
     !> \name
     !!#### BINARY specific variables
     LOGICAL                :: first       !< true if this is the first output
     !> \name
     !!#### VTK specific variables
     CHARACTER(LEN=64), DIMENSION(:), POINTER :: &
                               extent      !< extent of all pieces in VTK
     !!#### BINARY specific variables
     INTEGER                :: cfiletype   !< file data type for corner i/o
     INTEGER, DIMENSION(:), POINTER :: &
                               disp        !< array of displacements
#endif
  CONTAINS
    ! methods
    PROCEDURE :: InitFileio
    PROCEDURE (WriteHeader), DEFERRED :: WriteHeader
    PROCEDURE (WriteDataset), DEFERRED:: WriteDataset
    PROCEDURE :: AdjustTimestep
    PROCEDURE :: FinalizeFileio
    PROCEDURE :: GetFilename
    PROCEDURE :: GetBasename
    PROCEDURE :: MakeMultstr
    PROCEDURE :: IncTime
  END TYPE fileio_base

  ! Interfaces
  ABSTRACT INTERFACE
    SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
      IMPORT fileio_base,mesh_base,physics_base,fluxes_base,timedisc_base,Dict_TYP
      IMPLICIT NONE
      CLASS(fileio_base),   INTENT(INOUT) :: this
      CLASS(mesh_base),     INTENT(IN)    :: Mesh
      CLASS(physics_base),  INTENT(IN)    :: Physics
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
  END INTERFACE
  !--------------------------------------------------------------------------!

  !> \name Public Attributes
  !! #### file status and access modes
!  INTEGER, PARAMETER :: FILE_EXISTS = B'00000001' !< file status for existing files
  INTEGER, PARAMETER :: READONLY = 1       !< readonly access
  INTEGER, PARAMETER :: READEND  = 2       !< readonly access at end
  INTEGER, PARAMETER :: REPLACE  = 3       !< read/write access replacing file
  INTEGER, PARAMETER :: APPEND   = 4       !< read/write access at end
  !> \name
  !! #### file formats
  CHARACTER(LEN=9), PARAMETER :: ASCII = "formatted"   !< for ASCII data
  CHARACTER(LEN=11), PARAMETER :: BIN  = "unformatted" !< for BINARY data
  !> \}

  ! file formats
  INTEGER, PARAMETER :: BINARY  = 1
  INTEGER, PARAMETER :: GNUPLOT = 2
  INTEGER, PARAMETER :: VTK     = 4
  INTEGER, PARAMETER :: NPY     = 5
  INTEGER, PARAMETER :: HDF     = 6
  INTEGER, PARAMETER :: XDMF    = 7
  !--------------------------------------------------------------------------!
  INTEGER, SAVE :: lastunit = 10

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fileio_base, Output_TYP, TSOutput_TYP, ValPtr_TYP, &
       ! constants
       BINARY, GNUPLOT, VTK, NPY, HDF, XDMF, &
       READONLY, READEND, REPLACE, APPEND, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL,&
#endif
!       FILE_EXISTS, &
       ASCII, BIN


CONTAINS


  !> \public Generic constructor for file I/O
  !!
  SUBROUTINE InitFileIO(this,Mesh,Physics,Timedisc,Sources,config,IO, &
                        fmtname,fext)
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
    INTEGER                        :: fcycles    !< fcycles number of file cycles
    INTEGER                        :: unit       !<  unit fortran i/o unit number
    LOGICAL                        :: success
    CHARACTER(LEN=32)              :: timestamp
    INTEGER                        :: i,fstatus
    INTEGER                        :: count_def, fcycles_def, dtwall_def, ncfmt_def
    INTEGER                        :: cartcoords
    REAL                           :: stoptime_def
    REAL                           :: time,new_time
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
    CALL GetAttr(config, "unit"      , unit , lastunit+1)
    lastunit = unit

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
    ! check file name length
    IF (LEN_TRIM(fname).GT.FNAMLEN) &
       CALL this%Error("InitFileIO","file name too long")
    ! check file path length
    IF (LEN_TRIM(fpath).GT.FPATLEN) &
       CALL this%Error("InitFileIO","file path too long")
    ! check file name extension
    IF (LEN_TRIM(fext).GT.FEXTLEN-1) &
       CALL this%Error("InitFileIO","file name extension too long")

    ! format string for writing file names with explicit time step
    WRITE (cycfmt, "('(A,I',I1,'.',I1,',A)')") FCYCLEN-1,FCYCLEN-1

#ifdef PARALLEL
    ! turn on multiple file output if requested
    IF (this%multfiles) THEN
       ! check number of parallel processes
       IF (this%GetNumProcs().GT.MAXMLTFILES) &
          CALL this%Error("InitFileIO","number of processes for multiple file output exceeds limits")
       fmextstr = this%MakeMultstr()
    END IF
#endif

    this%path      = TRIM(fpath)
    this%filename  = TRIM(fname)
    this%extension = TRIM(fext)
    this%cycles    = fcycles_def
    this%unit      = unit
    this%stoptime  = stoptime_def
    this%dtwall    = dtwall_def
    this%time      = 0.
    this%count     = count_def
    this%step      = 0
    this%offset    = 0
    this%error_io  = 0

    ! compute the (actual) output time
    time = ABS(this%stoptime) / this%count
    this%time = time*FLOOR(Timedisc%time/time)

    ! print some information
    CALL this%Info(" FILEIO---> file name:         " // TRIM(this%GetFilename()))
    IF (success) THEN
       WRITE (timestamp,'(ES10.4)') Timedisc%time
       CALL this%Info("            time stamp:        " // TRIM(timestamp))
    END IF

    ! time for next output
    IF (Timedisc%time.GT.0.0) CALL this%IncTime()
    IF (ASSOCIATED(oldconfig)) CALL DeleteDict(oldconfig)
  END SUBROUTINE InitFileIO

  !> \public Increments the counter for timesteps and sets the time for next output
  !!
  PURE SUBROUTINE IncTime(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(INOUT) :: this!< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    this%time = this%time + ABS(this%stoptime) / this%count
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


  !> \public Get a file label (multiples files in parallel mode)
  !! without filenumber => use GetRank; with filenumber < 0 => empty label
  !! \result file label
  FUNCTION MakeMultstr(this,fn) RESULT (multstr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(IN) :: this !< \param [in,out] this fileio type
    INTEGER, OPTIONAL,  INTENT(IN) :: fn       !< \param [in] fn number of file
    CHARACTER(LEN=FMLTLEN)         :: multstr
    !------------------------------------------------------------------------!
    INTEGER                        :: fn_l
    CHARACTER(LEN=32)              :: mextfmt
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (PRESENT(fn)) THEN
      fn_l = fn
    ELSE
      fn_l = this%GetRank()
    END IF
    ! with fn < 0 you can suppress a label
    IF (this%multfiles .AND. fn_l .GE. 0) THEN
        WRITE (mextfmt, "('(A2,I',I1,'.',I1,')')") FMLTLEN-2,FMLTLEN-2

        WRITE (multstr, mextfmt) "-r", fn_l
    ELSE
      multstr = ""
    END IF
#else
    multstr = ""
#endif
  END FUNCTION MakeMultstr


  !> \public Get the current file name without path
  !! e.g. important for vtk (pvts files)
  !! \result current file name
  FUNCTION GetBasename(this,fn) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER, OPTIONAL, INTENT(IN)  :: fn   !< \param [in] fn number of file
    CHARACTER(LEN=256)             :: fname
    !------------------------------------------------------------------------!
    CHARACTER(LEN=FCYCLEN+2)       :: cycstr
    !------------------------------------------------------------------------!
    IF (this%cycles.GT.0) THEN
       ! generate a file name with time step
       WRITE (cycstr, FMT=TRIM(cycfmt)) "_", MODULO(this%step,this%cycles), "."
       WRITE (fname,"(A,A,A,A)")  TRIM(this%filename),&
              TRIM(MakeMultstr(this,fn)),TRIM(cycstr),TRIM(this%extension)
    ELSE
       ! file name + extension
       WRITE (fname,"(A,A,A,A)") TRIM(this%filename),&
              TRIM(MakeMultstr(this,fn)),".",TRIM(this%extension)
    END IF
  END FUNCTION GetBasename

  !> \public Get the current file name
  !! \result current file name
  FUNCTION GetFilename(this,fn) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(IN) :: this !< \param [in] this fileio type
    INTEGER, OPTIONAL, INTENT(IN)  :: fn   !< \param [in] fn number of file
    CHARACTER(LEN=256)             :: fname
    !------------------------------------------------------------------------!
    fname = TRIM(this%path) // TRIM(GetBasename(this,fn))
  END FUNCTION GetFilename

!  !> \public Generic routine to open a file
!  !!
!  SUBROUTINE OpenFile(this,action)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this      !< \param [in,out] this fileio type
!    INTEGER          :: action    !< \param [in] action
!    !------------------------------------------------------------------------!
!    INTENT(IN)       :: action
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!    IF (.NOT.Initialized(this)) &
!         CALL Error(this,"OpenFile","file module uninitialized")
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL OpenFile_binary(this,action)
!    CASE(GNUPLOT)
!       CALL OpenFile_gnuplot(this,action)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL OpenFile_netcdf(this,action)
!#endif
!    CASE(VTK)
!       CALL OpenFile_vtk(this,action)
!    CASE(NPY)
!       CALL OpenFile_npy(this,action)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL OpenFile_hdf5(this,action)
!#endif
!    CASE(XDMF)
!       CALL OpenFile_xdmf(this,action)
!    END SELECT
!  END SUBROUTINE OpenFile
!
!  !> \public Generic routine to close the a file
!  !!
!  SUBROUTINE CloseFile(this)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
!    !------------------------------------------------------------------------!
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!    IF (.NOT.Initialized(this)) &
!         CALL Error(this,"CloseFile","file module uninitialized")
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL CloseFile_binary(this)
!    CASE(GNUPLOT)
!       CALL CloseFile_gnuplot(this)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL CloseFile_netcdf(this)
!#endif
!    CASE(VTK)
!       CALL CloseFile_vtk(this)
!    CASE(NPY)
!       CALL CloseFile_npy(this)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL CloseFile_hdf5(this)
!#endif
!    CASE(XDMF)
!       CALL CloseFile_xdmf(this)
!    END SELECT
!  END SUBROUTINE CloseFile

!  !> \public Generic routine to write a header to a file
!  !!
!  SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(fileio_base),INTENT(INOUT)  :: this      !< \param [in,out] this fileio type
!    CLASS(mesh_base),INTENT(IN)  :: Mesh      !< \param [in] Mesh mesh type
!    CLASS(physics_base),INTENT(IN) :: Physics   !< \param [in] Physics physics type
!    !> \param [in] Header dict with header
!    TYPE(Dict_TYP),POINTER :: Header,&
!                              IO !< \param [in] IO dict with I/O arrays
!    !------------------------------------------------------------------------!
!    ! create new empty data file
!    CALL OpenFile(this,REPLACE)
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL WriteHeader_binary(this)
!    CASE(GNUPLOT)
!       CALL WriteHeader_gnuplot(this)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL WriteHeader_netcdf(this,Mesh,Physics)
!#endif
!    CASE(VTK)
!       CALL WriteHeader_vtk(this,Mesh,Physics,Header,IO)
!    CASE(NPY)
!       CALL WriteHeader_npy(this)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL WriteHeader_hdf5(this,Mesh,Header)
!#endif
!    CASE(XDMF)
!       CALL WriteHeader_xdmf(this)
!    END SELECT
!    CALL CloseFile(this)
!  END SUBROUTINE WriteHeader
!
!  !> \public Generic routine to read a header from a file
!  !!
!  SUBROUTINE ReadHeader(this,Mesh,Physics,Header,success)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
!    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
!    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
!    !> \param [in] Header dict with header
!    TYPE(Dict_TYP),POINTER :: Header
!    LOGICAL           :: success   !< \param [out] success
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Physics
!    INTENT(OUT)       :: success
!    INTENT(INOUT)     :: this
!    !------------------------------------------------------------------------!
!    CALL OpenFile(this,READONLY)
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL ReadHeader_binary(this,success)
!    CASE(GNUPLOT)
!       CALL ReadHeader_gnuplot(this,success)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL ReadHeader_netcdf(this,Mesh,Physics,success)
!#endif
!    CASE(VTK)
!       CALL ReadHeader_vtk(this,success)
!    CASE(NPY)
!       CALL ReadHeader_npy(this,success)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL ReadHeader_hdf5(this,Header,success)
!#endif
!    CASE(XDMF)
!       CALL ReadHeader_xdmf(this,success)
!    END SELECT
!    CALL CloseFile(this)
!  END SUBROUTINE ReadHeader
!
!  !> \public Generic routine to write a timestamp to a file
!  !!
!  SUBROUTINE WriteTimestamp(this,time)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
!    REAL             :: time    !< \param [in] time
!    !------------------------------------------------------------------------!
!    INTENT(IN)       :: time
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!    CALL OpenFile(this,APPEND)
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL WriteTimestamp_binary(this,time)
!    CASE(GNUPLOT)
!       CALL WriteTimestamp_gnuplot(this,time)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL WriteTimestamp_netcdf(this,time)
!#endif
!    CASE(VTK)
!       CALL WriteTimestamp_vtk(this,time)
!    CASE(NPY)
!       CALL WriteTimestamp_npy(this,time)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL WriteTimestamp_hdf5(this,time)
!#endif
!    CASE(XDMF)
!       CALL WriteTimestamp_xdmf(this,time)
!    END SELECT
!    CALL CloseFile(this)
!  END SUBROUTINE WriteTimestamp
!
!  !> \public Generic routine to read a timestamp from a file
!  !!
!  SUBROUTINE ReadTimestamp(this,time)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this    !< \param [in,out] this fileio type
!    REAL             :: time    !< \param [out] time
!    !------------------------------------------------------------------------!
!    INTENT(OUT)      :: time
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!    ! open datafile in r/w mode and seek to the end
!    CALL OpenFile(this,READEND)
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL ReadTimestamp_binary(this,time)
!    CASE(GNUPLOT)
!       CALL ReadTimestamp_gnuplot(this,time)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL ReadTimestamp_netcdf(this,time)
!#endif
!    CASE(VTK)
!       CALL ReadTimestamp_vtk(this,time)
!    CASE(NPY)
!       CALL ReadTimestamp_npy(this,time)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL ReadTimestamp_hdf5(this,time)
!#endif
!    CASE(XDMF)
!       CALL ReadTimestamp_xdmf(this,time)
!    END SELECT
!    CALL CloseFile(this)
!  END SUBROUTINE ReadTimestamp
!
!  !> \public Generic routine to write a dataset to a file
!  !!
!  SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
!    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
!    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
!    TYPE(Fluxes_TYP)  :: Fluxes    !< \param [in] Fluxes fluxes type
!    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] Timedisc timedisc type
!    !> \param [in] Header dict with header
!    TYPE(Dict_TYP),POINTER :: Header,&
!                              IO !< \param [in] IO dict with I/O arrays
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Physics,Fluxes
!    INTENT(INOUT)     :: this,Timedisc
!    !------------------------------------------------------------------------!
!    INTEGER           :: k, ierror
!    !------------------------------------------------------------------------!
!    !TODO: Wo gehört das hin?
!    ! calculate boundary fluxes, if they were requested for write
!    IF(ASSOCIATED(Timedisc%bflux)) THEN
!        DO k=1,4
!            Timedisc%bflux(:,k) = GetBoundaryFlux(Fluxes,Mesh,Physics,k)
!        END DO
!#ifdef PARALLEL
!        CALL MPI_BCAST(Timedisc%bflux, Physics%VNUM*4, DEFAULT_MPI_REAL, 0, &
!                Mesh%comm_cart, ierror)
!#endif
!    END IF
!    ! write the header if either this is the first data set we write or
!    ! each data set is written into a new file
!    IF ((this%step.EQ.0).OR.(this%cycles.GT.0)) THEN
!       CALL WriteHeader(this,Mesh,Physics,Header,IO)
!    END IF
!    CALL OpenFile(this,APPEND)
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL WriteDataset_binary(this,Mesh,Physics,Fluxes,Timedisc,IO)
!    CASE(GNUPLOT)
!       CALL WriteDataset_gnuplot(this,Mesh)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL WriteDataset_netcdf(this,Mesh,Physics,Timedisc)
!#endif
!    CASE(VTK)
!       CALL WriteDataset_vtk(this,Mesh,Physics,Fluxes,Timedisc,IO)
!    CASE(NPY)
!       CALL WriteDataset_npy(this,Mesh,Physics,Timedisc)
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL WriteDataset_hdf5(this,Mesh,Physics,Timedisc,Fluxes,IO)
!#endif
!    CASE(XDMF)
!       CALL WriteDataset_xdmf(this,Mesh,Physics,Fluxes,Timedisc,IO)
!    END SELECT
!    CALL CloseFile(this)
!    ! append the time stamp
!    CALL WriteTimestamp(this,Timedisc%time)
!    CALL IncTime(this)
!  END SUBROUTINE WriteDataset

!  !> \public Generic routine to read a dataset from a file
!  !!
!  SUBROUTINE ReadDataset(this,Mesh,Physics,Timedisc,IO)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
!    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] Mesh mesh type
!    TYPE(Physics_TYP) :: Physics   !< \param [in] Physics physics type
!    TYPE(Timedisc_TYP):: Timedisc  !< \param [in] Timedisc timedisc type
!    TYPE(Dict_TYP),POINTER :: IO   !< \param [in] IO dict with I/O arrays
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh
!    INTENT(INOUT)     :: this,Timedisc,Physics
!    !------------------------------------------------------------------------!
!    CALL OpenFile(this,APPEND)
!    SELECT CASE(GetType(this))
!    CASE(BINARY)
!       CALL ReadDataset_binary(this,Mesh,Physics,Timedisc)
!    CASE(GNUPLOT)
!       CALL ReadDataset_gnuplot(this,Mesh,Physics,Timedisc)
!#ifdef HAVE_NETCDF
!    CASE(NETCDF)
!       CALL ReadDataset_netcdf(this,Mesh,Physics,Timedisc)
!#endif
!    CASE(VTK)
!       CALL Error(this,"ReadDataset","function not supported")
!    CASE(NPY)
!       CALL Error(this,"ReadDataset","function not supported")
!#ifdef HAVE_HDF5_MOD
!    CASE(HDF)
!       CALL ReadDataset_hdf5(this,Mesh,Physics,Timedisc,IO)
!#endif
!    CASE(XDMF)
!       CALL ReadDataset_xdmf(this,Mesh,Physics,Timedisc)
!    END SELECT
!    CALL CloseFile(this)
!    ! calculate conservative variables
!    CALL Convert2Conservative(Physics,Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN,Mesh%JMAX, &
!         Timedisc%pvar,Timedisc%cvar)
!  END SUBROUTINE ReadDataset

  !> \public Generic deconstructor of the file I/O
  !!
  SUBROUTINE FinalizeFileio(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base),INTENT(INOUT) :: this !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("CloseFileIO","not initialized")
  END SUBROUTINE FinalizeFileIO

END MODULE fileio_base_mod
