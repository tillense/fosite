!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_gnuplot.f90                                                #
!#                                                                           #
!# Copyright (C) 2008-2017                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup fileio
!! \key{decimals,INTEGER,get number of decimal places; set to default if not given,5}
!! \key{cartcoords,INTEGER,check if cartesian coordinates are selected for output,0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Jannes Klee
!!
!! \brief I/O for GNUPLOT readable tabular files
!!
!! This module implements a file I/O, which files can be read by GNUPLOT.
!! It writes the configuration (dictionary) as header.
!! It is possible to select which data arrays should be written.
!!
!! \extends fileio_common
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_gnuplot_mod
  USE fileio_base_mod
  USE geometry_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE timedisc_base_mod
  USE fluxes_base_mod
  USE sources_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
!  INTERFACE Error
!     MODULE PROCEDURE Error_gnuplot, Error_common
!  END INTERFACE
  ! Private Attributes section starts here:
  !> \name some string lengths
  INTEGER, PARAMETER      :: HLEN = 10000         !< header length
  INTEGER, PARAMETER      :: DEFAULT_DECS = 5     !< default decimal places
  !> \name some special strings
  CHARACTER, PARAMETER    :: SP = ACHAR(32)       !< space
  CHARACTER, PARAMETER    :: LF = ACHAR(10)       !< line feed
  CHARACTER*2, PARAMETER  :: RECSEP = SP // SP    !< data record separator
  CHARACTER*2, PARAMETER  :: LINSEP = SP // LF    !< line separator
  CHARACTER*2, PARAMETER  :: BLKSEP = LF // LF    !< block separator
  !> the header string
  CHARACTER(LEN=30), PARAMETER :: &
         header_string = "# Data output of fosite" // LINSEP
  CHARACTER(LEN=HLEN)  :: header_buf              !< buffer of header
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(fileio_base) :: fileio_gnuplot

  CONTAINS
    ! methods
    PROCEDURE :: InitFileIO_gnuplot
    PROCEDURE :: Finalize
    PROCEDURE :: WriteHeader
!    PROCEDURE :: ReadHeader
!    PROCEDURE :: WriteTimestamp
!    PROCEDURE :: ReadTimestamp
    PROCEDURE :: WriteDataset
    PROCEDURE :: ReadDataset
    PROCEDURE :: OpenFile
    PROCEDURE :: CloseFile
    PROCEDURE :: AdjustTimestep
    PROCEDURE :: IncTime
    PROCEDURE :: GetFilename
    PROCEDURE :: GetFilestatus
  END TYPE

  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       FileIO_TYP, &
       ! constants
       READONLY, READEND, REPLACE, APPEND, &
       ASCII, BIN, &
       FILE_EXISTS
 !--------------------------------------------------------------------------!

CONTAINS

  !> \public Generic constructor for file I/O
  !!
  !! Initilizes the file I/O type, filename and extension, stoptime,
  !! number of outputs, number of files,
  !! mode for parallel output (separate files), unit number
  SUBROUTINE InitFileIO(this,Mesh,Physics,fmt,fmtname,fpath,filename,extension, &
       stoptime,dtwall,count,fcycles,sepfiles,unit)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base)  :: this            !< \param [in,out] this fileio type
    CLASS(mesh_base)    :: Mesh            !< \param [in] Mesh mesh type
    CLASS(physics_base) :: Physics         !< \param [in] Physics physics type
    INTEGER             :: fmt             !< \param [in] fmt fileio type number
    CHARACTER(LEN=*)    :: fmtname         !< \param [in] fmtname name of fileio
    CHARACTER(LEN=*)    :: fpath           !< \param [in] fpath
    CHARACTER(LEN=*)    :: filename        !< \param [in] filename
    CHARACTER(LEN=*)    :: extension       !< \param [in] extension file extension
    REAL                :: stoptime        !< \param [in] stoptime
    INTEGER             :: dtwall          !< \param [in] dtwall wall clock time
    INTEGER             :: count           !< \param [in] count number of outputs
    INTEGER             :: fcycles         !< \param [in] fcycles file cycle number
    LOGICAL             :: sepfiles        !< \param [in] sepfiles different files
    INTEGER, OPTIONAL   :: unit            !< \param [in] unit fileio unit number
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,fmtname,fpath,filename,extension,stoptime, &
         dtwall,count,fcycles,sepfiles,unit
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! basic FileIO initialization
    CALL this%InitFileIO_common(fmt,fmtname,fpath,filename,extension,fcycles,sepfiles,unit)

    this%stoptime = stoptime
    this%dtwall   = dtwall
    this%time     = 0.
    this%count    = count
    this%step     = 0
    ! count the number of output columns, i.e. fields per data point

#ifdef PARALLEL
    ! check data type extents in files
    ! first try to create a new dummy file
    CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(IOR(MPI_MODE_RDWR,&
         MPI_MODE_CREATE),IOR(MPI_MODE_EXCL,MPI_MODE_DELETE_ON_CLOSE)),&
         MPI_INFO_NULL,this%handle,this%error)
    ! maybe file exists
    IF (this%error.NE.0) CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),&
         MPI_MODE_RDONLY,MPI_INFO_NULL,this%handle,this%error)
    ! then check the data type sizes
    ! extent of integer in file
    IF (this%error.EQ.0) CALL MPI_File_get_type_extent(this%handle,MPI_INTEGER,&
         this%intext,this%error)
    ! extent of real in file
    IF (this%error.EQ.0) CALL MPI_File_get_type_extent(this%handle,DEFAULT_MPI_REAL,&
         this%realext,this%error)
    IF (this%error.NE.0) CALL Error(this,"InitFileIO","unable to check file properties")
    CALL MPI_File_close(this%handle,this%error)
#endif
  END SUBROUTINE InitFileIO


  !> \public Constructor for the GNUPLOT file I/O
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs,
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_gnuplot(this,Mesh,Physics,IO,fmt,fpath,filename,stoptime,dtwall,&
       count,fcycles,unit,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base)              :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base)                :: Mesh     !< \param [in] Mesh mesh type
    CLASS(physics_base)             :: Physics  !< \param [in] Physics Physics type
    TYPE(Dict_TYP),POINTER          :: IO       !< \param [in] IO Dictionary for I/O
    !> \param [in] config Dictionary with configuration
    TYPE(Dict_TYP),POINTER,OPTIONAL :: config
    INTEGER                         :: fmt      !< \param [in] fmt fileio type number
    CHARACTER(LEN=*)                :: fpath    !< \param [in] fpath
    CHARACTER(LEN=*)                :: filename !< \param [in] filename
    REAL                            :: stoptime !< \param [in] stoptime
    INTEGER                         :: dtwall   !< \param [in] dtwall wall clock time
    INTEGER                         :: count    !< \param [in] count number of outputs
    INTEGER                         :: fcycles  !< \param [in] fcycles file cycle number
    INTEGER                         :: unit     !< \param [in] unit fileio unit number
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER                   :: node
    REAL, DIMENSION(:,:), POINTER             :: dummy2
    REAL, DIMENSION(:,:,:), POINTER           :: dummy3
    INTEGER                                   :: cartcoords
    INTEGER                                   :: depth
    INTEGER                                   :: err
#ifdef PARALLEL
    INTEGER                                   :: i
    INTEGER, DIMENSION(Mesh%IMAX-Mesh%IMIN+1) :: blocklen,indices
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,fmt,fpath,filename,stoptime,dtwall,count,fcycles,&
                     unit
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%InitFileIO(Mesh,Physics,fmt,"GNUPLOT",fpath,filename,"dat",stoptime,&
         dtwall,count,fcycles,.FALSE.,unit)

    CALL GetAttr(config, "/datafile/decimals", this%DECS, DEFAULT_DECS)
    ! compute length of character field for real number output
    ! and check if linebuffer is large enough
    ! flen = 1 (sign) + 1 (one digit) + 1 (decimal point) + decs (decimal places)
    !      + 1 (E character) + 1 (sign of exponent) + 2 (exponent digits) + 2 (spaces)
    this%FLEN = this%DECS + 9
    this%maxcols = len(this%linebuf)/this%FLEN-1

    ALLOCATE(this%output(this%maxcols),STAT=err)
    IF (this%error.NE.0) &
       CALL this%Error("InitFileIO_gnuplot","memory allocation failed for this%output")

    ! check if cartesian coordinates are selected for output;
    ! default: curvilinear coordinates (0)
    CALL GetAttr(config, "/datafile/cartcoords", cartcoords, 0)
    depth = 1
    node => config
    CALL WriteHeaderString(header_buf,node,depth)
    IF (cartcoords.EQ.0) THEN
       CALL GetAttr(IO,"/mesh/bary_curv",dummy3)
    ELSE
       CALL GetAttr(IO,"/mesh/bary_centers",dummy3)
    END IF
    ! pointer to sub-array: set lower bounds
    ! use special function for bound remapping
    IF (Mesh%INUM.EQ.1) THEN
      ! use format string as temp
      WRITE (this%fmtstr,'(A5,I2,A1)')'(A1,A',this%FLEN-3,')'
      WRITE(this%linebuf,TRIM(this%fmtstr))'#','y'
      this%output(1)%val => dummy3(:,:,2)
      this%COLS = 1
    ELSE IF (Mesh%JNUM.EQ.1) THEN
      WRITE (this%fmtstr,'(A5,I2,A1)')'(A1,A',this%FLEN-3,')'
      WRITE(this%linebuf,TRIM(this%fmtstr))'#','x'
      this%output(1)%val => dummy3(:,:,1)
      this%COLS = 1
    ELSE
      WRITE (this%fmtstr,'(A5,I2,A2,I2,A1)')'(A1,A',this%FLEN-3,',A',this%FLEN-1,')'
      WRITE(this%linebuf,TRIM(this%fmtstr))'#','x','y'
      this%output(1)%val => dummy3(:,:,1)
      this%output(2)%val => dummy3(:,:,2)
      this%COLS = 2
    END IF

    ! set output-pointer and count the number of output columns
    WRITE (this%fmtstr,'(A,I2,A1)')'(A',this%FLEN-1,')'
    node => IO
    CALL this%GetOutputPointer(Mesh,node,this%COLS)

    ! length of one output line
    this%linelen = this%COLS * this%FLEN
    IF (this%linelen.GT.LEN(this%linebuf)) &
       CALL this%Error("InitFileIO_gnuplot", &
          "linebuffer to small; reducing decimals or number of output fields may help")

    header_buf = TRIM(header_string) // TRIM(header_buf) // this%linebuf(1:this%linelen)

    ! local domain size
    this%inum = Mesh%IMAX - Mesh%IMIN + 1
    this%jnum = Mesh%JMAX - Mesh%JMIN + 1

!#ifdef PARALLEL
!    ! create new data type handle for one line
!    CALL MPI_Type_contiguous(this%linelen,MPI_CHARACTER,this%basictype,this%error)
!    CALL MPI_Type_commit(this%basictype,this%error)
!
!    ! number of output blocks
!    this%blocknum = Mesh%IMAX - Mesh%IMIN + 1
!    ! size of the output buffer
!    this%bufsize  = Mesh%JMAX - Mesh%JMIN + 1
!
!    ! allocate memory for output buffer and displacement records
!    ALLOCATE(this%outbuf(this%linelen,Mesh%JMIN:Mesh%JMAX), &
!         STAT=err)
!    IF (this%error.NE.0) THEN
!       CALL Error(this,"InitFileIO_gnuplot","memory allocation failed for this%outbuf")
!    END IF
!
!    blocklen(:) = this%bufsize
!    DO i=Mesh%IMIN,Mesh%IMAX
!        indices(i-Mesh%IMIN+1) = (i-1)*Mesh%JNUM + Mesh%JMIN - 1
!    END DO
!
!    ! new file type for the staggered data
!    CALL MPI_Type_indexed(this%blocknum,blocklen,indices, &
!         this%basictype,this%filetype,this%error)
!    CALL MPI_Type_commit(this%filetype,this%error)
!#endif
    ! write the format string for one entry in the data file:
    ! FLEN-2 characters for the number and 2 for the separators
    WRITE (this%fmtstr,'(A3,I2,A,I2.2,A5)') '(ES', this%FLEN-2, '.', this%DECS,',A,A)'
    ! write format string for one output line
    WRITE (this%linefmt, '(A,I0,A)') "(A", this%linelen-1, ")"
  END SUBROUTINE InitFileIO_gnuplot


  !> Creates a string with the configuration (from the dictionary)
  !!
  RECURSIVE SUBROUTINE WriteHeaderString(string,root,k,prefix)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER    :: root,node,subnode
    CHARACTER(LEN=*)           :: string
    CHARACTER(LEN=*), OPTIONAL :: prefix
    CHARACTER(LEN=128)         :: buf
    !------------------------------------------------------------------------!
    INTEGER                    :: idummy, k
    LOGICAL                    :: ldummy
    CHARACTER(LEN=128)         :: cdummy
    REAL                       :: rdummy
    !------------------------------------------------------------------------!
    INTENT(INOUT)              :: string,k
    !------------------------------------------------------------------------!

    node => root
    DO WHILE(ASSOCIATED(node))
       SELECT CASE(GetDatatype(node))
       CASE(DICT_INT)
          CALL GetAttr(node,GetKey(node),idummy)
          WRITE(buf,'(A1,A25,I14,A)')'#',TRIM(GetKey(node))//": ",idummy, LINSEP
          WRITE(string(k:),'(A)')buf
          k = k + LEN(TRIM(buf))
       CASE(DICT_REAL)
          CALL GetAttr(node,GetKey(node),rdummy)
          WRITE(buf,'(A1,A25,ES14.5,A)')'#',TRIM(GetKey(node))//": ",rdummy, LINSEP
          WRITE(string(k:),'(A)')buf
          k = k + LEN(TRIM(buf))
       CASE(DICT_CHAR)
!          CALL GetAttr(node,GetKey(node),cdummy)
!          WRITE(buf,'(A1,A25,A,A)')'#',TRIM(GetKey(node))//": ",TRIM(cdummy), LINSEP
!          WRITE(string(k:),'(A)')buf
!          k = k + LEN(TRIM(buf))
       CASE(DICT_BOOL)
          CALL GetAttr(node,GetKey(node),ldummy)
          WRITE(buf,'(A1,A25,L14,A)')'#',TRIM(GetKey(node))//": ",ldummy, LINSEP
          WRITE(string(k:),'(A)')buf
          k = k + LEN(TRIM(buf))
       END SELECT
       IF(HasChild(node)) THEN
          IF (present(prefix)) THEN
             WRITE(buf,'(A)')'#  ['//TRIM(prefix)//'/'//TRIM(GetKey(node))//']' // LINSEP
          ELSE
             WRITE(buf,'(A)')'#  ['//TRIM(GetKey(node))//']' // LINSEP
          END IF
          WRITE(string(k:),'(A)')buf
          k = k + LEN(TRIM(buf))
          IF (present(prefix)) THEN
             buf = TRIM(prefix)//'/'//TRIM(GetKey(node))
          ELSE
             buf = TRIM(GetKey(node))
          END IF
          subnode => GetChild(node)
          CALL WriteHeaderString(string,subnode,k,TRIM(buf))
       END IF
       node => GetNext(node)
    END DO
  END SUBROUTINE WriteHeaderString

  !> Creates a list of all data arrays which will be written to file
  !!
  !! Therefore it ignores all arrays with coordinates and checks if the data
  !! arrays are of the dimension of the mesh.
  RECURSIVE SUBROUTINE GetOutputPointer(this,Mesh,node,k)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base)      :: this  !< \param [in,out] this fileio type
    CLASS(mesh_base)        :: Mesh  !< \param [in] mesh mesh type
    TYPE(Dict_TYP), POINTER :: node  !< \param [in,out] node pointer to (sub-)dict
    INTEGER                 :: k     !< \param [in,out] k number of data arrays
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER           :: dir
    REAL, DIMENSION(:,:), POINTER     :: dummy2
    REAL, DIMENSION(:,:,:), POINTER   :: dummy3
    REAL, DIMENSION(:,:,:,:), POINTER :: dummy4
    INTEGER, DIMENSION(4)             :: dims
    INTEGER                           :: dim3,dim4,i,j
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh
    INTENT(INOUT) :: this,k
    !------------------------------------------------------------------------!
    ! reset error code
    this%error = 0
    DO WHILE(ASSOCIATED(node))
      ! check for directory and exclude any coordinates (these are handled elsewhere)
      IF(HasChild(node)) THEN
      ! recursion
        CALL GetAttr(node,GetKey(node),dir)
        CALL GetOutputPointer(this,Mesh,dir,k)
      ELSE IF (HasData(node).AND. .NOT.(GetKey(node).EQ."bary_curv".OR. &
        GetKey(node).EQ."bary_centers".OR.GetKey(node).EQ."corners")) THEN
        WRITE(this%linebuf(k*this%FLEN:),TRIM(this%fmtstr))TRIM(GetKey(node))
      ! value found!
        SELECT CASE(GetDataType(node))
        CASE(DICT_REAL_TWOD)
          CALL GetAttr(node,GetKey(node),dummy2)
          dims(1:2) = SHAPE(dummy2)
          IF((dims(1).EQ.Mesh%IMAX-Mesh%IMIN+1).AND.&
             (dims(2).EQ.Mesh%JMAX-Mesh%JMIN+1)) THEN
            k = k+1
            IF (k .GT. this%maxcols) THEN
               this%error = 1
               EXIT
            END IF
            this%output(k)%val=> dummy2
          END IF
        CASE(DICT_REAL_THREED)
          CALL GetAttr(node,GetKey(node),dummy3)
          dims(1:3) = SHAPE(dummy3)
          IF((dims(1).EQ.Mesh%IMAX-Mesh%IMIN+1).AND.&
             (dims(2).EQ.Mesh%JMAX-Mesh%JMIN+1)) THEN
            dim3 = SIZE(dummy3, DIM = 3)
            IF (k+dim3 .GT. this%maxcols) THEN
               this%error = 1
               EXIT
            END IF
            DO i=k+1, k+dim3
              this%output(i)%val => dummy3(:,:,i-k)
            END DO
            k = k+dim3
          END IF
        CASE(DICT_REAL_FOURD)
          CALL GetAttr(node,GetKey(node),dummy4)
          dims(1:4) = SHAPE(dummy4)
          IF((dims(1).EQ.Mesh%IMAX-Mesh%IMIN+1).AND.&
             (dims(2).EQ.Mesh%JMAX-Mesh%JMIN+1)) THEN
            dim3 = SIZE(dummy4, DIM = 3)
            dim4 = SIZE(dummy4, DIM = 4)
            IF (k+dim3*dim4 .GT. this%maxcols) THEN
               this%error = 1
               EXIT
            END IF
            DO j=1,dim4
              DO i=1,dim3
              this%output(k+j+(i-1)*dim4)%val => dummy4(:,:,i,j)
              END DO
            END DO
            k = k+dim3*dim4
          END IF
        CASE DEFAULT
          !do nothing (wrong type)
        END SELECT
      END IF
      node=>GetNext(node)
    END DO
    IF (this%error.NE.0) &
         CALL Error(this,"GetOutputPointer_gnuplot","number of output fields exceeds upper limit")
  END SUBROUTINE GetOutputPointer


  !> \public Adjust the current timestep
  !!
  !! Last timestep before output must fit to desired time for output.
 PURE SUBROUTINE AdjustTimestep(this,time,dt,dtcause)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this    !< \param [in] this fileio type
    REAL               :: time    !< \param [in,out] time
    REAL               :: dt      !< \param [in,out] dt timestep
    INTEGER            :: dtcause !< \param [in,out] dtcause cause of smallest dt
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this
    INTENT(INOUT)      :: time,dt,dtcause
    !------------------------------------------------------------------------!
    IF ((time+dt)/this%time.GT.1.0) THEN
       dt = this%time - time
       dtcause = DTCAUSE_FILEIO
    ELSE IF((time+1.5*dt)/this%time.GT.1.0) THEN
       dt = 0.5*(this%time - time)
       dtcause = DTCAUSE_FILEIO
    END IF
  END SUBROUTINE AdjustTimestep

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

  !> \public Generic routine to open a file
  !!
  SUBROUTINE OpenFile(this,action,fformat)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this               !< \param [in,out] this fileio type
    INTEGER            :: action             !< \param [in] action mode of open
    CHARACTER(LEN=*)   :: fformat            !< \param [in] fformat file format
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset!< \param [in] offset offset for MPI
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action,fformat
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    SELECT CASE(action)
    CASE(READONLY)
#ifdef PARALLEL
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),MPI_MODE_RDONLY, &
            MPI_INFO_NULL,this%handle,this%error)
        this%offset = 0
        CALL MPI_File_seek(this%handle,this%offset,MPI_SEEK_SET,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="OLD", &
            ACTION="READ",POSITION="REWIND",IOSTAT=this%error)
       !REWIND (UNIT=this%unit,IOSTAT=this%error)
#endif
    CASE(READEND)
#ifdef PARALLEL
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_RDONLY,&
            MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)
       ! opening in append mode doesn't seem to work for pvfs2, hence ...
       offset = 0
       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
       CALL MPI_File_sync(this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="OLD", &
            ACTION="READ",POSITION="APPEND",IOSTAT=this%error)
#endif
    CASE(REPLACE)
#ifdef PARALLEL
       CALL MPI_File_delete(GetFilename(this),MPI_INFO_NULL,this%error)
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_WRONLY,&
            MPI_MODE_CREATE),MPI_INFO_NULL,this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="REPLACE",&
            ACTION="WRITE",POSITION="REWIND",IOSTAT=this%error)
#endif
    CASE(APPEND)
#ifdef PARALLEL
       CALL MPI_File_open(MPI_COMM_WORLD,GetFilename(this),IOR(MPI_MODE_RDWR,&
            MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)
       ! opening in append mode doesn't seem to work for pvfs2, hence ...
       offset = 0
       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
       CALL MPI_File_sync(this%handle,this%error)
#else
       OPEN(this%unit,FILE=GetFilename(this),FORM=fformat,STATUS="OLD",&
            ACTION="READWRITE",POSITION="APPEND",IOSTAT=this%error)
#endif
    CASE DEFAULT
       CALL Error(this,"OpenFile","Unknown access mode.")
    END SELECT
  END SUBROUTINE OpenFile


  !> \public Specific routine to open a file for gnuplot I/O
  !!
  SUBROUTINE OpenFile_gnuplot(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this   !< \param [in,out] this fileio type
    INTEGER            :: action !< \param [in] action mode of file access
    !------------------------------------------------------------------------!
    INTENT(IN)         :: action
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL this%OpenFile(action,ASCII)
  END SUBROUTINE OpenFile_gnuplot


  !> \public routine to close a file
  !!
  SUBROUTINE CloseFile_gnuplot(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_File_close(this%handle,this%error)
#else
    CLOSE(this%unit,IOSTAT=this%error)
#endif
  END SUBROUTINE CloseFile_gnuplot

  !> \public Writes the configuration as a header to the file
  !!
  SUBROUTINE WriteHeader(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
       CALL MPI_File_write(this%handle,TRIM(header_buf),LEN(TRIM(header_buf)), &
            MPI_CHARACTER,this%status,this%error)
#else
       WRITE (this%unit,FMT='(A)',IOSTAT=this%error) TRIM(header_buf) !(1:HLEN-1)
#endif
    END IF
  END SUBROUTINE WriteHeader

  !> \public Reads the header (not yet implemented)
  !!
  SUBROUTINE ReadHeader_gnuplot(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this     !< \param [in,out] this fileio type
    LOGICAL            :: success  !< \param [out] success
    !------------------------------------------------------------------------!
    INTENT(OUT)        :: success
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (GetRank(this).EQ.0) THEN
    END IF
#else
#endif
    success = .FALSE.
  END SUBROUTINE ReadHeader_gnuplot

  !> \public Writes the timestep (not yet implemented)
  !!
  SUBROUTINE WriteTimestamp_gnuplot(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this   !< \param [in,out] this fileio type
    REAL               :: time   !< \param [in] time
    !------------------------------------------------------------------------!
    INTENT(IN)         :: time
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
#else
#endif
    END IF
!!$    CALL Warning(this,"WriteTimestamp_gnuplot",&
!!$         "function is not implemented")
  END SUBROUTINE WriteTimestamp_gnuplot

  !> \public Reads the timestep (not yet implemented)
  !!
  SUBROUTINE ReadTimestamp_gnuplot(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this  !< \param [in,out] this fileio type
    REAL               :: time  !< \param [out] time
    !------------------------------------------------------------------------!
    INTENT(OUT)        :: time
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    IF (GetRank(this).EQ.0) THEN
#ifdef PARALLEL
#else
#endif
    END IF
    time = 0.0
!!$    CALL Warning(this,"ReadTimestamp_gnuplot",&
!!$         "function is not implemented")
  END SUBROUTINE ReadTimestamp_gnuplot

  !> \public Writes all desired data arrays to a file
  !!
  SUBROUTINE WriteDataset(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base) :: this !< \param [in,out] this fileio type
    CLASS(mesh_base)   :: Mesh !< \param [in] mesh mesh type
    !------------------------------------------------------------------------!
    INTEGER             :: i,j,k,l
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER                       :: request
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!

    IF (ASSOCIATED(Timedisc%w)) THEN
      IF (Mesh%FARGO.EQ.3.AND.Mesh%SN_shear) THEN
        CALL Physics%AddBackgroundVelocityX(Mesh,Timedisc%w,Timedisc%pvar,Timedisc%cvar)
      ELSE IF(Mesh%geometry%GetType().EQ.SPHERICAL) THEN
        CALL Physics%AddBackgroundVelocityZ(Mesh,Timedisc%w,Timedisc%pvar,Timedisc%cvar)
      ELSE
        CALL Physics%AddBackgroundVelocityY(Mesh,Timedisc%w,Timedisc%pvar,Timedisc%cvar)
      END IF
    END IF


#ifdef PARALLEL
    ! be sure to write at the end by getting the offset from the file's size
    CALL MPI_File_get_size(this%handle,offset,this%error)
    ! very importan
    CALL MPI_Barrier(MPI_COMM_WORLD,this%error)
    ! write _one_ line feed at the beginning of each time step
    IF (GetRank(this).EQ.0) THEN
       CALL MPI_File_write_at(this%handle, offset, LF, 1, MPI_CHARACTER, &
            this%status, this%error)
    END IF
    ! add the initial line feed and the general offset (depends on Mesh%IMIN)
    offset = offset + 1
    ! create the file view
    CALL MPI_File_set_view(this%handle,offset,this%basictype,this%filetype, &
         'native',MPI_INFO_NULL,this%error)
#else
    ! write _one_ line feed at the beginning of each time step
    WRITE (this%unit,FMT='(A)',ADVANCE='NO') LF
#endif

    DO k=1,this%cols
      ! trim the data for gnuplot output
      WHERE (ABS(this%output(k)%val(:,:)).LT.MAX(TINY(this%output(k)%val(:,:)),1.0D-99))
        this%output(k)%val(:,:) = 0.0E+00
      END WHERE
    END DO

    ! start i,j from 1 to rank local size, because RemapBounds has not been
    ! used/cannot be used.
    DO i=1,this%inum
       DO j=1,this%jnum
          ! write positions to line buffer
          DO k=1,this%COLS-1
             WRITE (this%linebuf((k-1)*this%FLEN+1:k*this%FLEN),TRIM(this%fmtstr)) &
                   this%output(k)%val(i,j), RECSEP
          END DO

          IF ((j.EQ.Mesh%JNUM).AND.((Mesh%JNUM.GT.1).OR.(Mesh%INUM.EQ.i))) THEN
             ! finish the block
             WRITE (this%linebuf((this%COLS-1)*this%FLEN+1:this%linelen),TRIM(this%fmtstr)) &
                  this%output(this%COLS)%val(i,j), BLKSEP
          ELSE
             ! finish the line
             WRITE (this%linebuf((this%COLS-1)*this%FLEN+1:this%linelen),TRIM(this%fmtstr)) &
                  this%output(this%COLS)%val(i,j), LINSEP
          END IF

#ifdef PARALLEL
          ! write line buffer to output buffer
          DO k=1,this%linelen
             this%outbuf(k,j-1+Mesh%JMIN) = this%linebuf(k:k)
          END DO
#else
          ! write line buffer to output file
          WRITE (this%unit,FMT=TRIM(this%linefmt),ADVANCE='YES') this%linebuf(1:this%linelen-1)
#endif
       END DO
#ifdef PARALLEL
       !*****************************************************************!
       ! This collective call doesn't work for pvfs2 -> bug in ROMIO ?
!!$       CALL MPI_File_write_all(this%handle,this%binout,this%bufsize,&
!!$            this%basictype, this%status, this%error)
       !*****************************************************************!
       ! so we use these two commands instead
       CALL MPI_File_iwrite(this%handle,this%outbuf,this%bufsize,this%basictype,&
            request,this%error)
       CALL MPI_Wait(request,this%status,this%error)
#endif
    END DO
  END SUBROUTINE WriteDataset

  !> \public Reads the data arrays from file (not yet implemented)
  !!
  SUBROUTINE ReadDataset_gnuplot(this,Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base)   :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base)     :: Mesh     !< \param [in] mesh mesh type
    CLASS(physics_base)  :: Physics  !< \param [in] physics physics type
    CLASS(timedisc_base) :: Timedisc !< \param [in] timedisc timedisc type
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics,Timedisc
    INTENT(INOUT)        :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE ReadDataset_gnuplot


  !> Closes the file I/O and calls a further error function
  !!
  SUBROUTINE Error(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), INTENT(INOUT) :: this    !< \param [in,out] this fileio type
    CHARACTER(LEN=*),  INTENT(IN)     :: modproc !< \param [in] modproc
    CHARACTER(LEN=*),  INTENT(IN)     :: msg     !< \param [in] msg error msg
    !------------------------------------------------------------------------!
    IF (Initialized(this)) &
         CALL this%CloseFile_gnuplot()
    CALL this%Error_fileio(modproc,msg)
  END SUBROUTINE Error

  !> \public Closes the file I/O
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    DEALLOCATE(this%outbuf)
#endif
    DEALLOCATE(this%output)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE fileio_gnuplot_mod
