!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_gnuplot.f90                                                #
!#                                                                           #
!# Copyright (C) 2008-2024                                                   #
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
!! \extends fileio_base
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
  ! Private Attributes section starts here:
  !> \name Private Attributes
  !>#### some default limits
  INTEGER, PARAMETER      :: HLEN = 10000         !< header length in bytes
  INTEGER, PARAMETER      :: DEFAULT_DECS = 5     !< default decimal places
  !>#### special strings
  CHARACTER, PARAMETER    :: SP = ACHAR(32)       !< space
  CHARACTER, PARAMETER    :: LF = ACHAR(10)       !< line feed
  CHARACTER(LEN=2), PARAMETER  :: RECSEP = SP // SP    !< data record separator
  CHARACTER(LEN=2), PARAMETER  :: LINSEP = SP // LF    !< line separator
  CHARACTER(LEN=2), PARAMETER  :: BLKSEP = LF // LF    !< block separator
  !> the header string
  CHARACTER(LEN=30), PARAMETER :: &
         header_string = "# Data output of fosite" // LINSEP
  CHARACTER(LEN=HLEN)  :: header_buf              !< buffer of header
  !--------------------------------------------------------------------------!
  !> output-pointer for 3D array data
  TYPE ValPtr_TYP
    REAL, DIMENSION(:,:,:), POINTER :: val
  END TYPE

  TYPE Output_TYP
    TYPE(ValPtr_TYP), DIMENSION(:), POINTER :: p
    CHARACTER(LEN=MAX_CHAR_LEN)             :: key
    CHARACTER(LEN=1024)                     :: path
    INTEGER                                 :: numbytes
  END TYPE Output_TYP

  !> output-pointer for time step scalar data
  TYPE TSOutput_TYP
    REAL, POINTER                           :: val
    CHARACTER(LEN=MAX_CHAR_LEN)             :: key
  END TYPE TSOutput_TYP

  !> FileIO gnuplot class
  TYPE, EXTENDS(fileio_base) :: fileio_gnuplot
    !> \name Variables
    TYPE(Output_TYP),DIMENSION(:), POINTER :: &     !< list of output fields
                              output => null()
    TYPE(TSOutput_TYP),DIMENSION(:), POINTER :: &   !< list of scalar time step output
                              tsoutput => null()
    CHARACTER(LEN=1), DIMENSION(:,:), POINTER :: &
                              outbuf => null()      !< output buffer
    CHARACTER(LEN=512)     :: heading     !< char buffer for heading (field data)
    CHARACTER(LEN=512)     :: tsheading   !< char buffer for heading (time step data)
    CHARACTER(LEN=512)     :: linebuf     !< char buffer fo field data
    CHARACTER(LEN=512)     :: tslinebuf   !< char buffer for time step data
    CHARACTER(LEN=64)      :: fmtstr      !< format string
    INTEGER                :: COLS        !< number of output columns
    INTEGER                :: TSCOLS      !< number of output columns for time step data
    INTEGER                :: MAXCOLS     !< upper limit for output cols
                                          !! (MAXCOLS < LEN(linebuf)/FLEN)
    INTEGER                :: DECS        !< decimal places for real number output
    INTEGER                :: FLEN        !< output field length
    INTEGER                :: linelen     !< length of a line
    INTEGER                :: tslinelen   !< length of a line for time step data
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFileIO_deferred => InitFileIO_gnuplot
    PROCEDURE :: WriteHeader
    PROCEDURE :: ReadHeader
    PROCEDURE :: WriteDataset_deferred => WriteDataset_gnuplot
!     PROCEDURE :: ReadDataset
    PROCEDURE :: GetOutputList
    FINAL :: Finalize
  END TYPE

  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fileio_gnuplot, Output_TYP, TSOutput_TYP, ValPtr_TYP
 !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for the GNUPLOT file I/O
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs,
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_gnuplot(this,Mesh,Physics,Timedisc,Sources,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot),INTENT(INOUT)      :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),    INTENT(IN)          :: Mesh     !< \param [in] Mesh mesh type
    CLASS(physics_base), INTENT(IN)          :: Physics  !< \param [in] Physics physics type
    CLASS(timedisc_base),INTENT(IN)          :: Timedisc !< \param [in] Timedisc timedisc type
    CLASS(sources_base), INTENT(IN), POINTER :: Sources  !< \param [in] Sources sources type
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: config   !< \param [in] config dict with I/O configuration
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: IO       !< \param [in] IO dict with pointers to I/O arrays
    !------------------------------------------------------------------------!
    TYPE(Output_TYP),DIMENSION(:),POINTER     :: poutput
    TYPE(Dict_TYP), POINTER                   :: node
    CHARACTER(LEN=MAX_CHAR_LEN),DIMENSION(4)  :: skip
    TYPE(real_t)                              :: dummy1
    REAL, DIMENSION(:,:,:,:), POINTER         :: dummy4
    INTEGER                                   :: cartcoords
    INTEGER                                   :: depth
    INTEGER                                   :: i,j,k,n
#ifdef PARALLEL
    INTEGER                                   :: blocknum
    INTEGER, DIMENSION(Mesh%JMAX-Mesh%JMIN+1,Mesh%KMAX-Mesh%KMIN+1) :: blocklen,indices
#endif
    !------------------------------------------------------------------------!
    CALL this%InitFileio(Mesh,Physics,Timedisc,Sources,config,IO,"gnuplot","dat",textfile=.FALSE.)

    CALL GetAttr(config, "decimals", this%DECS, DEFAULT_DECS)
    ! compute length of character field for real number output
    ! and check if linebuffer is large enough
    ! flen = 1 (sign) + 1 (one digit) + 1 (decimal point) + decs (decimal places)
    !      + 1 (E character) + 1 (sign of exponent) + 2 (exponent digits) + 2 (spaces)
    this%FLEN = this%DECS + 9
    this%MAXCOLS = LEN(this%linebuf)/this%FLEN-1

    ! this is for registering all arrays which are supposed to be written to the output file
    ALLOCATE(this%output(this%MAXCOLS),this%tsoutput(this%MAXCOLS),STAT=this%err)
    IF (this%err.NE.0) &
      CALL this%Error("fileio_gnuplot::InitFileIO","memory allocation failed for this%output, this%tsoutput")

    ! get pointer to simulation time to make it the first entry
    ! in the time step data set
    CALL GetAttr(IO,"/timedisc/time",dummy1)
    IF (ASSOCIATED(dummy1%p)) THEN
       this%tsoutput(1)%val => dummy1%p
       this%tsoutput(1)%key = "time"
       this%TSCOLS = 1
    ELSE
       this%TSCOLS = 0
    END IF

    ! check if cartesian coordinates are selected for output
    IF (this%cartcoords) THEN
       ! output cartesian coordinates
       CALL GetAttr(IO,"/mesh/bary_centers",dummy4)
       ! requires 3 pointers, one for each cartesian coordinate
       ALLOCATE(this%output(1)%p(3),STAT=this%err)
    ELSE
       ! output curvilinear coordinates
       CALL GetAttr(IO,"/mesh/bary_curv",dummy4)
       ! requires NDIMS pointers, one for each curvilinear coordinate, depending on the
       ! dimensionality of the mesh
       ALLOCATE(this%output(1)%p(3),STAT=this%err)
      !!!!! only output necessary curvilinear coordinates; needs some addional checks below
!        ALLOCATE(this%output(1)%p(Mesh%NDIMS),STAT=this%err)
    END IF
    IF (this%err.NE.0) &
       CALL this%Error("fileio_gnuplot::InitFileIO","memory allocation failed for this%output(1)%p")

    ! register coordinate array pointers for output
    DO n=1,3
      this%output(1)%p(n)%val => dummy4(:,:,:,n)
    END DO
    this%output(1)%key = '# ' // 'x' // REPEAT(' ',this%FLEN-3) // 'y' // REPEAT(' ',this%FLEN-1) &
                              // 'z' // REPEAT(' ',this%FLEN-1)
    
    !!!!!! ATTENTION old code suppresses coordinate output in 1D/2D simulations
    !!!!!! if the coordinate does not vary along a specific dimension
!     IF (Mesh%INUM.EQ.1) THEN
!       ! use format string as temp
!       WRITE (this%fmtstr,'(A5,I2,A1)')'(A1,A',this%FLEN-3,')'
!       WRITE(this%linebuf,TRIM(this%fmtstr))'#','y'
!       this%COLS = 1
!     ELSE IF (Mesh%JNUM.EQ.1) THEN
!       WRITE (this%fmtstr,'(A5,I2,A1)')'(A1,A',this%FLEN-3,')'
!       WRITE(this%linebuf,TRIM(this%fmtstr))'#','x'
!       this%COLS = 1
!     ELSE
!       WRITE (this%fmtstr,'(A5,I2,A2,I2,A1)')'(A1,A',this%FLEN-3,',A',this%FLEN-1,')'
!       WRITE(this%linebuf,TRIM(this%fmtstr))'#','x','y'
!       this%COLS = 2
!     END IF

    ! register scalars and arrays for output
    node => IO
    skip(1:4) = [CHARACTER(LEN=MAX_CHAR_LEN) :: "bary_centers", "bary_curv", "corners", "time"]
    k = 1
    CALL this%GetOutputList(Mesh,node,k,this%TSCOLS,skip)

    ! shrink this%output
    poutput => this%output
    ALLOCATE(this%output,SOURCE=poutput(1:k),STAT=this%err)
    IF (this%err.NE.0) &
      CALL this%Error("fileio_gnuplot::InitFileIO","memory allocation failed for this%output")
    DEALLOCATE(poutput)

    ! count number of output columns for array data
    this%COLS = 0
    DO n=1,SIZE(this%output)
       this%COLS = this%COLS + SIZE(this%output(n)%p)
    END DO

    ! length of one output line
    this%linelen = this%COLS * this%FLEN
    IF (this%linelen.GT.LEN(this%linebuf)) &
       CALL this%Error("fileio_gnuplot::InitFileIO", &
          "linebuffer to small; reducing decimals or number of output fields may help")

    ! create heading for field data
    this%linebuf = ""
    ! copy heading for coordinates first
    n = SIZE(this%output(1)%p)*this%FLEN
    this%linebuf(1:n-1) = this%output(1)%key(1:n-1)
    ! format string for writing one string of width FLEN
    WRITE(this%fmtstr,'(A,I2,A1)') '(A',this%FLEN,')'
    DO k=2,SIZE(this%output)  ! skip entry for coordinates, i.e. this%output(1)
       i = INDEX(this%output(k)%key,"/",BACK=.TRUE.) ! find last occurance of slash in key
       ! write key (without leading path)
       WRITE(this%linebuf(n:),TRIM(this%fmtstr)) this%output(k)%key(i+1:i+this%FLEN-1)
       ! append spaces
       DO i=2,SIZE(this%output(k)%p)
         WRITE(this%linebuf(n+(i-1)*this%FLEN:),TRIM(this%fmtstr)) REPEAT(' ',this%FLEN)
       END DO
       n = n + (i - 1)*this%FLEN
    END DO
    this%heading = TRIM(this%linebuf) // LF

    ! create heading for time step data
    ! prepend hash -> comment line in GNUPLOT
    this%linebuf = ""
    DO k=1,this%TSCOLS
       WRITE(this%linebuf(1+(k-1)*this%FLEN:),TRIM(this%fmtstr)) TRIM(this%tsoutput(k)%key(1:this%FLEN-1))
    END DO
    this%tsheading = REPEAT("#",this%FLEN-1) // TRIM(this%linebuf(1:)) // LF

    ! local domain size
    this%INUM = Mesh%IMAX - Mesh%IMIN + 1
    this%JNUM = Mesh%JMAX - Mesh%JMIN + 1
    this%KNUM = Mesh%KMAX - Mesh%KMIN + 1

   ! size of the output buffer (on this process)
   this%bufsize  = this%INUM

   ! allocate memory for output buffer and displacement records
   this%err = 0
   ALLOCATE(this%outbuf(this%linelen,this%bufsize),STAT=this%err)
   IF (this%err.NE.0) &
      CALL this%Error("fileio_gnuplot::InitFileIO_gnuplot","memory allocation failed for this%outbuf")

#ifdef PARALLEL
   ! number of output blocks (on this process)
   ! = number of cells (excluding ghost cells)
   blocknum = this%JNUM * this%KNUM
   blocklen(:,:) = this%bufsize
   DO k=1,this%KNUM
      DO j=1,this%JNUM
          indices(j,k) = ((k+Mesh%KMIN-2)*Mesh%JNUM + (j+Mesh%JMIN-2) )*Mesh%INUM + Mesh%IMIN - 1
      END DO
   END DO

   ! create new MPI data types
   SELECT TYPE(df=>this%datafile)
   CLASS IS(filehandle_mpi)
      ! basic type is one output line corresponding to one data point (coordinates + data)
      CALL MPI_Type_contiguous(this%linelen,MPI_CHARACTER,df%basictype,this%err)
      IF (this%err.EQ.0) CALL MPI_Type_commit(df%basictype,this%err)

      ! file type for blocknum indexed blocks of basic type
      IF (this%err.EQ.0) CALL MPI_Type_indexed(blocknum,RESHAPE(blocklen,(/blocknum/)), &
            RESHAPE(indices,(/blocknum/)),df%basictype,df%filetype,this%err)
      IF (this%err.EQ.0) CALL MPI_Type_commit(df%filetype,this%err)
   END SELECT
   IF (this%err.NE.0) &
      CALL this%Error("fileio_gnuplot::InitFileIO_gnuplot","creation of MPI file types failed")

#endif
    ! write the format string for one entry in the data file:
    ! FLEN-2 characters for the number and 2 for the separators
    WRITE (this%fmtstr,'(A3,I2,A,I2.2,A5)') '(ES', this%FLEN-2, '.', this%DECS,',A2)'

    ! print some information
    ! ...
  END SUBROUTINE InitFileIO_gnuplot


  !> Creates a string with the configuration (from the dictionary)
  !!
  RECURSIVE SUBROUTINE GetHeaderString(string,root,k,prefix)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER    :: root,node,subnode
    CHARACTER(LEN=*)           :: string
    CHARACTER(LEN=*), OPTIONAL :: prefix
    CHARACTER(LEN=128)         :: buf
    !------------------------------------------------------------------------!
    INTEGER                    :: idummy, k
    LOGICAL                    :: ldummy
    CHARACTER(LEN=MAX_CHAR_LEN):: cdummy
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
          CALL GetAttr(node,GetKey(node),cdummy)
          WRITE(buf,'(A1,A25,A,A)')'#',TRIM(GetKey(node))//": ",TRIM(cdummy), LINSEP
          WRITE(string(k:),'(A)')buf
          k = k + LEN(TRIM(buf))
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
          CALL GetHeaderString(string,subnode,k,TRIM(buf))
       END IF
       node => GetNext(node)
    END DO
  END SUBROUTINE GetHeaderString

  !> Creates a list of all data arrays which will be written to file
  !!
  !! Therefore it ignores all arrays with coordinates and checks if the data
  !! arrays are of the dimension of the mesh.
  RECURSIVE SUBROUTINE GetOutputList(this,Mesh,node,oarr,onum,skip,prefix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)         :: Mesh  !< \param [in] mesh mesh type
    TYPE(Dict_TYP), POINTER              :: node  !< \param [in,out] node pointer to (sub-)dict
    INTEGER, INTENT(INOUT)               :: oarr  !< \param [in,out] oarr number of output arrays
    INTEGER, INTENT(INOUT)               :: onum  !< \param [in,out] onum number of output numbers
    CHARACTER(LEN=MAX_CHAR_LEN), DIMENSION(:), OPTIONAL, INTENT(IN) &
                                         :: skip  !< \param [in] skip list of keys to skip
    CHARACTER(LEN=*), OPTIONAL, INTENT(INOUT) &
                                         :: prefix!<\param [in,out] prefix namespace (path) to sub-dict
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER              :: dir
    CHARACTER(LEN=MAX_CHAR_LEN)          :: key
    TYPE(real_t)                         :: dummy1
    REAL, DIMENSION(:,:), POINTER        :: dummy2
    REAL, DIMENSION(:,:,:), POINTER      :: dummy3
    REAL, DIMENSION(:,:,:,:), POINTER    :: dummy4
    REAL, DIMENSION(:,:,:,:,:), POINTER  :: dummy5
    INTEGER, DIMENSION(5)                :: dims
    INTEGER                              :: m,n
    !------------------------------------------------------------------------!
    ! reset error code
    DO WHILE(ASSOCIATED(node))
       ! check if node is directory
       IF(HasChild(node)) THEN
          ! recursion
          IF (PRESENT(prefix)) THEN
             ! add prefix to key
             key = TRIM(prefix)//'/'//TRIM(GetKey(node))
          ELSE
             key = '/'//TRIM(GetKey(node))
          END IF
          dir => GetChild(node)
          IF (PRESENT(skip)) THEN
             CALL GetOutputList(this,Mesh,dir,oarr,onum,skip,key)
          ELSE
             CALL GetOutputList(this,Mesh,dir,oarr,onum,prefix=key)
          END IF
       ELSE IF (HasData(node)) THEN
          ! node contains data
          IF (PRESENT(skip)) THEN
             IF (ANY(skip(:) == GetKey(node))) THEN
                ! skip entry and continue with next node
                node=>GetNext(node)
                CYCLE
             END IF
          END IF
          ! check shape of data
          dims = 0
          SELECT CASE(GetDatatype(node))
!           CASE(DICT_REAL_TWOD)
!               CALL GetAttr(node,TRIM(GetKey(node)),dummy2)
!               dims(1:2) = SHAPE(dummy2)
!               dims(3:5) = 1
          CASE(DICT_REAL_THREED)
              CALL GetAttr(node,TRIM(GetKey(node)),dummy3)
              dims(1:3) = SHAPE(dummy3)
              dims(4:5) = 1
          CASE(DICT_REAL_FOURD)
              CALL GetAttr(node,TRIM(GetKey(node)),dummy4)
              dims(1:4) = SHAPE(dummy4)
              dims(5)   = 1
          CASE(DICT_REAL_FIVED)
              CALL GetAttr(node,TRIM(GetKey(node)),dummy5)
              dims(1:5) = SHAPE(dummy5)
          CASE(DICT_REAL_P)
              CALL GetAttr(node,GetKey(node),dummy1)
              IF (ASSOCIATED(dummy1%p)) THEN
                onum=onum+1
                IF (onum.GT.this%MAXCOLS-1) &
                   CALL this%Error("fileio_gnuplot::GetOutputList", &
                            "number of scalar output fields exceeds upper limit")
                ! register number data for output
                this%tsoutput(onum)%val => dummy1%p
                this%tsoutput(onum)%key = TRIM(GetKey(node))
                IF (PRESENT(prefix)) THEN
                    this%tsoutput(onum)%key = TRIM(prefix) // '/' // this%tsoutput(onum)%key
                END IF
              END IF
          CASE DEFAULT
              CALL this%Warning("fileio_gnuplot::GetOutputList", &
                       "'" // GetKey(node) // "'" // " registered for output," &
                       // " but data type is currently not supported")
          END SELECT
          ! only register array data for output if it contains
          ! mesh data
          IF ((dims(1).EQ.(Mesh%IMAX-Mesh%IMIN+1)).AND.&
             (dims(2).EQ.(Mesh%JMAX-Mesh%JMIN+1)).AND.&
             (dims(3).EQ.(Mesh%KMAX-Mesh%KMIN+1))) THEN
             ! count number of output columns currently registered
             n=0
             DO m=1,oarr
                n = n + SIZE(this%output(m)%p)
             END DO
             ! check limit for output of array data
             IF (n+dims(4)*dims(5).GT.this%MAXCOLS) &
                CALL this%Error("fileio_gnuplot::GetOutputList", &
                         "number of array output fields exceeds upper limit")
             ! increase array output index
             oarr = oarr + 1
             ! store name of the output array
             this%output(oarr)%key = TRIM(GetKey(node))
             ! prepend prefix if present
             IF (PRESENT(prefix)) THEN
                this%output(oarr)%key = TRIM(prefix) // '/' // this%output(oarr)%key
             END IF
             ! allocate memory for pointer to output arrays
             ALLOCATE(this%output(oarr)%p(dims(4)*dims(5)),STAT=this%err)
             IF (this%err.NE.0) &
                CALL this%Error( "fileio_gnuplot::GetOutputList", "Unable to allocate memory.")
             ! set pointer to output arrays
             IF (dims(4).EQ.1) THEN
                ! 3D data
                this%output(oarr)%p(1)%val => dummy3
             ELSE IF (dims(5).EQ.1) THEN
                ! 4D data (for example pvars - density,xvel,yvel)
                DO n=1,dims(4)
                   this%output(oarr)%p(n)%val => dummy4(:,:,:,n)
                END DO
             ELSE
                ! 5D data, e.g., stress tensor components
                DO n=1,dims(5)
                   DO m=1,dims(4)
                      this%output(oarr)%p(n+(m-1)*dims(5))%val => dummy5(:,:,:,m,n)
                   END DO
                END DO
             END IF
          END IF

       END IF
       node=>GetNext(node)
    END DO
  END SUBROUTINE GetOutputList


  !> \public Writes the configuration as a header to the file
  !!
  SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base),      INTENT(IN)    :: Mesh      !< \param [in] Mesh mesh type
    CLASS(physics_base),   INTENT(IN)    :: Physics   !< \param [in] Physics physics type
    TYPE(Dict_TYP), POINTER              :: Header,IO
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER              :: node
    INTEGER                              :: depth
    !------------------------------------------------------------------------!
    IF (this%GetRank().EQ.0) THEN
       ! get settings from config dictionary and store in header_buf
       depth = 1
       node => Header
       CALL GetHeaderString(header_buf,node,depth)
       header_buf = TRIM(header_string) // TRIM(header_buf)

       SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
       CLASS IS(filehandle_fortran)
          WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) TRIM(header_buf) !(1:HLEN-1)
#else
       CLASS IS(filehandle_mpi)
          CALL MPI_File_write(df%GetUnitNumber(),TRIM(header_buf),LEN(TRIM(header_buf)), &
                MPI_CHARACTER,df%status,this%err)
#endif
       END SELECT
    END IF
  END SUBROUTINE WriteHeader

  !> \public Reads the header (not implemented)
  !!
  SUBROUTINE ReadHeader(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot) :: this     !< \param [in,out] this fileio type
    LOGICAL            :: success  !< \param [out] success
    !------------------------------------------------------------------------!
    INTENT(OUT)        :: success
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (this%GetRank().EQ.0) THEN
    END IF
#else
#endif
    CALL this%Warning("fileio_gnuplot::ReadHeader","reading file header not implemented yet")
    success = .FALSE.
  END SUBROUTINE ReadHeader

  !> \public Writes all desired data arrays to a file
  !!
  SUBROUTINE WriteDataset_gnuplot(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)        :: Mesh      !< \param [in] mesh mesh type
    CLASS(physics_base), INTENT(INOUT)  :: Physics   !< \param [in] physics physics type
    CLASS(fluxes_base), INTENT(IN)      :: Fluxes    !< \param [in] fluxes fluxes type
    CLASS(timedisc_base), INTENT(IN)    :: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP), POINTER             :: Header,IO !< \param [in,out] IO I/O dictionary
    INTEGER             :: i,j,k,m,l,n
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER                       :: request, status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
    ! some sanity checks:
    ! 1. no output data
    IF (.NOT.ASSOCIATED(this%output)) RETURN

    ! 2. uninitilized output data structure
    DO l=1,SIZE(this%output)
       IF (.NOT.ASSOCIATED(this%output(l)%p)) THEN
          CALL this%Error("fileio_gnuplot::WriteDataset", &
                          "this should not happen: output data pointer for " &
                          // TRIM(this%output(l)%key) // " not associated")
       ELSE
          DO m=1,SIZE(this%output(l)%p)
             IF (.NOT.ASSOCIATED(this%output(l)%p(m)%val)) &
                CALL this%Error("fileio_gnuplot::WriteDataset", &
                                "this should not happen: one of the data array pointers for " &
                                // TRIM(this%output(l)%key) // " not associated")
          END DO
       END IF
    END DO

    this%err = 0
    IF (this%GetRank().EQ.0) THEN
       ! generate string with time step data
       DO k=1,this%TSCOLS
          WRITE(this%tslinebuf(1+(k-1)*this%FLEN:),TRIM(this%fmtstr)) this%tsoutput(k)%val
       END DO
       this%tslinebuf = REPEAT("#",this%FLEN-1) // SP // TRIM(this%tslinebuf) // LF

       ! write time step data to file
       SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
       CLASS IS(filehandle_fortran)
          WRITE (df%GetUnitNumber(),IOSTAT=this%err) TRIM(this%tsheading) &
              // TRIM(this%tslinebuf) // TRIM(this%heading)
#else
       CLASS IS(filehandle_mpi)
          CALL MPI_File_write(df%GetUnitNumber(),TRIM(this%tsheading),LEN(TRIM(this%tsheading)), &
                MPI_CHARACTER,df%status,this%err)
          IF (this%err.EQ.0) CALL MPI_File_write(df%GetUnitNumber(),TRIM(this%tslinebuf), &
                LEN(TRIM(this%tslinebuf)),MPI_CHARACTER,df%status,this%err)
          IF (this%err.EQ.0) CALL MPI_File_write(df%GetUnitNumber(),TRIM(this%heading), &
                LEN(TRIM(this%heading)),MPI_CHARACTER,df%status,this%err)
#endif
       END SELECT
    END IF

    IF (this%err.NE.0) &
       CALL this%Error("fileio_gnuplot::WriteDataset","writing time step data to file failed")

    SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
    CLASS IS(filehandle_fortran)
      ! write _one_ line feed at the beginning of each time step
      WRITE (df%GetUnitNumber(),IOSTAT=this%err) LF
#else
    CLASS IS(filehandle_mpi)
      ! very important: wait for the header write command to finish
      CALL MPI_Barrier(MPI_COMM_WORLD,this%err)
      ! be sure to write at the end, i.e. after the time step data, by getting the offset from the file's size
      IF (this%err.EQ.0) CALL MPI_File_get_size(df%GetUnitNumber(),offset,this%err)
      ! write _one_ line feed at the beginning of each time step
      IF (this%GetRank().EQ.0) THEN
        IF (this%err.EQ.0) CALL MPI_File_write_at(df%GetUnitNumber(), offset, LF, 1, &
             MPI_CHARACTER, df%status, this%err)
      END IF
#endif
    END SELECT

    IF (this%err.NE.0) &
        CALL this%Error("fileio_gnuplot::WriteDataset","writing preceeding line feed to file failed")

#ifdef PARALLEL
    SELECT TYPE(df=>this%datafile)
    CLASS IS(filehandle_mpi)
      ! add the initial line feed and the general offset (depends on Mesh%IMIN)
      offset = offset + 1
      ! create the file view
      IF (this%err.EQ.0) CALL MPI_File_set_view(df%GetUnitNumber(),offset,df%basictype,df%filetype, &
          'native',MPI_INFO_NULL,this%err)
    END SELECT

    IF (this%err.NE.0) &
         CALL this%Error("fileio_gnuplot::WriteDataset","creating MPI file view failed")
#endif

    ! trim the floating point values for gnuplot output, i.e. set small numbers to 0
    DO l=1,SIZE(this%output)
      DO m=1,SIZE(this%output(l)%p)
        WHERE (ABS(this%output(l)%p(m)%val(:,:,:)).LT.MAX(TINY(this%output(l)%p(m)%val),1.0D-99))
          this%output(l)%p(m)%val(:,:,:) = 0.0E+00
        END WHERE
      END DO
    END DO

    ! write array data to file
    ! do not use Mesh%IMIN/MAX etc., because the indices of this%output()%p()%val start at 1
    this%err = 0
    DO k=Mesh%KMIN,Mesh%KMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          ! write array data to line buffer
          n = 1
          DO l=1,SIZE(this%output)
            DO m=1,SIZE(this%output(l)%p)
               ! fill output line buffer with data
               WRITE (this%linebuf((n-1)*this%FLEN+1:n*this%FLEN),TRIM(this%fmtstr)) &
                     this%output(l)%p(m)%val(i-Mesh%IMIN+1,j-Mesh%JMIN+1,k-Mesh%KMIN+1), RECSEP
               n = n + 1
            END DO
          END DO

          IF (Mesh%INUM.GT.1) THEN
             IF (i.EQ.Mesh%INUM) THEN
                ! finish the block
                this%linebuf(this%linelen-1:this%linelen) = BLKSEP
             ELSE
                ! finish the line
                this%linebuf(this%linelen-1:this%linelen) = LINSEP
             END IF
          ELSE IF (Mesh%JNUM.GT.1) THEN
             IF (j.EQ.Mesh%JNUM) THEN
                ! finish the block
                this%linebuf(this%linelen-1:this%linelen) = BLKSEP
             ELSE
                ! finish the line
                this%linebuf(this%linelen-1:this%linelen) = LINSEP
             END IF
          ELSE
             IF (k.EQ.Mesh%KNUM) THEN
                ! finish the block
                this%linebuf(this%linelen-1:this%linelen) = BLKSEP
             ELSE
                ! finish the line
                this%linebuf(this%linelen-1:this%linelen) = LINSEP
             END IF
          END IF

          ! write line buffer to output buffer
          this%outbuf(:,i-Mesh%IMIN+1) = TRANSFER(this%linebuf,this%outbuf(:,1),SIZE=this%linelen)
        END DO ! i-loop

        ! write output buffer to file
        SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
        CLASS IS(filehandle_fortran)
           ! write line buffer to file
           IF (this%err.EQ.0) &
              WRITE (df%GetUnitNumber(),IOSTAT=this%err) this%outbuf
#else
        CLASS IS(filehandle_mpi)
          !*****************************************************************!
          ! This collective call didn't work for pvfs2 -> bug in ROMIO ?
          IF (this%err.EQ.0) CALL MPI_File_write_all(df%GetUnitNumber(),this%outbuf, &
                                      this%bufsize,df%basictype,status,this%err)
          !*****************************************************************!
          ! If the collective call above fails try this:
!           IF (this%err.EQ.0) CALL MPI_File_iwrite(df%GetUnitNumber(),this%outbuf, &
!                                      this%bufsize,df%basictype,request,this%err)
!           IF (this%err.EQ.0) CALL MPI_Wait(request,df%status,this%err)
#endif
        END SELECT

      END DO ! j-loop
    END DO ! k-loop

    IF (this%err.NE.0) &
         CALL this%Error("fileio_gnuplot::WriteDataset","creating MPI file view failed")

  END SUBROUTINE WriteDataset_gnuplot

!   !> \public Reads the data arrays from file (not yet implemented)
!   !!
!   SUBROUTINE ReadDataset(this,Mesh,Physics,Timedisc)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(fileio_gnuplot)   :: this     !< \param [in,out] this fileio type
!     CLASS(mesh_base)     :: Mesh     !< \param [in] mesh mesh type
!     CLASS(physics_base)  :: Physics  !< \param [in] physics physics type
!     CLASS(timedisc_base) :: Timedisc !< \param [in] timedisc timedisc type
!     !------------------------------------------------------------------------!
!     INTENT(IN)           :: Mesh,Physics,Timedisc
!     INTENT(INOUT)        :: this
!     !------------------------------------------------------------------------!
!   END SUBROUTINE ReadDataset


  !> Closes the file I/O and calls a further error function
  !!
  SUBROUTINE Error(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot), INTENT(INOUT) :: this    !< \param [in,out] this fileio type
    CHARACTER(LEN=*),  INTENT(IN)     :: modproc !< \param [in] modproc
    CHARACTER(LEN=*),  INTENT(IN)     :: msg     !< \param [in] msg error msg
    !------------------------------------------------------------------------!
    IF (this%Initialized()) &
      CALL this%datafile%CloseFile(this%step)
    CALL this%Error(modproc,msg)
  END SUBROUTINE Error

  !> \public Closes the file I/O
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fileio_gnuplot),INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTEGER                          :: k
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(this%outbuf)) DEALLOCATE(this%outbuf)
    IF (ASSOCIATED(this%output)) THEN
      DO k=1,SIZE(this%output)
        IF (ASSOCIATED(this%output(k)%p)) DEALLOCATE(this%output(k)%p)
      END DO
      DEALLOCATE(this%output)
    END IF
    IF (ASSOCIATED(this%tsoutput)) DEALLOCATE(this%tsoutput)

    NULLIFY(this%outbuf,this%output,this%tsoutput)

    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE fileio_gnuplot_mod
