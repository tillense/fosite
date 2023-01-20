!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_gnuplot.f90                                                #
!#                                                                           #
!# Copyright (C) 2008-2023                                                   #
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
  ! Private Attributes section starts here:
  !> some string lengths
  INTEGER, PARAMETER      :: HLEN = 10000         !< header length
  INTEGER, PARAMETER      :: DEFAULT_DECS = 5     !< default decimal places
  !> some special strings
  CHARACTER, PARAMETER    :: SP = ACHAR(32)       !< space
  CHARACTER, PARAMETER    :: LF = ACHAR(10)       !< line feed
  CHARACTER*2, PARAMETER  :: RECSEP = SP // SP    !< data record separator
  CHARACTER*2, PARAMETER  :: LINSEP = SP // LF    !< line separator
  CHARACTER*2, PARAMETER  :: BLKSEP = LF // LF    !< block separator
  !> the header string
  CHARACTER(LEN=30), PARAMETER :: &
         header_string = "# Data output of fosite" // LINSEP
  CHARACTER(LEN=HLEN)  :: header_buf              !< buffer of header
  !--------------------------------------------------------------------------!
  !> FileIO gnuplot class
  TYPE, EXTENDS(fileio_base) :: fileio_gnuplot
    !> \name Variables
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
#ifdef PARALLEL
    !> \name Variables in parallel mode
    INTEGER                :: blocknum    !< number of output blocks
    INTEGER(KIND=MPI_ADDRESS_KIND) :: &
                              realext, &  !< real data type extent
                              intext      !< integer data type extent
    CHARACTER, DIMENSION(:,:), POINTER :: &
                              outbuf      !< output buffer
#endif
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFileIO_gnuplot
    PROCEDURE :: Finalize
    PROCEDURE :: WriteHeader
!    PROCEDURE :: ReadHeader
!    PROCEDURE :: WriteTimestamp
!    PROCEDURE :: ReadTimestamp
    PROCEDURE :: WriteDataset
!     PROCEDURE :: ReadDataset
    PROCEDURE :: GetOutputList
  END TYPE

  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fileio_gnuplot
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
    TYPE(Dict_TYP), POINTER                   :: node
    CHARACTER(LEN=128),DIMENSION(4)           :: skip
    REAL, DIMENSION(:,:,:,:), POINTER         :: dummy4
    INTEGER                                   :: cartcoords
    INTEGER                                   :: depth
    INTEGER                                   :: n
#ifdef PARALLEL
    INTEGER                                   :: i
    INTEGER, DIMENSION(Mesh%IMAX-Mesh%IMIN+1) :: blocklen,indices
#endif
    !------------------------------------------------------------------------!
    CALL this%InitFileio(Mesh,Physics,Timedisc,Sources,config,IO,"gnuplot","dat")

    CALL GetAttr(config, "/datafile/decimals", this%DECS, DEFAULT_DECS)
    ! compute length of character field for real number output
    ! and check if linebuffer is large enough
    ! flen = 1 (sign) + 1 (one digit) + 1 (decimal point) + decs (decimal places)
    !      + 1 (E character) + 1 (sign of exponent) + 2 (exponent digits) + 2 (spaces)
    this%FLEN = this%DECS + 9
    this%MAXCOLS = len(this%linebuf)/this%FLEN-1

    ALLOCATE(this%output(this%MAXCOLS),STAT=this%err)
    IF (this%err.EQ.0) &
      ! create pointer to coordinate arrays depending on mesh dimensionality
      ALLOCATE(this%output(1)%p(Mesh%NDIMS),STAT=this%err)
    IF (this%err.NE.0) &
      CALL this%Error("fileio_gnuplot::InitFileIO","memory allocation failed for this%output")

    ! check if cartesian coordinates are selected for output;
    ! default: curvilinear coordinates (0)
    CALL GetAttr(config, "/datafile/cartcoords", cartcoords, 0)
    depth = 1
    node => config
    CALL WriteHeaderString(header_buf,node,depth)
    IF (cartcoords.EQ.0) THEN
       CALL GetAttr(IO,"/mesh/bary_curv",dummy4)
    ELSE
       CALL GetAttr(IO,"/mesh/bary_centers",dummy4)
    END IF
    DO n=1,Mesh%NDIMS
       this%output(1)%p(n)%val => dummy4(:,:,:,n)
    END DO
    
    !!!!!! ATTENTION should be repaired
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

    ! set output-pointer and count the number of output columns
    WRITE (this%fmtstr,'(A,I2,A1)')'(A',this%FLEN-1,')'
    node => IO
    skip(1:4) = [CHARACTER(LEN=128) :: "bary_centers", "bary_curv", "corners", "time"]
    CALL this%GetOutputList(Mesh,skip,node,this%COLS,this%TSCOLS)

    ! length of one output line
    this%linelen = this%COLS * this%FLEN
    IF (this%linelen.GT.LEN(this%linebuf)) &
       CALL this%Error("fileio_gnuplot::InitFileIO", &
          "linebuffer to small; reducing decimals or number of output fields may help")

    header_buf = TRIM(header_string) // TRIM(header_buf) // this%linebuf(1:this%linelen)

    ! local domain size
    this%inum = Mesh%IMAX - Mesh%IMIN + 1
    this%jnum = Mesh%JMAX - Mesh%JMIN + 1
    this%knum = Mesh%KMAX - Mesh%KMIN + 1

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
  RECURSIVE SUBROUTINE GetOutputList(this,Mesh,skip,node,oarr,onum)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)         :: Mesh  !< \param [in] mesh mesh type
    CHARACTER(LEN=128), DIMENSION(:), INTENT(IN) &
                                         :: skip  !< \param [in] skip list of keys to skip
    TYPE(Dict_TYP), POINTER              :: node  !< \param [in,out] node pointer to (sub-)dict
    INTEGER, INTENT(INOUT)               :: oarr  !< \param [in,out] number of output arrays
    INTEGER, INTENT(INOUT)               :: onum  !< \param [in,out] number of output numbers
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER              :: dir
    TYPE(Real_t), POINTER                :: dummy1
    REAL, DIMENSION(:,:,:), POINTER      :: dummy3
    REAL, DIMENSION(:,:,:,:), POINTER    :: dummy4
    REAL, DIMENSION(:,:,:,:,:), POINTER  :: dummy5
    INTEGER, DIMENSION(5)                :: dims
    INTEGER                              :: m,n
    !------------------------------------------------------------------------!
    ! reset error code
    this%err = 0
    DO WHILE(ASSOCIATED(node).AND.this%err.EQ.0)
      ! check for directory and exclude any coordinates (these are handled elsewhere)
      IF(HasChild(node)) THEN
        ! recursion
        CALL GetAttr(node,GetKey(node),dir)
        CALL GetOutputList(this,Mesh,skip,dir,oarr,onum)
      ELSE IF (HasData(node).AND..NOT.ANY(skip(:) == GetKey(node))) THEN
        ! check data properties
!         SELECT CASE(GetDataType(node))
!         CASE(DICT_REAL_THREED)
!           ! real 3D array
!           CALL GetAttr(node,GetKey(node),dummy3)
!           dims(1:3) = SHAPE(dummy3)
!           ! check if it has mesh dimensions
!           IF((dims(1).EQ.Mesh%IMAX-Mesh%IMIN+1).AND.&
!              (dims(2).EQ.Mesh%JMAX-Mesh%JMIN+1).AND.&
!              (dims(3).EQ.Mesh%KMAX-Mesh%KMIN+1)) THEN
!             oarr = oarr + 1
!             IF (oarr .GT. this%MAXCOLS) THEN
!                this%err = 1
!                EXIT
!             END IF
!             this%output(oarr)%pval => dummy3
!             this%output(oarr)%key => GetKey(node)
!           END IF
!         CASE(DICT_REAL_FOURD)
!           ! real 4D array
!           CALL GetAttr(node,GetKey(node),dummy4)
!           dims(1:4) = SHAPE(dummy4)
!           IF((dims(1).EQ.Mesh%IMAX-Mesh%IMIN+1).AND.&
!              (dims(2).EQ.Mesh%JMAX-Mesh%JMIN+1).AND.&
!              (dims(3).EQ.Mesh%KMAX-Mesh%KMIN+1)) THEN
!             IF (oarr+dims(4) .GT. this%MAXCOLS) THEN
!                this%err = 1
!                EXIT
!             END IF
!             DO m=1,dims(4)
!               this%output(oarr+m)%val => dummy4(:,:,:,m)
!             END DO
!             this%output(oarr+1)%key = GetKey(node)
!             this%output(oarr+2:oarr+dims(4))%key = ""
!             oarr = oarr + dims(4)
!           END IF
!         CASE(DICT_REAL_FIVED)
!           ! real 5D array
!           CALL GetAttr(node,GetKey(node),dummy5)
!           dims(1:5) = SHAPE(dummy5)
!           IF((dims(1).EQ.Mesh%IMAX-Mesh%IMIN+1).AND.&
!              (dims(2).EQ.Mesh%JMAX-Mesh%JMIN+1).AND.&
!              (dims(3).EQ.Mesh%KMAX-Mesh%KMIN+1)) THEN
!             IF (oarr+dims(4)*dims(5) .GT. this%MAXCOLS) THEN
!                this%err = 1
!                EXIT
!             END IF
!             DO n=1,dims(5)
!               DO m=1,dims(4)
!                 this%output(oarr+n+(m-1)*dims(5))%val => dummy5(:,:,:,m,n)
!               END DO
!             END DO
!             this%output(oarr+1)%key = GetKey(node)
!             this%output(oarr+2:onum+dims(4)*dims(5))%key = ""
!             oarr = oarr + dims(4)*dims(5)
!           END IF          
!         CASE(DICT_REAL_P)
!           CALL GetAttr(node,GetKey(node),dummy1)
!           IF (ASSOCIATED(dummy1%p)) THEN
!              onum = onum + 1
!              this%tsoutput(onum)%val => dummy1%p
!              this%tsoutput(onum)%key = GetKey(node)
!           END IF
!         CASE DEFAULT
!           !do nothing (wrong type)
!         END SELECT
      END IF
      node=>GetNext(node)
    END DO
    IF (this%err.NE.0) &
         CALL Error(this,"fileio_gnuplot::GetOutputList","number of output fields exceeds upper limit")
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
    IF (this%GetRank().EQ.0) THEN
#ifdef PARALLEL
       CALL MPI_File_write(this%handle,TRIM(header_buf),LEN(TRIM(header_buf)), &
            MPI_CHARACTER,this%status,this%err)
#else
       WRITE (this%unit,FMT='(A)',IOSTAT=this%err) TRIM(header_buf) !(1:HLEN-1)
#endif
    END IF
  END SUBROUTINE WriteHeader

!   !> \public Reads the header (not yet implemented)
!   !!
!   SUBROUTINE ReadHeader_gnuplot(this,success)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(fileio_gnuplot) :: this     !< \param [in,out] this fileio type
!     LOGICAL            :: success  !< \param [out] success
!     !------------------------------------------------------------------------!
!     INTENT(OUT)        :: success
!     !------------------------------------------------------------------------!
! #ifdef PARALLEL
!     IF (GetRank(this).EQ.0) THEN
!     END IF
! #else
! #endif
!     success = .FALSE.
!   END SUBROUTINE ReadHeader_gnuplot

  !> \public Writes all desired data arrays to a file
  !!
  SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)        :: Mesh      !< \param [in] mesh mesh type
    CLASS(physics_base), INTENT(INOUT)  :: Physics   !< \param [in] physics physics type
    CLASS(fluxes_base), INTENT(IN)      :: Fluxes    !< \param [in] fluxes fluxes type
    CLASS(timedisc_base), INTENT(IN)    :: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP), POINTER             :: Header,IO !< \param [in,out] IO I/O dictionary
    INTEGER             :: i,j,k,m,l
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER                       :: request
#endif
    !------------------------------------------------------------------------!
    IF (ASSOCIATED(Timedisc%w)) THEN
      IF (Mesh%FARGO.EQ.3.AND.Mesh%shear_dir.EQ.1) THEN
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

    DO l=1,this%cols
      ! trim the data for gnuplot output
      DO m=1,SIZE(this%output(l)%p)
        WHERE (ABS(this%output(l)%p(m)%val(:,:,:)).LT.MAX(TINY(this%output(l)%p(m)%val),1.0D-99))
          this%output(l)%p(m)%val(:,:,:) = 0.0E+00
        END WHERE
      END DO
    END DO

    ! start i,j from 1 to rank local size, because RemapBounds has not been
    ! used/cannot be used.
    DO i=1,this%inum
      DO j=1,this%jnum
        DO k=1,this%knum   
          ! write positions to line buffer
          DO m=1,this%COLS-1
            DO l=1,SIZE(this%output(m)%p)
              WRITE (this%linebuf((m-1)*this%FLEN+1:m*this%FLEN),TRIM(this%fmtstr)) &
                     this%output(m)%p(l)%val(i,j,k), RECSEP
            END DO
          END DO
          
          m=this%COLS
          DO l=1,SIZE(this%output(m)%p)-1
            WRITE (this%linebuf((m-1)*this%FLEN+1:m*this%FLEN),TRIM(this%fmtstr)) &
                   this%output(m)%p(l)%val(i,j,k), RECSEP
          END DO

          l = SIZE(this%output(m)%p)
          !!!!! ATTENTION this is still 2D code should be different in 3D
          IF ((j.EQ.Mesh%JNUM).AND.((Mesh%JNUM.GT.1).OR.(Mesh%INUM.EQ.i))) THEN
             ! finish the block
             WRITE (this%linebuf((this%COLS-1)*this%FLEN+1:this%linelen),TRIM(this%fmtstr)) &
                  this%output(m)%p(l)%val(i,j,k), BLKSEP
          ELSE
             ! finish the line
             WRITE (this%linebuf((this%COLS-1)*this%FLEN+1:this%linelen),TRIM(this%fmtstr)) &
                  this%output(m)%p(l)%val(i,j,k), LINSEP
          END IF

#ifdef PARALLEL
!!!! ATTENTION parallel output broken
          ! write line buffer to output buffer
!          DO m=1,this%linelen
!             this%outbuf(m,j-1+Mesh%JMIN) = this%linebuf(m:m)
!          END DO
#else
          ! write line buffer to output file
          WRITE (this%unit,FMT=TRIM(this%linefmt),ADVANCE='YES') this%linebuf(1:this%linelen-1)
#endif
         END DO
       END DO
#ifdef PARALLEL
       !*****************************************************************!
       ! This collective call doesn't work for pvfs2 -> bug in ROMIO ?
!!$       CALL MPI_File_write_all(this%handle,this%binout,this%bufsize,&
!!$            this%basictype, this%status, this%error)
       !*****************************************************************!
       ! so we use these two commands instead
!        CALL MPI_File_iwrite(this%handle,this%outbuf,this%bufsize,this%basictype,&
!             request,this%error)
!        CALL MPI_Wait(request,this%status,this%error)
#endif
    END DO
  END SUBROUTINE WriteDataset

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
      CALL this%datafile%CloseFile()
    CALL this%Error(modproc,msg)
  END SUBROUTINE Error

  !> \public Closes the file I/O
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_gnuplot),INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    DEALLOCATE(this%outbuf)
#endif
    DEALLOCATE(this%output)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE fileio_gnuplot_mod
