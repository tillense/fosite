!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fileio_binary.f90                                                 #
!#                                                                           #
!# Copyright (C) 2015                                                        #
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
!----------------------------------------------------------------------------!
!> \author Manuel Jung
!!
!! \brief module for binary file I/O
!!
!! This module implements file I/O, which writes all of the output array in
!! a self describing binary data format including informations about
!! endianness and precision.
!!
!! Specification: [header],[data],[data],[data],..
!! - header: [magic],[endian],[version],[real size],[integer size]
!!              6   +   2    +    1    +     2     +      2       = 13 bytes
!!   These are all ASCII characters except the version, which is single byte
!!   unsigned integer.
!!
!! - data: [key length],[key],[type],[data length],[[dims]],[[data]]
!!              4         *     4          4          *        *       bytes
!!   If [type] indicates a 2D, 3D or 4D array, [data length] includes
!!   8, 12 or 16 bytes extra, for dimensional information. The ASCII [key]
!!   has the in [key length] specified size. The different [type]s are
!!   defined in common/common_dict.f90.
!!
!! To write the binary files, we need one of the following:
!! - a fortran compiler with f2003 Stream IO
!! - a MPI build
!! - on NEC sx9: Set the runtime enviroment variable F_NORCW=5555 (or to
!!   another value)
!!
!! \extends fileio_gnuplot
!! \ingroup fileio
!----------------------------------------------------------------------------!
!#ifdef FORTRAN_STREAMS
#define HAVE_VTK
!#elif defined(NECSXACE) || defined(NECSX9) || defined(NECSX8)
!#define HAVE_VTK
!#define NOSTREAM
!#endif
MODULE fileio_binary_mod
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
#ifdef PARALLEL
#define OFFSET_TYPE INTEGER(KIND=MPI_OFFSET_KIND)
#else
#define OFFSET_TYPE INTEGER
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER, PARAMETER    :: LF = ACHAR(10)        !< line feed
  TYPE, EXTENDS(fileio_base) :: fileio_binary
  CONTAINS
    PROCEDURE :: InitFileio_binary
    PROCEDURE :: WriteHeader
    !PROCEDURE :: ReadHeader
    !PROCEDURE :: WriteTimestamp
    !PROCEDURE :: ReadTimestamp
    PROCEDURE :: WriteDataset
    !PROCEDURE :: ReadDataset
    PROCEDURE :: SetMeshDims
    FINAL :: Finalize
    !PRIVATE
    PROCEDURE :: HasMeshDims
    PROCEDURE :: HasCornerDims
    PROCEDURE :: WriteNode
    PROCEDURE :: WriteKey
    PROCEDURE :: WriteDataAttributes
    PROCEDURE :: GetEndianness
    PROCEDURE :: OpenFile
    PROCEDURE :: CloseFile
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fileio_binary
  !--------------------------------------------------------------------------!

CONTAINS
  !> \public Constructor for the binary file I/O
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs,
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileio_binary(this,Mesh,Physics,Timedisc,Sources,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this          !< \param [in,out] this fileio type
    CLASS(mesh_base),INTENT(IN) :: Mesh          !< \param [in] Mesh mesh type
    CLASS(physics_base),INTENT(IN):: Physics       !< \param [in] Physics Physics type
    CLASS(timedisc_base),INTENT(IN) :: Timedisc        !< \param [in] Timedisc timedisc type
    CLASS(sources_base), POINTER :: Sources !< \param [in] Sources sources type
    TYPE(Dict_TYP),POINTER :: config,IO       !< \param [in] IO Dictionary for I/O
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, DIMENSION(3) :: gsizes,lsizes,indices
    INTEGER           :: lb,extent
#endif
    CHARACTER(LEN=1),DIMENSION(:),POINTER :: mold
    REAL              :: r
    INTEGER           :: i,err
    !------------------------------------------------------------------------!
    this%extension= 'bin'
    CALL this%InitFileio(Mesh,Physics,Timedisc,Sources,config,IO,"binary",this%extension)
    !! basic FileIO initialization
    !IF(PRESENT(fmtname)) THEN
    !  CALL InitFileIO(this,Mesh,Physics,fmt,fmtname,fpath,filename,"bin",&
    !                  stoptime,dtwall,count,fcycles,.FALSE.,unit)
    !ELSE
    !  CALL InitFileIO(this,Mesh,Physics,fmt,"binary",fpath,filename,"bin",&
    !                  stoptime,dtwall,count,fcycles,.FALSE.,unit)
    !END IF

    ! We mark the endianess similar to the TIFF format
    ! See: http://en.wikipedia.org/wiki/Endianness#Endianness_in_files_and_byte_swap
    ! II == Intel Format == Little Endian
    ! MM == Motorola Format == Big Endian
    CALL GetEndianness(this, this%endianness, 'II','MM')
    this%realsize = SIZE(TRANSFER(r, mold))
    this%intsize = SIZE(TRANSFER(i, mold))

    IF(this%realsize.GT.8) &
      CALL this%Error("InitFileIO_binary","Only single and double precision are allowed")

    WRITE(this%realfmt,'(A,I1,A)',IOSTAT=err) '"',this%realsize,'"'

#ifdef PARALLEL
    IF(Mesh%INUM.EQ.Mesh%IMAX) THEN
      this%inum = Mesh%IMAX-Mesh%IMIN+2
    ELSE
      this%inum = Mesh%IMAX-Mesh%IMIN+1
    END IF
    IF(Mesh%JNUM.EQ.Mesh%JMAX) THEN
      this%jnum = Mesh%JMAX-Mesh%JMIN+2
    ELSE
      this%jnum = Mesh%JMAX-Mesh%JMIN+1
    END IF
    IF(Mesh%KNUM.EQ.Mesh%KMAX) THEN
      this%knum = Mesh%KMAX-Mesh%KMIN+2
    ELSE
      this%knum = Mesh%KMAX-Mesh%KMIN+1
    END IF

    ! create the data type for the distributed array of
    ! coordinates and simulation data
    gsizes(1) = Mesh%INUM
    gsizes(2) = Mesh%JNUM
    gsizes(3) = Mesh%KNUM
    lsizes(1) = Mesh%IMAX-Mesh%IMIN+1
    lsizes(2) = Mesh%JMAX-Mesh%JMIN+1
    lsizes(3) = Mesh%KMAX-Mesh%KMIN+1
    indices(1)= Mesh%IMIN-1
    indices(2)= Mesh%JMIN-1
    indices(3)= Mesh%KMIN-1
    this%bufsize = PRODUCT(lsizes)
    CALL MPI_Type_create_subarray(3, gsizes, lsizes, indices, MPI_ORDER_FORTRAN,&
         DEFAULT_MPI_REAL,this%filetype,this%error)
    CALL MPI_Type_commit(this%filetype,this%error)

    ! create the data type for the distributed array of
    ! mesh corner positions
    gsizes(1) = Mesh%INUM+1
    gsizes(2) = Mesh%JNUM+1
    gsizes(3) = Mesh%KNUM+1
    lsizes(1) = this%inum
    lsizes(2) = this%jnum
    lsizes(3) = this%knum
    indices(1)= Mesh%IMIN-1
    indices(2)= Mesh%JMIN-1
    indices(3)= Mesh%KMIN-1
    this%cbufsize = PRODUCT(lsizes)
    CALL MPI_Type_create_subarray(3, gsizes, lsizes, indices, MPI_ORDER_FORTRAN,&
         DEFAULT_MPI_REAL,this%cfiletype,this%error)
    CALL MPI_Type_commit(this%cfiletype,this%error)

    this%first = .TRUE.
#endif
  END SUBROUTINE InitFileio_binary


  !> \public Specific routine to open a file for binary I/O
  !!
  SUBROUTINE OpenFile(this,action)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this    !< \param [in,out] this fileio type
    INTEGER          :: action  !< \param [in] action mode of file access
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset!< \param [in] offset offset for MPI
#endif
    CHARACTER(LEN=3) :: fformat            !< \param [in] fformat file format
    INTEGER :: err
    !------------------------------------------------------------------------!
    INTENT(IN)       :: action
    !------------------------------------------------------------------------!
!#ifdef PARALLEL
!    fformat = 'BIN'
!    SELECT CASE(action)
!    CASE(READONLY)
!#ifdef PARALLEL
!       CALL MPI_File_open(MPI_COMM_WORLD,this%GetFilename(),MPI_MODE_RDONLY, &
!            MPI_INFO_NULL,this%handle,this%error)
!        this%offset = 0
!        CALL MPI_File_seek(this%handle,this%offset,MPI_SEEK_SET,this%error)
!#else
!       OPEN(this%unit,FILE=this%GetFilename(),FORM=fformat,STATUS="OLD", &
!            ACTION="READ",POSITION="REWIND",IOSTAT=this%error)
!       !REWIND (UNIT=this%unit,IOSTAT=this%error)
!#endif
!    CASE(READEND)
!#ifdef PARALLEL
!       CALL MPI_File_open(MPI_COMM_WORLD,this%GetFilename(),IOR(MPI_MODE_RDONLY,&
!            MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)
!       ! opening in append mode doesn't seem to work for pvfs2, hence ...
!       offset = 0
!       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
!       CALL MPI_File_sync(this%handle,this%error)
!#else
!       OPEN(this%unit,FILE=this%GetFilename(),FORM=fformat,STATUS="OLD", &
!            ACTION="READ",POSITION="APPEND",IOSTAT=this%error)
!#endif
!    CASE(REPLACE)
!#ifdef PARALLEL
!       CALL MPI_File_delete(this%GetFilename(),MPI_INFO_NULL,this%error)
!       CALL MPI_File_open(MPI_COMM_WORLD,this%GetFilename(),IOR(MPI_MODE_WRONLY,&
!            MPI_MODE_CREATE),MPI_INFO_NULL,this%handle,this%error)
!#else
!       OPEN(this%unit,FILE=this%GetFilename(),FORM=fformat,STATUS="REPLACE",&
!            ACTION="WRITE",POSITION="REWIND",IOSTAT=this%error)
!#endif
!    CASE(APPEND)
!#ifdef PARALLEL
!       CALL MPI_File_open(MPI_COMM_WORLD,this%GetFilename(),IOR(MPI_MODE_RDWR,&
!            MPI_MODE_APPEND),MPI_INFO_NULL,this%handle,this%error)       
!       ! opening in append mode doesn't seem to work for pvfs2, hence ...
!       offset = 0
!       CALL MPI_File_seek(this%handle,offset,MPI_SEEK_END,this%error)
!       CALL MPI_File_sync(this%handle,this%error)
!#else
!       OPEN(this%unit,FILE=this%GetFilename(),FORM=fformat,STATUS="OLD",&
!            ACTION="READWRITE",POSITION="APPEND",IOSTAT=this%error)
!#endif
!    CASE DEFAULT
!       CALL this%Error("OpenFile","Unknown access mode.")
!    END SELECT
!#else
!#ifdef HAVE_VTK
!    SELECT CASE(action)
!    CASE(READONLY)
!       OPEN(this%unit, FILE=this%GetFilename(), &
!         STATUS     = 'OLD',          &
!#ifndef NOSTREAM 
!         ACCESS     = 'STREAM' ,   &
!#else
!         FORM='UNFORMATTED',&
!#endif
!         action     = 'READ',         &
!         POSITION   = 'REWIND',       &
!         iostat     = this%error_code)
!    CASE(READEND)
!       open(this%unit, FILE=this%GetFilename(), &
!         STATUS     = 'OLD',          &
!#ifndef NOSTREAM
!         ACCESS     = 'STREAM' ,   &
!#else
!         FORM='UNFORMATTED',&
!#endif
!         action     = 'READ',         &
!         POSITION   = 'APPEND',       &
!         iostat     = this%error_code)
!    CASE(REPLACE)
!       open(this%unit, FILE=this%GetFilename(), &
!         STATUS     = 'REPLACE',      &
!#ifndef NOSTREAM
!         ACCESS     = 'STREAM' ,   &
!#else
!         FORM='UNFORMATTED',&
!#endif
!         action     = 'WRITE',        &
!         POSITION   = 'REWIND',       &
!         iostat     = this%error_code)
!#ifdef PARALLEL
!    ! open pvts-file
!      IF (this%GetRank().EQ.0) THEN
!        this%extension='pvts'
!        OPEN(this%unit+100, FILE=this%GetFilename(-1), &
!           STATUS     = 'REPLACE',      &
!#ifndef NOSTREAM
!           ACCESS     = 'STREAM' ,   &
!#else
!           FORM='UNFORMATTED',&
!#endif
!           ACTION     = 'WRITE',        &
!           POSITION   = 'REWIND',       &
!           IOSTAT     = this%error_code)
!        this%extension='vts'
!      END IF
!#endif
!    CASE(APPEND)
!       open(this%unit, FILE=this%GetFilename(), &
!         STATUS     = 'OLD',          &
!#ifndef NOSTREAM
!         ACCESS     = 'STREAM' ,   &
!#else
!         FORM='UNFORMATTED',&
!#endif
!         action     = 'READWRITE',    &
!         POSITION   = 'APPEND',       &
!         iostat     = this%error_code)
!#ifdef PARALLEL
!    ! open pvts-file
!      IF (GetRank(this).EQ.0) THEN
!        this%extension='pvts'
!        OPEN(this%unit+100, FILE=this%GetFilename(-1), &
!           STATUS     = 'OLD',      &
!#ifndef NOSTREAM
!           ACCESS     = 'STREAM' ,   &
!#else
!           FORM='UNFORMATTED',&
!#endif
!           ACTION     = 'READWRITE',        &
!           POSITION   = 'APPEND',       &
!           IOSTAT     = this%error)
!        this%extension='vts'
!      END IF
!#endif
!
!    CASE DEFAULT
!       CALL this%Error("OpenFile","Unknown access mode.")
!    END SELECT
!    IF (this%error_code.NE. 0) CALL this%Error("OpenFile_vtk","Can't open file")
!#endif
!#endif
    OPEN(this%unit,FILE=this%GetFilename(),ACCESS='STREAM',POSITION='APPEND',IOSTAT=err)
    IF (err.NE. 0) CALL this%Error("OpenFile_vtk","Can't open file")
  END SUBROUTINE OpenFile

  !> \public Write the file header
  !! The header is written in ASCII and is 13 Byte long.
  !! First a "magic" identifier is written, than the endianness (II=little,
  !! MM=big), a single byte file format version number, two char realsize,
  !! two char intsize. This would result for example to (\0=hex 0):
  !! "FOSITEII\0 8 4"
  SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT)  :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base),INTENT(IN)  :: Mesh      !< \param [in] Mesh mesh type
    CLASS(physics_base),INTENT(IN) :: Physics   !< \param [in] Physics physics type
    TYPE(Dict_TYP),POINTER :: Header,IO
    !------------------------------------------------------------------------!
    CHARACTER(LEN=6) :: magic = "FOSITE"
    CHARACTER(LEN=1) :: version = ACHAR(0)
    CHARACTER(LEN=4) :: sizes
    CHARACTER(LEN=13):: sheader
    !------------------------------------------------------------------------!
    WRITE(sizes, '(I2,I2)') this%realsize, this%intsize
    sheader = magic // this%endianness(1:2) // version // sizes
    this%offset = 0
#ifndef PARALLEL
    WRITE(this%unit) sheader
#else
    CALL MPI_File_set_view(this%handle,this%offset,MPI_BYTE,&
         MPI_BYTE, 'native', MPI_INFO_NULL, this%error)
    IF(GetRank(this).EQ.0) &
      CALL MPI_File_write(this%handle,sheader,LEN(sheader),MPI_BYTE, &
        this%status,this%error)
#endif
    this%offset = this%offset + LEN(sheader)
  END SUBROUTINE WriteHeader

!  !> \public Reads the header of a file (not yet implemented)
!  !!
!  SUBROUTINE ReadHeader_binary(this,success)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
!    LOGICAL          :: success  !< \param [out] success
!    !------------------------------------------------------------------------!
!    INTENT(OUT)      :: success
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!    success = .FALSE.
!  END SUBROUTINE ReadHeader_binary
!
!  !> \public Writes the timestep (not yet implemented)
!  !!
!  SUBROUTINE WriteTimestamp_binary(this,time)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
!    REAL             :: time     !< \param [in] time
!    !------------------------------------------------------------------------!
!    INTENT(IN)       :: time
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!  END SUBROUTINE WriteTimestamp_binary
!
!  !> \public Writes the timestep (not yet implemented)
!  !!
!  SUBROUTINE ReadTimestamp_binary(this,time)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP) :: this     !< \param [in,out] this fileio type
!    REAL             :: time     !< \param [out] time
!    !------------------------------------------------------------------------!
!    INTENT(OUT)      :: time
!    INTENT(INOUT)    :: this
!    !------------------------------------------------------------------------!
!  END SUBROUTINE ReadTimestamp_binary


  !> Writes key structure
  !! This subroutine writes the the key, data type and data sizes. 
  !! It is defined as following (suppose 4B Integer)
  !! | 4B length of key | *B key | 4B data type | 4B data size in bytes |
  !! If the data type is a 2D,3D or 4D array, 2, 3 or 4 (4 Byte) integers
  !! are appended with the shape information. There storage is included in
  !! the data size field. Therefore without knowing the data types, one can
  !! jump over the data to the next key structure.
  SUBROUTINE WriteKey(this,key,type,bytes,dims)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this         !< \param [in,out] this fileio type
    CHARACTER(LEN=*)  :: key
    INTEGER           :: type,bytes
    INTEGER,DIMENSION(5),OPTIONAL :: dims
    INTEGER           :: bufsize
    CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE :: buf
    !------------------------------------------------------------------------!
    INTEGER           :: keylen,l,b,o
    !------------------------------------------------------------------------!
    INTENT(IN)        :: key,type,bytes
    !------------------------------------------------------------------------!
    keylen = LEN_TRIM(key)
    IF(PRESENT(dims)) THEN
      l = 3
      IF(dims(4).GT.1) l = 4
      IF(dims(5).GT.1) l = 5
    ELSE
      l = 0
    END IF
    b = bytes + l * this%intsize
    bufsize = (3+l) * this%intsize + keylen
    ALLOCATE(buf(bufsize))
    o = 1
    CALL Append(buf,o,TRANSFER(keylen,buf))
    CALL Append(buf,o,TRANSFER(TRIM(key),buf))
    CALL Append(buf,o,TRANSFER(type,buf))
    CALL Append(buf,o,TRANSFER(b,buf))
    IF(l.GT.0) THEN
      CALL Append(buf,o,TRANSFER(dims(1:l),buf))
    END IF
#ifndef PARALLEL
     WRITE(this%unit) buf
#else
    CALL MPI_File_set_view(this%handle,this%offset,MPI_BYTE,&
           MPI_BYTE, 'native', MPI_INFO_NULL, this%error)
    IF(GetRank(this).EQ.0) &
      CALL MPI_File_write(this%handle,buf,bufsize,MPI_BYTE, &
        this%status,this%error)
#endif
    DEALLOCATE(buf)
    this%offset = this%offset + bufsize

    CONTAINS
      SUBROUTINE Append(buffer,i,d)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        CHARACTER(LEN=1),DIMENSION(:) :: buffer,d
        INTEGER :: i,s
        !--------------------------------------------------------------------!
        s = SIZE(d)
        buffer(i:i+s-1) = d
        o = o + s
      END SUBROUTINE
  END SUBROUTINE WriteKey


  !> Writes data attributes to a file
  !!
  RECURSIVE SUBROUTINE WriteDataAttributes(this,Mesh,config,path)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this         !< \param [in,out] this fileio type
    CLASS(mesh_base),INTENT(IN) :: Mesh         !< \param [in] mesh mesh type
    TYPE(Dict_TYP),POINTER :: config  !< \param [in] config dict of configuration
    CHARACTER(LEN=*),OPTIONAL &
                      :: path    !< \param [in] path
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: str, key
    TYPE(dict_TYP),POINTER   :: node
    REAL,DIMENSION(2,Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) :: buf
    CHARACTER(LEN=1),DIMENSION(:),POINTER :: val
    !------------------------------------------------------------------------!
    INTENT(IN)        :: path
    !------------------------------------------------------------------------!
    IF(PRESENT(path)) THEN
        str = path
    ELSE
        str = ""
    ENDIF
    node => config
    DO WHILE(ASSOCIATED(node))
      key = TRIM(str)//"/"//TRIM(GetKey(node))

      CALL this%WriteNode(Mesh,key,node)

      IF(HasChild(node)) THEN
        CALL this%WriteDataAttributes(Mesh, GetChild(node), key)
      END IF
      node=>GetNext(node)
    END DO
  END SUBROUTINE WriteDataAttributes

  SUBROUTINE WriteNode(this,Mesh,key,node)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this         !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN) ::  Mesh         !< \param [in] mesh mesh type
    CHARACTER(LEN=MAX_CHAR_LEN) :: key
    TYPE(Dict_TYP),POINTER :: node    !< \param [in] data node
    !------------------------------------------------------------------------!
    TYPE(real_t) :: ptr0
    TYPE(int_t) :: ptrint
    REAL,DIMENSION(:,:),POINTER :: ptr2
    REAL,DIMENSION(:,:,:),POINTER :: ptr3
    REAL,DIMENSION(:,:,:,:),POINTER :: ptr4
    REAL,DIMENSION(:,:,:,:,:),POINTER :: ptr5
    INTEGER,DIMENSION(5) :: dims,dims2
    INTEGER                :: bytes,k,l,pos,imax,jmax
    CHARACTER(LEN=1),DIMENSION(:),POINTER :: val
    REAL,DIMENSION(:,:),ALLOCATABLE :: buf
#ifdef PARALLEL
    OFFSET_TYPE       :: omax,omin
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)        :: key
    !------------------------------------------------------------------------!
    dims(:) = 1
    NULLIFY(ptr2,ptr3,ptr4,ptr5)

    SELECT CASE(GetDataType(node))
    CASE(DICT_REAL_THREED)
      CALL GetAttr(node,key,ptr3)
      dims(1:3) = SHAPE(ptr3)
    CASE(DICT_REAL_FOURD)
      CALL GetAttr(node,key,ptr4)
      dims(1:4) = SHAPE(ptr4)
    CASE(DICT_REAL_FIVED)
      CALL GetAttr(node,key,ptr5)
      dims(1:5) = SHAPE(ptr5)
    END SELECT

    IF(PRODUCT(dims).GT.1) THEN
      CALL this%SetMeshDims(Mesh,dims)
      bytes = PRODUCT(dims(1:3)) * this%realsize
      CALL this%WriteKey(key,GetDataType(node),&
        PRODUCT(dims)*this%realsize,dims)
      DO l=1,dims(5)
        DO k=1,dims(4)
          IF(ASSOCIATED(ptr4)) THEN
            ptr3 => ptr4(:,:,:,k)
          ELSE IF(ASSOCIATED(ptr5)) THEN
            ptr3 => ptr5(:,:,:,k,l)
          END IF
#ifdef PARALLEL
          IF(this%HasMeshDims(Mesh,SHAPE(ptr3))) THEN
            CALL MPI_File_set_view(this%handle,this%offset,DEFAULT_MPI_REAL,&
                   this%filetype, 'native', MPI_INFO_NULL, this%error)
            CALL MPI_File_write_all(this%handle,ptr3,this%bufsize,&
                   DEFAULT_MPI_REAL,this%status,this%error)

          ELSE IF(this%HasCornerDims(Mesh,SHAPE(ptr3))) THEN
            CALL MPI_File_set_view(this%handle,this%offset,DEFAULT_MPI_REAL,&
                   this%cfiletype, 'native', MPI_INFO_NULL, this%error)
            CALL MPI_File_write_all(this%handle,ptr3(1:this%inum,1:this%jnum,1:this%knum),&
                   this%cbufsize,DEFAULT_MPI_REAL,this%status,this%error)
          ELSE

            CALL MPI_File_set_view(this%handle,this%offset,MPI_BYTE,&
                 MPI_BYTE, 'native', MPI_INFO_NULL, this%error)
            IF(this%GetRank().EQ.0) THEN
              CALL MPI_File_write(this%handle,ptr3,bytes,MPI_BYTE, &
                   this%status,this%error)
            END IF
          END IF
#else
          WRITE(this%unit) ptr3
#endif
          this%offset = this%offset + bytes
        END DO
      END DO
    ELSE
      SELECT CASE(GetDatatype(node))
      CASE(DICT_REAL_P)
        CALL GetAttr(node,key,ptr0)
        bytes = this%realsize
        CALL this%WriteKey(key,DICT_REAL,bytes)
#ifndef PARALLEL
          WRITE(this%unit) ptr0%p
#else
          CALL MPI_File_set_view(this%handle,this%offset,MPI_BYTE,&
               MPI_BYTE, 'native', MPI_INFO_NULL, this%error)
          IF(this%GetRank().EQ.0) THEN
            CALL MPI_File_write(this%handle,ptr0%p,bytes,MPI_BYTE, &
                 this%status,this%error)
          END IF
#endif
      CASE(DICT_INT_P)
        CALL GetAttr(node,key,ptrint)
        bytes = this%intsize
        CALL this%WriteKey(key,DICT_INT,bytes)
#ifndef PARALLEL
          WRITE(this%unit) ptrint%p
#else
          CALL MPI_File_set_view(this%handle,this%offset,MPI_BYTE,&
               MPI_BYTE, 'native', MPI_INFO_NULL, this%error)
          IF(this%GetRank().EQ.0) THEN
            CALL MPI_File_write(this%handle,ptrint%p,bytes,MPI_BYTE, &
                 this%status,this%error)
          END IF
#endif
      CASE DEFAULT
        val => GetData(node)
        IF(ASSOCIATED(val)) THEN
          bytes = SIZE(val)
        ELSE
          bytes = 0
        END IF
        CALL this%WriteKey(key,GetDataType(node),bytes)
        IF(bytes.GT.0) THEN
#ifndef PARALLEL
          WRITE(this%unit) val
#else
          CALL MPI_File_set_view(this%handle,this%offset,MPI_BYTE,&
               MPI_BYTE, 'native', MPI_INFO_NULL, this%error)
          IF(this%GetRank().EQ.0) THEN
            CALL MPI_File_write(this%handle,val,bytes,MPI_BYTE, &
                 this%status,this%error)
          END IF
#endif
        END IF
      END SELECT
      this%offset = this%offset + bytes
    END IF
#ifdef PARALLEL
    IF(this%first) THEN
      ! If MPI is used and this is the first output of the run,
      ! it is checked, if the offsets of the different nodes are still
      ! in sync. They can get easily out of sync is an output array is not
      ! a mesh or corner array, but still has a different size on some nodes.
      ! This easily happens, e.g. if one forgets to specify the subarray limits.
      ! output_array(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
      ! The array would without them be interpreted as generic 2D array in
      ! this case and therefore the calculated size could be different on
      ! different nodes.
      CALL MPI_Allreduce(this%offset,omax,1,MPI_OFFSET,MPI_MAX,&
        Mesh%comm_cart,this%error)
      CALL MPI_Allreduce(this%offset,omin,1,MPI_OFFSET,MPI_MIN,&
        Mesh%comm_cart,this%error)
      IF((this%offset.NE.omax).OR.(this%offset.NE.omin)) &
        CALL this%Error("WriteNode_binary",&
          "The offsets on different nodes are not in sync anymore." // LF &
        //"The last key was '" // TRIM(key) // "'.")
    END IF
#endif
  END SUBROUTINE WriteNode

  FUNCTION HasMeshDims(this,Mesh,dims) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT)  :: this      !< \param [in] this type
    CLASS(mesh_base),INTENT(IN)  :: Mesh      !< \param [in] mesh mesh type
    INTEGER,DIMENSION(:) :: dims
    LOGICAL           :: res
    !------------------------------------------------------------------------!
    INTENT(IN)        :: dims
    !------------------------------------------------------------------------!
    IF(SIZE(dims).GE.3) THEN
      res = (dims(1).EQ.(Mesh%IMAX-Mesh%IMIN+1)) &
            .AND.(dims(2).EQ.(Mesh%JMAX-Mesh%JMIN+1)) &
            .AND.(dims(3).EQ.(Mesh%KMAX-Mesh%KMIN+1))
    ELSE
      res = .FALSE.
    END IF
  END FUNCTION HasMeshDims


  FUNCTION HasCornerDims(this,Mesh,dims) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT)  :: this      !< \param [in] this type
    CLASS(mesh_base),INTENT(IN)    :: Mesh      !< \param [in] mesh mesh type
    INTEGER,DIMENSION(:) :: dims
    LOGICAL           :: res
    !------------------------------------------------------------------------!
    INTENT(IN)        :: dims
    !------------------------------------------------------------------------!
    IF(SIZE(dims).GE.3) THEN
      res = (dims(1).EQ.(Mesh%IMAX-Mesh%IMIN+2)) &
            .AND.(dims(2).EQ.(Mesh%JMAX-Mesh%JMIN+2)) &
            .AND.(dims(3).EQ.(Mesh%KMAX-Mesh%KMIN+2))
    ELSE
      res = .FALSE.
    END IF
  END FUNCTION HasCornerDims


  SUBROUTINE SetMeshDims(this,Mesh,dims)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN) :: Mesh
    INTEGER,DIMENSION(:)   :: dims
    !------------------------------------------------------------------------!
    INTENT(INOUT)          :: dims
    !------------------------------------------------------------------------!
    IF(this%HasMeshDims(Mesh,dims)) THEN
      dims(1) = Mesh%INUM
      dims(2) = Mesh%JNUM
      dims(3) = Mesh%KNUM
    ELSE IF(this%HasCornerDims(Mesh,dims)) THEN
      dims(1) = Mesh%INUM+1
      dims(2) = Mesh%JNUM+1
      dims(3) = Mesh%KNUM+1
    END IF
  END SUBROUTINE SetMeshDims


  !> \public Writes all desired data arrays to a file 
  !!
  SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base),INTENT(IN) :: Mesh      !< \param [in] mesh mesh type
    CLASS(physics_base),INTENT(IN) :: Physics   !< \param [in] physics physics type
    CLASS(fluxes_base),INTENT(IN) :: Fluxes    !< \param [in] fluxes fluxes type
    CLASS(timedisc_base),INTENT(IN) :: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP),POINTER :: Header,IO   !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: meshIO
    !------------------------------------------------------------------------!
    ! write data
    CALL this%OpenFile(APPEND)
    CALL this%WriteHeader(Mesh,Physics,Header,IO)
    CALL this%WriteDataAttributes(Mesh,IO)
    CALL this%CloseFile()
    CALL this%IncTime()

#ifdef PARALLEL
    this%first = .FALSE.
#endif

  END SUBROUTINE WriteDataset

!  !> \public Reads the data arrays from file (not yet implemented)
!  !!
!  SUBROUTINE ReadDataset_binary(this,Mesh,Physics,Timedisc)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(FileIO_TYP)  :: this      !< \param [in,out] this fileio type
!    TYPE(Mesh_TYP)    :: Mesh      !< \param [in] mesh mesh type
!    TYPE(Physics_TYP) :: Physics   !< \param [in] physics physics type
!    TYPE(Timedisc_TYP):: Timedisc  !< \param [in,out] timedisc timedisc type
!    !------------------------------------------------------------------------!
!    INTENT(IN)        :: Mesh,Physics
!    INTENT(INOUT)     :: this,Timedisc
!    !------------------------------------------------------------------------!
!  END SUBROUTINE ReadDataset_binary

  !> \public Determines the endianness of the system
  !!
  !! Determines the the endianess of the system (big or little endian)
  SUBROUTINE GetEndianness(this, res, littlestr, bigstr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this         !< \param [in,out] this fileio type
    CHARACTER(LEN=*)  :: res          !< \param [out] res result string
    CHARACTER(LEN=*)  :: littlestr    !< \param [in] littlestr little endian str
    CHARACTER(LEN=*)  :: bigstr       !< \param [in] bigstr big endian str
    !------------------------------------------------------------------------!
    INTEGER           :: k,err,iTIPO
    CHARACTER, POINTER:: cTIPO(:)
    !------------------------------------------------------------------------!
    INTENT(OUT)       :: res
    !------------------------------------------------------------------------!

    !endianness
    k = BIT_SIZE(iTIPO)/8
    ALLOCATE(cTIPO(k),STAT = err)
       IF (err.NE.0) &
         CALL this%Error("GetEndianness_binary", "Unable to allocate memory.")
    cTIPO(1)='A'
    !cTIPO(2:k-1) = That's of no importance.
    cTIPO(k)='B'

    iTIPO = transfer(cTIPO, iTIPO)
    DEALLOCATE(cTIPO)
    !Test of 'B'=b'01000010' ('A'=b'01000001')
    IF (BTEST(iTIPO,1)) THEN ! big endian
       write(res,'(A)',IOSTAT=err)bigstr
    ELSE ! little endian
       write(res,'(A)',IOSTAT=err)littlestr
    END IF
  END SUBROUTINE GetEndianness

  !> \public routine to close a file
  !!
  SUBROUTINE CloseFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary),INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    INTEGER :: err
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_File_close(this%handle,this%error_code)
#else
    CLOSE(this%unit,IOSTAT=err)
#endif
  END SUBROUTINE CloseFile

  !> \public Closes the file I/O
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fileio_binary),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_Type_free(this%cfiletype,this%error_code)
    CALL MPI_Type_free(this%filetype,this%error_code)
#endif
  END SUBROUTINE Finalize

END MODULE fileio_binary_mod
