!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_binary.f90                                                 #
!#                                                                           #
!# Copyright (C) 2015-2024                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> \author Manuel Jung
!! \author Jannes Klee
!! \author Tobias Illenseer
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
!!   If [type] indicates a 2D, 3D, 4D or 5D array, [data length] includes
!!   8, 12, 16 or 24 bytes extra, for dimensional information. The ASCII [key]
!!   has the in [key length] specified size. The different [type]s are
!!   defined in common/common_dict.f90.
!!
!! To write the binary files, we need one of the following:
!! - a fortran compiler with f2003 Stream IO
!! - a MPI build
!! - on NEC sx9: Set the runtime enviroment variable F_NORCW=5555 (or to
!!   another value)
!!
!! \extends fileio_base
!! \ingroup fileio
!----------------------------------------------------------------------------!
#define HAVE_VTK
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
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER, PARAMETER       :: LF = ACHAR(10)        !< line feed
  !--------------------------------------------------------------------------!
  !> FileIO binary class
  TYPE, EXTENDS(fileio_base) :: fileio_binary
    !> \name Variables
    CHARACTER(LEN=12)      :: realfmt     !< real format string
    CHARACTER(LEN=14)      :: endianness  !< endianness string
    INTEGER                :: realsize    !< byte size of real numbers
    INTEGER                :: intsize     !< byte size of integer numbers
#ifdef PARALLEL
    LOGICAL                :: first       !< true if this is the first output
    INTEGER                :: cbufsize, & !< size of corner output
                              clsizes(3)  !< local array sizes for corner output
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset !< skip header bytes
#else
    INTEGER                :: offset
#endif
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFileio_deferred => InitFileio_binary
    PROCEDURE :: WriteHeader
    !PROCEDURE :: ReadHeader
    !PROCEDURE :: WriteTimestamp
    !PROCEDURE :: ReadTimestamp
    PROCEDURE :: WriteDataset_deferred => WriteDataset_binary
    !PROCEDURE :: ReadDataset
    FINAL :: Finalize
    !PRIVATE
    PROCEDURE :: HasMeshDims
    PROCEDURE :: HasCornerDims
    PROCEDURE :: SetOutputDims
    PROCEDURE :: WriteNode
    PROCEDURE :: WriteKey
    PROCEDURE :: WriteDataAttributes
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
    CLASS(fileio_binary), INTENT(INOUT) :: this       !< \param [in,out] this fileio type
    CLASS(mesh_base),     INTENT(IN)    :: Mesh       !< \param [in] Mesh mesh type
    CLASS(physics_base),  INTENT(IN)    :: Physics    !< \param [in] Physics Physics type
    CLASS(timedisc_base), INTENT(IN)    :: Timedisc   !< \param [in] Timedisc timedisc type
    CLASS(sources_base),  INTENT(IN), POINTER :: Sources    !< \param [in] Sources sources type
    TYPE(Dict_TYP),       INTENT(IN), POINTER :: config,IO  !< \param [in] IO Dictionary for I/O
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, DIMENSION(3)                 :: gsizes,lsizes,indices
#endif
    CHARACTER(LEN=1),DIMENSION(:),POINTER :: mold
    REAL                                  :: r
    INTEGER                               :: i,err
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
      CALL this%InitFileio(Mesh,Physics,Timedisc,Sources,config,IO,"binary","bin",textfile=.FALSE.)

    ! We mark the endianess similar to the TIFF format
    ! See: http://en.wikipedia.org/wiki/Endianness#Endianness_in_files_and_byte_swap
    ! II == Intel Format == Little Endian
    ! MM == Motorola Format == Big Endian
    CALL this%GetEndianness(this%endianness, 'II','MM')
    this%realsize = SIZE(TRANSFER(r, mold))
    this%intsize = SIZE(TRANSFER(i, mold))

    SELECT CASE(this%realsize)
    CASE(4,8)
      WRITE(this%realfmt,'(A,I1,A)') '"',this%realsize,'"'
    CASE DEFAULT
      CALL this%Error("fileio_binary::InitFileIO_binary","Only single and double precision are allowed")
    END SELECT

    ! set local array dimensions (frequently used)
    this%INUM = Mesh%IMAX-Mesh%IMIN+1
    this%JNUM = Mesh%JMAX-Mesh%JMIN+1
    this%KNUM = Mesh%KMAX-Mesh%KMIN+1

#ifdef PARALLEL
    ! create the data type for the distributed array of
    ! coordinates and simulation data
    gsizes(1:3)  = (/ Mesh%INUM, Mesh%JNUM, Mesh%KNUM /)
    lsizes(1:3)  = (/ this%INUM, this%JNUM, this%KNUM /)
    indices(1:3) = (/ Mesh%IMIN-1, Mesh%JMIN-1, Mesh%KMIN-1 /)
    this%bufsize = PRODUCT(lsizes)
    this%err = 0
    SELECT TYPE(df=>this%datafile)
    CLASS IS(filehandle_mpi)
      CALL MPI_Type_create_subarray(3, gsizes, lsizes, indices, MPI_ORDER_FORTRAN,&
           DEFAULT_MPI_REAL,df%filetype,this%err)
      IF (this%err.EQ.0) CALL MPI_Type_commit(df%filetype,this%err)
    END SELECT

    ! create the data type for the distributed array of
    ! mesh corner positions
    this%clsizes(:) = lsizes(:)
    WHERE ((/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/).GT.1)
      ! add +1 to the global size if dimension has real extend,i.e. not one cell
      gsizes(:) = gsizes(:)+1
      WHERE ((/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/).EQ.(/Mesh%IMAX,Mesh%JMAX,Mesh%KMAX/))
        ! add +1 to the local size if process contains the outermost domain
        this%clsizes(:) = lsizes(:)+1
      END WHERE
    END WHERE
    this%cbufsize = PRODUCT(this%clsizes)
    SELECT TYPE(df=>this%datafile)
    CLASS IS(filehandle_mpi)
      IF (this%err.EQ.0) CALL MPI_Type_create_subarray(3, gsizes, this%clsizes, indices, MPI_ORDER_FORTRAN,&
          DEFAULT_MPI_REAL,df%cfiletype,this%err)
      IF (this%err.EQ.0) CALL MPI_Type_commit(df%cfiletype,this%err)
    END SELECT

    IF (this%err.NE.0) &
       CALL this%Error("fileio_binary::InitFileIO_binary","creating MPI data types failed")
    this%first = .TRUE.
#endif
  END SUBROUTINE InitFileio_binary

  !> \public Write the file header
  !! The header is written in ASCII and is 13 Byte long.
  !! First a "magic" identifier is written, than the endianness (II=little,
  !! MM=big), a single byte file format version number, two char realsize,
  !! two char intsize. This would result for example to (\0=hex 0):
  !! "FOSITEII\0 8 4"
  SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)        :: Mesh      !< \param [in] Mesh mesh type
    CLASS(physics_base), INTENT(IN)     :: Physics   !< \param [in] Physics physics type
    TYPE(Dict_TYP), POINTER             :: Header,IO
    !------------------------------------------------------------------------!
    CHARACTER(LEN=6)  :: magic = "FOSITE"
    CHARACTER(LEN=1)  :: version = ACHAR(0)
    CHARACTER(LEN=4)  :: sizes
    CHARACTER(LEN=13) :: sheader
    !------------------------------------------------------------------------!
    WRITE(sizes, '(I2,I2)') this%realsize, this%intsize
    sheader = magic // this%endianness(1:2) // version // sizes
    this%offset = 0
    SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
    CLASS IS(filehandle_fortran)
      WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) sheader
#else
    CLASS IS(filehandle_mpi)
      CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,MPI_BYTE,&
               MPI_BYTE, 'native', MPI_INFO_NULL, this%err)
      IF (this%GetRank().EQ.0) &
        CALL MPI_File_write(df%GetUnitNumber(),sheader,LEN(sheader),MPI_BYTE, &
                 df%status,this%err)
#endif
    CLASS DEFAULT
      CALL this%Error("fileio_binary::WriteHeader","unknown file handle")
    END SELECT
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
    CLASS(fileio_binary), INTENT(INOUT)         :: this !< \param [in,out] this fileio type
    CHARACTER(LEN=*)                            :: key
    INTEGER                                     :: type,bytes
    INTEGER, DIMENSION(5), OPTIONAL             :: dims
    INTEGER                                     :: bufsize
    CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: buf
    !------------------------------------------------------------------------!
    INTEGER                                     :: keylen,l,b,o
    !------------------------------------------------------------------------!
    INTENT(IN)                                  :: key,type,bytes
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

    SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
    CLASS IS(filehandle_fortran)
      WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) buf
#else
    CLASS IS(filehandle_mpi)
      CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,MPI_BYTE,&
            MPI_BYTE, 'native', MPI_INFO_NULL, this%err)
      IF (this%GetRank().EQ.0) &
        CALL MPI_File_write(df%GetUnitNumber(),buf,bufsize,MPI_BYTE, &
          df%status,this%err)
#endif
    CLASS DEFAULT
      CALL this%Error("fileio_binary::WriteKey","unknown file handle")
    END SELECT

    DEALLOCATE(buf)
    this%offset = this%offset + bufsize

    CONTAINS
      SUBROUTINE Append(buffer,i,d)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        CHARACTER(LEN=1), DIMENSION(:) :: buffer,d
        INTEGER                        :: i,s
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
    CLASS(fileio_binary), INTENT(INOUT)     :: this   !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)            :: Mesh   !< \param [in] mesh mesh type
    TYPE(Dict_TYP), POINTER                 :: config !< \param [in] config dict of configuration
    CHARACTER(LEN=*), OPTIONAL              :: path   !< \param [in] path
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN)             :: str, key
    TYPE(dict_TYP), POINTER                 :: node
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: path
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
    CLASS(fileio_binary), INTENT(INOUT) :: this !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)        :: Mesh !< \param [in] mesh mesh type
    CHARACTER(LEN=MAX_CHAR_LEN)         :: key
    TYPE(Dict_TYP), POINTER             :: node !< \param [in] data node
    !------------------------------------------------------------------------!
    TYPE(real_t)                            :: ptr0
    TYPE(int_t)                             :: ptrint
    REAL, DIMENSION(:,:), POINTER, CONTIGUOUS :: ptr2
    REAL, DIMENSION(:,:,:), POINTER         :: ptr3
    REAL, DIMENSION(:,:,:,:), POINTER       :: ptr4
    REAL, DIMENSION(:,:,:,:,:), POINTER     :: ptr5
    INTEGER, DIMENSION(5)                   :: dims
    INTEGER                                 :: bytes,k,l
    CHARACTER(LEN=1), DIMENSION(:), POINTER :: val
#ifdef PARALLEL
    INTEGER(KIND=MPI_OFFSET_KIND)           :: omax,omin
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: key
    !------------------------------------------------------------------------!
    dims(:) = 1
    NULLIFY(ptr2,ptr3,ptr4,ptr5)

    SELECT CASE(GetDataType(node))
    CASE(DICT_REAL_TWOD)
      CALL GetAttr(node,key,ptr2)
      dims(1:2) = SHAPE(ptr2)
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
      CALL this%SetOutputDims(Mesh,dims)
      bytes = PRODUCT(dims(1:3)) * this%realsize
      CALL this%WriteKey(key,GetDataType(node),&
        PRODUCT(dims)*this%realsize,dims)
      DO l=1,dims(5)
        DO k=1,dims(4)
          IF(ASSOCIATED(ptr4)) THEN
            ptr3 => ptr4(:,:,:,k)
          ELSE IF(ASSOCIATED(ptr5)) THEN
            ptr3 => ptr5(:,:,:,k,l)
          ELSE IF(ASSOCIATED(ptr2)) THEN
            ptr3(1:dims(1),1:dims(2),1:1) => ptr2
          END IF
          SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
          CLASS IS(filehandle_fortran)
            WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) ptr3
#else
          CLASS IS(filehandle_mpi)
            IF(this%HasMeshDims(Mesh,SHAPE(ptr3))) THEN
              CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,DEFAULT_MPI_REAL,&
                    df%filetype, 'native', MPI_INFO_NULL, this%err)
              CALL MPI_File_write_all(df%GetUnitNumber(),ptr3,this%bufsize,&
                    DEFAULT_MPI_REAL,df%status,this%err)

            ELSE IF(this%HasCornerDims(Mesh,SHAPE(ptr3))) THEN
              CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,DEFAULT_MPI_REAL,&
                    df%cfiletype, 'native', MPI_INFO_NULL, this%err)
              CALL MPI_File_write_all(df%GetUnitNumber(),ptr3(1:this%clsizes(1),1:this%clsizes(2),1:this%clsizes(3)),&
                    this%cbufsize,DEFAULT_MPI_REAL,df%status,this%err)
            ELSE
              CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,MPI_BYTE,&
                  MPI_BYTE, 'native', MPI_INFO_NULL, this%err)
              IF(this%GetRank().EQ.0) THEN
                CALL MPI_File_write(df%GetUnitNumber(),ptr3,bytes,MPI_BYTE, &
                    df%status,this%err)
              END IF
            END IF
#endif
          CLASS DEFAULT
            CALL this%Error("fileio_binary::WriteNode","unknown file handle")
          END SELECT
          this%offset = this%offset + bytes
        END DO
      END DO
    ELSE
      SELECT CASE(GetDatatype(node))
      CASE(DICT_REAL_P)
        CALL GetAttr(node,key,ptr0)
        bytes = this%realsize
        CALL this%WriteKey(key,DICT_REAL,bytes)
          SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
          CLASS IS(filehandle_fortran)
            WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) ptr0%p
#else
          CLASS IS(filehandle_mpi)
            CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,MPI_BYTE,&
                MPI_BYTE, 'native', MPI_INFO_NULL, this%err)
            IF(this%GetRank().EQ.0) &
              CALL MPI_File_write(df%GetUnitNumber(),ptr0%p,bytes,MPI_BYTE, &
                  df%status,this%err)
#endif
          CLASS DEFAULT
            CALL this%Error("fileio_binary::WriteNode","unknown file handle")
          END SELECT
      CASE(DICT_INT_P)
        CALL GetAttr(node,key,ptrint)
        bytes = this%intsize
        CALL this%WriteKey(key,DICT_INT,bytes)
          SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
          CLASS IS(filehandle_fortran)
            WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) ptrint%p
#else
          CLASS IS(filehandle_mpi)
            CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,MPI_BYTE,&
                MPI_BYTE, 'native', MPI_INFO_NULL, this%err)
            IF(this%GetRank().EQ.0) &
              CALL MPI_File_write(df%GetUnitNumber(),ptrint%p,bytes,MPI_BYTE, &
                  df%status,this%err)
#endif
          CLASS DEFAULT
            CALL this%Error("fileio_binary::WriteNode","unknown file handle")
          END SELECT
      CASE DEFAULT
        val => GetData(node)
        IF(ASSOCIATED(val)) THEN
          bytes = SIZE(val)
        ELSE
          bytes = 0
        END IF
        CALL this%WriteKey(key,GetDataType(node),bytes)
        IF(bytes.GT.0) THEN
          SELECT TYPE(df=>this%datafile)
#ifndef PARALLEL
          CLASS IS(filehandle_fortran)
            WRITE (UNIT=df%GetUnitNumber(),IOSTAT=this%err) val
#else
          CLASS IS(filehandle_mpi)
            CALL MPI_File_set_view(df%GetUnitNumber(),this%offset,MPI_BYTE,&
                MPI_BYTE, 'native', MPI_INFO_NULL, this%err)
            IF(this%GetRank().EQ.0) THEN
              CALL MPI_File_write(df%GetUnitNumber(),val,bytes,MPI_BYTE, &
                  df%status,this%err)
            END IF
#endif
          CLASS DEFAULT
            CALL this%Error("fileio_binary::WriteNode","unknown file handle")
          END SELECT
        END IF
      END SELECT
      this%offset = this%offset + bytes
    END IF
#ifdef PARALLEL
    IF(this%first) THEN
      ! If MPI is used and this is the first output of the run,
      ! it is checked, if the offsets of the different nodes are still
      ! in sync. They can get easily out of sync if an output array is not
      ! a mesh or corner array, but still has a different size on some nodes.
      ! This easily happens, e.g. if one forgets to specify the subarray limits.
      ! output_array(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
      ! The array would without them be interpreted as generic 2D array in
      ! this case and therefore the calculated size could be different on
      ! different nodes.
      CALL MPI_Allreduce(this%offset,omax,1,MPI_OFFSET,MPI_MAX,&
        Mesh%comm_cart,this%err)
      CALL MPI_Allreduce(this%offset,omin,1,MPI_OFFSET,MPI_MIN,&
        Mesh%comm_cart,this%err)
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
    CLASS(fileio_binary), INTENT(INOUT) :: this !< \param [in] this type
    CLASS(mesh_base), INTENT(IN)        :: Mesh !< \param [in] mesh mesh type
    INTEGER, DIMENSION(:), INTENT(IN)   :: dims !< \param [in] dims array dimensions
    LOGICAL                             :: res
    !------------------------------------------------------------------------!
    IF(SIZE(dims).GE.3) THEN
      res = ALL(dims(1:3).EQ.(/Mesh%IMAX-Mesh%IMIN+1, &
                               Mesh%JMAX-Mesh%JMIN+1, &
                               Mesh%KMAX-Mesh%KMIN+1/))
    ELSE
      res = .FALSE.
    END IF
  END FUNCTION HasMeshDims


  FUNCTION HasCornerDims(this,Mesh,dims) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary), INTENT(INOUT) :: this !< \param [in] this type
    CLASS(mesh_base), INTENT(IN)        :: Mesh !< \param [in] mesh mesh type
    INTEGER, DIMENSION(:)               :: dims
    LOGICAL                             :: res
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: dims
    !------------------------------------------------------------------------!
    IF(SIZE(dims).GE.3) THEN
      ! remark: xP1 = 1 if xNUM > 1 and 0 if xNUM = 1 for x=I,J,K
      res = ALL(dims(1:3).EQ.(/Mesh%IMAX-Mesh%IMIN+Mesh%IP1+1, &
                               Mesh%JMAX-Mesh%JMIN+Mesh%JP1+1, &
                               Mesh%KMAX-Mesh%KMIN+Mesh%KP1+1/))
    ELSE
      res = .FALSE.
    END IF
  END FUNCTION HasCornerDims


  SUBROUTINE SetOutputDims(this,Mesh,dims)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary), INTENT(INOUT) :: this
    CLASS(mesh_base), INTENT(IN)        :: Mesh
    INTEGER, DIMENSION(:), INTENT(INOUT):: dims
    !------------------------------------------------------------------------!
    IF(this%HasMeshDims(Mesh,dims)) THEN
      dims(1:3) = (/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/)
    ELSE IF(this%HasCornerDims(Mesh,dims)) THEN
      dims(1:3) = (/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/)
      WHERE ((/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/).GT.1)
        ! add +1 to the global size to write one additional data point
        ! in each dimension which has extend larger than 1
        dims(:) = dims(:) + 1
      END WHERE
    END IF
  END SUBROUTINE SetOutputDims


  !> \public Writes all desired data arrays to a file
  !!
  SUBROUTINE WriteDataset_binary(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_binary), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CLASS(mesh_base), INTENT(IN)        :: Mesh      !< \param [in] mesh mesh type
    CLASS(physics_base), INTENT(INOUT)  :: Physics   !< \param [in] physics physics type
    CLASS(fluxes_base), INTENT(IN)      :: Fluxes    !< \param [in] fluxes fluxes type
    CLASS(timedisc_base), INTENT(IN)    :: Timedisc  !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP), POINTER             :: Header,IO !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    CALL this%WriteDataAttributes(Mesh,IO)

#ifdef PARALLEL
    this%first = .FALSE.
#endif
  END SUBROUTINE WriteDataset_binary

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

  !> \public Closes the file I/O
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fileio_binary), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    SELECT TYPE(df=>this%datafile)
    CLASS IS(filehandle_mpi)
      CALL MPI_Type_free(df%cfiletype,this%err)
      CALL MPI_Type_free(df%filetype,this%err)
    END SELECT
#endif
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE fileio_binary_mod
