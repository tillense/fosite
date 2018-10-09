!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_xdmf.f90                                                   #
!#                                                                           #
!# Copyright (C) 2013-2015                                                   #
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
!----------------------------------------------------------------------------!
!> \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief module for XDMF file I/O
!!
!! The xdmf file format [1] carries light data in a xml file and can store
!! heavy data in binary or hdf5 files. This implementation stores heavy data
!! in binary files (See fileio_binary.f90).
!! \todo check if fortran streams are not possible at nec aurora
!! XDMF file I/O needs Fortran Stream IO (F2003 standard). The exception are
!! NEC SX Vector machines. Here can the unformatted output be used, if the
!! record control parts are switched off. To switch it off, choose a high
!! unit number, e.g. 5555, and export the following variable in the batch
!! script:
!!
!! export F_NORCW=5555
!!
!! If it is a MPI calculation, you have addionally to define
!!
!! export MPIEXPORT="F_NORCW"
!!
!! Otherwise the variable will not be exported to the eviornment of the
!! fosite runtime. Now choose this unit number for the XDMF output. Set the
!! key "/datafile/unit" to 5555.
!!
!! [1]: http://www.xdmf.org/index.php/XDMF_Model_and_Format
!!
!! \extends fileio_gnuplot
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_xdmf_mod
  USE fileio_base_mod
  USE fileio_binary_mod
  USE geometry_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE timedisc_base_mod
  USE fluxes_base_mod
  USE sources_base_mod
  USE common_dict
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER, PARAMETER         :: LF = ACHAR(10)        !< line feed

  TYPE, EXTENDS(fileio_binary) :: fileio_xdmf
    CHARACTER(LEN=14)          :: endian_xdmf      !< endianness string
    INTEGER                    :: unit_xdmf
  CONTAINS
    PROCEDURE :: InitFileIO_xdmf
    PROCEDURE :: IterateDict
    PROCEDURE :: WriteXMF
    PROCEDURE :: WriteNode_xdmf
    PROCEDURE :: WriteKey_xdmf
    PROCEDURE :: WriteAttribute
    PROCEDURE :: WriteMeshXML
    PROCEDURE :: WriteVector
    PROCEDURE :: WriteDataItem
    PROCEDURE :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fileio_xdmf
       ! constants
       ! methods
  !--------------------------------------------------------------------------!

CONTAINS
  !> \public Constructor for the xdmf file I/O
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs,
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_xdmf(this,Mesh,Physics,Timedisc,Sources,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf),   INTENT(INOUT) :: this     !< \param [in,out] fileio class
    CLASS(mesh_base),     INTENT(IN)    :: Mesh     !< \param [in] Mesh class
    CLASS(physics_base),  INTENT(IN)    :: Physics  !< \param [in] Physics class
    CLASS(timedisc_base), INTENT(IN)    :: Timedisc !< \param [in] Timedisc class
    CLASS(sources_base),  INTENT(IN), POINTER &
                                        :: Sources  !< \param [in] Sources class
    TYPE(Dict_TYP),       INTENT(IN), POINTER &
                                        :: config
    TYPE(Dict_TYP),       INTENT(IN), POINTER &
                                        :: IO       !< \param [in] IO Dictionary for I/O
    !------------------------------------------------------------------------!
    INTEGER               :: err
    !------------------------------------------------------------------------!
    !\attention interface changed, however "xdmf" is not passed anymore
    CALL this%InitFileIO_binary(Mesh,Physics,Timedisc,Sources,config,IO)
!    CALL InitFileIO_binary(this,Mesh,Physics,IO,fmt,fpath,filename,stoptime,dtwall,&
!       count,fcycles,unit,"xdmf")

    CALL this%GetEndianness(this%endian_xdmf, '"Little"', '"Big"')
    IF(this%realsize.GT.8) &
      CALL this%Error("WriteXMF","Only single and double precision are allowed")

    WRITE(this%realfmt,'(A,I1,A)',IOSTAT=err) '"',this%realsize,'"'

#ifndef PARALLEL
#if defined(NECSXAURORA)
    CALL this%Warning("InitFileIO_xdmf","On NEC SX the Environment variable " &
      // "'F_NORCW' and")
    CALL this%Warning("InitFileIO_xdmf","the key 'datafile/unit' have to " &
      // "be set to the same (high) unit number (e.g. 5555).")
    CALL this%Warning("InitFileIO_xdmf","Otherwise there will be record " &
      // "format indentifier in the xmf file.")
#endif
#endif
    CALL this%WriteXMF(Mesh,IO)

  END SUBROUTINE InitFileIO_xdmf


  !> Iterate the dictionary and run a Subroutine on every node
  RECURSIVE SUBROUTINE IterateDict(this,Mesh,config,offset,filename,path)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT)  :: this     !< \param [in,out] fileio class
    CLASS(mesh_base),   INTENT(IN)     :: Mesh     !< \param [in] mesh class
    TYPE(Dict_TYP),     INTENT(IN), POINTER  &
                                       :: config   !< \param [in] config dict
    CHARACTER(LEN=*),   INTENT(IN)     :: filename
    CHARACTER(LEN=*),   INTENT(IN), OPTIONAL &
                                       :: path     !< \param [in] path
    INTEGER, INTENT(INOUT)             :: offset   !< \param [in,out] offset
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN)        :: str, key
    TYPE(Dict_TYP),POINTER             :: node
    !------------------------------------------------------------------------!
    IF(PRESENT(path)) THEN
        str = path
    ELSE
        str = ""
    ENDIF
    node => config
    DO WHILE(ASSOCIATED(node))
      key = TRIM(str)//"/"//TRIM(GetKey(node))

      CALL this%WriteNode_xdmf(Mesh,key,node,offset,filename)

      IF(HasChild(node)) THEN
        CALL this%IterateDict(Mesh, GetChild(node), offset, filename, key)
      END IF
      node=>GetNext(node)
    END DO
  END SUBROUTINE IterateDict

  !> Write the xdmf key
  !!
  !! \attention There is also a WriteKey function in the parent class
  !! (binary) which should not be inherited
  SUBROUTINE WriteKey_xdmf(this,offset,key,type,bytes,dims)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this     !< \param [in,out]
    CHARACTER(LEN=*),   INTENT(IN)    :: key
    INTEGER,            INTENT(INOUT) :: offset
    INTEGER,            INTENT(IN)    :: type
    INTEGER,            INTENT(IN)    :: bytes
    INTEGER,DIMENSION(5), OPTIONAL, INTENT(IN) &
                                      :: dims
    !------------------------------------------------------------------------!
    INTEGER :: keylen,l
    !------------------------------------------------------------------------!
    keylen = LEN_TRIM(key)
    IF(PRESENT(dims)) THEN
      l = 3
      IF(dims(4).GT.1) l = 4
      IF(dims(5).GT.1) l = 5
    ELSE
      l = 0
    END IF
    offset = offset + (3+l) * this%intsize + keylen
  END SUBROUTINE WriteKey_xdmf


  !> Write the xdmf node
  !!
  !! \attention There is also a WriteNode function in the parent class
  !! (binary) which should not be inherited
  SUBROUTINE WriteNode_xdmf(this,Mesh,key,node,offset,filename)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh     !< \param [in] mesh mesh type
    CHARACTER(LEN=MAX_CHAR_LEN), INTENT(IN) &
                                      :: key
    TYPE(Dict_TYP),     INTENT(IN), POINTER &
                                      :: node     !< \param [in] data node
    INTEGER ,           INTENT(INOUT) :: offset   !< \param [in,out] offset
    CHARACTER(LEN=*),   INTENT(IN)    :: filename
    !------------------------------------------------------------------------!
    REAL,DIMENSION(:,:,:),    POINTER :: ptr3 => null()
    REAL,DIMENSION(:,:,:,:),  POINTER :: ptr4 => null()
    REAL,DIMENSION(:,:,:,:,:),POINTER :: ptr5 => null()
    INTEGER,DIMENSION(5)              :: dims
    INTEGER                           :: bytes
    CHARACTER(LEN=1),DIMENSION(:),POINTER  &
                                      :: val
    LOGICAL                           :: ref
    INTEGER                           :: type
    !------------------------------------------------------------------------!
    ref = .FALSE.
    dims(:) = 1

    type = GetDataType(node)
    SELECT CASE(type)
    CASE(DICT_REAL_P)
      bytes = this%realsize
      type = DICT_REAL
    CASE(DICT_REAL_THREED)
      CALL GetAttr(node,key,ptr3)
      dims(1:3) = SHAPE(ptr3)
    CASE(DICT_REAL_FOURD)
      CALL GetAttr(node,key,ptr4)
      dims(1:4) = SHAPE(ptr4)
    CASE(DICT_REAL_FIVED)
      CALL GetAttr(node,key,ptr5)
      dims(1:5) = SHAPE(ptr5)
    CASE DEFAULT
      val => GetData(node)
      IF(ASSOCIATED(val)) THEN
        bytes = SIZE(val)
      ELSE
        bytes = 0
      END IF
    END SELECT

    IF(PRODUCT(dims).GT.1) THEN
      CALL this%SetMeshDims(Mesh,dims)
      bytes = PRODUCT(dims) * this%realsize
      CALL this%WriteKey_xdmf(offset,key,type,bytes,dims)
    ELSE
      CALL this%WriteKey_xdmf(offset,key,type,bytes)
    END IF

    SELECT CASE(GetDataType(node))
    CASE(DICT_REAL_THREED)
      CALL this%WriteAttribute(Mesh,TRIM(key),dims(1:3),filename,offset,ref)
! It is not clear, how other datatypes than DICT_REAL_TWOD should be mapped in
! XDMF.
!    CASE(DICT_REAL_THREED)
!      CALL this%WriteAttribute(Mesh,TRIM(key),dims(1:3),filename,offset,ref)
!    CASE(DICT_REAL_FOURD)
!      CALL this%WriteAttribute(Mesh,TRIM(key),dims(1:4),filename,offset,ref)
    END SELECT

    offset = offset + bytes

  END SUBROUTINE WriteNode_xdmf


  !> Writes description of data item in xml syntax
  !!
  SUBROUTINE WriteDataItem(this,dims,filename,offset)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CHARACTER(LEN=*),   INTENT(IN)    :: dims      !< \param [in] dims
    CHARACTER(LEN=*),   INTENT(IN)    :: filename  !< \param [in] filename
    INTEGER,            INTENT(IN)    :: offset    !< \param [in] offset
    !------------------------------------------------------------------------!
    CHARACTER(LEN=16)                 :: seek
    INTEGER                           :: err
    !------------------------------------------------------------------------!
    WRITE(this%unit, IOSTAT=err)&
           '<DataItem Dimensions=' // TRIM(dims) // ' ' &
        // 'NumberType="Float" Precision=' // TRIM(this%realfmt) // LF &
        // 'Format="Binary" Endian=' // TRIM(this%endian_xdmf)
    WRITE(seek,'(I16)',IOSTAT=err) offset
    WRITE(this%unit, IOSTAT=err)&
       LF // 'Seek="' // TRIM(ADJUSTL(seek)) // '"'
    WRITE(this%unit, IOSTAT=err)&
           '>' // LF &
        // TRIM(filename) // LF &
        // '</DataItem>' // LF

  END SUBROUTINE WriteDataItem

  !> Writes description of data item in xml syntax
  !!
  SUBROUTINE WriteAttribute(this,Mesh,name,dims,filename,offset,ref)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this       !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh       !< \param [in] mesh mesh type
    CHARACTER(LEN=*),   INTENT(IN)    :: name       !< \param [in] name
    CHARACTER(LEN=*),   INTENT(IN)    :: filename   !< \param [in] filename
    INTEGER,            INTENT(IN), DIMENSION(:) &
                                      :: dims       !< \param [in] dims
    INTEGER,            INTENT(IN)    :: offset     !< \param [in,out] offset
    LOGICAL,            INTENT(IN)    :: ref        !< \param [in] ref
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)                 :: type,center
    INTEGER                           :: dsize, err
    !------------------------------------------------------------------------!
    !data size
    dsize = SIZE(dims)

    ! TODO CHECK THIS - 3D HACK
    ! If the data has mesh sizes, it is cell data
    IF((dims(1).EQ.Mesh%INUM).AND.&
       (dims(2).EQ.Mesh%JNUM).AND.&
       (dims(3).EQ.Mesh%KNUM)) THEN
      center = "Cell"
      dsize = dsize - 3
    ELSE IF((dims(1).EQ.Mesh%INUM+1).AND.&
            (dims(2).EQ.Mesh%JNUM+1).AND.&
            (dims(3).EQ.Mesh%KNUM+1)) THEN
      center = "Node"
      dsize = dsize - 3
    ELSE ! otherwise data for the whole grid
      center = "Grid"
    END IF

    ! The remaining dimensions define the data type
    SELECT CASE(dsize)
    CASE(0)
        type = "Scalar"
    CASE(1)
        type = "Vector"
    CASE DEFAULT
        type = "Matrix"
    END SELECT

    WRITE(this%unit, IOSTAT=err)&
         '<Attribute Name="' // TRIM(name) // '" '&
      // 'AttributeType="' // TRIM(type) // '" ' &
      // 'Center="' // TRIM(center) // '">' // LF

    CALL this%WriteDataItem(GetDimsStr(Mesh,dims), TRIM(filename), offset)

    WRITE(this%unit, IOSTAT=err)&
         '</Attribute>' // LF
  END SUBROUTINE WriteAttribute

  !> Writes the mesh to file
  !!
  SUBROUTINE WriteVector(this,Mesh,name,dims,ref1,ref2,ref3,step)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this       !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh       !< \param [in] mesh mesh type
    CHARACTER(LEN=*),   INTENT(IN)    :: name       !< \param [in] name
    INTEGER,            INTENT(IN), DIMENSION(3) &
                                      :: dims       !< \param [in] dims
    CHARACTER(LEN=*),   INTENT(IN)    :: ref1,ref2,ref3  !< \param [in] name
    CHARACTER(LEN=*),   INTENT(IN)    :: step       !< \param [in]
    !------------------------------------------------------------------------!
    CHARACTER(LEN=128)                :: dstr
    INTEGER                           :: err
    !------------------------------------------------------------------------!
    dstr = GetDimsStr(Mesh,dims)
    WRITE(this%unit, IOSTAT=err)&
         '<Attribute Name="' // TRIM(name) // '" ' &
      // 'AttributeType="Vector" Center="Cell">' // LF&
      // '<DataItem ItemType="Function" Function="JOIN($0, $1, $2)" '&
      // 'Dimensions=' // TRIM(dstr) // '>' // LF &
      ! first component
      // '<DataItem Reference="/Xdmf/Domain/Grid'&
      // "/Grid[@Name='step" // step // "']" &
      // '/Attribute[@Name=' &
      // "'" // TRIM(ref1) // "'" // ']/DataItem[1]"/>' // LF &
      ! second component
      // '<DataItem Reference="/Xdmf/Domain/Grid'&
      // "/Grid[@Name='step" // step // "']" &
      // '/Attribute[@Name=' &
      // "'" // TRIM(ref2) // "'" // ']/DataItem[1]"/>' // LF &
      ! third component
      // '<DataItem Reference="/Xdmf/Domain/Grid'&
      // "/Grid[@Name='step" // step // "']" &
      // '/Attribute[@Name=' &
      // "'" // TRIM(ref3) // "'" // ']/DataItem[1]"/>' // LF &
      ! end
      // '</DataItem>' // LF &
      // '</Attribute>' // LF
  END SUBROUTINE WriteVector

  !> Writes the mesh to file
  !!
  SUBROUTINE WriteMeshXML(this,Mesh,filename,offset)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh     !< \param [in] mesh mesh type
    CHARACTER(LEN=*),   INTENT(IN)    :: filename !< \param [in] filename
    INTEGER,            INTENT(INOUT) :: offset   !< \param [in,out] offset
    !------------------------------------------------------------------------!
    CHARACTER(LEN=64) :: dims
    CHARACTER(LEN=8)  :: inum, jnum, knum
    INTEGER           :: err
    !------------------------------------------------------------------------!
    WRITE(inum,'(I8)') Mesh%INUM+1
    WRITE(jnum,'(I8)') Mesh%JNUM+1
    WRITE(knum,'(I8)') Mesh%KNUM+1
    WRITE(dims,'(A,A,A,A,A)') TRIM(ADJUSTL(knum)), ' ',TRIM(ADJUSTL(jnum)), &
                              ' ',TRIM(ADJUSTL(inum))

! TODO 3D HACK should not be needed anymore, since now everything is 3D now
! \attention here was a bianglespherical hack
    WRITE(this%unit, IOSTAT=err) &
            '<Topology TopologyType="3DSMesh" ' &
         // 'NumberOfElements="' // TRIM(dims) // '"/>' // LF &
         // '<Geometry GeometryType="X_Y_Z">' // LF
    WRITE(this%unit, IOSTAT=err) &
         '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name='//"'/mesh/grid_x'"//']/DataItem[1]"/>' // LF &
      // '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name='//"'/mesh/grid_y'"//']/DataItem[1]"/>' // LF &
      // '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name='//"'/mesh/grid_z'"//']/DataItem[1]"/>' // LF

    WRITE(this%unit, IOSTAT=err) &
         '</Geometry>' // LF

  END SUBROUTINE WriteMeshXML

  !> Main routine to write all data to xmf file
  !!
  SUBROUTINE WriteXMF(this,Mesh,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this   !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh   !< \param [in] mesh mesh type
    TYPE(Dict_TYP),     INTENT(IN), POINTER &
                                      :: IO     !< \param [in] io I/O dictionary
    !------------------------------------------------------------------------!
    INTEGER                 :: i, offset, err
    REAL                    :: ftime
    CHARACTER(LEN=4)        :: step
    CHARACTER(LEN=32)       :: time
    CHARACTER(LEN=256)      :: filename
    TYPE(Dict_TYP), POINTER :: meshIO => Null()
    !------------------------------------------------------------------------!

    IF (this%GetRank().EQ.0) THEN
      ! write a xmf-file
      OPEN(this%unit, FILE=TRIM(this%filename)//'.xmf', &
           STATUS     = 'REPLACE',      &
           ACTION     = 'WRITE',        &
           ACCESS     = 'STREAM',       &
           POSITION   = 'REWIND',       &
           IOSTAT     = err)
       IF (err.NE. 0) CALL this%Error("WriteXMF","Can't open xmf-file")

      WRITE(this%unit, IOSTAT=err)&
           '<?xml version="1.0" ?>' // LF &
        // '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' // LF &
        // '<Xdmf Version="2.0">' // LF &
        // '<Domain>' // LF &
        // '<Grid Name="mesh" GridType="Collection" CollectionType="Temporal">' // LF

! this works in paraview, but not in visit :(. So we have to use the simpler method for now
!        // '<Time TimeType="List">' // LF &
!        // '<DataItem Format="XML" NumberType="Float" Dimensions="11">' // LF &
!        // '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0' // LF &
!        // '</DataItem>' // LF &
!        // '</Time>' // LF
!      CALL WriteDataItem_xdmf(this, '"1"', TRIM(filename)//'_'//step//'.bin')
!      WRITE(this%unit, IOSTAT=err)&

      CALL GetAttr(IO,"mesh",meshIO)

      DO i=0,this%count
        offset = 13
        ftime = this%stoptime/(this%count)*i
        WRITE(step,'(I4.4)') i
        WRITE(time,'(ES16.9)') ftime
        WRITE(filename, '(A)') TRIM(this%filename)//"_"//step//".bin"

        WRITE(this%unit, IOSTAT=err)&
             '<Grid Name="step' // step // '" GridType="Uniform">' // LF &
          // '<Time Value="' // TRIM(time) // '" />' // LF

        CALL this%WriteMeshXML(Mesh,filename,offset)
        CALL this%WriteVector(Mesh,"/timedisc/velocity", &
                              (/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/), &
                              "/timedisc/xvelocity","/timedisc/yvelocity", &
                              "/timedisc/zvelocity", step)
        CALL this%IterateDict(Mesh, IO, offset, filename)

        WRITE(this%unit, IOSTAT=err) &
             '</Grid>' // LF
      END DO

      WRITE(this%unit, IOSTAT=err) &
           '</Grid>' // LF


      WRITE(this%unit, IOSTAT=err) &
           '</Domain>' // LF &
        // '</Xdmf>' // LF

      IF(err.NE.0) CALL this%Error("WriteXMF", "Can't write xmf-file")
      CLOSE(this%unit,IOSTAT=err)
      IF(err.NE.0) CALL this%Error("WriteXMF", "Can't close xmf-file")
   END IF

  END SUBROUTINE WriteXMF


  FUNCTION GetDimsStr(Mesh,dims) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    INTEGER,DIMENSION(:), INTENT(IN) :: dims
    CHARACTER(LEN=128)               :: res
    !------------------------------------------------------------------------!
    INTEGER                          :: i, l
    CHARACTER(LEN=128), DIMENSION(:), ALLOCATABLE &
                                     :: buf
    !------------------------------------------------------------------------!
    l = SIZE(dims)
    ALLOCATE(buf(l))
    IF(l.GE.3) THEN
      WRITE(buf(1),'(I8)') dims(3)
      WRITE(buf(2),'(I8)') dims(2)
      WRITE(buf(3),'(I8)') dims(1)
      DO i=4,l
        WRITE(buf(i),'(I8)') dims(i)
      END DO
      SELECT CASE(l)
      CASE(2)
        WRITE(res,'(A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),'"'
      CASE(3)
        WRITE(res,'(A,A,A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),' ',TRIM(ADJUSTL(buf(3))),'"'
      CASE(4)
        WRITE(res,'(A,A,A,A,A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),' ',TRIM(ADJUSTL(buf(3))), &
          ' ',TRIM(ADJUSTL(buf(4))),'"'
      CASE(5)
        WRITE(res,'(A,A,A,A,A,A,A,A,A,A,A)') '"',TRIM(ADJUSTL(buf(1))),&
          ' ',TRIM(ADJUSTL(buf(2))),' ',TRIM(ADJUSTL(buf(3))), &
          ' ',TRIM(ADJUSTL(buf(4))),' ',TRIM(ADJUSTL(buf(5))),'"'
      END SELECT
    ELSE
      WRITE(buf(1),'(I8)') dims(1)
      WRITE(res,'(A,A,A)') '"',TRIM(ADJUSTL(buf(1))),'"'
    END IF
    DEALLOCATE(buf)
  END FUNCTION GetDimsStr

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this
    !-----------------------------------------------------------------------!
    CALL this%fileio_binary%Finalize()
  END SUBROUTINE

END MODULE fileio_xdmf_mod
