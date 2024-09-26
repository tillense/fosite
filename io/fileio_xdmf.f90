!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_xdmf.f90                                                   #
!#                                                                           #
!# Copyright (C) 2013-2024                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!! \brief module for XDMF file I/O
!!
!! The xdmf file format [1,4] carries light data in a xml file and can store
!! heavy data in binary or hdf5 files. This implementation stores heavy data
!! in binary files (See fileio_binary.f90).
!!
!! If the xdmf file fails loading in paraview [2] try validation using
!!
!! > `xmllint --noout --dtdvalid Xdmf.dtd <example_data_file.xmf>`
!!
!! with the dtd file from [3] .
!!
!! [1]: http://www.xdmf.org/index.php/XDMF_Model_and_Format
!! [2]: https://www.paraview.org
!! [3]: https://gitlab.kitware.com/xdmf/xdmf/blob/master/Xdmf.dtd
!! [4]: https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/data_into_visit/XdmfFormat.html#
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
  INTEGER, PARAMETER           :: BLK_INDENT = 2        !< block indentation
  CHARACTER, PARAMETER         :: SP = ACHAR(32)        !< space
  CHARACTER, PARAMETER         :: LF = ACHAR(10)        !< line feed
  CHARACTER(LEN=2), PARAMETER  :: TB = REPEAT(SP,BLK_INDENT)
  !--------------------------------------------------------------------------!
  !> FileIO gnuplot class
  TYPE, EXTENDS(fileio_binary) :: fileio_xdmf
    TYPE(filehandle_fortran)   :: xmffile          !< paraview master file
    CHARACTER(LEN=14)          :: endian_xdmf      !< endianness string
  CONTAINS
    PROCEDURE :: InitFileIO_deferred => InitFileIO_xdmf
    PROCEDURE :: IterateDict
    PROCEDURE :: WriteXMF
    PROCEDURE :: WriteNode_xdmf
    PROCEDURE :: WriteKey_xdmf
    PROCEDURE :: WriteAttribute
    PROCEDURE :: WriteMeshXML
    PROCEDURE :: WriteVector
    PROCEDURE :: WriteDataItem
    FINAL :: Finalize
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
    CLASS(sources_base),  INTENT(IN)    :: Sources  !< \param [in] Sources class
    TYPE(Dict_TYP),       INTENT(IN), POINTER &
                                        :: config
    TYPE(Dict_TYP),       INTENT(IN), POINTER &
                                        :: IO       !< \param [in] IO Dictionary for I/O
    !------------------------------------------------------------------------!
    INTEGER :: idx(1)
    !------------------------------------------------------------------------!
    !\attention interface changed, however "xdmf" is not passed anymore
    IF (.NOT.this%Initialized()) &
      CALL this%InitFileio(Mesh,Physics,Timedisc,Sources,config,IO,"xdmf","bin",textfile=.FALSE.)

    CALL this%fileio_binary%InitFileIO(Mesh,Physics,Timedisc,Sources,config,IO)

    CALL this%GetEndianness(this%endian_xdmf, '"Little"', '"Big"')

    SELECT CASE (Mesh%NDIMS)
    CASE(2,3)
      ! 2D & 3D -> do nothing
    CASE DEFAULT
      CALL this%Error("fileio_xdmf::InitFileIO_xdmf","only 2D and 3D mesh is supported")
    END SELECT

    ! create file handle for the xmf file; only on rank 0 in parallel mode
#ifdef PARALLEL
    IF (this%GetRank().EQ.0) &
#endif
      CALL this%xmffile%InitFilehandle(this%datafile%filename,this%datafile%path,"xmf",textfile=.FALSE.,onefile=.TRUE.,cycles=1)

    CALL this%WriteXMF(Mesh,IO)
  END SUBROUTINE InitFileIO_xdmf


  !> Iterate the dictionary and run a Subroutine on every node
  RECURSIVE SUBROUTINE IterateDict(this,Mesh,config,offset,filename,path,indent)
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
    INTEGER                            :: indent   !< \param [in] indent indentation level
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

      CALL this%WriteNode_xdmf(Mesh,key,node,offset,filename,indent)

      IF(HasChild(node)) THEN
        CALL this%IterateDict(Mesh, GetChild(node), offset, filename, key, indent)
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
  SUBROUTINE WriteNode_xdmf(this,Mesh,key,node,offset,filename,indent)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh     !< \param [in] mesh mesh type
    CHARACTER(LEN=MAX_CHAR_LEN), INTENT(IN) &
                                      :: key
    TYPE(Dict_TYP),     INTENT(IN), POINTER &
                                      :: node     !< \param [in] data node
    INTEGER ,           INTENT(INOUT) :: offset   !< \param [in,out] offset
    CHARACTER(LEN=*),   INTENT(IN)    :: filename !< \param [in] filename
    INTEGER,            INTENT(IN)    :: indent   !< \param [in] indentation level
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
      CALL this%SetOutputDims(Mesh,dims)
      bytes = PRODUCT(dims) * this%realsize
      CALL this%WriteKey_xdmf(offset,key,type,bytes,dims)
    ELSE
      CALL this%WriteKey_xdmf(offset,key,type,bytes)
    END IF

    SELECT CASE(GetDataType(node))
    CASE(DICT_REAL_THREED)
      CALL this%WriteAttribute(Mesh,TRIM(key),dims(1:3),filename,offset,ref,indent)
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
  SUBROUTINE WriteDataItem(this,dims,filename,offset,indent)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CHARACTER(LEN=*),   INTENT(IN)    :: dims      !< \param [in] dims
    CHARACTER(LEN=*),   INTENT(IN)    :: filename  !< \param [in] filename
    INTEGER,            INTENT(IN)    :: offset    !< \param [in] offset
    INTEGER,            INTENT(IN)    :: indent     !< \param [in] indentation level
    !------------------------------------------------------------------------!
    CHARACTER(LEN=16)                 :: seek
    !------------------------------------------------------------------------!
    WRITE(seek,'(I16)',IOSTAT=this%err) offset

    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        REPEAT(TB,indent) // '<DataItem Dimensions=' // TRIM(dims) // SP &
        // 'NumberType="Float" Precision=' // TRIM(this%realfmt)
    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        SP // 'Format="Binary" Endian=' // TRIM(this%endian_xdmf) // &
        SP // 'Seek="' // TRIM(ADJUSTL(seek)) // '">' // LF
    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        REPEAT(TB,indent+1) // TRIM(filename) // LF // &
        REPEAT(TB,indent) // '</DataItem>' // LF
  END SUBROUTINE WriteDataItem

  !> Writes description of data item in xml syntax
  !!
  SUBROUTINE WriteAttribute(this,Mesh,name,dims,filename,offset,ref,indent)
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
    INTEGER,            INTENT(IN)    :: indent     !< \param [in] indentation level
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)                 :: dtype,center,dstr
    INTEGER                           :: dsize, err
    !------------------------------------------------------------------------!
    !data size
    dsize = SIZE(dims)
    ! default (no mesh data)
    center = "Grid"

    ! check for mesh data
    IF (dsize.GE.3) THEN
      IF (ALL(dims(1:3).EQ.(/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/))) &
         center = "Cell"
      IF (ALL(dims(1:3).EQ.(/Mesh%INUM+Mesh%IP1,Mesh%JNUM+Mesh%JP1,Mesh%KNUM+Mesh%KP1/))) &
         center = "Node"
      IF (TRIM(center).NE."Grid") &
         dsize = dsize - 3
    END IF

    SELECT CASE(dsize)
    CASE(0)
        dtype = "Scalar"
    CASE(1)
        dtype = "Vector"
    CASE(2)
        dtype = "Matrix"
    CASE DEFAULT
        CALL this%Error("fileio_xdmf::WriteAttribute","data type currently not supported")
    END SELECT

    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
      REPEAT(TB,indent) // '<Attribute Name="' // TRIM(name) // '" ' &
      // 'AttributeType="' // TRIM(dtype) // '" ' &
      // 'Center="' // TRIM(center) // '">' // LF

    CALL GetDimsStr(dims,dstr)
    CALL this%WriteDataItem(TRIM(dstr), TRIM(filename), offset,indent+1)

    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
      REPEAT(TB,indent) // '</Attribute>' // LF
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
    CHARACTER(LEN=32)                 :: dstr
    INTEGER                           :: err
    !------------------------------------------------------------------------!
    CALL GetDimsStr(dims,dstr)
    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
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
  SUBROUTINE WriteMeshXML(this,Mesh,filename,indent)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_xdmf), INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),   INTENT(IN)    :: Mesh     !< \param [in] mesh mesh type
    CHARACTER(LEN=*),   INTENT(IN)    :: filename !< \param [in] filename
    INTEGER                           :: indent   !< \param [in] indent indentation level
    !------------------------------------------------------------------------!
    INTEGER           :: dims(3)
    CHARACTER(LEN=32) :: dstr
    CHARACTER(LEN=14) :: gstr(3) = (/"'/mesh/grid_x'","'/mesh/grid_y'","'/mesh/grid_z'"/)
    INTEGER           :: i
    !------------------------------------------------------------------------!
    dims = (/ Mesh%INUM+Mesh%IP1, Mesh%JNUM+Mesh%JP1, Mesh%KNUM+Mesh%KP1 /)
    CALL GetDimsStr(dims(:),dstr)

    SELECT CASE(Mesh%NDIMS)
    CASE(2)
      WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        REPEAT(TB,indent) // '<Topology TopologyType="2DSMesh" ' // 'NumberOfElements=' // TRIM(dstr) // '/>' // LF // &
        REPEAT(TB,indent) // '<Geometry GeometryType="X_Y_Z">' // LF
    CASE(3)
      WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        REPEAT(TB,indent) // '<Topology TopologyType="3DSMesh" ' // 'NumberOfElements=' // TRIM(dstr) // '/>' // LF // &
        REPEAT(TB,indent) // '<Geometry GeometryType="X_Y_Z">' // LF
    CASE DEFAULT
      CALL this%Error("fileio_xdmf::WriteMeshXML","only 2D and 3D mesh is currently supported")
    END SELECT

    ! register all 3 coordinates even for 2D grids, because 2D curved surfaces may have 3D cartesian coordinates
    ! e.g., surface of a sphere
    DO i=1,3
        WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
          REPEAT(TB,indent+1) // '<DataItem Reference="/Xdmf/Domain/Grid/Grid/Attribute[@Name='//gstr(i)//']/DataItem[1]"/>' // LF
    END DO

    WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
         REPEAT(TB,indent) // '</Geometry>' // LF
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
    this%err = 0
    IF (this%GetRank().EQ.0) THEN
      ! write a xmf-file
      CALL this%xmffile%OpenFile(REPLACE,this%step)

      WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
           '<?xml version="1.0" ?>' // LF &
        // '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' // LF &
        // '<Xdmf Version="2.0">' // LF &
        // TB // '<Domain>' // LF &
        // TB // TB // '<Grid Name="mesh" GridType="Collection" CollectionType="Temporal">' // LF

! this works in paraview, but not in visit :(. So we have to use the simpler method for now
!        // '<Time TimeType="List">' // LF &
!        // '<DataItem Format="XML" NumberType="Float" Dimensions="11">' // LF &
!        // '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0' // LF &
!        // '</DataItem>' // LF &
!        // '</Time>' // LF
!      CALL WriteDataItem_xdmf(this, '"1"', TRIM(filename)//'_'//step//'.bin')

      CALL GetAttr(IO,"mesh",meshIO)

      DO i=0,this%count
        offset = 13
        ftime = this%stoptime/(this%count)*i
        WRITE(step,'(I4.4)') i
        WRITE(time,'(ES16.9)') ftime
        WRITE(filename, '(A)') TRIM(this%datafile%GetFilename(i))

        WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
          REPEAT(TB,3) // '<Grid Name="step' // step // '" GridType="Uniform">' // LF // &
          REPEAT(TB,4) // '<Time Value="' // TRIM(time) // '" />' // LF

        CALL this%WriteMeshXML(Mesh,filename,indent=5)
! quick HACK to bind velocity components into vector structure
!         CALL this%WriteVector(Mesh,"/timedisc/velocity", &
!                               (/Mesh%INUM,Mesh%JNUM,Mesh%KNUM/), &
!                               "/timedisc/xvelocity","/timedisc/yvelocity", &
!                               "/timedisc/zvelocity", step)
        CALL this%IterateDict(Mesh, IO, offset, filename,indent=5)

        WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
          REPEAT(TB,3) // '</Grid>' // LF
      END DO

      WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        TB // TB // '</Grid>' // LF


      WRITE(UNIT=this%xmffile%GetUnitNumber(),IOSTAT=this%err) &
        TB // '</Domain>' // LF &
        // '</Xdmf>' // LF

      IF(this%err.NE.0) CALL this%Error("fileio_xdmf::WriteXMF", "Can't write xmf-file")
      CALL this%xmffile%CloseFile(this%step)
    END IF
  END SUBROUTINE WriteXMF

  !> \public determines the string for the dimension attribute, i.e. array dimensions of mesh data
  PURE SUBROUTINE GetDimsStr(dims,res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER,DIMENSION(3), INTENT(IN) :: dims
    CHARACTER(LEN=*), INTENT(OUT)    :: res
    !------------------------------------------------------------------------!
    INTEGER                          :: i,idx(3)
    CHARACTER(LEN=8)                 :: dim_str
    !------------------------------------------------------------------------!
    ! reverse dimensional order (probably because of Fortran vs. C array order)
    idx(1:3) = dims(3:1:-1)

    res = ''
    DO i=1,3
      IF (idx(i).GT.1) THEN ! skip suppressed dimension in 2D case
        WRITE(dim_str,'(I8)') idx(i)
        ! append to resulting string
        res = TRIM(res) // ' ' // TRIM(ADJUSTL(dim_str))
      END IF
    END DO
    res = '"' // TRIM(ADJUSTL(res)) // '"'
  END SUBROUTINE GetDimsStr

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    TYPE(fileio_xdmf), INTENT(INOUT) :: this
    !-----------------------------------------------------------------------!
    ! do nothing specific for xdmf
    ! finalizer of fileio_binary is called automatically
  END SUBROUTINE

END MODULE fileio_xdmf_mod
