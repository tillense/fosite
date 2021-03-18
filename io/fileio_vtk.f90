!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_vtk.f90                                                    #
!#                                                                           #
!# Copyright (C) 2010-2021                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief I/O for VTK files in XML format (vtkStructuredGrid)
!!
!! This module implements VTK file I/O to write vts files (vtkStructuredGrid)
!! in XML syntax. Each output file contains a single timestep.
!! In case of parallel computation there is one file per job (and timestep)
!! and only one global container file (pvts) which groups all vts files.
!! Additionally a PVD file is generated, which can be used to load
!! the data in Paraview.
!!
!! References:
!! - Visualization Toolkit (VTK)
!!   * Wiki: http://www.vtk.org/Wiki/VTK_XML_Formats
!!   * file format specifications: http://www.vtk.org/VTK/img/file-formats.pdf
!!   * VTK XML Reader/Writer: http://www.vtk.org/doc/nightly/html/IOXMLInformationFormat.html
!! - Paraview, http://www.paraview.org/Wiki/ParaView/Data_formats
!! - Visit, https://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
!!
!! \extends fileio_common
!! \ingroup fileio
!----------------------------------------------------------------------------!
MODULE fileio_vtk_mod
  USE fileio_base_mod
  USE geometry_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE timedisc_base_mod
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
  !--------------------------------------------------------------------------!
   INTEGER, PARAMETER      :: MAXCOMP  = 9   !< max. of allowed components
                                             !! 9 is a tensor (rank 2, dim 3)
   INTEGER, PARAMETER      :: MAXCOLS  = 40  !< max of different output arrays
   INTEGER, PARAMETER      :: MAXKEY   = 64  !< max length of keyname
   CHARACTER, PARAMETER    :: LF = ACHAR(10) !< line feed
   !> names of fluxes
!   CHARACTER(LEN=16),DIMENSION(6),PARAMETER  :: fluxkey = (/'bflux_WEST  ', &
!                                                            'bflux_EAST  ', &
!                                                            'bflux_SOUTH ', &
!                                                            'bflux_NORTH ', &
!                                                            'bflux_BOTTOM', &
!                                                            'bflux_TOP   ' /)

  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(fileio_base) :: fileio_vtk
    PRIVATE
    INTEGER   :: IP1,JP1,KP1
    CONTAINS
    ! methods
    PROCEDURE :: InitFileio_vtk
    PROCEDURE :: OpenFile
    PROCEDURE :: CloseFile
    PROCEDURE :: WriteHeader
    PROCEDURE :: ReadHeader
    PROCEDURE :: WriteTimestamp
    PROCEDURE :: WriteParaviewFile
    PROCEDURE :: ReadTimestamp
    PROCEDURE :: WriteDataset
    PROCEDURE :: GetPrecision
    PROCEDURE :: GetOutputlist
    PROCEDURE :: GetEndianness
    PROCEDURE :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    fileio_vtk

  !--------------------------------------------------------------------------!
CONTAINS
  !> \public Constructor for the VTK file I/O
  !!
  !! Initilizes the file I/O type, filename, stoptime, number of outputs,
  !! number of files, unit number, config as a dict
  SUBROUTINE InitFileIO_vtk(this,Mesh,Physics,Timedisc,Sources,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk),   INTENT(INOUT)       :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),    INTENT(IN)          :: Mesh     !< \param [in] Mesh mesh type
    CLASS(physics_base), INTENT(IN)          :: Physics  !< \param [in] Physics Physics type
    CLASS(timedisc_base),INTENT(IN)          :: Timedisc !< \param [in] Physics Physics type
    CLASS(sources_base), INTENT(IN), POINTER :: Sources  !< \param [in] Physics Physics type
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: config   !< \param [in] IO Dictionary for I/O
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: IO       !< \param [in] IO Dictionary for I/O
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER             :: node
    INTEGER                             :: k,i,j,err
#ifdef PARALLEL
    INTEGER, DIMENSION(:), POINTER      :: sendbuf,recvbuf
    INTEGER                             :: n
#endif
    REAL, DIMENSION(:,:,:,:,:), POINTER :: corners
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    this%multfiles = .TRUE.
#endif
    ! start init form base class in beginning
    CALL this%InitFileIO(Mesh,Physics,Timedisc,Sources,config,IO,"VTK","vts")

    ! some sanity checks
    IF (this%cycles .NE. this%count+1) CALL this%Error('InitFileIO','VTK need filecycles = count+1')

    ! determine dimensions for the coordinate array
    IF (ABS(Mesh%xmax-Mesh%xmin).LT.2*EPSILON(Mesh%xmin) &
       .OR.(Mesh%INUM.EQ.1.AND.Mesh%ROTSYM.EQ.1)) THEN
      ! suppress x-direction in output
      this%IP1 = 0
    ELSE
      this%IP1 = 1
    END IF
    IF (ABS(Mesh%ymax-Mesh%ymin).LT.2*EPSILON(Mesh%ymin) &
       .OR.(Mesh%JNUM.EQ.1.AND.Mesh%ROTSYM.EQ.2)) THEN
      ! suppress y-direction in output
      this%JP1  = 0
    ELSE
      this%JP1  = 1
    END IF
    IF (ABS(Mesh%zmax-Mesh%zmin).LT.2*EPSILON(Mesh%zmin) &
       .OR.(Mesh%KNUM.EQ.1.AND.Mesh%ROTSYM.EQ.3)) THEN
      ! suppress z-direction in output
      this%KP1  = 0
    ELSE
      this%KP1  = 1
    END IF

    ! allocate memory for auxilliary data fields
    ALLOCATE(this%binout(MAXCOMP,                &
                         Mesh%IMIN:Mesh%IMAX,    &
                         Mesh%JMIN:Mesh%JMAX,    &
                         Mesh%KMIN:Mesh%KMAX),   &
             this%vtkcoords(3,                   &
                         Mesh%IMIN:Mesh%IMAX+this%IP1,  &
                         Mesh%JMIN:Mesh%JMAX+this%JP1,  &
                         Mesh%KMIN:Mesh%KMAX+this%KP1), &
             this%output(MAXCOLS),               &
             this%tsoutput(MAXCOLS),             &
             this%bflux(Physics%VNUM,6),         &
             STAT = err)
    IF (err.NE.0) &
         CALL this%Error("InitFileio_vtk", "Unable to allocate memory.")

    ! determine extend and endianness of real numbers
    CALL this%GetPrecision(this%realsize)
    ! save size of floats to string: realfmt
    WRITE(this%linebuf,FMT='(I4)',IOSTAT=this%error_io) 8*this%realsize
    WRITE(this%realfmt,FMT='(A)',IOSTAT=this%error_io) '"Float' &
          // TRIM(AdjustL(this%linebuf)) // '"'
    CALL this%GetEndianness(this%endianness,'"LittleEndian"','"BigEndian"')

    ! sanity check in case someone ignores the comment in
    ! fileio_common.f90 regarding the KIND type of this variable
    IF (BIT_SIZE(this%output(1)%numbytes) .NE. 4*8) &
       CALL this%Error("InitFileio_vtk", "output%numbytes must be a 4 byte integer !!!")

#ifdef PARALLEL
    ! allocate memory for this%extent, send & receive buffers;
    ! recvbuf and this%extent are only needed on the receiver node (rank 0)
    ! but parallel profiling using scalasca with mpich3 crashes
    ! if recvbuf is not allocated on each node hence we allocate a minimal
    ! array of size 1 on all nodes except for the receiver node;
    ! remark: according to the MPI-2 standard recvbuf must only be allocated
    !         on the receiver node
    IF (this%GetRank() .EQ. 0 ) THEN
       n = this%GetNumProcs()
       ALLOCATE(this%extent(0:(n-1)),sendbuf(6),recvbuf(6*n),STAT = this%error_io)
    ELSE
       ! see comment above
       ALLOCATE(sendbuf(6),recvbuf(1),STAT = this%error_io)
    END IF
    IF (this%error_io.NE.0) &
       CALL this%Error( "InitFileio_vtk", "Unable to allocate memory.")

    ! store information about grid extent on each node
    !> \bug This is most probably broken !
    sendbuf(1) = Mesh%IMIN
    sendbuf(2) = Mesh%IMAX+this%IP1
    sendbuf(3) = Mesh%JMIN
    sendbuf(4) = Mesh%JMAX+this%JP1
    sendbuf(5) = Mesh%KMIN
    sendbuf(6) = Mesh%KMAX+this%KP1

    ! gather information about grid extent on each node at the rank 0 node
    CALL MPI_Gather(sendbuf, 6, MPI_INTEGER, recvbuf, 6, MPI_INTEGER, &
               0, MPI_COMM_WORLD, this%error_io)

    ! store information about grid extent at the rank 0 node
    IF (this%GetRank() .EQ. 0 ) THEN
       DO i=0,n-1
          WRITE(this%extent(i),FMT='(6(I7))',IOSTAT=this%error_io)&
          recvbuf(i*6+1),recvbuf(i*6+2),recvbuf(i*6+3),recvbuf(i*6+4),recvbuf(i*6+5),recvbuf(i*6+6)
       END DO
    END IF
    ! free buffer memory
    DEALLOCATE(sendbuf,recvbuf,STAT=this%error_io)
#endif

    ! check whether mesh coordinates should be cartesian or curvilinear
    IF (this%cartcoords) THEN
       corners => Mesh%RemapBounds(Mesh%cart%corners)
    ELSE
       corners => Mesh%RemapBounds(Mesh%curv%corners)
    END IF

    ! store mesh corners in appropriate array for VTK coordinate output
    DO k=Mesh%KMIN,Mesh%KMAX+Mesh%KP1
      DO j=Mesh%JMIN,Mesh%JMAX+Mesh%JP1
         DO i=Mesh%IMIN,Mesh%IMAX+Mesh%IP1
            this%vtkcoords(1:3,i,j,k) = corners(i,j,k,1,1:3)
         END DO
      END DO
    END DO

    ! add coordinates in 1D/2D simulations with non-vanishing extent
    IF (Mesh%INUM.EQ.1.AND.this%IP1.EQ.1) THEN
      ! no transport in x-direction, but x-extent > 0
      i = Mesh%IMIN
      DO k=Mesh%KMIN,Mesh%KMAX+Mesh%KP1
         DO j=Mesh%JMIN,Mesh%JMAX+Mesh%JP1
            this%vtkcoords(1:3,i+1,j,k) = corners(i,j,k,2,1:3)
         END DO
      END DO
    END IF

    IF (Mesh%JNUM.EQ.1.AND.this%JP1.EQ.1) THEN
      ! no transport in y-direction, but y-extent > 0
      j = Mesh%JMIN
      DO k=Mesh%KMIN,Mesh%KMAX+Mesh%KP1
         DO i=Mesh%IMIN,Mesh%IMAX+Mesh%IP1
            this%vtkcoords(1:3,i,j+1,k) = corners(i,j,k,3,1:3)
         END DO
      END DO
    END IF

    IF (Mesh%KNUM.EQ.1.AND.this%KP1.EQ.1) THEN
      ! no transport in z-direction, but z-extent > 0
      k = Mesh%KMIN
      DO j=Mesh%JMIN,Mesh%JMAX+Mesh%JP1
         DO i=Mesh%IMIN,Mesh%IMAX+Mesh%IP1
            this%vtkcoords(1:3,i,j,k+1) = corners(i,j,k,5,1:3)
         END DO
      END DO
    END IF

    ! byte offset, i.e. extent of coordinate data + 4 bytes for storing
    ! the size information which must be a 4 byte integer number
    this%offset = SIZE(this%vtkcoords)*this%realsize + 4

    ! create list of output data arrays
    node => IO
    this%cols = 0
    this%tscols = 0
    CALL this%GetOutputlist(Mesh,node,this%cols,this%tscols)

    IF (this%GetRank().EQ.0) CALL this%WriteParaviewFile()
  END SUBROUTINE InitFileIO_vtk

  SUBROUTINE WriteParaviewFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    INTEGER                          :: k
    REAL                             :: ftime
    CHARACTER(LEN=256)               :: basename
#ifdef PARALLEL
    INTEGER                          :: i
#endif
    !------------------------------------------------------------------------!
    ! write a pvd-file: this is a "master" file of all timesteps of all ranks
    CALL this%OpenFile(REPLACE,'pvd')

    ! write vtk header
    WRITE(this%unit+200,FMT='(A)',IOSTAT=this%error_io) &
          '<?xml version="1.0"?>'//LF &
          //'<VTKFile type="Collection" version="0.1" byte_order=' &
          //this%endianness//'>'//LF &
          //REPEAT(' ',2)//'<Collection>'

    ! write entries for each time step
    DO k=0,this%cycles-1
       this%step=k  ! we need this to get the correct file name with time step appended
                    ! ATTENTION: remember to reset the this%step (see below)
       ftime = this%stoptime*k/(this%cycles-1)
#ifdef PARALLEL
       ! in parallel mode each node generates its own file
       DO i=0,this%GetNumProcs()-1
          basename=this%GetBasename(i)
          IF (this%error_io.EQ. 0) WRITE(this%unit+200,FMT='(A,E11.5,A,I4.4,A)',IOSTAT=this%error_io) &
             REPEAT(' ',4)//'<DataSet timestep="',&
             ftime,'" part="', i ,'" file="'//TRIM(basename)//'"/>'
       END DO
#else
       basename=this%GetBasename()
       IF (this%error_io.EQ. 0) WRITE(this%unit+200,FMT='(A,E11.5,A)',IOSTAT=this%error_io) &
          REPEAT(' ',4)//'<DataSet timestep="',ftime,'" part="0" file="' &
                       //TRIM(basename)//'"/>'
#endif
    END DO
    ! reset the step (see comment above)
    this%step=0

    IF (this%error_io.EQ. 0) WRITE(this%unit+200,FMT='(A)',IOSTAT=this%error_io) &
             REPEAT(' ',2)//'</Collection>'//LF//'</VTKFile>'
    IF (this%error_io.GT.0) CALL this%Error("InitFileIO_vtk","cannot write pvd-file")

    CALL this%CloseFile('pvd')
  END SUBROUTINE WriteParaviewFile


  !> \public Determines precision of real numbers in bytes
  !!
  !! Determines the precision (aka size) of a real number.
  !! Single precision (4 bytes), double precision (8 bytes) and
  !! quad precision (16 bytes) are possible results.
  SUBROUTINE GetPrecision(this,realsize)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    INTEGER                          :: realsize !< \param [out] realsize size of real (byte)
    !------------------------------------------------------------------------!
    REAL                             :: real_number
    !------------------------------------------------------------------------!
    INTENT(OUT)                      :: realsize
    !------------------------------------------------------------------------!
    ! size of default real numbers in bytes
    ! Fortran 2008 standard function!
    realsize = STORAGE_SIZE(real_number)/8
  END SUBROUTINE GetPrecision

  !> \public Determines the endianness of the system
  !!
  !! Determines the the endianess of the system (big or little endian)
  SUBROUTINE GetEndianness(this, res, littlestr, bigstr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this      !< \param [in,out] this fileio type
    CHARACTER(LEN=*)                 :: res       !< \param [out] res result string
    CHARACTER(LEN=*)                 :: littlestr !< \param [in] littlestr little endian str
    CHARACTER(LEN=*)                 :: bigstr    !< \param [in] bigstr big endian str
    !------------------------------------------------------------------------!
    INTEGER                          :: k,iTIPO
    CHARACTER, POINTER               :: cTIPO(:)
    !------------------------------------------------------------------------!
    INTENT(IN)                       :: littlestr, bigstr
    INTENT(OUT)                      :: res
    !------------------------------------------------------------------------!

    !endianness
    k = BIT_SIZE(iTIPO)/8
    ALLOCATE(cTIPO(k),STAT = this%error_io)
       IF (this%error_io.NE.0) &
         CALL this%Error("GetEndianness", "Unable to allocate memory.")
    cTIPO(1)='A'
    !cTIPO(2:k-1) = That's of no importance.
    cTIPO(k)='B'

    iTIPO = transfer(cTIPO, iTIPO)
    DEALLOCATE(cTIPO)
    !Test of 'B'=b'01000010' ('A'=b'01000001')
    IF (BTEST(iTIPO,1)) THEN
       write(res,'(A)',IOSTAT=this%error_io)bigstr
    ELSE
       write(res,'(A)',IOSTAT=this%error_io)littlestr
    END IF
  END SUBROUTINE GetEndianness


  !> Creates a list of all data arrays which will be written to file
  !!
  !! Therefore it ignores all arrays with coordinates and checks if the data
  !! arrays are of the dimension of the mesh.
  RECURSIVE SUBROUTINE GetOutputlist(this,Mesh,node,k,l,prefix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this !< \param [in,out] this fileio type
    CLASS(mesh_base),  INTENT(IN)    :: Mesh !< \param [in] mesh mesh type
    TYPE(Dict_TYP),    POINTER       :: node !< \param [in,out] node pointer to (sub-)dict
    INTEGER,           INTENT(INOUT) :: k    !< \param [in,out] k number of data arrays
    INTEGER,           INTENT(INOUT) :: l    !< \param [in,out] l number of data output scalars
    !> \param [in,out] prefix namespace (path) to sub-dict
    CHARACTER(LEN=*), OPTIONAL, INTENT(INOUT) :: prefix
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER           :: dir
    CHARACTER(LEN=MAXKEY)             :: key
    TYPE(real_t)                      :: dummy1
    REAL, DIMENSION(:,:), POINTER     :: dummy2
    REAL, DIMENSION(:,:,:), POINTER   :: dummy3
    REAL, DIMENSION(:,:,:,:), POINTER :: dummy4
    INTEGER, DIMENSION(4)             :: dims
    INTEGER                           :: n
    !------------------------------------------------------------------------!
    DO WHILE(ASSOCIATED(node))
      IF(HasChild(node)) THEN
      ! recursion
        IF(PRESENT(prefix)) THEN
           key = TRIM(prefix)//'/'//TRIM(GetKey(node))
        ELSE
           key = '/'//TRIM(GetKey(node))
        END IF
        dir => GetChild(node)
        CALL this%GetOutputlist(Mesh,dir,k,l,key)
      ELSE IF (HasData(node).AND.&
               .NOT.(GetKey(node).EQ."time")) THEN ! skip time, this is handled in WriteTimestep
         IF (k+1 .GT. MAXCOLS) CALL this%Error("GetOutputlist","reached MAXCOLS")
         ! check shape of output
         dims = 0
         SELECT CASE(GetDatatype(node))
         CASE(DICT_REAL_TWOD)
            CALL GetAttr(node,TRIM(GetKey(node)),dummy2)
            dims(1:2) = SHAPE(dummy2)
            dims(3:4) = 1
         CASE(DICT_REAL_THREED)
            CALL GetAttr(node,TRIM(GetKey(node)),dummy3)
            dims(1:3) = SHAPE(dummy3)
            dims(4) = 1
         CASE(DICT_REAL_FOURD)
            CALL GetAttr(node,TRIM(GetKey(node)),dummy4)
            dims(1:4) = SHAPE(dummy4)
         CASE(DICT_REAL_P)
            CALL GetAttr(node,GetKey(node),dummy1)
            IF (ASSOCIATED(dummy1%p)) THEN
               l=l+1
               this%tsoutput(l)%val => dummy1%p
               this%tsoutput(l)%key = TRIM(GetKey(node))
               IF (PRESENT(prefix)) THEN
                  this%tsoutput(l)%key = TRIM(prefix) // '/' // this%tsoutput(l)%key
               END IF
            END IF
         END SELECT
         ! only register output if it contains mesh data
         ! if not => reset k
         ! TODO: HACKATHON 3D - Dritte Meshdimension abfragen.
         IF ((dims(1).EQ.(Mesh%IMAX-Mesh%IMIN+1)).AND.&
             (dims(2).EQ.(Mesh%JMAX-Mesh%JMIN+1)).AND.&
             (dims(3).EQ.(Mesh%KMAX-Mesh%KMIN+1))) THEN
            ! increase output index
            k = k + 1
            ! store name of the output array
            this%output(k)%key = TRIM(GetKey(node))
            ! prepend prefix if present
            IF (PRESENT(prefix)) THEN
               this%output(k)%key = TRIM(prefix) // '/' // this%output(k)%key
            END IF
            ! allocate memory for pointer to output arrays
            ALLOCATE(this%output(k)%p(dims(4)),STAT=this%error_io)
            IF (this%error_io.NE.0) &
               CALL this%Error( "GetOutputlist", "Unable to allocate memory.")
            ! set pointer to output arrays
            IF (dims(4).EQ.1) THEN
               ! 3D data
               this%output(k)%p(1)%val => dummy3
               this%output(k)%numbytes = SIZE(dummy3)*this%realsize
            ELSE
               ! 4D data (for example pvars - density,xvel,yvel)
               DO n=1,dims(4)
                  this%output(k)%p(n)%val => dummy4(:,:,:,n)
                  this%output(k)%numbytes = SIZE(dummy4)*this%realsize
               END DO
            END IF
         END IF
      END IF
      node=>GetNext(node)
    END DO
  END SUBROUTINE GetOutputlist

  !> \public Specific routine to open a file for vtk I/O
  !!
  SUBROUTINE OpenFile(this,action,ftype)
     IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT)        :: this   !< \param [in,out] this fileio type
    INTEGER,           INTENT(IN)           :: action !< \param [in] action mode of file access
    CHARACTER(LEN=*),  INTENT(IN), OPTIONAL :: ftype  !< \param [in] file type
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)  :: sta,act,pos,fext
    !------------------------------------------------------------------------!
    this%error_io=1
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
       CALL this%Error("OpenFile","Unknown access mode.")
    END SELECT

    ! check file type to open
    IF (PRESENT(ftype)) THEN
       fext=TRIM(ftype)
    ELSE
       ! default
       fext='vts'
    END IF

    SELECT CASE(TRIM(fext))
    CASE('vts')
       ! open the VTS file
       OPEN(this%unit, FILE=this%GetFilename(),STATUS=TRIM(sta), &
            ACCESS = 'STREAM' ,   &
            ACTION=TRIM(act),POSITION=TRIM(pos),IOSTAT=this%error_io)
#ifdef PARALLEL
    CASE('pvts')
       ! open PVTS file (only parallel mode)
       IF (this%GetRank().EQ.0) THEN
          this%extension='pvts' ! change file name extension
          OPEN(this%unit+100, FILE=this%GetFilename(-1),STATUS=TRIM(sta), &
               FORM='FORMATTED',ACTION=TRIM(act),POSITION=TRIM(pos),IOSTAT=this%error_io)
          this%extension='vts' ! reset default extension
       END IF
#endif
    CASE('pvd')
       OPEN(this%unit+200, FILE=TRIM(this%filename)//'.'//TRIM(fext),STATUS=TRIM(sta), &
            FORM='FORMATTED',ACTION=TRIM(act),POSITION=TRIM(pos),IOSTAT=this%error_io)
    CASE DEFAULT
       CALL this%Error("OpenFile","Unknown file type")
    END SELECT

    IF (this%error_io.GT.0) &
       CALL this%Error("OpenFile","cannot open " // TRIM(this%extension) // "-file")
  END SUBROUTINE OpenFile

  !> \public Specific routine to close a file for vtk I/O
  !!
  SUBROUTINE CloseFile(this,ftype)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
    CHARACTER(LEN=*), OPTIONAL       :: ftype !< \param [in] file type
    !------------------------------------------------------------------------!
    CHARACTER(LEN=4)                 :: fext
    !------------------------------------------------------------------------!
    INTENT(IN)                       :: ftype
    !------------------------------------------------------------------------!
    this%error_io=1
    ! check file type to open
    IF (PRESENT(ftype)) THEN
       fext=TRIM(ftype)
    ELSE
       ! default
       fext='vts'
    END IF

    ! close file depending on the file type
    SELECT CASE(TRIM(fext))
    CASE('vts')
       CLOSE(this%unit,IOSTAT=this%error_io)
#ifdef PARALLEL
    CASE('pvts')
        IF (this%GetRank() .EQ. 0 ) THEN
           CLOSE(this%unit+100,IOSTAT=this%error_io)
        END IF
#endif
    CASE('pvd')
        CLOSE(this%unit+200,IOSTAT=this%error_io)
    CASE DEFAULT
       CALL this%Error("OpenFile","Unknown file type")
    END SELECT
    ! set default extension

    IF(this%error_io .NE. 0) &
       CALL this%Error( "CloseFileIO", "cannot close "//TRIM(fext)//"-file")
  END SUBROUTINE CloseFile

!  !> Extract a subkey of a string (obsolete)
!  !!
!  !! Extract a subkey from a string (key) of type "abcd/SUBKEY/efgh"
!  FUNCTION GetSubKey(key) RESULT(subkey)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CHARACTER(LEN=*)      :: key
!    !------------------------------------------------------------------------!
!    CHARACTER(LEN=MAXKEY) :: subkey
!    !------------------------------------------------------------------------!
!    INTENT(IN)            :: key
!    !------------------------------------------------------------------------!
!    ! format of key like: "abcd/SUBKEY/efgh"
!    subkey = key(SCAN(key,"/",.TRUE.)+1:)
!  END FUNCTION GetSubKey

  !> \public Writes XML header to file
  !!
  SUBROUTINE WriteHeader(this,Mesh,Physics,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk),   INTENT(INOUT) :: this    !< \param [in,out] this fileio type
    CLASS(mesh_base),    INTENT(IN)    :: Mesh    !< \param [in] mesh mesh type
    CLASS(physics_base), INTENT(IN)    :: Physics !< \param [in] physics physics type
    TYPE(Dict_TYP),      POINTER       :: IO      !< \param [in,out] IO I/O dictionary
    TYPE(Dict_TYP),      POINTER       :: Header  !< \param [in,out] config config dictionary
    !------------------------------------------------------------------------!
    CALL this%OpenFile(REPLACE)

    this%error_io = 1
    ! write header in vts files
    WRITE(this%linebuf,FMT='(A,6(I7),A)') '<?xml version="1.0"?>' //LF// &
        '<VTKFile type="StructuredGrid" version="0.1" byte_order='&
                   //this%endianness//'>'//LF//&
        '  <StructuredGrid WholeExtent="',&
                  Mesh%IMIN,Mesh%IMAX+this%IP1,Mesh%JMIN,Mesh%JMAX+this%JP1, &
                  Mesh%KMIN,Mesh%KMAX+this%KP1,'">'//LF// &
!                   Mesh%IMIN,Mesh%IMAX+1,Mesh%JMIN,Mesh%JMAX+1,Mesh%KMIN,Mesh%KMAX+1,'">'//LF// &
        '    <FieldData>'//LF// &
        '      <InformationKey type="String" Name="fosite version" format="ascii">'//LF// &
        '        '//TRIM(VERSION)//LF// &
        '      </InformationKey>'//LF
    WRITE(this%unit,IOSTAT=this%error_io) TRIM(this%linebuf)

#ifdef PARALLEL
    ! write header in pvts file (only rank 0)
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%OpenFile(REPLACE,'pvts')
       WRITE(this%unit+100,FMT='(A,6(I7),A)',IOSTAT=this%error_io) '<?xml version="1.0"?>' //LF// &
           '<VTKFile type="PStructuredGrid" version="0.1" byte_order='&
                   //this%endianness//'>'//LF//&
           '  <PStructuredGrid WholeExtent="',&
                   1,Mesh%INUM+this%IP1,1,Mesh%JNUM+this%JP1,1,Mesh%KNUM+this%KP1,'">'//LF//&
           '    <PFieldData>'//LF// &
           '      <InformationKey type="String" Name="fosite version" format="ascii">'//LF// &
           '        '//TRIM(VERSION)//LF//&
           '      </InformationKey>'
       CALL this%CloseFile('pvts')
    END IF
#endif
    IF (this%error_io.GT.0) CALL this%Error("WriteHeader","cannot write (P)VTS file header")

    CALL this%CloseFile()
  END SUBROUTINE WriteHeader

  !> \public Reads the header (not yet implemented)
  !!
  SUBROUTINE ReadHeader(this,success)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this
    LOGICAL                          :: success
    !------------------------------------------------------------------------!
    INTENT(OUT)                      :: success
    !------------------------------------------------------------------------!
    CALL this%OpenFile(READONLY)
    success = .FALSE.
    CALL this%CloseFile()
  END SUBROUTINE ReadHeader

  !> \public Writes the timestep
  !!
  SUBROUTINE WriteTimestamp(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk) :: this
    REAL              :: time
    !------------------------------------------------------------------------!
    INTENT(IN)        :: time
    !------------------------------------------------------------------------!
    CALL this%OpenFile(APPEND)

    WRITE(this%linebuf,FMT='(A,ES14.7,A)') &
        '      <DataArray type="Float64" Name="TIME" NumberOfTuples="1" format="ascii">'//LF// &
                   REPEAT(' ',8),time,LF//&
        '      </DataArray>'//LF
    WRITE(this%unit,IOSTAT=this%error_io) TRIM(this%linebuf)
#ifdef PARALLEL
    ! write time stamp in PVTS file in parallel mode
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%OpenFile(APPEND,'pvts')
       WRITE(this%unit+100,FMT='(A)',ADVANCE='NO',IOSTAT=this%error_io) TRIM(this%linebuf)
       CALL this%CloseFile('pvts')
    END IF

    CALL this%CloseFile()
#endif

  END SUBROUTINE WriteTimestamp

  !> \public Reads the timestep (not yet implemented)
  !!
  SUBROUTINE ReadTimestamp(this,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this
    REAL                             :: time
    !------------------------------------------------------------------------!
    INTENT(OUT)                      :: time
    !------------------------------------------------------------------------!
    CALL this%OpenFile(READONLY)
    time = 0.0
    CALL this%CloseFile()
  END SUBROUTINE ReadTimestamp

  !> \public Writes all desired data arrays to a file
  !!
  SUBROUTINE WriteDataset(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk),    INTENT(INOUT) :: this     !< \param [in,out] this fileio type
    CLASS(mesh_base),     INTENT(IN)    :: Mesh     !< \param [in] mesh mesh type
    CLASS(physics_base),  INTENT(INOUT) :: Physics  !< \param [in] physics physics type
    CLASS(fluxes_base),   INTENT(IN)    :: Fluxes   !< \param [in] fluxes fluxes type
    CLASS(timedisc_base), INTENT(IN)    :: Timedisc !< \param [in] timedisc timedisc type
    TYPE(Dict_TYP),       POINTER       :: Header   !< \param [in,out] IO I/O dictionary
    TYPE(Dict_TYP),       POINTER       :: IO       !< \param [in,out] IO I/O dictionary
    !------------------------------------------------------------------------!
    INTEGER                             :: k,m
    INTEGER                             :: n, offset
#ifdef PARALLEL
    CHARACTER(LEN=256)                  :: basename
    INTEGER                             :: i
#endif
    !------------------------------------------------------------------------!
    ! this part is formerly from fileio_generic.f90
    ! calculate boundary fluxes, if they were requested for write
!    IF(ASSOCIATED(Timedisc%bflux)) THEN
!        print *, "test"
!        DO k=1,4
!            Timedisc%bflux(:,k) = Fluxes%GetBoundaryFlux(Mesh,Physics,k)
!        END DO
!#ifdef PARALLEL
!        CALL MPI_BCAST(Timedisc%bflux, Physics%VNUM*4, DEFAULT_MPI_REAL, 0, &
!                Mesh%comm_cart, ierror)
!#endif
!    END IF

    IF (ASSOCIATED(Timedisc%w)) THEN
      IF (Mesh%FARGO.EQ.3.AND.Mesh%shear_dir.EQ.1) THEN
        CALL Physics%AddBackgroundVelocityX(Mesh,Timedisc%w,Timedisc%pvar,Timedisc%cvar)
      ELSE
        CALL Physics%AddBackgroundVelocityY(Mesh,Timedisc%w,Timedisc%pvar,Timedisc%cvar)
      END IF
    END IF


    ! write the header if either this is the first data set we write or
    ! each data set is written into a new file
    IF ((this%step.EQ.0).OR.(this%cycles.GT.0)) THEN
       CALL this%WriteHeader(Mesh,Physics,Header,IO)
    END IF
    ! write the time stamp
    CALL this%WriteTimestamp(Timedisc%time)
    CALL this%OpenFile(APPEND)

    !------------------------------------------------------------------------!
    ! REMARK: VTS/PVTS data files are opened in unformated or stream mode,
    !         hence we always write the formated data into a line buffer
    !         and then write the line buffer to the data file
    ! -----------------------------------------------------------------------!
    ! 1. META data
    ! -----------------------------------------------------------------------!
    ! a) global information; currently only scalar real numbers supported
#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%OpenFile(APPEND,'pvts')
       IF(this%error_io .NE. 0) CALL this%Error( "WriteDataset", "Can't open pvts file")
    END IF
#endif
    ! loop over all time step data registred in GetOutputlist
    DO k=1,this%TSCOLS
       WRITE(this%linebuf,FMT='(A,ES14.7,A)') &
        '      <DataArray type="Float64" Name="'//TRIM(this%tsoutput(k)%key)// &
                       '" NumberOfTuples="1" format="ascii">'//LF// &
                   REPEAT(' ',8),this%tsoutput(k)%val,LF//&
        '      </DataArray>'//LF
       WRITE(this%unit,IOSTAT=this%error_io) TRIM(this%linebuf)
#ifdef PARALLEL
       ! write time step data in PVTS file in parallel mode
       IF (this%GetRank() .EQ. 0 ) THEN
          WRITE(this%unit+100,FMT='(A)',ADVANCE='NO',IOSTAT=this%error_io) TRIM(this%linebuf)
       END IF
#endif
    END DO

!> \todo : implement generic output routine for non-field/non-scalar data

! simple hack to write positions of the binary star system
!     IF (HasKey(IO, "/sources/grav/binary/binpos/value"))THEN
!       CALL GetAttr(IO, "/sources/grav/binary/binpos", dummy2)
!
!       WRITE(this%tslinebuf,fmt='(ES14.7,A,ES14.7)',IOSTAT=this%error_io)&
!             dummy2(1,1), repeat(' ',4),dummy2(2,1)
!
!       WRITE(this%unit, IOSTAT=this%error_io)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position primary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
!
! #ifdef PARALLEL
!       IF (GetRank(this) .EQ. 0 ) &
!          WRITE(this%unit+100, IOSTAT=this%error_io)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position primary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
!
! #endif
!
!       WRITE(this%tslinebuf,fmt='(ES14.7,A,ES14.7)',IOSTAT=this%error_io)&
!             dummy2(1,2), repeat(' ',4),dummy2(2,2)
!
!       WRITE(this%unit, IOSTAT=this%error_io)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position secondary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
!
! #ifdef PARALLEL
!       IF (GetRank(this) .EQ. 0 ) &
!          WRITE(this%unit+100, IOSTAT=this%error_io)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position secondary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
! #endif
!
!     END IF

    ! -----------------------------------------------------------------------!
    ! b) coordinates
    WRITE(this%linebuf,FMT='(A,(6(I7)),A)')&
         '    </FieldData>'//LF// &
         '    <Piece Extent="', &
                    Mesh%IMIN,Mesh%IMAX+this%IP1,Mesh%JMIN,Mesh%JMAX+this%JP1,&
                    Mesh%KMIN,Mesh%KMAX+this%KP1,'">'//LF// &
         '    <Points>'//LF//&
         '      <DataArray type='//this%realfmt//' NumberOfComponents="3"'//&
                         ' Name="Point" format="appended" offset="0">'//LF// &
         '      </DataArray>'//LF// &
         '    </Points>'//LF// &
         '    <CellData>'//LF
    WRITE(this%unit,IOSTAT=this%error_io) TRIM(this%linebuf)

#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
       WRITE(this%unit+100,FMT='(A)',IOSTAT=this%error_io)&
              repeat(' ',4)//'</PFieldData>' &
        //LF//repeat(' ',4)//'<PPoints>' &
        //LF//repeat(' ',6)//'<DataArray type='//this%realfmt &
                //' NumberOfComponents="3" Name="Point"/>' &
        //LF//repeat(' ',4)//'</PPoints>' &
        //LF//repeat(' ',4)//'<PCellData>'
    END IF
#endif

    ! -----------------------------------------------------------------------!
    ! c) output data fields
    offset = 0
    DO k = 1, this%cols
       WRITE(this%linebuf,FMT='(A,I3,A,I8,A)')&
           REPEAT(' ',6)//'<DataArray type='//this%realfmt//&
                          ' NumberOfComponents="', SIZE(this%output(k)%p), &
              '" Name="'//TRIM(this%output(k)%key)//&
              '" format="appended" offset="',this%offset+offset,'"/>'//LF
       WRITE(this%unit,IOSTAT=this%error_io) TRIM(this%linebuf)

#ifdef PARALLEL
       IF (this%GetRank() .EQ. 0 ) THEN
           WRITE(this%unit+100,FMT='(A,I3,A)',IOSTAT=this%error_io)&
               REPEAT(' ',6)//'<DataArray type='//this%realfmt//&
                              ' NumberOfComponents="', SIZE(this%output(k)%p),&
               '" Name="'//TRIM(this%output(k)%key)//'"/>'
       END IF
#endif
       ! add size of data set (in bytes) + 4 (4 byte integer for size information)
       ! -> "jump address" for next data set in the file in bytes
       offset = offset + this%output(k)%numbytes + 4
    END DO

    WRITE(this%unit,IOSTAT=this%error_io)&
             repeat(' ',6)//'</CellData>' &
       //LF//repeat(' ',4)//'</Piece>' &
       //LF//repeat(' ',2)//'</StructuredGrid>' &
       //LF//repeat(' ',2)//'<GlobalData>'

#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
      WRITE(this%unit+100,FMT='(A)',IOSTAT=this%error_io) REPEAT(' ',4)//'</PCellData>'

      DO i=0,this%GetNumProcs()-1
         basename=this%GetBasename(i)
         WRITE(this%unit+100,FMT='(A)',IOSTAT=this%error_io) &
            REPEAT(' ',4)//'<Piece Extent="'//TRIM(this%extent(i))&
            //'" Source="'//TRIM(basename)//'"/>'
      END DO

      WRITE(this%unit+100,FMT='(A)',IOSTAT=this%error_io)&
        repeat(' ',2)//'</PStructuredGrid>' &
        //LF//'</VTKFile>'

    END IF
#endif

    ! -----------------------------------------------------------------------!
    ! d) boundary fluxes
!     DO i=1,4
!        ! in PARALLEL mode the result is sent to process with rank 0
!        this%bflux(:,i) = GetBoundaryFlux(Fluxes,Mesh,Physics,i)
!
!        WRITE(this%unit,IOSTAT=this%error_io)&
!               LF//repeat(' ',4)//'<DataArray type='//this%realfmt//&
!               'Name="'//TRIM(fluxkey(i))//'" format="ascii" >' &
!               //LF//repeat(' ',6)
!        DO k=1,Physics%VNUM
!            WRITE(this%linebuf((k-1)*20+1:k*20),fmt='(E16.9,A)',IOSTAT=this%error_io)&
!                 this%bflux(k,i),repeat(' ',4)
!        END DO
!
!        WRITE(this%unit,IOSTAT=this%error_io)TRIM(this%linebuf(1:Physics%VNUM*20)) &
!               //LF//repeat(' ',4)//'</DataArray>'
!     END DO

    ! end META data
    ! -----------------------------------------------------------------------------------!

    WRITE(this%unit,IOSTAT=this%error_io)&
       LF//repeat(' ',2)//'</GlobalData>' &
       //LF//'<AppendedData encoding="raw">' &
       //LF//'_'

#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%CloseFile('pvts')
       IF(this%error_io .NE. 0) CALL this%Error( "WriteDataset", "Can't close pvts file")
    END IF
#endif

    ! -----------------------------------------------------------------------!
    ! 2. REAL data
    ! -----------------------------------------------------------------------!
    ! TODO: Hier nochmal rübergucken für 3D
    ! a) coordinates:
    !    (i)  size of data array in bytes (must be a 4 byte integer)
    !    (ii) coordinate data array (generated in InitFileio)
    WRITE(this%unit,IOSTAT=this%error_io) INT(this%offset,4),this%vtkcoords(:,:,:,:)
    IF (this%error_io.GT.0) &
       CALL this%Error( "WriteDataset", "cannot write coordinates")

    ! -----------------------------------------------------------------------!
    ! b) data fields
    !    (i)  size of data array in bytes (must be a 4 byte integer);
    !         this is determined in GetOutputlist and storen in
    !         this%output(:)%numbytes
    !    (ii) the data given by a list of pointers which is also
    !         generated in GetOutputlist
    ! loop over all output data fields
    DO k = 1, this%cols
         ! reorganize the data to comply with VTK style ordering
         n = SIZE(this%output(k)%p)  ! get first dimension
         DO m=1,n
            this%binout(m,:,:,:)=this%output(k)%p(m)%val(:,:,:)
         END DO

         WRITE(this%unit,IOSTAT=this%error_io) this%output(k)%numbytes, &
            this%binout(1:n,:,:,:)

         IF (this%error_io.GT.0) &
            CALL this%Error( "WriteDataset", "cannot write data")
    END DO
    WRITE(this%unit,IOSTAT=this%error_io)LF//repeat(' ',2)//&
         '</AppendedData>'//LF//'</VTKFile>'//LF

    CALL this%CloseFile()
    CALL this%IncTime()
  END SUBROUTINE WriteDataset

  !> \public Closes the file I/O
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this   !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    INTEGER                          :: k
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    DO k=1,this%cols
       IF (ASSOCIATED(this%output(k)%p)) DEALLOCATE(this%output(k)%p)
    END DO
    DEALLOCATE(this%binout,this%vtkcoords,this%output,this%tsoutput,this%bflux)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE fileio_vtk_mod
