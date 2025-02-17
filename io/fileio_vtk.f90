!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_vtk.f90                                                    #
!#                                                                           #
!# Copyright (C) 2010-2024                                                   #
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
  USE fileio_gnuplot_mod
  USE geometry_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE timedisc_base_mod
  USE sources_generic_mod
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
  !>#### some restrictions to shape & amount of output data
  INTEGER, PARAMETER      :: MAXCOMP  = 9   !< max. of allowed components
                                            !! 9 is a tensor (rank 2, dim 3)
  INTEGER, PARAMETER      :: MAXCOLS  = 64  !< max of different output arrays
  CHARACTER, PARAMETER    :: LF = ACHAR(10) !< line feed
  !> \name
  !!#### handling multiple data files with process rank in their names
#ifdef PARALLEL
  INTEGER, PARAMETER      :: RANK_STR_LEN = 5     !< length of multi process rank string
#endif
   !> names of fluxes
!   CHARACTER(LEN=16),DIMENSION(6),PARAMETER  :: fluxkey = (/'bflux_WEST  ', &
!                                                            'bflux_EAST  ', &
!                                                            'bflux_SOUTH ', &
!                                                            'bflux_NORTH ', &
!                                                            'bflux_BOTTOM', &
!                                                            'bflux_TOP   ' /)
  !--------------------------------------------------------------------------!
#ifdef PARALLEL
  !> class for Fortran file handle with process rank in file name
  TYPE, EXTENDS(filehandle_fortran) :: filehandle_vts
    !> \name Variables
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFilehandle
    PROCEDURE :: GetBasename
    PROCEDURE :: GetRankString
!     FINAL :: Finalize_vts
  END TYPE filehandle_vts
#endif
  !> FileIO class for VTK output
  TYPE, EXTENDS(fileio_gnuplot) :: fileio_vtk
    !> \name Variables
    TYPE(filehandle_fortran) :: pvdfile   !< paraview master file
#ifdef PARALLEL
    TYPE(filehandle_fortran) :: pvtsfile  !< parallel vts file
#endif
    CHARACTER(LEN=32)      :: buf         !< buffer for character i/o
    CHARACTER(LEN=32)      :: endianness  !< big/little endian
    CHARACTER(LEN=12)      :: realfmt     !< real number format as string
    INTEGER                :: realsize    !< size of real numbers (bits)
    INTEGER                :: IP1,JP1,KP1
    REAL, DIMENSION(:,:,:,:), POINTER :: &
                              binout      !< binary data output buffer
    REAL, DIMENSION(:,:) , POINTER  :: &
                              bflux       !< boundary flux output buffer
    REAL, DIMENSION(:,:,:,:), POINTER :: &
                              vtkcoords   !< VTK coordinate data
#ifdef PARALLEL
    CHARACTER(LEN=64), DIMENSION(:), POINTER :: &
                              extent      !< extent of all pieces in VTK
#endif
  CONTAINS
    !> \name Methods
    PROCEDURE :: InitFileio_deferred => InitFileio_vtk
    PROCEDURE :: WriteHeader
!     PROCEDURE :: ReadHeader
    PROCEDURE :: WriteDataset_deferred => WriteDataset_vtk
    PROCEDURE :: WriteParaviewFile
    FINAL :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    fileio_vtk


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
    CLASS(sources_list), ALLOCATABLE, INTENT(IN) :: Sources !< \param [in] Sources sources type
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: config   !< \param [in] IO Dictionary for I/O
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: IO       !< \param [in] IO Dictionary for I/O
    !------------------------------------------------------------------------!
    TYPE(Output_TYP),DIMENSION(:),POINTER    :: poutput
    TYPE(real_t)                             :: dummy1
    TYPE(Dict_TYP), POINTER                  :: node
    CHARACTER(LEN=MAX_CHAR_LEN),DIMENSION(2) :: skip
    INTEGER                             :: k,i,j
#ifdef PARALLEL
    INTEGER, DIMENSION(:), POINTER      :: sendbuf,recvbuf
    INTEGER                             :: n
#endif
    REAL, DIMENSION(:,:,:,:,:), POINTER :: corners
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    ALLOCATE(filehandle_vts::this%datafile)
#endif
    ! start init form base class in beginning
    CALL this%InitFileIO(Mesh,Physics,Timedisc,Sources,config,IO,"VTK","vts",textfile=.FALSE.)

    ! some sanity checks
    IF (this%datafile%cycles.NE.this%count+1) &
      CALL this%Error('fileio_vtk::InitFileIO','filecycles = count+1 required for VTK output')

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

    this%MAXCOLS = MAXCOLS

    ! allocate memory for auxilliary data fields
    ALLOCATE(this%binout(MAXCOMP,                &
                         Mesh%IMIN:Mesh%IMAX,    &
                         Mesh%JMIN:Mesh%JMAX,    &
                         Mesh%KMIN:Mesh%KMAX),   &
             this%vtkcoords(3,                   &
                         Mesh%IMIN:Mesh%IMAX+this%IP1,  &
                         Mesh%JMIN:Mesh%JMAX+this%JP1,  &
                         Mesh%KMIN:Mesh%KMAX+this%KP1), &
             this%output(this%MAXCOLS),               &
             this%tsoutput(this%MAXCOLS),             &
             this%bflux(Physics%VNUM,6),         &
             STAT = this%err)
    IF (this%err.NE.0) &
         CALL this%Error("fileio_vtk::InitFileio_vtk", "Unable to allocate memory.")

    ! determine size in bits of default real numbers
    ! Fortran 2008 standard function!
    this%realsize = STORAGE_SIZE(this%binout)
    SELECT CASE(this%realsize)
    CASE(32,64,128)
      ! generate string with default real number format (size in bits)
      WRITE(this%linebuf,FMT='(I4)',IOSTAT=this%err) this%realsize
      WRITE(this%realfmt,FMT='(A)',IOSTAT=this%err) '"Float' // TRIM(AdjustL(this%linebuf)) // '"'
      ! size of default real numbers in byte
      this%realsize = this%realsize/8
    CASE DEFAULT
      CALL this%Error("fileio_vtk::InitFileio_vtk", &
        "only single, double or quadruple precision real numbers (4,8,16 bytes) supported")
    END SELECT

    ! determine the system endianness
    CALL this%GetEndianness(this%endianness,'"LittleEndian"','"BigEndian"')

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
       ALLOCATE(this%extent(0:(n-1)),sendbuf(6),recvbuf(6*n),STAT = this%err)
    ELSE
       ! see comment above
       ALLOCATE(sendbuf(6),recvbuf(1),STAT = this%err)
    END IF
    IF (this%err.NE.0) &
       CALL this%Error( "fileio_vtk::InitFileio_vtk", "Unable to allocate memory.")

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
               0, MPI_COMM_WORLD, this%err)

    ! store information about grid extent at the rank 0 node
    IF (this%GetRank() .EQ. 0 ) THEN
       DO i=0,n-1
          WRITE(this%extent(i),FMT='(6(I7))',IOSTAT=this%err)&
          recvbuf(i*6+1),recvbuf(i*6+2),recvbuf(i*6+3),recvbuf(i*6+4),recvbuf(i*6+5),recvbuf(i*6+6)
       END DO
    END IF
    ! free buffer memory
    DEALLOCATE(sendbuf,recvbuf,STAT=this%err)
    IF (this%err.NE.0) &
       CALL this%Error( "fileio_vtk::InitFileio_vtk", "Deallocation of sendbuf/recvbuf failed.")
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

    ! get pointer to simulation time to make it the first entry
    ! in the time step data set
    CALL GetAttr(IO,"/timedisc/time",dummy1)
    IF (ASSOCIATED(dummy1%p)) THEN
       this%tsoutput(1)%val => dummy1%p
       this%tsoutput(1)%key = "TIME"
       this%TSCOLS = 1
    ELSE
       this%TSCOLS = 0
    END IF

    ! create list of output data arrays
    node => IO
    k = 0
    skip(1:2) = [CHARACTER(LEN=MAX_CHAR_LEN) :: "corners", "time"]
    CALL this%GetOutputlist(Mesh,node,k,this%TSCOLS,skip)

    ! shrink this%output
    poutput => this%output
    ALLOCATE(this%output,SOURCE=poutput(1:k),STAT=this%err)
    IF (this%err.NE.0) &
      CALL this%Error("fileio_gnuplot::InitFileIO","memory allocation failed for this%output")
    DEALLOCATE(poutput)

    ! set size of output arrays in bytes
    DO k = 1, SIZE(this%output)
       this%output(k)%numbytes = 0
       DO i = 1, SIZE(this%output(k)%p)
          this%output(k)%numbytes = this%output(k)%numbytes + SIZE(this%output(k)%p(i)%val)*this%realsize
       END DO
    END DO

    ! create file handle for the pvd file; only on rank 0 in parallel mode
#ifdef PARALLEL
    IF (this%GetRank().EQ.0) THEN
      CALL this%pvtsfile%InitFilehandle(this%datafile%filename,this%datafile%path,"pvts",&
        textfile=.TRUE.,onefile=.FALSE.,cycles=this%count+1)
#endif
      CALL this%pvdfile%InitFilehandle(this%datafile%filename,this%datafile%path,"pvd",&
        textfile=.TRUE.,onefile=.TRUE.,cycles=1)
#ifdef PARALLEL
    END IF
#endif
  END SUBROUTINE InitFileIO_vtk

  !> \public Write the Paraview global description file
  !!
  SUBROUTINE WriteParaviewFile(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_vtk), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    INTEGER                          :: k
    REAL                             :: ftime
#ifdef PARALLEL
    INTEGER                          :: i
#endif
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (this%GetRank().EQ.0) THEN
#endif
      CALL this%pvdfile%OpenFile(REPLACE,this%step)

      ! write vtk header
      WRITE(UNIT=this%pvdfile%GetUnitNumber(),FMT='(A)',IOSTAT=this%err) &
            '<?xml version="1.0"?>'//LF &
            //'<VTKFile type="Collection" version="0.1" byte_order=' &
            //TRIM(this%endianness)//'>'//LF &
            //REPEAT(' ',2)//'<Collection>'

      ! write entries for each time step
      DO k=0,this%step
         ftime = this%stoptime*k/this%count
         IF (this%err.EQ. 0) &
            WRITE(UNIT=this%pvdfile%GetUnitNumber(),FMT='(A,ES12.5,A)',IOSTAT=this%err) &
                  REPEAT(' ',4)//'<DataSet timestep="',ftime,'" part="0" file="' // &
#ifdef PARALLEL
                       TRIM(this%pvtsfile%GetBasename(k))//&
                       "." // TRIM(this%pvtsfile%extension)//'"/>'
#else
                       TRIM(this%datafile%GetBasename(k))//&
                       "." // TRIM(this%datafile%extension)//'"/>'
#endif
#ifdef PARALLEL
         ! in parallel mode each node generates its own file, but the rank 0 node
         ! creates the pvd and pvts files
!          DO i=0,this%GetNumProcs()-1
!             IF (this%err.EQ. 0) WRITE(this%pvtsfile%GetUnitNumber(),FMT='(A,E11.5,A,I4.4,A)',IOSTAT=this%err) &
!                REPEAT(' ',4)//'<DataSet timestep="',&
!                ftime,'" part="', i ,'" file="'//TRIM(this%pvtsfile%GetBasename(i))//'"/>'
!          END DO
#else
#endif
      END DO

      IF (this%err.EQ. 0) &
         WRITE(UNIT=this%pvdfile%GetUnitNumber(),FMT='(A)',IOSTAT=this%err) &
               REPEAT(' ',2)//'</Collection>'//LF//'</VTKFile>'

      CALL this%pvdfile%CloseFile(this%step)

      IF (this%err.NE.0) CALL this%Error("fileio_vtk::WriteParaviewFile","cannot write pvd-file")

#ifdef PARALLEL
    END IF
#endif
  END SUBROUTINE WriteParaviewFile


!   !> \public Specific routine to open a file for vtk I/O
!   !!
!   SUBROUTINE OpenFile(this,action,ftype)
!      IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(fileio_vtk), INTENT(INOUT)        :: this   !< \param [in,out] this fileio type
!     INTEGER,           INTENT(IN)           :: action !< \param [in] action mode of file access
!     CHARACTER(LEN=*),  INTENT(IN), OPTIONAL :: ftype  !< \param [in] file type
!     !------------------------------------------------------------------------!
!     CHARACTER(LEN=32)  :: sta,act,pos,fext
!     !------------------------------------------------------------------------!
!     this%err=1
!     SELECT CASE(action)
!     CASE(READONLY)
!        sta = 'OLD'
!        act = 'READ'
!        pos = 'REWIND'
!     CASE(READEND)
!        sta = 'OLD'
!        act = 'READ'
!        pos = 'APPEND'
!     CASE(REPLACE)
!        sta = 'REPLACE'
!        act = 'WRITE'
!        pos = 'REWIND'
!     CASE(APPEND)
!        sta = 'OLD'
!        act = 'READWRITE'
!        pos = 'APPEND'
!     CASE DEFAULT
!        CALL this%Error("OpenFile","Unknown access mode.")
!     END SELECT
! 
!     ! check file type to open
!     IF (PRESENT(ftype)) THEN
!        fext=TRIM(ftype)
!     ELSE
!        ! default
!        fext='vts'
!     END IF
! 
!     SELECT CASE(TRIM(fext))
!     CASE('vts')
!        ! open the VTS file
!        OPEN(this%unit, FILE=this%GetFilename(),STATUS=TRIM(sta), &
!             ACCESS = 'STREAM' ,   &
!             ACTION=TRIM(act),POSITION=TRIM(pos),IOSTAT=this%err)
! #ifdef PARALLEL
!     CASE('pvts')
!        ! open PVTS file (only parallel mode)
!        IF (this%GetRank().EQ.0) THEN
!           this%extension='pvts' ! change file name extension
!           OPEN(this%unit+100, FILE=this%GetFilename(-1),STATUS=TRIM(sta), &
!                FORM='FORMATTED',ACTION=TRIM(act),POSITION=TRIM(pos),IOSTAT=this%err)
!           this%extension='vts' ! reset default extension
!        END IF
! #endif
!     CASE('pvd')
!        OPEN(this%unit+200, FILE=TRIM(this%filename)//'.'//TRIM(fext),STATUS=TRIM(sta), &
!             FORM='FORMATTED',ACTION=TRIM(act),POSITION=TRIM(pos),IOSTAT=this%err)
!     CASE DEFAULT
!        CALL this%Error("OpenFile","Unknown file type")
!     END SELECT
! 
!     IF (this%err.GT.0) &
!        CALL this%Error("OpenFile","cannot open " // TRIM(this%extension) // "-file")
!   END SUBROUTINE OpenFile

!   !> \public Specific routine to close a file for vtk I/O
!   !!
!   SUBROUTINE CloseFile(this,ftype)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(fileio_vtk), INTENT(INOUT) :: this  !< \param [in,out] this fileio type
!     CHARACTER(LEN=*), OPTIONAL       :: ftype !< \param [in] file type
!     !------------------------------------------------------------------------!
!     CHARACTER(LEN=4)                 :: fext
!     !------------------------------------------------------------------------!
!     INTENT(IN)                       :: ftype
!     !------------------------------------------------------------------------!
!     this%err=1
!     ! check file type to open
!     IF (PRESENT(ftype)) THEN
!        fext=TRIM(ftype)
!     ELSE
!        ! default
!        fext='vts'
!     END IF
! 
!     ! close file depending on the file type
!     SELECT CASE(TRIM(fext))
!     CASE('vts')
!        CLOSE(this%unit,IOSTAT=this%err)
! #ifdef PARALLEL
!     CASE('pvts')
!         IF (this%GetRank() .EQ. 0 ) THEN
!            CLOSE(this%unit+100,IOSTAT=this%err)
!         END IF
! #endif
!     CASE('pvd')
!         CLOSE(this%unit+200,IOSTAT=this%err)
!     CASE DEFAULT
!        CALL this%Error("OpenFile","Unknown file type")
!     END SELECT
!     ! set default extension
! 
!     IF(this%err .NE. 0) &
!        CALL this%Error( "CloseFileIO", "cannot close "//TRIM(fext)//"-file")
!   END SUBROUTINE CloseFile

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
    INTEGER :: k
    !------------------------------------------------------------------------!
    ! write header in vts files
    WRITE(this%linebuf,FMT='(A,6(I7),A)') '<?xml version="1.0"?>' //LF// &
        '<VTKFile type="StructuredGrid" version="0.1" byte_order='&
                   //TRIM(this%endianness)//'>'//LF//&
        '  <StructuredGrid WholeExtent="',&
                  Mesh%IMIN,Mesh%IMAX+this%IP1,Mesh%JMIN,Mesh%JMAX+this%JP1, &
                  Mesh%KMIN,Mesh%KMAX+this%KP1,'">'//LF// &
        '    <FieldData>'//LF// &
        '      <InformationKey type="String" Name="fosite version" format="ascii">'//LF// &
        '        '//TRIM(VERSION)//LF// &
        '      </InformationKey>'//LF
    WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) TRIM(this%linebuf)

#ifdef PARALLEL
    ! write header in pvts file (only rank 0)
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%pvtsfile%OpenFile(REPLACE,this%step)
       WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A,6(I7),A)',IOSTAT=this%err) '<?xml version="1.0"?>' //LF// &
           '<VTKFile type="PStructuredGrid" version="0.1" byte_order='&
                   //TRIM(this%endianness)//'>'//LF//&
           '  <PStructuredGrid WholeExtent="',&
                   1,Mesh%INUM+this%IP1,1,Mesh%JNUM+this%JP1,1,Mesh%KNUM+this%KP1,'">'//LF//&
           '    <PFieldData>'//LF// &
           '      <InformationKey type="String" Name="fosite version" format="ascii">'//LF// &
           '        '//TRIM(VERSION)//LF//&
           '      </InformationKey>'
       IF (this%err.NE.0) &
          CALL this%Error("fileio_vtk::WriteHeader","Writing pvts file header failed")
    END IF
#endif

    ! loop over all time step data registred in GetOutputlist
    DO k=1,this%TSCOLS
       WRITE(this%linebuf,FMT='(A,ES14.7,A)') &
        '      <DataArray type='//TRIM(this%realfmt)//' Name="'//TRIM(this%tsoutput(k)%key)// &
                       '" NumberOfTuples="1" format="ascii">'//LF// &
                   REPEAT(' ',8),this%tsoutput(k)%val,LF//&
        '      </DataArray>'
       IF (this%err.EQ.0) THEN
#ifndef PARALLEL
          ! write time step data to VTS file in serial mode
          WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) TRIM(this%linebuf) // LF
#else
          ! write time step data to PVTS file in parallel mode
          IF (this%GetRank() .EQ. 0 ) THEN
             WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A)',ADVANCE='NO',IOSTAT=this%err) TRIM(this%linebuf)
             CALL this%pvtsfile%CloseFile(this%step)
          END IF
#endif
       END IF
    END DO

    IF (this%err.NE.0) CALL this%Error("fileio_vtk::WriteHeader","cannot write (P)VTS file header")
  END SUBROUTINE WriteHeader

  !> \public Reads the header (not yet implemented)
  !!
!   SUBROUTINE ReadHeader(this,success)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     CLASS(fileio_vtk), INTENT(INOUT) :: this
!     LOGICAL                          :: success
!     !------------------------------------------------------------------------!
!     INTENT(OUT)                      :: success
!     !------------------------------------------------------------------------!
!     CALL this%OpenFile(READONLY,this%step)
!     success = .FALSE.
!     CALL this%CloseFile(this%step)
!   END SUBROUTINE ReadHeader

  !> \public Writes all desired data arrays to a file
  !!
  SUBROUTINE WriteDataset_vtk(this,Mesh,Physics,Fluxes,Timedisc,Header,IO)
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
    INTEGER                             :: i,j,k,m,l
    INTEGER                             :: n, offset
#ifdef PARALLEL
    CHARACTER(LEN=256)                  :: filename
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

    CALL this%WriteParaviewFile()

    !------------------------------------------------------------------------!
    ! REMARK: VTS data files are opened in unformated or stream mode,
    !         hence we always write the formated data into a line buffer
    !         and then write the line buffer to the data file.
    !         PVTS and PVD files are opened in formated mode.
    ! -----------------------------------------------------------------------!
    ! 1. META data
    ! -----------------------------------------------------------------------!
    ! a) global information; currently only scalar real numbers supported
#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%pvtsfile%OpenFile(APPEND,this%step)
       IF(this%err .NE. 0) CALL this%Error( "fileio_vtk::WriteDataset_vtk", "Can't open pvts file")
    END IF
#endif

!> \todo : implement generic output routine for non-field/non-scalar data

! simple hack to write positions of the binary star system
!     IF (HasKey(IO, "/sources/grav/binary/binpos/value"))THEN
!       CALL GetAttr(IO, "/sources/grav/binary/binpos", dummy2)
!
!       WRITE(this%tslinebuf,fmt='(ES14.7,A,ES14.7)',IOSTAT=this%err)&
!             dummy2(1,1), repeat(' ',4),dummy2(2,1)
!
!       WRITE(this%unit, IOSTAT=this%err)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position primary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
!
! #ifdef PARALLEL
!       IF (GetRank(this) .EQ. 0 ) &
!          WRITE(this%unit+100, IOSTAT=this%err)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position primary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
!
! #endif
!
!       WRITE(this%tslinebuf,fmt='(ES14.7,A,ES14.7)',IOSTAT=this%err)&
!             dummy2(1,2), repeat(' ',4),dummy2(2,2)
!
!       WRITE(this%unit, IOSTAT=this%err)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position secondary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
!
! #ifdef PARALLEL
!       IF (GetRank(this) .EQ. 0 ) &
!          WRITE(this%unit+100, IOSTAT=this%err)&
!              LF//repeat(' ',6)//'<DataArray type="Float64" Name="position secondary"' &
!                            //' NumberOfTuples="2" format="ascii">' &
!            //LF//repeat(' ',8)//TRIM(this%tslinebuf) &
!            //LF//repeat(' ',6)//'</DataArray>'
! #endif
!
!     END IF

    ! -----------------------------------------------------------------------!
    ! b) coordinates
    offset = 0
    WRITE(this%linebuf,FMT='(A,(6(I7)),A,I8,A)')&
         '    </FieldData>'//LF// &
         '    <Piece Extent="', &
                    Mesh%IMIN,Mesh%IMAX+this%IP1,Mesh%JMIN,Mesh%JMAX+this%JP1,&
                    Mesh%KMIN,Mesh%KMAX+this%KP1,'">'//LF// &
         '      <Points>'//LF//&
         '        <DataArray type='//TRIM(this%realfmt)//' NumberOfComponents="3"'//&
                         ' Name="Point" format="appended" offset="',offset,'">'//LF// &
         '        </DataArray>'//LF// &
         '      </Points>'//LF// &
         '      <CellData>'//LF
    IF (this%err.EQ.0) &
      WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) TRIM(this%linebuf)
    ! add size of data set (in bytes) + 4 (4 byte integer for size information)
    ! -> "jump address" for next data set in the file in bytes
    offset = offset + SIZE(this%vtkcoords)*this%realsize + 4

#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
       WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A)',IOSTAT=this%err)&
              repeat(' ',4)//'</PFieldData>' &
        //LF//repeat(' ',4)//'<PPoints>' &
        //LF//repeat(' ',6)//'<DataArray type='//TRIM(this%realfmt) &
                //' NumberOfComponents="3" Name="Point"/>' &
        //LF//repeat(' ',4)//'</PPoints>' &
        //LF//repeat(' ',4)//'<PCellData>'
    END IF
#endif

    ! -----------------------------------------------------------------------!
    ! c) output data fields
    DO k = 1, SIZE(this%output)
       WRITE(this%linebuf,FMT='(A,I3,A,I8,A)')&
           REPEAT(' ',8)//'<DataArray type='//this%realfmt//&
                          ' NumberOfComponents="', SIZE(this%output(k)%p), &
              '" Name="'//TRIM(this%output(k)%key)//&
              '" format="appended" offset="',offset,'"/>'//LF
    IF (this%err.EQ.0) &
      WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) TRIM(this%linebuf)

#ifdef PARALLEL
       IF (this%GetRank() .EQ. 0 ) THEN
           WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A,I3,A)',IOSTAT=this%err)&
               REPEAT(' ',6)//'<DataArray type='//this%realfmt//&
                              ' NumberOfComponents="', SIZE(this%output(k)%p),&
               '" Name="'//TRIM(this%output(k)%key)//'"/>'
       END IF
#endif
       ! add size of data set (in bytes) + 4 (4 byte integer for size information)
       ! -> "jump address" for next data set in the file in bytes
       offset = offset + this%output(k)%numbytes + 4
    END DO

    IF (this%err.EQ.0) &
      WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) &
             repeat(' ',6)//'</CellData>' &
       //LF//repeat(' ',4)//'</Piece>' &
       //LF//repeat(' ',2)//'</StructuredGrid>' &
       //LF//repeat(' ',2)//'<GlobalData>'

#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
      WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A)',IOSTAT=this%err) REPEAT(' ',4)//'</PCellData>'

      SELECT TYPE(df=>this%datafile)
      CLASS IS(filehandle_vts)
         DO i=0,this%GetNumProcs()-1
            filename = TRIM(df%path) // &
                     TRIM(df%filename) // &
                     TRIM(df%GetRankString(i)) // &
                     TRIM(df%GetStepString(this%step)) // &
                     "." // TRIM(this%datafile%extension)
            WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A)',IOSTAT=this%err) &
               REPEAT(' ',4)//'<Piece Extent="'//TRIM(this%extent(i))&
               //'" Source="'//TRIM(filename)//'"/>'
         END DO
      CLASS DEFAULT
        CALL this%Error("fileio_vtk::WriteDataset","datafile must be of type filehandle_vts in parallel mode")
      END SELECT

      WRITE(UNIT=this%pvtsfile%GetUnitNumber(),FMT='(A)',IOSTAT=this%err)&
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
!        WRITE(this%unit,IOSTAT=this%err)&
!               LF//repeat(' ',4)//'<DataArray type='//this%realfmt//&
!               'Name="'//TRIM(fluxkey(i))//'" format="ascii" >' &
!               //LF//repeat(' ',6)
!        DO k=1,Physics%VNUM
!            WRITE(this%linebuf((k-1)*20+1:k*20),fmt='(E16.9,A)',IOSTAT=this%err)&
!                 this%bflux(k,i),repeat(' ',4)
!        END DO
!
!        WRITE(this%unit,IOSTAT=this%err)TRIM(this%linebuf(1:Physics%VNUM*20)) &
!               //LF//repeat(' ',4)//'</DataArray>'
!     END DO

    ! end META data
    ! -----------------------------------------------------------------------------------!

    IF (this%err.EQ.0) &
      WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err)&
        LF//repeat(' ',2)//'</GlobalData>' &
        //LF//'  <AppendedData encoding="raw">' &
        //LF//'_'
    IF (this%err.NE.0) &
       CALL this%Error("fileio_vtk::WriteDataset", "writing data array specs failed")

#ifdef PARALLEL
    IF (this%GetRank() .EQ. 0 ) THEN
       CALL this%pvtsfile%CloseFile(this%step)
       IF(this%err .NE. 0) CALL this%Error( "fileio_vtk::WriteDataset_vtk", "Can't close pvts file")
    END IF
#endif

    ! -----------------------------------------------------------------------!
    ! 2. REAL data
    ! -----------------------------------------------------------------------!
    ! a) coordinates:
    !    (i)  size of data array in bytes (must be a 4 byte integer)
    !    (ii) coordinate data array (generated in InitFileio)
    IF (this%err.EQ.0) &
       WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) &
             INT(SIZE(this%vtkcoords)*this%realsize + 4,KIND=4),this%vtkcoords(:,:,:,:)
    IF (this%err.NE.0) &
       CALL this%Error("fileio_vtk::WriteDataset", "writing coordinates failed")

    ! -----------------------------------------------------------------------!
    ! b) data fields
    !    (i)  size of data array in bytes (must be a 4 byte integer);
    !         this is determined in GetOutputlist and stored in
    !         this%output(:)%numbytes
    !    (ii) the data given by a list of pointers which is also
    !         generated in GetOutputlist
    ! loop over all output data fields
    DO l = 1, SIZE(this%output)
         ! reorganize the data to comply with VTK style ordering
         n = SIZE(this%output(l)%p)
         DO m=1,n
           DO k=1,Mesh%KMAX-Mesh%KMIN+1
             DO j=1,Mesh%JMAX-Mesh%JMIN+1
               DO i=1,Mesh%IMAX-Mesh%IMIN+1
                 this%binout(m,Mesh%IMIN+i-1,Mesh%JMIN+j-1,Mesh%KMIN+k-1) &
                   = this%output(l)%p(m)%val(i,j,k)
               END DO
             END DO
           END DO
         END DO
         IF (this%err.EQ.0) &
            WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) &
                INT(this%output(l)%numbytes,KIND=4), &
!               this%binout(1:n,:,:,:)
                ((((this%binout(m,i,j,k),m=1,n),i=Mesh%IMIN,Mesh%IMAX), &
                    j=Mesh%JMIN,Mesh%JMAX),k=Mesh%KMIN,Mesh%KMAX)

         IF (this%err.NE.0) &
            CALL this%Error("fileio_vtk::WriteDataset", "writing data")
    END DO
    WRITE(UNIT=this%datafile%GetUnitNumber(),IOSTAT=this%err) LF//repeat(' ',2)//&
         '</AppendedData>'//LF//'</VTKFile>'//LF
  END SUBROUTINE WriteDataset_vtk

  !> \public Destructor of VTK file I/O class
  !!
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fileio_vtk), INTENT(INOUT) :: this   !< \param [in,out] this fileio type
    !------------------------------------------------------------------------!
    ! REMARK: this%output and this%tsoutput are deallocated in the gnuplot finalizer
    !         which is called automatically
    IF (ASSOCIATED(this%binout)) DEALLOCATE(this%binout)
    IF (ASSOCIATED(this%vtkcoords)) DEALLOCATE(this%vtkcoords)
    IF (ASSOCIATED(this%bflux)) DEALLOCATE(this%bflux)
    NULLIFY(this%binout,this%vtkcoords,this%bflux)
  END SUBROUTINE Finalize

#ifdef PARALLEL
  !> basic initialization of Fortran file handle with extension for parallel vts files
  SUBROUTINE InitFilehandle(this,filename,path,extension,textfile,onefile,cycles,unit)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(filehandle_vts), INTENT(INOUT) :: this !< \param [inout] this file handle class
    CHARACTER(LEN=*), INTENT(IN) :: filename !< \param [in] filename file name without extension
    CHARACTER(LEN=*), INTENT(IN) :: path     !< \param [in] path file path without filename
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: extension !< \param [in] extension file name extension
    LOGICAL, OPTIONAL, INTENT(IN)      :: textfile !< \param [in] textfile true for text data
    LOGICAL, OPTIONAL, INTENT(IN)      :: onefile  !< \param [in] onefile true if all data goes into one file
    INTEGER, OPTIONAL, INTENT(IN)      :: cycles   !< \parma [in] cycles max number of files
    INTEGER, OPTIONAL, INTENT(IN)      :: unit     !< \parma [in] unit force fortran i/o unit number
    !-------------------------------------------------------------------!
    CALL this%filehandle_fortran%InitFilehandle(filename,path,extension,textfile,onefile,cycles,unit)
    ! check number of parallel processes
    IF (RANK_STR_LEN.LT.3) &
      CALL this%Error("filehandle_vts::InitFilehandle", &
        "we need at least 3 bytes for rank string: 2 for '-r' + 1 for process number (0..9)")
    IF (this%GetNumProcs().GT.10**(RANK_STR_LEN-3)) &
      CALL this%Error("filehandle_vts::InitFilehandle", &
        "number of processes for multiple file output exceeds limits")
  END SUBROUTINE InitFilehandle

  !> \public get file name of Fortran stream including rank and time step
  FUNCTION GetBasename(this,step) RESULT (fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_vts), INTENT(IN)     :: this   !< \param [in] this file handle
    INTEGER, INTENT(IN)                   :: step   !< \param [in] step time step
    CHARACTER(LEN=256)                    :: fname
    !------------------------------------------------------------------------!
    fname = TRIM(this%path) // TRIM(this%filename) // TRIM(this%GetRankString()) &
      // TRIM(this%GetStepString(step))
  END FUNCTION GetBasename

  !> \public convert process rank to string
  FUNCTION GetRankString(this,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(filehandle_vts), INTENT(IN) :: this   !< \param [in] this file handle
    INTEGER, OPTIONAL, INTENT(IN)     :: rank   !< \param [in] rank process rank
    CHARACTER(LEN=RANK_STR_LEN)       :: GetRankString
    !------------------------------------------------------------------------!
    CHARACTER(LEN=16)                 :: fmtstr
    !-------------------------------------------------------------------!
    ! process rank as string with "-r" + leading zeros
    WRITE (fmtstr ,'(A,I1,A)') "(A2,I0.",RANK_STR_LEN-2,")"
    IF (PRESENT(rank)) THEN
      WRITE (GetRankString,FMT=TRIM(fmtstr)) "-r",rank
    ELSE
      WRITE (GetRankString,FMT=TRIM(fmtstr)) "-r",this%GetRank()
    END IF
  END FUNCTION GetRankString

!   !> \public Destructor of VTS file handle
!   !!
!   SUBROUTINE Finalize_vts(this)
!     IMPLICIT NONE
!     !------------------------------------------------------------------------!
!     TYPE(filehandle_vts), INTENT(INOUT) :: this   !< \param [in,out] this filehandle
!     !------------------------------------------------------------------------!
!   END SUBROUTINE Finalize_vts
#endif
END MODULE fileio_vtk_mod
