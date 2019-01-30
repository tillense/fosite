!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: restart.f90                                                       #
!#                                                                           #
!# Copyright (C) 2015-2019                                                   #
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
!> Restart fosite using a binary output file
!!
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \example restart.f90
!!
!! The file uses the binary input and copies its config to the config of
!! fosite (corresponds to Makeconfig in other intializations). The data of
!! the file is afterwards used for initialization of the fields (corresponds
!! to Initdata).
!!
!! \warning restart module does not work with farfield boundaries (tested on gresho)
!!
!! \todo At the moment this file inhibits all the source code in order to read
!!       in a file properly. However, this code should later be moved to the
!!       ReadHeader and ReadData routines in fileio_binary.f90.
!!
!! Usage:
!! ./restart filename
!----------------------------------------------------------------------------!
PROGRAM Restart
  USE fosite_mod
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
  CLASS(fosite), ALLOCATABLE :: Sim
  TYPE(Dict_TYP), POINTER :: input
  REAL                    :: time
  INTEGER                 :: i, step
  CHARACTER(LEN=100)      :: filename, filename_tmp
#ifdef PARALLEL
  INTEGER, DIMENSION(3)   :: decomposition
#endif
  LOGICAL                 :: file_exists
  !--------------------------------------------------------------------------!

  ! load file
  CALL GET_COMMAND_ARGUMENT(1, filename)
  INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
  IF(.NOT.file_exists) THEN
    print *, "restart, ", "Input file does not exist!"
    STOP
  END IF

  ALLOCATE(Sim)

  CALL Sim%InitFosite()

  !-------------------- load configuration-----------------------------------!
  input => LoadConfig(TRIM(filename))
  CALL GetAttr(input, "config", Sim%config)

  ! set starting time in fosite
  CALL GetAttr(input, '/timedisc/time', time)
  time = 0.0
  CALL SetAttr(Sim%config,"/timedisc/starttime", time)

#ifdef PARALLEL
  decomposition(1) = 1
  decomposition(2) = -1
  decomposition(3) = 1
  CALL SetAttr(Sim%config, "/mesh/decomposition", decomposition)
#endif

  ! get step number by stripping from filename
  filename_tmp = filename
  i = SCAN(filename_tmp,'_',back=.TRUE.)
  IF(i.GT.0) THEN
    filename_tmp = filename_tmp(i+1:i+4)
    IF(filename_tmp(1:1) .EQ. "0") filename_tmp = filename_tmp(2:)
    IF(filename_tmp(1:1) .EQ. "0") filename_tmp = filename_tmp(2:)
    IF(filename_tmp(1:1) .EQ. "0") filename_tmp = filename_tmp(2:)
    READ(filename_tmp,*) step
    CALL SetAttr(Sim%config,"/datafile/step", step)
  ELSE
    CALL Sim%Error("Restart", "Unable to find starting step in filename_tmp.")
  END IF


  !-------------------- setup & data initialization -------------------------!
  CALL Sim%Setup()

  CALL LoadData(Sim,TRIM(filename))

  CALL Sim%Info(" DATA-----> initial condition: restarted from data file - " // &
  TRIM(filename))

  IF(HasKey(input,"/timedisc/xmomentum")) THEN
    CALL Sim%Physics%Convert2Primitive(Sim%Timedisc%cvar,Sim%Timedisc%pvar)
  ELSE
    CALL Sim%Physics%Convert2Conservative(Sim%Timedisc%pvar,Sim%Timedisc%cvar)
  END IF

  CALL Sim%Run()
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

CONTAINS


!> Loads the whole dictionary from binary input AND scalar data
!!
!! Load the config PLUS the scalar data written down in Fosite binary format. The
!! structure of the binary format can be found in \link common_dict.f90 \endlink .
FUNCTION LoadConfig(filename) RESULT(res)
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  CHARACTER(LEN=*)        :: filename
  TYPE(Dict_TYP),POINTER  :: res
  !--------------------------------------------------------------------------!
  INTEGER                 :: file, error, offset, keylen, intsize, realsize,  &
                             type, bytes, l, dims(5)
  CHARACTER(LEN=6)        :: magic
  CHARACTER(LEN=2)        :: endian
  CHARACTER(LEN=1)        :: version
  CHARACTER(LEN=4)        :: sizestr
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: keybuf
  CHARACTER(LEN=MAX_CHAR_LEN)           :: key,kf
  REAL,DIMENSION(:,:,:),        POINTER :: ptr3 => null()
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: val

  !--------------------------------------------------------------------------!
  INTENT(IN)         :: filename
  !--------------------------------------------------------------------------!
  NULLIFY(res)
  offset = 1
  file = 5555
  OPEN(file, &
       FILE       =TRIM(filename), &
       STATUS     = 'OLD',      &
       ACCESS     = 'STREAM' ,   &
       ACTION     = 'READ', &
       POSITION   = 'REWIND', &
       IOSTAT     = error)
  READ(file) magic, endian, version, sizestr
  offset = offset + 13
  READ(sizestr, '(I2,I2)') realsize, intsize
  READ(file, IOSTAT=error) keylen
  offset = offset + intsize
  DO WHILE(error.EQ.0)
    key = ''
    ALLOCATE(keybuf(keylen))
    READ(file) keybuf,type,bytes
    WRITE(kf,'(A, I4, A)') '(',keylen,'(A))'
    WRITE(key,kf)keybuf
    DEALLOCATE(keybuf)
    key = TRIM(key)
    offset = offset + keylen + 2*intsize
    dims(:) = 1
    SELECT CASE(type)
    CASE(DICT_REAL_TWOD)
      l = 2
    CASE(DICT_REAL_THREED)
      l = 3
    CASE(DICT_REAL_FOURD)
      l = 4
    CASE(DICT_REAL_FIVED)
      l = 5
    CASE DEFAULT
      l = 0
    END SELECT
    IF(l.GE.2) THEN
      READ(file) dims(1:l)
      bytes = bytes - l*intsize
      offset = offset + l*intsize
    END IF

    ! Here nothing is read in, the data parts are skipped
    ! TODO: The allocation of the field is not good here and just a workaround
    !   - Problem: The data parts in the file need to be skipped somewhow,
    !     but the POS argument is not available in F95, which can be used
    !     on NEC/SX-Ace. At the field decomposition by MPI is done after
    !     reading in the dictionary in SetupFosite.
    SELECT CASE(l)
    CASE(3)
      ALLOCATE(ptr3(dims(1), dims(2), dims(3)))
      READ(file) ptr3
      DEALLOCATE(ptr3)
    CASE DEFAULT
      IF(bytes.GT.0) THEN
        ALLOCATE(val(bytes))
        READ(file) val
        IF(key.EQ.'/config/mesh/decomposition')THEN
          ! Do not set the key for composition. Use new one.
        ELSE
          CALL SetAttr(res,key,val,type)
        END IF
        DEALLOCATE(val)
      END IF
    END SELECT

    offset = offset + bytes

    READ(file, IOSTAT=error) keylen
    offset = offset + intsize
  END DO

  CLOSE(file)
END FUNCTION LoadConfig

!> Loads the datafields from binary input that are necessary for initialization
!!
!! Loads density, xvelocity, yvelocity and pressure. Since datafields can
!! become very large not very process is loading the whole file into memory.
!! Instead MPI routines are used in order to directly load the fields in the
!! according arrays.
SUBROUTINE LoadData(this,filename)
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  CLASS(fosite)                         :: this
  CHARACTER(LEN=*)                      :: filename
  !--------------------------------------------------------------------------!
  INTEGER                               :: unit, error, keylen, &
                                           intsize, realsize, type, bytes, &
                                           l, dims(5)
!  CHARACTER(LEN=64)                     :: keybufsize
  CHARACTER(LEN=13)                     :: header
  CHARACTER(LEN=6)                      :: magic
  CHARACTER(LEN=2)                      :: endian
  CHARACTER(LEN=1),DIMENSION(50)        :: buffer
#ifndef PARALLEL
  INTEGER                               :: offset
  CHARACTER(LEN=1)                      :: version
#else
  INTEGER(KIND=MPI_OFFSET_KIND)         :: offset, offset_0, filesize
  INTEGER                               :: version
#endif
  CHARACTER(LEN=4)                      :: sizestr
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: keybuf
  CHARACTER(LEN=MAX_CHAR_LEN)           :: key,kf
  REAL,DIMENSION(:,:),          POINTER :: ptr2 => null()
  REAL,DIMENSION(:,:,:),        POINTER :: ptr3 => null()
  REAL,DIMENSION(:,:,:,:),      POINTER :: ptr4 => null()
  REAL,DIMENSION(:,:,:,:,:),    POINTER :: ptr5 => null()
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: val
  INTEGER                               :: counter
#ifdef PARALLEL
  INTEGER                               :: handle      !< MPI file handle
  INTEGER                               :: filetype    !< data type for data i/o
  INTEGER                               :: ierror      !< MPI error output
  INTEGER                               :: bufsize     !< output data buffer size
  INTEGER                               :: position
  INTEGER, DIMENSION(2)                 :: gsizes,lsizes,indices,memsizes
  INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status      !< MPI i/o status record
#endif
  !--------------------------------------------------------------------------!
  INTENT(IN)                            :: filename
  INTENT(INOUT)                         :: this
  !--------------------------------------------------------------------------!
  offset = 1
  counter = 0 !temporary

#ifndef PARALLEL
  unit = 5555
  OPEN(unit, &
       FILE       = TRIM(filename), &
       STATUS     = 'OLD', &
       ACCESS     = 'STREAM', &
       ACTION     = 'READ', &
       POSITION   = 'REWIND', &
       IOSTAT     = error)
#else
  ! First create Datatype that will be read in
  gsizes(1) = Sim%Mesh%INUM
  gsizes(2) = Sim%Mesh%JNUM
  lsizes(1) = Sim%Mesh%IMAX-Sim%Mesh%IMIN+1
  lsizes(2) = Sim%Mesh%JMAX-Sim%Mesh%JMIN+1
  indices(1)= Sim%Mesh%IMIN-1
  indices(2)= Sim%Mesh%JMIN-1
  bufsize = PRODUCT(lsizes)
  CALL MPI_Type_create_subarray(2, gsizes, lsizes, indices, MPI_ORDER_FORTRAN,&
       DEFAULT_MPI_REAL,filetype,ierror)
  CALL MPI_Type_commit(filetype,ierror)

  ! Open File
  CALL MPI_File_open(MPI_COMM_WORLD,TRIM(filename),MPI_MODE_RDONLY, &
       MPI_INFO_NULL,handle,error)
!  CALL MPI_File_seek(handle,offset-1,MPI_SEEK_SET,ierror)
#endif

#ifndef PARALLEL
  READ(unit) magic, endian, version, sizestr
#else
  CALL MPI_File_read_all(handle, header, LEN(header), MPI_BYTE, &
    status,ierror)
  WRITE (magic, '(A6)')  header(1:6)
  WRITE (endian, '(A2)')  header(7:8)
  version = IACHAR(header(9:9))
  WRITE (sizestr, '(A4)')  header(10:13)
#endif
  offset = offset + 13
  READ(sizestr, '(I2,I2)') realsize, intsize
#ifndef PARALLEL
  READ(unit, IOSTAT=error) keylen
#else
! TODO: keylen from data does not need to have the same intsize like from the actual run with Fosite
  buffer = ''
  CALL MPI_File_read_all(handle, buffer, intsize, MPI_BYTE, status,ierror)
  keylen = TRANSFER(buffer(1:intsize),keylen)
!  print *, keylen, magic, endian, version, sizestr
#endif
  offset = offset + intsize

!------------------- loop over data ------------------------------------------!
  DO WHILE(error.EQ.0)
    NULLIFY(ptr2)
    key = ''
    buffer = ''
    ALLOCATE(keybuf(keylen))
#ifndef PARALLEL
    READ(unit) keybuf,type,bytes
#else
    CALL MPI_File_read_all(handle, buffer, keylen+2*intsize, MPI_BYTE, &
      status,ierror)
    keybuf = TRANSFER(buffer(1:keylen), keybuf)
    type = TRANSFER(buffer(keylen+1:keylen+intsize), type)
    bytes = TRANSFER(buffer(keylen+intsize+1:keylen+2*intsize), bytes)
#endif
    WRITE(kf,'(A, I4, A)') '(',keylen,'(A))'
    WRITE(key,kf) keybuf
    DEALLOCATE(keybuf)
    offset = offset + keylen + 2*intsize
    dims(:) = 1
    SELECT CASE(type)
    CASE(DICT_REAL_TWOD)
      l = 2
    CASE(DICT_REAL_THREED)
      l = 3
    CASE(DICT_REAL_FOURD)
      l = 4
    CASE(DICT_REAL_FIVED)
      l = 5
    CASE DEFAULT
      l = 0
    END SELECT
    IF(l.GE.2) THEN
#ifndef PARALLEL
      READ(unit) dims(1:l)
#else
      buffer = ''
      CALL MPI_File_read_all(handle,buffer(1:l*intsize),l*intsize,MPI_BYTE,status,error)
      dims = TRANSFER(buffer(1:l*intsize),dims)
#endif
      bytes = bytes - l*intsize
      offset = offset + l*intsize
    END IF

    SELECT CASE(TRIM(key))
    CASE('/timedisc/density')
      Sim%Timedisc%pvar%data4d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
        Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%DENSITY) = 0.0
#ifndef PARALLEL
      READ(unit) Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%DENSITY)
#else
      CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
      CALL MPI_File_read_all(handle, &
        Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%DENSITY), &
        bufsize, DEFAULT_MPI_REAL, status, ierror)
      offset_0 = 0
      CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
    CASE('/timedisc/xvelocity')
      Sim%Timedisc%pvar%data4d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
        Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%XVELOCITY) = 0.0
#ifndef PARALLEL
      READ(unit) Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%XVELOCITY)
#else
      CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
      CALL MPI_File_read_all(handle, &
        Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%XVELOCITY),&
        bufsize, DEFAULT_MPI_REAL, status, ierror)
      offset_0 = 0
      CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
    CASE('/timedisc/yvelocity')
      Sim%Timedisc%pvar%data4d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
        Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%YVELOCITY) = 0.0
#ifndef PARALLEL
      READ(unit) Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%YVELOCITY)
#else
      CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
      CALL MPI_File_read_all(handle, &
        Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%YVELOCITY),&
        bufsize, DEFAULT_MPI_REAL, status, ierror)
      offset_0 = 0
      CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
    CASE('/timedisc/zvelocity')
      Sim%Timedisc%pvar%data4d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
        Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%ZVELOCITY) = 0.0
#ifndef PARALLEL
      READ(unit) Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%ZVELOCITY)
#else
      CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
      CALL MPI_File_read_all(handle, &
        Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,Sim%Physics%ZVELOCITY),&
        bufsize, DEFAULT_MPI_REAL, status, ierror)
      offset_0 = 0
      CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
    CASE('/timedisc/pressure')
      SELECT TYPE (phys => Sim%Physics)
      TYPE IS(physics_euler)
      Sim%Timedisc%pvar%data4d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
        Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,phys%PRESSURE) = 0.0
#ifndef PARALLEL
      READ(unit) Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,phys%PRESSURE)
#else
      CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
      CALL MPI_File_read_all(handle, &
        Sim%Timedisc%pvar%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
        Sim%Mesh%KMIN:Sim%Mesh%KMAX,phys%PRESSURE),&
        bufsize, DEFAULT_MPI_REAL, status, ierror)
      offset_0 = 0
      CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
      END SELECT
    CASE('/physics/bccsound')
      SELECT TYPE (phys => Sim%Physics)
      TYPE IS(physics_eulerisotherm)
        phys%bccsound%data3d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
          Sim%Mesh%KGMIN:Sim%Mesh%KGMAX) = 0.0
#ifndef PARALLEL
        READ(unit) phys%bccsound%data3d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
          Sim%Mesh%KMIN:Sim%Mesh%KMAX)
#else
        CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
        CALL MPI_File_read_all(handle, &
          phys%bccsound%data3d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
          Sim%Mesh%KMIN:Sim%Mesh%KMAX),&
          bufsize, DEFAULT_MPI_REAL, status, ierror)
        offset_0 = 0
        CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
      END SELECT
    CASE('/physics/fccsound')
      SELECT TYPE (phys => Sim%Physics)
      TYPE IS(physics_eulerisotherm)
        phys%fcsound%data4d(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX, &
          Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,1:Sim%Mesh%NFACES) = 0.0
#ifndef PARALLEL
        READ(unit) phys%fcsound%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
          Sim%Mesh%KMIN:Sim%Mesh%KMAX,1:Sim%Mesh%NFACES)
#else
        CALL MPI_File_set_view(handle, offset-1,DEFAULT_MPI_REAL,filetype, 'native', MPI_INFO_NULL, ierror)
        CALL MPI_File_read_all(handle, &
          phys%fcsound%data4d(Sim%Mesh%IMIN:Sim%Mesh%IMAX,Sim%Mesh%JMIN:Sim%Mesh%JMAX, &
          Sim%Mesh%KMIN:Sim%Mesh%KMAX,1:Sim%Mesh%NFACES),&
          bufsize, DEFAULT_MPI_REAL, status, ierror)
        offset_0 = 0
        CALL MPI_File_set_view(handle,offset_0,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierror)
#endif
      END SELECT
    CASE DEFAULT
      SELECT CASE(l)
      CASE(2)
#ifndef PARALLEL
        ALLOCATE(ptr2(dims(1),dims(2)))
        READ(unit) ptr2
        DEALLOCATE(ptr2)
#endif
      CASE(3)
#ifndef PARALLEL
        ALLOCATE(ptr3(dims(1),dims(2),dims(3)))
        READ(unit) ptr3
        DEALLOCATE(ptr3)
#endif
      CASE(4)
#ifndef PARALLEL
        ALLOCATE(ptr4(dims(1),dims(2),dims(3),dims(4)))
        READ(unit) ptr4
        DEALLOCATE(ptr4)
#endif
      CASE(5)
#ifndef PARALLEL
        ALLOCATE(ptr5(dims(1),dims(2),dims(3),dims(4),dims(5)))
        READ(unit) ptr5
        DEALLOCATE(ptr5)
#endif
      CASE DEFAULT
        IF(bytes.GT.0) THEN
#ifndef PARALLEL
          ALLOCATE(val(bytes))
          READ(unit) val
          DEALLOCATE(val)
#endif
        END IF
      END SELECT
    END SELECT
    offset = offset + bytes

#ifndef PARALLEL
    READ(unit, IOSTAT=error) keylen
#else
    buffer = ''
    CALL MPI_File_seek(handle,offset-1,MPI_SEEK_SET,ierror)
    CALL MPI_File_read_all(handle,buffer(1:intsize),intsize,MPI_BYTE,status,ierror)
    keylen = TRANSFER(buffer(1:intsize),keylen)

    CALL MPI_File_get_position(handle,offset_0,ierror)
    CALL MPI_File_get_size(handle,filesize,ierror)
    IF (filesize.LE.offset_0) THEN
      ierror=1
    END IF
    error = ierror
#endif
    offset = offset + intsize
  END DO
#ifndef PARALLEL
  CLOSE(unit)
#else
  CALL MPI_File_close(handle,ierror)
#endif
END SUBROUTINE LoadData

END PROGRAM Restart
