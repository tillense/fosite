!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: restart.f90                                                       #
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
!> Restart fosite using a binary output file
!----------------------------------------------------------------------------!
PROGRAM Restart
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE sources_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  TYPE(Dict_TYP),POINTER :: input
  CHARACTER(LEN=100) :: filename
  LOGICAL :: file_exists
  !--------------------------------------------------------------------------!

  CALL GETARG(1, filename)

  INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
  IF(.NOT.file_exists) THEN
    print *, "restart, ", "Input file does not exist!"
    STOP
  END IF

  CALL InitFosite(Sim)

  input => LoadConfig(TRIM(filename))
  CALL PrintDict(input)
  CALL GetAttr(input, "config", Sim%config)

  CALL PrintDict(Sim%config)

  CALL SetupFosite(Sim)

  CALL CopyDictValues(Sim, Sim%Mesh,input,Sim%IO)

  IF(HasKey(input,"/timedisc/xmomentum")) THEN
    CALL Convert2Primitive(Sim%Physics,Sim%Mesh,Sim%Timedisc%cvar,Sim%Timedisc%pvar)
  ELSE
    CALL Convert2Conservative(Sim%Physics,Sim%Mesh,Sim%Timedisc%pvar,Sim%Timedisc%cvar)
  END IF
  !CALL DeleteDict(input)

  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

FUNCTION LoadConfig(filename) RESULT(res)
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  CHARACTER(LEN=*)   :: filename
  TYPE(Dict_TYP),POINTER :: res
  !--------------------------------------------------------------------------!
  TYPE(Dict_TYP),POINTER :: node
  INTEGER            :: file, error, offset, keylen, intsize, realsize, type, &
                        bytes, l, dims(4)
  CHARACTER(LEN=6)   :: magic
  CHARACTER(LEN=2)   :: endian
  CHARACTER(LEN=1)   :: version
  CHARACTER(LEN=4)   :: sizestr
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: keybuf
  CHARACTER(LEN=MAX_CHAR_LEN) :: key,kf
  REAL,DIMENSION(:,:),POINTER :: ptr2 => null()
  REAL,DIMENSION(:,:,:),POINTER :: ptr3 => null()
  REAL,DIMENSION(:,:,:,:),POINTER :: ptr4 => null()
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: val

  !--------------------------------------------------------------------------!
  INTENT(IN)         :: filename
  !--------------------------------------------------------------------------!
#ifdef FORTRAN_STREAM
    NULLIFY(res)
    offset = 1
    OPEN(file, FILE=TRIM(filename), &
         STATUS     = 'OLD',      &
         ACCESS     = 'STREAM' ,   &
         IOSTAT     = error)
    READ(file, POS=offset) magic, endian, version, sizestr
    offset = offset + 13
    READ(sizestr, '(I2,I2)') realsize,intsize
    !print *,magic,endian,version,realsize,intsize
    READ(file, POS=offset, IOSTAT=error) keylen
    offset = offset + intsize
    DO WHILE(error.EQ.0)
      key = ''
      ALLOCATE(keybuf(keylen))
      READ(file, POS=offset) keybuf,type,bytes
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
      CASE DEFAULT
        l = 0
      END SELECT
      IF(l.GE.2) THEN
        READ(file, POS=offset) dims(1:l)
        bytes = bytes - l*intsize
        offset = offset + l*intsize
      END IF

      SELECT CASE(l)
      CASE(2)
        ALLOCATE(ptr2(dims(1),dims(2)))
        READ(file,POS=offset) ptr2
        CALL SetAttr(res, key, ptr2)
!        DEALLOCATE(ptr2)
      CASE(3)
        ALLOCATE(ptr3(dims(1),dims(2),dims(3)))
        READ(file,POS=offset) ptr3
        node => Dict(key / ptr3)
!        DEALLOCATE(ptr3)
      CASE(4)
        ALLOCATE(ptr4(dims(1),dims(2),dims(3),dims(4)))
        READ(file,POS=offset) ptr4
        node => Dict(key / ptr4)
!        DEALLOCATE(ptr4)
      CASE DEFAULT
        IF(bytes.GT.0) THEN
          ALLOCATE(val(bytes))
          READ(file,POS=offset) val
          CALL SetAttr(res,key,val,type)
          DEALLOCATE(val)
        END IF
      END SELECT
      offset = offset + bytes

!      print *, key, type, bytes, dims
      READ(file, POS=offset, IOSTAT=error) keylen
      offset = offset + intsize
    END DO
    CLOSE(file)
#endif

END FUNCTION LoadConfig

!> Copy complete Dictionary
RECURSIVE SUBROUTINE CopyDictValues(this, Mesh, source, sink, key_)
IMPLICIT NONE
  !------------------------------------------------------------------------!
  TYPE(Fosite_TYP) :: this
  TYPE(Mesh_TYP) :: Mesh
  TYPE(Dict_TYP),POINTER &
                  :: source, sink, node, tmp
  CHARACTER(LEN=*),OPTIONAL :: key_
  CHARACTER(LEN=MAX_CHAR_LEN) :: key,k
  REAL,DIMENSION(:,:),POINTER :: ptr2a,ptr2b
  REAL,DIMENSION(:,:,:),POINTER :: ptr3a,ptr3b
  REAL,DIMENSION(:,:,:,:),POINTER :: ptr4a,ptr4b
  INTEGER,DIMENSION(4) :: dims
  CHARACTER(LEN=1),DIMENSION(:),POINTER :: val
  INTEGER :: i,j
  !------------------------------------------------------------------------!
  key = ''
  IF(PRESENT(key_)) &
    key = key_
  node => source
  DO WHILE(ASSOCIATED(node))
    k = TRIM(key)//'/'//TRIM(GetKey(node))
    IF(HasData(node)) THEN
      NULLIFY(ptr2a,ptr2b,ptr3a,ptr3b,ptr4a,ptr4b)
      dims(:) = 1
      SELECT CASE(GetDatatype(node))
      CASE(DICT_REAL_TWOD)
        CALL GetAttr(sink,k,ptr2a)
        CALL GetAttr(node,GetKey(node),ptr2b)
        dims(1:2) = SHAPE(ptr2b)
      CASE(DICT_REAL_THREED)
        CALL GetAttr(sink,k,ptr3a)
        CALL GetAttr(node,GetKey(node),ptr3b)
        dims(1:3) = SHAPE(ptr3b)
      CASE(DICT_REAL_FOURD)
        CALL GetAttr(sink,k,ptr4a)
        CALL GetAttr(node,GetKey(node),ptr4b)
        dims(1:4) = SHAPE(ptr4b)
      CASE DEFAULT
        CALL GetAttr(sink,TRIM(k),tmp)
        IF(ASSOCIATED(tmp)) THEN
          IF(GetDatatype(tmp).NE.GetDatatype(node)) &
            CALL Error(this,"CopyDictValues","Types of key '" &
              //TRIM(k)//"' are not matching.")
          CALL SetData(tmp,GetData(node))
        ELSE
          CALL SetAttr(sink,TRIM(key),GetData(node),GetDatatype(node))
        END IF
      END SELECT
      IF(PRODUCT(dims).GT.1) THEN
        DO j = 1,dims(4)
          DO i = 1,dims(3)
            IF(ASSOCIATED(ptr3b)) THEN
              ptr2a => ptr3a(:,:,i)
              ptr2b => ptr3b(:,:,i)
            ELSE IF(ASSOCIATED(ptr4b)) THEN
              ptr2a => ptr4a(:,:,i,j)
              ptr2b => ptr4b(:,:,i,j)
            END IF
            IF((dims(1).EQ.Mesh%INUM).AND.&
               (dims(2).EQ.Mesh%JNUM)) THEN
              ptr2a = ptr2b(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
            ELSE IF((dims(1).EQ.Mesh%INUM+1).AND.&
                    (dims(2).EQ.Mesh%JNUM+1)) THEN
              ptr2a = ptr2b(Mesh%IMIN:Mesh%IMAX+1,Mesh%JMIN:Mesh%JMAX+1)
            ELSE
              ptr2a = ptr2b
            END IF
          END DO
        END DO
      END IF
    END IF
    IF(HasChild(node)) THEN
      CALL CopyDictValues(Sim, Mesh, GetChild(node), sink, TRIM(k))
    END IF
    node => GetNext(node)
  END DO
END SUBROUTINE CopyDictValues

END PROGRAM Restart
