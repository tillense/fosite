!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: common_dict.f90                                                   #
!#                                                                           #
!# Copyright (C) 2015 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!> \author Manuel
!!
!! \brief Dictionary for generic data types.
!!
!! This module defines a dictionary (key/value storage) of arbitrary data
!! types. At the moment the following are defined:
!!  1. Integer
!!  2. Real
!!  3. Character(len=MAX_CHAR_LEN)
!!  4. Logical
!!  5. 1D Real array
!!  6. Pointer to 2D Real array
!!  7. Pointer to 3D Real array
!!  8. Pointer to 4D Real array
!!  9. 1D Integer array
!! 10. Pointer to Real
!! 11. Pointer to Integer
!! 12. Pointer to 5D Real array
!!
!! The dictionary is implemented as Trie(1) with comprepressed keys, which is
!! called radix tree(2). Each node can hold data with a defined type, which
!! is stored in a generic container with the TRANSFER intrinsic(3). This
!! technique has been explained in (4), see also (5). Each node can have
!! severall children. 'child' points to the first child and other can be
!! found by iterating to the 'next' one.
!!
!! (1): http://en.wikipedia.org/wiki/Trie
!! (2): http://en.wikipedia.org/wiki/Radix_tree
!! (3): http://fortranwiki.org/fortran/show/transfer
!! (4): \cite blevins2009
!! (5): http://fortranwiki.org/fortran/show/gen_list
!!
!! \extends common_types
!! \ingroup dict
!----------------------------------------------------------------------------!
MODULE common_dict
  !USE dict_common, InitDict_common => InitDict
  USE logging_base_mod
  !--------------------------------------------------------------------------!
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE SetAttr
    MODULE PROCEDURE SetAttr0a, SetAttr0b,                 &
        SetAttr1, SetAttr2, SetAttr3, SetAttr4, SetAttr5,  &
        SetAttr6, SetAttr7, SetAttr8, SetAttr9, SetAttr10, &
        SetAttr11, SetAttr12
  END INTERFACE SetAttr
  INTERFACE GetAttr
    MODULE PROCEDURE GetAttr0,                             &
        GetAttr1, GetAttr2, GetAttr3, GetAttr4, GetAttr5,  &
        GetAttr6, GetAttr7, GetAttr8, GetAttr9, GetAttr10, &
        GetAttr11, GetAttr12
  END INTERFACE GetAttr
  INTERFACE Ref
    MODULE PROCEDURE Ref1, Ref2
  END INTERFACE Ref
  INTERFACE OPERATOR (/)
    MODULE PROCEDURE Assign0, Assign1, Assign2, Assign3, Assign4, &
        Assign5, Assign6, Assign7, Assign8, Assign9, Assign10,    &
        Assign11, Assign12
  END INTERFACE
  !> \endcond
  ! constants
  INTEGER, PARAMETER :: MAX_CHAR_LEN     = 128
  INTEGER, PARAMETER :: DICT_NONE        = 0
  INTEGER, PARAMETER :: DICT_INT         = 1
  INTEGER, PARAMETER :: DICT_REAL        = 2
  INTEGER, PARAMETER :: DICT_CHAR        = 3
  INTEGER, PARAMETER :: DICT_BOOL        = 4
  INTEGER, PARAMETER :: DICT_REAL_ONED   = 5
  INTEGER, PARAMETER :: DICT_REAL_TWOD   = 6
  INTEGER, PARAMETER :: DICT_REAL_THREED = 7
  INTEGER, PARAMETER :: DICT_REAL_FOURD  = 8
  INTEGER, PARAMETER :: DICT_INT_ONED    = 9
  INTEGER, PARAMETER :: DICT_REAL_P      = 10
  INTEGER, PARAMETER :: DICT_INT_P       = 11
  INTEGER, PARAMETER :: DICT_REAL_FIVED  = 12
#define TYPE_DICT_KEY CHARACTER(LEN=MAX_CHAR_LEN)
#define TYPE_DICT_INT INTEGER
#define TYPE_DICT_REAL REAL
#define TYPE_DICT_CHAR CHARACTER(LEN=MAX_CHAR_LEN)
#define TYPE_DICT_BOOL LOGICAL
#define TYPE_DICT_REAL_ONED REAL, DIMENSION(:)
#define TYPE_DICT_REAL_TWOD REAL, DIMENSION(:,:), POINTER
#define TYPE_DICT_REAL_THREED REAL, DIMENSION(:,:,:), POINTER
#define TYPE_DICT_REAL_FOURD REAL, DIMENSION(:,:,:,:), POINTER
#define TYPE_DICT_REAL_FIVED REAL, DIMENSION(:,:,:,:,:), POINTER
#define TYPE_DICT_INT_ONED INTEGER, DIMENSION(:)
#define TYPE_DICT_REAL_P REAL, POINTER
#define TYPE_DICT_INT_P INTEGER, POINTER
! Type of the mold, which holds all the data
#define TYPE_DICT_MOLD CHARACTER(LEN=1), DIMENSION(:)
  ! common data structure
  TYPE_DICT_MOLD, ALLOCATABLE :: mold
  TYPE Dict_TYP
     PRIVATE
     CHARACTER(LEN=MAX_CHAR_LEN) :: key = ""
     INTEGER                     :: type = DICT_NONE
     TYPE_DICT_MOLD, POINTER     :: value => null()
     TYPE(Dict_TYP), POINTER     :: child => null()
     TYPE(Dict_TYP), POINTER     :: next => null()
  END TYPE Dict_TYP
  TYPE real_t
    TYPE_DICT_REAL_P       :: p
  END TYPE
  TYPE int_t
    TYPE_DICT_INT_P        :: p
  END TYPE
  TYPE real_twod_t
    TYPE_DICT_REAL_TWOD    :: p
  END TYPE
  TYPE real_threed_t
    TYPE_DICT_REAL_THREED  :: p
  END TYPE
  TYPE real_fourd_t
    TYPE_DICT_REAL_FOURD   :: p
  END TYPE
  TYPE real_fived_t
    TYPE_DICT_REAL_FIVED   :: p
  END TYPE
  TYPE(logging_base), SAVE :: this
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Dict_TYP, real_t, int_t, &
       ! constants
       MAX_CHAR_LEN, DICT_NONE, DICT_INT, DICT_REAL, DICT_CHAR, DICT_BOOL, &
       DICT_REAL_ONED, DICT_REAL_TWOD, DICT_REAL_THREED, DICT_REAL_FOURD,  &
       DICT_INT_ONED, DICT_REAL_P, DICT_INT_P, DICT_REAL_FIVED,            &
       ! methods
       SetAttr,       &
       GetAttr,       &
       GetNext,       &
       GetChild,      &
       GetKey,        &
       GetDataSize,   &
       GetDatatype,   &
       GetData,       &
       SetData,       &
       HasChild,      &
       HasData,       &
       HasKey,        &
       Ref,           &
       OPERATOR(/),   &
       Dict,          &
       CopyHierarchy, &
       CopyDict,      &
       PrintDict,     &
       InitDict,      &
       DeleteDict
  !--------------------------------------------------------------------------!

CONTAINS
  SUBROUTINE InitDict()
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, PARAMETER :: type = 0
    CHARACTER(LEN=32), PARAMETER :: name = "Dictionary"
    !------------------------------------------------------------------------!
    !CALL InitDict_common(this, type, name)
  END SUBROUTINE InitDict


  !> Search for the path in 'key' beginning at root and return a pointer to
  !! this node in 'res'. If create is set to true (default is false), the
  !! path is created if it is not existing.
  !! If create=false and the path does not exist, res is null().
  FUNCTION FindPath(root, key, create) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: root
    CHARACTER(LEN=*)       :: key
    LOGICAL, OPTIONAL      :: create
    TYPE(Dict_TYP),POINTER :: res
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: parent,node
    LOGICAL                :: c
    TYPE_DICT_CHAR         :: k, key_
    !------------------------------------------------------------------------!
    INTENT(IN)             :: key, create
    !------------------------------------------------------------------------!
    c = .FALSE.
    IF(PRESENT(create)) &
      c = create
    key_ = TRIM(key)
    NULLIFY(parent)
    node => root
    k = Tokenize(key_)
    DO WHILE(LEN_TRIM(k).GT.0)
      node => FindChild(node,TRIM(k))
      IF(ASSOCIATED(node)) THEN
        parent => node
        node => node%child
      ELSE
        IF(c) THEN
          ALLOCATE(node)
          node%key = k
          IF(ASSOCIATED(parent)) THEN
            IF(ASSOCIATED(parent%child)) THEN
              parent => GetLast(parent%child)
              parent%next => node
            ELSE
              parent%child => node
            END IF
          ELSE
            IF(ASSOCIATED(root)) THEN
              parent => GetLast(root)
              parent%next => node
            ELSE
              root => node
            END IF
          END IF
          parent => node
        ELSE
          NULLIFY(parent)
          k = ''
          key_ = ''
        END IF
      END IF
      k = Tokenize(key_)
    END DO
    res => parent

    IF(c.EQV..TRUE..AND..NOT.ASSOCIATED(res)) &
      CALL this%Error("FindPath","Create was activated, so res should be associated.")
  END FUNCTION FindPath

  !> Set the dictionary 'value' as child at the path 'key' relative to 'root'.
  !! If a child at this path is already existing, it will be deleted.
  RECURSIVE SUBROUTINE SetAttr0a(root, key, value)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE(Dict_TYP), TARGET  :: value
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: node, parent, val
    TYPE_DICT_CHAR          :: k, str
    INTEGER                 :: status
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key
    !------------------------------------------------------------------------!
    node => value
    DO WHILE(ASSOCIATED(node))
      WRITE(k,'(A,A,A)') TRIM(key),'/',TRIM(node%key)
      parent => node

      IF(ASSOCIATED(node%child)) THEN
        CALL SetAttr0a(root, TRIM(k), node%child)
      ELSE
        CALL SetAttr0b(root, TRIM(k), node%value, node%type)
      END IF
      node => node%next
      CALL DeleteNode(parent,k)
    END DO
  END SUBROUTINE SetAttr0a

  !> Create an empty node at path 'key' relative to 'root', if value and type
  !! are not defined. If they are defined, also fill the node with data.
  SUBROUTINE SetAttr0b(root, key, value, type)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: root
    CHARACTER(LEN=*)       :: key
    !INTEGER, OPTIONAL :: type
    !TYPE_DICT_MOLD, OPTIONAL :: value
    INTEGER                :: type
    TYPE_DICT_MOLD         :: value
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: node, last
    TYPE_DICT_CHAR         :: k, key_
    INTEGER                :: i
    !------------------------------------------------------------------------!
    INTENT(IN)             :: key, type, value
    !------------------------------------------------------------------------!
    node => FindPath(root,key,.TRUE.)
    node%type = type
    CALL SetData(node,value)
  END SUBROUTINE SetAttr0b

  !> Retrieve the node at path 'key' relative to 'root'. The result will be
  !! given as third argument 'parent'. If the path can not be found, the
  !! result is null().
  SUBROUTINE GetAttr0a(root, key, parent)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: root, parent, child
    CHARACTER(LEN=*)       :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_CHAR         :: k, key_
    !------------------------------------------------------------------------!
    INTENT(IN)             :: key
    !------------------------------------------------------------------------!
    parent => root
    child => root
    key_ = TRIM(key)
    DO WHILE(ASSOCIATED(child).AND.LEN_TRIM(key_).GT.0)
      k = Tokenize(key_)
      parent => FindChild(child,k)
      IF(ASSOCIATED(parent)) &
        child => parent%child
    END DO
    IF(LEN_TRIM(key_).GT.0) &
      NULLIFY(parent)
  END SUBROUTINE GetAttr0a

  !> Retrieve the data 'value' of kind 'type' at path 'key' relative to
  !! 'root'. If the path can not be found and default is not defined, an
  !! error is raised.
  !! If default is present, it is set and returned instead.
  SUBROUTINE GetAttr0b(root, key, type, value, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: root
    CHARACTER(LEN=*)         :: key
    INTEGER                  :: type
    TYPE_DICT_MOLD, POINTER  :: value
    TYPE_DICT_MOLD, OPTIONAL :: default
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: node
    CHARACTER(LEN=10)        :: b1,b2
    !------------------------------------------------------------------------!
    INTENT(IN)               :: key, type, default
    !------------------------------------------------------------------------!
    CALL GetAttr0a(root, key, node)
    IF(.NOT.ASSOCIATED(node)) THEN
      IF(PRESENT(default)) THEN
        CALL SetAttr0b(root,key,default,type)
        CALL GetAttr0a(root, key, node)
        IF(.NOT.ASSOCIATED(node)) &
          CALL this%Error("GetAttr", "Setting a default value has gone wrong.")
      ELSE
        CALL this%Error("GetAttr", "Couldn't find key '"//TRIM(key)//"'.")
      END IF
    END IF
    IF(type.NE.node%type) THEN
      WRITE(b1,"(I4.4)") node%type
      WRITE(b2,"(I4.4)") type
      CALL this%Error("GetAttr", "Key '"//TRIM(key)//"' is of type "//TRIM(b1) &
        //", but type "//TRIM(b2)//" requested.")
    END IF
    value => node%value
  END SUBROUTINE GetAttr0b

  !> Find the direct child with key 'key' in a list of childs. 'root' points
  !! to the first child. If 'key' is not found, null() is returned.
  FUNCTION FindChild(root,key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*), INTENT(IN) :: key
    TYPE(Dict_TYP), POINTER      :: root,res
    !------------------------------------------------------------------------!
    res => root
    DO WHILE(ASSOCIATED(res))
      IF(TRIM(res%key).EQ.TRIM(key)) EXIT
      res => res%next
    END DO
  END FUNCTION FindChild

  !> Get the pointer to the next child.
  FUNCTION GetNext(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root, res
    !------------------------------------------------------------------------!
    res => root%next
  END FUNCTION GetNext

  !> Get the pointer to the last child.
  FUNCTION GetLast(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root, res
    !------------------------------------------------------------------------!
    res => root
    IF(ASSOCIATED(res)) THEN
      DO WHILE(ASSOCIATED(res%next))
        res => res%next
      END DO
    END IF
  END FUNCTION GetLast

  !> Get the pointer to a direct child of the pointer 'root'.
  FUNCTION GetChild(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root, res
    !------------------------------------------------------------------------!
    res => root%child
  END FUNCTION GetChild

  !> Get the key of pointer 'root'.
  FUNCTION GetKey(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    TYPE_DICT_CHAR          :: res
    !------------------------------------------------------------------------!
    res = root%key
  END FUNCTION GetKey

  !> Get the size of the data in node 'root'. If there is no data 0 is
  !! returned.
  !! note: This is also the byte size of the data.
  FUNCTION GetDatasize(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    INTEGER                 :: res
    !------------------------------------------------------------------------!
    IF(ASSOCIATED(root%value)) THEN
      res = SIZE(root%value)
    ELSE
      res = 0
    END IF
  END FUNCTION GetDatasize

  !> Return the datatype of node 'root'
  FUNCTION GetDatatype(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    INTEGER                 :: res
    !------------------------------------------------------------------------!
    res = root%type
  END FUNCTION GetDatatype

  !> Return the datatype of node 'root'
  FUNCTION GetData(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    TYPE_DICT_MOLD, POINTER :: res
    !------------------------------------------------------------------------!
    res => root%value
  END FUNCTION GetData

  !> Check if the node 'root' has one or more children.
  FUNCTION HasChild(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    LOGICAL                 :: res
    !------------------------------------------------------------------------!
    res = ASSOCIATED(root%child)
  END FUNCTION HasChild

  !> Checks if the node 'root' has data associated.
  FUNCTION HasData(root) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    LOGICAL                 :: res
    !------------------------------------------------------------------------!
    res = ASSOCIATED(root%value)
  END FUNCTION HasData

  !> Checks if a node with key 'key' exists.
  FUNCTION HasKey(root, key) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    LOGICAL                 :: res
    !------------------------------------------------------------------------!
    res = ASSOCIATED(FindPath(root,TRIM(key)))
  END FUNCTION HasKey

  !> Set data of 'node' ot 'val'
  SUBROUTINE SetData(node, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: node
    TYPE_DICT_MOLD          :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: val
    !------------------------------------------------------------------------!
    IF(SIZE(val).LE.0) &
      CALL this%Error("SetData", "Array size smaller 0 is not possible.")
    IF(ASSOCIATED(node%value)) &
      DEALLOCATE(node%value)
    ALLOCATE(node%value(SIZE(val)))
    node%value = val
  END SUBROUTINE SetData

  !> Cuts a path into two tokens, which is explained best with an example:
  !! back=.FALSE.: key = /sources/grav/mass => res = sources, key = /grav/mass
  !! back=.TRUE.:  key = /sources/grav/mass => res = mass, key = /sources/grav
  FUNCTION Tokenize(key,back) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=MAX_CHAR_LEN) :: key
    LOGICAL,OPTIONAL            :: back
    TYPE_DICT_CHAR              :: res
    !------------------------------------------------------------------------!
    LOGICAL                     :: back_
    INTEGER                     :: i
    !------------------------------------------------------------------------!
    INTENT(IN)                  :: back
    INTENT(INOUT)               :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(back)) THEN
      back_ = back
    ELSE
      back_ = .FALSE.
    END IF
    res = ''
    IF(key(1:1).EQ.'/') key = key(2:)
    i = SCAN(key,'/',back_)
    IF(i.GT.0) THEN
      IF(.NOT.back_) THEN
        res = key(1:i-1)
        key = key(i+1:)
      ELSE
        res = key(i+1:)
        key = key(1:i-1)
      END IF
    ELSE
      res = key
      key = ''
    END IF
  END FUNCTION Tokenize

  !> Construct a new dictionary from several key/value pairs.
  !! Together with the Assign subroutine and overloading the '/'-operator,
  !! this makes the syntactic sugar available to write the following syntax:
  !!   physics => Dict( "problem" / EULER2D, &
  !!                    "gamma"   / GAMMA)
  FUNCTION  Dict(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,&
                 n18,n19,n20) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER           :: res,n1
    TYPE(Dict_TYP), POINTER, OPTIONAL :: n2,n3,n4,n5,n6,n7,n8,n9,n10,n11, &
                                         n12,n13,n14,n15,n16,n17,n18,n19,n20
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr0a(res,'',n1)
    IF(PRESENT( n2)) CALL SetAttr0a(res,'', n2)
    IF(PRESENT( n3)) CALL SetAttr0a(res,'', n3)
    IF(PRESENT( n4)) CALL SetAttr0a(res,'', n4)
    IF(PRESENT( n5)) CALL SetAttr0a(res,'', n5)
    IF(PRESENT( n6)) CALL SetAttr0a(res,'', n6)
    IF(PRESENT( n7)) CALL SetAttr0a(res,'', n7)
    IF(PRESENT( n8)) CALL SetAttr0a(res,'', n8)
    IF(PRESENT( n9)) CALL SetAttr0a(res,'', n9)
    IF(PRESENT(n10)) CALL SetAttr0a(res,'',n10)
    IF(PRESENT(n11)) CALL SetAttr0a(res,'',n11)
    IF(PRESENT(n12)) CALL SetAttr0a(res,'',n12)
    IF(PRESENT(n13)) CALL SetAttr0a(res,'',n13)
    IF(PRESENT(n14)) CALL SetAttr0a(res,'',n14)
    IF(PRESENT(n15)) CALL SetAttr0a(res,'',n15)
    IF(PRESENT(n16)) CALL SetAttr0a(res,'',n16)
    IF(PRESENT(n17)) CALL SetAttr0a(res,'',n17)
    IF(PRESENT(n18)) CALL SetAttr0a(res,'',n18)
    IF(PRESENT(n19)) CALL SetAttr0a(res,'',n19)
    IF(PRESENT(n20)) CALL SetAttr0a(res,'',n20)
  END FUNCTION Dict

  SUBROUTINE SetAttr1(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_INT           :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_INT)
  END SUBROUTINE SetAttr1

  SUBROUTINE SetAttr2(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_REAL          :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_REAL)
  END SUBROUTINE SetAttr2

  SUBROUTINE SetAttr3(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key, val
    TYPE_DICT_CHAR          :: val_
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    val_ = TRIM(val)
    CALL SetAttr0b(root,key,TRANSFER(val_,mold),DICT_CHAR)
  END SUBROUTINE SetAttr3

  SUBROUTINE SetAttr4(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_BOOL          :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_BOOL)
  END SUBROUTINE SetAttr4

  SUBROUTINE SetAttr5(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_REAL_ONED     :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_REAL_ONED)
  END SUBROUTINE SetAttr5

  SUBROUTINE SetAttr6(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER      :: root
    CHARACTER(LEN=*)             :: key
    REAL, DIMENSION(:,:), TARGET :: val
    TYPE(real_twod_t)            :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                   :: key,val
    !------------------------------------------------------------------------!
    c%p => val
    CALL SetAttr0b(root,key,TRANSFER(c,mold),DICT_REAL_TWOD)
  END SUBROUTINE SetAttr6

  SUBROUTINE SetAttr7(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER        :: root
    CHARACTER(LEN=*)               :: key
    REAL, DIMENSION(:,:,:), TARGET :: val
    TYPE(real_threed_t)            :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                     :: key, val
    !------------------------------------------------------------------------!
    c%p => val
    CALL SetAttr0b(root,key,TRANSFER(c,mold),DICT_REAL_THREED)
  END SUBROUTINE SetAttr7

  SUBROUTINE SetAttr8(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER          :: root
    CHARACTER(LEN=*)                 :: key
    REAL, DIMENSION(:,:,:,:), TARGET :: val
    TYPE(real_fourd_t)               :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                       :: key, val
    !------------------------------------------------------------------------!
    c%p => val
    CALL SetAttr0b(root,key,TRANSFER(c,mold),DICT_REAL_FOURD)
  END SUBROUTINE SetAttr8

  SUBROUTINE SetAttr9(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_INT_ONED      :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_INT_ONED)
  END SUBROUTINE SetAttr9

  SUBROUTINE SetAttr10(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE(real_t)            :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_REAL_P)
  END SUBROUTINE SetAttr10

  SUBROUTINE SetAttr11(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    CHARACTER(LEN=*)        :: key
    TYPE(int_t)             :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    CALL SetAttr0b(root,key,TRANSFER(val,mold),DICT_INT_P)
  END SUBROUTINE SetAttr11

  SUBROUTINE SetAttr12(root, key, val)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER            :: root
    CHARACTER(LEN=*)                   :: key
    REAL, DIMENSION(:,:,:,:,:), TARGET :: val
    TYPE(real_fived_t)                 :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: key, val
    !------------------------------------------------------------------------!
    c%p => val
    CALL SetAttr0b(root,key,TRANSFER(c,mold),DICT_REAL_FIVED)
  END SUBROUTINE SetAttr12

  RECURSIVE SUBROUTINE PrintDict(root, prefix)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), TARGET     ::root
    CHARACTER(LEN=*), OPTIONAL :: prefix
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: node
    TYPE_DICT_CHAR          :: prefix_
    TYPE_DICT_CHAR          :: s,str
    !------------------------------------------------------------------------!
    node => root
    prefix_ = ''
    IF(PRESENT(prefix)) &
      prefix_ = TRIM(prefix)
    DO WHILE(ASSOCIATED(node))
      WRITE(s,'(A,A,A)')TRIM(prefix_),'/',TRIM(node%key)
      IF(ASSOCIATED(node%value)) THEN
        WRITE(str,'(A,I2,A,A,A,I4)') "type=",node%type,", key=",TRIM(s),&
          ", size=",SIZE(node%value)
      ELSE
        WRITE(str,'(A,I2,A,A)') "type=",node%type,", key=",TRIM(s)
      END IF
      !CALL Info(this,str)
      IF(ASSOCIATED(node%child)) &
        CALL PrintDict(node%child, s)
      node => node%next
    END DO
  END SUBROUTINE PrintDict

  !> Copy complete Dictionary
  RECURSIVE SUBROUTINE CopyDict(root, outdir)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER     :: root, outdir, dir, odir, tmp
    CHARACTER(LEN=MAX_CHAR_LEN) :: key
    !------------------------------------------------------------------------!
    dir => root
    NULLIFY(outdir)
    NULLIFY(odir)
    DO WHILE(ASSOCIATED(dir))
      IF(ASSOCIATED(odir)) THEN
        ALLOCATE(odir%next)
        odir => odir%next
      ELSE
        ALLOCATE(odir)
        outdir => odir
      END IF
      odir%type = dir%type
      odir%key = dir%key
      IF(ASSOCIATED(dir%value)) THEN
        CALL SetData(odir,dir%value)
      END IF

      IF(ASSOCIATED(dir%child)) THEN
        CALL CopyDict(dir%child, tmp)
        odir%child => tmp
      END IF
      dir => dir%next
    END DO
  END SUBROUTINE CopyDict


  !> Copy all nodes, which have children from 'root' to 'outdir'.
  RECURSIVE SUBROUTINE CopyHierarchy(root, outdir)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER     :: root, outdir, dir, odir
    CHARACTER(LEN=MAX_CHAR_LEN) :: key
    !------------------------------------------------------------------------!
    dir => root
    NULLIFY(outdir)
    NULLIFY(odir)
    DO WHILE(ASSOCIATED(dir))
       IF(ASSOCIATED(dir%child)) THEN
         IF(ASSOCIATED(odir)) THEN
           ALLOCATE(odir%next)
           odir => odir%next
         ELSE
           ALLOCATE(odir)
           outdir => odir
         END IF
         odir%key = dir%key
         CALL CopyHierarchy(dir%child,odir%child)
       END IF
       dir => dir%next
    END DO
  END SUBROUTINE CopyHierarchy

  !> Return the node at path 'key' relative to 'root' in 'res'. If this node
  !! has no data, but a child (e.g. a directory), return the child instead.
  SUBROUTINE GetAttr0(root, key, res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root, res
    CHARACTER(LEN=*)        :: key
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key
    !------------------------------------------------------------------------!
    CALL GetAttr0a(root, key, res)
    IF(ASSOCIATED(res)) THEN
      IF(.NOT.HasData(res).AND.HasChild(res)) THEN
        res => res%child
      END IF
    END IF
  END SUBROUTINE GetAttr0

  SUBROUTINE GetAttr1(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    TYPE_DICT_INT           :: res
    TYPE_DICT_INT, OPTIONAL :: default
    CHARACTER(LEN=*)        :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER :: value
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, default
    INTENT(INOUT)           :: res
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_INT, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_INT, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr1

  SUBROUTINE GetAttr2(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: root
    TYPE_DICT_REAL           :: res
    TYPE_DICT_REAL, OPTIONAL :: default
    CHARACTER(LEN=*)         :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER  :: value
    !------------------------------------------------------------------------!
    INTENT(IN)               :: key, default
    INTENT(INOUT)            :: res
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_REAL, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr2

  SUBROUTINE GetAttr3(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: root
    TYPE_DICT_CHAR           :: res
    TYPE_DICT_CHAR, OPTIONAL :: default
    CHARACTER(LEN=*)         :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER  :: value
    !------------------------------------------------------------------------!
    INTENT(IN)               :: key
    INTENT(INOUT)            :: res
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_CHAR, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_CHAR, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr3

  SUBROUTINE GetAttr4(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: root
    TYPE_DICT_BOOL           :: res
    TYPE_DICT_BOOL, OPTIONAL :: default
    CHARACTER(LEN=*)         :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER  :: value
    !------------------------------------------------------------------------!
    INTENT(IN)               :: key
    INTENT(INOUT)            :: res
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_BOOL, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_BOOL, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr4

  SUBROUTINE GetAttr5(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER       :: root
    TYPE_DICT_REAL_ONED           :: res
    TYPE_DICT_REAL_ONED, OPTIONAL :: default
    CHARACTER(LEN=*)              :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER       :: value
    !------------------------------------------------------------------------!
    INTENT(IN)                    :: key
    INTENT(INOUT)                 :: res
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_REAL_ONED, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL_ONED, value)
    END IF
    IF(SIZE(TRANSFER(value,res)).NE.SIZE(res)) &
      CALL this%Error("GetAttr5","1D array with key '" // TRIM(key) &
        // "' has the wrong size.")

    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr5

  SUBROUTINE GetAttr6(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER       :: root
    TYPE_DICT_REAL_TWOD           :: res
    TYPE_DICT_REAL_TWOD, OPTIONAL :: default
    CHARACTER(LEN=*)              :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER       :: value
    TYPE(real_twod_t)             :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                    :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      c%p => default
      CALL GetAttr0b(root, key, DICT_REAL_TWOD, value, TRANSFER(c,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL_TWOD, value)
    END IF
    c = TRANSFER(value,c)
    res => c%p
  END SUBROUTINE GetAttr6

  SUBROUTINE GetAttr7(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER         :: root
    TYPE_DICT_REAL_THREED           :: res
    TYPE_DICT_REAL_THREED, OPTIONAL :: default
    CHARACTER(LEN=*)                :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER         :: value
    TYPE(real_threed_t)             :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      c%p => default
      CALL GetAttr0b(root, key, DICT_REAL_THREED, value, TRANSFER(c,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL_THREED, value)
    END IF
    c = TRANSFER(value,c)
    res => c%p
  END SUBROUTINE GetAttr7

  SUBROUTINE GetAttr8(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER        :: root
    TYPE_DICT_REAL_FOURD           :: res
    TYPE_DICT_REAL_FOURD, OPTIONAL :: default
    CHARACTER(LEN=*)               :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER        :: value
    TYPE(real_fourd_t)             :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                     :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      c%p => default
      CALL GetAttr0b(root, key, DICT_REAL_FOURD, value, TRANSFER(c,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL_FOURD, value)
    END IF
    c = TRANSFER(value,c)
    res => c%p
  END SUBROUTINE GetAttr8

  SUBROUTINE GetAttr9(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER      :: root
    TYPE_DICT_INT_ONED           :: res
    TYPE_DICT_INT_ONED, OPTIONAL :: default
    CHARACTER(LEN=*)             :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER      :: value
    !------------------------------------------------------------------------!
    INTENT(IN)                   :: key
    INTENT(INOUT)                :: res
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_INT_ONED, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_INT_ONED, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr9

  SUBROUTINE GetAttr10(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    TYPE(real_t)            :: res
    TYPE(real_t), OPTIONAL  :: default
    CHARACTER(LEN=*)        :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER :: value
    TYPE(real_t)            :: c
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_REAL_P, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL_P, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr10

  SUBROUTINE GetAttr11(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root
    TYPE(int_t)             :: res
    TYPE(int_t), OPTIONAL   :: default
    CHARACTER(LEN=*)        :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER :: value
    TYPE(int_t)             :: c
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      CALL GetAttr0b(root, key, DICT_INT_P, value, TRANSFER(default,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_INT_P, value)
    END IF
    res = TRANSFER(value,res)
  END SUBROUTINE GetAttr11

  SUBROUTINE GetAttr12(root, key, res, default)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER        :: root
    TYPE_DICT_REAL_FIVED           :: res
    TYPE_DICT_REAL_FIVED, OPTIONAL :: default
    CHARACTER(LEN=*)               :: key
    !------------------------------------------------------------------------!
    TYPE_DICT_MOLD, POINTER        :: value
    TYPE(real_fived_t)             :: c
    !------------------------------------------------------------------------!
    INTENT(IN)                     :: key
    !------------------------------------------------------------------------!
    IF(PRESENT(default)) THEN
      c%p => default
      CALL GetAttr0b(root, key, DICT_REAL_FIVED, value, TRANSFER(c,mold))
    ELSE
      CALL GetAttr0b(root, key, DICT_REAL_FIVED, value)
    END IF
    c = TRANSFER(value,c)
    res => c%p
  END SUBROUTINE GetAttr12

  SUBROUTINE DeleteNode(node,k)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: node
    CHARACTER(LEN=*)        :: k
    LOGICAL, SAVE           :: first=.TRUE.
    !------------------------------------------------------------------------!
    INTEGER            :: status
    CHARACTER(LEN=512) :: str
    !------------------------------------------------------------------------!
    NULLIFY(node%child,node%next)
    IF(ASSOCIATED(node%value)) &
      DEALLOCATE(node%value)
    DEALLOCATE(node,stat=status)
    IF(status.NE.0.AND.first) THEN
      ! This warning (no. 195) does occur on the SX ACE. The definite reason
      ! is unknown. It is probably related to deallocating a pointer, which
      ! has been associated to a TARGET subroutine parameter. But this cannot
      ! be change, since a TARGET instead of a POINTER is required for
      ! overloading the '/' operator.
      WRITE(str,'(A,A,A,I3)')&
        "Deallocating key '",TRIM(k),&
        "' throws the error nio.: ",status
      !CALL Warning(this,"SetAttr0a",TRIM(str))
      WRITE(str,'(A,A,A,I3)')&
        "More invocations of this error will be suppressed."
      !CALL Warning(this,"SetAttr0a",TRIM(str))
      first = .FALSE.
    END IF
  END SUBROUTINE DeleteNode

  !> Delete the dictionary 'root' and all subnodes.
  RECURSIVE SUBROUTINE DeleteDict(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: root, node, next, child
    !------------------------------------------------------------------------!
    INTEGER                 :: status
    TYPE_DICT_CHAR          :: k,str
    !------------------------------------------------------------------------!
    node => root
    DO WHILE(ASSOCIATED(node))
        next => node%next
        child => node%child
        k = TRIM(node%key)
        CALL DeleteNode(node,TRIM(k))
        IF (ASSOCIATED(child)) &
           CALL DeleteDict(child)
        node => next
    END DO
    NULLIFY(root)
  END SUBROUTINE DeleteDict

  FUNCTION Assign0(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE(Dict_TYP), TARGET  :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key,val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign0

  FUNCTION Assign1(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_INT           :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign1

  FUNCTION Assign2(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_REAL          :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign2

  FUNCTION Assign3(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key, val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign3

  FUNCTION Assign4(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_BOOL          :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign4

  FUNCTION Assign5(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_REAL_ONED     :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign5

  FUNCTION Assign6(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER      :: res
    CHARACTER(LEN=*)             :: key
    REAL, DIMENSION(:,:), TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)                   :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign6

  FUNCTION Assign7(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER        :: res
    CHARACTER(LEN=*)               :: key
    REAL, DIMENSION(:,:,:), TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)                     :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign7

  FUNCTION Assign8(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER          :: res, tmp
    CHARACTER(LEN=*)                 :: key
    REAL, DIMENSION(:,:,:,:), TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)                       :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign8

  FUNCTION Assign9(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE_DICT_INT_ONED      :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign9

  FUNCTION Assign10(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE(real_t)            :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign10

  FUNCTION Assign11(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER :: res
    CHARACTER(LEN=*)        :: key
    TYPE(int_t)             :: val
    !------------------------------------------------------------------------!
    INTENT(IN)              :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign11

  FUNCTION Assign12(key, val) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER            :: res, tmp
    CHARACTER(LEN=*)                   :: key
    REAL, DIMENSION(:,:,:,:,:), TARGET :: val
    !------------------------------------------------------------------------!
    INTENT(IN)                         :: key, val
    !------------------------------------------------------------------------!
    NULLIFY(res)
    CALL SetAttr(res, key, val)
  END FUNCTION Assign12

  FUNCTION Ref1(p) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE_DICT_REAL_P :: p
    TYPE(real_t)     :: res
    !------------------------------------------------------------------------!
    res%p => p
  END FUNCTION Ref1

  FUNCTION Ref2(p) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE_DICT_INT_P :: p
    TYPE(int_t)     :: res
    !------------------------------------------------------------------------!
    res%p => p
  END FUNCTION Ref2

END MODULE common_dict
