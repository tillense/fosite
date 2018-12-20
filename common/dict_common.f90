!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: dict_common.f90                                                   #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de<                               #
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
!> basic dict module
!----------------------------------------------------------------------------!
MODULE dict_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetDictType, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetDictName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetDictRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetDictNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE DictInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE DictInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE DictWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE DictError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  TYPE Dict2_TYP
     TYPE(Common_TYP)       :: parent              
  END TYPE Dict2_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Dict2_TYP, &
       ! methods
       InitDict, &
       CloseDict, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       !GetErrorMap, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitDict(this,type,name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP)  :: this
    INTEGER           :: type
    CHARACTER(LEN=32) :: name
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: type,name
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%parent,type,name)
  END SUBROUTINE InitDict


  PURE FUNCTION GetDictType(this) RESULT(ap)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN) :: this
    INTEGER :: ap
    !------------------------------------------------------------------------!
    ap = GetType_common(this%parent)
  END FUNCTION GetDictType


  PURE FUNCTION GetDictName(this) RESULT(an)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: an
    !------------------------------------------------------------------------!
    an = GetName_common(this%parent)
  END FUNCTION GetDictName


  PURE FUNCTION GetDictRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%parent)
  END FUNCTION GetDictRank

  PURE FUNCTION GetDictNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%parent)
  END FUNCTION GetDictNumProcs


!  PURE FUNCTION GetErrorMap(this, error) RESULT(c)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(Fosite_TYP), INTENT(IN) :: this
!    INTEGER, INTENT(IN):: error
!    CHARACTER(LEN=1) :: c
!    !------------------------------------------------------------------------!
!    c = this%errormap(error)
!  END FUNCTION GetErrorMap


  PURE FUNCTION DictInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%parent)
  END FUNCTION DictInitialized

 
  SUBROUTINE DictInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN)  :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%parent,msg)
  END SUBROUTINE DictInfo


  !> \public Print warning message on standard error;
  !! overloads \b Warning from \link common_types::Warning \endlink
  SUBROUTINE DictWarning(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN)   :: this    !< \param [in] this dictionary type
    CHARACTER(LEN=*),  INTENT(IN) :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*),  INTENT(IN) :: msg     !< \param [in] msg warning message
    INTEGER, OPTIONAL, INTENT(IN) :: rank    !< \param [in] rank MPI rank
    !------------------------------------------------------------------------!
    CALL Warning_common(this%parent,modproc,msg,rank)
  END SUBROUTINE DictWarning


  SUBROUTINE DictError(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP), INTENT(IN)  :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, OPTIONAL, INTENT(IN) :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%parent,modproc,msg,rank)
  END SUBROUTINE DictError


  SUBROUTINE CloseDict(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict2_TYP) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%parent)
  END SUBROUTINE CloseDict


END MODULE dict_common
