!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: constants_common.f90                                              #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> \defgroup constants constants
!! \{
!! \brief Family of constants modules
!!
!! This is the family of constants modules.The generic interface routines are
!! defined in the module \link constants_generic \endlink. The basic constants
!! data type and common basic subroutines and functions are defined in
!! \link constants_common \endlink. Any other module of the family defines
!! physical constants and conversion factors for specific unit systems.
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief Basic physical constants module
!!
!! \extends common_types
!----------------------------------------------------------------------------!
MODULE constants_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName,         &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info,       &
       Warning_common => Warning, Error_common => Error
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetUnits, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetUnitsName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetConstantsRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetConstantsNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE ConstantsInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE ConstantsInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE ConstantsWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE ConstantsError, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> \name Public Attributes
  !!#### Physical constants in SI units
  DOUBLE PRECISION, PARAMETER :: &
     C  = 2.99792458E+08, & !< vacuum speed of light [m/s]
     GN = 6.6742E-11,     & !< gravitational constant [m^3/kg/s^2]
     KB = 1.3806505E-23,  & !< Boltzmann constant[J/K]
     NA = 6.022E+23,      & !< Avogadro constant [1/mol]
     SB = 5.6704E-8,      & !< Stefan-Boltzmann constant[W/m^2/K^4
     KE = 3.48E-02          !< electron scattering opacity [m^2/kg]
  !--------------------------------------------------------------------------!
  !> physical constants data structure
  !!
  !! This data type stores information on the unit system, the value
  !! of important physical constants and conversion factors to SI units.
  !!
  TYPE Constants_TYP
     TYPE(Common_TYP) :: units                      !< SI, natural, etc.
     !> \name Variables
     !!#### Named physical constants
     DOUBLE PRECISION :: C                          !< vacuum speed of light
     DOUBLE PRECISION :: GN                         !< gravitational constant
     DOUBLE PRECISION :: KB                         !< Boltzmann constant
     DOUBLE PRECISION :: NA                         !< Avogadro constant
     DOUBLE PRECISION :: SB                         !< Stefan-Boltzmann constant
     DOUBLE PRECISION :: RG                         !< universal gas constant
     DOUBLE PRECISION :: KE                         !< electron scattering opacity
     !> \name
     !!#### Conversion factors from SI to other units
     DOUBLE PRECISION :: cf_time                    !< time scale
     DOUBLE PRECISION :: cf_mass                    !< mass scale
     DOUBLE PRECISION :: cf_momentum                !< momentum scale
     DOUBLE PRECISION :: cf_energy                  !< energy scale
     DOUBLE PRECISION :: cf_power                   !< power scale
     DOUBLE PRECISION :: cf_temperature             !< temperature scale
     DOUBLE PRECISION :: cf_density                 !< density scale
     DOUBLE PRECISION :: cf_opacity                 !< opacity scale
  END TYPE Constants_TYP
  SAVE
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC ::                   &
       ! types
       Constants_TYP,         &
       ! constants
       C, GN, KB, NA, SB, KE, &
       ! methods
       InitConstants,         &
       GetType,               &
       GetName,               &
       GetRank,               &
       GetNumProcs,           &
       Initialized,           &
       Info,                  &
       Warning,               &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of common physical constants module
  !!
  !! Initializes the physical constants module with the given unit system.
  SUBROUTINE InitConstants(this,ut,un)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this !< \param [inout] this unit system
    INTEGER             :: ut   !< \param [in] ut unit system number
    CHARACTER(LEN=32)   :: un   !< \param [in] un unit system name
    !------------------------------------------------------------------------!
    INTENT(IN)          :: ut,un
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%units,ut,un)
  END SUBROUTINE InitConstants


  !> \public Destructor of common physical constants module
  SUBROUTINE CloseConstants(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(INOUT) :: this !< \param [inout] this unit system
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%units)
  END SUBROUTINE CloseConstants


  !> \public Get the physical unit system number;
  !! overloads \b GetType from \link common_types::GetType \endlink
  !! \return physical units number
  PURE FUNCTION GetUnits(this) RESULT(ut)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this !< \param [in] this unit system
    INTEGER                         :: ut
    !------------------------------------------------------------------------!
    ut = GetType_common(this%units)
  END FUNCTION GetUnits


  !> \public Get the physical unit system name;
  !! overloads \b GetName from \link common_types::GetName \endlink
  !! \return physical unit system name
  PURE FUNCTION GetUnitsName(this) RESULT(un)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this !< \param [in] this unit system
    CHARACTER(LEN=32)               :: un
    !------------------------------------------------------------------------!
    un = GetName_common(this%units)
  END FUNCTION GetUnitsName


  !> \public Get the MPI rank;
  !! overloads \b GetRank from \link common_types::GetRank \endlink
  !! \return MPI rank
  PURE FUNCTION GetConstantsRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this !< \param [in] this unit system
    INTEGER                         :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%units)
  END FUNCTION GetConstantsRank


  !> \public Get the total number of MPI processes;
  !! overloads \b GetNumProcs from \link common_types::GetNumProcs \endlink
  !! \return number of MPI processes
  PURE FUNCTION GetConstantsNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this !< \param [in] this unit system
    INTEGER                         :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%units)
  END FUNCTION GetConstantsNumProcs


  !> \public Query initialization status;
  !! overloads \b Initialized from \link common_types::Initialized \endlink
  PURE FUNCTION ConstantsInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP), INTENT(IN) :: this !< \param [in] this unit system
    LOGICAL                         :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%units)
  END FUNCTION ConstantsInitialized

 
  !> \public Print information on standard output;
  !! overloads \b Info from \link common_types::Info \endlink
  SUBROUTINE ConstantsInfo(this,msg,rank,node_info,tostderr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this      !< \param [in] this unit system
    CHARACTER(LEN=*)    :: msg       !< \param [in] msg info message
    INTEGER, OPTIONAL   :: rank      !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL   :: node_info !< \param [in] node_info enable rank output
    LOGICAL, OPTIONAL   :: tostderr  !< \param [in] tostderr enable STDERR output
    !------------------------------------------------------------------------!
    INTENT(IN)          :: this,msg,rank,node_info,tostderr
    !------------------------------------------------------------------------!
    CALL Info_common(this%units,msg)
  END SUBROUTINE ConstantsInfo


  !> \public Print warning message on standard error;
  !! overloads \b Warning from \link common_types::Warning \endlink
  SUBROUTINE ConstantsWarning(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this    !< \param [in] this unit system
    CHARACTER(LEN=*)    :: modproc !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*)    :: msg     !< \param [in] msg warning message
    INTEGER, OPTIONAL   :: rank    !< \param [in] rank MPI rank
    !------------------------------------------------------------------------!
    INTENT(IN)          :: this,modproc,msg,rank
    !------------------------------------------------------------------------!
    CALL Warning_common(this%units,modproc,msg,rank)
  END SUBROUTINE ConstantsWarning


  !> \public Print error message on standard error and terminate the program;
  !! overloads \b Error from \link common_types::Error \endlink
  SUBROUTINE ConstantsError(this,modproc,msg,rank,node_info)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Constants_TYP) :: this      !< \param [in] this unit system
    CHARACTER(LEN=*)    :: modproc   !< \param [in] modproc name of module procedure
    CHARACTER(LEN=*)    :: msg       !< \param [in] msg warning message
    INTEGER, OPTIONAL   :: rank      !< \param [in] rank MPI rank
    LOGICAL, OPTIONAL   :: node_info !< \param [in] node_info enable rank output
    !------------------------------------------------------------------------!
    INTENT(IN)          :: this,modproc,msg,rank,node_info
    !------------------------------------------------------------------------!
    CALL Error_common(this%units,modproc,msg,rank,node_info)
  END SUBROUTINE ConstantsError

END MODULE constants_common
