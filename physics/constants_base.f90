!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: constants_generic.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2008,2011                                              #
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
!> \defgroup constants Constants
!! \{
!! \brief Module class for system of units and natural constants
!! \}
!----------------------------------------------------------------------------!
!> \addtogroup constants
!! \key{units,INTEGER,physical unit system}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief generic module for units and physical constants
!!
!! \ingroup constants
!----------------------------------------------------------------------------!
MODULE constants_base_mod
  USE logging_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: SI          = 1
  INTEGER, PARAMETER :: CGS         = 2
  INTEGER, PARAMETER :: GEOMETRICAL = 3
  !> physical constants data structure
  !!
  !! This data type stores information on the unit system, the value
  !! of important physical constants and conversion factors to SI units.
  !!
  TYPE, EXTENDS(logging_base) ::  constants_base
     !> \name Variables
     !!#### Named physical constants
     REAL :: C                          !< vacuum speed of light
     REAL :: GN                         !< gravitational constant
     REAL :: KB                         !< Boltzmann constant
     REAL :: NA                         !< Avogadro constant
     REAL :: SB                         !< Stefan-Boltzmann constant
     REAL :: RG                         !< universal gas constant
     REAL :: KE                         !< electron scattering opacity
     !> \name
     !!#### Conversion factors from SI to other units
     REAL :: cf_time                    !< time scale
     REAL :: cf_mass                    !< mass scale
     REAL :: cf_momentum                !< momentum scale
     REAL :: cf_energy                  !< energy scale
     REAL :: cf_power                   !< power scale
     REAL :: cf_temperature             !< temperature scale
     REAL :: cf_density                 !< density scale
     REAL :: cf_opacity                 !< opacity scale
  CONTAINS
    PROCEDURE :: InitConstants
  END TYPE constants_base
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       constants_base, &
       ! constant flags
       SI, CGS, GEOMETRICAL
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitConstants(this,units,units_name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(constants_base),INTENT(INOUT) :: this
    INTEGER                             :: units
    CHARACTER(LEN=*)                    :: units_name
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: units
    !------------------------------------------------------------------------!
    CALL this%InitLogging(units,units_name)

    ! derived constants
    this%RG = this%KB * this%NA         ![J/mol/K]    universal gas constant !

    ! print some information
    CALL this%Info(" CONSTANTS> physical units:    " // TRIM(this%GetName()))
  END SUBROUTINE InitConstants

END MODULE constants_base_mod
