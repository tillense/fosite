!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: constants_SI.f90                                                  #
!#                                                                           #
!# Copyright (C) 2007-2016                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief module for SI units and physical constants
!!
!! \extends constants_common
!! \ingroup constants
!----------------------------------------------------------------------------!
MODULE constants_SI_mod
  USE constants_base_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(constants_base) :: constants_SI
  CONTAINS
    PROCEDURE :: InitConstants_SI
  END TYPE
  CHARACTER(LEN=32), PARAMETER :: units_name = 'SI'
  !--------------------------------------------------------------------------!
  ! basic numerical constants
  !!#### Physical constants in SI units
  REAL, PARAMETER :: &
     C  = 2.99792458E+08, &   !< vacuum speed of light       [m/s]
     GN = 6.6742E-11, &       !< gravitational constant      [m^3/kg/s^2]
     KB = 1.3806505E-23, &    !< Boltzmann constant          [J/K]
     NA = 6.022E+23, &        !< Avogadro constant           [1/mol]
     SB = 5.6704E-8, &        !< Stefan-Boltzmann constant   [W/m^2/K^4]
     KE = 3.48E-02, &         !< electron scattering opacity [m^2/kg]
     AU = 1.49597870691D+11, &!< astronomical unit           [m]
     MSUN     = 1.9885D+30, & !< sun mass                    [kg]
     MJUPITER = 1.89819D+27, &!< jupiter mass                [kg]
     MEARTH   = 5.9723D+24, & !< earth mass                  [kg]
     RSUN     = 6.957E+8, &   !< sun radius                  [m]
     RJUPITER = 6.9911E+7, &  !< jupiter radius              [m]
     REARTH   = 6.371E+6, &   !< earth radius                [m]
     DAY      = 8.64E+4, &    !< length of a day             [s]
     PI = 3.1415926535897932384626433832795028842
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       constants_SI, &
       C, GN, KB, NA, SB, KE, AU, &
       MSUN, MJUPITER, MEARTH, RSUN, RJUPITER, REARTH, DAY, PI
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of physical constants module using SI units
  SUBROUTINE InitConstants_SI(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(constants_SI), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! assign numerical values of physical constants in SI units;
    ! (C, GN, etc. are defined in constants_common)
    this%C  = C
    this%GN = GN
    this%KB = KB
    this%NA = NA
    this%SB = SB
    this%KE = KE
    this%AU = AU
    this%MSUN     = MSUN
    this%MJUPITER = MJUPITER
    this%MEARTH   = MEARTH
    this%RSUN     = RSUN
    this%RJUPITER = RJUPITER
    this%REARTH   = REARTH
    this%DAY      = DAY
    this%PI       = PI
    ! conversion factors to SI units are unity
    this%cf_time        = 1.0
    this%cf_mass        = 1.0
    this%cf_momentum    = 1.0
    this%cf_energy      = 1.0
    this%cf_power       = 1.0
    this%cf_temperature = 1.0
    this%cf_density     = 1.0
    this%cf_opacity     = 1.0

    CALL this%InitConstants(SI,units_name)
  END SUBROUTINE InitConstants_SI

END MODULE constants_SI_mod
