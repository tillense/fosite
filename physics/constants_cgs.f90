!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: constants_cgs.f03                                                 #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
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
!> \author Jannes Klee
!!
!! \brief module for cgs units and physical constants
!!
!! \extends constants_common
!! \ingroup constants
!----------------------------------------------------------------------------!
MODULE constants_cgs_mod
  USE constants_base_mod
  USE constants_SI_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(constants_SI) :: constants_cgs
  CONTAINS
    PROCEDURE :: InitConstants_cgs
  END TYPE
  CHARACTER(LEN=32), PARAMETER :: units_name = 'cgs'
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! classes
       constants_cgs
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of physical constants module using cgs units
  SUBROUTINE InitConstants_cgs(this,units)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(constants_cgs), INTENT(INOUT) :: this
    INTEGER,              INTENT(IN)    :: units
    !------------------------------------------------------------------------!
    ! assign numerical values of physical constants in cgs units;
    ! (C, GN, etc. are defined in constants_common)
    this%C  = C*1.0e2           !< vacuum speed of light [cm/s]
    this%GN = GN*1.0e3          !< gravitational constant [cm^3/g/s^2]
    this%KB = KB*1.0e7          !< Boltzmann constant [erg/K]
    this%NA = NA                !< Avogadro constant [1/mol]
    this%SB = SB*1.0e3          !< Stefan-Boltzmann constant [erg/s/cm^2/K^4]
    this%KE = KE*1.0e1          !< electron scattering opacity [cm^2/g]
    ! conversion factors to from SI to cgs
    this%cf_time = 1.0          !< [s] -> [s]
    this%cf_mass = 1.0e3        !< [kg] -> [g]
    this%cf_momentum = 1.0e5    !< [kg*m/s] -> [g*cm/s]
    this%cf_energy = 1.0e7      !< [kg*m/s] -> [g*cm/s]
    this%cf_power = 1.0e7       !< [kg*m/s^2] -> [g*cm/s^2]
    this%cf_temperature = 1.0   !< [K] -> [K]
    this%cf_density = 1.0e-3    !< [kg/m^3] -> [g/cm^3]
    this%cf_opacity = 1.0e1     !< [m^2/kg] -> [cm^2/kg]

    CALL this%InitConstants(cgs,units_name)
  END SUBROUTINE InitConstants_cgs

END MODULE constants_cgs_mod
