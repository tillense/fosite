!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: constants_geometrical.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2016                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \author Tobias Illenseer
!! \author Manuel Jung
!!
!! \brief module for geometrical units and physical constants
!!
!! \extends constants_common
!! \ingroup constants
!----------------------------------------------------------------------------!
MODULE constants_geometrical_mod
  USE constants_base_mod
  USE constants_SI_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(constants_SI) :: constants_geometrical
  CONTAINS
    PROCEDURE :: InitConstants_geometrical
  END TYPE
  CHARACTER(LEN=32), PARAMETER :: units_name = 'geometrical'
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       constants_geometrical
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of physical constants module using geometrical units
  SUBROUTINE InitConstants_geometrical(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(constants_geometrical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    REAl :: C, GN, KB, NA, SB, KE
    !------------------------------------------------------------------------!
    CALL this%InitConstants(GEOMETRICAL,units_name)
    ! numerical values of physical constants in geometrical units (c = G = 1)
    this%C  = 1.0                                   ! lightspeed             !
    this%GN = 1.0                                   ! Newtons grav. constant !
    this%KB = 1.0                                   ! Boltzmann constant     !
    this%NA = NA  ![1/mol]                            Avogadro constant      !
    ! factors for conversion to SI units
    ! time and mass have unit of length i.e. metre

    C  = this%constants_SI%C
    GN = this%constants_SI%GN
    KB = this%constants_SI%KB
    NA = this%constants_SI%NA
    SB = this%constants_SI%SB
    KE = this%constants_SI%KE
    this%cf_time = C                 ! i.e. 1 sec [SI] = 3e8 m [geometrical] !
    this%cf_mass = (GN/C)/C
    this%cf_momentum = this%cf_mass/C
    this%cf_energy = this%cf_momentum/C
    this%cf_power = this%cf_energy/C
    this%cf_temperature = KB * this%cf_energy
    this%cf_density = this%cf_mass
    ! some scaled constants
    this%KE = KE*C/(GN/C)               !  [m]        electron scat. opacity !
  END SUBROUTINE InitConstants_geometrical

END MODULE constants_geometrical_mod
