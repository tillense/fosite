!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_dormand_prince.f90                                       #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!!
!! \brief subroutines for Dormand-Prince method
!!
!! Reference:
!! - \cite dormand1980 Dormand, J. R.; Prince, P. J. (1980),
!!            "A family of embedded Runge-Kutta formulae",
!!            Journal of Computational and Applied Mathematics 6 (1): 19–26,
!!            doi:10.1016/0771-050X(80)90013-3
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_dormand_prince_mod
  USE timedisc_base_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE sources_base_mod
  USE timedisc_rkfehlberg_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(timedisc_rkfehlberg) :: timedisc_dormand_prince
  CONTAINS
    PROCEDURE :: InitTimedisc_dormand_prince
    PROCEDURE :: SetButcherTableau
    PROCEDURE :: Finalize
  END TYPE timedisc_dormand_prince
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Dormand-Prince method"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       timedisc_dormand_prince
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_dormand_prince(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_dormand_prince), INTENT(INOUT) :: this
    CLASS(mesh_base),               INTENT(INOUT) :: Mesh
    CLASS(physics_base),            INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER :: config,IO
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 5)

    ! set number of coefficients
    SELECT CASE(this%GetOrder())
    CASE(5)
       this%m = 7
    CASE DEFAULT
       CALL this%Error("InitTimedisc_dormand_prince","only order 5 supported")
    END SELECT

    CALL this%InitTimedisc(Mesh,Physics,config,IO,DORMAND_PRINCE,ODEsolver_name)
  END SUBROUTINE InitTimedisc_dormand_prince

  !> set coefficients for dormand_prince scheme
  SUBROUTINE SetButcherTableau(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_dormand_prince)   :: this
    !------------------------------------------------------------------------!
    SELECT CASE(this%GetOrder())
    CASE(5)
       this%b_high = (/ 35./384.,0.,500./1113.,125./192.,-2187./6784.,11./84.,0. /)
       this%b_low  = (/ 5179./57600., 0., 7571./16695., 393./640., &
                       -92097./339200., 187./2100., 1./40. /)
       this%c  = (/ 0., .2, 3./10., 4./5., 8./9., 1., 1. /)
       this%a  = TRANSPOSE(RESHAPE((/ &
                 0., 0., 0., 0., 0., 0., 0., &
                 .2, 0., 0., 0., 0., 0., 0., &
                 3./40., 9./40., 0., 0., 0., 0., 0., &
                 44./45., -56./15., 32./9., 0., 0., 0., 0., &
                 19372./6561., -25360./2187., 64448./6561., -212./729., 0., 0., 0.,&
                 9017./3168.,-355./33., 46732./5247.,49./176.,-5103./18656.,0.,0., &
                 35./384.,0.,500./1113.,125./192.,-2187./6784.,11./84.,0. /),&
                 (/this%m,this%m/)))
    CASE DEFAULT
      CALL this%Error("timedisc_dormand_prince::SetButcherTableau","only order 5 supported")
    END SELECT
  END SUBROUTINE SetButcherTableau

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    CLASS(timedisc_dormand_prince) :: this
    !-----------------------------------------------------------------------!
    CALL this%timedisc_rkfehlberg%Finalize()
  END SUBROUTINE
END MODULE timedisc_dormand_prince_mod
