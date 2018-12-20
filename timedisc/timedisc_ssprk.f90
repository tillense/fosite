!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_ssprk.f90                                                #
!#                                                                           #
!# Copyright (C) 2013                                                        #
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
!> \author Manuel Jung
!!
!! \brief subroutines for strong stability preserving (SSP) Runge Kutta methods
!!
!! References:
!! - \cite macdonald2003 : Constructing High-Order Runge-Kutta Methods with Embedded
!!   Strong-Stability-Preserving Pairs (PhD thesis)
!! - \cite gottlieb2011 Gottlieb et. al (2011): Strong stability preserving runge-kutta
!!   and multistep time discretization (book)
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_ssprk_mod
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
  TYPE, EXTENDS (timedisc_rkfehlberg) :: timedisc_ssprk
  CONTAINS
    PROCEDURE :: InitTimedisc_ssprk
    PROCEDURE :: SetButcherTableau
    PROCEDURE :: Finalize
  END TYPE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "SSP Runge-Kutta method"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       timedisc_ssprk
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_ssprk(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_ssprk), INTENT(INOUT) :: this
    CLASS(Mesh_base),      INTENT(INOUT) :: Mesh
    CLASS(Physics_base),   INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER :: config,IO
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 5)

    ! set number of coefficients
    SELECT CASE(this%GetOrder())
    CASE(3)
       this%m = 3
    CASE(5)
       this%m = 5
    CASE DEFAULT
       CALL this%Error("timedisc_ssprk::InitTimedisc_ssprk","time order must be 3 or 5")
    END SELECT

    CALL this%InitTimedisc(Mesh,Physics,config,IO,SSPRK,ODEsolver_name)
  END SUBROUTINE InitTimedisc_ssprk

  !> set coefficients for strong stability preserving schemes
  SUBROUTINE SetButcherTableau(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_ssprk)   :: this
    !------------------------------------------------------------------------!
    SELECT CASE(this%GetOrder())
    CASE(3) ! SSPRK order 3
      CALL this%Warning("timedisc_ssprk", &
        "This 3rd order embedded SSP RK scheme has been constructed from a second " // &
        "third order SSP RK scheme by hand! It seems to work, but I am not sure, that one " // &
        "is allowed to do so. An optimal embedded third order SSP RK scheme is described " // &
        "in chapter 6.3 of the main reference (see above), but still needs to be translated into " // &
        "a butchers tableau. => Better use the 5th order scheme!")
      this%b_high = (/ 1./6., 1./6., 2./3. /)
      this%b_low  = (/ 0.5,   0.5,   0.    /)
      this%c  = (/ 0.,    1.,    0.5   /)
      this%a  = TRANSPOSE(RESHAPE((/ &
                0.0,  0.0,  0.0, &
                1.0,  0.0,  0.0, &
                0.25, 0.25, 0.0 /),(/this%m,this%m/)))
    CASE(5) ! SSPRK order 5
      this%b_high = (/ 0.17279, 0.094505, 0.12947, 0.29899, 0.30424 /)
      this%b_low  = (/ 0.12293, 0.31981, -0.15316, 0.31887, 0.39155 /)
      this%c  = (/ 0., 0.36717, 0.58522, 0.44156, 0.8464 /)
      this%a  = TRANSPOSE(RESHAPE((/ &
                0., 0., 0., 0., 0., &
                0.36717, 0., 0., 0., 0., &
                0.26802, 0.3172, 0., 0., 0., &
                0.11606, 0.13735, 0.18816, 0., 0., &
                0.11212, 0.13269, 0.18178, 0.4198, 0. /),(/this%m,this%m/)))
    CASE DEFAULT
      CALL this%Error("timedisc_ssprk::SetButcherTableau","only order 3 or 5 supported")
    END SELECT
  END SUBROUTINE SetButcherTableau

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_ssprk) :: this
    !------------------------------------------------------------------------!
    CALL this%timedisc_rkfehlberg%Finalize()
  END SUBROUTINE

END MODULE timedisc_ssprk_mod
