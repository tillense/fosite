!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
    TYPE(Dict_TYP), POINTER &
                       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER            :: err,method,ShowBut
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 5)

    CALL GetAttr(config, "method", method)
    CALL this%InitTimedisc(Mesh,Physics,config,IO,method,ODEsolver_name)

!CDIR IEXPAND
    SELECT CASE(this%GetOrder())
    CASE(5)
       !set number of coefficients
       this%m = 7
       ! allocate memory
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM,this%m), &
                this%b_high(this%m),&
                this%b_low(this%m),&
                this%c(this%m),&
                this%a(this%m,this%m), &
                STAT = err)
       IF (err.NE.0) THEN
          CALL this%Error("timedisc_dormand_prince", "Unable to allocate memory.")
       END IF
       !set coefficient scheme of dormand_prince
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
       CALL this%Error("timedisc_dormand_prince","time order must be 5")
    END SELECT

    ! set RHS pointer to the first entry of the coeff field
    this%rhs => Mesh%RemapBounds(this%coeff(:,:,:,:,1))

    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL this%Error("timedisc_dormand_prince", &
         "error tolerance levels must be greater than 0")
    IF (this%tol_rel.GT.1.0) THEN
         CALL this%Warning("timedisc_dormand_prince", &
            "adaptive step size control disabled (tol_rel>1)",0)
    ELSE IF(this%tol_rel.GE.0.01) THEN
         CALL this%Warning("timedisc_dormand_prince", &
             "You chose a relatively high tol_rel (in comparison to order)",0)
    END IF
    CALL GetAttr(config, "ShowButcherTableau", ShowBut, 0)
    IF (ShowBut .EQ. 1) CALL this%ShowButcherTableau()
  END SUBROUTINE InitTimedisc_dormand_prince

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    CLASS(timedisc_dormand_prince) :: this
    !-----------------------------------------------------------------------!
!    DEALLOCATE(this%coeff,this%b_high,this%b_low,this%c,this%a)
    CALL this%timedisc_rkfehlberg%Finalize()
  END SUBROUTINE
END MODULE timedisc_dormand_prince_mod
