!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
!! Reference:
!! - \cite gottlieb2011 Gottlieb et. al (2011): Strong stability preserving runge-kutta
!! and multistep time discretization (book)
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
    TYPE(Dict_TYP), POINTER &
                       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER            :: err,method,ShowBut
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 5)

    CALL GetAttr(config, "method", method)
    CALL this%InitTimedisc(Mesh,Physics,config,IO,method,ODEsolver_name)

!NEC$ IEXPAND
    SELECT CASE(this%GetOrder())
    CASE(3)
       CALL this%Warning("timedisc_ssprk", &
         "This 3rd order embedded SSP RK scheme has been constructed from a second " // &
         "third order SSP RK scheme by hand! It seems to work, but i am not sure, that one " // &
         "is allowed to do so. An optimal embedded third order SSP RK scheme is described " // &
         "in chapter 6.3 of the main reference (see above), but still needs to be translated into " // &
         "a butchers tableau. => Better use the 5th order scheme!")
       !set number of coefficients
       this%m = 3
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM,this%m), &
                this%b_high(this%m),&
                this%b_low(this%m),&
                this%c(this%m),&
                this%a(this%m,this%m), &
                STAT = err)
       IF (err.NE.0) THEN
          CALL this%Error("timedisc_ssprk", "Unable to allocate memory.")
       END IF
       this%b_high = (/ 1./6., 1./6., 2./3. /)
       this%b_low  = (/ 0.5,   0.5,   0.    /)
       this%c  = (/ 0.,    1.,    0.5   /)
       this%a  = TRANSPOSE(RESHAPE((/ &
                 0.0,  0.0,  0.0, &
                 1.0,  0.0,  0.0, &
                 0.25, 0.25, 0.0 /),(/this%m,this%m/)))
    CASE(5) 
       ! RK(5,4) SSP(5,3) 
       ! Reference: Colin Barr Macdonald: Constructing High-Order Runge_lutta Methods with Embedded 
       ! Strong-Stability-PReserving Pairs (2003)
       !set number of coefficients
       this%m = 5
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM,this%m), &
                this%b_high(this%m),&
                this%b_low(this%m),&
                this%c(this%m),&
                this%a(this%m,this%m), &
                STAT = err)
       IF (err.NE.0) THEN
          CALL this%Error("timedisc_ssprk", "Unable to allocate memory.")
       END IF
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
       CALL this%Error("timedisc_ssprk","time order must be 3 or 5")
    END SELECT

    ! set RHS pointer to the first entry of the coeff field
    this%rhs => Mesh%RemapBounds(this%coeff(:,:,:,:,1))

    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL this%Error("timedisc_ssprk", &
         "error tolerance levels must be greater than 0")
    IF (this%tol_rel.GT.1.0) THEN
         CALL this%Warning("timedisc_ssprk", &
            "adaptive step size control disabled (tol_rel>1)",0)
    ELSE IF(this%tol_rel.GE.0.01) THEN
         CALL this%Warning("timedisc_ssprk", &
             "You chose a relatively high tol_rel (in comparison to order)",0)
    END IF
    CALL GetAttr(config, "ShowButcherTableau", ShowBut, 0)
    IF (ShowBut .EQ. 1) CALL this%ShowButcherTableau()
  END SUBROUTINE InitTimedisc_ssprk

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_ssprk) :: this
    !------------------------------------------------------------------------!
    CALL this%timedisc_rkfehlberg%Finalize()
  END SUBROUTINE

END MODULE timedisc_ssprk_mod
