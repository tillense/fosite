!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_rkfehlberg .f90                                          #
!#                                                                           #
!# Copyright (C) 2011,2014                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Lars Boesch
!!
!! \brief subroutines for Runge-Kutta Fehlberg method
!!
!! Reference:
!! - \cite engeln2006 G.Engeln-Müllges & F.Reutter (Book)
!! - \cite fehlberg1969 Fehlberg, E. (1969). Low-order classical Runge-Kutta
!!      formulas with stepsize control and their application to some heat
!!      transfer problems.
!!
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_rkfehlberg_mod
  USE timedisc_base_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE timedisc_modeuler_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE sources_base_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE coeff_type
    CLASS(marray_compound), POINTER :: p => null()
  END TYPE coeff_type
  TYPE, EXTENDS (timedisc_modeuler) :: timedisc_rkfehlberg
     TYPE(coeff_type), DIMENSION(:), ALLOCATABLE :: coeff
     REAL, DIMENSION(:), POINTER        :: b_low,b_high,c   !<    needed by
     REAL, DIMENSION(:,:), POINTER      :: a                !<    embedded RK
  CONTAINS
    PROCEDURE :: InitTimedisc
    PROCEDURE :: InitTimedisc_rkfehlberg
    PROCEDURE :: Finalize
    PROCEDURE :: SolveODE
    PROCEDURE :: ComputeCVar_rkfehlberg
    GENERIC   :: ComputeCVar => ComputeCVar_rkfehlberg
    PROCEDURE :: SetButcherTableau
    PROCEDURE :: ShowButcherTableau
  END TYPE timedisc_rkfehlberg
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Runge-Kutta Fehlberg"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       coeff_type, &
       timedisc_rkfehlberg
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc(this,Mesh,Physics,config,IO,ttype,tname)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_rkfehlberg), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(INOUT) :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    TYPE(Dict_TYP),       POINTER       :: config,IO
    INTEGER,              INTENT(IN)    :: ttype
    CHARACTER(LEN=32),    INTENT(IN)    :: tname
    !------------------------------------------------------------------------!
    INTEGER :: err,k,ShowBT
    !------------------------------------------------------------------------!
    CALL this%timedisc_modeuler%InitTimedisc(Mesh,Physics,config,IO,ttype,tname)

    ALLOCATE(this%b_high(this%m),this%b_low(this%m),&
             this%c(this%m),this%a(this%m,this%m), &
             this%coeff(this%m), &
             STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("timedisc_rkfehlberg::InitTimedisc", "Unable to allocate memory.")
    END IF

    ! init Butcher tableau
    CALL this%SetButcherTableau()
    CALL GetAttr(config, "ShowButcherTableau", ShowBT, 0)
    IF (ShowBT.EQ.1) CALL this%ShowButcherTableau()

    ! init coefficient compounds
    this%coeff(1)%p => this%rhs
    DO k=2,this%m
      CALL Physics%new_statevector(this%coeff(k)%p,CONSERVATIVE)
      this%coeff(k)%p%data1d(:) = 0.0
    END DO

    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL this%Error("timedisc_rkfehlberg::InitTimedisc", &
         "error tolerance levels must be greater than 0")

    IF (this%tol_rel.GT.1.0) THEN
         CALL this%Warning("timedisc_rkfehlberg::InitTimedisc", &
            "adaptive step size control disabled (tol_rel>1)",0)
    ELSE IF(this%tol_rel.GE.0.01 .AND. this%order .GE. 5) THEN
         CALL this%Warning("timedisc_rkfehlberg::InitTimedisc", &
             "You chose a relatively high tol_rel (in comparison to order)",0)
    END IF

  END SUBROUTINE

  SUBROUTINE InitTimedisc_rkfehlberg(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg), INTENT(INOUT) :: this
    CLASS(Mesh_base),           INTENT(INOUT) :: Mesh
    CLASS(Physics_base),        INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER                   :: config,IO
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 5)

    ! set number of coefficients
    SELECT CASE(this%GetOrder())
    CASE(3)
      this%m = 3
    CASE(5)
      this%m = 6
    CASE DEFAULT
       CALL this%Error("InitTimedisc_rkfehlberg","time order must be 3 or 5")
    END SELECT

    CALL this%InitTimedisc(Mesh,Physics,config,IO,RK_FEHLBERG,ODEsolver_name)
  END SUBROUTINE InitTimedisc_rkfehlberg


  SUBROUTINE SolveODE(this,Mesh,Physics,Sources,Fluxes,time,dt,err)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg), INTENT(INOUT) :: this
    CLASS(Mesh_base),           INTENT(IN)    :: Mesh
    CLASS(Physics_base),        INTENT(INOUT) :: Physics
    CLASS(Fluxes_base),         INTENT(INOUT) :: Fluxes
    CLASS(sources_base),        POINTER       :: Sources
    REAL,                       INTENT(IN)    :: time
    REAL,                       INTENT(INOUT) :: dt,err
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,l,m
    REAL               :: t
    !------------------------------------------------------------------------!
    t = time
    ! ATTENTION: rhs and therefore coeff(1)%p should already be up to date.
    ! In the very first step this is done in FirstStepFosite, any later step
    ! must update coeff(1)%p at the end of the whole time step update. This is
    ! done in timedisc_base::AcceptTimestep or timedisc_base::RejectTimestep.
    ! Don't forget: this%rhs => this%coeff(1)%p.
!NEC$ SHORTLOOP
    DO m=2,this%m
       ! time step update of cell mean values
       CALL this%ComputeCVar_rkfehlberg(Mesh,Physics,Fluxes,dt,m,this%coeff,this%cvar,this%ctmp)
       ! compute right-hand-side
       ! coeff_m is k_m/dt from Butcher tableau
       CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,t+this%c(m)*dt,dt,&
                           this%ptmp,this%ctmp,CHECK_NOTHING,this%coeff(m)%p)
    END DO

    !reset ctmp
    this%ctmp = this%cvar
!NEC$ SHORTLOOP
    DO m=1,this%m
        ! compute two solutions with different numerical orders
        ! y_n+1 = y_n + SUM(b_i*k_i) = y_n + SUM(b_i*dt*coeff_i)
        this%ctmp%data1d(:) = this%ctmp%data1d(:) &
                           - this%b_low(m)*dt*this%coeff(m)%p%data1d(:)
        this%cvar%data1d(:) = this%cvar%data1d(:) &
                           - this%b_high(m)*dt*this%coeff(m)%p%data1d(:)
    END DO

    ! at the boundary the this%rhs contains the boundary fluxes
    ! (see subroutine ComputeRHS_modeuler )
!NEC$ SHORTLOOP
    DO m=1,this%m
!NEC$ SHORTLOOP
      DO l=1,Physics%VNUM
      IF(Mesh%INUM.GT.1) THEN
        ! western and eastern
!NEC$ IVDEP
        DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
          DO j=Mesh%JMIN,Mesh%JMAX
             Fluxes%bxflux(j,k,1,l) = Fluxes%bxflux(j,k,1,l) &
                               - dt*this%b_high(m)*this%coeff(m)%p%data4d(Mesh%IMIN-Mesh%IP1,j,k,l)
             Fluxes%bxflux(j,k,2,l) = Fluxes%bxflux(j,k,2,l) &
                               - dt*this%b_high(m)*this%coeff(m)%p%data4d(Mesh%IMAX+Mesh%IP1,j,k,l)
          END DO
        END DO
      END IF
        ! southern and northern
      IF(Mesh%JNUM.GT.1) THEN
!NEC$ IVDEP
        DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN,Mesh%IMAX
            Fluxes%byflux(k,i,1,l) = Fluxes%byflux(k,i,1,l) &
                              - dt*this%b_high(m)*this%coeff(m)%p%data4d(i,Mesh%JMIN-Mesh%JP1,k,l)
            Fluxes%byflux(k,i,2,l) = Fluxes%byflux(k,i,2,l) &
                              - dt*this%b_high(m)*this%coeff(m)%p%data4d(i,Mesh%JMAX+Mesh%JP1,k,l)
          END DO
        END DO
      END IF
        ! bottom and top
      IF(Mesh%KNUM.GT.1) THEN
!NEC$ IVDEP
        DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN,Mesh%IMAX
            Fluxes%bzflux(i,j,1,l) = Fluxes%bzflux(i,j,1,l) &
                              - dt*this%b_high(m)*this%coeff(m)%p%data4d(i,j,Mesh%KMIN-Mesh%KP1,l)
            Fluxes%bzflux(i,j,2,l) = Fluxes%bzflux(i,j,2,l) &
                              - dt*this%b_high(m)*this%coeff(m)%p%data4d(i,j,Mesh%KMAX+Mesh%KP1,l)
          END DO
        END DO
      END IF
      END DO
    END DO

    ! compute an error estimate based on the two independent numerical
    ! solutions err and dt
    IF (this%tol_rel.LE.1.0) THEN
      err = this%ComputeError(Mesh,Physics,this%cvar,this%ctmp)
      dt = this%AdjustTimestep(err,dt)
    END IF

  END SUBROUTINE SolveODE

  !> \private perfroms the time step update using the RHS
  !!
  !! This subroutine computes new conservative variables
  !! according to the update formula:
  !! \f[
  !!     y_n^{(i)} = y_n + dt \sum_{j=1}^{m-1} a_{ij} coeff_j
  !! \f]
  SUBROUTINE ComputeCVar_rkfehlberg(this,Mesh,Physics,Fluxes,dt,m,coeff,cvar,cnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg),INTENT(INOUT) :: this
    CLASS(Mesh_base)          ,INTENT(IN)    :: Mesh
    CLASS(Physics_base)       ,INTENT(INOUT) :: Physics
    CLASS(Fluxes_base)        ,INTENT(INOUT) :: Fluxes
    REAL                      ,INTENT(IN)    :: dt
    INTEGER                   ,INTENT(IN)    :: m
    TYPE(coeff_type), DIMENSION(m), INTENT(IN) :: coeff
    CLASS(marray_compound)    ,INTENT(INOUT) :: cvar,cnew
    !------------------------------------------------------------------------!
    INTEGER :: mm
    !------------------------------------------------------------------------!
    ! compute y + SUM(a_mi*k_i) = y + SUM(a_mi*dt*rhs_i)  : i in [1,m-1]
    ! see Runge-Kutta method (Butcher tableau)
    IF (m.EQ.1) THEN
      cnew = cvar
    ELSE
      cnew%data1d(:) = cvar%data1d(:) - this%a(m,1)*dt*coeff(1)%p%data1d(:)
!NEC$ SHORTLOOP
      DO mm=2,m-1
        cnew%data1d(:) = cnew%data1d(:) - this%a(m,mm)*dt*coeff(mm)%p%data1d(:)
      END DO
    END IF
  END SUBROUTINE ComputeCVar_rkfehlberg

  !> set coefficients for RK-Fehlberg schemes
  SUBROUTINE SetButcherTableau(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_rkfehlberg)   :: this
    !------------------------------------------------------------------------!
    SELECT CASE(this%GetOrder())
    CASE(3) ! RK-Fehlberg rkf23
      this%b_high = (/ 1.0/6.0, 2.0/3.0, 1.0/6.0 /)
      this%b_low  = (/ 0.0, 1.0, 0.0 /)
      this%c  = (/ 0.0, 0.5, 1.0 /)
      this%a  = TRANSPOSE(RESHAPE((/ &
                0.0, 0.0, 0.0, &
                0.5, 0.0, 0.0, &
              -1.0, 2.0, 0.0 /),(/this%m,this%m/)))
    CASE(5) ! RK-Fehlberg rkf45
      this%b_high = (/ 16.0/135.0, 0.0, 6656.0/12825.0, &
                      28561.0/56430.0, -9.0/50.0, 2.0/55.0 /)
      this%b_low  = (/ 25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -0.2, 0.0 /)
      this%c  = (/ 0.0, 0.25, 3.0/8.0, 12.0/13.0, 1.0, 0.5 /)
      this%a  = TRANSPOSE(RESHAPE((/ &
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                0.25, 0.0, 0.0, 0.0, 0.0, 0.0, &
                3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, 0.0, &
                1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0, 0.0, &
                439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0, 0.0, &
                -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0.0/),&
                (/this%m,this%m/)))
    CASE DEFAULT
       CALL this%Error("timedisc_rkfehlberg::SetButcherTableau","only order 3 or 5 supported")
    END SELECT
  END SUBROUTINE SetButcherTableau


  SUBROUTINE ShowButcherTableau(this)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg) :: this
    !------------------------------------------------------------------------!
    INTEGER            :: i
    CHARACTER(LEN=64)  :: sformat
    CHARACTER(LEN=128)  :: buffer
    !------------------------------------------------------------------------!
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL this%Info(repeat("+",(this%m+1)*12+4))
    CALL this%Info("Butcher tableau")
    DO i=1,this%m
      IF (i == 1) THEN
        WRITE (sformat,'(A)') '(ES12.3,A)'
        WRITE (buffer,fmt=sformat) this%c(i), " | "
      ELSE
        WRITE (sformat,'(A,I1,A)') '(ES12.3,A,',i-1,'(ES12.3))'
        WRITE (buffer,fmt=sformat) this%c(i), " | " , this%a(i,1:i-1)
      END IF
      CALL this%Info(buffer)
    END DO
    CALL this%Info(repeat("-",(this%m+1)*12+4))
    WRITE (sformat,'(A,I1,A)') '(A,',this%m,'(ES12.3))'
    WRITE (buffer,fmt=sformat) " high order  | ", this%b_high(:)
    CALL this%Info(buffer)
    WRITE (buffer,fmt=sformat) "  low order  | ", this%b_low(:)
    CALL this%Info(buffer)
    CALL this%Info(repeat("+",(this%m+1)*12+4))
    WRITE (buffer,fmt=sformat) "c_i-SUM(a_ij)| ", this%c(:)-SUM(this%a(:,:),DIM=2)
    CALL this%Info(buffer)
  END SUBROUTINE ShowButcherTableau


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg)   :: this
    !------------------------------------------------------------------------!
    INTEGER :: k
    !------------------------------------------------------------------------!
    DEALLOCATE(this%b_high,this%b_low,this%c,this%a)
    DO k=2,this%m
      DEALLOCATE(this%coeff(k)%p)
    END DO
    DEALLOCATE(this%coeff)
    CALL this%timedisc_modeuler%Finalize()
  END SUBROUTINE Finalize

END MODULE timedisc_rkfehlberg_mod
