!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_rkfehlberg_mod
  USE timedisc_base_mod
  USE mesh_base_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE timedisc_modeuler_mod!, &
!          AdjustTimestep_rkfehlberg => AdjustTimestep_modeuler, &
!          CalcTimestep_rkfehlberg => CalcTimestep_modeuler, &
!          ComputeError_rkfehlberg => ComputeError_modeuler, &
!          ComputeRHS_rkfehlberg => ComputeRHS_modeuler
  USE sources_base_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS (timedisc_modeuler) :: timedisc_rkfehlberg
     REAL, DIMENSION(:,:,:,:,:),POINTER :: coeff            !< coefficents
     REAL, DIMENSION(:), POINTER       :: b_low,b_high,c   !<    needed by
     REAL, DIMENSION(:,:), POINTER     :: a                !<    embedded RK
  CONTAINS
    PROCEDURE :: InitTimedisc_rkfehlberg
    PROCEDURE :: Finalize
    PROCEDURE :: SolveODE
    PROCEDURE :: ComputeCVar_rkfehlberg
    PROCEDURE :: ShowButcherTableau
  END TYPE timedisc_rkfehlberg
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Runge-Kutta Fehlberg"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       timedisc_rkfehlberg!, &
       ! methods 
!       InitTimedisc_rkfehlberg, &
!       CloseTimedisc_rkfehlberg, &
!       SolveODE_rkfehlberg, &
!       CalcTimestep_rkfehlberg, &
!       ComputeCVar_rkfehlberg, &
!       ComputeRHS_rkfehlberg, &
!       ComputeError_rkfehlberg, &
!       ShowButcherTableau, &
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_rkfehlberg(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg), INTENT(INOUT) :: this
    CLASS(Mesh_base),           INTENT(INOUT) :: Mesh
    CLASS(Physics_base),        INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER                   :: config,IO
    !------------------------------------------------------------------------!
    INTEGER            :: err,method,ShowBut
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 5)

 !   CALL GetAttr(config, "method", method)
    CALL this%InitTimedisc(Mesh,Physics,config,IO,RK_FEHLBERG,ODEsolver_name)

!NEC$ IEXPAND
    ! set number of coefficients
    SELECT CASE(this%GetOrder())
    CASE(3)
       this%m = 3
    CASE(5)
       this%m = 6
    CASE DEFAULT
       CALL this%Error("InitTimedisc_rkfehlberg","time order must be 3 or 5")
    END SELECT

    ! allocate memory
    ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM,this%m), &
             this%b_high(this%m),&
             this%b_low(this%m),&
             this%c(this%m),&
             this%a(this%m,this%m), &
             STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("InitTimedisc_rkfehlberg", "Unable to allocate memory.")
    END IF

    ! set RHS pointer to the first entry of the coeff field
    this%rhs => Mesh%RemapBounds(this%coeff(:,:,:,:,1))

!NEC$ IEXPAND
    SELECT CASE(this%GetOrder())
    CASE(3)
       !set coefficient scheme of RK-Fehlberg rkf23
       this%b_high = (/ 1.0/6.0, 2.0/3.0, 1.0/6.0 /)
       this%b_low  = (/ 0.0, 1.0, 0.0 /)
       this%c  = (/ 0.0, 0.5, 1.0 /)
       this%a  = TRANSPOSE(RESHAPE((/ &
                 0.0, 0.0, 0.0, &
                 0.5, 0.0, 0.0, &
                -1.0, 2.0, 0.0 /),(/this%m,this%m/)))
    CASE(5)
       !set coefficient scheme of RK-Fehlberg rkf45
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
    END SELECT

    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL this%Error("InitTimedisc_rkfehlberg", &
         "error tolerance levels must be greater than 0")

    IF (this%tol_rel.GT.1.0) THEN
         CALL this%Warning("InitTimedisc_rkfehlberg", &
            "adaptive step size control disabled (tol_rel>1)",0)
    ELSE IF(this%tol_rel.GE.0.01 .AND. this%order .GE. 5) THEN
         CALL this%Warning("InitTimedisc_rkfehlberg", &
             "You chose a relatively high tol_rel (in comparison to order)",0)
    END IF

    CALL GetAttr(config, "ShowButcherTableau", ShowBut, 0)
    IF (ShowBut .EQ. 1) CALL this%ShowButcherTableau()
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
    INTEGER            :: n,i,j,k,l,m
    REAL               :: t
    !------------------------------------------------------------------------!
    t = time
    ! ATTENTION: coeff(:,:,:,:,1) should already be up to date. In the very first
    ! step this is done in FirstStepFosite, any later step must update coeff(:,:,:,1)
    ! at the end of the whole time step update. This is done in AcceptTimestep
    ! or RejectTimestep. Don't forget: this%rhs => this%coeff(:,:,:,1).
    DO m=2,this%m
       ! time step update of cell mean values
!NEC$ IEXPAND
       CALL this%ComputeCVar_rkfehlberg(Mesh,Physics,Fluxes,dt,m,this%coeff,this%cvar,this%ctmp)
       ! compute right-hand-side
       ! coeff_m is k_m/dt from Butcher tableau
       CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,t+this%c(m)*dt,dt,&
                           this%ptmp,this%ctmp,CHECK_NOTHING,this%coeff(:,:,:,:,m))
    END DO

    !reset ctmp
    this%ctmp(:,:,:,:) = this%cvar(:,:,:,:)
!NEC$ NOVECTOR
    DO m=1,this%m
!NEC$ NOVECTOR
      DO l=1,Physics%VNUM
!NEC$ COLLAPSE
       DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ COLLAPSE
        DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
             ! compute two solutions with different numerical orders
             ! y_n+1 = y_n + SUM(b_i*k_i) = y_n + SUM(b_i*dt*coeff_i)
             this%ctmp(i,j,k,l) = this%ctmp(i,j,k,l) &
                                - this%b_low(m)*dt*this%coeff(i,j,k,l,m)
             this%cvar(i,j,k,l) = this%cvar(i,j,k,l) &
                                - this%b_high(m)*dt*this%coeff(i,j,k,l,m)
          END DO
        END DO
      END DO
     END DO
    END DO

    ! at the boundary the this%rhs contains the boundary fluxes
    ! (see subroutine ComputeRHS_modeuler )
    DO m=1,this%m
      DO l=1,Physics%VNUM
        ! western and eastern
        DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
          DO j=Mesh%JMIN,Mesh%JMAX
             Fluxes%bxflux(j,k,1,l) = Fluxes%bxflux(j,k,1,l) &
                               - dt*this%b_high(m)*this%coeff(Mesh%IMIN-Mesh%IP1,j,k,l,m)
             Fluxes%bxflux(j,k,2,l) = Fluxes%bxflux(j,k,2,l) &
                               - dt*this%b_high(m)*this%coeff(Mesh%IMAX+Mesh%IP1,j,k,l,m)
          END DO
        ! southern and northern
!NEC$ IVDEP
          DO i=Mesh%IMIN,Mesh%IMAX
            Fluxes%byflux(k,i,1,l) = Fluxes%byflux(k,i,1,l) &
                              - dt*this%b_high(m)*this%coeff(i,Mesh%JMIN-Mesh%JP1,k,l,m)
            Fluxes%byflux(k,i,2,l) = Fluxes%byflux(k,i,2,l) &
                              - dt*this%b_high(m)*this%coeff(i,Mesh%JMAX+Mesh%JP1,k,l,m)
          END DO
        END DO
        DO j=Mesh%JMIN,Mesh%JMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            Fluxes%bzflux(i,j,1,l) = Fluxes%bzflux(i,j,1,l) &
                              - dt*this%b_high(m)*this%coeff(i,j,Mesh%KMIN-Mesh%KP1,l,m)
            Fluxes%bzflux(i,j,2,l) = Fluxes%bzflux(i,j,2,l) &
                              - dt*this%b_high(m)*this%coeff(i,j,Mesh%KMAX+Mesh%KP1,l,m)
          END DO
        END DO
      END DO
    END DO

    ! compute an error estimate based on the two independent numerical
    ! solutions err and dt
    err = this%ComputeError(Mesh,Physics,this%cvar,this%ctmp)
    dt = this%AdjustTimestep(err,dt)

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
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
             INTENT(IN)          :: cvar
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
             INTENT(OUT)         :: cnew
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM,this%m), &
             INTENT(IN)          :: coeff
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,mm
    !------------------------------------------------------------------------!
    ! compute y + SUM(a_mi*k_i) = y + SUM(a_mi*dt*rhs_i)  : i in [1,m-1] 
    ! see Runge-Kutta method (Butcher tableau)
!     cnew(:,:,:) = cvar(:,:,:) - this%a(m,1)*dt*coeff(:,:,:,1)
! 
!     ! time step update of boundary fluxes
!     DO j=Mesh%JMIN,Mesh%JMAX
!        ! western and eastern boundary
!        Fluxes%bxflux(j,1,:) = Fluxes%bxflux(j,1,:) - this%a(m,1)*dt*coeff(Mesh%IMIN-1,j,:,1)
!        Fluxes%bxflux(j,2,:) = Fluxes%bxflux(j,2,:) - this%a(m,1)*dt*coeff(Mesh%IMAX+1,j,:,1)
!     END DO
! 
!    ! southern and northern boundary fluxes
! !NEC$ IVDEP
!    DO i=Mesh%IMIN,Mesh%IMAX
!       ! time step update of boundary fluxes
!       Fluxes%byflux(i,1,:) = Fluxes%byflux(i,1,:) - this%a(m,1)*dt*coeff(i,Mesh%JMIN-1,:,1)
!       Fluxes%byflux(i,2,:) = Fluxes%byflux(i,2,:) - this%a(m,1)*dt*coeff(i,Mesh%JMAX+1,:,1)
!    END DO

   cnew(:,:,:,:) = cvar(:,:,:,:)
!NEC$ NOVECTOR   
    DO mm=1,m-1
       cnew(:,:,:,:) = cnew(:,:,:,:) - this%a(m,mm)*dt*coeff(:,:,:,:,mm)

!        ! time step update of boundary fluxes
!        DO j=Mesh%JMIN,Mesh%JMAX
!           ! western and eastern boundary
!           Fluxes%bxflux(j,1,k) = Fluxes%bxflux(j,1,k) - this%a(m,mm)*dt*coeff(Mesh%IMIN-1,j,k,mm)
!           Fluxes%bxflux(j,2,k) = Fluxes%bxflux(j,2,k) - this%a(m,mm)*dt*coeff(Mesh%IMAX+1,j,k,mm)
!        END DO
! 
!        ! southern and northern boundary fluxes
! !NEC IVDEP
!        DO i=Mesh%IMIN,Mesh%IMAX
!           ! time step update of boundary fluxes
!           Fluxes%byflux(i,1,k) = Fluxes%byflux(i,1,k) - this%a(m,mm)*dt*coeff(i,Mesh%JMIN-1,k,mm)
!           Fluxes%byflux(i,2,k) = Fluxes%byflux(i,2,k) - this%a(m,mm)*dt*coeff(i,Mesh%JMAX+1,k,mm)
!        END DO
    END DO
  END SUBROUTINE ComputeCVar_rkfehlberg

  SUBROUTINE ShowButcherTableau(this)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Timedisc_rkfehlberg) :: this
    !------------------------------------------------------------------------!
    INTEGER            :: i,j
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
      ELSE
        WRITE (sformat,'(A,I1,A)') '(ES12.3,A,',i-1,'ES12.3)'
      END IF
      WRITE (buffer,fmt=sformat) this%c(i), " | " , this%a(i,1:i-1)
      CALL this%Info(buffer)
    END DO
    CALL this%Info(repeat("-",(this%m+1)*12+4))
    WRITE (sformat,'(A,I1,A)') '(A,',this%m,'ES12.3)'
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
     DEALLOCATE(this%b_high,this%b_low,this%c,this%a)

     IF(ASSOCIATED(this%coeff)) &
      DEALLOCATE(this%coeff)

    CALL this%timedisc_modeuler%Finalize()
  END SUBROUTINE Finalize

END MODULE timedisc_rkfehlberg_mod
