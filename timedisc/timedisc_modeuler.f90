!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f03                                             #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Jannes Klee
!!
!! \brief subroutines for modified Euler i.e. Runge-Kutta methods
!!
!! \cite shu1988
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_modeuler_mod
  USE timedisc_base_mod
  USE mesh_base_mod
  USE marray_compound_mod
  USE marray_base_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE physics_base_mod
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
  TYPE, EXTENDS (timedisc_base) :: timedisc_modeuler
  CONTAINS
    PROCEDURE :: InitTimedisc_modeuler
    PROCEDURE :: ComputeCVar_modeuler
    PROCEDURE :: SolveODE
    PROCEDURE :: Finalize
  END TYPE timedisc_modeuler
  !--------------------------------------------------------------------------!
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "modified Euler"
  REAL, PARAMETER :: eta(3,3) = &
       RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.75, 1.0/3.0 /),(/3,3/))
  REAL, PARAMETER :: zeta(3,3) = &
       RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5 /),(/3,3/))
 !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       timedisc_modeuler, &
       ! constants
       DTCAUSE_CFL,DTCAUSE_ERRADJ,DTCAUSE_SMALLERR, &
       CHECK_ALL, CHECK_NOTHING, CHECK_CSOUND, CHECK_PMIN, CHECK_RHOMIN, &
       CHECK_INVALID, CHECK_TMIN
       ! methods
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_modeuler(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(INOUT) :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    TYPE(Dict_TYP),           POINTER       :: config, IO
    !------------------------------------------------------------------------!
!     INTEGER                                 :: method
    !------------------------------------------------------------------------!
    ! set default order
    CALL GetAttr(config, "order", this%order, 3)
!    CALL GetAttr(config, "method", method)

    CALL this%InitTimedisc(Mesh,Physics,config,IO,MODIFIED_EULER,ODEsolver_name)

    SELECT CASE(this%GetOrder())
    CASE(1)
       ! set relative error tolarance to value > 1.0
       ! to disable adaptive step size control
       ! (not available for 1st order scheme)
       this%tol_rel    = 10.0
       this%tol_abs(:) = 1.0
      CALL this%Warning("InitTimedisc_modeuler", &
            "adaptive step size control not supported in 1st order scheme",0)
    CASE(2,3)
       IF (this%tol_rel.GT.1.0) &
            CALL this%Warning("InitTimedisc_modeuler", &
            "adaptive step size control disabled (tol_rel>1)",0)
    CASE DEFAULT
       CALL this%Error("InitTimedisc_modeuler","time order must be one of 1,2,3")
    END SELECT
    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL this%Error("InitTimedisc_modeuler", &
         "error tolerance levels must be greater than 0")

  END SUBROUTINE InitTimedisc_modeuler


  SUBROUTINE SolveODE(this,Mesh,Physics,Sources,Fluxes,time,dt,err)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(INOUT) :: Physics
    CLASS(sources_base),      POINTER       :: Sources
    CLASS(fluxes_base),       INTENT(INOUT) :: Fluxes
    REAL                                    :: time,dt,err
    !------------------------------------------------------------------------!
    INTEGER                                 :: n
    INTEGER                                 :: order
    REAL                                    :: t
    TYPE var_typ
      CLASS(marray_compound), POINTER :: var
    END TYPE var_typ
    TYPE(var_typ)                           :: p(4),c(4)
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: time
    INTENT(INOUT)                           :: dt,err
    !------------------------------------------------------------------------!
    t = time
    order = this%GetOrder()
    ! check if adaptive step size control is enabled
    IF (this%tol_rel.GE.1.0) THEN
       ! no adaptive step size control
!NEC$ UNROLL(3)
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! time step update of cvar and bfluxes
          CALL this%ComputeCVar_modeuler(Mesh,Physics,Fluxes,eta(n,order), &
               t,dt,this%cold,this%pvar,this%cvar,this%rhs,this%cvar)
          ! compute right hand side for next time step update
          CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,t,dt,&
               this%pvar,this%cvar,CHECK_NOTHING,this%rhs)
       END DO
       err = 0.0
       dt = HUGE(dt)
    ELSE
      ! with adaptive step size control (2nd / 3rd order scheme)
       p(1)%var => this%pvar
       p(2)%var => this%ptmp  ! store intermediate result for error control
       p(3)%var => this%pvar
       p(4)%var => this%pvar
       c(1)%var => this%cvar
       c(2)%var => this%ctmp  ! store intermediate result for error control
       c(3)%var => this%cvar
       c(4)%var => this%cvar
!NEC$ UNROLL(3)
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! time step update of cvar and bfluxes
          CALL this%ComputeCVar_modeuler(Mesh,Physics,Fluxes,eta(n,order), &
               t,dt,this%cold,p(n)%var,c(n)%var,this%rhs,c(n+1)%var)
          ! for 3rd order scheme compute the 2nd order result with the same RHS
          ! and store it in this%ctmp, bfluxes are not required
          IF (n.EQ.2.AND.order.EQ.3) &
             this%ctmp%data4d(:,:,:,:) = UpdateTimestep_modeuler(eta(2,2),dt,this%cold%data4d(:,:,:,:), &
                                this%ctmp%data4d(:,:,:,:),this%rhs(:,:,:,:))
          ! compute right hand side for next time step update
          IF (n.LT.order) &
             CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,t,dt,p(n+1)%var,c(n+1)%var,&
               CHECK_NOTHING,this%rhs)
       END DO

       err = this%ComputeError(Mesh,Physics,this%cvar,this%ctmp)
       dt = this%AdjustTimestep(err,dt)

    END IF

  END SUBROUTINE SolveODE


  !> \private performs the time step update using the RHS
  !!
  SUBROUTINE ComputeCVar_modeuler(this,Mesh,Physics,Fluxes,eta,time,dt, &
                                  cold,pvar,cvar,rhs,cnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(INOUT) :: Physics
    CLASS(fluxes_base),       INTENT(INOUT) :: Fluxes
    CLASS(marray_compound),   INTENT(INOUT) :: cold,pvar,cvar,cnew
    REAL                                    :: eta,time,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                   Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM)          &
                                            :: rhs
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k,l
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: eta,time,dt,rhs
    !------------------------------------------------------------------------!
!NEC$ NOVECTOR
    DO l=1,Physics%VNUM
!NEC$ OUTERLOOP_UNROLL(8)
      DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN,Mesh%IMAX
            ! time step update of conservative variables
            cnew%data4d(i,j,k,l) = UpdateTimestep_modeuler(eta,dt, &
               cold%data4d(i,j,k,l),cvar%data4d(i,j,k,l),rhs(i,j,k,l))
          END DO
        END DO
      END DO

      ! western and eastern boundary fluxes
!NEC$ IVDEP
      DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j=Mesh%JMIN,Mesh%JMAX
          ! time step update of boundary fluxes
          Fluxes%bxflux(j,k,1,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,k,1,l), &
               Fluxes%bxflux(j,k,1,l),rhs(Mesh%IMIN-Mesh%Ip1,j,k,l))
          Fluxes%bxflux(j,k,2,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,k,2,l), &
               Fluxes%bxflux(j,k,2,l),rhs(Mesh%IMAX+Mesh%Ip1,j,k,l))
        END DO
      END DO

      ! southern and northern boundary fluxes
!NEC$ IVDEP
      DO i=Mesh%IMIN,Mesh%IMAX
!NEC$ IVDEP
        DO k=Mesh%KMIN,Mesh%KMAX
          ! time step update of boundary fluxes
          Fluxes%byflux(k,i,1,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(k,i,1,l), &
               Fluxes%byflux(k,i,1,l),rhs(i,Mesh%JMIN-Mesh%Jp1,k,l))
          Fluxes%byflux(k,i,2,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(k,i,2,l), &
               Fluxes%byflux(k,i,2,l),rhs(i,Mesh%JMAX+Mesh%Jp1,k,l))
        END DO
      END DO

      ! bottom and top boundary fluxes
!NEC$ IVDEP
      DO j=Mesh%JMIN,Mesh%JMAX
        ! time step update of boundary fluxes
!NEC$ IVDEP
        DO i=Mesh%IMIN,Mesh%IMAX
          Fluxes%bzflux(i,j,1,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bzfold(i,j,1,l), &
               Fluxes%bzflux(i,j,1,l),rhs(i,j,Mesh%KMIN-Mesh%Kp1,l))
          Fluxes%bzflux(i,j,2,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bzfold(i,j,2,l), &
               Fluxes%bzflux(i,j,2,l),rhs(i,j,Mesh%KMAX+Mesh%Kp1,l))
        END DO
      END DO
    END DO
  END SUBROUTINE ComputeCVar_modeuler


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize


  ELEMENTAL FUNCTION UpdateTimestep_modeuler(eta_n,dt,y0,yn,rhs) RESULT(y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: eta_n,dt,y0,yn,rhs
    REAL             :: y
    !------------------------------------------------------------------------!
    ! ATTENTION:
    ! The time step update is computed according to:
    !    y = eta_n * y0 + (1.0 - eta_n) * (yn - dt * rhs)
    ! but to minimize the truncation error it is essential to sort the terms
    ! in this way:
    y = yn-dt*rhs+eta_n*(y0-yn+dt*rhs)
  END FUNCTION UpdateTimestep_modeuler

END MODULE timedisc_modeuler_mod
