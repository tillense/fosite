!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
  USE fluxes_base_mod
  USE boundary_base_mod
  USE physics_base_mod
  USE sources_base_mod
!  USE gravity_base
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
    PROCEDURE :: ComputeCVar
    PROCEDURE :: SolveODE
    PROCEDURE :: FinalizeTimedisc_modeuler
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
!       FargoAddVelocity, &
!       FargoSubstractVelocity, &
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_modeuler(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)        :: Mesh
    CLASS(physics_base),  INTENT(IN)        :: Physics
    TYPE(Dict_TYP), POINTER                 :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                                 :: method
    !------------------------------------------------------------------------!
    CALL this%InitTimedisc(Mesh,Physics,config,IO,MODIFIED_EULER,ODEsolver_name)
    ! set default order
    CALL GetAttr(config, "order", this%order, 3)

    CALL GetAttr(config, "method", method)
!    CALL this%InitTimedisc(method,ODEsolver_name)


!CDIR IEXPAND
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


  SUBROUTINE SolveODE(this,Mesh,Physics,Fluxes,time,dt,err)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(INOUT) :: Physics
    CLASS(fluxes_base),       INTENT(INOUT) :: Fluxes
    REAL                                    :: time,dt,err
    !------------------------------------------------------------------------!
    INTEGER                                 :: n
    INTEGER                                 :: order
    REAL                                    :: t
    TYPE var_typ
       REAL, DIMENSION(:,:,:,:), POINTER    :: var
    END TYPE var_typ
    TYPE(var_typ)                           :: p(4),c(4)
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: time
    INTENT(INOUT)                           :: dt,err
    !------------------------------------------------------------------------!
    t = time
!CDIR IEXPAND
    order = this%GetOrder()
    ! check if adaptive step size control is enabled
    IF (this%tol_rel.GE.1.0) THEN
       ! no adaptive step size control
!CDIR UNROLL=3
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! time step update of cvar and bfluxes
          CALL this%ComputeCVar(Mesh,Physics,Fluxes,eta(n,order), &
               t,dt,this%cold,this%pvar,this%cvar,this%rhs,this%cvar)
          ! compute right hand side for next time step update
          CALL this%ComputeRHS(Mesh,Physics,Fluxes,t,dt,&
               this%pvar,this%cvar,CHECK_NOTHING,this%rhs)
       END DO
       err = 0.0
       dt = HUGE(dt)
    ELSE
       p(1)%var => Mesh%RemapBounds(this%pvar)
       p(2)%var => Mesh%RemapBounds(this%ptmp)  ! store intermediate result for error control
       p(3)%var => Mesh%RemapBounds(this%pvar)
       p(4)%var => Mesh%RemapBounds(this%pvar)
       c(1)%var => Mesh%RemapBounds(this%cvar)
       c(2)%var => Mesh%RemapBounds(this%ctmp)  ! store intermediate result for error control
       c(3)%var => Mesh%RemapBounds(this%cvar)
       c(4)%var => Mesh%RemapBounds(this%cvar)
!CDIR UNROLL=3
       DO n=1,order
          ! update time variable
          t = time+zeta(n,order)*dt
          ! time step update of cvar and bfluxes
          CALL this%ComputeCVar(Mesh,Physics,Fluxes,eta(n,order), &
               t,dt,this%cold,p(n)%var,c(n)%var,this%rhs,c(n+1)%var)
          ! for 3rd order scheme compute the 2nd order result with the same RHS
          ! and store it in this%ctmp, bfluxes are not required
          IF (n.EQ.2.AND.order.EQ.3) &
!CDIR IEXPAND
             this%ctmp(:,:,:,:) = UpdateTimestep_modeuler(eta(2,2),dt,this%cold(:,:,:,:), &
                                this%ctmp(:,:,:,:),this%rhs(:,:,:,:))
          ! compute right hand side for next time step update
          IF (n.LT.order) &
             CALL this%ComputeRHS(Mesh,Physics,Fluxes,t,dt,p(n+1)%var,c(n+1)%var,&
               CHECK_NOTHING,this%rhs)
       END DO

       err = this%ComputeError(Mesh,Physics,this%cvar,this%ctmp)
       dt = this%AdjustTimestep(err,dt)

    END IF

  END SUBROUTINE SolveODE


  !> \private performs the time step update using the RHS
  !!
  SUBROUTINE ComputeCVar(this,Mesh,Physics,Fluxes,eta,time,dt, &
                                  cold,pvar,cvar,rhs,cnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(INOUT) :: Physics
    CLASS(fluxes_base),       INTENT(INOUT) :: Fluxes
    REAL                                    :: eta,time,dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                   Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM)          &
                                            :: cold,pvar,cvar,cnew,rhs
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k,l
    !------------------------------------------------------------------------!
    INTENT(IN)                              :: eta,time,dt,cold,pvar,cvar,rhs
    INTENT(OUT)                             :: cnew
    !------------------------------------------------------------------------!
!CDIR NOVECTOR
    DO l=1,Physics%VNUM
!CDIR OUTERUNROLL=8
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
          DO i=Mesh%IMIN,Mesh%IMAX
            ! time step update of conservative variables
            cnew(i,j,k,l) = UpdateTimestep_modeuler(eta,dt,cold(i,j,k,l),cvar(i,j,k,l),rhs(i,j,k,l))
          END DO
        END DO
      END DO

      ! western and eastern boundary fluxes
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          ! time step update of boundary fluxes
!CDIR IEXPAND
          Fluxes%bxflux(j,k,1,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,k,1,l), &
               Fluxes%bxflux(j,k,1,l),rhs(Mesh%IMIN-Mesh%Ip1,j,k,l))
!CDIR IEXPAND
          Fluxes%bxflux(j,k,2,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bxfold(j,k,2,l), &
               Fluxes%bxflux(j,k,2,l),rhs(Mesh%IMAX+Mesh%Ip1,j,k,l))
        END DO
      END DO

      ! southern and northern boundary fluxes
      DO i=Mesh%IMIN,Mesh%IMAX
!CDIR NODEP
        DO k=Mesh%KMIN,Mesh%KMAX
          ! time step update of boundary fluxes
!CDIR IEXPAND
          Fluxes%byflux(k,i,1,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(k,i,1,l), &
               Fluxes%byflux(k,i,1,l),rhs(i,Mesh%JMIN-Mesh%Jp1,k,l))
!CDIR IEXPAND
          Fluxes%byflux(k,i,2,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%byfold(k,i,2,l), &
               Fluxes%byflux(k,i,2,l),rhs(i,Mesh%JMAX+Mesh%Jp1,k,l))
        END DO
      END DO

      ! bottom and top boundary fluxes
      DO j=Mesh%JMIN,Mesh%JMAX
        ! time step update of boundary fluxes
        DO i=Mesh%IMIN,Mesh%IMAX
!CDIR IEXPAND
          Fluxes%bzflux(i,j,1,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bzfold(i,j,1,l), &
               Fluxes%bzflux(i,j,1,l),rhs(i,j,Mesh%KMIN-Mesh%Kp1,l))
!CDIR IEXPAND
          Fluxes%bzflux(i,j,2,l) = UpdateTimestep_modeuler(eta,dt,Fluxes%bzfold(i,j,2,l), &
               Fluxes%bzflux(i,j,2,l),rhs(i,j,Mesh%KMAX+Mesh%Kp1,l))
        END DO
      END DO
    END DO
  END SUBROUTINE ComputeCVar


  SUBROUTINE FinalizeTimedisc_modeuler(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_modeuler) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%rhs)
  END SUBROUTINE FinalizeTimedisc_modeuler


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


!  SUBROUTINE FargoAddVelocity(this,Mesh,Physics)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(timedisc_modeuler)   :: this
!    TYPE(mesh_base)       :: Mesh
!    TYPE(physics_base)    :: Physics
!    !------------------------------------------------------------------------!
!    INTEGER              :: i,j
!    !------------------------------------------------------------------------!
!    INTENT(IN)           :: Mesh,Physics
!    INTENT(INOUT)        :: this
!    !------------------------------------------------------------------------!
!    DO j=Mesh%JGMIN,Mesh%JGMAX
!      DO i=Mesh%IGMIN,Mesh%IGMAX
!        IF(Physics%PRESSURE.GT.0) &
!          this%cvar(i,j,Physics%ENERGY) = &
!            this%cvar(i,j,Physics%ENERGY) &
!            + 0.5*this%cvar(i,j,Physics%DENSITY)&
!              *(2.*this%pvar(i,j,Physics%YVELOCITY)+this%w(i))*this%w(i)
!        this%pvar(i,j,Physics%YVELOCITY) = &
!          this%pvar(i,j,Physics%YVELOCITY) &
!          + this%w(i)
!        this%cvar(i,j,Physics%YVELOCITY) = &
!          this%cvar(i,j,Physics%YVELOCITY) &
!          + this%w(i)*this%cvar(i,j,Physics%DENSITY)
!      END DO
!    END DO
!  END SUBROUTINE FargoAddVelocity
!
!
!  SUBROUTINE FargoSubstractVelocity(this,Mesh,Physics)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    TYPE(timedisc_modeuler)   :: this
!    TYPE(mesh_base)       :: Mesh
!    TYPE(physics_base)    :: Physics
!    !------------------------------------------------------------------------!
!    INTEGER              :: i,j
!    !------------------------------------------------------------------------!
!    INTENT(IN)           :: Mesh,Physics
!    INTENT(INOUT)        :: this
!    !------------------------------------------------------------------------!
!    DO j=Mesh%JGMIN,Mesh%JGMAX
!      DO i=Mesh%IGMIN,Mesh%IGMAX
!        this%pvar(i,j,Physics%YVELOCITY) = &
!          this%pvar(i,j,Physics%YVELOCITY) &
!          - this%w(i)
!        this%cvar(i,j,Physics%YVELOCITY) = &
!          this%cvar(i,j,Physics%YVELOCITY) &
!          - this%w(i)*this%cvar(i,j,Physics%DENSITY)
!        IF(Physics%PRESSURE.GT.0) &
!          this%cvar(i,j,Physics%ENERGY) = &
!            this%cvar(i,j,Physics%ENERGY) &
!            - 0.5*this%cvar(i,j,Physics%DENSITY)&
!              *(2.*this%pvar(i,j,Physics%YVELOCITY)+this%w(i))*this%w(i)
!      END DO
!    END DO
!  END SUBROUTINE FargoSubstractVelocity


END MODULE timedisc_modeuler_mod
