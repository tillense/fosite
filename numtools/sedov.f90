!#############################################################################
!#                                                                           #
!# sedov - blast wave solutions                                              #
!# module: sedov.f90                                                         #
!#                                                                           #
!# Copyright (C) 2019                                                        #
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
!> \author Tobias Illenseer
!!
!! \brief module providing sedov class with basic init and solve subroutines
!!
!----------------------------------------------------------------------------!
MODULE sedov_mod
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
! #define DEBUG_OUTPUT 1
#undef DEBUG_OUTPUT
#ifdef f2003
  USE, INTRINSIC :: iso_fortran_env, ONLY : STDIN=>input_unit, &
                                          STDOUT=>output_unit, &
                                          STDERR=>error_unit
#else
#define STDIN  5
#define STDOUT 6
#define STDERR 0
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: PI  = 3.1415926535897932384626433832795028842
  !--------------------------------------------------------------------------!
  TYPE sedov_typ
    LOGICAL :: initialized = .FALSE.
    INTEGER :: ndim  = 3
    REAL    :: gamma = 1.4
    REAL    :: omega = 0.0
    REAL    :: Vmin  = 2./7.   ! = 2./(ndim+2-omega)/gamma
    REAL    :: Vmax  = 1./3.   ! = 4./(ndim+2-omega)/(gamma+1.)
    REAL    :: rho0  = 1.0
    REAL    :: E0    = 1.0
    REAL    :: p0    = 1.0E-5
    REAL    :: alpha = 0.8511
    REAL    :: eps   = 1.0E-12 ! tolerance for integration and root finding
                               ! smaller values yield better results, but
                               ! may run into convergence problems
    REAL, ALLOCATABLE, DIMENSION(:) :: plist
  CONTAINS
    PROCEDURE :: InitParams
    PROCEDURE :: ComputeEnergyIntegral
    PROCEDURE :: PrintConfiguration
    PROCEDURE :: ComputeSolution
    PROCEDURE :: ShockPosition
    PROCEDURE :: Abort
    PROCEDURE :: Destroy
  END TYPE sedov_typ
  !--------------------------------------------------------------------------!
  INTERFACE sedov_typ
    MODULE PROCEDURE CreateSedov
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: sedov_typ

CONTAINS

  !> constructor for sedov class
  FUNCTION CreateSedov(ndim,gamma,omega,rho0,E0,p0) RESULT(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    TYPE(sedov_typ) :: this
    INTEGER, OPTIONAL, INTENT(IN) :: ndim
    REAL, OPTIONAL, INTENT(IN)    :: gamma,omega,rho0,E0,p0
    !-------------------------------------------------------------------!
    ! perform sanity checks and set basic parameters
    IF (PRESENT(ndim)) THEN
      IF (ndim.LT.1.AND.ndim.GT.3) &
        CALL this%Abort("CreateSedov: ndim must be one of {1,2,3}")
      this%ndim = ndim
    END IF
    IF (PRESENT(gamma)) THEN
      IF (gamma.LE.1.0.OR.gamma.GT.5./3.) &
        CALL this%Abort("CreateSedov: gamma not in ( 1 , 5/3 )")
      this%gamma = gamma
    END IF
    IF (PRESENT(omega)) THEN
      IF (omega.LT.0.0.OR.omega.GE.this%ndim) &
        CALL this%Abort("CreateSedov: omega not in [ 0 , ndim )")
      this%omega = omega
      ! default precision does not work for omega > 0
      IF (this%omega.GT.0.0) this%eps = 1.0E-10
    END IF
    IF (PRESENT(rho0)) THEN
      IF (rho0.LE.TINY(rho0).OR.rho0.GT.HUGE(rho0)) &
        CALL this%Abort("CreateSedov: rho0 not in ( 0 , huge )")
      this%rho0 = rho0
    END IF
    IF (PRESENT(E0)) THEN
      IF (E0.LE.TINY(E0).OR.E0.GT.HUGE(E0)) &
        CALL this%Abort("CreateSedov: E0 not in ( 0 , huge )")
      this%E0 = E0
    END IF
    IF (PRESENT(p0)) THEN
      IF (p0.LE.TINY(p0).OR.p0.GT.HUGE(p0)) &
        CALL this%Abort("CreateSedov: p0 not in ( 0 , huge )")
      this%p0 = p0
    END IF
    CALL this%InitParams()
    this%initialized = .TRUE.
  END FUNCTION CreateSedov

  !> sets constants according to Kamm & Timmes
  !! "On Efficient Generation of Numerically Robust Sedov Solutions", Tech report, 2007
  !! see http://cococubed.asu.edu/papers/la-ur-07-2849.pdf
  SUBROUTINE InitParams(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    !-------------------------------------------------------------------!
    INTEGER           :: err
    REAL              :: gam,gam_p1,gam_m1,j,w,w2,w3,J1,J2
#ifdef DEBUG_OUTPUT
    INTEGER           :: i,inum
    REAL              :: V,dV,x1,x2,x3,x4,f,g,h,lambda,dlambda
#endif
    !-------------------------------------------------------------------!
    ! some frequently used constants
    gam = this%gamma
    gam_p1 = gam + 1.0D+00
    gam_m1 = gam - 1.0D+00
    j  = REAL(this%ndim)
    w  = this%omega
    w2 = (2*gam_m1 + j)/gam
    w3 = j*(2-gam)
    this%Vmin = 2./(j+2-w)/gam
    this%Vmax = 4./(j+2-w)/gam_p1

    IF (ALLOCATED(this%plist)) DEALLOCATE(this%plist)
    ALLOCATE(this%plist(19),STAT=err)
    IF (err.NE.0) CALL this%Abort("InitParams: memory allocation failed")

    this%plist(1)  = gam                                   ! = gamma
    this%plist(2)  = 2.0/(j+2.0-w)                         ! = alpha_0
    this%plist(4)  = gam_m1/(gam*(w-w2))                   ! = alpha_2
    this%plist(3)  = (j+2.0-w)*gam/(2.0+j*gam_m1) &        ! = alpha_1
                       *(2*(j*(2.0-gam)-w)/(gam*(j+2.0-w)**2)-this%plist(4))
    this%plist(5)  = (j-w)/(gam*(w2-w))                    ! = alpha_3
    this%plist(6)  = (j+2-w)*(j-w)/(w3-w)*this%plist(3)    ! = alpha_4
    this%plist(7)  = (w*gam_p1-2*j)/(w3-w)                 ! = alpha_5
    this%plist(8)  = 0.25*(j+2-w)*gam_p1                   ! = a = dx1/dV
    this%plist(9)  = gam_p1/gam_m1                         ! = b
    this%plist(10) = 0.5*gam*(j+2.0-w)                     ! = c = 1/Vmin
    this%plist(12) = 0.5*(2+j*gam_m1)                      ! = e
    this%plist(11) = 1.0 &                    ! = d = a/(a-e) = 1/(1-e/a)
                        /(1.0-this%plist(12)/this%plist(8))
    this%plist(13) = this%plist(9)*this%plist(10)    ! = dx2/dV
    this%plist(14) = -this%plist(11)*this%plist(12)  ! = dx3/dV
    this%plist(15) = -this%plist(9)*this%plist(10) & ! = dx4/dV
                        /this%plist(1)
    this%plist(16) = j                               ! = j = ndim : {1,2,3}
    this%plist(17) = w                               ! = omega
    ! these are set later, see ComputeSolution
    this%plist(18) = 1.0
    this%plist(19) = 1.0

    this%Vmin = this%Vmin + 2*EPSILON(this%Vmin) ! avoid lambda(Vmin) -> infty

    CALL this%ComputeEnergyIntegral(J1,J2)
#ifdef DEBUG_OUTPUT
    PRINT '(A)',"# " // REPEAT("=",77)
    PRINT '(A)',"# Parameter:"
    PRINT '(A)',"#       gamma      omega      dim        alpha       J1         J2"
    PRINT '(A,6(F11.4))',"# ",this%gamma,this%omega,REAL(this%ndim),this%alpha,J1,J2
    PRINT '(A)',"# " // REPEAT("-",77)
    PRINT '(A,10(F10.4))',"# p1-p7  ",this%plist(1:7)
    PRINT '(A,10(F10.4))',"# p8-p13 ",this%plist(8:13)
    PRINT '(A,10(F10.4))',"# p14-p19",this%plist(14:19)
    PRINT '(A)',"# " // REPEAT("=",77)
    inum = 20
    dV = (this%Vmax-this%Vmin)/inum
    PRINT '(A)',"#      lambda       V          f          g          h      dJ1/dlam   dJ2/dlam"
    DO i=0,inum
      V = this%Vmin+i*dV
      CALL funcd_lambda(V,this%plist,x1,x2,x3,lambda,dlambda)
      x4 = this%plist(9)*(1.0-this%plist(10)*V/this%plist(1))
      f = lambda*x1
      g = func_g(x1,x2,x3,x4,this%plist)
      h = func_h(x1,x2,x3,x4,this%plist)
      PRINT '(A,7(F11.4))','# ',lambda,V,f,g,h,func_J1(V,this%plist),func_J2(V,this%plist)
    END DO
#endif
  END SUBROUTINE InitParams

  !> compute the energy integral and set parameter alpha
  SUBROUTINE ComputeEnergyIntegral(this,J1,J2)
    USE integration
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    !-------------------------------------------------------------------!
    REAL, OPTIONAL, INTENT(OUT)     :: J1,J2
    !-------------------------------------------------------------------!
    ! compute energy integral
    J1 = integrate(func_J1,this%Vmin,this%Vmax,this%eps,this%plist,1)
    J2 = integrate(func_J2,this%Vmin,this%Vmax,this%eps,this%plist,1)
    IF (this%ndim.EQ.1) THEN
      ! 1D
      this%alpha = 0.5*(J1+1./(this%gamma-1.0)*J2)
    ELSE ! 2D and 3D
      this%alpha = (this%ndim-1)*PI*(J1+2./(this%gamma-1.0)*J2)
    END IF
  END SUBROUTINE ComputeEnergyIntegral

  !> compute the energy integral and set parameter alpha
  SUBROUTINE ComputeSolution(this,t,r,rho,v,p)
    USE roots
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    REAL, INTENT(IN)  :: t,r
    REAL, INTENT(OUT) :: rho,v,p
    !-------------------------------------------------------------------!
    INTEGER :: error
    REAL    :: Rshock,Vshock,Vstar,rho2,v2,p2,x1,x2,x3,x4,lambda,dlambda
    REAL    :: lnx,lnlambda,lnx1,lnx2,lnx3,lnx4
    !-------------------------------------------------------------------!
    ! compute current shock position
    Rshock = this%ShockPosition(t)
    IF (r.GT.0.0.AND.r.LT.Rshock) THEN
      ! estimate for Vstar valid for 0< Vstar/Vmin - 1 << 1
      lnx = -(LOG(this%plist(9)) + (LOG(r/Rshock) + this%plist(2) &
        *LOG(this%plist(8)/this%plist(10)) + this%plist(3) &
        *LOG(this%plist(11)*(1.0-this%plist(12)/this%plist(10))))/this%plist(4))
      Vstar = (1.0+EXP(lnx)) / this%plist(10)
      IF (lnx.LT.LOG(2*EPSILON(lnx))) THEN
        ! Vstar = Vmin*(1.0+EXP(lnx)) with EXP(lnx) << epsilon
        ! => computation of Vstar and other quantities depending on Vstar fails
        ! => use logarithmic quantities and approximate for Vstar ~ Vmin
        Vshock = this%plist(2) * Rshock / t
        rho2   = this%plist(9) * this%rho0
        v2     = 2./(this%gamma+1.0) * Vshock
        p2     = this%rho0 * v2 * Vshock
        lnlambda = LOG(r/Rshock)
        lnx1     = LOG(this%plist(8)*this%Vmin)      ! = LOG(a*Vmin)
        lnx2     = LOG(this%plist(9)) + lnx          ! = LOG(b*x)
        lnx3     = LOG(this%plist(11)*(1.0-this%plist(12)/this%plist(10))) ! = LOG(d*(1-e/c))
        lnx4     = LOG(this%plist(9)*(1.-this%plist(10)*this%Vmin/this%gamma)) ! = LOG(b*(1-c*Vmin/gamma))
        rho = rho2 * EXP(func_lng(lnx1,lnx2,lnx3,lnx4,this%plist))
        v   = v2 * EXP(lnx1+lnlambda)
        p   = p2 * EXP(func_lnh(lnx1,lnx2,lnx3,lnx4,this%plist))
      ELSE
        ! Newton-Raphson to solve the implicit equation
        this%plist(18) = Rshock
        this%plist(19) = r
        ! Newton-Raphson to solve the implicit equation for Vstar
        ! see Kamm & Timmes, 2007 eq. (68)
        CALL GetRoot_Newton(funcd,this%Vmin,this%Vmax,Vstar,error,this%plist)
!         CALL GetRoot(func,this%Vmin,this%Vmax,Vstar,error,this%plist)
        CALL funcd_lambda(Vstar,this%plist,x1,x2,x3,lambda,dlambda)
        x4 = this%plist(9)*(1.0-this%plist(10)*Vstar/this%gamma)
        Vshock = this%plist(2) * Rshock / t
        rho2   = this%plist(9) * this%rho0
        v2     = 2./(this%gamma+1.0) * Vshock
        p2     = this%rho0 * v2 * Vshock
        rho    = rho2 * func_g(x1,x2,x3,x4,this%plist)
        v      = v2 * lambda * x1
        p      = p2 * func_h(x1,x2,x3,x4,this%plist)
      END IF
    ELSE
      ! return unperturbed state
      rho = this%rho0/r**this%omega
      v   = 0.0
      p   = this%p0
    END IF
  END SUBROUTINE ComputeSolution

  FUNCTION ShockPosition(this,time) RESULT(Rshock)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    REAL, INTENT(IN)  :: time
    REAL              :: Rshock
    !-------------------------------------------------------------------!
    IF (.NOT.this%initialized) CALL this%Abort("sedov uninitialzed, aborting ...")
    Rshock = (this%E0/this%rho0*time*time/this%alpha)**(0.5*this%plist(2))
  END FUNCTION ShockPosition

  !> prints configuration on the screen (stdout)
  SUBROUTINE PrintConfiguration(this)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    !-------------------------------------------------------------------!
    WRITE (STDOUT,"(A)") "# ==========================================="
    WRITE (STDOUT,"(A)") "#  << Solution of Sedov explosion problem >> "
    WRITE (STDOUT,"(A)") "# ==========================================="
    WRITE (STDOUT,"(A,I3)")     "# Geometry / Dimensions       ", this%ndim
    WRITE (STDOUT,"(A,F10.6)") "# Ratio of specific heat cap. ", this%gamma
    WRITE (STDOUT,"(A,ES10.2)") "# Initial energy input        ", this%E0
    WRITE (STDOUT,"(A,ES10.2)") "# Ambient density             ", this%rho0
    WRITE (STDOUT,"(A,ES10.2)") "# Peak density                ", this%rho0*this%plist(9)
    WRITE (STDOUT,"(A,F10.6)") "# Density power law exponent  ", this%omega
    WRITE (STDOUT,"(A,F10.6)") "# Shock position parameter    ", this%alpha
  END SUBROUTINE PrintConfiguration

  SUBROUTINE Abort(this,msg)
    IMPLICIT NONE
    !-------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    !-------------------------------------------------------------------!
    CHARACTER(LEN=*) :: msg
    !-------------------------------------------------------------------!
    WRITE (STDERR,'(A)') "ERROR in " // TRIM(msg)
    STOP
  END SUBROUTINE Abort

  SUBROUTINE Destroy(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sedov_typ), INTENT(INOUT) :: this
    !-------------------------------------------------------------------!
    DEALLOCATE(this%plist)
    this%initialized = .FALSE.
  END SUBROUTINE Destroy

  PURE SUBROUTINE funcd(y,fy,dfy,plist)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: y
    REAL, INTENT(OUT) :: fy,dfy
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    !------------------------------------------------------------------------!
    REAL :: x1,x2,x3,lambda,dlambda
    !------------------------------------------------------------------------!
    CALL funcd_lambda(y,plist,x1,x2,x3,lambda,dlambda)
    fy  = plist(18)*lambda - plist(19)
    dfy = plist(18)*dlambda
  END SUBROUTINE funcd

  PURE SUBROUTINE func(y,fy,plist)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: y
    REAL, INTENT(OUT) :: fy
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    !------------------------------------------------------------------------!
    REAL :: x1,x2,x3,lambda,dlambda
    !------------------------------------------------------------------------!
    CALL funcd_lambda(y,plist,x1,x2,x3,lambda,dlambda)
    fy  = plist(18)*lambda - plist(19)
  END SUBROUTINE func

  FUNCTION func_J1(x,plist) RESULT(fx)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x
    REAL, INTENT(INOUT), DIMENSION(:), OPTIONAL :: plist
    REAL :: fx
    !------------------------------------------------------------------------!
    REAL :: x1,x2,x3,x4,lambda,dlambda,g
    !------------------------------------------------------------------------!
    CALL funcd_lambda(x,plist,x1,x2,x3,lambda,dlambda)
    x4 = plist(9)*(1.0-plist(10)*x/plist(1))
    g  = func_g(x1,x2,x3,x4,plist)
    fx = plist(9) * lambda**(plist(16)+1.0) * g * x**2 * dlambda
  END FUNCTION func_J1

  FUNCTION func_J2(x,plist) RESULT(fx)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x
    REAL, INTENT(INOUT), DIMENSION(:), OPTIONAL :: plist
    REAL :: fx
    !------------------------------------------------------------------------!
    REAL :: x1,x2,x3,x4,lambda,dlambda,h
    !------------------------------------------------------------------------!
    CALL funcd_lambda(x,plist,x1,x2,x3,lambda,dlambda)
    x4 = plist(9)*(1.0-plist(10)*x/plist(1))
    h  = func_h(x1,x2,x3,x4,plist)
    !!! ATTENTION: There is an error in eq. 56 of Kamm & Timmes (2007);
    !!!            it must be lambda^(gamma-1) and NOT ..^(gamma+1)
    fx = 0.5*(plist(1)+1.0)/plist(8)**2 * lambda**(plist(16)-1.0) * h * dlambda
  END FUNCTION func_J2

  PURE SUBROUTINE funcd_lambda(V,plist,x1,x2,x3,lambda,dlambda)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: V
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL, INTENT(OUT) :: x1,x2,x3,lambda,dlambda
    !------------------------------------------------------------------------!
    x1 = plist(8)*V
    x2 = plist(9)*MAX(0.0,plist(10)*V-1.0)  ! should be >=0, since V*p10 >= 1
    x3 = plist(11)*(1.0-plist(12)*V)
    lambda  = x1**(-plist(2)) * x2**(-plist(4)) * x3**(-plist(3))
!    dlambda = -lambda*(plist(2)/x1*plist(8) + plist(4)/x2*plist(13) + plist(3)/x3*plist(14))
    dlambda = -plist(2)*plist(8) * x1**(-plist(2)-1.0) * x2**(-plist(4)) * x3**(-plist(3)) &
              -plist(4)*plist(13) * x1**(-plist(2)) * x2**(-plist(4)-1.0) * x3**(-plist(3)) &
              -plist(3)*plist(14) * x1**(-plist(2)) * x2**(-plist(4)) * x3**(-plist(3)-1.0)
  END SUBROUTINE funcd_lambda

  PURE FUNCTION func_g(x1,x2,x3,x4,plist) RESULT(g)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2,x3,x4
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: g
    !------------------------------------------------------------------------!
    g = x1**(plist(2)*plist(17)) * x2**(plist(5)+plist(4)*plist(17)) &
        *x3**(plist(6)+plist(3)*plist(17)) * x4**plist(7)
  END FUNCTION func_g

  PURE FUNCTION func_lng(lnx1,lnx2,lnx3,lnx4,plist) RESULT(lng)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: lnx1,lnx2,lnx3,lnx4
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: lng
    !------------------------------------------------------------------------!
    lng = (plist(2)*plist(17))*lnx1 + (plist(5)+plist(4)*plist(17))*lnx2 &
        + (plist(6)+plist(3)*plist(17))*lnx3 + plist(7)*lnx4
  END FUNCTION func_lng

  PURE FUNCTION func_h(x1,x2,x3,x4,plist) RESULT(h)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2,x3,x4
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: h
    !------------------------------------------------------------------------!
    h = x1**(plist(2)*plist(16)) * x3**(plist(6)+plist(3)*(plist(17)-2.0)) &
        *x4**(1.0+plist(7))
  END FUNCTION func_h

  PURE FUNCTION func_lnh(lnx1,lnx2,lnx3,lnx4,plist) RESULT(lnh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: lnx1,lnx2,lnx3,lnx4
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    REAL :: lnh
    !------------------------------------------------------------------------!
    lnh = (plist(2)*plist(16))*lnx1 + (plist(6)+plist(3)*(plist(17)-2.0))*lnx3 &
        + (1.0+plist(7))*lnx4
  END FUNCTION func_lnh

END MODULE sedov_mod


  
