!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: roots.f90                                                         #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
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
!! \brief root finding subroutines
!!
!! 1. Newton's method combined with bisection, for ill-posed problems
!! 2. Regula falsi
!! 3. Bisection
!! 4. Pegasus
!! 5. King
!! 6. Anderson Bjoerk
!! 7. Ridder
!! 8. Brent Dekker (default)
!!
!! The algorithms can be found in standard textbooks on numerical analysis.
!! See e. g. \cite press2007 and \cite engeln2011 and references therein.
!!
!! \ingroup numtools
!----------------------------------------------------------------------------!
MODULE roots
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !> basic type for root finding functions
  TYPE Roots_TYP
    REAL    :: eps,tol,root,dx,dxold,df,d,e
    REAL    :: x(3),f(3)
    INTEGER :: iter,max_iterations,error
    LOGICAL :: iterate,do_modified_step
  END TYPE Roots_TYP
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetRoot
     MODULE PROCEDURE  GetRoot_BrentDekker
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  ! some constants
! #define DEBUG_OUTPUT 1
#undef DEBUG_OUTPUT
  INTEGER, PARAMETER :: DEFAULT_MAX_ITERATIONS = 1000
  REAL, PARAMETER    :: DEFAULT_ACCURACY       = 4*EPSILON(DEFAULT_ACCURACY)
  CHARACTER(LEN=64), PARAMETER, DIMENSION(0:4) ::  ERROR_MESSAGE = (/ &
     "unknown error                                                   ", &
     "root not bracketed                                              ", &
     "iteration exceeds maximum                                       ", &
     "requested accuracy smaller than machine precission              ", &
     "upper limit for iterations should be larger than 1              " /)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! constants
       DEFAULT_ACCURACY, &
       DEFAULT_MAX_ITERATIONS, &
       ! types
       Roots_TYP, &
       ! methods
       InitRoots, &
       GetErrorMessage, &
       GetRoot, &
       GetRoot_Newton, &
       GetRoot_RegulaFalsi, &
       GetRoot_Bisection, &
       GetRoot_Pegasus, &
       GetRoot_King, &
       GetRoot_AndersonBjoerk, &
       GetRoot_Ridder, &
       GetRoot_BrentDekker
  !--------------------------------------------------------------------------!

CONTAINS

  PURE SUBROUTINE InitRoots(this,x1,x2,f1,f2,dxacc,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    REAL, INTENT(IN)                :: x1,x2
    REAL, INTENT(IN)                :: f1,f2
    REAL,OPTIONAL, INTENT(IN)       :: dxacc
    INTEGER,OPTIONAL, INTENT(IN)    :: maxiter
    !------------------------------------------------------------------------!
    ! check if root is bracketed
    IF (f1*f2.GT.0.0) THEN
        ! f1 and f2 should have opposite signs
        this%error = 1
    ELSE
        IF (x1.LT.x2) THEN
            this%x(1) = x1
            this%x(2) = x2
            this%f(1) = f1
            this%f(2) = f2
        ELSE
            this%x(1) = x2
            this%x(2) = x1
            this%f(1) = f2
            this%f(2) = f1
        END IF

        this%dx = ABS(this%x(2)-this%x(1))
        this%dxold = this%dx

        IF(PRESENT(dxacc)) THEN
            IF (dxacc.GT.EPSILON(dxacc)) THEN
                this%eps = dxacc
            ELSE
                this%error = 3
            END IF
        ELSE
            this%eps = DEFAULT_ACCURACY
        END IF

        IF(PRESENT(maxiter)) THEN
            IF (maxiter.GT.1) THEN
                this%max_iterations = maxiter
            ELSE
                this%error = 4
            END IF
        ELSE
            this%max_iterations = DEFAULT_MAX_ITERATIONS
        END IF

        ! only used in Brent-Dekker method
        this%d = this%x(2)-this%x(1)
        this%e = this%d

        this%error = 0
        this%iter  = 0
        this%iterate = .TRUE.
        this%do_modified_step = .TRUE.
    END IF

  END SUBROUTINE InitRoots


  FUNCTION GetErrorMessage(error) RESULT(msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: error
    CHARACTER(LEN=64)   :: msg
    !------------------------------------------------------------------------!
    SELECT CASE (error)
    CASE(1,2)
        msg = ERROR_MESSAGE(error)
    CASE DEFAULT
        msg = ERROR_MESSAGE(0)
    END SELECT
  END FUNCTION GetErrorMessage


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_generic(this,Stepper,func,x1,x2,dxacc,maxiter,plist,xm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP)   :: this
    REAL              :: x1,x2
    REAL, OPTIONAL    :: dxacc, plist(:), xm
    INTEGER, OPTIONAL :: maxiter
    !------------------------------------------------------------------------!
    INTENT(IN)        :: x1,x2,dxacc,maxiter,xm,plist
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE Stepper(root)
         IMPORT roots_typ
         IMPLICIT NONE
         TYPE(Roots_TYP), INTENT(INOUT) :: root
       END SUBROUTINE Stepper
    END INTERFACE
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    REAL            :: f1,f2
    !------------------------------------------------------------------------!
    ! compute function values at the boundaries
    CALL func(x1,f1,plist)
    CALL func(x2,f2,plist)

    ! initialize root finding algorithm
    CALL InitRoots(this,x1,x2,f1,f2)

    IF (this%error.EQ.0) THEN
        ! initial guess for the root
        IF (PRESENT(xm)) THEN
            this%x(3) = xm
        ELSE
            this%x(3) = ArithmeticMean(this%x(1),this%x(2))
        END IF
#ifdef DEBUG_OUTPUT
        PRINT '(A3,A18,A10,A18,A18,2(A13,A10))', &
             "#","x3","f3","dx","tol","x1","f1","x2","f2"
#endif
        ! main loop
        DO
            ! compute function value at x3
            CALL func(this%x(3),this%f(3),plist)
            ! print debug information if requested
#ifdef DEBUG_OUTPUT
            PRINT '(I3,ES18.10,ES10.2,2(ES18.10),2(ES13.5,ES10.2))', &
                this%iter,this%x(3),this%f(3),this%dx,this%tol,this%x(1),this%f(1),this%x(2),this%f(2)
#endif
            ! check convergence
            CALL TestConvergence(this)
            IF (.NOT.this%iterate) EXIT

            ! determine new estimate for x3
            CALL Stepper(this)
        END DO
    END IF
  END SUBROUTINE GetRoot_generic

#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_Newton(funcd,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE funcd(x,fx,dfx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx,dfx
       END SUBROUTINE funcd
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    REAL    :: f1,f2,df1,df2
    !------------------------------------------------------------------------!
    ! compute function values and derivatives at the boundaries
    CALL funcd(x1,f1,df1,plist)
    CALL funcd(x2,f2,df2,plist)

    ! initialize root finding algorithm
    CALL InitRoots(this,x1,x2,f1,f2)

    IF (this%error.EQ.0) THEN
        ! initial guess for the root
        IF (PRESENT(xm)) THEN
            this%x(3) = xm
        ELSE
            this%x(3) = ArithmeticMean(this%x(1),this%x(2))
        END IF
#ifdef DEBUG_OUTPUT
        PRINT '(A3,A18,A10,A18,A18,2(A13,A10))', &
             "#","x3","f3","dx","tol","x1","f1","x2","f2"
#endif
        ! main loop
        DO
            ! compute function value and derivative at x3
            CALL funcd(this%x(3),this%f(3),this%df,plist)
            ! print debug information if requested
#ifdef DEBUG_OUTPUT
            PRINT '(I3,ES18.10,ES10.2,2(ES18.10),2(ES13.5,ES10.2))', &
                this%iter,this%x(3),this%f(3),this%dx,this%tol,this%x(1),this%f(1),this%x(2),this%f(2)
#endif
            ! check convergence
            CALL TestConvergence(this)
            IF (.NOT.this%iterate) EXIT

            ! determine new interval for enclosure and estimate x3
            CALL Step_Newton(this)
       END DO
    END IF

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_Newton


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_RegulaFalsi(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_RegulaFalsi,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_RegulaFalsi


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_Pegasus(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_Pegasus,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_Pegasus


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_King(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_King,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_King


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_AndersonBjoerk(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_AndersonBjoerk,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_AndersonBjoerk


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_Ridder(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_Ridder,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_Ridder


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_BrentDekker(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_BrentDekker,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_BrentDekker


#ifndef DEBUG_OUTPUT
  PURE &
#endif
  SUBROUTINE GetRoot_Bisection(func,x1,x2,root,error,plist,xm,iterations)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2
    REAL, INTENT(OUT) :: root
    INTEGER, INTENT(OUT) :: error
    REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: plist
    REAL, INTENT(IN), OPTIONAL :: xm
    INTEGER, INTENT(OUT), OPTIONAL :: iterations
    !------------------------------------------------------------------------!
    INTERFACE
       PURE SUBROUTINE func(x,fx,plist)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
         REAL, INTENT(OUT) :: fx
       END SUBROUTINE func
    END INTERFACE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP) :: this
    !------------------------------------------------------------------------!
    CALL GetRoot_generic(this,Step_Bisection,func,x1,x2, &
        DEFAULT_ACCURACY,DEFAULT_MAX_ITERATIONS,plist,xm)

    root  = this%root
    error = this%error
    IF (PRESENT(iterations)) iterations = this%iter

  END SUBROUTINE GetRoot_Bisection


  PURE SUBROUTINE UpdateBounds_Bisection(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    ! check sign
    IF (this%f(1)*this%f(3).GT.0.0) THEN
        ! sign(f1) == sign(f3) => x3 < root < x2
        this%x(1)=this%x(3)
        this%f(1)=this%f(3)
    ELSE
        ! sign(f1) != sign(f3) => x1 < root < x3
        this%x(2)=this%x(3)
        this%f(2)=this%f(3)
    END IF
  END SUBROUTINE UpdateBounds_Bisection


  PURE SUBROUTINE Step_Bisection(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    CALL UpdateBounds_Bisection(this)
    CALL Step_ArithmeticMean(this)
  END SUBROUTINE Step_Bisection


  PURE SUBROUTINE Step_Newton(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    ! determine new brackets for the root
!CDIR IEXPAND
    CALL UpdateBounds_Bisection(this)
    this%dxold = this%dx
    ! check if we are out of bounds
    ! or if the convergence is to slow
    IF ( ( (this%df*(this%x(3)-this%x(2))-this%f(3).GE.0.0) &
        .OR. (this%df*(this%x(3)-this%x(1))-this%f(3).LE.0.0) ) &
        .OR. (ABS(2*this%f(3)).GT.ABS(this%dxold*this%df)) ) THEN
        ! compute new estimate for the route using arithmetic mean value
!CDIR IEXPAND
        CALL Step_ArithmeticMean(this)
    ELSE
        ! compute new estimate for the root using Newton iteration
        this%dx   = this%f(3) / this%df
        this%x(3) = this%x(3) - this%dx
        this%dx   = ABS(this%dx)
    END IF
  END SUBROUTINE Step_Newton


  PURE SUBROUTINE Step_RegulaFalsi(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    ! calculate new interval boundaries
!CDIR IEXPAND
    CALL UpdateBounds_Bisection(this)
    ! determine new approximation for the root
!CDIR IEXPAND
    CALL Step_SecantMethod(this)
  END SUBROUTINE Step_RegulaFalsi


  PURE SUBROUTINE Step_Pegasus(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    REAL             :: g
    !------------------------------------------------------------------------!
    ! calculate new interval boundaries
    IF (this%f(2)*this%f(3).LT.0.0) THEN
       this%x(1) = this%x(2)
       this%x(2) = this%x(3)
       this%f(1) = this%f(2)
       this%f(2) = this%f(3)
    ELSE
       g = this%f(2)/(this%f(3)+this%f(2))
       this%x(2) = this%x(3)
       this%f(1) = g*this%f(1)
       this%f(2) = this%f(3)
    END IF
    ! determine new approximation for the root
!CDIR IEXPAND
    CALL Step_SecantMethod(this)
  END SUBROUTINE Step_Pegasus


  PURE SUBROUTINE Step_King(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    REAL :: g
    !------------------------------------------------------------------------!
    IF (this%do_modified_step) THEN
        IF (this%f(2)*this%f(3).LT.0.0) THEN
            ! interchange (x1,f1) and (x2,f2)
            g       = this%x(1)
            this%x(1) = this%x(2)
            this%x(2) = g
            g       = this%f(1)
            this%f(1) = this%f(2)
            this%f(2) = g
        END IF
        ! Pegasus step
        IF (this%f(2)*this%f(3).GT.0.0) THEN
            g = this%f(2)/(this%f(3)+this%f(2))
            this%x(2) = this%x(3)
            this%f(1) = g*this%f(1)
            this%f(2) = this%f(3)
        END IF
        this%do_modified_step = .FALSE.
    ELSE
        IF (this%f(2)*this%f(3).LT.0.0) THEN
            this%x(1) = this%x(2)
            this%f(1) = this%f(2)
            this%x(2) = this%x(3)
            this%f(2) = this%f(3)
            this%do_modified_step = .TRUE.
        ELSE
            ! Pegasus step
            IF (this%f(2)*this%f(3).GT.0.0) THEN
                g = this%f(2)/(this%f(3)+this%f(2))
                this%x(2) = this%x(3)
                this%f(1) = g*this%f(1)
                this%f(2) = this%f(3)
            END IF
            this%do_modified_step = .FALSE.
        END IF
    END IF
    ! determine new approximation for the root
!CDIR IEXPAND
    CALL Step_SecantMethod(this)
  END SUBROUTINE Step_King


  PURE SUBROUTINE Step_AndersonBjoerk(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    REAL             :: g
    !------------------------------------------------------------------------!
    ! calculate new interval boundaries
    IF (this%f(2)*this%f(3).LT.0.0) THEN
       this%x(1) = this%x(2)
       this%x(2) = this%x(3)
       this%f(1) = this%f(2)
       this%f(2) = this%f(3)
    ELSE
       g = 1.-this%f(3)/this%f(2)
       IF(g.LE.0.) g = 0.5
       this%x(2) = this%x(3)
       this%f(1) = g*this%f(1)
       this%f(2) = this%f(3)
    END IF
    ! determine new approximation for the root
!CDIR IEXPAND
    CALL Step_SecantMethod(this)
  END SUBROUTINE Step_AndersonBjoerk
  
  
  PURE SUBROUTINE Step_Ridder(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    REAL :: denom,dx
    !------------------------------------------------------------------------!
    IF (this%do_modified_step) THEN
        ! compute arithmetic mean of the new boundaries
        CALL Step_Bisection(this)
        ! next update uses Ridders 3 point formula
        this%do_modified_step = .FALSE.
    ELSE
        ! estimate the new step from the arithmetic mean value
        denom = SQRT(this%f(3)*this%f(3) - this%f(1)*this%f(2)) + TINY(denom)
        dx = (this%x(3)-this%x(1)) * SIGN(1.0,this%f(1)-this%f(2))*this%f(3) / denom
        IF (dx.GT.0.0) THEN
            this%x(1) = this%x(3)
            this%f(1) = this%f(3)
        ELSE
            this%x(2) = this%x(3)
            this%f(2) = this%f(3)
        END IF
        this%dx = 0.5*this%dx
        this%x(3) = this%x(3) + dx
        ! next step is a bisection step
        this%do_modified_step = .TRUE.
    END IF
  END SUBROUTINE Step_Ridder


  PURE SUBROUTINE Step_BrentDekker(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    REAL :: m,p,q
    !------------------------------------------------------------------------!
    ! bracket the root
    IF (this%f(3)*this%f(2).GT.0.0) THEN
        ! x3, x2 on the same side of the root
        this%x(2) = this%x(1)
        this%f(2) = this%f(1)
        this%d  = this%x(1)-this%x(3)
        this%e  = this%d
    END IF
    ! root is now bracketed between x3 and x2

    ! x3 should be closer to the root than x2, i.e. |f3|<|f2|
    IF (ABS(this%f(2)).LT.ABS(this%f(3))) THEN
        ! interchange (x3,f3) and (x2,f2) using (x1,f1) as temporary space
        this%x(1) = this%x(3)
        this%x(3) = this%x(2)
        this%x(2) = this%x(1)
        this%f(1) = this%f(3)
        this%f(3) = this%f(2)
        this%f(2) = this%f(1)
    END IF

    m = 0.5*(this%x(2)-this%x(3))

    IF ((ABS(this%e).LT.this%tol) .OR. ABS(this%f(1)).LE.ABS(this%f(3))) THEN
        ! convergence too slow -> do bisection step
        this%d = m
        this%e = m
    ELSE
        ! try interpolation
        CALL SaveInverseQuadraticInterpolation(this%x(1),this%x(3),this%x(2), &
               this%f(1),this%f(3),this%f(2),p,q)

        IF (2*p.LT. MIN(3*m*q-ABS(this%tol*q), ABS(this%e*q))) THEN
            ! interpolation worked store old correction d in e
            ! and compute new correction d
            this%e = this%d
            this%d = p / q
        ELSE
            ! interpolation failed, do bisection
            this%e = m
            this%d = m
        END IF
    END IF

    this%x(1) = this%x(3)
    this%f(1) = this%f(3)

    IF (ABS(this%d).GT.this%tol) THEN
        this%x(3) = this%x(3) + this%d
    ELSE
        this%x(3) = this%x(3) + SIGN(this%tol,m)
    END IF
    
    this%dx = ABS(this%x(2)-this%x(3))
  END SUBROUTINE Step_BrentDekker


  PURE SUBROUTINE Step_SecantMethod(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    INTEGER :: ONE,TWO
    !------------------------------------------------------------------------!
    IF (ABS(this%f(1)).GT.ABS(this%f(2))) THEN
      ! x2 is probably the better approximation for the root
      ! set indices as usual
      ONE = 1
      TWO = 2
    ELSE
      ! x1 is probably the better approximation for the root
      ! interchange indices
      ONE = 2
      TWO = 1
    END IF
!CDIR IEXPAND
    this%dx = LinearInterpolation(this%x(ONE),this%x(TWO),this%f(ONE),this%f(TWO))
    ! improves convergence
    IF (ABS(this%dx).LE.this%tol) THEN
      this%dx = 0.8*this%tol*SIGN(1.0,this%x(ONE)-this%x(TWO))
    END IF
    ! new estimate for the root
    this%x(3) = this%x(TWO) + this%dx
    ! store this for error control
    this%dx = ABS(this%x(ONE)-this%x(TWO))
  END SUBROUTINE Step_SecantMethod


  PURE SUBROUTINE SaveInverseQuadraticInterpolation(a,b,c,fa,fb,fc,db_numer,db_denom)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,b,c,fa,fb,fc
    REAL, INTENT(OUT) :: db_numer,db_denom
    !------------------------------------------------------------------------!
    REAL :: fbdivfa,fadivfc,fbdivfc,p,q
    !------------------------------------------------------------------------!
    fbdivfa = fb / fa
    IF ((fa.EQ.fc).OR.(fb.EQ.fc)) THEN
        ! fallback: linear interpolation
        p = (c-b) * fbdivfa
        q = 1.0 - fbdivfa
    ELSE
        ! try inverse quadratic interpolation
        fadivfc = fa / fc
        fbdivfc = fb / fc
        p = fbdivfa * ((c-b)*fadivfc*(fadivfc-fbdivfc) - (b-a)*(fbdivfc-1.0) )
        q = (fadivfc - 1.0) * (fbdivfc - 1.0) * (fbdivfa - 1.0)
    END IF
    ! move the sign of the ratio p/q to the denominator q;
    ! thus the return value of the numerator is always > 0
    IF (p.GT.0.0) THEN
        db_numer = p
        db_denom = -q
    ELSE
        db_numer = -p
        db_denom = q
    END IF
  END SUBROUTINE SaveInverseQuadraticInterpolation


  PURE SUBROUTINE Step_ArithmeticMean(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    ! new estimate for the root
!CDIR IEXPAND
    this%x(3) = ArithmeticMean(this%x(1),this%x(2))
    ! store this for error control
    this%dx = ABS(this%x(2)-this%x(1))
  END SUBROUTINE Step_ArithmeticMean


  PURE FUNCTION ArithmeticMean(a,b) RESULT(s)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: a,b
    REAL             :: s
    !------------------------------------------------------------------------!
    s = 0.5*(a+b)
  END FUNCTION ArithmeticMean


  PURE FUNCTION LinearInterpolation(x1,x2,f1,f2) RESULT(dx)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2,f1,f2
    REAL             :: dx
    !------------------------------------------------------------------------!
    dx = (x1-x2)*f2 / (f2 - f1 + TINY(f1))
  END FUNCTION LinearInterpolation


  PURE SUBROUTINE TestConvergence(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Roots_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (this%error.NE.0) THEN
        ! error occurred => abort iteration
        this%iterate = .FALSE.
    ELSE
        ! check convergence
        this%tol = this%eps*ABS(this%x(3))
        IF (ABS(this%f(3)).LE.TINY(this%f(3)) .OR. this%dx.LE.this%tol) THEN
            ! converged => store current estimate xm in root and finish iteration
            this%iterate = .FALSE.
            this%root = this%x(3)
        ELSE
            ! not converged => check iteration limit
            IF (this%iter.GE.this%MAX_ITERATIONS) THEN
                ! exceeded => store error code and abort iteration
                this%error = 2
                this%iterate = .FALSE.
            ELSE
                ! continue iteration
                this%iter = this%iter + 1
                this%iterate = .TRUE.
            END IF
        END IF
    END IF
  END SUBROUTINE TestConvergence

END MODULE Roots
