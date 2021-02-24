!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: rootstest.f90                                                     #
!#                                                                           #
!# Copyright (C) 2016                                                        #
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
!> \test Check solutions of root finding algorithms
!! \author Tobias Illenseer
!!
!! The examples are taken from the book \cite engeln2011 .
!----------------------------------------------------------------------------!
PROGRAM rootstest
  USE roots
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE
    PURE SUBROUTINE func(x,fx,plist)
      IMPLICIT NONE
      REAL, INTENT(IN)  :: x
      REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
      REAL, INTENT(OUT) :: fx
    END SUBROUTINE func
  END INTERFACE
  INTERFACE
    PURE SUBROUTINE funcd(x,fx,dfx,plist)
      IMPLICIT NONE
      REAL, INTENT(IN)  :: x
      REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
      REAL, INTENT(OUT) :: fx,dfx
    END SUBROUTINE funcd
  END INTERFACE
  COMMON /test_num/ k
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: NUM_TESTS = 13
  INTEGER, PARAMETER :: NUM_METHODS = 8
  CHARACTER(LEN=16), PARAMETER, DIMENSION(NUM_METHODS) :: method_name = (/ &
                                                       "       Bisection", &
                                                       "    Regula Falsi", &
                                                       "         Pegasus", &
                                                       "            King", &
                                                       " Anderson Bjoerk", &
                                                       "          Ridder", &
                                                       "    Brent Dekker", &
                                                       "          Newton" /)
  REAL, PARAMETER, DIMENSION(2,NUM_TESTS) :: bounds = RESHAPE((/ 0.0, 1.2, &
                                                                 0.4, 1.6, &
                                                                -0.5, 1.9, &
                                                                -0.5, 0.7, &
                                                                -1.4, 1.0, &
                                                                -0.8, 1.6, &
                                                                -0.5, 1.9, &
                                                               0.001, 1.201, &
                                                                -0.9, 1.5, &
                                                                 0.4, 1.0, &
                                                                -1.2, 0.0, &
                                                                 1.0, 3.4, &
                                                                 0.0, 5.0 /), &
                                             (/2,NUM_TESTS/))
  DOUBLE PRECISION, PARAMETER, DIMENSION(NUM_TESTS) :: ref_roots = &
     (/ 3.9942229171096819451D-01, 8.0413309750366432374D-01, 9.0340766319186021294D-01, &
      7.7014241346192677110D-02, 2.5920449372984746773D-01, 5.3674166257799978186D-01, &
      4.4754176206055907112D-01, 1.1111111111111111111D-01, 5.0000003403025908310D-01, &
      6.7980892150470050192D-01, -3.5938136638046273022D-01, 1.6487212707001281468D-00, &
      1.0000000000000000000D-00 &
     /)
  !--------------------------------------------------------------------------!
  REAL              :: root, xm, dx_rel, dx_acc, plist(1)
  INTEGER           :: i,k, iter, error
  CHARACTER(LEN=64) :: tap_message
  LOGICAL           :: verbose_results = .FALSE.
  !--------------------------------------------------------------------------!
  ! some information
  WRITE (*,"(A)") "====================================================================="
  WRITE (*,"(A)") "                << Testing root finding algorithms >>                "
  WRITE (*,"(A)") "====================================================================="
  
  
  TAP_PLAN((NUM_METHODS-1)*NUM_TESTS)

  DO k=1,NUM_TESTS
    ! plist is an array of real parameters that can be passed to the function func
    ! so far just a dummy, but it can be used if necessary;
    plist(1) = 1.0
    ! reset error code; if GetRoot returns error .ne. 0 something bad happend
    error = 0
    ! initial guess: arithmetic mean of boundary values
    xm = 0.5*(bounds(1,k)+bounds(2,k))
    PRINT '(A,I2,A,2(ES10.2))', "Test #",k,"  Search Interval: ",bounds(1,k),bounds(2,k)
    PRINT '(A)', "---------------------------------------------------------------------"
    IF (verbose_results) PRINT '(A20,ES27.19)', "reference result", ref_roots(k)
    DO i=1,NUM_METHODS
      SELECT CASE(i)
      CASE(1)
        CALL GetRoot_Bisection(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(2)
        ! don't use Regula Falsi because of slow convergence in some cases
        CALL GetRoot_RegulaFalsi(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(3)
        CALL GetRoot_Pegasus(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(4)
        CALL GetRoot_King(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(5)
        CALL GetRoot_AndersonBjoerk(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(6)
        CALL GetRoot_Ridder(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(7)
        CALL GetRoot_BrentDekker(func,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      CASE(8)
        CALL GetRoot_Newton(funcd,bounds(1,k),bounds(2,k),root,error,plist,xm,iter)
      END SELECT

      ! compute the relative  error using reference root
      dx_rel = ABS(1.0-root/ref_roots(k))
      IF (verbose_results) CALL PrintResults(method_name(i),root,ref_roots(k),error,iter)
      ! exclude Regular Falsi from check
      IF (i.NE.2) THEN
      
        dx_acc = DEFAULT_ACCURACY/ABS(ref_roots(k))
        IF (k.EQ.12) THEN
          ! this test never yields the required accuracy, see \cite engeln2011
          dx_acc = 1.0E-5
        END IF
        WRITE (tap_message,'(A16,A11,ES9.2,A3,ES9.2,A7,I4)') method_name(i), ": dx_rel = ", &
                dx_rel, " < ", dx_acc, " iter = ", iter
        TAP_CHECK(dx_rel.LE.dx_acc,tap_message)
        
      END IF
    END DO
    PRINT '(A)', "---------------------------------------------------------------------"
  END DO

  TAP_DONE

CONTAINS


  SUBROUTINE PrintResults(method,root,ref_root,error,iter)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: method
    REAL, INTENT(IN) :: root
    DOUBLE PRECISION, INTENT(IN) :: ref_root
    INTEGER, INTENT(IN) :: error,iter
    IF (error.NE.0) THEN
        PRINT *,TRIM(GetErrorMessage(error))
    ELSE
        PRINT '(A20,ES23.15,ES9.1,I4)',TRIM(method),root,ABS(1.0-root/ref_root),iter
    END IF
  END SUBROUTINE PrintResults


!   TAP_PLAN(7)
! 
!   
! 
!   ! Check if the random numbers are in (0,1)
!   ! and if the average is near 0.5
!   TAP_CHECK_CLOSE(r/imax,0.5,1.E-4,"Average close to 0.5.")
!   TAP_CHECK_GE(rmin,0.,"All are bigger (or equal) than 0.")
!   TAP_CHECK_LE(root,1.,"All are smaller (or equal) than 1.")
!   TAP_CHECK_CLOSE(rmin,0.,1.E-4,"Lower limit is close to 0.")
!   TAP_CHECK_CLOSE(rmax,1.,1.E-4,"Upper limit is close to 1.")
! 
!   TAP_CHECK(x.EQ.x0,"SuperKiss64")
! 

END PROGRAM rootstest


! function definitions
PURE SUBROUTINE func(x,fx,plist)
IMPLICIT NONE
  COMMON /test_num/ k
  REAL, INTENT(IN)  :: x
  REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
  REAL, INTENT(OUT) :: fx
  INTEGER :: k
  SELECT CASE(k)
  CASE(1)
    fx = x*x*(x*x/3. + SQRT(2.)*SIN(x)) - SQRT(3.)/18.
  CASE(2)
    fx = 11*x**11 - 1.0
  CASE(3)
    fx = 35*x**35 - 1.0
  CASE(4)
    fx = 2*(x*EXP(-9.)-EXP(-9*x)) + 1.
  CASE(5)
    fx = x*x - (1.-x)**9
  CASE(6)
    fx = (x-1.)*EXP(-9*x) + x**9
  CASE(7)
    fx = x*x + SIN(x/9.) - 0.25
  CASE(8)
    fx = 0.125*(9.-1./x)
  CASE(9)
    fx = TAN(x) - x - 0.0463025
  CASE(10)
    fx = x*(x+SIN(SQRT(75.)*x)) - 0.2
  CASE(11)
    fx = x**9 + 1e-4
  CASE(12)
    fx = LOG(x) + 0.5*x*x/EXP(1.) -2*x/SQRT(EXP(1.)) + 1.
  CASE DEFAULT
    fx = x - 1.0
  END SELECT
END SUBROUTINE func


PURE SUBROUTINE funcd(x,fx,dfx,plist)
  IMPLICIT NONE
  COMMON /test_num/ k
  REAL, INTENT(IN)  :: x
  REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
  REAL, INTENT(OUT) :: fx,dfx
  INTEGER :: k
  SELECT CASE(k)
  CASE(1)
    fx = x*x*(x*x/3. + SQRT(2.)*SIN(x)) - SQRT(3.)/18.
    dfx= 2*x*(SQRT(2.)*SIN(x)+x*x/3.)+x*x*(SQRT(2.)*COS(x)+(2*x)/3.)
  CASE(2)
    fx = 11*x**11 - 1.0
    dfx = 121*x**10
  CASE(3)
    fx = 35*x**35 - 1.0
    dfx= 35*35*x**34 - 1.0
  CASE(4)
    fx = 2*(x*EXP(-9.)-EXP(-9*x)) + 1.
    dfx= 2*(9*EXP(-9*x)+EXP(-9.))
  CASE(5)
    fx = x*x - (1.-x)**9
    dfx= 2*x + 9*(1.-x)**8
  CASE(6)
    fx = (x-1.)*EXP(-9*x) + x**9
    dfx= (10.-9*x)*EXP(-9*x) + 9*x**8
  CASE(7)
    fx = x*x + SIN(x/9.) - 0.25
    dfx= 2*x + COS(x/9.)/9.
  CASE(8)
    fx = 0.125*(9.-1./x)
    dfx= 0.125/x**2
  CASE(9)
    fx = TAN(x) - x - 0.0463025
    dfx= 1./COS(x)**2 - 1.
  CASE(10)
    fx = x*(x+SIN(SQRT(75.)*x)) - 0.2
    dfx= SIN(SQRT(75.)*x) + x*(2. + SQRT(75.)*COS(SQRT(75.)*x))
  CASE(11)
    fx = x**9 + 1e-4
    dfx= 9*x**8
  CASE(12)
    fx = LOG(x) + 0.5*x*x/EXP(1.) - 2*x/SQRT(EXP(1.)) + 1.
    dfx= 1./x + x/EXP(1.) - 2./SQRT(EXP(1.))
  CASE DEFAULT
    fx = x - 1.0
    dfx= 1.
  END SELECT
END SUBROUTINE funcd
  
