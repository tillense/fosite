!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: tap.f90                                                           #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Manuel Jung <mjunf@astrophysik.uni-kiel.de>                               #
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
!> TAP: Test Anything Protocoll
!! For a short definition visit
!! http://podwiki.hexten.net/TAP/TAP13.html?page=TAP13
!----------------------------------------------------------------------------!

MODULE tap
  USE logging_base_mod, ONLY : SetPrefix
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
    INTEGER, PARAMETER :: NO_PLAN  = -1
    INTEGER, PARAMETER :: SKIP_ALL = -2
    INTEGER, SAVE :: expected_tests = NO_PLAN
    INTEGER, SAVE :: failed_tests = 0
    INTEGER, SAVE :: current_test = 0
  !--------------------------------------------------------------------------!
  PUBLIC :: &
        tap_plan, &
        tap_done, &
        tap_check_at_loc, &
        tap_check_op_at_loc, &
        tap_check_close_at_loc, &
        tap_check_small_at_loc
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE tap_diag(str)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)   :: str
    !------------------------------------------------------------------------!
    INTENT(IN)         :: str
    !------------------------------------------------------------------------!
    WRITE(*,'(A,A)') "#  ",TRIM(str)
  END SUBROUTINE tap_diag


  FUNCTION itoa(i) RESULT(str)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER         :: i
    CHARACTER(LEN=CEILING(LOG10(REAL(i)+1))):: str
    !------------------------------------------------------------------------!
    INTEGER         :: n
    CHARACTER(LEN=16):: tmp
    !------------------------------------------------------------------------!
    INTENT(IN)      :: i
    !------------------------------------------------------------------------!
    n = CEILING(LOG10(REAL(i)+1))
    WRITE(tmp,'(I16)') i
    str = tmp(16-n+1:16)
  END FUNCTION itoa


  FUNCTION rtoa(r) RESULT(str)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL            :: r
    CHARACTER(LEN=12):: str
    !------------------------------------------------------------------------!
    INTENT(IN)      :: r
    !------------------------------------------------------------------------!
    WRITE(str,'(ES12.3)') r
  END FUNCTION rtoa


  SUBROUTINE tap_plan(tests, why)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER     :: tests
    CHARACTER(len=128),OPTIONAL :: why
    !------------------------------------------------------------------------!
    INTENT(IN)  :: tests, why
    !------------------------------------------------------------------------!
    CALL SetPrefix('#')
    expected_tests = tests
    failed_tests = 0
    current_test = 0
    WRITE(*, '(A)') "TAP version 13"
    SELECT CASE(tests)
    CASE(SKIP_ALL)
      IF(PRESENT(why)) THEN
         WRITE(*, '(A,A)') "1..0 SKIP ", why
      ELSE
         WRITE(*, '(A)') "1..0 SKIP"
      END IF
    CASE(NO_PLAN)
    CASE DEFAULT
      WRITE(*, '(A,A)') "1..", itoa(tests)
    END SELECT
  END SUBROUTINE tap_plan


  SUBROUTINE tap_done()
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    IF(expected_tests.EQ.NO_PLAN) THEN
      WRITE(*,'(A,A)') "1..",itoa(current_test)
      STOP 0
    ELSE IF(current_test.NE.expected_tests) THEN
      CALL tap_diag("Looks like " // itoa(expected_tests) //&
        " tests were planned but " // itoa(current_test) // " ran.")
      STOP 255
    END IF
    IF(failed_tests.GT.0) THEN
      CALL tap_diag("Looks like " // itoa(failed_tests) // &
        " tests of " // itoa(current_test) // " have failed.")
      !IF(expected_tests.EQ.NO_PLAN) THEN
      !  retval = failed_tests
      !ELSE
      !  retval = expected_tests - current_test + failed_tests
      !END IF
      ! STOP retval ...not possible in FORTRAN
      STOP 255
    END IF
  END SUBROUTINE tap_done


  SUBROUTINE tap_check(file,line,test,name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)    :: file
    INTEGER             :: line
    LOGICAL             :: test
    CHARACTER(LEN=*)    :: name
    !------------------------------------------------------------------------!
    INTENT(IN)          :: file, line, test, name
    !------------------------------------------------------------------------!
    current_test = current_test + 1
    IF(test) THEN
      WRITE(*,'(A,A,A,A)') "ok ",itoa(current_test)," - ",TRIM(name)
    ELSE
      failed_tests = failed_tests + 1
      WRITE(*,'(A,A,A,A)') "not ok ",itoa(current_test)," - ",TRIM(name)
    END IF
  END SUBROUTINE tap_check


  SUBROUTINE tap_check_at_loc(file,line,test,name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)    :: file
    INTEGER             :: line
    LOGICAL             :: test
    CHARACTER(LEN=*)    :: name
    !------------------------------------------------------------------------!
    INTENT(IN)          :: file, line, test, name
    !------------------------------------------------------------------------!
    CALL tap_check(file,line,test,name)
    IF(.NOT.test) THEN
      CALL tap_diag("Failed test '"//TRIM(name)//"'")
      CALL tap_diag("at "//TRIM(file)//":"//itoa(line)//".")
    END IF
  END SUBROUTINE tap_check_at_loc


  SUBROUTINE tap_check_op_at_loc(file,line,op,test,a,b,name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)    :: file
    INTEGER             :: line
    CHARACTER(LEN=*)    :: op
    LOGICAL             :: test
    REAL                :: a, b
    CHARACTER(LEN=*)    :: name
    !------------------------------------------------------------------------!
    INTENT(IN)          :: file, line, op, test, a, b, name
    !------------------------------------------------------------------------!
    CALL tap_check(file,line,test,name)
    IF(.NOT.test) THEN
      CALL tap_diag("Failed test '"//TRIM(name)//"'")
      CALL tap_diag(rtoa(a)//" "//TRIM(op)//" "//rtoa(b))
      CALL tap_diag("at "//TRIM(file)//":"//itoa(line)//".")
    END IF
  END SUBROUTINE tap_check_op_at_loc


  SUBROUTINE tap_check_close_at_loc(file,line,a,b,eps,name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)    :: file
    INTEGER             :: line
    REAL                :: a, b, eps
    CHARACTER(LEN=*)    :: name
    !------------------------------------------------------------------------!
    LOGICAL             :: test
    !------------------------------------------------------------------------!
    INTENT(IN)          :: file, line, a, b, eps, name
    !------------------------------------------------------------------------!
    test = ABS(a-b).LT.eps
    CALL tap_check(file,line,test,name)
    IF(.NOT.test) THEN
      CALL tap_diag("Failed test '"//TRIM(name)//"'")
      CALL tap_diag("ABS("//rtoa(a)//"-"//rtoa(b)//") .LT. "//rtoa(eps))
      CALL tap_diag("at "//TRIM(file)//":"//itoa(line)//".")
    END IF
  END SUBROUTINE tap_check_close_at_loc


  SUBROUTINE tap_check_small_at_loc(file,line,a,eps,name)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CHARACTER(LEN=*)    :: file
    INTEGER             :: line
    REAL                :: a, eps
    CHARACTER(LEN=*)    :: name
    !------------------------------------------------------------------------!
    LOGICAL             :: test
    !------------------------------------------------------------------------!
    INTENT(IN)          :: file, line, a, eps, name
    !------------------------------------------------------------------------!
    test = ABS(a).LT.eps
    CALL tap_check(file,line,test,name)
    IF(.NOT.test) THEN
      CALL tap_diag("Failed test '"//TRIM(name)//"'")
      CALL tap_diag("ABS("//rtoa(a)//") .LT. "//rtoa(eps))
      CALL tap_diag("at "//TRIM(file)//":"//itoa(line)//".")
    END IF
  END SUBROUTINE tap_check_small_at_loc

END MODULE Tap
