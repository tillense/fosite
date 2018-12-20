!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: factors.f90                                                       #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
!> module for prime factorization
!----------------------------------------------------------------------------!
MODULE Factors
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! number of prim numbers
  INTEGER, PARAMETER :: NUMPRIMS = 100
  ! array of first "NUMPRIMS" prime numbers;
  ! see http://primes.utm.edu/
  INTEGER, PARAMETER, DIMENSION(NUMPRIMS) :: PRIMS = &
       (/ 2,     3,     5,     7,    11,    13,    17,    19,    23,    29, &
         31,    37,    41,    43,    47,    53,    59,    61,    67,    71, &
         73,    79,    83,    89,    97,   101,   103,   107,   109,   113, &
        127,   131,   137,   139,   149,   151,   157,   163,   167,   173, &
        179,   181,   191,   193,   197,   199,   211,   223,   227,   229, &
        233,   239,   241,   251,   257,   263,   269,   271,   277,   281, &
        283,   293,   307,   311,   313,   317,   331,   337,   347,   349, &
        353,   359,   367,   373,   379,   383,   389,   397,   401,   409, &
        419,   421,   431,   433,   439,   443,   449,   457,   461,   463, &
        467,   479,   487,   491,   499,   503,   509,   521,   523,   541 /)
  ! allows for prime factorization up to 299208
  INTEGER, PARAMETER :: MAXNUM = 299208
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! some constants
       NUMPRIMS, MAXNUM, &
       ! methods
       GetFactor
  !--------------------------------------------------------------------------!

CONTAINS

  ! return the smallest prime factor of "number"
  ! p=0 for erroneous input
  PURE FUNCTION GetFactor(number) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN) :: number
    INTEGER :: p
    !------------------------------------------------------------------------!
    INTEGER :: i
    !------------------------------------------------------------------------!
    p = 0 ! for errornous input return 0
    IF ((number.GT.MAXNUM).OR.(number.LT.1)) RETURN
    DO i=1,NUMPRIMS
       p = MOD(number,PRIMS(i))
       IF (p.EQ.0) THEN
          p = PRIMS(i)
          RETURN
       END IF
    END DO
    p = number ! if there is no prim factor in the list, the number is
               ! (a pretty large) prime number itself
  END FUNCTION GetFactor

END MODULE Factors

