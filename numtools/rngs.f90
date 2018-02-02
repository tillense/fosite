!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: rngs.f90                                                          #
!#                                                                           #
!# Copyright (C) 2015 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!! random number generators
!! Kiss64: Period: >2^124 \approx 2.1*10^37
!!   Source: https://de.wikipedia.org/wiki/KISS_(Zufallszahlengenerator)
!!           http://fortranwiki.org/fortran/show/kiss64
!! DKiss64: Kiss64 converted to double intervall [0,1]
!! SuperKiss64: A Super KISS. Period of 5*2^1320480*(2^64-1)
!!   A KISS (Keep-It-Simple-Stupid) RNG combining,
!!   by addition mod 2^32, three simple RNGs:
!!     CMWC (Complementary-Multiply-With-Carry)
!!   + CNG (Congruential)
!!   + XS(Xorshift)
!!   with resulting period greater than 10^402575
!!   Source: http://mathforum.org/kb/message.jspa?messageID=6914945
!----------------------------------------------------------------------------!
MODULE rngs
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
#if defined(NECSX8) || defined(NECSX9) || defined(NECSXACE)
  INTEGER, PARAMETER :: I8 = 8
#else
  INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND(18)
#endif
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       I8, &
       Kiss64, &
       SuperKiss64, &
       DKiss64
  !--------------------------------------------------------------------------!

CONTAINS

  FUNCTION Kiss64(seed) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER(KIND=I8)           :: res
    INTEGER(KIND=I8), OPTIONAL :: seed
    !------------------------------------------------------------------------!
    INTEGER(KIND=I8), save     :: x, y, z, c
    INTEGER(KIND=I8)           :: t, k, m, s
    data x, y, z, c &
      / 1234567890987654321_I8, &
        362436362436362436_I8, &
        1066149217761810_I8, &
        123456123456123456_I8 /
    !------------------------------------------------------------------------!
    IF(PRESENT(seed)) THEN
      z = seed
      res = 0
    ELSE
      t = ISHFT(x, 58) + c
      IF (ISHFT(x,-63_I8) .EQ. ISHFT(t,-63_I8)) THEN
        c = ISHFT(x, -6) + ISHFT(x, -63_I8)
      ELSE
        c = ISHFT(x, -6) + 1 - ISHFT(x+t, -63_I8)
      END IF
      x = t + x
      y = IEOR(y, ISHFT(y,13_I8))
      y = IEOR(y, ISHFT(y,-17_I8))
      y = IEOR(y, ISHFT(y,43_I8))

      z = 6906969069_I8 * z + 1234567
      res = x + y + z
    END IF
  END FUNCTION Kiss64

  FUNCTION SuperKiss64() RESULT(x)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER(KIND=I8)       :: x
    !------------------------------------------------------------------------!
    INTEGER(KIND=I8), save :: i,Q(20632),carry,xcng,xs,indx
    data i,carry,xcng,xs,indx &
      / 0_I8, &
        36243678541_I8, &
        12367890123456_I8,&
        521288629546311_I8,&
        20633_I8 /
    !------------------------------------------------------------------------!
! Nec SX internal! compiler error:
! "f90 fatal: Internal error in optimization phase."
#if !(defined(NECSX8) || defined(NECSX9) || defined(NECSXACE))
    !Initalize if i=0
    IF(i.EQ.0) THEN
      DO i=1,20632
        !fill Q with Congruential+Xorshift
        xcng = xcng*6906969069_I8+123
        xs = IEOR(xs,ISHFT(xs,13))
        xs = IEOR(xs,ISHFT(xs,-17))
        xs = IEOR(xs,ISHFT(xs,43))
        Q(i) = xcng + xs
      END DO
    END IF

    IF(indx <= 20632) THEN
      x=Q(indx)
      indx=indx+1
    ELSE
      x=refill()
    END IF

    xcng=xcng*6906969069_I8+123
    xs=IEOR(xs,ISHFT(xs,13))
    xs=IEOR(xs,ISHFT(xs,-17))
    xs=IEOR(xs,ISHFT(xs,43))
    x=x+xcng+xs
  CONTAINS

    FUNCTION refill() RESULT(s)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      INTEGER(KIND=I8) :: i,s,z,h
      !------------------------------------------------------------------------!
      DO i=1,20632
        h=IAND(carry,1_I8)
        z = ISHFT(ISHFT(Q(i),41),-1) &
            + ISHFT(ISHFT(Q(i),39),-1) &
            + ISHFT(carry,-1)
        carry=ISHFT(Q(i),-23)+ISHFT(Q(i),-25)+ISHFT(z,-63)
        Q(i)=NOT(ISHFT(z,1)+h)
      END DO
      indx=2
      s=Q(1)
      RETURN
    END FUNCTION refill
#else
    x = 0_I8
#endif
  END FUNCTION SuperKiss64

  !> Convert integer I8 random numbers to double precision real values in [0,1]
  FUNCTION I8toD(x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER(KIND=I8) :: x
    REAL             :: res
#if defined(NECSX8) || defined(NECSX9) || defined(NECSXACE)
    INTEGER(KIND=8)  :: long
#endif
    !------------------------------------------------------------------------!
#if defined(NECSX8) || defined(NECSX9) || defined(NECSXACE)
    ! No idea what is happeing on the SX.., but HUGE(x)!=HUGE(long)
    ! The 2048? This is 2.*1024. No idea, why the 1024 is needed. This shouldn't
    ! be required.
    res = (REAL(x) / HUGE(long)) / 2048. + 0.5
#else
    res = 0.5E+0 * (REAL(x) / REAL(HUGE(x))) + 0.5E+0
#endif
  END FUNCTION I8toD


  FUNCTION DKiss64() RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL             :: res
    !------------------------------------------------------------------------!
    res = I8toD(Kiss64())
  END FUNCTION DKiss64


END MODULE rngs
