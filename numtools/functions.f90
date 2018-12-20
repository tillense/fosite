!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: functions.f90                                                     #
!#                                                                           #
!# Copyright (C) 2006-2013                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Tobias Illenseer, Manuel Jung, Björn Sperling
!!
!! \brief Mathematical functions.
!!
!! This module implements several mathematical functions, including
!! - inverse hyperbolic sine,cosine and tangens
!! - complete elliptic integrals of first and second kind
!! - error and inverse error function
!! - some Bessel functions of the first and second kind
!! - exponential integral function Ei(x)
!! - Gamma function
!----------------------------------------------------------------------------!
MODULE functions
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER    :: PI = 3.14159265358979323846
  REAL, PARAMETER    :: SQRT_TWO = 1.41421356237309504880
  REAL, PARAMETER    :: EPS_AGM = EPSILON(EPS_AGM)   ! precision of the AGM
  INTEGER, PARAMETER :: MAX_AGM_ITERATIONS = 20      ! limit iteration steps
  ! store the first 20 coefficients used in computing the Legendre polynomials
  INTEGER :: iii
  INTEGER, PARAMETER :: MAX_PL_COEFF = 20
  REAL, DIMENSION(MAX_PL_COEFF), PARAMETER &
       :: PL_COEFF = (/ (iii/(iii+1.0), iii=1,MAX_PL_COEFF) /)
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE LegendrePolynomial
     MODULE PROCEDURE LegendrePolynomial_one, LegendrePolynomial_all
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       Asinh, Acosh, Atanh, &
       Heaviside, &
       Barrier, &
       NormalRand, &
       EllipticIntegrals, &
       EllipticIntegral_K, &
       EllipticIntegral_E, &
       IncompleteEllipticIntegral_F, &
       IncompleteEllipticIntegral_E, &
       ErrorFunc,&
       CompErrorFunc,&
       InvErf,&
       LegendrePolynomials, &
       LegendrePolynomial, &
       LegendreFunction_QminHalf, &
       Bessel_I0, &
       Bessel_I1, &
       Bessel_K0, &
       Bessel_K0e, &
       Bessel_K1, &
       LnGamma, &
       Ei, &
       Gamma
  !--------------------------------------------------------------------------!

CONTAINS

  !> inverse hyperbolic sine function
  ! not included in FORTRAN before 2008 standard
  ELEMENTAL FUNCTION Asinh(x) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL :: fx
    fx = LOG(x+SQRT(x*x+1.0))
  END FUNCTION Asinh
  

  !> inverse hyperbolic cosine function
  ! not included in FORTRAN before 2008 standard
  ELEMENTAL FUNCTION Acosh(x) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL :: fx
    fx = LOG(x+SQRT(x*x-1.0))
  END FUNCTION Acosh
  

  !> inverse hyperbolic tangens function
  ! not included in FORTRAN before 2008 standard
  ELEMENTAL FUNCTION Atanh(x) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL :: fx
    fx = 0.5*LOG((1.0+x)/(1.0-x))
  END FUNCTION Atanh
  

  !> step function
  !      returns:   0   for x<a
  !                 1/2 for x=a
  !                 1   for x>a
  ELEMENTAL FUNCTION Heaviside(x,a) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x,a
    REAL :: fx
    fx = 0.5*(SIGN(1.0,x-a)+1.0)
  END FUNCTION Heaviside

  !> barrier function
  !         returns    0   for x<a and x>b
  !                    1/2 for x=a and x=b
  !                    1   for a<x<b
  ELEMENTAL FUNCTION Barrier(x,a,b) RESULT(fx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x,a,b
    REAL :: fx
    fx = Heaviside(x,a)-Heaviside(x,b)
  END FUNCTION Barrier

  !> Box-Muller alg. for gaussian distributed random variables
  !! input are uniform distributed random variables [0,1)
  ELEMENTAL SUBROUTINE NormalRand(u1,u2,x,y) 
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: u1,u2 
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    REAL    :: r,t
    !------------------------------------------------------------------------!
    r = SQRT(-2.0*LOG(u1))
    t = 2.0*PI*u2
    x = r*COS(t)
    y = r*SIN(t)
  END SUBROUTINE NormalRand


  !> compute the complete elliptic integrals of first and second kind
  !! using the AGM method (arithmetic-geometric-mean);
  !! 
  !! Function definitions and algorithm are given in
  !! [1] M. Abramowitz and I. A. Stegun: Handbook of Mathematical Functions,
  !!     Applied Mathematics Series. National Bureau of Standards, Vol. 55, 1964
  !!     online resource: http://people.math.sfu.ca/~cbm/aands/
  !! 
  !! check the value of K(k):
  !! (a) K((SQRT(6)-SQRT(2))/4) = 2**(-7/3) * 3**(1/4) * Gamma(1/3)**3 / PI
  !! (b) K(1/SQRT(2)) = SQRT(PI)/4 * Gamma(1/4)**2
  !! (c) K((SQRT(6)+SQRT(2))/4) = 2**(-7/3) * 3**(3/4) * Gamma(1/3)**3 / PI
  !!
  !! check value of E(k):
  !! (a) E((SQRT(6)-SQRT(2))/4) = SQRT(PI)*3**(-1/4) * (
  !!        2**(1/3) / SQRT(3) * (SQRT(PI)/Gamma(1/3))**3 
  !!      + 0.125*(SQRT(3)+1) / (2**(1/3)) * (Gamma(1/3)/SQRT(PI))**3)
  !! (b) E(1/SQRT(2)) = SQRT(PI)*(PI/(Gamma(1/4)**2)+0.125*(Gamma(1/4)**2)/PI)
  !! (a) E((SQRR(6)+SQRT(2))/4) = SQRT(PI)*3**(1/4) * (
  !!        2**(1/3) / SQRT(3) * (SQRT(PI)/Gamma(1/3))**3 
  !!      + 0.125*(SQRT(3)-1) / (2**(1/3)) * (Gamma(1/3)/SQRT(PI))**3)
  !!
  !! with Gamma(1/3) = 2.6789385347 ... and Gamma(1/4) = 3.6256099082 ...
  !!
  !! returns NaN on error, i.e. for |k| > 1
  ELEMENTAL SUBROUTINE EllipticIntegrals(k,Kell,Eell)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: k ! elliptic modulus
    REAL, INTENT(OUT) :: Eell,Kell
    !------------------------------------------------------------------------!
    INTEGER :: i,n
    REAL    :: an,bn,cn,tmp
    !------------------------------------------------------------------------!
    tmp = 0.0
    IF (ABS(k).LT.1.0) THEN
       i   = 1
       an  = 1.0
       bn  = SQRT(1.0-k*k)
       cn  = k
!CDIR UNROLL=20
       DO n=1,MAX_AGM_ITERATIONS ! do not loop forever
          tmp = tmp + cn*cn*i
          IF (ABS(cn).GT.0.5*EPS_AGM*ABS(an)) THEN
             cn = 0.5*(an-bn)
             bn = SQRT(an*bn)
             an = an-cn     ! = 0.5*(an+bn_old)
             i = ISHFT(i,1) ! = 2*i
#if !(defined(NECSXACE) || defined(NECSX9) || defined(NECSX8))
          ELSE
             EXIT  ! exit prohibits vectorization
#endif
          END IF
       END DO
       ! Kell = 0.5*PI / an better: Kell = 0.5*PI/0.5*(an+bn)
       Kell = PI / (an + bn + TINY(Kell)) ! avoid division by 0
    ELSE
       Kell = SQRT(-1.0*ABS(k)) ! return NaN
    END IF
    Eell = (1.0-0.5*tmp)*Kell
  END SUBROUTINE EllipticIntegrals
  

  !> compute the complete elliptic integral of the first kind
  !! using the AGM method (arithmetic-geometric-mean)
  !!
  !! returns NaN on error, i.e. for |k| > 1
  !! compute the complete elliptic integral of the first kind
  !! using the AGM method (arithmetic-geometric-mean)
  !!
  !! returns NaN on error, i.e. for |k| > 1
  ELEMENTAL FUNCTION EllipticIntegral_K(k) RESULT(Kell)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: k ! elliptic modulus
    REAL :: Kell
    !------------------------------------------------------------------------!
    REAL :: an,bn,tmp
    INTEGER :: i
    !------------------------------------------------------------------------!
    IF (ABS(k).LT.1.0) THEN
       an = 1.0
       bn = SQRT(1.0-k*k) ! = NaN if |k| > 1
!CDIR UNROLL=20
       DO i=1,MAX_AGM_ITERATIONS ! do not loop forever
          IF (ABS(an-bn).GT.EPS_AGM*ABS(an)) THEN
             tmp = 0.5*(an + bn)
             bn  = SQRT(an*bn)
             an  = tmp
#if !(defined(NECSXACE) || defined(NECSX9) || defined(NECSX8))
          ELSE
             EXIT  ! exit prohibits vectorization
#endif
          END IF
       END DO
       ! Kell = 0.5*PI / an -> next iteration would compute an = 0.5*(an+bn)
       ! hence return an even better approximation:
       Kell = PI / ( an + bn + TINY(Kell) ) ! avoid division by 0
    ELSE
       Kell = SQRT(-1.0*ABS(k)) ! return NaN       
    END IF
  END FUNCTION EllipticIntegral_K


  !> compute the complete elliptic integral of the second kind
  !! using the AGM method (arithmetic-geometric-mean)
  !!
  !! returns NaN on error, i.e. for |k| > 1
  ELEMENTAL FUNCTION EllipticIntegral_E(k) RESULT(Eell)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: k ! elliptic modulus 
    REAL :: Eell
    !------------------------------------------------------------------------!
    REAL :: dummy
    !------------------------------------------------------------------------!
    CALL EllipticIntegrals(k,dummy,Eell)
  END FUNCTION EllipticIntegral_E


  !> Computes Carlson's elliptic integral of the first kind R_F(x,y,z),
  !! y,y,z > 0. One can be zero.
  !! TINY must be at least 5x the machine underflow limit, BIG at most one 
  !! fifth of the machine overflow limit.
  !! similar to Numerical recipes 2nd edition
  ELEMENTAL FUNCTION rf(x,y,z) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL            :: x,y,z,res
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: ERRTOL = .0025
    REAL, PARAMETER :: THIRD = 1./3.
    REAL, PARAMETER :: TIN = 1.5E-38
    REAL, PARAMETER :: BIG = 3.E37
    REAL, PARAMETER :: C1 = 1./24.
    REAL, PARAMETER :: C2 = 0.1
    REAL, PARAMETER :: C3 = 3./44.
    REAL, PARAMETER :: C4 = 1./14.
    REAL            :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,&
                       yt,zt
    !------------------------------------------------------------------------!
    INTENT(IN)      :: x,y,z
    !------------------------------------------------------------------------!
    IF((MIN(x,y,z).LT.0.).OR.&
       (MIN(x+y,x+z,y+z).LT.TIN).OR.&
       (MAX(x,y,z).GT.BIG)) THEN
      res = SQRT(-1.*ABS(x)) ! return NAN
    ELSE
      xt = x
      yt = y
      zt = z
      DO
        sqrtx = SQRT(xt)
        sqrty = SQRT(yt)
        sqrtz = SQRT(zt)
        alamb = sqrtx * (sqrty+sqrtz) + sqrty*sqrtz
        xt = .25*(xt+alamb)
        yt = .25*(yt+alamb)
        zt = .25*(zt+alamb)
        ave = THIRD*(xt+yt+zt)
        delx = (ave-xt)/ave
        dely = (ave-yt)/ave
        delz = (ave-zt)/ave
        IF(MAX(ABS(delx),ABS(dely),ABS(delz)).LE.ERRTOL) &
          EXIT
      END DO
      
      e2 = delx*dely-delz**2
      e3 = delx*dely*delz
      res = (1. + (C1*e2-C2-C3*e3)*e2 + C4*e3)/SQRT(ave)
    END IF
  END FUNCTION rf
  

  !> Computes Carlson's elliptic integral of the second kind R_D(x,y,z),
  !! y,y > 0. One can be zero. z > 0.
  !! TINY must be at least twice -2/3 power of the machine underflow limit, 
  !! BIG at most 0.1 x ERRTOL times the -2/3 power of the machine overflow 
  !! limit.
  !! similar to Numerical recipes 2nd edition
  ELEMENTAL FUNCTION rd(x,y,z) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL            :: x,y,z,res
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: ERRTOL = .0015
    REAL, PARAMETER :: THIRD = 1./3.
    REAL, PARAMETER :: TIN = 1.5E-25
    REAL, PARAMETER :: BIG = 4.5E+21
    REAL, PARAMETER :: C1 = 3./14.
    REAL, PARAMETER :: C2 = 1./6.
    REAL, PARAMETER :: C3 = 9./22.
    REAL, PARAMETER :: C4 = 3./26.
    REAL, PARAMETER :: C5 = .25*C3
    REAL, PARAMETER :: C6 = 1.5*C4
    REAL            :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,&
                       sqrty,sqrtz,summ,xt,yt,zt
    !------------------------------------------------------------------------!
    INTENT(IN)      :: x,y,z
    !------------------------------------------------------------------------!
    IF((MIN(x,y).LT.0.).OR.&
       (MIN(x+y,z).LT.TIN).OR.&
       (MAX(x,y,z).GT.BIG)) THEN
      res = SQRT(-1.*ABS(x)) ! return NAN
    ELSE
      xt = x
      yt = y
      zt = z
      summ = 0.
      fac = 1.
      DO
        sqrtx = SQRT(xt)
        sqrty = SQRT(yt)
        sqrtz = SQRT(zt)
        alamb = sqrtx * (sqrty+sqrtz) + sqrty*sqrtz
        summ = summ + fac/(sqrtz*(zt+alamb))
        fac = .25*fac
        xt = .25*(xt+alamb)
        yt = .25*(yt+alamb)
        zt = .25*(zt+alamb)
        ave = .2*(xt+yt+3.*zt)
        delx = (ave-xt)/ave
        dely = (ave-yt)/ave
        delz = (ave-zt)/ave
        IF(MAX(ABS(delx),ABS(dely),ABS(delz)).LE.ERRTOL) &
          EXIT
      END DO
      ea = delx*dely
      eb = delz*delz
      ec = ea - eb
      ed = ea - 6. * eb
      ee = ed + ec + ec
      res = 3. * summ + fac*(1.+ed*(-C1+C5*ed-C6*delz*ee) &
        + delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*SQRT(ave))
    END IF
  END FUNCTION rd

  ! Legendre elliptic integral of the first kind F(phi,k), evaluated using
  ! Carlson's function R_F. The argument ranges are 0 <= phi <= pi/2, 
  ! 0 <= k*sin(phi) <= 1
  ! similar to Numerical recipes 2nd edition
  ELEMENTAL FUNCTION IncompleteEllipticIntegral_F(phi,ak) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL            :: phi,ak,res
    !------------------------------------------------------------------------!
    REAL            :: s
    !------------------------------------------------------------------------!
    INTENT(IN)      :: phi,ak
    !------------------------------------------------------------------------!
    s = SIN(phi)
    res = s*rf(COS(phi)**2,(1.-s*ak)*(1.+s*ak),1.)
  END FUNCTION IncompleteEllipticIntegral_F


  !> Legendre elliptic integral of the second kind E(phi,k), evaluated using
  !! Carlson's function R_D and R_F. The argument ranges are 0 <= phi <= pi/2, 
  !! 0 <= k*sin(phi) <= 1
  !! similar to Numerical recipes 2nd edition
  ELEMENTAL FUNCTION IncompleteEllipticIntegral_E(phi,ak) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL            :: phi,ak,res
    !------------------------------------------------------------------------!
    REAL            :: cc,q,s
    !------------------------------------------------------------------------!
    INTENT(IN)      :: phi,ak
    !------------------------------------------------------------------------!
    s = SIN(phi)
    cc = COS(phi)**2
    q = (1.-s*ak)*(1.+s*ak)
    res = s*(rf(cc,q,1.)-((s*ak)**2)*rd(cc,q,1.)/3.)
  END FUNCTION IncompleteEllipticIntegral_E


  PURE SUBROUTINE LegendrePolynomials(l,x,P)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER              :: l
    REAL                 :: x
    REAL, DIMENSION(0:l) :: P
    !------------------------------------------------------------------------!
    INTEGER              :: i
    !------------------------------------------------------------------------!
    INTENT(IN)           :: l,x
    INTENT(OUT)          :: P
    !------------------------------------------------------------------------!
    IF(l.LT.0) THEN
       P(:) = SQRT(1.0*l) ! return NaN
    ELSE
       P(0) = 1.0
       IF (l.GE.1) THEN
          P(1) = x
          ! use the constant coefficients for the first 
          ! MAX_PL_COEFF Legendre Polynomials
          DO i=2,MIN(l,MAX_PL_COEFF)
             P(i) = (1.0+PL_COEFF(i))*x*P(i-1) - PL_COEFF(i)*P(i-2)
          END DO
          ! from MAX_PL_COEFF+1 to l compute the coefficients i/(i+1);
          ! this is probably slower, because of the division
          DO i=MAX_PL_COEFF+1,l
             ! recurrence formula for the Legendre polynomials:
             ! P(i) = ( (2*i+1)*x*P(i-1) - i*P(i-2) ) / (i+1)
             !      = (1+i/(i+1))*x*P(i-1) - i/(i+1)*P(i-2)
             P(i) = i/(i+1.0) ! temporary (is a number close to 1)
             P(i) = (1.0+P(i))*x*P(i-1) - P(i)*P(i-2)
          END DO
       END IF
    END IF
  END SUBROUTINE LegendrePolynomials


  ELEMENTAL FUNCTION LegendrePolynomial_one(l,x,Plminus1,Plminus2) RESULT (Pl)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER    :: l
    REAL       :: x,Plminus1,Plminus2,Pl
    !------------------------------------------------------------------------!
    REAL       :: c
    !------------------------------------------------------------------------!
    INTENT(IN) :: l,x,Plminus1,Plminus2
    !------------------------------------------------------------------------!
    IF(l.LT.0) THEN
       Pl = SQRT(1.0*l) ! return NaN
    ELSE
       SELECT CASE(l)
       CASE(0)
          Pl = 1.0
       CASE(1)
          Pl = x
       CASE DEFAULT
          ! speed up things a little with predefined coefficients
          IF (l.LE.MAX_PL_COEFF) THEN
             c = PL_COEFF(l) 
          ELSE
             c = l/(l+1.0) ! temporary storage
          END IF
          ! do recursion step
          Pl = (1.0+c)*x*Plminus1 - c*Plminus2
       END SELECT
    END IF
  END FUNCTION LegendrePolynomial_one


  ELEMENTAL FUNCTION LegendrePolynomial_all(l,x) RESULT (Pl)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER        :: l
    REAL           :: x,Pl
    !------------------------------------------------------------------------!
    INTEGER        :: i
    REAL           :: Plminus1,Plminus2
    !------------------------------------------------------------------------!
    INTENT(IN)     :: l,x
    !------------------------------------------------------------------------!
    IF(l.LT.0) THEN
        Pl = SQRT(1.0*l) ! return NaN
    ELSE
!CDIR UNROLL=20
        DO i=0,l
!CDIR IEXPAND
           Pl = LegendrePolynomial_one(i,x,Plminus1,Plminus2)
           Plminus2 = Plminus1
           Plminus1 = Pl
        END DO
    END IF
  END FUNCTION LegendrePolynomial_all

  !> Computation of the order 0 and -1/2 degree
  !! Legendre function of the second kind Q_{-1/2}
  !! using the AGM method (arithmetic-geometric-mean)
  !!
  !! returns NaN on error, i.e. for x <= 1
  ELEMENTAL FUNCTION LegendreFunction_QminHalf(x) RESULT (q)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL       :: x,q
    !------------------------------------------------------------------------!
    INTEGER    :: i
    REAL       :: an,bn,tmp
    !------------------------------------------------------------------------!
    INTENT(IN) :: x
    !------------------------------------------------------------------------!
    IF (x.GE.1.0) THEN
       an = SQRT(x-1.0)
       bn = SQRT(x+1.0) ! if x < 1 bn = NaN, hence q = NaN
!CDIR UNROLL=20
       DO i=1,MAX_AGM_ITERATIONS ! do not loop forever
          IF (ABS(an-bn).GT.EPS_AGM*ABS(an)) THEN
             tmp = 0.5*(an + bn)
             bn = SQRT(an*bn)
             an = tmp
#if !(defined(NECSXACE) || defined(NECSX9) || defined(NECSX8))
          ELSE
             EXIT  ! exit prohibits vectorization
#endif
          END IF
       END DO
!       q = 0.5*PI*SQRT_TWO / an
       q = PI*SQRT_TWO / (an+bn+TINY(q))
    ELSE
       q = SQRT(-1.0*ABS(x)) ! return NaN
    END IF
  END FUNCTION LegendreFunction_QminHalf



!> Compute the modified Bessel function of the first kind as polynomial
!! expansion, which has been proposed by Numerical Recipes in fortran, Second
!! Edition on page 229ff. Nontheless the implementation is different and only the
!! idea is used.
  ELEMENTAL FUNCTION Bessel_I0(x) RESULT(I0)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: I0
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ 1.0d0, 3.5156229d0, 3.0899424d0, &
                                  1.2067492d0, 0.2659732d0, 0.360768d-1, &
                                  0.45813d-2 /)
    REAL, DIMENSION(9), PARAMETER &
                        :: q = (/ 0.39894228d0, 0.1328592d-1, 0.225319d-2, &
                                  -0.157565d-2, 0.916281d-2, -0.2057706d-1, &
                                  0.2635537d-1, -0.1647633d-1, 0.392377d-2 /)
    REAL                :: t, absx
    !------------------------------------------------------------------------!

    IF(ABS(x).LT.3.75) THEN
        t = (x/3.75)**2
        I0 = p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7))))))
    ELSE
        absx = ABS(x)
        t = 3.75/absx
        I0 = (EXP(absx)/sqrt(absx)) &
             * (q(1)+t*(q(2)+t*(q(3)+t*(q(4) &
                + t*(q(5)+t*(q(6)+t*(q(7)+t*(q(8)+t*q(9)))))))))
    ENDIF

    END FUNCTION Bessel_I0



!> Compute the modified Bessel function of the first kind as polynomial
!! expansion, which has been proposed by Numerical Recipes in fortran, Second
!! Edition on page 229ff. Nontheless the implementation is different and only the
!! idea is used.
  ELEMENTAL FUNCTION Bessel_I1(x) RESULT(I1)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: I1
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ 0.5d0, 0.87890594d0, 0.51498869d0, &
                                  0.15084934d0, 0.2658733d-1, 0.301532d-2, &
                                  0.32411d-3 /)
    REAL, DIMENSION(9), PARAMETER &
                        :: q = (/ 0.39894228d0, -0.3988024d-1, -0.362018d-2, &
                                  0.163801d-2, -0.1031555d-1, 0.2282967d-1, &
                                  -0.2895312d-1, 0.1787654d-1, -0.420059d-2 /)
    REAL                :: t, absx
    !------------------------------------------------------------------------!

    IF(ABS(x).LT.3.75) THEN
        t = (x/3.75)**2
        I1 = x*(p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7)))))))
    ELSE
        absx = ABS(x)
        t = 3.75/absx
        I1 = (EXP(absx)/sqrt(absx)) &
             * (q(1)+t*(q(2)+t*(q(3)+t*(q(4) &
                + t*(q(5)+t*(q(6)+t*(q(7)+t*(q(8)+t*q(9)))))))))
        IF(x.LT.0.) THEN
            I1 = -I1
        END IF

    ENDIF

    END FUNCTION Bessel_I1



!> Compute the modified Bessel function of the second kind as polynomial
!! expansion, which has been proposed by Numerical Recipes in fortran, Second
!! Edition on page 229ff. Nontheless the implementation is different and only
!! the idea is used.
  ELEMENTAL FUNCTION Bessel_K0(x) RESULT(K0)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: K0
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ -0.57721566d0, 0.42278420d0, 0.23069756d0, &
                                  0.3488590d-1, 0.262698d-2, 0.10750d-3, &
                                  0.74d-5 /)
    REAL, DIMENSION(7), PARAMETER &
                        :: q = (/ 1.25331414d0, -0.7832358d-1, 0.2189568d-1, &
                                  -0.1062446d-1, 0.587872d-2, -0.251540d-2, &
                                  0.53208d-3 /)
    REAL                :: t
    !------------------------------------------------------------------------!

    IF(x.LE.2.0) THEN
            t = x*x / 4.0
            K0 = (-LOG(x/2.0)*Bessel_I0(x)) &
                 + (p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7)))))))
    ELSE
            t = (2.0/x)
            K0 = (EXP(-x)/SQRT(x)) &
                 * (q(1)+t*(q(2)+t*(q(3)+t*(q(4)+t*(q(5)+t*(q(6)+t*q(7)))))))
    ENDIF

    END FUNCTION Bessel_K0

!> Compute the exponential scaled modified Bessel function of the second kind 
!! e.g. K0e = EXP(x) * K0(x) as polynomial expansion, using coefficients from
!! Abramowitz p.379 (http://people.math.sfu.ca/~cbm/aands/page_379.htm) for
!! x >= 2, and exp(x)*K0(x) directly for 0<x<2.
  ELEMENTAL FUNCTION Bessel_K0e(x) RESULT(K0e)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: K0e
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: q = (/ 1.25331414, -0.07832358, 0.02189568,&
                                  -0.01062446, 0.00587872, -0.00251540,&
                                  0.00053208 /)
    REAL                :: t
    !------------------------------------------------------------------------!

    IF(x.LT.2.0) THEN
            K0e = EXP(x) * Bessel_K0(x)
    ELSE
            t = (2.0/x)
            K0e = (1.0/SQRT(x)) &
                 * (q(1)+t*(q(2)+t*(q(3)+t*(q(4)+t*(q(5)+t*(q(6)+t*q(7)))))))
    ENDIF

    END FUNCTION Bessel_K0e

!> Compute the modified Bessel function of the second kind as polynomial
!! expansion, which has been proposed by Numerical Recipes in fortran, Second
!! Edition on page 229ff. Nontheless the implementation is different and only
!! the idea is used.
  ELEMENTAL FUNCTION Bessel_K1(x) RESULT(K1)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: x
    REAL                :: K1
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    REAL, DIMENSION(7), PARAMETER &
                        :: p = (/ 1.0d0, 0.15443144d0, -0.67278579d0, &
                                  -0.18156897d0, -0.1919402d-1, -0.110404d-2, &
                                  -0.4686d-4 /)
    REAL, DIMENSION(7), PARAMETER &
                        :: q = (/ 1.25331414d0, 0.23498619d0, -0.3655620d-1, &
                                  0.1504268d-1, -0.780353d-2, 0.325614d-2, &
                                  -0.68245d-3 /)
    REAL                :: t
    !------------------------------------------------------------------------!

    IF(x.LE.2.0) THEN
            t = x*x / 4.0
            K1 = (LOG(x/2.0)*Bessel_I1(x)) &
                 + (1.0/x) &
                    * (p(1)+t*(p(2)+t*(p(3)+t*(p(4)+t*(p(5)+t*(p(6)+t*p(7)))))))
    ELSE
            t = (2.0/x)
            K1 = (EXP(-x)/SQRT(x)) &
                 * (q(1)+t*(q(2)+t*(q(3)+t*(q(4)+t*(q(5)+t*(q(6)+t*q(7)))))))
    ENDIF

    END FUNCTION Bessel_K1


!> Computes the exponential integral Ei(x) for x != 0
!! Ei(x) := int_-oo^x exp(t)/t dt
!! Implementation similar to:
!! http://www.mymathlib.com/functions/exponential_integrals.html
  ELEMENTAL FUNCTION Ei(x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN):: x
    REAL            :: res
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: EPS=EPSILON(x), EULER=.577215664902, FPMIN=TINY(x)/EPS
    !------------------------------------------------------------------------!
    IF(x.LT.-5.) THEN
      res = Continued_Fraction_Ei(x)
    ELSE IF(x.EQ.0.) THEN
      res = -HUGE(res)
    ELSE IF(x.LT.6.8) THEN
      res = Power_Series_Ei(x)
    ELSE IF(x.LT.50.) THEN
      res = Argument_Addition_Series_Ei(x)
    ELSE
      res = Continued_Fraction_Ei(x)
    END IF
    END FUNCTION Ei


  ELEMENTAL FUNCTION Continued_Fraction_Ei(x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN):: x
    REAL            :: res
    !------------------------------------------------------------------------!
    REAL            :: Am1,A0,Bm1,B0,a,b,Ap1,Bp1,EPS
    INTEGER         :: j
    !------------------------------------------------------------------------!
    Am1 = 1.
    A0 = 0.
    Bm1 = 0.
    B0 = 1.
    a = EXP(x)
    b = -x + 1.
    Ap1 = b* A0 + a * Am1
    Bp1 = b* B0 + a * Bm1
    EPS = 10. * EPSILON(x) 
    j = 1
    a = 1.
    DO WHILE(ABS(Ap1 * B0 - A0 * Bp1).GT.EPS*ABS(A0*Bp1))
      IF(ABS(Bp1).GT.1.) THEN
        Am1 = A0 / Bp1
        A0 = Ap1 / Bp1
        Bm1 = B0 / Bp1
        B0 = 1.
      ELSE
        Am1 = A0
        A0 = Ap1
        Bm1 = B0
        B0 = Bp1
      END IF
      a = -j*j
      b = b + 2.
      Ap1 = b * A0 + a * Am1
      Bp1 = b * B0 + a * Bm1
      j = j + 1
    END DO
    res = -Ap1 / Bp1
  END FUNCTION Continued_Fraction_Ei


  ELEMENTAL FUNCTION Power_Series_Ei(x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN):: x
    REAL            :: res
    !------------------------------------------------------------------------!
    REAL            :: xn,Sn,Sm1,hsum,g,y,factorial,EPS
    !------------------------------------------------------------------------!
    xn = -x
    Sn = -x
    Sm1 = 0.
    hsum = 1.
    g = 0.5772156649015328606065121
    y = 1.
    factorial = 1.
    EPS = 10.*EPSILON(x)
    IF(x.EQ.0.) &
      res = -HUGE(res)
    DO WHILE(ABS(Sn-Sm1).GT.EPS*ABS(Sm1))
      Sm1 = Sn
      y = y + 1.
      xn = xn * (-x)
      factorial = factorial * y
      hsum = hsum + 1. / y
      Sn = Sn + hsum * xn / factorial
    END DO
    res = g + LOG(ABS(x)) - EXP(x) * Sn
  END FUNCTION Power_Series_Ei

  ELEMENTAL FUNCTION Argument_Addition_Series_Ei(x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN):: x
    REAL            :: res
    !------------------------------------------------------------------------!
    REAL, DIMENSION(44), PARAMETER :: eic = &
    (/1.915047433355013959531e2,  4.403798995348382689974e2, &
      1.037878290717089587658e3,  2.492228976241877759138e3, &
      6.071406374098611507965e3,  1.495953266639752885229e4, &
      3.719768849068903560439e4,  9.319251363396537129882e4, &
      2.349558524907683035782e5,  5.955609986708370018502e5, &
      1.516637894042516884433e6,  3.877904330597443502996e6, &
      9.950907251046844760026e6,  2.561565266405658882048e7, &
      6.612718635548492136250e7,  1.711446713003636684975e8, &
      4.439663698302712208698e8,  1.154115391849182948287e9, &
      3.005950906525548689841e9,  7.842940991898186370453e9, &
      2.049649711988081236484e10, 5.364511859231469415605e10,&
      1.405991957584069047340e11, 3.689732094072741970640e11,&
      9.694555759683939661662e11, 2.550043566357786926147e12,&
      6.714640184076497558707e12, 1.769803724411626854310e13,&
      4.669055014466159544500e13, 1.232852079912097685431e14,&
      3.257988998672263996790e14, 8.616388199965786544948e14,&
      2.280446200301902595341e15, 6.039718263611241578359e15,&
      1.600664914324504111070e16, 4.244796092136850759368e16,&
      1.126348290166966760275e17, 2.990444718632336675058e17,&
      7.943916035704453771510e17, 2.111342388647824195000e18,&
      5.614329680810343111535e18, 1.493630213112993142255e19,&
      3.975442747903744836007e19, 1.058563689713169096306e20 /)
    INTEGER         :: k, j
    REAL            :: xx,dx,xxj,edx,Sm,Sn,term,factorial,dxj,EPS
    !------------------------------------------------------------------------!
    k = NINT(x + 0.5)
    j = 0
    xx = k
    dx = x - xx
    xxj = xx
    edx = EXP(dx)
    Sm = 1.
    Sn = (edx - 1.) / xxj
    term = HUGE(x)
    factorial = 1.
    dxj = 1.
    EPS = 10.*EPSILON(x)
    
    DO WHILE(ABS(term).GT.EPS*ABS(Sn))
      j = j + 1
      factorial = factorial * j
      xxj = xxj * xx
      dxj = dxj * (-dx)
      Sm = Sm + dxj / factorial
      term = (factorial * (edx * Sm - 1.)) / xxj
      Sn = Sn + term
    END DO
    res = eic(k-6) + Sn * EXP(xx)
  END FUNCTION Argument_Addition_Series_Ei


  !> returns the error function of argument x 
  !! calculation of error fuction and complementary errorfunction are based on
  !! M. Abramowitz and I. A. Stegun: Handbook of Mathematical Functions,
  !! Applied Mathematics Series. National Bureau of Standards, Vol. 55, 1964
  !! online resource: http://people.math.sfu.ca/~cbm/aands/
  !! and 
  !! "Rational Chebyshev approximations for the error function"
  !! by W. J. Cody, Math. Comp., 1969, PP. 631-638.
  !! Coefficients for the Polynoms are taken from the library routine
  !! http://www.netlib.org/specfun/erf written by W. J. Cody
   ELEMENTAL FUNCTION ErrorFunc(x) RESULT(erf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x  
    REAL :: erf
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: XERF  = 0.46875
    REAL, PARAMETER :: XERFC = 4.0
    !------------------------------------------------------------------------!
    !Cody error function approximation
    IF (ABS(x) .LT. TINY(1.) ) THEN !!0 besser abfangen
            erf = 0.0
    ELSE IF(x.GT.0.0)THEN 
           IF (x .LT. XERF) THEN
               erf = CodyerfApprox(x)
           ELSE IF (x .LT. XERFC) THEN
               erf = 1. - Codycerf1Approx(x)
           ELSE
               erf = 1. - Codycerf2Approx(x)
           END IF
    !use erf,erfc symmetry for negative arguments
    ELSE 
           IF (-x .LT. XERF) THEN
               erf = - CodyerfApprox(-x)
           ELSE IF (-x .LT. XERFC) THEN
               erf = Codycerf1Approx(-x) - 1.
           ELSE
               erf = Codycerf2Approx(-x) - 1.
           END IF
    END IF

  END FUNCTION ErrorFunc

  ! returns the complementary error function 
  ELEMENTAL FUNCTION CompErrorFunc(x) RESULT(cerf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x  
    REAL :: cerf
    !------------------------------------------------------------------------!
    REAL, PARAMETER :: XERF  = 0.46875
    REAL, PARAMETER :: XERFC = 4.0
    !------------------------------------------------------------------------!
    !Cody complementary error function approximation
    IF(x.GE.0.0) THEN  
        IF (x .LT.XERF) THEN
            cerf = 1. - CodyerfApprox(x)
        ELSE IF (x .LT. XERFC) THEN
            cerf = Codycerf1Approx(x)
        ELSE
            cerf = Codycerf2Approx(x)
        END IF
    !use erf,erfc symmetry for negative arguments
    ELSE 
        IF (-x .LT. XERF) THEN
            cerf = 1. + CodyerfApprox(-x)
        ELSE IF (-x .LT. XERFC) THEN 
            cerf = 2. - Codycerf1Approx(-x)
        ELSE
            cerf = 2. - Codycerf2Approx(-x)
        END IF
    END IF

  END FUNCTION CompErrorFunc


  ELEMENTAL FUNCTION CodyerfApprox(x) RESULT(erf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x  
    REAL :: erf
    !------------------------------------------------------------------------!
    INTEGER :: i
    REAL :: denom, nom, x2, x4, x6
    REAL,DIMENSION(4) :: P_denom
    REAL,DIMENSION(5) :: P_nom
    !------------------------------------------------------------------------!
    P_nom   = (/ 3.20937758913846947E03 , 3.77485237685302021E02 ,&
                 1.13864154151050156E02 , 3.16112374387056560    ,&
                 1.85777706184603153E-1   /)
    P_denom = (/ 2.84423683343917062E03 , 1.28261652607737228E03 ,&
                 2.44024637934444173E02 , 2.36012909523441209E01 /)

    x2 = x  * x
    x4 = x2 * x2
    x6 = x4 * x2
    
    nom   = P_nom(1)  + x2 * (P_nom(2)   + P_nom(3)* x2   + P_nom(4) *   x4 &
            + P_nom(5)*x6) 
    denom = P_denom(1)+ x2 * (P_denom(2) + P_denom(3)* x2 + P_denom(4) * x4 &
            + x6)
     
    erf   = x * nom / denom

  END FUNCTION CodyerfApprox

  !> Cody approximation for complementary error function for 0.46875 < |x| < 4
  ELEMENTAL FUNCTION Codycerf1Approx(x) RESULT(cerf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x  
    REAL :: cerf
    !------------------------------------------------------------------------!
    INTEGER :: i
    REAL :: denom, nom
    REAL,DIMENSION(9) :: P_nom , P_denom
    !------------------------------------------------------------------------!
     P_nom   = (/ 1.23033935479799725E03 , 2.05107837782607147E03 ,&
                  1.71204761263407058E03 , 8.81952221241769090E02 ,&
                  2.98635138197400131E02 , 6.61191906371416295E01 ,&
                  8.88314979438837594E0  , 5.64188496988670089E-1 ,&
                  2.15311535474403846E-8 /)
     P_denom = (/ 1.23033935480374942E03 , 3.43936767414372164E03 ,&
                  4.36261909014324716E03 , 3.29079923573345963E03 ,&
                  1.62138957456669019E03 , 5.37181101862009858E02 ,&
                  1.17693950891312499E02 , 1.57449261107098347E01 ,&
                  1.0 /)

     nom   = 0.0
     denom = 0.0

     DO i= 0,7
        nom   = (nom   + P_nom(9-i))   * x
        denom = (denom + P_denom(9-i)) * x
     END DO
     
     nom   = nom   + P_nom(1)  
     denom = denom + P_denom(1)
    
     cerf = exp (-x*x)* nom / denom

  END FUNCTION Codycerf1Approx
  
  !> Cody approximation for complementary error function for |x| > 4
  ELEMENTAL FUNCTION Codycerf2Approx(x) RESULT(cerf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x  
    REAL :: cerf
    !------------------------------------------------------------------------!
    INTEGER :: i
    REAL :: denom, nom, x2, x4, x6, x8
    REAL,DIMENSION(5) :: P_denom
    REAL,DIMENSION(6) :: P_nom 
    !------------------------------------------------------------------------!
    P_nom   = (/ 6.58749161529837803E-4 , 1.60837851487422766E-2 ,&
                 1.25781726111229246E-1 , 3.60344899949804439E-1 ,&
                 3.05326634961232344E-1 , 1.63153871373020978E-2 /)    
    P_denom = (/ 2.33520497626869185E-3 , 6.05183413124413191E-2 ,&
                 5.27905102951428412E-1 , 1.87295284992346047E00 ,&
                 2.56852019228982242E00 /)

    x2 = x  * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4

    nom   = -x2 * P_nom(1) - P_nom(2) - P_nom(3)/x2 - P_nom(4)/(x4) - &
            P_nom(5)/(x6) - P_nom(6)/(x8)
    denom = x2 * P_denom(1) + P_denom(2) + P_denom(3)/x2 + P_denom(4)/(x4) + &
            P_denom(5)/(x6) + 1./(x8)

    cerf  = EXP(-x*x)/x * (1./SQRT(PI) +1./x2 * nom/denom)

  END FUNCTION Codycerf2Approx

  !> Inverse of the error function
  !! Argument of this function has to be in ]-1,1[, returns huge(1.0) for 
  !! values x with abs(x)>=1 without an error !!!! 
  !! uses and idea of P.J. Acklam 
  !! http://home.online.no/~pjacklam/notes/invnorm/ to calculate the inverse
  !! using the Halley root finding algorithm 
  !! Halley is preferable to Newton- Raphson as the second derivative can 
  !! be computed easily f''(y) = -2y*f'(y) for f(y)=erf(y)
  ELEMENTAL FUNCTION InvErf(x) RESULT(ierf)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x  
    REAL :: ierf
    !------------------------------------------------------------------------!
! FIXME : find global error variable, depending on precision(single,double,..)
    REAL, PARAMETER :: error = 1e-14
    REAL            :: q,y ,SQRTPI,SQRTPI_0_5,fy
    !------------------------------------------------------------------------!
     
     IF(x.LT.1.0)THEN
       SQRTPI = SQRT(PI)
       SQRTPI_0_5 = SQRTPI*0.5 

       ! initial guess: 
       y  = SQRTPI*(0.5*x + PI/24.*x**3 + 7.*PI**2/960.*x**5) 
       q  = SQRTPI_0_5 * EXP(y*y)*(ErrorFunc(y)-x)   !f(y)/f'(y)   
       fy = ErrorFunc(y) - x

       DO WHILE (ABS(fy) .GT. error)
          y = y - q/(1.+q*y)
          fy = ErrorFunc(y) - x
          q = SQRTPI_0_5 * EXP(y*y)*(ErrorFunc(y)-x)
       END DO
     ierf = y
     ELSE
     ! infinity
        ierf = HUGE(1.)
     END IF
     
   END FUNCTION InvErf

!> Compute the Gamma function
!! for arguments != 0,-1,-2,..
!!   x           Gamma(x)
!! ----------------------------
!!  1/3       2.678938534708
!!  0.5       1.772453850906
!! -0.5      -3.544907701811
!! -1.5       2.363271801207
!!  5.0      24.000000000000
  ELEMENTAL FUNCTION Gamma(x) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN):: x
    REAL            :: res
    !------------------------------------------------------------------------!
    REAL, DIMENSION(26), PARAMETER :: G = &
    (/ 1.0E0,                 0.5772156649015329E0,&
      -0.6558780715202538E0, -0.420026350340952E-1, &
       0.1665386113822915E0, -0.421977345555443E-1, &
      -0.96219715278770E-2,   0.72189432466630E-2, &
      -0.11651675918591E-2,  -0.2152416741149E-3, &
       0.1280502823882E-3,   -0.201348547807E-4, &
      -0.12504934821E-5,      0.11330272320E-5, &
      -0.2056338417E-6,       0.61160950E-8, &
       0.50020075E-8,        -0.11812746E-8, &
       0.1043427E-9,          0.77823E-11, &
      -0.36968E-11,           0.51E-12, &
      -0.206E-13,            -0.54E-14, &
       0.14E-14,              0.1E-15/)
    INTEGER         :: M1,K,M
    REAL            :: Z,R,GR
    !------------------------------------------------------------------------!
    IF (x.EQ.INT(x)) THEN
      IF (x.GT.0.) THEN
        res=1.0
        M1=x-1
        DO K=2,M1
          res=res*K
        END DO
      ELSE
        res=1.0E+300
      END IF
    ELSE
      IF (ABS(X).GT.1.) THEN
        Z=ABS(X)
        M=INT(Z)
        R=1.
        DO K=1,M
          R=R*(Z-K)
        END DO
        Z=Z-M
      ELSE
        Z=X
      ENDIF
      GR=G(26)
      DO K=25,1,-1
        GR=GR*Z+G(K)
      END DO
      res=1.0/(GR*Z)
      IF (ABS(X).GT.1.) THEN
        res=res*R
        IF (X.LT.0.) &
          res=-PI/(X*res*SIN(PI*X))
      ENDIF
    ENDIF
  END FUNCTION Gamma

!> computes the logarithm of the gammafunction with Lanczos approximation
!! found in Numerical Recipes for C++ p. 257 (different implementation)
!! coefficiants from Paul Godfrey for g=607/128 and N=15
!! http://my.fit.edu/~gabdo/gamma.txt
!! Gamma(x+1)=x!
 ELEMENTAL  FUNCTION LnGamma(x) RESULT(lg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    DOUBLE PRECISION    :: x
    REAL                :: lg
    DOUBLE PRECISION    :: alpha,beta,Ag,con
    INTEGER             :: i
    !------------------------------------------------------------------------!
    INTENT(IN)          :: x
    !------------------------------------------------------------------------!
    DOUBLE PRECISION, DIMENSION(15), PARAMETER &
                        :: p = &
          (/ .99999999999999709182d0,      57.156235665862923517d0,&
            -59.597960355475491248d0,      14.136097974741747174d0,&
            -.49191381609762019978d0,     .33994649984811888699d-4,&
            .46523628927048575665d-4,    -.98374475304879564677d-4,&
            .15808870322491248884d-3,    -.21026444172410488319d-3,&
            .21743961811521264320d-3,    -.16431810653676389022d-3,&
            .84418223983852743293d-4,    -.26190838401581408670d-4,&
            .36899182659531622704d-5                                /)

    ! prefactor 
    alpha = 2.506628274631000502415d0
    beta  = x + 5.24218750d0
    con = LOG(alpha)+LOG(beta)*(x+0.5d0) - beta

    ! sum
    IF (x.GE.0.0) THEN
      Ag = p(1)
      DO i=2,15
        Ag = Ag + p(i)/(x+i-1.d0)
      END DO
      lg  = con + LOG(Ag)
    ELSE
!      CALL Error(this,"LnGamma","bad value for for gammafunction") 
    END IF
  END FUNCTION LnGamma
END MODULE functions
