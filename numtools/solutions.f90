!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: solutions.f90                                                     #
!#                                                                           #
!# Copyright (C) 2013 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!! This program solves the Riemann-Problem. For details please read
!! chapter 4 of "Riemann Solvers and Numerical Methods for Fluid Dynamics"
!! from E.F. Toro (Springer 1999, 2nd Ed.).
!----------------------------------------------------------------------------!

MODULE solutions
  USE roots, ONLY: GetRoot
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
    REAL :: rho_L,u_L,p_L,A_L,B_L,c_L,&
            rho_R,u_R,p_R,A_R,B_R,c_R,gamma
    REAL :: gam_p1, gam_m1, R0, Rt, Vmin, Vmax, &
            Vshock, n1, n2, n3, n4, n5
    REAL :: r, xi, vxi, gxi, zxi, zgxi, xacc
    REAL :: E0, rho0, P1
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    riemann, &
    sedov
  !--------------------------------------------------------------------------!
CONTAINS
  SUBROUTINE riemann(x0,gamma_,rho_l_,u_l_,p_l_,rho_r_,u_r_,p_r_,t,x,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL                :: gamma_,t,x0,rho_l_,u_l_,p_l_,rho_r_,u_r_,p_r_
    REAL,DIMENSION(:)   :: x
    REAL,DIMENSION(:,:) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i
    Real    :: rho_Lstar, rho_Rstar, p_star, u_star, &
               c_Lstar, c_Rstar,           &
               S_L, S_HL, S_TL, S_R, S_HR, S_TR
    INTEGER :: error
    !------------------------------------------------------------------------!
    gamma = gamma_ !1.4                  !heat capacity ratio

    rho_L = rho_l_
    rho_R = rho_R_

    u_L = u_l_
    u_R = u_r_

    p_L = p_l_
    p_R = p_r_

    ! sound speeds of far regions
    c_L = SQRT(gamma*p_L/rho_L)
    c_R = SQRT(gamma*p_R/rho_R)

    A_L = 2.0/((gamma+1.0)*rho_L)
    A_R = 2.0/((gamma+1.0)*rho_R)

    B_L = p_L*(gamma-1.0)/(gamma+1.0)
    B_R = p_R*(gamma-1.0)/(gamma+1.0)

    ! calculate p_star by findung the root of function f
    ! upper limit is just a somewhat big number. Is this really always
    ! sufficient?
    CALL GetRoot(f, 0., 2.E+3, p_star, error)
    !p_star = GetRoot_andbjo(f, 0., 1.E+4,1.E-6)

    u_star = 0.5*((u_L + u_R) + ( f_x(p_star, rho_R, u_R, p_R, A_R, B_R, c_R) &
                                  - f_x(p_star, rho_L, u_L, p_L, A_L, B_L, c_L)))

    IF (p_star > p_L) THEN
      !left shock
      rho_Lstar = rho_L * ( ((p_star/p_L) + (gamma-1.0)/(gamma+1.0)) &
                            / (((gamma-1.0)/(gamma+1.0))*(p_star/p_L) + 1.0))
    ELSE
      !left rarefaction
      rho_Lstar = rho_L * (p_star/p_L)**(1.0/gamma)
    END IF

    IF (p_star > p_R) THEN
      !right shock
      rho_Rstar = rho_R * ( ((p_star/p_R) + (gamma-1.0)/(gamma+1.0)) &
                           /(((gamma-1.0)/(gamma+1.0))*(p_star/p_R) + 1.0))
    ELSE
      !right rarefaction
      rho_Rstar = rho_R * (p_star/p_R)**(1.0/gamma)
    END IF

    ! sound speeds of star regions
    c_Lstar = c_L * (p_star/p_L)**((gamma-1.0)/(2.0*gamma))
    c_Rstar = c_R * (p_star/p_R)**((gamma-1.0)/(2.0*gamma))

    ! shock speeds
    S_L = u_L - c_L * SQRT( ((gamma+1.0)/(2.0*gamma))*(p_star/p_L) + (gamma-1.0)/(2.0*gamma))
    S_R = u_R + c_R * SQRT( ((gamma+1.0)/(2.0*gamma))*(p_star/p_R) + (gamma-1.0)/(2.0*gamma))

    ! head speed of rarefraction waves
    S_HL = u_L - c_L
    S_HR = u_R + c_R
    ! tail speed of rarefraction waves
    S_TL = u_star - c_Lstar
    S_TR = u_star + c_Rstar

    DO i=1,SIZE(x)
      CALL sample(x(i),pvar(i,1),pvar(i,2),pvar(i,3))
    END DO

    CONTAINS
      SUBROUTINE sample(x,rho,u,p)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL        :: x, rho, u, p
        !----------------------------------------------------------------------!
        REAL        :: S
        !----------------------------------------------------------------------!
        INTENT(IN)  :: x
        INTENT(OUT) :: rho, u, p
        !----------------------------------------------------------------------!

        !Sampling the solution
        S = (x-x0)/t
        IF (S <= u_star) THEN
          !left side of contact
          IF (p_star > p_L) THEN
            !left shock
            IF (S < S_L) THEN
              !far left region
              rho = rho_L
              u   = u_L
              p   = p_L
            ELSE
              !left star region
              rho = rho_Lstar
              u   = u_star
              p   = p_star
            END IF
          ELSE
            !left fan
            IF (S < S_HL) THEN
              !far left region
              rho = rho_L
              u   = u_L
              p   = p_L
            ELSE IF (S > S_TL) THEN
              !left star region
              rho = rho_Lstar
              u   = u_star
              p   = p_star
            ELSE
              !left fan region
              rho = rho_L * ( 2.0/(gamma+1.0) + ((gamma-1.0)/((gamma+1.0)*c_L)) &
                *(u_L - (x-x0)/t) ) ** (2.0/(gamma-1.0))
              u   = 2.0/(gamma+1.0) * (c_L + (gamma-1.0)*u_L/2.0 + (x-x0)/t)
              p   = p_L * ( 2.0/(gamma+1.0) + ((gamma-1.0)/((gamma+1.0)*c_L)) &
                * (u_L - (x-x0)/t) ) ** (2.0*gamma/(gamma-1.0))
            END IF
          END IF
        ELSE
          !right side of contact
          IF (p_star > p_R) THEN
            !right shock
            IF (S > S_R) THEN
              rho = rho_R
              u   = u_R
              p   = p_R
            ELSE
              rho = rho_Rstar
              u   = u_star
              p   = p_star
            END IF
          ELSE
            !right fan
            IF (S > S_HR) THEN
              rho = rho_R
              u   = u_R
              p   = p_R
            ELSE IF (S < S_TR) THEN
              rho = rho_Rstar
              u   = u_star
              p   = p_star
            ELSE
              rho = rho_R * ( 2.0/(gamma+1.0) - ((gamma-1.0)/((gamma+1.0)*c_R)) &
                * (u_R - (x-x0)/t) ) ** (2.0/(gamma-1.0))
              u = 2.0/(gamma+1.0) * ( -c_R + (gamma-1.0)*u_R/2.0 + (x-x0)/t)
              p = p_R * ( 2.0/(gamma+1.0) - ((gamma-1.0)/((gamma+1.0)*c_R)) &
                * (u_R - (x-x0)/t) ) ** (2.0*gamma/(gamma-1.0))
            END IF
          END IF
        END IF
      END SUBROUTINE sample
  END SUBROUTINE riemann

  SUBROUTINE sedov(gamma_, E0_, rho0_, P1_, time, dim, x, pvar)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    REAL, INTENT(IN)                    :: gamma_, E0_, rho0_, P1_, time
    INTEGER, INTENT(IN)                 :: dim
    REAL, DIMENSION(:), INTENT(IN)      :: x
    REAL, DIMENSION(:,:), INTENT(INOUT) :: pvar
    !----------------------------------------------------------------------!
    REAL, PARAMETER                     :: EPS = 1.0D-09
    INTEGER, PARAMETER                  :: MAXIT = 1.0D+08
    INTEGER                             :: i
    !----------------------------------------------------------------------!
    gamma = gamma_
    E0    = E0_
    rho0  = rho0_
    P1    = P1_
    gam_p1 = gamma + 1.0
    gam_m1 = gamma - 1.0

    SELECT CASE(dim)
    CASE(2)
      R0 = 1.0
      Rt = R0*(E0*time**2/rho0)**0.25
      Vmin = 1.0 / gamma
      Vmax = 2.0 / gamma
      Vshock = 0.5 *Rt/time
      n1 = -2.0
      n2 = 2.0*gam_m1/gamma
      n3 = 1.0/gamma
      n4 = -n1 / (2.0-gamma)
      n5 = -2.0 / (2.0-gamma)
    CASE(3)
      R0 = 1.033 ! for gamma=7/5, see Padmanabhan: Theo. Astro., Vol 1,  p.409
      Rt = R0*(E0*time**2/rho0)**0.2
      Vmin = 1./gamma
      Vmax = 5./(3.*gamma-1.)
      Vshock = 2.*Rt/(5.*time)
      n1 = -(13.*gamma**2.-7.*gamma+12.) / ((3.*gamma-1.)*(2.*gamma+1.))
      n2 = 5.*gam_m1 / (2.*gamma+1.)
      n3 = 3. / (2.*gamma+1.)
      n4 = -n1 / (2.-gamma)
      n5 = -2. / (2.-gamma)
    END SELECT
  ! accuracy
  xacc = EPS
  
  ! main loop
  DO i=1, SIZE(x)
     r = x(i)
     xi = r/Rt
     IF (xi.LE.1.0) THEN
        ! Newton-Raphson to solve the implicit equation
        IF (dim.EQ.2) THEN
           vxi = GetRoot_test(funcd2D,Vmin,Vmax*0.99,xacc)
           gxi = gam_p1/gam_m1*(gam_p1/gam_m1*(gamma*vxi-1.0))**n3 &
                *(gam_p1/2.0*(2.0-gamma*vxi))**n4 &
                *(gam_p1/gam_m1*(1.0-vxi))**n5
        ELSE
           vxi = GetRoot_test(funcd3D,Vmin,Vmax*0.99,xacc)
           gxi = gam_p1/gam_m1*(gam_p1/gam_m1*(gamma*vxi-1.))**n3 &
                *(gam_p1/(7.-gamma)*(5.-(3.*gamma-1.)*vxi))**n4 &
                *(gam_p1/gam_m1*(1.-vxi))**n5
        END IF
        zxi = gamma*gam_m1*(1.-vxi)*vxi*vxi/(2.*(gamma*vxi-1.))
!        IF (xi.GT.0.01) THEN
!           zgxi = 0.5*gam*gam_p1*(gam_p1/gam_m1)**(n3+n5) * vxi**2 &
!                *(gam*vxi-1.)**(n3-1.) * (1.-vxi)**(n5+1.) &
!                *(gam_p1/(7.-gam)*(5.-(3.*gam-1.)*vxi))**n4
!        ELSE
!           zgxi = 0.5*gam*gam_m1*(1.-vxi)*vxi**2 * xi**(5/n2*(n3-1.))
!        END IF
!        CALL funcd2D(vxi,fv,dfv)
!        WRITE (*,"(5(ES15.7))") xi,vxi,gxi,zxi,fv
!        STOP
     ELSE
        vxi = 0.
        gxi = 1.
        zxi = ((2.+dim)*time/(2.*r))**2 * gamma*P1/rho0
        zgxi= zxi
     END IF
     ! set primitive variables
     pvar(i,1) = rho0*gxi
     pvar(i,2) = Vshock*xi*vxi
     pvar(i,3) = pvar(i,1)/gamma*zxi*(2.*r/((2.+dim)*time))**2
!     pvar(i,3) = Control%rho0/gam*(2.*r/(5.*Control%time))**2 * zgxi
  END DO

  CONTAINS

  FUNCTION GetRoot_test(funcd,x1,x2,xacc) RESULT(root)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: x1,x2,xacc
    REAL :: root
    !------------------------------------------------------------------------!
    INTERFACE
       SUBROUTINE funcd(x,fx,dfx)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: x
         REAL, INTENT(OUT) :: fx,dfx
       END SUBROUTINE funcd
    END INTERFACE
    !------------------------------------------------------------------------!
    REAL    :: fm,dfm,fl,dfl,fr,dfr
    REAL    :: xm,xl,xr,dx
    INTEGER :: i
    !------------------------------------------------------------------------!
    ! compute left and right function values
    xl = MIN(x1,x2)
    xr = MAX(x1,x2)    
    CALL funcd(xl,fl,dfl)
    CALL funcd(xr,fr,dfr)
    ! check if root is within the interval [x1,x2]
    IF (fl*fr.GT.0.0) THEN
       WRITE (*,*) xl, xr,fl, fr,gamma, "Error: f(x1)*f(x2) should be < 0, aborting!"
       STOP
    END IF
    ! main loop
    DO i=1,MAXIT
  !     WRITE (*,"(4(ES20.12))") xl,fl,xr,fr
       ! regular falsi
       dx = fl*(xl-xr)/(fl-fr)
       xm = xl - dx
       root = xm
       CALL funcd(xm,fm,dfm)
       ! check abort criteron
       IF (ABS(fm).LT.xacc) THEN
          EXIT
       END IF
       IF (fm*fl.GT.0.0) THEN
          xl=xm
          fl=fm
       ELSE
          xr=xm
          fr=fm
       END IF
    END DO
    IF (i.GT.MAXIT) THEN
       WRITE (*,*) "WARNING: limit of iterations exceeded!"
    END IF
  END FUNCTION GetRoot_test

  END SUBROUTINE sedov


  ! the root of this function gives p_star
  !REAL FUNCTION f(p)
  PURE SUBROUTINE f(p,fx,plist)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    REAL, INTENT(IN)                         :: p
    REAL, INTENT(OUT)                        :: fx
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    !----------------------------------------------------------------------!
    fx = f_x(p, rho_L, u_L, p_L, A_L, B_L, c_L) &
       + f_x(p, rho_R, u_R, p_R, A_R, B_R, c_R) &
       + (u_R - u_L)
  END SUBROUTINE f

  !left side or right side function
  REAL PURE FUNCTION f_x(p, rho_x, u_x, p_x, A_x, B_x, c_x)
    IMPLICIT NONE
    !----------------------------------------------------------------------!
    REAL, INTENT(IN) :: p, rho_x, u_x, p_x, A_x, B_x, c_x
    !----------------------------------------------------------------------!
    IF (p > p_x) THEN
      f_x = (p-p_x) * SQRT(A_x/(p+B_x))
    ELSE
      f_x = (2.0*c_x/(gamma-1.0)) * ((p/p_x)**((gamma-1.0)/(2.0*gamma)) - 1.0)
    END IF
  END FUNCTION f_x

SUBROUTINE funcd2D(y,fy,dfy)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
  REAL, INTENT(IN)  :: y
  REAL, INTENT(OUT) :: fy,dfy
  !------------------------------------------------------------------------!
  REAL :: Ay,dAy,By,dBy,Cy,dCy
  !------------------------------------------------------------------------!
  
  Ay = 2./(gam_p1*y)
  dAy= -Ay/y
  By = gam_p1/2. * (2.-gamma*y)
  dBy= -gam_p1/2. * gamma
  Cy = gam_p1/gam_m1*(gamma*y-1.0)
  dCy= gamma*gam_p1/gam_m1
  fy = Ay**2 * By**n1 * Cy**n2 - xi**4
  dfy= dAy*By*Cy + Ay*dBy*Cy + Ay*By*dCy
!  PRINT "(A,4(ES14.6))", "xi,n1,n2       = ", xi,n1,n2
!  PRINT "(A,3(ES14.6))", "A(y),B(y),C(y) = ", Ay,By,Cy
!  PRINT "(A,3(ES14.6))", "y,f(y),df(y)   = ", y,fy,dfy
!  PRINT "(A,3(ES14.6))", "gamma,gam_p1,gam_m1   = ", gamma,gam_p1,gam_m1
END SUBROUTINE funcd2D


SUBROUTINE funcd3D(y,fy,dfy)
  IMPLICIT NONE
  !------------------------------------------------------------------------!
  REAL, INTENT(IN)  :: y
  REAL, INTENT(OUT) :: fy,dfy
  !------------------------------------------------------------------------!
  REAL :: Ay,dAy,By,dBy,Cy,dCy
  !------------------------------------------------------------------------!

    Ay = 2./(gam_p1*y)
    dAy= -Ay/y
    By = ABS(gam_p1/(7.-gamma)*(5.-(3*gamma-1.)*y))
    dBy= -gam_p1/(7.-gamma)*(3*gamma-1.)
    Cy = gam_p1/gam_m1*(gamma*y-1.)
    dCy= gamma*gam_p1/gam_m1
    IF (ABS(gamma*y-1.).LT.1.0D-20) THEN
       fy = -xi**5
    ELSE
       fy = Ay**2 * By**n1 * Cy**n2 - xi**5
    END IF
    fy = Ay**2 * By**n1 * Cy**n2 - xi**5
    dfy= dAy*By*Cy + Ay*dBy*Cy + Ay*By*dCy
!    PRINT "(A,ES14.6)","xi = ", xi
!    PRINT "(A,3(ES14.6))", "A(y),B(y),C(y) = ", Ay,By,Cy
!    PRINT "(A,3(ES14.6))", "y,f(y),df(y)   = ", y,fy,dfy
END SUBROUTINE funcd3D


END MODULE

