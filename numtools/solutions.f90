!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: riemannsolver.f90                                                 #
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
  !--------------------------------------------------------------------------!
  PUBLIC :: &
    riemann
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
END MODULE

