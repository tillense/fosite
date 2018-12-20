!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: specialfunctions.f90                                              #
!#                                                                           #
!# Copyright (C) 2015                                                        #
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
!> \test Check some known values of special functions
!! \author Manuel Jung
!!
!! The results have been generated with functions from scipy.special
!----------------------------------------------------------------------------!
PROGRAM specialfunctions
  USE functions
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  REAL, PARAMETER :: eps = 1.E-4
  REAL :: val
  !--------------------------------------------------------------------------!

  TAP_PLAN(10)

  ! scipy.special.iv(0,1.)
  val = Bessel_I0(1.)
  TAP_CHECK_CLOSE(val,1.2660658777520084,eps,"Bessel_I0(1.)")

  ! scipy.special.iv(0,10.)
  val = Bessel_I0(10.)
  TAP_CHECK_CLOSE(val,2815.7166284662544,eps,"Bessel_I0(10.)")

  ! scipy.special.iv(1,1.)
  val = Bessel_I1(1.)
  TAP_CHECK_CLOSE(val,0.56515910399248503,eps,"Bessel_I1(1.)")

  ! scipy.special.iv(1,10.)
  val = Bessel_I1(10.)
  TAP_CHECK_CLOSE(val,2670.9883037012542,eps,"Bessel_I1(10.)")

  ! scipy.special.kv(0,1.)
  val = Bessel_K0(1.)
  TAP_CHECK_CLOSE(val,0.42102443824070834,eps,"Bessel_K0(1.)")

  ! scipy.special.kv(0,10.)
  val = Bessel_K0(10.)
  TAP_CHECK_CLOSE(val,1.778006231616765e-05,eps,"Bessel_K0(10.)")

  ! scipy.special.kv(1,1.)
  val = Bessel_K1(1.)
  TAP_CHECK_CLOSE(val,0.60190723019723458,eps,"Bessel_K1(1.)")

  ! scipy.special.kv(1,10.)
  val = Bessel_K1(10.)
  TAP_CHECK_CLOSE(val,1.8648773453825585e-05,eps,"Bessel_K1(10.)")

  ! scipy.special.kve(0,1.)
  val = Bessel_K0e(1.)
  TAP_CHECK_CLOSE(val,1.1444630798068949,eps,"Bessel_K0e(1.)")

  ! scipy.special.kve(0,10.)
  val = Bessel_K0e(10.)
  TAP_CHECK_CLOSE(val,0.39163193443659866,eps,"Bessel_K0e(10.)")

  TAP_DONE

END PROGRAM specialfunctions
