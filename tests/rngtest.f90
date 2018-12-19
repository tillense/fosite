!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: rngtest.f03                                                       #
!#                                                                           #
!# Copyright (C) 2015-2018                                                   #
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
!> \test Check some known values of the random numbers
!! \author Manuel Jung
!!
!! We also test the basic quality of the real values from the DKiss64
!! generator.
!----------------------------------------------------------------------------!
PROGRAM rngtest
  USE rngs
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTEGER(KIND=I8)   :: i,x,x0,imax
  REAL               :: r,rmax,rmin,rnew
  !--------------------------------------------------------------------------!

  TAP_PLAN(6)

  ! Check Kiss64
  DO i=1, 100000000
    x = Kiss64()
  END DO

  x0 = 1666297717051644203_I8
  TAP_CHECK(x.EQ.x0,"Kiss64")

  imax = 100000000_I8
  r = 0.
  rmin = 1.
  rmax = 0.
  DO i=1, imax
   rnew = DKiss64()
   r = r + rnew
   rmin = MIN(rmin,rnew)
   rmax = MAX(rmax,rnew)
  END DO

  ! Check if the random numbers are in (0,1)
  ! and if the average is near 0.5
  TAP_CHECK_CLOSE(r/imax,0.5,1.E-4,"Average close to 0.5.")
  TAP_CHECK_GE(rmin,0.,"All are bigger (or equal) than 0.")
  TAP_CHECK_LE(rmax,1.,"All are smaller (or equal) than 1.")
  TAP_CHECK_CLOSE(rmin,0.,1.E-4,"Lower limit is close to 0.")
  TAP_CHECK_CLOSE(rmax,1.,1.E-4,"Upper limit is close to 1.")

  ! Check SuperKiss64
  TAP_DONE

END PROGRAM rngtest
