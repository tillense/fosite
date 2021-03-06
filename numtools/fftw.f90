!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fftw.f90                                                          #
!#                                                                           #
!# Copyright (C) 2011 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!> fftw module
!----------------------------------------------------------------------------!

MODULE fftw
  USE, INTRINSIC :: iso_c_binding
#if defined(HAVE_FFTW) && !defined(PARALLEL)
#ifdef NECSXAURORA
  INCLUDE 'aslfftw3.f03'
#else
  INCLUDE 'fftw3.f03'
#endif
#elif defined(HAVE_FFTW) && defined(PARALLEL)
#ifdef NECSXAURORA
  INCLUDE 'aslfftw3-mpi.f03'
#else
  INCLUDE 'fftw3-mpi.f03'
#endif
#endif
END MODULE fftw
