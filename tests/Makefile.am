#                                               -*- Automake -*-
# Process this file with automake --add-missing to produce a Makefile.
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# configuration file: Makefile.am                                           #
#                                                                           #
# Copyright (C) 2012                                                        #
# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
#                                                                           #
# This program is free software; you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 2 of the License, or (at     #
# your option) any later version.                                           #
#                                                                           #
# This program is distributed in the hope that it will be useful, but       #
# WITHOUT ANY WARRANTY; without even the implied warranty of                #
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
# NON INFRINGEMENT.  See the GNU General Public License for more            #
# details.                                                                  #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program; if not, write to the Free Software               #
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                           #
#############################################################################

include ../common.am
#generate  with "find tests/ -name "*.f90" -exec basename \{\} .f90 \; | xargs" in basedir
#SOURCES=bondi3d bondi2d gauss3d noh3d noh2d KHI3 vortex3d sedov2d riemann3d RTI sedov3d vortex2d pringle riemann1d riemann2d collapse gauss2d ode specialfunctions rngtest
SOURCES=vortex2d sedov2d

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
						 $(top_srcdir)/tap-driver.sh
TESTS = $(SOURCES)
