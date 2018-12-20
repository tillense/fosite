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
#SOURCES=boundary_axis.f90 boundary_nogradients.f90 boundary_farfield.f90 boundary_custom.f90 boundary_noh.f90 boundary_fixed.f90 boundary_noslip.f90 boundary_extrapolation.f90 boundary_reflecting.f90 boundary_absorbing.f90 boundary_folded.f90 boundary_periodic.f90 boundary_generic.f90 boundary_dmr.f90 boundary_shearing.f90
SOURCES=boundary_base.f90 boundary_reflecting.f90 boundary_generic.f90 boundary_periodic.f90
OBJECTS = $(SOURCES:.f90=.o)

noinst_LIBRARIES = libboundary.a
libboundary_a_SOURCES = $(SOURCES)

$(DEPFILE) : $(SOURCES)
	$(AWK) -f $(MAKEDEPS) -v objfiles="$(OBJECTS)" $(addprefix $(srcdir)/,$(SOURCES)) > $@

-include $(DEPFILE)                                                                                                                                                                 
