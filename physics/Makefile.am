#                                               -*- Automake -*-
# Process this file with automake --add-missing to produce a Makefile.
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# configuration file: Makefile.am                                           #
#                                                                           #
# Copyright (C) 2012-2013                                                   #
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
#SOURCES=physics_euler3Drotamt.f90 physics_generic.f90 constants_generic.f90 physics_euler2Diamt.f90 constants_SI.f90 physics_euler2Disoiamt.f90 physics_euler2Disothm.f90 physics_euler3Drotsym.f90 physics_euler2D.f90 constants_geometrical.f90 physics_euler3Drotsymsgs.f90 physics_euler2Dsgs.f90 physics_euler3Drotamtsgs.f90 physics_euler2Diamrot.f90 physics_euler2Disoiamrot.f90
SOURCES=physics_base.f90 constants_base.f90 constants_SI.f90 physics_euler2Disothm.f90 constants_geometrical.f90 physics_generic.f90 constants_generic.f90
OBJECTS = $(SOURCES:.f90=.o)

noinst_LIBRARIES = libphysics.a
libphysics_a_SOURCES = $(SOURCES)

$(DEPFILE) : $(SOURCES)
	$(AWK) -f $(MAKEDEPS) -v objfiles="$(OBJECTS)" $(addprefix $(srcdir)/,$(SOURCES)) > $@

-include $(DEPFILE)                                                                                                                                                                 
