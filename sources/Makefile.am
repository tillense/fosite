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
#SOURCES=sources_cooling.f90 gravity_monopol.f90 gravity_binary.f90 multipole_cylindrical.f90 sources_generic.f90 gravity_multigrid.f90 sources_diskthomson.f90 multipole_spherical.f90 gravity_pointmass.f90 sources_rotframe.f90 sources_c_accel.f90 sources_sgs.f90 sources_diskcooling.f90 sources_viscosity.f90 gravity_spectral.f90 multipole_generic.f90 sources_wave_damping.f90 sources_forcing.f90 gravity_potential.f90 gravity_generic.f90 sources_stellarheating.f90 sources_planetheating.f90 sources_planetcooling.f90 sources_shearbox.f90 gravity_sboxspectral.f90
SOURCES=sources_base.f90 sources_generic.f90
OBJECTS = $(SOURCES:.f90=.o)

noinst_LIBRARIES = libsources.a
libsources_a_SOURCES = $(SOURCES)

$(DEPFILE) : $(SOURCES)
	$(AWK) -f $(MAKEDEPS) -v objfiles="$(OBJECTS)" $(addprefix $(srcdir)/,$(SOURCES)) > $@

-include $(DEPFILE)                                                                                                                                                                 
