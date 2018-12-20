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
SOURCES=geometry_base.f90 mesh_base.f90 mesh_generic.f90 mesh_trapezoidal.f90 mesh_midpoint.f90 geometry_cartesian.f90 geometry_polar.f90 geometry_generic.f90
#SOURCES=mesh_trapezoidal.f90 geometry_channel.f90 geometry_cylindrical.f90 geometry_spherical.f90 mesh_generic.f90 geometry_oblatespheroidal.f90 geometry_polar.f90 geometry_cartesian.f90 geometry_sinhpolar.f90 mesh_midpoint.f90 geometry_tanpolar.f90 geometry_generic.f90 geometry_logpolar.f90 geometry_bipolar.f90 geometry_polypolar.f90 geometry_sinhtanh.f90 geometry_lncoshcylindrical.f90 geometry_sinhspherical.f90 geometry_tancylindrical.f90 geometry_elliptic.f90 geometry_sinhcartesian.f90 geometry_bianglespherical.f90
OBJECTS = $(SOURCES:.f90=.o)

noinst_LIBRARIES = libmesh.a
libmesh_a_SOURCES = $(SOURCES)

$(DEPFILE) : $(SOURCES)
	$(AWK) -f $(MAKEDEPS) -v objfiles="$(OBJECTS)" $(addprefix $(srcdir)/,$(SOURCES)) > $@

-include $(DEPFILE)                                                                                                                                                                 
