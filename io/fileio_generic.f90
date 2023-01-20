!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fileio_generic.f90                                                #
!#                                                                           #
!# Copyright (C) 2016 Manuel Jung <mjung@astrophysik.uni-kiel.de>            #
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
!! \author Jannes Klee
!!
!! \brief constructor for fileio class
!!
!! This module allocates the mesh class and decides which specific
!! mesh to use from the config.
!----------------------------------------------------------------------------!
MODULE fileio_generic_mod
  USE fileio_base_mod
  USE fileio_gnuplot_mod
!   USE fileio_vtk_mod
!   USE fileio_binary_mod
!   USE fileio_xdmf_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE timedisc_base_mod
  USE sources_base_mod
  USE common_dict

CONTAINS

  SUBROUTINE new_fileio(Fileio,Mesh,Physics,Timedisc,Sources,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fileio_base), ALLOCATABLE           :: Fileio
    CLASS(mesh_base),     INTENT(IN)          :: Mesh
    CLASS(physics_base),  INTENT(IN)          :: Physics
    CLASS(timedisc_base), INTENT(IN)          :: Timedisc
    CLASS(sources_base),  INTENT(IN), POINTER :: Sources
    TYPE(DICT_TYP),       INTENT(IN), POINTER :: config
    TYPE(DICT_TYP),       INTENT(IN), POINTER ::  IO
    !------------------------------------------------------------------------!
    INTEGER                                   :: fileformat
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"fileformat",fileformat)

    ! allocate data
    SELECT CASE(fileformat)
    CASE(GNUPLOT)
      ALLOCATE(fileio_gnuplot::Fileio)
!     CASE(VTK)
!       ALLOCATE(fileio_vtk::Fileio)
!     CASE(BINARY)
!       ALLOCATE(fileio_binary::Fileio)
!     CASE(XDMF)
!       ALLOCATE(fileio_xdmf::Fileio)
    CASE DEFAULT
      CALL Fileio%Error("fileio_generic::new_fileio","Unknown file format.")
    END SELECT

    ! call initialization
    SELECT TYPE(obj => Fileio)
    TYPE IS (fileio_gnuplot)
      CALL obj%InitFileio_gnuplot(Mesh,Physics,Timedisc,Sources,config,IO)
!     TYPE IS (fileio_vtk)
!       CALL obj%InitFileio_vtk(Mesh,Physics,Timedisc,Sources,config,IO)
!     TYPE IS (fileio_binary)
!       CALL obj%InitFileio_binary(Mesh,Physics,Timedisc,Sources,config,IO)
!     TYPE IS (fileio_xdmf)
!       CALL obj%InitFileio_xdmf(Mesh,Physics,Timedisc,Sources,config,IO)
    END SELECT
  END SUBROUTINE
END MODULE fileio_generic_mod
