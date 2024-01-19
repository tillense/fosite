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

  !> constructor for FileIO class
  FUNCTION new_fileio(Mesh,Physics,Timedisc,Sources,config,IO) RESULT(new_fio)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),     INTENT(IN)          :: Mesh
    CLASS(physics_base),  INTENT(IN)          :: Physics
    CLASS(timedisc_base), INTENT(IN)          :: Timedisc
    CLASS(sources_base),  INTENT(IN), POINTER :: Sources
    TYPE(DICT_TYP),       INTENT(IN), POINTER :: config
    TYPE(DICT_TYP),       INTENT(IN), POINTER ::  IO
    !------------------------------------------------------------------------!
    CLASS(fileio_base), ALLOCATABLE           :: new_fio
    INTEGER                                   :: fileformat
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"fileformat",fileformat)

    ! allocate data
    SELECT CASE(fileformat)
    CASE(GNUPLOT)
      ALLOCATE(fileio_gnuplot::new_fio)
!     CASE(VTK)
!       ALLOCATE(fileio_vtk::new_fio)
!     CASE(BINARY)
!       ALLOCATE(fileio_binary::new_fio)
!     CASE(XDMF)
!       ALLOCATE(fileio_xdmf::new_fio)
    CASE DEFAULT
      CALL Mesh%Error("fileio_base::CreateFileIO","Unknown file format.")
    END SELECT

    ! call initialization
    CALL new_fio%InitFileIO(Mesh,Physics,Timedisc,Sources,config,IO)
  END FUNCTION new_fileio

END MODULE fileio_generic_mod
