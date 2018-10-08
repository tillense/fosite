!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: reconstruction_constant.f90                                       #
!#                                                                           #
!# Copyright (C) 2007-2017                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> \author Tobias Illenseer
!> \author Jannes Klee
!!
!! \brief basic module for constant (zero order) reconstruction
!!
!! \extends reconstruction_common
!! \ingroup reconstruction
!----------------------------------------------------------------------------!
MODULE reconstruction_constant_mod
  USE logging_base_mod
  USE reconstruction_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS (reconstruction_base) :: reconstruction_constant
  ! no data declarations
  CONTAINS
    PROCEDURE             :: InitReconstruction_constant
    PROCEDURE             :: CalculateStates
    PROCEDURE             :: Finalize
  END TYPE reconstruction_constant
  !--------------------------------------------------------------------------!
  CHARACTER(LEN=32), PARAMETER    :: recontype_name = "constant"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       reconstruction_constant
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitReconstruction_constant(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_constant), INTENT(INOUT) :: this
    CLASS(mesh_base),               INTENT(IN)    :: Mesh
    CLASS(physics_base),            INTENT(IN)    :: Physics
    TYPE(Dict_TYP),POINTER                        :: config,IO
    !-----------------------------------------------------------------------!
    CALL this%InitReconstruction(Mesh,Physics,config,IO,CONSTANT,recontype_name)
  END SUBROUTINE InitReconstruction_constant


  PURE SUBROUTINE CalculateStates(this,Mesh,Physics,npos,dx,dy,dz,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_constant), INTENT(INOUT) :: this
    CLASS(mesh_base),               INTENT(IN)    :: Mesh
    CLASS(physics_base),            INTENT(IN)    :: Physics
    INTEGER                                       :: npos
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,npos)   :: dx,dy,dz
    REAL :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,    &
                 Mesh%KGMIN:Mesh%KGMAX,Physics%vnum)
    REAL :: rstates(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,npos,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER                                       :: n
    !------------------------------------------------------------------------!
    INTENT(IN)                                    :: npos,rvar,dx,dy,dz
    INTENT(OUT)                                   :: rstates
    !------------------------------------------------------------------------!

    ! reconstruct boundary states
    FORALL (n=1:npos)
       rstates(:,:,:,n,:) = rvar(:,:,:,:)
    END FORALL
  END SUBROUTINE CalculateStates


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_constant), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE reconstruction_constant_mod
