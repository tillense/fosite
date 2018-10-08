!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: boundary_nogradients.f90                                          #
!#                                                                           #
!# Copyright (C) 2006-2017                                                   #
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
!! \author Jannes Klee
!!
!! \brief Boundary module for reflecting boundaries
!!
!! \extends boundary_nogradients
!! \ingroup boundary
!----------------------------------------------------------------------------!
MODULE boundary_nogradients_mod
  USE boundary_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS(boundary_base) :: boundary_nogradients
  CONTAINS
    PROCEDURE :: InitBoundary_nogradients
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE boundary_nogradients
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "nogradients"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       boundary_nogradients
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for nogradients boundary conditions
  SUBROUTINE InitBoundary_nogradients(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_nogradients), INTENT(INOUT) :: this
    CLASS(physics_base), INTENT(IN)            :: Physics
    CLASS(mesh_base), INTENT(IN)               :: Mesh
    TYPE(Dict_TYP), POINTER                    :: config
    INTEGER                                    :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)                                 :: dir
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,NO_GRADIENTS,boundcond_name,dir,config)
  END SUBROUTINE InitBoundary_nogradients

  !> \public Applies the nogradients boundary condition
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_nogradients), INTENT(INOUT) :: this
    CLASS(mesh_base), INTENT(IN)            :: Mesh
    CLASS(physics_base), INTENT(IN)         :: Physics
    REAL,                 INTENT(IN)        :: time
    REAL :: pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                 Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM)
    !------------------------------------------------------------------------!
    INTEGER                                 :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(INOUT)                           :: pvar
    !------------------------------------------------------------------------!
    SELECT CASE(this%direction%GetType())
    CASE(WEST)
       ! UNROLL=Mesh%GNUM would be sufficient, but the compiler does
       ! not know the value of Mesh%GNUM, hence we set UNROLL=4 and
       ! hope that nobody sets Mesh%GNUM to a value greater than 4
!NEC$ UNROLL(4)
       DO i=1,Mesh%GINUM
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) = &
            pvar(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
       END DO
    CASE(EAST)
!NEC$ UNROLL(4)
       DO i=1,Mesh%GINUM
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:) = &
            pvar(Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
       END DO
    CASE(SOUTH)
!NEC$ UNROLL(4)
       DO j=1,Mesh%GJNUM
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,:) = &
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,Mesh%KMIN:Mesh%KMAX,:)
       END DO
    CASE(NORTH)
!NEC$ UNROLL(4)
       DO j=1,Mesh%GJNUM
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,:) = &
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:)
       END DO
    CASE(BOTTOM)
!NEC$ UNROLL(4)
       DO k=1,Mesh%GKNUM
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-k,:) = &
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN,:)
       END DO
    CASE(TOP)
!NEC$ UNROLL(4)
       DO k=1,Mesh%GKNUM
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX+k,:) = &
            pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMAX,:)
       END DO
    END SELECT

  END SUBROUTINE SetBoundaryData

  !> \public Destructor for nogradients boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_nogradients), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE boundary_nogradients_mod
