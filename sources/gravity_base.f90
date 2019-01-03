!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_base.f90                                                  #
!#                                                                           #
!# Copyright (C) 2014-2018                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup gravity
!! - general parameters of gravity group as key-values
!! \key{gtype,INTEGER,Type of gravity source}
!! \key{energy,INTEGER,Add source terms to energy equation?}
!! \key{output/accel,INTEGER,enable(=1) output of acceleration}
!! \key{output/height,INTEGER,enable(=1) output of disc height}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief generic gravity terms module providing functionaly common to all
!! gravity terms
!----------------------------------------------------------------------------!
MODULE gravity_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, ABSTRACT, EXTENDS(logging_base) :: gravity_base
    !> \name Variables
    CLASS(logging_base), ALLOCATABLE   :: gravitytype  !< type of gravity term
    CLASS(gravity_base), POINTER       :: next => null() !< next gravity in list
    REAL, DIMENSION(:,:,:,:), POINTER  :: accel        !< acceleration
    REAL, DIMENSION(:,:,:,:), POINTER  :: pot          !< general potential
  CONTAINS
    PROCEDURE :: InitGravity
    PROCEDURE :: SetOutput
    PROCEDURE (InfoGravity),     DEFERRED :: InfoGravity
    PROCEDURE (UpdateGravity_single), DEFERRED :: UpdateGravity_single
    PROCEDURE (CalcDiskHeight_single), DEFERRED :: CalcDiskHeight_single
    PROCEDURE :: GetGravityPointer
    PROCEDURE :: Finalize_base
    PROCEDURE (Finalize), DEFERRED :: Finalize
  END TYPE gravity_base

  ABSTRACT INTERFACE
    SUBROUTINE InfoGravity(this,Mesh)
      IMPORT gravity_base, mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(IN)    :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
    END SUBROUTINE
    SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
      IMPORT gravity_base, mesh_base, physics_base, fluxes_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                           INTENT(IN)    :: pvar
      REAL,                INTENT(IN)    :: time,dt
    END SUBROUTINE
    SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
      IMPORT gravity_base, mesh_base, physics_base, marray_compound, marray_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      CLASS(marray_compound), INTENT(INOUT) :: pvar
      TYPE(marray_base),      INTENT(INOUT) :: bccsound,h_ext,height
    END SUBROUTINE
    SUBROUTINE Finalize(this)
      IMPORT gravity_base
      IMPLICIT NONE
      CLASS(gravity_base), INTENT(INOUT) :: this
    END SUBROUTINE
  END INTERFACE
  ! flags for source terms
  INTEGER, PARAMETER :: POINTMASS        = 1
  INTEGER, PARAMETER :: POINTMASS_BINARY = 2
!  INTEGER, PARAMETER :: MONOPOL          = 3
!  INTEGER, PARAMETER :: MULTIGRID        = 4
  INTEGER, PARAMETER :: SPECTRAL         = 5
!  INTEGER, PARAMETER :: POTENTIAL        = 6
  INTEGER, PARAMETER :: SBOXSPECTRAL     = 7
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       gravity_base, &
       ! constants
       POINTMASS, POINTMASS_BINARY, &
       SPECTRAL, SBOXSPECTRAL
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGravity(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base),  INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    TYPE(Dict_TYP),POINTER              :: config,IO
    !------------------------------------------------------------------------!
    CALL this%Info(" GRAVITY--> gravity term:      " // this%GetName())
    CALL this%InfoGravity(Mesh)
  END SUBROUTINE InitGravity

  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    TYPE(Dict_TYP),      POINTER       :: config,IO
    !------------------------------------------------------------------------!
    CHARACTER(LEN=1) :: xyz(3) = (/"x","y","z"/)
    INTEGER          :: valwrite,k
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "output/accel", valwrite, 1)
    IF (valwrite .EQ. 1) THEN
       DO k=1,SIZE(this%accel,4)
          CALL SetAttr(IO, ("accel_" // xyz(k)),&
             this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,k))
       END DO
    END IF

    CALL GetAttr(config, "output/potential", valwrite, 1)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "potential",&
         this%pot(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1))

  END SUBROUTINE SetOutput

  FUNCTION GetGravityPointer(list,stype) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), TARGET  :: list
    INTEGER, INTENT(IN)          :: stype
    CLASS(gravity_base), POINTER :: gp
    !------------------------------------------------------------------------!
    gp => list
    DO
       IF (ASSOCIATED(gp).EQV..FALSE.) EXIT
!CDIR IEXPAND
       IF (gp%GetType().EQ.stype) RETURN
       gp => gp%next
    END DO
  END FUNCTION GetGravityPointer

  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base) :: this
    !------------------------------------------------------------------------!
    ! nothing intializaed

  END SUBROUTINE Finalize_base

END MODULE gravity_base_mod
