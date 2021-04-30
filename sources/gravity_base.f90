!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_base.f90                                                  #
!#                                                                           #
!# Copyright (C) 2014-2019                                                   #
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
!> \defgroup gravity Gravity
!! \{
!! \brief Source terms for gravitational forces
!! \}
!----------------------------------------------------------------------------!
!> \addtogroup gravity
!! - general parameters of gravity group as key-values
!! \key{gtype,INTEGER,Type of gravity source}
!! \key{output/accel,INTEGER,enable(=1) output of acceleration,0}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief generic gravity terms module providing functionaly common to all
!! gravity terms
!!
!! \ingroup gravity
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
    CLASS(marray_base), ALLOCATABLE    :: accel        !< acceleration
    REAL, DIMENSION(:,:,:,:), POINTER  :: pot          !< general potential
  CONTAINS
    PROCEDURE :: InitGravity
    PROCEDURE (SetOutput),             DEFERRED :: SetOutput
    PROCEDURE (UpdateGravity_single),  DEFERRED :: UpdateGravity_single
    PROCEDURE (CalcDiskHeight_single), DEFERRED :: CalcDiskHeight_single
    PROCEDURE :: GetGravityPointer
    PROCEDURE :: Finalize_base
  END TYPE gravity_base

  ABSTRACT INTERFACE
    SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
      IMPORT gravity_base, mesh_base, physics_base, dict_typ
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      TYPE(Dict_TYP),      POINTER       :: config,IO
    END SUBROUTINE
    SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
      IMPORT gravity_base, mesh_base, physics_base, fluxes_base, marray_compound
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
      REAL,                INTENT(IN)    :: time,dt
      CLASS(marray_compound), INTENT(INOUT) :: pvar
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

  SUBROUTINE InitGravity(this,Mesh,Physics,gravity_name,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base),  INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CHARACTER(LEN=*),     INTENT(IN)    :: gravity_name
    TYPE(Dict_TYP),POINTER              :: config,IO
    !------------------------------------------------------------------------!
    CHARACTER(LEN=1) :: xyz(3) = (/"x","y","z"/)
    INTEGER :: gtype,k,valwrite
    !------------------------------------------------------------------------!
    ! basic initialization of gravity module
    CALL GetAttr(config, "gtype", gtype)
    ! allocate memory for new gravity term
    CALL this%InitLogging(gtype,gravity_name)
    ALLOCATE(this%accel)
    this%accel = marray_base(Physics%VDIM)
    ! reset acceleration
    this%accel%data1d(:) = 0.0
    ! check whether gravitational acceleration of a particular gravity
    ! module should be written into the data file
    CALL GetAttr(config, "output/accel", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       DO k=1,SIZE(this%accel%data4d,4)
          CALL SetAttr(IO, ("accel_" // xyz(k)),&
             this%accel%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,k))
       END DO
    END IF
    CALL this%Info(" GRAVITY--> gravity term:      " // this%GetName())
  END SUBROUTINE InitGravity

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
    CLASS(gravity_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (ALLOCATED(this%accel)) DEALLOCATE(this%accel)
  END SUBROUTINE Finalize_base

END MODULE gravity_base_mod
