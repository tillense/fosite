!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_gravity.f90                                               #
!#                                                                           #
!# Copyright (C) 2014-2024                                                   #
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
!> \addtogroup sources
!! - general parameters of gravity group as key-values
!! \key{energy,INTEGER,Add source terms to energy equation?}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief generic gravity terms module providing functionaly common to all
!! gravity terms
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_gravity_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE physics_base_mod
  USE sources_base_mod
  USE sources_c_accel_mod
  USE gravity_base_mod
  USE gravity_generic_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: source_name = "gravity"
  TYPE, EXTENDS(sources_c_accel) :: sources_gravity
    CLASS(gravity_base),    POINTER :: glist => null() !< list of gravity terms
     !> 0: no src term in energy equation
     !! 1: src term in energy equation
     LOGICAL                        :: addtoenergy
  CONTAINS
    PROCEDURE :: InitSources_gravity
    PROCEDURE :: InfoSources
    PROCEDURE :: UpdateGravity
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcDiskHeight
    FINAL :: Finalize
  END TYPE sources_gravity
  ABSTRACT INTERFACE
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       sources_gravity
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_gravity(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_gravity), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(fluxes_base),     INTENT(IN)    :: Fluxes
    CLASS(physics_base),    INTENT(IN)    :: Physics
    TYPE(Dict_TYP),         POINTER       :: config, IO
    !------------------------------------------------------------------------!
    INTEGER :: err, stype,i, valwrite
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"stype",stype)
    CALL this%InitLogging(stype,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    ALLOCATE(this%accel)
    this%accel = marray_base(Physics%VDIM)
    this%accel%data1d(:) = 0.
    this%pot = marray_base(4)
    this%pot%data1d(:) = 0.

    ! Add source terms to energy equation?
    ! Set this to zero, if a potential is defined in physics_euler2Diamt
    CALL GetAttr(config, "energy", i, 1)
    IF(i.EQ.0) THEN
      this%addtoenergy = .FALSE.
    ELSE
      this%addtoenergy = .TRUE.
    END IF

    ! initialize gravity
    CALL new_gravity(this%glist,Mesh,Fluxes,Physics,config,IO)

    ! activate output of scale height if requested (default: disabled)
    CALL GetAttr(config, "output/height", valwrite, 0)

    ! this parameter is usually set automatically depending on other
    ! settings (see sources_generic), but if output of height
    ! is enabled, computation of disk height is activated as well
    ! (works only 1D/2D accretion disk simulations with flat geometry)
    IF (valwrite.EQ.1) CALL SetAttr(config, "update_disk_height", 1)
    CALL GetAttr(config, "update_disk_height",valwrite, 0)

    IF (valwrite.EQ.1) THEN
      !> \todo check if this is sufficient, what we really need is
      !! to check whether the geometry is flat or not
      IF (Physics%VDIM.EQ.2) THEN
        this%update_disk_height = .TRUE.
        this%height = marray_base()
        this%h_ext= marray_base()
        this%invheight2 = marray_base()
        this%height%data1d = 0.0
        this%h_ext%data1d = 0.0
        this%invheight2%data1d = 0.0

        IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, "height", &
             this%height%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
        END IF
      ELSE
        CALL this%Error("InitGravity", "Computation of disk scale height is only supported in 1D/2D flat geometries")
      END IF
    ELSE
       this%update_disk_height = .FALSE.
    END IF

    ! activate output of gravitational acceleration
    CALL GetAttr(config, "output/accel", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       CALL SetAttr(IO, "accel", &
          this%accel%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1:Physics%VDIM))
    END IF

  END SUBROUTINE InitSources_gravity

  !> Evaluates source-terms by gravitation
  !!
  !! The gravitational source term evaluates all forces that are produced by
  !! gravitational participants.
  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE physics_euler_mod, ONLY : statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_gravity), INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(INOUT) :: Physics
    CLASS(fluxes_base),     INTENT(IN)    :: Fluxes
    CLASS(sources_base),    INTENT(INOUT) :: Sources
    REAL,                   INTENT(IN)    :: time, dt
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    ! go through all gravity terms in the list
    CALL this%UpdateGravity(Mesh,Physics,Fluxes,pvar,time,dt)

    ! update disk scale height if requested
    IF (this%update_disk_height) THEN
      SELECT TYPE(phys => Physics)
      CLASS IS (physics_eulerisotherm)
        CALL this%CalcDiskHeight(Mesh,phys,pvar)
      END SELECT
    END IF

    ! gravitational source terms
    CALL Physics%ExternalSources(this%accel,pvar,cvar,sterm)

    !> \todo The treatment of energy sources should be handled in the physics
    !! module and not here!
    ! Set src term in energy equation to zero, if it is handeled in the physics
    ! module
    SELECT TYPE(s => sterm)
    TYPE IS (statevector_euler)
      IF (.NOT.this%addtoenergy) s%energy%data1d(:) = 0.
    END SELECT

  END SUBROUTINE ExternalSources_single

  !> Updates gravity of all gravity source modules
  SUBROUTINE UpdateGravity(this,Mesh,Physics,Fluxes,pvar,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_gravity), TARGET, INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    REAL,                INTENT(IN)    :: time, dt
    !------------------------------------------------------------------------!
    INTEGER                            :: i,j,k,l
    CLASS(gravity_base), POINTER       :: gravptr
    !------------------------------------------------------------------------!
    ! update acceleration of all gravity sources
    ! reset gterm
    this%pot%data1d(:) = 0.
    this%accel%data1d(:) = 0.

    gravptr => this%glist
    DO WHILE (ASSOCIATED(gravptr))
      ! get gravitational acceleration of each gravity module
      CALL gravptr%UpdateGravity_single(Mesh,Physics,Fluxes,pvar,time,dt)

      ! add contribution to the overall gravitational acceleration
      this%accel = this%accel + gravptr%accel

      ! add potential if available
      IF (ALLOCATED(gravptr%pot)) THEN
        this%pot = this%pot + gravptr%pot
      END IF

      ! next source term
      gravptr => gravptr%next
    END DO
  END SUBROUTINE

  SUBROUTINE CalcDiskHeight(this,Mesh,Physics,pvar)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE gravity_pointmass_mod, ONLY : gravity_pointmass
    USE gravity_spectral_mod, ONLY : gravity_spectral
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_gravity),TARGET,INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_eulerisotherm), INTENT(INOUT) :: Physics
    CLASS(marray_compound),       INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER :: grav_ptr,selfgrav_ptr => null()
    LOGICAL                      :: has_external_potential = .FALSE.
    !------------------------------------------------------------------------!
    ! reset inverse scale height^2
    this%invheight2%data1d(:) = 0.0
    ! go through all gravity terms in the list
    grav_ptr => this%glist
    DO WHILE(ASSOCIATED(grav_ptr))
      SELECT TYPE (grav => grav_ptr)
      CLASS IS (gravity_pointmass)
        CALL grav%CalcDiskHeight_single(Mesh,Physics,pvar,Physics%bccsound,this%h_ext,this%height)
        this%invheight2%data1d(:) = this%invheight2%data1d(:) + 1./this%height%data1d(:)**2
        has_external_potential = .TRUE.
      CLASS IS (gravity_spectral)
        selfgrav_ptr => grav_ptr
      END SELECT

      ! next gravity term
      grav_ptr => grav_ptr%next
    END DO

    ! self-gravity of the disk needs special treatment
    IF (ASSOCIATED(selfgrav_ptr)) THEN
      IF (has_external_potential) THEN
        ! compute the resultant height due to all external gravitational forces
        this%h_ext%data1d(:) = 1./SQRT(this%invheight2%data1d(:))
      END IF
      CALL selfgrav_ptr%CalcDiskHeight_single(Mesh,Physics,pvar,Physics%bccsound,this%h_ext,this%height)
    ELSE
      ! non-selfgravitating disk
      this%height%data1d(:) = 1./SQRT(this%invheight2%data1d(:))
    END IF

  END SUBROUTINE CalcDiskHeight



  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_gravity), INTENT(IN) :: this
    CLASS(mesh_base),       INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
  END SUBROUTINE InfoSources

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(sources_gravity),  INTENT(INOUT) :: this
    CLASS(gravity_base),    POINTER       :: gravptr,tmpptr
    !------------------------------------------------------------------------!
    gravptr => this%glist
    DO WHILE (ASSOCIATED(gravptr))
      tmpptr => gravptr
      gravptr => gravptr%next
      DEALLOCATE(tmpptr)
    END DO
    ! deallocation of this%accel is done in inherited destructor
    ! which is called automatically
  END SUBROUTINE Finalize

END MODULE sources_gravity_mod
