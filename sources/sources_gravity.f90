!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_gravity.f03                                               #
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
!> \addtogroup sources
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
MODULE sources_gravity_mod
  USE logging_base_mod
  USE mesh_base_mod
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
  CONTAINS
    PROCEDURE :: InitSources_gravity
    PROCEDURE :: InfoSources
    PROCEDURE :: UpdateGravity
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcDiskHeight
    PROCEDURE :: Finalize
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

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VDIM), &
             this%pot(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4), &
             STAT=err)
    IF (err.NE.0) CALL this%Error("InitGravity", "Unable allocate memory!")
    this%accel(:,:,:,:) = 0.
    NULLIFY(this%height,this%h_ext,this%invheight2)

    ! Add source terms to energy equation?
    ! Set this to zero, if a potential is defined in physics_euler2Diamt
    CALL GetAttr(config, "energy", i, 1)
    IF(i.EQ.0) THEN
      this%addtoenergy = .FALSE.
    ELSE
      this%addtoenergy = .TRUE.
    END IF

    ! enable update of disk scale height if requested
    CALL GetAttr(config, "update_disk_height", i, 0)
    IF (i.EQ.1) THEN
      !> \todo check if this is really sufficient, what we really need is
      !! to check whether the geometry is flat or not
      IF (Physics%VDIM.EQ.2) THEN
        this%update_disk_height = .TRUE.
        IF (.NOT.ASSOCIATED(this%height)) &
          ALLOCATE(this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                   STAT=err)
          this%height = 0.0
        IF (err.EQ.0) &
          ALLOCATE(this%h_ext(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                   this%invheight2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                   STAT=err)
          this%h_ext = 0.0
          this%invheight2 = 0.0
        IF (err.NE.0) CALL this%Error("InitGravity","Memory allocation failed!")
      ELSE
         CALL this%Error("InitGravity", "DiskHeight is only supported in 2D")
      END IF
    ELSE
       this%update_disk_height = .FALSE.
    END IF

    ! initialize gravity
    CALL new_gravity(this%glist,Mesh,Fluxes,Physics,config,IO)

    CALL GetAttr(config, "output/height", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       CALL SetAttr(config, "update_disk_height", 1)
       IF (.NOT.ASSOCIATED(this%height)) THEN
          ALLOCATE(this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                   STAT=err)
          IF (err.NE.0) CALL this%Error("SetOutput", &
                                   "Memory allocation failed for this%height!")
       END IF
       CALL SetAttr(IO, "height", &
         this%height(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
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
    CLASS(marray_compound),INTENT(INOUT)  :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    ! update acceleration of all gravity sources
    ! reset gterm
    this%pot(:,:,:,:) = 0.
    this%accel(:,:,:,:) = 0.

    ! go through all gravity terms in the list
    CALL this%UpdateGravity(Mesh,Physics,Fluxes,pvar%data4d,time,dt)

    ! update disk scale height if requested
    IF (this%glist%update_disk_height) THEN
      SELECT TYPE(phys => Physics)
      CLASS IS (physics_eulerisotherm)
        CALL this%CalcDiskHeight(Mesh,phys,pvar%data4d)
      END SELECT
    END IF

    ! gravitational source terms
    CALL Physics%ExternalSources(Mesh,this%accel,pvar,cvar,sterm)

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
    REAL,                INTENT(IN)    :: time, dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER       :: gravptr
    !------------------------------------------------------------------------!

    gravptr => this%glist
    DO WHILE (ASSOCIATED(gravptr))
      ! call specific subroutine
!CDIR IEXPAND
      CALL gravptr%UpdateGravity_single(Mesh,Physics,Fluxes,pvar,time,dt)

      ! add to the sources
      this%accel(:,:,:,:) = this%accel(:,:,:,:) + gravptr%accel(:,:,:,:)
      IF(ASSOCIATED(gravptr%pot)) &
        this%pot(:,:,:,:) = this%pot(:,:,:,:) + gravptr%pot(:,:,:,:)

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
    CLASS(sources_gravity), TARGET, INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)            :: Mesh
    CLASS(physics_eulerisotherm), INTENT(IN)   :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER :: grav_ptr,selfgrav_ptr => null()
    LOGICAL                      :: has_external_potential = .FALSE.
    !------------------------------------------------------------------------!
    ! reset inverse scale height^2
    this%invheight2(:,:,:) = 0.0
    ! go through all gravity terms in the list
    grav_ptr => this%glist
    DO WHILE(ASSOCIATED(grav_ptr))
      SELECT TYPE (grav => grav_ptr)
      CLASS IS (gravity_pointmass)
!CDIR IEXPAND
        CALL grav%CalcDiskHeight_single(Mesh,Physics,pvar,Physics%bccsound%data3d,this%h_ext,this%height)
        this%invheight2(:,:,:) = this%invheight2(:,:,:) + 1./this%h_ext(:,:,:)**2
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
        this%h_ext(:,:,:) = 1./SQRT(this%invheight2(:,:,:))
      END IF
      CALL selfgrav_ptr%CalcDiskHeight_single(Mesh,Physics,pvar,Physics%bccsound%data3d, &
                                  this%h_ext,this%height)
    ELSE
      ! non-selfgravitating disk
      this%height(:,:,:) = 1./SQRT(this%invheight2(:,:,:))
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
    CLASS(sources_gravity), INTENT(INOUT) :: this
    CLASS(gravity_base),    POINTER       :: gravptr
    !------------------------------------------------------------------------!
    DEALLOCATE(this%pot,this%accel)
    IF(ASSOCIATED(this%h_ext)) DEALLOCATE(this%h_ext)
    IF(ASSOCIATED(this%height)) DEALLOCATE(this%height)

    gravptr => this%glist
    DO WHILE (ASSOCIATED(gravptr))

      CALL gravptr%Finalize()

      gravptr => gravptr%next
    END DO

    CALL this%Finalize_base()
    IF(ASSOCIATED(this%next)) CALL this%next%Finalize()
  END SUBROUTINE Finalize

END MODULE sources_gravity_mod
