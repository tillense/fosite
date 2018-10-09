!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_shearbox.f90                                              #
!#                                                                           #
!# Copyright (C) 2010-2018                                                   #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup sources
!!
!! \brief source term module for shearing forces in a shearingsheet
!!
!! \todo This module can only be used with fixed 2D plane in x-y direction.
!!       It should be no problem to make it also run in 3D.
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!!
!! \brief Source terms module for fictious forces in a shearingsheet.
!!
!! \extends sources_c_accel
!! \ingroup sources
!----------------------------------------------------------------------------!
MODULE sources_shearbox_mod
  USE sources_c_accel_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
    CHARACTER(LEN=32) :: source_name = "forces in shearing-box"
  !--------------------------------------------------------------------------!

  TYPE, EXTENDS(sources_c_accel) :: sources_shearbox
 CONTAINS
    PROCEDURE :: InitSources_shearbox
    PROCEDURE :: InfoSources
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: Finalize
  END TYPE

  PUBLIC :: &
       ! classes
       sources_shearbox

CONTAINS

  !> \public Constructor of sources shearbox module.
  !!
  !! This subroutine reads the necessary config data for the source module
  !! for a shearingsheet. It initializes the sources type and the
  !! accelerations array.
  SUBROUTINE InitSources_shearbox(this,Mesh,Physics,Fluxes,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox)          :: this
    CLASS(mesh_base),     INTENT(IN) :: Mesh
    CLASS(physics_base),  INTENT(IN) :: Physics
    CLASS(fluxes_base),   INTENT(IN) :: Fluxes
    TYPE(Dict_TYP),          POINTER :: config, IO
    !------------------------------------------------------------------------!
    INTEGER                          :: err, valwrite, stype
    !------------------------------------------------------------------------!
    CALL GetAttr(config,"stype",stype)
    CALL this%InitLogging(stype,source_name)
    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

    SELECT CASE(Physics%GetType())
    CASE(EULER2D,EULER2D_ISOTHERM)
      ! do nothing
    CASE DEFAULT
      CALL this%Error("InitSources_shearbox","physics not supported")
    END SELECT

    SELECT CASE(Mesh%Geometry%GetType())
    CASE(CARTESIAN)
      ! do nothing
    CASE DEFAULT
      CALL this%Error("InitSources_shearbox","mesh not supported")
    END SELECT

    ALLOCATE(&
       this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                  Mesh%KGMIN:Mesh%KGMAX,Physics%DIM), &
       STAT = err)

    IF (err.NE.0) CALL this%Error("InitSources_shearbox", "Unable allocate memory!")

    ! reset acceleration term
    this%accel(:,:,:,:) = 0.0
  END SUBROUTINE InitSources_shearbox

  !> \public Write shearbox fictious forces parameters to screen.
  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox), INTENT(IN) :: this
    CLASS(mesh_base),        INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)       :: omega_str,q_str
    !------------------------------------------------------------------------!
    WRITE (omega_str,'(ES8.2)') Mesh%OMEGA
    WRITE (q_str,'(ES8.2)')     Mesh%Q
    CALL this%Info("            angular velocity:  " // TRIM(omega_str))
    CALL this%Info("            shearing parameter:" // TRIM(q_str))
  END SUBROUTINE InfoSources

  !> \public Compute fictious forces within a shearingsheet.
  !!
  !! Computation is done via
  !! \f[
  !!    \mathbf{a} = - 2 \mathbf{\Omega} \times \mathbf{v} + 2 q \Omega^2 x
  !!                    \mathbf{\hat{e}_x}.
  !! \f]
  !! See for example \cite gammie2001 or \cite hawley1995 .
  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox), INTENT(INOUT) :: this
    CLASS(mesh_base),        INTENT(IN)    :: Mesh
    CLASS(physics_base),     INTENT(INOUT) :: Physics
    CLASS(fluxes_base),      INTENT(IN)    :: Fluxes
    REAL,                    INTENT(IN)    :: time, dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                             INTENT(IN)    :: cvar,pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                             INTENT(OUT)   :: sterm
    !------------------------------------------------------------------------!
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------!
    IF (Mesh%FARGO.EQ.0) THEN
      sterm(:,:,:,:) = 0.0
!NEC$ IVDEP
      FORALL(i=Mesh%IMIN:Mesh%IMAX,j=Mesh%JMIN:Mesh%JMAX,k=Mesh%KMIN:Mesh%KMAX)
            this%accel(i,j,k,1) = Mesh%OMEGA*2.0*(Mesh%Q*Mesh%OMEGA* &
              Mesh%bcenter(i,j,k,1) + pvar(i,j,k,Physics%YVELOCITY))
            this%accel(i,j,k,2) = -Mesh%OMEGA*2.0*pvar(i,j,k,Physics%XVELOCITY)
      END FORALL
      ! shearingsheet inertial forces source terms
      CALL Physics%ExternalSources(Mesh,this%accel,pvar,cvar,sterm)
    ELSE IF (Mesh%FARGO.EQ.3) THEN
      sterm(:,:,:,Physics%DENSITY) = 0.0
!NEC$ IVDEP
      FORALL(i=Mesh%IMIN:Mesh%IMAX,j=Mesh%JMIN:Mesh%JMAX,k=Mesh%KMIN:Mesh%KMAX)
            sterm(i,j,k,Physics%XMOMENTUM) = &
                pvar(i,j,k,Physics%DENSITY)*Mesh%OMEGA*2.0*pvar(i,j,k,Physics%YVELOCITY)
            sterm(i,j,k,Physics%YMOMENTUM) = &
                pvar(i,j,k,Physics%DENSITY)*Mesh%OMEGA*(Mesh%Q-2.0)*pvar(i,j,k,Physics%XVELOCITY)
      END FORALL
      IF (Physics%PRESSURE .GT. 0) THEN
        FORALL(i=Mesh%IMIN:Mesh%IMAX,j=Mesh%JMIN:Mesh%JMAX,k=Mesh%KMIN:Mesh%KMAX)
            sterm(i,j,k,Physics%ENERGY) = &
                 pvar(i,j,k,Physics%DENSITY)*Mesh%Q*Mesh%OMEGA* &
                 pvar(i,j,k,Physics%XVELOCITY)*pvar(i,j,k,Physics%YVELOCITY)
        END FORALL
      END IF
    END IF
  END SUBROUTINE ExternalSources_single

  !> \public Closes the shearingsheet source term.
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_shearbox), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%sources_c_accel%Finalize()
    CALL this%next%Finalize()
  END SUBROUTINE Finalize

END MODULE sources_shearbox_mod
