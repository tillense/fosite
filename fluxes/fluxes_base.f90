!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fluxes_base.f90                                                   #
!#                                                                           #
!# Copyright (C) 2007-2017                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief base module for numerical flux functions
!!
!! \todo implement constructor and other fluxes, e.g. HLLC
!!
!! \ingroup fluxes
!----------------------------------------------------------------------------!
MODULE fluxes_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE reconstruction_generic_mod
  USE physics_base_mod
  USE physics_generic_mod
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
  TYPE, ABSTRACT, EXTENDS (logging_base) ::  fluxes_base
     !> \name Classes
     CLASS(reconstruction_base), ALLOCATABLE &
                                  :: Reconstruction  !< reconstruction method
     CLASS(marray_base), ALLOCATABLE &
                                  :: minwav,maxwav   !< min/max wave speeds
     CLASS(marray_compound), POINTER :: prim, &      !< primitive/conservative
                                        cons         !< state vectors on cell faces
     !> \name
     !! #### various data fields
     REAL, DIMENSION(:,:,:,:), POINTER &
                                  :: bxflux,byflux, &
                                     bzflux,bxfold, &
                                     byfold,bzfold   !< boundary fluxes
     REAL, DIMENSION(:,:,:,:,:), POINTER &
                                  :: pfluxes         !< physical fluxes
  CONTAINS
    PROCEDURE                               :: InitFluxes
    PROCEDURE (CalculateFluxes), DEFERRED   :: CalculateFluxes
    PROCEDURE                               :: CalculateFaceData
    PROCEDURE                               :: GetBoundaryFlux
    PROCEDURE (Finalize), DEFERRED          :: Finalize
    PROCEDURE                               :: Finalize_base
  END TYPE fluxes_base

  ABSTRACT INTERFACE
    SUBROUTINE CalculateFluxes(this,Mesh,Physics,pvar,cvar, &
                xfluxdydz,yfluxdzdx,zfluxdxdy)
      IMPORT fluxes_base,mesh_base,physics_base,marray_compound
      IMPLICIT NONE
      CLASS(fluxes_base),   INTENT(INOUT) :: this
      CLASS(mesh_base),     INTENT(IN)    :: Mesh
      CLASS(physics_base),  INTENT(INOUT) :: Physics
      CLASS(marray_compound),INTENT(INOUT):: pvar,cvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                            INTENT(OUT)   :: xfluxdydz,yfluxdzdx,zfluxdxdy
    END SUBROUTINE
    SUBROUTINE Finalize(this)
      IMPORT fluxes_base
      IMPLICIT NONE
      CLASS(fluxes_base), INTENT(INOUT)   :: this
    END SUBROUTINE
  END INTERFACE
  INTEGER, PARAMETER :: KT       = 1
  INTEGER, PARAMETER :: HLL      = 2
  INTEGER, PARAMETER :: HLLC     = 3
  INTEGER, PARAMETER :: EXACT    = 4
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fluxes_base, &
       ! constants
       KT, HLL, HLLC, EXACT
  !--------------------------------------------------------------------------!

CONTAINS

  !> Initialize Fluxes
  SUBROUTINE InitFluxes(this,Mesh,Physics,config,IO,ftype,fname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_base),  INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    TYPE(Dict_TYP),POINTER             :: config,IO
    INTEGER                            :: err
    INTEGER                            :: ftype
    CHARACTER(LEN=*)                   :: fname
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER             :: IOrec => null()
    INTEGER                            :: valwrite,i
    CHARACTER(LEN=60)                  :: key
    !------------------------------------------------------------------------!
    CALL this%InitLogging(ftype,fname)

    ! check initialization of Mesh and Physics
    IF (.NOT.Mesh%Initialized().OR..NOT.Physics%Initialized()) &
         CALL this%Error("InitFluxes","mesh and/or physics module uninitialized")

    ! allocate memory for all arrays used in fluxes
    ALLOCATE(this%minwav,this%maxwav, &
      this%pfluxes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Mesh%NFACES,Physics%VNUM), &
      this%bxflux(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,2,Physics%VNUM),      &
      this%byflux(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,2,Physics%VNUM),      &
      this%bzflux(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2,Physics%VNUM),      &
      this%bxfold(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,2,Physics%VNUM),      &
      this%byfold(Mesh%KGMIN:Mesh%KGMAX,Mesh%IGMIN:Mesh%IGMAX,2,Physics%VNUM),      &
      this%bzfold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2,Physics%VNUM),      &
      STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("InitFluxes", "Unable to allocate memory.")
    END IF

    ! create RANK 1 mesh arrays
    this%minwav = marray_base(Mesh%NDIMS)
    this%maxwav = marray_base(Mesh%NDIMS)

    ! initialize state vectors on cell interfaces
    CALL new_statevector(Physics,this%prim,PRIMITIVE,Mesh%NFACES)
    CALL new_statevector(Physics,this%cons,CONSERVATIVE,Mesh%NFACES)

    ! print some information
    CALL this%Info(" FLUXES---> fluxes type        " // TRIM(this%GetName()))

    ! enable output of reconstructed cell face data if requested
    CALL GetAttr(config, "output/pstates", valwrite,0)
    IF(valwrite.EQ.1) THEN
      DO i=1, Physics%VNUM
        key = TRIM(Physics%pvarname(i)) // "_pstates"
        CALL SetAttr(IO, TRIM(key), &
                     this%prim%data5d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:,i))
      END DO
    END IF
    CALL GetAttr(config, "output/cstates", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      DO i=1, Physics%VNUM
        key = TRIM(Physics%cvarname(i)) // "_cstates"
        CALL SetAttr(IO, TRIM(key), &
                     this%cons%data5d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:,i))
      END DO
    END IF
    CALL GetAttr(config, "output/pfluxes", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      DO i=1, Physics%VNUM
        key = TRIM(Physics%pvarname(i)) // "_pfluxes"
        CALL SetAttr(IO, TRIM(key), &
                     this%pfluxes(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,:,i))
      END DO
    END IF

    CALL GetAttr(config, "output/wave_speeds", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      CALL SetAttr(IO, "minwav", &
                   this%minwav%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1:Mesh%NDIMS))
      CALL SetAttr(IO, "maxwav", &
                   this%maxwav%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1:Mesh%NDIMS))
      END IF

    ! initialize reconstruction modules
    CALL new_reconstruction(this%Reconstruction,Mesh,Physics,config,IOrec)

    ! add output data of reconstruction modules to IO dictionary
    IF(ASSOCIATED(IOrec)) CALL SetAttr(IO,"reconstruction",IOrec)


    ! initialize boundary fluxes
    this%bxflux(:,:,:,:) = 0.
    this%byflux(:,:,:,:) = 0.
    this%bzflux(:,:,:,:) = 0.
    this%bxfold(:,:,:,:) = 0.
    this%byfold(:,:,:,:) = 0.
    this%bzfold(:,:,:,:) = 0.
  END SUBROUTINE InitFluxes

  !> Get fluxes at boundaries
  !!
  !! \todo MPI communication
  FUNCTION GetBoundaryFlux(this,Mesh,Physics,direction,comm) RESULT(bflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_base),  INTENT(IN) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    INTEGER                         :: direction
    REAL, DIMENSION(Physics%VNUM)   :: bflux
    INTEGER, OPTIONAL               :: comm        ! communicator for MPI
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    REAL, DIMENSION(Physics%VNUM) :: bflux_local
    INTEGER                       :: sender_rank(1),dest_ranks(1),rank0(1)
    INTEGER                       :: dest_comm,union_comm
    INTEGER                       :: world_group,dest_group,union_group,sender_group
    INTEGER                       :: ierror
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)                    :: direction
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    IF (PRESENT(comm)) THEN
       dest_comm = comm
    ELSE
       dest_comm = MPI_COMM_NULL
    END IF
#endif
    SELECT CASE(direction)
    CASE(WEST) ! western boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(1).EQ.0) THEN
          bflux_local(:) = SUM(SUM(this%bxflux(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1,:),1),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(SUM(this%bxflux(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1,:),1),1)
#endif
    CASE(EAST) ! eastern boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(1).EQ.Mesh%dims(1)-1) THEN
          bflux_local(:) = SUM(SUM(this%bxflux(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2,:),1),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(SUM(this%bxflux(Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2,:),1),1)
#endif
    CASE(SOUTH) ! southern boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(2).EQ.0) THEN
          bflux_local(:) = SUM(SUM(this%byflux(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1,:),1),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(SUM(this%byflux(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,1,:),1),1)
#endif
    CASE (NORTH) ! northern boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(2).EQ.Mesh%dims(2)-1) THEN
          bflux_local(:) = SUM(SUM(this%byflux(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2,:),1),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(SUM(this%byflux(Mesh%KMIN:Mesh%KMAX,Mesh%IMIN:Mesh%IMAX,2,:),1),1)
#endif
    CASE (BOTTOM) ! bottom boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(3).EQ.0) THEN
          bflux_local(:) = SUM(SUM(this%bzflux(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1,:),1),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(SUM(this%bzflux(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1,:),1),1)
#endif
    CASE (TOP) ! topper boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(3).EQ.Mesh%dims(3)-1) THEN
          bflux_local(:) = SUM(SUM(this%bzflux(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2,:),1),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(SUM(this%bzflux(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2,:),1),1)
#endif
    CASE DEFAULT
       CALL this%Error("GetBoundaryFlux","wrong direction")
    END SELECT

#ifdef PARALLEL
    ! if dest_comm is the world comm, use simpler and faster AllReduce
    IF(dest_comm.EQ.MPI_COMM_WORLD) THEN
      CALL MPI_AllReduce(bflux_local,bflux,Physics%VNUM,DEFAULT_MPI_REAL, &
           MPI_SUM,MPI_COMM_WORLD,ierror)
    ELSE
      ! create world group
      CALL MPI_Comm_group(MPI_COMM_WORLD,world_group,ierror)
      ! determine the destination for the result
      IF(dest_comm.NE.MPI_COMM_NULL) THEN
        dest_comm = comm
        ! set destination group
        CALL MPI_Comm_group(dest_comm,dest_group,ierror)
      ELSE
        ! if no communicator is given, send the result
        ! to the rank0 process
        dest_ranks(1) = 0
        CALL MPI_Group_incl(world_group,1,dest_ranks,dest_group,ierror)
      END IF
      ! collect and sum up the result in process with rank 0 with respect to the
      ! subset of MPI processes at the boundary under consideration
      IF (Mesh%comm_boundaries(direction).NE.MPI_COMM_NULL) THEN
         CALL MPI_Reduce(bflux_local,bflux,Physics%VNUM,DEFAULT_MPI_REAL, &
              MPI_SUM,0,Mesh%comm_boundaries(direction),ierror)
      ELSE
         bflux(:) = 0.0
      END IF
      ! get sender group
      rank0(1) = Mesh%rank0_boundaries(direction)
      CALL MPI_Group_incl(world_group,1,rank0,sender_group,ierror)
      ! merge sender with destination
      CALL MPI_Group_union(sender_group,dest_group,union_group,ierror)
      ! create a communicator for the union group
      CALL MPI_Comm_create(MPI_COMM_WORLD,union_group,union_comm,ierror)
      IF (union_comm.NE.MPI_COMM_NULL) THEN
         ! get rank of sender in union group
         rank0(1) = 0
         CALL MPI_Group_translate_ranks(sender_group,1,rank0,union_group,sender_rank,ierror)
         IF (sender_rank(1).EQ.MPI_UNDEFINED) &
             CALL this%Error("GetBoundaryFlux","sender rank undefined")
         ! send result to all processes in communicator 'union_comm'
         CALL MPI_Bcast(bflux,Physics%VNUM,DEFAULT_MPI_REAL,sender_rank(1),union_comm,ierror)
         ! free union communicator
         CALL MPI_Comm_free(union_comm,ierror)
      END IF
      ! free all groups
      CALL MPI_Group_free(union_group,ierror)
      CALL MPI_Group_free(sender_group,ierror)
      CALL MPI_Group_free(world_group,ierror)
      CALL MPI_Group_free(dest_group,ierror)
    END IF

#endif
  END FUNCTION GetBoundaryFlux

  !> Calcualtes face data with reconstruction methods (e. g. limiters)
  PURE SUBROUTINE CalculateFaceData(this,Mesh,Physics,pvar,cvar)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_base),  INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar
    !------------------------------------------------------------------------!
    ! reconstruct data on cell faces
    IF (this%Reconstruction%PrimRecon()) THEN
       CALL this%Reconstruction%CalculateStates(Mesh,Physics,pvar,this%prim%data5d)
       CALL Physics%Convert2Conservative(this%prim,this%cons)
!        CALL Physics%Convert2Conservative(Mesh,this%prim%data5d,this%cons%data5d)
    ELSE
       CALL this%Reconstruction%CalculateStates(Mesh,Physics,cvar,this%cons%data5d)
       CALL Physics%Convert2Primitive(this%cons,this%prim)
!        CALL Physics%Convert2Primitive(Mesh,this%cons%data5d,this%prim%data5d)
    END IF

    ! update the speed of sound on cell faces (non-isotherml physics only)
    SELECT TYPE(phys => Physics)
    CLASS IS(physics_euler)
      SELECT TYPE(prim => this%prim)
      CLASS IS(statevector_euler)
        CALL phys%UpdateSoundSpeed(prim)
      END SELECT
    END SELECT

    ! get minimal & maximal wave speeds on cell interfaces
    CALL Physics%CalculateWaveSpeeds(Mesh,this%prim%data5d,this%cons%data5d,this%minwav,this%maxwav)
  END SUBROUTINE CalculateFaceData

  !> Destructor
  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("CloseFluxes","not initialized")
    CALL this%minwav%Destroy()
    CALL this%maxwav%Destroy()
    CALL this%prim%Destroy()
    CALL this%cons%Destroy()
    DEALLOCATE(this%cons,this%prim,this%pfluxes,this%minwav,this%maxwav, &
         this%bxflux,this%byflux,this%bzflux,this%bxfold,this%byfold,this%bzfold)
    CALL this%Reconstruction%Finalize()
    DEALLOCATE(this%Reconstruction)
  END SUBROUTINE Finalize_base

END MODULE fluxes_base_mod
