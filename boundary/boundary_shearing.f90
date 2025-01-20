  !#############################################################################
  !#                                                                           #
  !# fosite - 3D hydrodynamical simulation program                             #
  !# module: boundary_shearing.f90                                             #
  !#                                                                           #
  !# Copyright (C) 2006-2024                                                   #
  !# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
  !> \author Jannes Klee
  !!
  !! \brief Boundary module for a shearingsheet/shearingbox.
  !!        (see e.g. the \link shearingsheet.f90 standard run \endlink)
  !!
  !! \todo Implementation for shearing in third dimension is not included, yet.
  !!
  !! Implementation of the boundaries in a shearing box. For details have a
  !! look at \cite hawley1995 or \cite gammie2001 .
  !!
  !! The dependent variables \f$ f = (\Sigma, E, v_{x}, v_{y}) \f$ are sheared-
  !! periodic, thus for shearing in the western/eastern boundaries:
  !! \f[
  !!      f(x,y) = f(x,y+L) \\
  !!      f(x,y) = f(x+L,y-q \Omega L t).
  !! \f]
  !!
  !! There is an integral shift and a residual one. The latter is done by
  !! linear interpolation between two cells.
  !!
  !! \attention The shearing boundaries are either applying periodic boundaries
  !!            first or assuming that periodic boundaries are applied (parallel).
  !!
  !! <div class="row"> <div class="col-md-6">
  !! \image html boundaries_shift.png "illustration of boundary mapping"
  !! </div> <div class="col-md-6">
  !! \image html boundaries_velshift.png "velocity substraction for the v_y component"
  !! </div></div>
  !!
  !> \ingroup boundary
  !----------------------------------------------------------------------------!
MODULE boundary_shearing_mod
  USE boundary_base_mod
  USE boundary_periodic_mod
  USE marray_compound_mod
  USE mesh_base_mod
  USE physics_base_mod
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
  CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "shearing"

  TYPE, EXTENDS(boundary_periodic) :: boundary_shearing
    REAL                        :: velocity_offset !< velocity-shift distance L
    REAL, DIMENSION(:), POINTER :: velocity_shift  !< extra shift in one direction
  CONTAINS
    PROCEDURE :: InitBoundary_shearing
    PROCEDURE :: SetBoundaryData
    PROCEDURE :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: boundary_shearing
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor for shearing boundary conditions.
  SUBROUTINE InitBoundary_shearing(this,Mesh,Physics,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_shearing), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    TYPE(Dict_TYP), POINTER                 :: config
    INTEGER,                  INTENT(IN)    :: dir
    !------------------------------------------------------------------------!
    INTEGER            :: err, l
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,SHEARING,boundcond_name,dir,config)

    ALLOCATE( &
             this%velocity_shift(Physics%VNUM+Physics%PNUM), &
             STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("boundary_shearing::InitBoundaryg", "Unable to allocate memory.")
    END IF

    DO l=1,Physics%VNUM + Physics%PNUM
      ! this part for WEST-EAST shear
      IF (Mesh%shear_dir.EQ.2) THEN
        IF (l.EQ.Physics%YVELOCITY) THEN
          SELECT CASE(Mesh%fargo%GetType())
          CASE(0) ! fargo disabled
            this%velocity_shift(l) = Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)
          CASE(3) ! shearing sheet/box fargo enabled
            this%velocity_shift(l) = 0.0
          CASE DEFAULT
            CALL this%Error("boundary_shearing::InitBoundary", &
              "Shearing boundaries are only compatible without Fargo or with Fargo type 3 (shearing box).")
          END SELECT
        ELSE
            this%velocity_shift(l) = 0.0
        END IF
      ! this part for SOUTH-NORTH shear
      ELSE IF (Mesh%shear_dir.EQ.1) THEN
        IF (l.EQ.Physics%XVELOCITY) THEN
          SELECT CASE(Mesh%fargo%GetType())
          CASE(0) ! fargo disabled
            this%velocity_shift(l) = Mesh%Q*Mesh%OMEGA*(Mesh%ymax-Mesh%ymin)
          CASE(3) ! shearing sheet/box fargo enabled
            this%velocity_shift(l) = 0.0
          CASE DEFAULT
            CALL this%Error("boundary_shearing::InitBoundary", &
              "Shearing boundaries are only compatible without Fargo or with Fargo type 3 (shearing box).")
          END SELECT
        ELSE
            this%velocity_shift(l) = 0.0
        END IF
      ELSE
        CALL this%Error("boundary_shearing::InitBoundary", &
          "Shearing boundaries in top/south direction not allowed, yet.")
      END IF
    END DO
    this%velocity_offset = Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)/Mesh%dy

  END SUBROUTINE InitBoundary_shearing

  !> \public Applies the shearing boundary conditions.
  !!
  !! The physical meaning of the implementation is explained in the
  !! overall module descriptions. Some of the data required for setting
  !! this boundary condition is calculated during initialization
  !! (see \ref boundary_shearing_mod.initboundary_shearing ) and depends
  !! on whether fast advection (FARGO) (see \ref timedisc_base_mod )
  !! is enabled or not.
  !!
  !! The implemenation is done in a way, that periodic boundaries are applied at
  !! the very beginning (in serial) or are already applied (in parallel). Thus,
  !! only the shifting on the same side is applied. This gives us the benefit to be able
  !! to reuse the periodic boundaries, that are automatically applied by MPI from
  !! of domain decomposition. When running on a single core the periodic
  !! boundaries are applied by hand at the beginning.
  !!
  !! For parallel usage: currently, this implementation only works
  !! in parallel with processes along the y-axis.
  !!
  !!                                 NB
  !!                          ----------------
  !!                                 p1
  !!                          ----------------
  !!                                 p2
  !!                          ----------------
  !!                                 p3
  !!                          ----------------
  !!                                 p4
  !!                          ----------------
  !!                                 SB
  !!
  !! Here, NB, SB are northern and southern boundaries, respectively.
  !! p1,p2,... are the used processes. This kind of parallization goes in
  !! concordance with the parallelization used for the gravitation fourier
  !! solvers, which need full strides in one direction. We use solely the
  !! x-direction, because it is the first dimension and lies coherently behind
  !! each other for vectorization and generally fast access.
  SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_shearing), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    REAL,                     INTENT(IN)    :: time
    CLASS(marray_compound),   INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,l,intshift
    REAL               :: offremain,offset,offset_tmp
#ifdef PARALLEL
    INTEGER            :: status(MPI_STATUS_SIZE)
    INTEGER            :: ierror,req(4)
    CHARACTER(LEN=80)  :: str
    REAL               :: mpi_buf(2*Mesh%GNUM)
#endif
    !------------------------------------------------------------------------!
    ! the routine below expect that periodic boundaries are already applied
#ifndef PARALLEL
    CALL this%boundary_periodic%SetBoundaryData(Mesh,Physics,time,pvar)
#endif

    ! make sure all MPI processes use the same step if domain is decomposed
    ! along the y-direction (can be different due to round-off errors)
    offset = -this%velocity_offset*time
    DO WHILE(offset<0)
      offset = offset + Mesh%JNUM
    END DO

    DO l=1,Physics%VNUM+Physics%PNUM
      ! chose velocity offset for fargo and considered pvar
      SELECT CASE(this%GetDirection())
      CASE(WEST)
        offset_tmp = offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO i=1,Mesh%GNUM
          !------- integral shift ---------------------------------------------!
          pvar%data4d(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar%data4d(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l),intshift)
          !------- residual shift ---------------------------------------------!
          pvar%data4d(Mesh%IMIN-i,Mesh%JMAX+1,:,l) = pvar%data4d(Mesh%IMIN-i,Mesh%JMIN,:,l)
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
              pvar%data4d(Mesh%IMIN-i,j,k,l) = (1.0 - offremain)*pvar%data4d(Mesh%IMIN-i,j,k,l) + &
                offremain*pvar%data4d(Mesh%IMIN-i,j+1,k,l) + this%velocity_shift(l)
            END DO
          END DO
        END DO
      CASE(EAST)
        offset_tmp = -offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO i=1,Mesh%GNUM
          !------- integral shift ---------------------------------------------!
          pvar%data4d(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar%data4d(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l),intshift)
          !------- residual shift ---------------------------------------------!
          pvar%data4d(Mesh%IMAX+i,Mesh%JMAX+1,:,l) = pvar%data4d(Mesh%IMAX+i,Mesh%JMIN,:,l)
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
              pvar%data4d(Mesh%IMAX+i,j,k,l) = (1.0 - offremain)*pvar%data4d(Mesh%IMAX+i,j,k,l) + &
                offremain*pvar%data4d(Mesh%IMAX+i,j+1,k,l) - this%velocity_shift(l)
            END DO
          END DO
        END DO
      CASE(SOUTH)
        offset_tmp = -offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO j=1,Mesh%GNUM
          !------- integral shift ---------------------------------------------!
          pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,l),intshift)
          !------- residual shift ---------------------------------------------!
          pvar%data4d(Mesh%IMAX+1,Mesh%JMIN-j,:,l) = pvar%data4d(Mesh%IMIN,Mesh%JMIN-j,:,l)
          DO k=Mesh%KMIN,Mesh%KMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              pvar%data4d(i,Mesh%JMIN-j,k,l) = (1.0 - offremain)*pvar%data4d(i,Mesh%JMIN-j,k,l) + &
                offremain*pvar%data4d(i+1,Mesh%JMIN-j,k,l) - this%velocity_shift(l)
            END DO
          END DO
        END DO
      CASE(NORTH)
        offset_tmp = offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO j=1,Mesh%GNUM
          !------- integral shift ---------------------------------------------!
          pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,l),intshift)
          !------- residual shift ---------------------------------------------!
          pvar%data4d(Mesh%IMAX+1,Mesh%JMAX+j,:,l) = pvar%data4d(Mesh%IMIN,Mesh%JMAX+j,:,l)
          DO k=Mesh%KMIN,Mesh%KMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              pvar%data4d(i,Mesh%JMAX+j,k,l) = (1.0 - offremain)*pvar%data4d(i,Mesh%JMAX+j,k,l) + &
                offremain*pvar%data4d(i+1,Mesh%JMAX+j,k,l) + this%velocity_shift(l)
            END DO
          END DO
        END DO
      CASE(BOTTOM)
       !CALL this%Error("SetBoundary", "Shearing not supported in BOTTOM direction.")
      CASE(TOP)
       !CALL this%Error("SetBoundary", "Shearing not supported in TOP direction.")
      END SELECT
    END DO
  END SUBROUTINE SetBoundaryData


  !> \public Destructor for periodic boundary conditions
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_shearing), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%velocity_shift)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize



END MODULE boundary_shearing_mod
