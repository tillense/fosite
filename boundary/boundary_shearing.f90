  !#############################################################################
  !#                                                                           #
  !# fosite - 3D hydrodynamical simulation program                             #
  !# module: boundary_shearing.f03                                             #
  !#                                                                           #
  !# Copyright (C) 2006-2018                                                   #
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
  !> \addtogroup boundary
  !----------------------------------------------------------------------------!
  !> \author Jannes Klee
  !!
  !! \brief Boundary module for a shearing box
  !!        (see e.g. the \link shearingsheet.f90 standard run \endlink)
  !!
  !! Implementation of the boundaries in a shearing box. For details have a
  !! look at \cite hawley1995 or \cite gammie2001 . The northern and southern boundary
  !! conditions are periodic boundaries.
  !!
  !! The dependent variables \f$ f = (\Sigma, E, v_{x}, v_{y}) \f$ are sheared-
  !! periodic, thus:
  !! \f[
  !!      f(x,y) = f(x,y+L) \\
  !!      f(x,y) = f(x+L,y-q \Omega L t).
  !! \f]
  !!
  !! Since the sheared cells move continuously in time, linear interpolation
  !! is used, in order to set the boundary in the fixed ghost cells.
  !! For implementation it is important to notice that there is an integral
  !! shift and a rest one. The latter is done by linear interpolation between
  !! two cells.
  !!
  !! \attention For a good numerical setup, the shearing boundary should be applied
  !!            after the other boundaries, because the shearing accesses the edge cells.
  !!
  !! <div class="row"> <div class="col-md-6">
  !! \image html boundaries_shift.png "illustration of boundary mapping"
  !! </div> <div class="col-md-6">
  !! \image html boundaries_velshift.png "velocity substraction for the v_y component"
  !! </div></div>
  !----------------------------------------------------------------------------!
MODULE boundary_shearing_mod
  USE boundary_base_mod
  USE boundary_periodic_mod
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
    REAL, DIMENSION(:), POINTER :: velocity_shift  !< shear box parameter
    REAL                        :: velocity_offset !< shear box parameter
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
    INTEGER            :: btype, err, l
    !------------------------------------------------------------------------!
    CALL this%InitBoundary(Mesh,Physics,SHEARING,boundcond_name,dir,config)

    ALLOCATE( &
             this%velocity_shift(Physics%VNUM+Physics%PNUM), &
             STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("InitBoundary_shearing", "Unable to allocate memory.")
    END IF

    DO l=1,Physics%VNUM + Physics%PNUM
      ! this part for WEST-EAST shear
      IF (Mesh%WE_shear) THEN
        IF (l.EQ.Physics%YVELOCITY) THEN
          IF (Mesh%FARGO.EQ.0) THEN
            this%velocity_shift(l) = Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)
          ELSE IF (Mesh%FARGO.EQ.3) THEN
            this%velocity_shift(l) = 0.0
          ELSE
            CALL this%Error("InitTimedisc", &
            "Shearing boundaries are only compatible without Fargo or Fargo type 3 (shearing box).")
          END IF
        ELSE
            this%velocity_shift(l) = 0.0
        END IF
      ! this part for SOUTH-NORTH shear
      ELSE IF (Mesh%SN_shear) THEN
        IF (l.EQ.Physics%XVELOCITY) THEN
          IF (Mesh%FARGO.EQ.0) THEN
            this%velocity_shift(l) = Mesh%Q*Mesh%OMEGA*(Mesh%ymax-Mesh%ymin)
          ELSE IF (Mesh%FARGO.EQ.3) THEN
            this%velocity_shift(l) = 0.0
          ELSE
            CALL this%Error("InitTimedisc", &
            "Shearing boundaries are only compatible without Fargo or Fargo type 3 (shearing box).")
          END IF
        ELSE
            this%velocity_shift(l) = 0.0
        END IF
      ELSE
        CALL this%Error("InitBoundary", &
        "Shearing boundaries in top/south direction not allowed, yet.")
      END IF
    END DO
    this%velocity_offset = Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)/Mesh%dy

  END SUBROUTINE InitBoundary_shearing

  !> \public Applies the shearing boundary conditions.
  !!
  !! The physical meaning of the implementation is explained in the
  !! overall module descriptions. To the implementation used values which
  !! are calculated during initialization of the \link timedisc_generic \endlink
  !! module.
  !!
  !! \attention The implemenation is done in a way, that periodic boundaries
  !! are already applied. Thus, only the shifting on the **same** side is
  !! applied to the primitive variables. This gives us the benefit to be able
  !! to reuse the periodic boundaries automatically applied by MPI because
  !! of domain decomposition. When running on a single core the periodic
  !! boundaries are applied by hand at the beginning.
  !!
  !! \attention Currently, this implementation only works parallelized with processes
  !!  along the y-axis, which means that
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
  !! \attention Here, NB, SB are northern and southern boundaries, respectively.
  !! p1,p2,... are the used processes. This kind of parallization goes in
  !! concordance with the parallelization used for the gravitation fourier
  !! solvers, which need full strides in one direction. We use solely the
  !! x-direction, because it is the first dimension and lies coherently behind
  !! each other for vectorization and generally fast access.
  PURE SUBROUTINE SetBoundaryData(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(boundary_shearing), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    REAL,                     INTENT(IN)    :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX, &
                    Physics%VNUM+Physics%PNUM), &
                              INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,l,intshift
    REAL               :: velocity_shift,offremain,offset,offset_tmp
    REAL               :: pvar_old,pvar_old2
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
          !------- residual shift ---------------------------------------------!
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
              pvar(Mesh%IMIN-i,j,k,l) = (1.0 - offremain)*pvar(Mesh%IMIN-i,j,k,l) + &
                offremain*pvar(Mesh%IMIN-i,j+1,k,l) + this%velocity_shift(l)
            END DO
          END DO
          !------- integral shift ---------------------------------------------!
          pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l),intshift)
        END DO
      CASE(EAST)
        offset_tmp = -offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO i=1,Mesh%GNUM
          !------- residual shift ---------------------------------------------!
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
              pvar(Mesh%IMAX+i,j,k,l) = (1.0 - offremain)*pvar(Mesh%IMAX+i,j,k,l) + &
                offremain*pvar(Mesh%IMAX+i,j+1,k,l) - this%velocity_shift(l)
            END DO
          END DO
          !------- integral shift ---------------------------------------------!
          pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l),intshift)
        END DO
      CASE(SOUTH)
        offset_tmp = -offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO j=1,Mesh%GNUM
          !------- residual shift ---------------------------------------------!
          DO k=Mesh%KMIN,Mesh%KMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              pvar(i,Mesh%JMIN-j,k,l) = (1.0 - offremain)*pvar(i,Mesh%JMIN-j,k,l) + &
                offremain*pvar(i+1,Mesh%JMIN-j,k,l) - this%velocity_shift(l)
            END DO
          END DO
          !------- integral shift ---------------------------------------------!
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX,l),intshift)
        END DO
      CASE(NORTH)
        offset_tmp = offset
        intshift = FLOOR(offset_tmp)
        offremain = offset_tmp - intshift
        DO j=1,Mesh%GNUM
          !------- residual shift ---------------------------------------------!
          DO k=Mesh%KMIN,Mesh%KMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              pvar(i,Mesh%JMAX+j,k,l) = (1.0 - offremain)*pvar(i,Mesh%JMAX+j,k,l) + &
                offremain*pvar(i+1,Mesh%JMAX+j,k,l) + this%velocity_shift(l)
            END DO
          END DO
          !------- integral shift ---------------------------------------------!
          pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,l)  = &
            CSHIFT(pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX,l),intshift)
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
