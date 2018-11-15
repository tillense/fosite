!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: fluxes_kt.f90                                                     #
!#                                                                           #
!# Copyright (C) 2007-2017                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
!# Jubin Lirawi <jlirawi@astrophysik.uni-kiel.de>                            #
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
!! \author Jubin Lirawi
!!
!! \brief numerical fluxes module for midpoint quadrature rule
!!
!! \extends fluxes_common
!! \ingroup fluxes
!----------------------------------------------------------------------------!
MODULE fluxes_kt_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE marray_compound_mod
!#ifdef PARALLEL
!       DEFAULT_MPI_REAL
!#endif
  USE physics_base_mod
  USE reconstruction_base_mod
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
  TYPE, EXTENDS (fluxes_base) :: fluxes_kt
  CONTAINS
    PROCEDURE         :: InitFluxes_kt
    PROCEDURE         :: CalculateFluxes
    PROCEDURE         :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       fluxes_kt
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes_kt(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_kt),    INTENT(INOUT)       :: this
    CLASS(mesh_base),    INTENT(IN)          :: Mesh
    CLASS(physics_base), INTENT(IN)          :: Physics
    TYPE(Dict_TYP),      INTENT(IN), POINTER :: config,IO
    !------------------------------------------------------------------------!
    CALL this%InitFluxes(Mesh,Physics,config,IO,KT,"KT")

    ! set relative positions for reconstruction points:
    ! cell face positions
    SELECT CASE(Mesh%GetType())
    CASE(MIDPOINT)
         this%dx(:,:,:,1) = Mesh%curv%faces(:,:,:,1,1) - Mesh%curv%bcenter(:,:,:,1)
         this%dx(:,:,:,2) = Mesh%curv%faces(:,:,:,2,1) - Mesh%curv%bcenter(:,:,:,1)
         this%dy(:,:,:,1:2) = 0.0
         this%dz(:,:,:,1:2) = 0.0
         this%dx(:,:,:,3:4) = 0.0
         this%dy(:,:,:,3) = Mesh%curv%faces(:,:,:,3,2) - Mesh%curv%bcenter(:,:,:,2)
         this%dy(:,:,:,4) = Mesh%curv%faces(:,:,:,4,2) - Mesh%curv%bcenter(:,:,:,2)
         this%dz(:,:,:,3:4) = 0.0
         this%dx(:,:,:,5:6) = 0.0
         this%dy(:,:,:,5:6) = 0.0
         this%dz(:,:,:,5) = Mesh%curv%faces(:,:,:,5,3) - Mesh%curv%bcenter(:,:,:,3)
         this%dz(:,:,:,6) = Mesh%curv%faces(:,:,:,6,3) - Mesh%curv%bcenter(:,:,:,3)
    CASE(TRAPEZOIDAL)
      CALL this%Error("InitFluxes_kt", "Trapezoidal not supported from > 0.6")
    CASE DEFAULT
      CALL this%Error("InitFluxes_kt", "Unknown mesh type.")
    END SELECT
  END SUBROUTINE InitFluxes_kt


  !> Calculates the fluxes with the midpoint method
  !!
  !! \warning{The return values xfluxdydz, yfluxdzdx, zfluxdxdy are the numerical fluxes
  !!          devided by dy or dx respectively. This reduces numerical errors
  !!          because otherwise we would multiply the fluxes by dy (or dx) here and
  !!          devide by dy (or dx) later when computing flux differences.}
  !! \warning{In an older version of Fosite, there was still the possibility to
  !!          use the trapezoidal rule. This is not supported anymore from version
  !!          <=0.6.}
  PURE SUBROUTINE CalculateFluxes(this,Mesh,Physics,pvar,cvar, &
                  xfluxdydz,yfluxdzdx,zfluxdxdy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_kt),    INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%vnum), &
                         INTENT(OUT)   :: xfluxdydz,yfluxdzdx,zfluxdxdy
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,l,m
    !------------------------------------------------------------------------!
    ! execute generic tasks common to all flux types
    CALL this%CalculateFaceData(Mesh,Physics,pvar,cvar)

    SELECT CASE(Mesh%GetType())
    CASE(MIDPOINT)
      ! compute numerical fluxes along x-direction (west and east) devided by dy
      IF (Mesh%INUM.GT.1) THEN
        m = 1 ! index counting the number of spatial directions for which transport is enabled
        ! physical fluxes
        CALL Physics%CalculateFluxesX(Mesh,2*m-1,2*m,this%prim%data5d,this%cons%data5d,this%pfluxes)
!NEC$ UNROLL(8)
        DO l=1,Physics%VNUM
!NEC$ UNROLL(8)
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IMIN+Mesh%IM1,Mesh%IMAX
                  xfluxdydz(i,j,k,l) = Mesh%dAxdydz(i+1,j,k,1) / &
                         (this%maxwav%data4d(i,j,k,m) - this%minwav%data4d(i,j,k,m)) * &
                         (this%maxwav%data4d(i,j,k,m) * this%pfluxes(i,j,k,m+1,l) - &
                          this%minwav%data4d(i,j,k,m) * this%pfluxes(i+1,j,k,m,l) + &
                          this%minwav%data4d(i,j,k,m) * this%maxwav%data4d(i,j,k,m) * &
                      (this%cons%data5d(i+1,j,k,m,l) - this%cons%data5d(i,j,k,m+1,l)))
              END DO
            END DO
          END DO
        END DO
      ELSE
         m = 0
         xfluxdydz(:,:,:,:) = 0.0
      END IF

      ! compute numerical fluxes along y-direction (south and north) devided by dx
      IF (Mesh%JNUM.GT.1) THEN
        m = m + 1 ! increase wave speed index, may be 1 or 2 now depending
                  ! on whether there was transport in x-direction or not
        ! physical fluxes
        CALL Physics%CalculateFluxesY(Mesh,2*m-1,2*m,this%prim%data5d,this%cons%data5d,this%pfluxes)
!NEC$ UNROLL(8)
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JMIN+Mesh%JM1,Mesh%JMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                  yfluxdzdx(i,j,k,l) = Mesh%dAydzdx(i,j+1,k,1) / &
                         (this%maxwav%data4d(i,j,k,m) - this%minwav%data4d(i,j,k,m)) * &
                         (this%maxwav%data4d(i,j,k,m) * this%pfluxes(i,j,k,2*m,l) - &
                          this%minwav%data4d(i,j,k,m) * this%pfluxes(i,j+1,k,2*m-1,l) + &
                          this%minwav%data4d(i,j,k,m) * this%maxwav%data4d(i,j,k,m) * &
                      (this%cons%data5d(i,j+1,k,2*m-1,l) - this%cons%data5d(i,j,k,2*m,l)))
              END DO
            END DO
          END DO
        END DO
      ELSE
         yfluxdzdx(:,:,:,:) = 0.0
      END IF

      IF (Mesh%KNUM.GT.1) THEN
        m = m + 1 ! increase wave speed index, may be 1, 2 or 3 now depending
                  ! on whether there was transport in x- and/or y-direction or not
         ! physical fluxes
        CALL Physics%CalculateFluxesZ(Mesh,2*m-1,2*m,this%prim%data5d,this%cons%data5d,this%pfluxes)
!NEC$ UNROLL(8)
        DO l=1,Physics%VNUM
!NEC$ UNROLL(8)
          DO k=Mesh%KMIN+Mesh%KM1,Mesh%KMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN,Mesh%IGMAX
                  zfluxdxdy(i,j,k,l) = Mesh%dAzdxdy(i,j,k+1,1) / &
                         (this%maxwav%data4d(i,j,k,m) - this%minwav%data4d(i,j,k,m)) * &
                         (this%maxwav%data4d(i,j,k,m) * this%pfluxes(i,j,k,2*m,l) - &
                          this%minwav%data4d(i,j,k,m) * this%pfluxes(i,j,k+1,2*m-1,l) + &
                          this%minwav%data4d(i,j,k,m) * this%maxwav%data4d(i,j,k,m) * &
                      (this%cons%data5d(i,j,k+1,2*m-1,l) - this%cons%data5d(i,j,k,2*m,l)))
             END DO
            END DO
          END DO
        END DO
      ELSE
         zfluxdxdy(:,:,:,:) = 0.0
      END IF
    CASE(TRAPEZOIDAL)
      ! Here was once the trapezoidal rule implemented. See releases (<= 0.6) in
      ! order to see this. Since 3D fosite it is not supported anymore.
    END SELECT
  END SUBROUTINE CalculateFluxes

  ! Here was once a routine for BilinearInterpolation, which was needed for
  ! trapezoidal rule. See (<= 0.6) version of fosite.

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fluxes_kt), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE fluxes_kt_mod
