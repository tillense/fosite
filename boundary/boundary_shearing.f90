  !#############################################################################
  !#                                                                           #
  !# fosite - 2D hydrodynamical simulation program                             #
  !# module: boundary_shearing.f90                                             #
  !#                                                                           #
  !# Copyright (C) 2006-2015                                                   #
  !# Jannes Klee <jklee@astrophysik.uni-kiel.de                                #
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
  !! <div class="row"> <div class="col-md-6">
  !! \image html boundaries_shift.png "illustration of boundary mapping"
  !! </div> <div class="col-md-6">
  !! \image html boundaries_velshift.png "velocity substraction for the v_y component"
  !! </div></div>
  !!
  !! \extends boundary_common
  !! \ingroup boundary
  !----------------------------------------------------------------------------!
  MODULE boundary_shearing
    USE mesh_common, ONLY : Mesh_TYP
    USE physics_common, ONLY : Physics_TYP
    USE common_dict
    USE boundary_nogradients
    USE boundary_common
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
    CHARACTER(LEN=32), PARAMETER  :: boundcond_name = "shearing"
    !--------------------------------------------------------------------------!
    PUBLIC :: &
         ! types
         Boundary_TYP, &
         ! constants
         WEST, EAST, SOUTH, NORTH, &
         ! methods
         InitBoundary_shearing, &
         CenterBoundary_shearing
    !--------------------------------------------------------------------------!

  CONTAINS

  !> \public Constructor for shearing boundary conditions.
  SUBROUTINE InitBoundary_shearing(this,btype,dir,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Dict_TYP),POINTER &
                       :: config
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
    INTENT(IN)         :: btype,dir
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitBoundary(this,btype,boundcond_name,dir)
  END SUBROUTINE InitBoundary_shearing

  !> \public Applies the shearing boundary conditions.
  !!
  !! The physical meaning of the implementation is explained in the
  !! overall module descriptions. To the implementation used values which
  !! are calculated during initialization of the \link timedisc_generic \endlink
  !! module.
  SUBROUTINE CenterBoundary_shearing(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    REAL               :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                       :: pvar
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k
    REAL               :: velocity_shift,offremain,offset
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,time
    INTENT(INOUT)      :: pvar
    !------------------------------------------------------------------------!
    DO k=1,Physics%VNUM
      ! chose velocity offset for fargo and considered pvar
      IF (k == Physics%YVELOCITY) THEN
        velocity_shift = this%velocity_shift
      ELSE
        velocity_shift = 0.0
      END IF

      SELECT CASE(GetDirection(this))
      CASE(WEST)
      offset = -this%velocity_offset*time
      offremain = offset - FLOOR(offset)
!CDIR NODEP
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=1,Mesh%GNUM
            pvar(Mesh%IMIN-i,j,k) = &
              (1.0 - offremain) * &
              pvar(Mesh%IMAX-i+1,1+MODULO(j-1+FLOOR(offset),Mesh%JNUM),k) + &
              offremain*pvar(Mesh%IMAX-i+1,1+MODULO(j+FLOOR(offset),Mesh%JNUM),k) + &
              velocity_shift
          END DO
        END DO
      CASE(EAST)
      offset = this%velocity_offset*time
      offremain = offset - FLOOR(offset)
!CDIR NODEP
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=1,Mesh%GNUM
            pvar(Mesh%IMAX+i,j,k) = &
              (1.0 - offremain) * &
              pvar(Mesh%IMIN+i-1,1+MODULO(j-1+FLOOR(offset),Mesh%JNUM),k) + &
              offremain*pvar(Mesh%IMIN+i-1,1+MODULO(j+FLOOR(offset),Mesh%JNUM),k) - &
              velocity_shift
          END DO
        END DO
      END SELECT
    END DO

    ! for southern and northern boundaries use periodic
    SELECT CASE(GetDirection(this))
    CASE(SOUTH)
      CALL Error(this,"boundary_shearing", &
                 "Shearing boundaries not distinguishable from periodic" // &
                 "ones in southern direction. Use periodic boundaries!")
    CASE(NORTH)
      CALL Error(this,"boundary_shearing", &
                 "Shearing boundaries not distinguishable from periodic" // &
                 "ones in northern direction. Use periodic boundaries!")
    END SELECT
  END SUBROUTINE CenterBoundary_shearing
END MODULE boundary_shearing
