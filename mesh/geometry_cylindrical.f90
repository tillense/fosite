!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_cylindrical.f90                                          #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
!#                    Jubin Lirawi     <jlirawi@astrophysik.uni-kiel.de>     #
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
!! \author Jubin Lirawi
!! \author Jannes Klee
!!
!! \brief defines properties of a 3D cylindrical mesh
!!
!! \extends geometry_cylindrical
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_cylindrical_mod
  USE geometry_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS (geometry_base) :: geometry_cylindrical
  CONTAINS
    PROCEDURE :: InitGeometry_cylindrical
    PROCEDURE :: ScaleFactors_1
    PROCEDURE :: ScaleFactors_2
    PROCEDURE :: ScaleFactors_3
    PROCEDURE :: ScaleFactors_4
    PROCEDURE :: Radius_1
    PROCEDURE :: Radius_2
    PROCEDURE :: Radius_3
    PROCEDURE :: Radius_4
    PROCEDURE :: PositionVector_1
    PROCEDURE :: PositionVector_2
    PROCEDURE :: PositionVector_3
    PROCEDURE :: PositionVector_4
    PROCEDURE :: Convert2Cartesian_coords_1
    PROCEDURE :: Convert2Cartesian_coords_2
    PROCEDURE :: Convert2Cartesian_coords_3
    PROCEDURE :: Convert2Cartesian_coords_4
    PROCEDURE :: Convert2Curvilinear_coords_1
    PROCEDURE :: Convert2Curvilinear_coords_2
    PROCEDURE :: Convert2Curvilinear_coords_3
    PROCEDURE :: Convert2Curvilinear_coords_4
    PROCEDURE :: Convert2Cartesian_vectors_1
    PROCEDURE :: Convert2Cartesian_vectors_2
    PROCEDURE :: Convert2Cartesian_vectors_3
    PROCEDURE :: Convert2Cartesian_vectors_4
    PROCEDURE :: Convert2Curvilinear_vectors_1
    PROCEDURE :: Convert2Curvilinear_vectors_2
    PROCEDURE :: Convert2Curvilinear_vectors_3
    PROCEDURE :: Convert2Curvilinear_vectors_4
    PROCEDURE :: Finalize
  END TYPE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "cylindrical"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_cylindrical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_cylindrical(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(INOUT) :: this
    TYPE(DICT_TYP), POINTER                    :: config
    !------------------------------------------------------------------------!
    CALL this%InitGeometry(CYLINDRICAL,geometry_name,config)
    CALL this%SetAzimuthIndex(2)
  END SUBROUTINE InitGeometry_cylindrical

  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,1),hx(:),hy(:),hz(:))
  END SUBROUTINE ScaleFactors_1

  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,1),hx(:,:),hy(:,:),hz(:,:))
  END SUBROUTINE ScaleFactors_2

  PURE SUBROUTINE ScaleFactors_3(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,:,1),hx(:,:,:),hy(:,:,:),hz(:,:,:))
  END SUBROUTINE ScaleFactors_3

  PURE SUBROUTINE ScaleFactors_4(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,:,:,1),hx(:,:,:,:),hy(:,:,:,:),hz(:,:,:,:))
  END SUBROUTINE ScaleFactors_4

  PURE SUBROUTINE Radius_1(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,1),coords(:,3))
  END SUBROUTINE Radius_1

  PURE SUBROUTINE Radius_2(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,1),coords(:,:,3))
  END SUBROUTINE Radius_2

  PURE SUBROUTINE Radius_3(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,:,1),coords(:,:,:,3))
  END SUBROUTINE Radius_3

  PURE SUBROUTINE Radius_4(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,:,:,1),coords(:,:,:,:,3))
  END SUBROUTINE Radius_4

  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,1),coords(:,3),posvec(:,1),posvec(:,2),posvec(:,3))
  END SUBROUTINE PositionVector_1

  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,:,1),coords(:,:,3),posvec(:,:,1), &
                        posvec(:,:,2),posvec(:,:,3))
  END SUBROUTINE PositionVector_2

  PURE SUBROUTINE PositionVector_3(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,:,:,1),coords(:,:,:,3),posvec(:,:,:,1), &
                        posvec(:,:,:,2),posvec(:,:,:,3))
  END SUBROUTINE PositionVector_3

  PURE SUBROUTINE PositionVector_4(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,:,:,:,1),coords(:,:,:,:,3),posvec(:,:,:,:,1), &
                        posvec(:,:,:,:,2),posvec(:,:,:,:,3))
  END SUBROUTINE PositionVector_4

  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,1),curv(:,2),curv(:,3), &
                                  cart(:,1),cart(:,2),cart(:,3))
  END SUBROUTINE Convert2Cartesian_coords_1

  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,:,1),curv(:,:,2),curv(:,:,3), &
                                  cart(:,:,1),cart(:,:,2),cart(:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_2

  PURE SUBROUTINE Convert2Cartesian_coords_3(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                  cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_3

  PURE SUBROUTINE Convert2Cartesian_coords_4(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                  cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_4

  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,1),cart(:,2),cart(:,3), &
                                    curv(:,1),curv(:,2),curv(:,3))
  END SUBROUTINE Convert2Curvilinear_coords_1

  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,:,1),cart(:,:,2),cart(:,:,3), &
                                    curv(:,:,1),curv(:,:,2),curv(:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_2

  PURE SUBROUTINE Convert2Curvilinear_coords_3(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3), &
                                    curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_3

  PURE SUBROUTINE Convert2Curvilinear_coords_4(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3), &
                                    curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_4

  PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_vectors(curv(:,1),curv(:,2),curv(:,3), &
                                  v_curv(:,1),v_curv(:,2),v_curv(:,3), &
                                  v_cart(:,1),v_cart(:,2),v_cart(:,3))
  END SUBROUTINE Convert2Cartesian_vectors_1

  PURE SUBROUTINE Convert2Cartesian_vectors_2(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_vectors(curv(:,:,1),curv(:,:,2),curv(:,:,3), &
                                  v_curv(:,:,1),v_curv(:,:,2),v_curv(:,:,3), &
                                  v_cart(:,:,1),v_cart(:,:,2),v_cart(:,:,3))
  END SUBROUTINE Convert2Cartesian_vectors_2

  PURE SUBROUTINE Convert2Cartesian_vectors_3(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_vectors(curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                  v_curv(:,:,:,1),v_curv(:,:,:,2),v_curv(:,:,:,3), &
                                  v_cart(:,:,:,1),v_cart(:,:,:,2),v_cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_vectors_3

  PURE SUBROUTINE Convert2Cartesian_vectors_4(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_vectors(curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                  v_curv(:,:,:,:,1),v_curv(:,:,:,:,2),v_curv(:,:,:,:,3), &
                                  v_cart(:,:,:,:,1),v_cart(:,:,:,:,2),v_cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_vectors_4

  PURE SUBROUTINE Convert2Curvilinear_vectors_1(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_vectors(curv(:,1),curv(:,2),curv(:,3), &
                                  v_cart(:,1),v_cart(:,2),v_cart(:,3), &
                                  v_curv(:,1),v_curv(:,2),v_curv(:,3))
  END SUBROUTINE Convert2Curvilinear_vectors_1

  PURE SUBROUTINE Convert2Curvilinear_vectors_2(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_vectors(curv(:,:,1),curv(:,:,2),curv(:,:,3), &
                                  v_cart(:,:,1),v_cart(:,:,2),v_cart(:,:,3), &
                                  v_curv(:,:,1),v_curv(:,:,2),v_curv(:,:,3))
  END SUBROUTINE Convert2Curvilinear_vectors_2

  PURE SUBROUTINE Convert2Curvilinear_vectors_3(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_vectors(curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                  v_cart(:,:,:,1),v_cart(:,:,:,2),v_cart(:,:,:,3), &
                                  v_curv(:,:,:,1),v_curv(:,:,:,2),v_curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_vectors_3

  PURE SUBROUTINE Convert2Curvilinear_vectors_4(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_vectors(curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                  v_cart(:,:,:,:,1),v_cart(:,:,:,:,2),v_cart(:,:,:,:,3), &
                                  v_curv(:,:,:,:,1),v_curv(:,:,:,:,2),v_curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_vectors_4

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

  ELEMENTAL SUBROUTINE ScaleFactors(r,hr,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL, INTENT(OUT) :: hr,hphi,hz
    !------------------------------------------------------------------------!
    hr   = 1.
    hphi = r
    hz   = 1.
  END SUBROUTINE ScaleFactors

  ELEMENTAL FUNCTION Radius(r,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,z
    REAL :: Radius
    !------------------------------------------------------------------------!
    Radius = SQRT(r*r+z*z)
  END FUNCTION Radius

  ELEMENTAL SUBROUTINE PositionVector(r,z,rx,ry,rz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,z
    REAL, INTENT(OUT) :: rx,ry,rz
    !------------------------------------------------------------------------!
    rx = Radius(r,z)
    ry = 0.0
    rz = z
  END SUBROUTINE PositionVector

  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(r,phi,zz,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,phi,zz
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = r*COS(phi)
    y = r*SIN(phi)
    z = zz
  END SUBROUTINE Convert2Cartesian_coords

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(x,y,z,r,phi,zz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y,z
    REAL, INTENT(OUT) :: r,phi,zz
    !------------------------------------------------------------------------!
    r = SQRT(x*x+y*y)
    phi = ATAN2(y,x)
    IF (phi.LT.0.0) THEN
      phi = phi + 2.0*PI
    END IF
    zz = z
  END SUBROUTINE Convert2Curvilinear_coords

  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors(r,phi,z,vr,vphi,vzz,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)                        :: r,phi,z,vr,vphi,vzz
    REAL, INTENT(OUT)                       :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vr*COS(phi) - vphi*SIN(phi)
    vy = vr*SIN(phi) + vphi*COS(phi)
    vz = vzz
  END SUBROUTINE Convert2Cartesian_vectors


  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors(r,phi,z,vx,vy,vz,vr,vphi,vzz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)                        :: r,phi,z,vx,vy,vz
    REAL, INTENT(OUT)                       :: vr,vphi,vzz
    !------------------------------------------------------------------------!
    vr   = vx*COS(phi) + vy*SIN(phi)
    vphi = -vx*SIN(phi) + vy*COS(phi)
    vzz  = vz
  END SUBROUTINE Convert2Curvilinear_vectors


END MODULE geometry_cylindrical_mod
