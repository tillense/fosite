!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_spherical_planet.f90                                     #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Manuel Jung                                                               #
!# Lars Bösch                                                                #
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
!! \author Lars Bösch
!! \author Jannes Klee
!!
!! \brief defines properties of a 3D spherical mesh
!!
!!
!!
!! \extends geometry_base
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_spherical_planet_mod
  USE geometry_base_mod
  USE marray_cellscalar_mod
  USE marray_cellvector_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(geometry_base) :: geometry_spherical_planet
  CONTAINS
    PROCEDURE :: InitGeometry_spherical_planet
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
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "spherical_planet"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_spherical_planet
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_spherical_planet(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(INOUT) :: this
    TYPE(DICT_TYP),POINTER                   :: config
    !------------------------------------------------------------------------!
    CALL this%InitGeometry(SPHERICAL_PLANET,geometry_name,config)
    CALL this%SetAzimuthIndex(3)
  END SUBROUTINE InitGeometry_spherical_planet

  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,1),coords(:,3),hx(:),hy(:),hz(:))
  END SUBROUTINE ScaleFactors_1

  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,1),coords(:,:,3),hx(:,:),hy(:,:),hz(:,:))
  END SUBROUTINE ScaleFactors_2

  PURE SUBROUTINE ScaleFactors_3(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,:,1),coords(:,:,:,3),hx(:,:,:), &
                      hy(:,:,:),hz(:,:,:))
  END SUBROUTINE ScaleFactors_3

  PURE SUBROUTINE ScaleFactors_4(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,:,:,1),coords(:,:,:,:,3),hx(:,:,:,:), &
                      hy(:,:,:,:),hz(:,:,:,:))
  END SUBROUTINE ScaleFactors_4

  PURE SUBROUTINE Radius_1(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,3))
  END SUBROUTINE Radius_1

  PURE SUBROUTINE Radius_2(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,3))
  END SUBROUTINE Radius_2

  PURE SUBROUTINE Radius_3(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,:,3))
  END SUBROUTINE Radius_3

  PURE SUBROUTINE Radius_4(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,:,:,3))
  END SUBROUTINE Radius_4

  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN)  :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,3),posvec(:,1),posvec(:,2),posvec(:,3))
  END SUBROUTINE PositionVector_1

  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,:,3),posvec(:,:,1),posvec(:,:,2),posvec(:,:,3))
  END SUBROUTINE PositionVector_2

  PURE SUBROUTINE PositionVector_3(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,:,:,3),posvec(:,:,:,1), &
                        posvec(:,:,:,2),posvec(:,:,:,3))
  END SUBROUTINE PositionVector_3

  PURE SUBROUTINE PositionVector_4(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(coords(:,:,:,:,3),posvec(:,:,:,:,1), &
                        posvec(:,:,:,:,2),posvec(:,:,:,:,3))
  END SUBROUTINE PositionVector_4

  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,1),curv(:,2),curv(:,3), &
                                  cart(:,1),cart(:,2),cart(:,3))
  END SUBROUTINE Convert2Cartesian_coords_1

  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,:,1),curv(:,:,2),curv(:,:,3), &
                                  cart(:,:,1),cart(:,:,2),cart(:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_2

  PURE SUBROUTINE Convert2Cartesian_coords_3(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                  cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_3

  PURE SUBROUTINE Convert2Cartesian_coords_4(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                  cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_4

  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,1),cart(:,2),cart(:,3), &
                                    curv(:,1),curv(:,2),curv(:,3))
  END SUBROUTINE Convert2Curvilinear_coords_1

  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,:,1),cart(:,:,2),cart(:,:,3), &
                                    curv(:,:,1),curv(:,:,2),curv(:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_2

  PURE SUBROUTINE Convert2Curvilinear_coords_3(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3), &
                                    curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_3

  PURE SUBROUTINE Convert2Curvilinear_coords_4(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3), &
                                    curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_4

  PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(IN) :: this
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
    CLASS(geometry_spherical_planet), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

  ELEMENTAL SUBROUTINE ScaleFactors(theta,r,htheta,hphi,hr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta
    REAL, INTENT(OUT) :: hr,htheta,hphi
    !------------------------------------------------------------------------!
    htheta = r
    hphi   = r * SIN(theta)
    hr     = 1.
  END SUBROUTINE ScaleFactors

  ELEMENTAL FUNCTION Radius(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL :: Radius
    !------------------------------------------------------------------------!
    Radius = r
  END FUNCTION Radius

  ELEMENTAL SUBROUTINE PositionVector(r,rx,ry,rz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r
    REAL, INTENT(OUT) :: rx,ry,rz
    !------------------------------------------------------------------------!
    rx = 0.0
    ry = 0.0
    rz = r
  END SUBROUTINE PositionVector

  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(theta,phi,r,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: r,theta,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = r*SIN(theta)*COS(phi)
    y = r*SIN(theta)*SIN(phi)
    z = r*COS(theta)
  END SUBROUTINE Convert2Cartesian_coords

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(x,y,z,theta,phi,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y,z
    REAL, INTENT(OUT) :: r,theta,phi
    !------------------------------------------------------------------------!
    r = SQRT(x*x+y*y+z*z)
    theta = ACOS(z/r)
    phi = ATAN2(y,x)
    IF (phi.LT.0.0) THEN
      phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Convert2Curvilinear_coords

  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors(theta,phi,r,vtheta,vphi,vr,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)                      :: r,theta,phi,vr,vtheta,vphi
    REAL, INTENT(OUT)                     :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vr*SIN(theta)*COS(phi) + vtheta*COS(theta)*COS(phi) - vphi*SIN(phi)
    vy = vr*SIN(theta)*SIN(phi) + vtheta*COS(theta)*SIN(phi) + vphi*COS(phi)
    vz = vr*COS(theta) - vtheta*SIN(theta)
  END SUBROUTINE Convert2Cartesian_vectors

  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors(theta,phi,r,vx,vy,vz,vtheta,vphi,vr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)                      :: r,theta,phi,vx,vy,vz
    REAL, INTENT(OUT)                     :: vr,vtheta,vphi
    !------------------------------------------------------------------------!
    vr     = vx*SIN(theta)*COS(phi) + vy*SIN(theta)*SIN(phi) + vz*COS(theta)
    vtheta = vx*COS(theta)*COS(phi) + vy*COS(theta)*SIN(phi) - vz*SIN(theta)
    vphi   = -vx*SIN(phi) + vy*COS(phi)
  END SUBROUTINE Convert2Curvilinear_vectors

END MODULE geometry_spherical_planet_mod
