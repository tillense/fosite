!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_cartesian.f90                                            #
!#                                                                           #
!# Copyright (C) 2007-2010                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!! \author Jannes Klee
!!
!! \brief defines properties of a 3D cartesian mesh
!!
!! \extends geometry_common
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_cartesian_mod
  USE geometry_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS (geometry_base) ::  geometry_cartesian
  CONTAINS
    PROCEDURE :: InitGeometry_cartesian
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
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "cartesian"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_cartesian
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_cartesian(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(INOUT) :: this
    TYPE(DICT_TYP),POINTER                   :: config
    !------------------------------------------------------------------------!
    REAL                                     :: dz
    !------------------------------------------------------------------------!
    CALL this%InitGeometry(CARTESIAN,geometry_name,config)
    CALL GetAttr(config, "dz", dz, 1.0)
  END SUBROUTINE InitGeometry_cartesian

  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,1),coords(:,2),coords(:,3),hx(:),hy(:),hz(:))
  END SUBROUTINE ScaleFactors_1

  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,1),coords(:,:,2),coords(:,:,3), &
                      hx(:,:),hy(:,:),hz(:,:))
  END SUBROUTINE ScaleFactors_2

  PURE SUBROUTINE ScaleFactors_3(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,:,1),coords(:,:,:,2),coords(:,:,:,3), &
                      hx(:,:,:),hy(:,:,:),hz(:,:,:))
  END SUBROUTINE ScaleFactors_3

  PURE SUBROUTINE ScaleFactors_4(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(coords(:,:,:,:,1),coords(:,:,:,:,2),coords(:,:,:,:,3), &
                      hx(:,:,:,:),hy(:,:,:,:),hz(:,:,:,:))
  END SUBROUTINE ScaleFactors_4

  PURE SUBROUTINE Radius_1(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,1),coords(:,2),coords(:,3))
  END SUBROUTINE Radius_1

  PURE SUBROUTINE Radius_2(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,1),coords(:,:,2),coords(:,:,3))
  END SUBROUTINE Radius_2

  PURE SUBROUTINE Radius_3(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,:,1),coords(:,:,:,2),coords(:,:,:,3))
  END SUBROUTINE Radius_3

  PURE SUBROUTINE Radius_4(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(coords(:,:,:,:,1),coords(:,:,:,:,2),coords(:,:,:,:,3))
  END SUBROUTINE Radius_4

  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN)  :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    posvec = coords
  END SUBROUTINE PositionVector_1

  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    posvec = coords
  END SUBROUTINE PositionVector_2

  PURE SUBROUTINE PositionVector_3(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    posvec = coords
  END SUBROUTINE PositionVector_3

  PURE SUBROUTINE PositionVector_4(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    posvec = coords
  END SUBROUTINE PositionVector_4

  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    cart = curv
  END SUBROUTINE Convert2Cartesian_coords_1

  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    cart = curv
  END SUBROUTINE Convert2Cartesian_coords_2

  PURE SUBROUTINE Convert2Cartesian_coords_3(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    cart = curv
  END SUBROUTINE Convert2Cartesian_coords_3

  PURE SUBROUTINE Convert2Cartesian_coords_4(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    cart = curv
  END SUBROUTINE Convert2Cartesian_coords_4

  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    curv = cart
  END SUBROUTINE Convert2Curvilinear_coords_1

  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    curv = cart
  END SUBROUTINE Convert2Curvilinear_coords_2

  PURE SUBROUTINE Convert2Curvilinear_coords_3(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    curv = cart
  END SUBROUTINE Convert2Curvilinear_coords_3

  PURE SUBROUTINE Convert2Curvilinear_coords_4(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    curv = cart
  END SUBROUTINE Convert2Curvilinear_coords_4

  ! vector transformations
  PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    v_cart = v_curv
  END SUBROUTINE Convert2Cartesian_vectors_1

  PURE SUBROUTINE Convert2Cartesian_vectors_2(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    v_cart = v_curv
  END SUBROUTINE Convert2Cartesian_vectors_2

  PURE SUBROUTINE Convert2Cartesian_vectors_3(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    v_cart = v_curv
  END SUBROUTINE Convert2Cartesian_vectors_3

  PURE SUBROUTINE Convert2Cartesian_vectors_4(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: v_curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: v_cart
    !------------------------------------------------------------------------!
    v_cart = v_curv
  END SUBROUTINE Convert2Cartesian_vectors_4

  PURE SUBROUTINE Convert2Curvilinear_vectors_1(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    v_curv = v_cart
  END SUBROUTINE Convert2Curvilinear_vectors_1

  PURE SUBROUTINE Convert2Curvilinear_vectors_2(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    v_curv = v_cart
  END SUBROUTINE Convert2Curvilinear_vectors_2

  PURE SUBROUTINE Convert2Curvilinear_vectors_3(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    v_curv = v_cart
  END SUBROUTINE Convert2Curvilinear_vectors_3

  PURE SUBROUTINE Convert2Curvilinear_vectors_4(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: v_cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: v_curv
    !------------------------------------------------------------------------!
    v_curv = v_cart
  END SUBROUTINE Convert2Curvilinear_vectors_4

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

  ELEMENTAL SUBROUTINE ScaleFactors(x,y,z,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y,z
    REAL, INTENT(OUT) :: hx,hy,hz
    !------------------------------------------------------------------------!
    hx = 1.
    hy = 1.
    hz = 1.
  END SUBROUTINE ScaleFactors

  ELEMENTAL FUNCTION Radius(x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y,z
    REAL :: Radius
    !------------------------------------------------------------------------!
    Radius = SQRT(x*x+y*y+z*z)
  END FUNCTION Radius

END MODULE geometry_cartesian_mod
