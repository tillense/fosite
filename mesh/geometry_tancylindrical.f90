!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_tancylindrical.f90                                       #
!#                                                                           #
!# Copyright (C) 2009 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!! \author Lars Boesch
!!
!! \brief define properties of a 3D tancylindrical mesh
!!
!! dimensionless vertical coordinate zeta (with -pi/2 < zeta < +pi/2)
!! according to:
!!   x = r
!!   z = z0 * tan(zeta)
!!
!! \extends geometry_cylindrical
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_tancylindrical_mod
  USE geometry_base_mod
  USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS (geometry_cylindrical) :: geometry_tancylindrical
  CONTAINS
    PROCEDURE :: InitGeometry_tancylindrical
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
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "tancylindrical"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_tancylindrical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_tancylindrical(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(INOUT) :: this
    TYPE(DICT_TYP), POINTER                      :: config
    !------------------------------------------------------------------------!
    REAL :: gparam
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gparam", gparam, 1.0)

    CALL this%SetScale(gparam)
    CALL this%InitGeometry(TANCYLINDRICAL,geometry_name,config)
    CALL this%SetAzimuthIndex(3)
  END SUBROUTINE InitGeometry_tancylindrical


  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,1),coords(:,2),hx(:),hy(:),hz(:))
  END SUBROUTINE ScaleFactors_1

  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,:,1),coords(:,:,2),hx(:,:),hy(:,:),hz(:,:))
  END SUBROUTINE ScaleFactors_2

  PURE SUBROUTINE ScaleFactors_3(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,:,:,1),coords(:,:,:,2),hx(:,:,:),hy(:,:,:),hz(:,:,:))
  END SUBROUTINE ScaleFactors_3

  PURE SUBROUTINE ScaleFactors_4(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,:,:,:,1),coords(:,:,:,:,2),hx(:,:,:,:),hy(:,:,:,:),hz(:,:,:,:))
  END SUBROUTINE ScaleFactors_4

  PURE SUBROUTINE Radius_1(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,1),coords(:,2))
  END SUBROUTINE Radius_1

  PURE SUBROUTINE Radius_2(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,:,1),coords(:,:,2))
  END SUBROUTINE Radius_2

  PURE SUBROUTINE Radius_3(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,:,:,1),coords(:,:,:,2))
  END SUBROUTINE Radius_3

  PURE SUBROUTINE Radius_4(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,:,:,:,1),coords(:,:,:,:,2))
  END SUBROUTINE Radius_4

  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,1),coords(:,2),posvec(:,1),posvec(:,2),posvec(:,3))
  END SUBROUTINE PositionVector_1

  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,:,1),coords(:,:,2),posvec(:,:,1), &
                        posvec(:,:,2),posvec(:,:,3))
  END SUBROUTINE PositionVector_2

  PURE SUBROUTINE PositionVector_3(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,:,:,1),coords(:,:,:,2),posvec(:,:,:,1), &
                        posvec(:,:,:,2),posvec(:,:,:,3))
  END SUBROUTINE PositionVector_3

  PURE SUBROUTINE PositionVector_4(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,:,:,:,1),coords(:,:,:,:,2),posvec(:,:,:,:,1), &
                        posvec(:,:,:,:,2),posvec(:,:,:,:,3))
  END SUBROUTINE PositionVector_4

  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,1),curv(:,2),curv(:,3), &
                                  cart(:,1),cart(:,2),cart(:,3))
  END SUBROUTINE Convert2Cartesian_coords_1

  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,:,1),curv(:,:,2),curv(:,:,3), &
                                  cart(:,:,1),cart(:,:,2),cart(:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_2

  PURE SUBROUTINE Convert2Cartesian_coords_3(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                  cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_3

  PURE SUBROUTINE Convert2Cartesian_coords_4(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                  cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_4

  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,1),cart(:,2),cart(:,3), &
                                    curv(:,1),curv(:,2),curv(:,3))
  END SUBROUTINE Convert2Curvilinear_coords_1

  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,:,1),cart(:,:,2),cart(:,:,3), &
                                    curv(:,:,1),curv(:,:,2),curv(:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_2

  PURE SUBROUTINE Convert2Curvilinear_coords_3(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3), &
                                    curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_3

  PURE SUBROUTINE Convert2Curvilinear_coords_4(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3), &
                                    curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_4

  PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(IN) :: this
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
    CLASS(geometry_tancylindrical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize



  ELEMENTAL SUBROUTINE ScaleFactors(gp,zeta,r,hzeta,hr,hphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL, INTENT(OUT) :: hzeta,hr,hphi
    !------------------------------------------------------------------------!
    hzeta = gp/COS(zeta)**2
    hr = 1.
    hphi = r
  END SUBROUTINE ScaleFactors


  ELEMENTAL FUNCTION Radius(gp,zeta,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL              :: Radius
    !------------------------------------------------------------------------!
    REAL :: z
    !------------------------------------------------------------------------!
    z = gp*TAN(zeta)
    Radius = SQRT(z*z+r*r)
  END FUNCTION Radius


  ELEMENTAL SUBROUTINE PositionVector(gp,zeta,r,rx,ry,rz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r
    REAL, INTENT(OUT) :: rx,ry,rz
    !------------------------------------------------------------------------!
    rx = gp*TAN(zeta)
    ry = r
    rz = 0.0
  END SUBROUTINE PositionVector


  ! coordinate transformations
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(gp,zeta,r,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,zeta,r,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = r*COS(phi)
    y = r*SIN(phi)
    z = gp*TAN(zeta)
  END SUBROUTINE Convert2Cartesian_coords

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(gp,x,y,z,zeta,r,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y,z
    REAL, INTENT(OUT) :: zeta,r,phi
    !------------------------------------------------------------------------!
    zeta = ATAN(z/gp)
    r = SQRT(x*x+y*y)
    phi = ATAN2(y,x)
    IF (phi.LT.0.0) THEN
      phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Convert2Curvilinear_coords


!  ! vector transformations
!  ELEMENTAL SUBROUTINE Tancyl2Cartesian_vectors(vzeta,vr,vx,vy)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: vzeta,vr
!    REAL, INTENT(OUT) :: vx,vy
!    !------------------------------------------------------------------------!
!    vx = vr
!    vy = vzeta
!  END SUBROUTINE Tancyl2Cartesian_vectors
!
!  ! cartesian -> tancylindrical
!  ELEMENTAL SUBROUTINE Cartesian2Tancyl_vectors(vx,vy,vzeta,vr)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    REAL, INTENT(IN)  :: vx,vy
!    REAL, INTENT(OUT) :: vzeta,vr
!    !------------------------------------------------------------------------!
!    vzeta = vy
!    vr = vx
!  END SUBROUTINE Cartesian2Tancyl_vectors

  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors(zeta,r,phi,vzeta,vr,vphi,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)                        :: zeta,r,phi,vzeta,vr,vphi
    REAL, INTENT(OUT)                       :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vr*COS(phi) - vphi*SIN(phi)
    vy = vr*SIN(phi) + vphi*COS(phi)
    vz = vzeta
  END SUBROUTINE Convert2Cartesian_vectors


  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors(zeta,r,phi,vx,vy,vz,vzeta,vr,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)                        :: zeta,r,phi,vx,vy,vz
    REAL, INTENT(OUT)                       :: vzeta,vr,vphi
    !------------------------------------------------------------------------!
    vzeta  = vz
    vr   = vx*COS(phi) + vy*SIN(phi)
    vphi = -vx*SIN(phi) + vy*COS(phi)
  END SUBROUTINE Convert2Curvilinear_vectors


END MODULE geometry_tancylindrical_mod
