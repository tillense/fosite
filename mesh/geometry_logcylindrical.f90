!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_logcylindrical.f90                                       #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
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
!! \brief defines properties of a 3D logcylindrical mesh
!----------------------------------------------------------------------------!
MODULE geometry_logcylindrical_mod
  USE geometry_base_mod
  USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS (geometry_cylindrical) :: geometry_logcylindrical
  CONTAINS
    PROCEDURE :: InitGeometry_logcylindrical
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
    PROCEDURE :: Convert2Cartesian_vectors_0
    PROCEDURE :: Convert2Curvilinear_vectors_0
    PROCEDURE :: Finalize
  END TYPE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "logcylindrical"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_logcylindrical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_logcylindrical(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(INOUT) :: this
    TYPE(DICT_TYP), POINTER                    :: config
    !------------------------------------------------------------------------!
    REAL :: gparam
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gparam", gparam, 1.0)

    CALL this%SetScale(gparam)
    CALL this%InitGeometry(LOGCYLINDRICAL,geometry_name,config)
    CALL this%SetAzimuthIndex(2)
  END SUBROUTINE InitGeometry_logcylindrical

  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,1),coords(:,2),coords(:,3),hx(:),hy(:),hz(:))
  END SUBROUTINE ScaleFactors_1

  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,:,1),coords(:,:,2),coords(:,:,3), &
                      hx(:,:),hy(:,:),hz(:,:))
  END SUBROUTINE ScaleFactors_2

  PURE SUBROUTINE ScaleFactors_3(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,:,:,1),coords(:,:,:,2),coords(:,:,:,3), &
                      hx(:,:,:),hy(:,:,:),hz(:,:,:))
  END SUBROUTINE ScaleFactors_3

  PURE SUBROUTINE ScaleFactors_4(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL ScaleFactors(this%geoparam(1),coords(:,:,:,:,1),coords(:,:,:,:,2),coords(:,:,:,:,3), &
                      hx(:,:,:,:),hy(:,:,:,:),hz(:,:,:,:))
  END SUBROUTINE ScaleFactors_4

  PURE SUBROUTINE Radius_1(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,1),coords(:,3))
  END SUBROUTINE Radius_1

  PURE SUBROUTINE Radius_2(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,:,1),coords(:,:,3))
  END SUBROUTINE Radius_2

  PURE SUBROUTINE Radius_3(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,:,:,1),coords(:,:,:,3))
  END SUBROUTINE Radius_3

  PURE SUBROUTINE Radius_4(this,coords,r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
    REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: r
    !------------------------------------------------------------------------!
    r = Radius(this%geoparam(1),coords(:,:,:,:,1),coords(:,:,:,:,3))
  END SUBROUTINE Radius_4

  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,1),coords(:,3),posvec(:,1), &
                        posvec(:,2),posvec(:,3))
  END SUBROUTINE PositionVector_1

  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,:,1),coords(:,:,3), &
                        posvec(:,:,1),posvec(:,:,2),posvec(:,:,3))
  END SUBROUTINE PositionVector_2

  PURE SUBROUTINE PositionVector_3(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,:,:,1),coords(:,:,:,3), &
                        posvec(:,:,:,1),posvec(:,:,:,2),posvec(:,:,:,3))
  END SUBROUTINE PositionVector_3

  PURE SUBROUTINE PositionVector_4(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN)  :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: coords
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: posvec
    !------------------------------------------------------------------------!
    CALL PositionVector(this%geoparam(1),coords(:,:,:,:,1),coords(:,:,:,:,3), &
                        posvec(:,:,:,:,1),posvec(:,:,:,:,2),posvec(:,:,:,:,3))
  END SUBROUTINE PositionVector_4

  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,1),curv(:,2),curv(:,3), &
                                  cart(:,1),cart(:,2),cart(:,3))
  END SUBROUTINE Convert2Cartesian_coords_1

  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,:,1),curv(:,:,2),curv(:,:,3), &
                                  cart(:,:,1),cart(:,:,2),cart(:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_2

  PURE SUBROUTINE Convert2Cartesian_coords_3(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                  cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_3

  PURE SUBROUTINE Convert2Cartesian_coords_4(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: cart
    !------------------------------------------------------------------------!
    CALL Convert2Cartesian_coords(this%geoparam(1),curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                  cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_4

  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,1),cart(:,2),cart(:,3), &
                                    curv(:,1),curv(:,2),curv(:,3))
  END SUBROUTINE Convert2Curvilinear_coords_1

  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,:,1),cart(:,:,2),cart(:,:,3), &
                                    curv(:,:,1),curv(:,:,2),curv(:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_2

  PURE SUBROUTINE Convert2Curvilinear_coords_3(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3), &
                                    curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_3

  PURE SUBROUTINE Convert2Curvilinear_coords_4(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: cart
    REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: curv
    !------------------------------------------------------------------------!
    CALL Convert2Curvilinear_coords(this%geoparam(1),cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3), &
                                    curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_4

  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors_0(this,xi,eta,phi,vxi,veta,vphi,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi,vxi,veta,vphi
    REAL, INTENT(OUT)                       :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vxi*COS(eta) - veta*SIN(eta)
    vy = vxi*SIN(eta) + veta*COS(eta)
    vz = vphi
  END SUBROUTINE Convert2Cartesian_vectors_0


  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vx,vy,vz,vxi,veta,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi,vx,vy,vz
    REAL, INTENT(OUT)                       :: vxi,veta,vphi
    !------------------------------------------------------------------------!
    vxi  = vx*COS(eta) + vy*SIN(eta)
    veta = -vx*SIN(eta) + vy*COS(eta)
    vphi = vz
  END SUBROUTINE Convert2Curvilinear_vectors_0

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logcylindrical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

  ELEMENTAL SUBROUTINE ScaleFactors(gp,logr,phi,z,hlogr,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,logr,phi,z
    REAL, INTENT(OUT) :: hlogr,hphi,hz
    !------------------------------------------------------------------------!
    hlogr = gp*EXP(logr)
    hphi  = hlogr
    hz    = 1.
  END SUBROUTINE ScaleFactors

  ELEMENTAL FUNCTION Radius(gp,logr,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,logr,z
    REAL :: Radius
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = gp*EXP(logr)
    Radius = SQRT(r*r+z*z)
  END FUNCTION Radius

  ELEMENTAL SUBROUTINE PositionVector(gp,logr,z,rx,ry,rz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,logr,z
    REAL, INTENT(OUT) :: rx,ry,rz
    !------------------------------------------------------------------------!
    rx = Radius(gp,logr,z)
    ry = 0.0
    rz = z
  END SUBROUTINE PositionVector

  ELEMENTAL SUBROUTINE Convert2Cartesian_coords(gp,logr,phi,zz,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,logr,phi,zz
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = gp*EXP(logr)
    x = r*COS(phi)
    y = r*SIN(phi)
    z = zz
  END SUBROUTINE Convert2Cartesian_coords

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords(gp,x,y,z,logr,phi,zz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: gp,x,y,z
    REAL, INTENT(OUT) :: logr,phi,zz
    !------------------------------------------------------------------------!
    REAL :: x_,y_
    !------------------------------------------------------------------------!
    x_ = x/gp
    y_ = y/gp
    logr = 0.5*LOG(x_*x_+y_*y_)
    phi = ATAN2(y_,x_)
    IF (phi.LT.0.0) THEN
      phi = phi + 2.0*PI
    END IF
    zz = z
  END SUBROUTINE Convert2Curvilinear_coords

END MODULE geometry_logcylindrical_mod
