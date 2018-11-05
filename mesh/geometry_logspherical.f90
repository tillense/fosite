!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_logspherical.f03                                         #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!! \brief defines properties of a 3D logspherical mesh
!----------------------------------------------------------------------------!
MODULE geometry_logspherical_mod
  USE geometry_base_mod
  USE geometry_spherical_mod, ONLY: geometry_spherical
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(geometry_spherical) :: geometry_logspherical
  CONTAINS
    PROCEDURE :: InitGeometry_logspherical
    PROCEDURE :: ScaleFactors_0
    PROCEDURE :: Radius_0
    PROCEDURE :: PositionVector_0
    PROCEDURE :: Convert2Cartesian_coords_0
    PROCEDURE :: Convert2Cartesian_vectors_0
    PROCEDURE :: Convert2Curvilinear_coords_0
    PROCEDURE :: Convert2Curvilinear_vectors_0
    PROCEDURE :: Finalize
  END TYPE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "logspherical"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_logspherical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_logspherical(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(INOUT) :: this
    TYPE(DICT_TYP),POINTER                   :: config
    !------------------------------------------------------------------------!
    REAL                                     :: gparam
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gparam", gparam, 1.0)

    CALL this%SetScale(gparam)
    CALL this%InitGeometry(LOGSPHERICAL,geometry_name,config)
  END SUBROUTINE InitGeometry_logspherical

  ELEMENTAL SUBROUTINE ScaleFactors_0(this,xi,eta,phi,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: xi,eta,phi
    REAL, INTENT(OUT)                     :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL this%Radius_0(xi,eta,phi,hx)
    hy = hx
    hz = hx * SIN(eta)
  END SUBROUTINE ScaleFactors_0

  ELEMENTAL SUBROUTINE Radius_0(this,xi,eta,phi,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: xi,eta,phi
    REAL, INTENT(OUT)                     :: radius
    !------------------------------------------------------------------------!
    radius = this%GetScale() * EXP(xi)
  END SUBROUTINE Radius_0

  ELEMENTAL SUBROUTINE PositionVector_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: xi,eta,phi
    REAL, INTENT(OUT)                     :: x,y,z
    !------------------------------------------------------------------------!
    CALL this%Radius_0(xi,eta,phi,x)
    y = 0.0
    z = 0.0
  END SUBROUTINE PositionVector_0

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: xi,eta,phi
    REAL, INTENT(OUT)                     :: x,y,z
    !------------------------------------------------------------------------!
    REAL :: r
    !------------------------------------------------------------------------!
    r = this%GetScale()*EXP(xi)
    x = r*SIN(eta)*COS(phi)
    y = r*SIN(eta)*SIN(phi)
    z = r*COS(eta)
  END SUBROUTINE Convert2Cartesian_coords_0

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords_0(this,x,y,z,xi,eta,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: x,y,z
    REAL, INTENT(OUT)                     :: xi,eta,phi
    !------------------------------------------------------------------------!
    REAL :: x1,y1,z1
    !------------------------------------------------------------------------!
    x1 = x/this%GetScale()
    y1 = y/this%GetScale()
    z1 = z/this%GetScale()
    xi  = SQRT(x1*x1+y1*y1+z1*z1)
    eta = ACOS(z1/xi)
    phi = ATAN2(y1,x1)
    IF(phi.LT.0.0) THEN
        phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Convert2Curvilinear_coords_0

  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors_0(this,xi,eta,phi,vxi,veta,vphi,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: xi,eta,phi,vxi,veta,vphi
    REAL, INTENT(OUT)                     :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vxi*SIN(eta)*COS(phi) + veta*COS(eta)*COS(phi) - vphi*SIN(phi)
    vy = vxi*SIN(eta)*SIN(phi) + veta*COS(eta)*SIN(phi) + vphi*COS(phi)
    vz = vxi*COS(eta) - veta*SIN(eta)
  END SUBROUTINE Convert2Cartesian_vectors_0

  !> Reference: \cite bronstein2008 , Tabelle 13.1
  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vx,vy,vz,vxi,veta,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(IN) :: this
    REAL, INTENT(IN)                      :: xi,eta,phi,vx,vy,vz
    REAL, INTENT(OUT)                     :: vxi,veta,vphi
    !------------------------------------------------------------------------!
    vxi  = vx*SIN(eta)*COS(phi) + vy*SIN(eta)*SIN(phi) + vz*COS(eta)
    veta = vx*COS(eta)*COS(phi) + vy*COS(eta)*SIN(phi) - vz*SIN(eta)
    vphi = -vx*SIN(phi) + vy*COS(phi)
  END SUBROUTINE Convert2Curvilinear_vectors_0

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_logspherical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

END MODULE geometry_logspherical_mod
