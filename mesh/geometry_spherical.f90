!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_polar.f90                                                #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>    #
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
!!
!! \brief define properties of a 2D polar mesh
!!
!!
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_spherical_mod
  USE geometry_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE, EXTENDS(geometry_base) :: geometry_spherical
  CONTAINS
    PROCEDURE :: InitGeometry_spherical
    PROCEDURE :: ScaleFactors_0
    PROCEDURE :: Radius_0
    PROCEDURE :: PositionVector_0
    PROCEDURE :: Convert2Cartesian_coords_0
    PROCEDURE :: Convert2Cartesian_vectors_0
    PROCEDURE :: Convert2Curvilinear_coords_0
    PROCEDURE :: Convert2Curvilinear_vectors_0
    FINAL     :: Finalize
  END TYPE
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "spherical"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_spherical
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_spherical(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(INOUT) :: this
    TYPE(DICT_TYP),POINTER            :: config
    !------------------------------------------------------------------------!
    REAL                              :: dz
    !------------------------------------------------------------------------!
    CALL this%InitGeometry(SPHERICAL,geometry_name,config)
    CALL GetAttr(config, "dz", dz, 1.0)
  END SUBROUTINE InitGeometry_spherical

  ELEMENTAL SUBROUTINE ScaleFactors_0(this,xi,eta,phi,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: hx,hy,hz
    !------------------------------------------------------------------------!
    hx = 1.
    hy = xi
    hz = xi * SIN(eta)
  END SUBROUTINE ScaleFactors_0

  ELEMENTAL SUBROUTINE Radius_0(this,xi,eta,phi,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: radius
    !------------------------------------------------------------------------!
    radius = xi
  END SUBROUTINE Radius_0

  ELEMENTAL SUBROUTINE PositionVector_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    CALL this%Radius_0(xi,eta,phi,x)
    y = 0.0
    z = 0.0
  END SUBROUTINE PositionVector_0

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = xi*SIN(eta)*COS(phi)
    y = xi*SIN(eta)*SIN(phi)
    z = xi*COS(eta)
  END SUBROUTINE Convert2Cartesian_coords_0

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords_0(this,x,y,z,xi,eta,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: x,y,z
    REAL, INTENT(OUT) :: xi,eta,phi
    !------------------------------------------------------------------------!
    xi  = SQRT(x*x+y*y+z*z)
    eta = ACOS(z/xi)
    phi = ATAN2(y,x)
    IF(phi.LT.0.0) THEN
        phi = phi + 2.0*PI
    END IF
  END SUBROUTINE Convert2Curvilinear_coords_0

  ! vector transformations
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors_0(this,xi,eta,phi,vxi,veta,vphi,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi,vxi,veta,vphi
    REAL, INTENT(OUT) :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = (vxi * COS(eta) - veta * SIN(eta)) * COS(phi) - vphi * SIN(phi)
    vy = (vxi * COS(eta) - veta * SIN(eta)) * SIN(phi) + vphi * COS(phi)
    vz =  vxi * SIN(eta) + veta * COS(eta)
  END SUBROUTINE Convert2Cartesian_vectors_0

  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vx,vy,vz,vxi,veta,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_spherical), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi,vx,vy,vz
    REAL, INTENT(OUT) :: vxi,veta,vphi
    !------------------------------------------------------------------------!
!    vxi  =  (vx * COS(eta) + vy * SIN(eta)) * COS(phi) + vz * SIN(phi)
!    veta = -(vx * COS(eta) + vy * SIN(eta)) * SIN(phi) + vz * COS(phi)
!    vphi = - vx * SIN(eta) + vy * COS(eta)
    vxi  =  (vx * COS(phi) + vy * SIN(phi)) * SIN(eta) + vz * COS(eta)
    veta =  (vx * COS(phi) + vy * SIN(phi)) * COS(eta) - vz * SIN(eta)
    vphi = - vx * SIN(phi) + vy * COS(phi)
  END SUBROUTINE Convert2Curvilinear_vectors_0

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(geometry_spherical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizeGeometry()
  END SUBROUTINE Finalize

END MODULE geometry_spherical_mod
