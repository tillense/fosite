!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_cylindrical.f90                                          #
!#                                                                           #
!# Copyright (C) 2007                                                        #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jubin Lirawi     <jlirawi@astrophysik.uni-kiel.de>                        #
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
!> \author Jubin Lirawi
!!
!! \brief define properties of a 3D cylindrical mesh
!!
!! \extends geometry_cartesian
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
    PROCEDURE :: ScaleFactors_cylindrical
    PROCEDURE :: ScaleFactors_0
!    PROCEDURE :: Radius_cylindrical
    PROCEDURE :: Radius_0
    PROCEDURE :: PositionVector_0
    PROCEDURE :: Convert2Cartesian_cylindrical
    PROCEDURE :: Convert2Cartesian_coords_0
    PROCEDURE :: Convert2Cartesian_vectors_0
    PROCEDURE :: Convert2Curvilinear_cylindrical
    PROCEDURE :: Convert2Curvilinear_coords_0
    PROCEDURE :: Convert2Curvilinear_vectors_0
    FINAL     :: Finalize
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
  END SUBROUTINE InitGeometry_cylindrical


  ELEMENTAL SUBROUTINE ScaleFactors_cylindrical(this,r,hr,hphi,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: r
    REAL, INTENT(OUT)                       :: hr,hphi,hz
    !------------------------------------------------------------------------!
    hr   = 1.
    hphi = r
    hz   = 1.
  END SUBROUTINE ScaleFactors_cylindrical


  ELEMENTAL SUBROUTINE ScaleFactors_0(this,xi,eta,phi,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi
    REAL, INTENT(OUT)                       :: hx,hy,hz
    !------------------------------------------------------------------------!
    hx   = 1.
    hy = xi
    hz   = 1.
  END SUBROUTINE ScaleFactors_0


  ELEMENTAL SUBROUTINE Radius_0(this,xi,eta,phi,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi
    REAL, INTENT(OUT)                       :: radius
    !------------------------------------------------------------------------!
    radius = xi
  END SUBROUTINE Radius_0



  ELEMENTAL SUBROUTINE PositionVector_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi
    REAL, INTENT(OUT)                       :: x,y,z
    !------------------------------------------------------------------------!
    x = xi
    y = eta
    z = phi
  END SUBROUTINE PositionVector_0


  ! coordinate/vector transformation
  ! cylindrical -> cartesian
  ! this works for both, vectors and coords
  ELEMENTAL SUBROUTINE Convert2Cartesian_cylindrical(this,r,phi,z1,x,y,z2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: r,phi,z2
    REAL, INTENT(OUT)                       :: x,y,z1
    !------------------------------------------------------------------------!
    x = r * COS(phi)
    y = r * SIN(phi)
    z1 = z2
  END SUBROUTINE Convert2Cartesian_cylindrical


  ELEMENTAL SUBROUTINE Convert2Cartesian_coords_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi
    REAL, INTENT(OUT)                       :: x,y,z
    !------------------------------------------------------------------------!
    x = xi * COS(eta)
    y = xi * SIN(eta)
    z = phi
  END SUBROUTINE Convert2Cartesian_coords_0


  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors_0(this,xi,eta,phi,vxi,veta,vphi,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi,vxi,veta,vphi
    REAL, INTENT(OUT)                       :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vxi * xi * COS(eta)
    vy = vxi * xi * SIN(eta)
    vz = vphi * phi
  END SUBROUTINE Convert2Cartesian_vectors_0


  ! coordinate/vector transformation
  ! cartesian -> cylindrical
  ! this works for both, vectors and coords
  ELEMENTAL SUBROUTINE Convert2Curvilinear_cylindrical(this,x,y,z1,r,phi,z2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: x,y,z1
    REAL, INTENT(OUT)                       :: r,phi,z2
    !------------------------------------------------------------------------!
    r = SQRT(x*x + y*y)
    phi = ATAN2(y, x)
    z2 = z1
  END SUBROUTINE Convert2Curvilinear_cylindrical

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords_0(this,x,y,z,xi,eta,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: x,y,z
    REAL, INTENT(OUT)                       :: xi,eta,phi
    !------------------------------------------------------------------------!
    xi = SQRT(x*x + y*y)
    eta = ATAN2(y, x)
    phi = z
  END SUBROUTINE Convert2Curvilinear_coords_0

!  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vxi,phi1,z2)
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    CLASS(geometry_cylindrical), INTENT(IN) :: this
!    REAL, INTENT(IN)                        :: xi,eta,phi
!    REAL, INTENT(OUT)                       :: vxi,phi1,z2
!    !------------------------------------------------------------------------!
!    vxi  = SQRT(xi*xi + eta*eta)
!    phi1 = ATAN2(eta, xi)
!    z2   = phi
!  END SUBROUTINE Convert2Curvilinear_vectors_0

  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vx,vy,vz,vxi,veta,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cylindrical), INTENT(IN) :: this
    REAL, INTENT(IN)                        :: xi,eta,phi,vx,vy,vz
    REAL, INTENT(OUT)                       :: vxi,veta,vphi
    !------------------------------------------------------------------------!
    vxi  = vx*COS(eta) + vy*SIN(eta)
    veta = - vx*SIN(eta) + vy*COS(eta)
    vphi = vz
  END SUBROUTINE Convert2Curvilinear_vectors_0

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(geometry_cylindrical), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizeGeometry()
  END SUBROUTINE Finalize

END MODULE geometry_cylindrical_mod
