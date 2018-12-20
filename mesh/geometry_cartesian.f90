!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
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
!! \brief define properties of a 2D cartesian mesh
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
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "cartesian"
  !--------------------------------------------------------------------------!
  PUBLIC :: geometry_cartesian
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_cartesian(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(INOUT) :: this
    TYPE(DICT_TYP),POINTER            :: config
    !------------------------------------------------------------------------!
    REAL                              :: dz
    !------------------------------------------------------------------------!
    CALL this%InitGeometry(CARTESIAN,geometry_name,config)
    CALL GetAttr(config, "dz", dz, 1.0)
  END SUBROUTINE InitGeometry_cartesian

  ELEMENTAL SUBROUTINE ScaleFactors_0(this,xi,eta,phi,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: hx,hy,hz
    !------------------------------------------------------------------------!
    ! scale factors are unity
    hx = 1.
    hy = 1.
    hz = 1.
  END SUBROUTINE ScaleFactors_0

  ELEMENTAL SUBROUTINE Radius_0(this,xi,eta,phi,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: radius
    !------------------------------------------------------------------------!
    radius = SQRT(xi*xi+eta*eta+phi*phi)
  END SUBROUTINE Radius_0

  ELEMENTAL SUBROUTINE PositionVector_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = xi
    y = eta
    z = phi
  END SUBROUTINE PositionVector_0

  ! coordinate transformations
  ELEMENTAL SUBROUTINE Convert2Cartesian_coords_0(this,xi,eta,phi,x,y,z)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi
    REAL, INTENT(OUT) :: x,y,z
    !------------------------------------------------------------------------!
    x = xi
    y = eta
    z = phi
  END SUBROUTINE Convert2Cartesian_coords_0

  ELEMENTAL SUBROUTINE Convert2Curvilinear_coords_0(this,x,y,z,xi,eta,phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: x,y,z
    REAL, INTENT(OUT) :: xi,eta,phi
    !------------------------------------------------------------------------!
    xi  = x
    eta = y
    phi = z
  END SUBROUTINE Convert2Curvilinear_coords_0

  ! vector transformations
  ELEMENTAL SUBROUTINE Convert2Cartesian_vectors_0(this,xi,eta,phi,vxi,veta,vphi,vx,vy,vz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi,vxi,veta,vphi
    REAL, INTENT(OUT) :: vx,vy,vz
    !------------------------------------------------------------------------!
    vx = vxi  * xi
    vy = veta * eta
    vz = vphi * phi
  END SUBROUTINE Convert2Cartesian_vectors_0

  ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vx,vy,vz,vxi,veta,vphi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_cartesian), INTENT(IN) :: this
    REAL, INTENT(IN)  :: xi,eta,phi,vx,vy,vz
    REAL, INTENT(OUT) :: vxi,veta,vphi
    !------------------------------------------------------------------------!
    vxi  = vx * xi
    veta = vy * eta
    vphi = vz * phi
  END SUBROUTINE Convert2Curvilinear_vectors_0

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(geometry_cartesian), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL this%FinalizeGeometry()
  END SUBROUTINE Finalize
  
END MODULE geometry_cartesian_mod
