!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: geometry_generic.f90                                              #
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
!! \author Manuel Jung
!!
!! \brief base class for geometrical properties
!!
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_base_mod
  USE logging_base_mod
  USE marray_cellscalar_mod
  USE marray_cellvector_mod
  USE common_dict
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: PI = 3.1415926535897932384626433832795028842
!  INTERFACE GetScale
!    MODULE PROCEDURE GetScale1, GetScale2
!  END INTERFACE
!  INTERFACE SetScale
!    MODULE PROCEDURE SetScale1, SetScale2, SetScale3
!  END INTERFACE

  TYPE, ABSTRACT, EXTENDS (logging_base) ::  geometry_base
    REAL,DIMENSION(3) :: geoparam                  !< geometry parameter
    INTEGER, PRIVATE  :: azimuthIndex = 0          !< index of azimuthal angle
  CONTAINS

!    PRIVATE
    PROCEDURE :: ScaleFactors_0
    PROCEDURE (ScaleFactors_1), DEFERRED :: ScaleFactors_1
    PROCEDURE (ScaleFactors_2), DEFERRED :: ScaleFactors_2
    PROCEDURE (ScaleFactors_3), DEFERRED :: ScaleFactors_3
    PROCEDURE (ScaleFactors_4), DEFERRED :: ScaleFactors_4
    PROCEDURE :: Radius_0
    PROCEDURE (Radius_1), DEFERRED :: Radius_1
    PROCEDURE (Radius_2), DEFERRED :: Radius_2
    PROCEDURE (Radius_3), DEFERRED :: Radius_3
    PROCEDURE (Radius_4), DEFERRED :: Radius_4
    PROCEDURE :: PositionVector_0
    PROCEDURE (PositionVector_1), DEFERRED :: PositionVector_1
    PROCEDURE (PositionVector_2), DEFERRED :: PositionVector_2
    PROCEDURE (PositionVector_3), DEFERRED :: PositionVector_3
    PROCEDURE (PositionVector_4), DEFERRED :: PositionVector_4
    PROCEDURE :: Convert2Cartesian_coords
    PROCEDURE (Convert2Cartesian_coords_1), DEFERRED :: Convert2Cartesian_coords_1
    PROCEDURE (Convert2Cartesian_coords_2), DEFERRED :: Convert2Cartesian_coords_2
    PROCEDURE (Convert2Cartesian_coords_3), DEFERRED :: Convert2Cartesian_coords_3
    PROCEDURE (Convert2Cartesian_coords_4), DEFERRED :: Convert2Cartesian_coords_4
    PROCEDURE :: Convert2Curvilinear_coords
    PROCEDURE (Convert2Curvilinear_coords_1), DEFERRED :: Convert2Curvilinear_coords_1
    PROCEDURE (Convert2Curvilinear_coords_2), DEFERRED :: Convert2Curvilinear_coords_2
    PROCEDURE (Convert2Curvilinear_coords_3), DEFERRED :: Convert2Curvilinear_coords_3
    PROCEDURE (Convert2Curvilinear_coords_4), DEFERRED :: Convert2Curvilinear_coords_4
    PROCEDURE :: Convert2Cartesian_vectors
    PROCEDURE (Convert2Cartesian_vectors_1), DEFERRED :: Convert2Cartesian_vectors_1
    PROCEDURE (Convert2Cartesian_vectors_2), DEFERRED :: Convert2Cartesian_vectors_2
    PROCEDURE (Convert2Cartesian_vectors_3), DEFERRED :: Convert2Cartesian_vectors_3
    PROCEDURE (Convert2Cartesian_vectors_4), DEFERRED :: Convert2Cartesian_vectors_4
    PROCEDURE :: Convert2Curvilinear_vectors
    PROCEDURE (Convert2Curvilinear_vectors_1), DEFERRED :: Convert2Curvilinear_vectors_1
    PROCEDURE (Convert2Curvilinear_vectors_2), DEFERRED :: Convert2Curvilinear_vectors_2
    PROCEDURE (Convert2Curvilinear_vectors_3), DEFERRED :: Convert2Curvilinear_vectors_3
    PROCEDURE (Convert2Curvilinear_vectors_4), DEFERRED :: Convert2Curvilinear_vectors_4
    PROCEDURE :: SetScale1
    PROCEDURE :: SetScale2
    PROCEDURE :: SetScale3
    PROCEDURE :: GetScale1
    PROCEDURE :: GetScale2
    PROCEDURE :: SetAzimuthIndex
    PROCEDURE :: GetAzimuthIndex

!    PUBLIC
    PROCEDURE, PUBLIC   :: InitGeometry
    PROCEDURE, PUBLIC   :: Finalize_base
    PROCEDURE (Finalize), DEFERRED    :: Finalize

    !> \public Compute scale factors for the given geometry
    GENERIC, PUBLIC :: ScaleFactors => ScaleFactors_0, ScaleFactors_1, ScaleFactors_2, &
                                       ScaleFactors_3, ScaleFactors_4
    !> \public Compute radial distances to the origin
    GENERIC, PUBLIC :: Radius => Radius_0, Radius_1, Radius_2, Radius_3, Radius_4
    !> \public Compute position vector components
    GENERIC, PUBLIC :: PositionVector => PositionVector_0, PositionVector_1, PositionVector_2, &
                                         PositionVector_3, PositionVector_4
    !> \public Convert curvilinear to cartesian coordinates
    GENERIC, PUBLIC :: Convert2Cartesian => Convert2Cartesian_coords, Convert2Cartesian_vectors, &
                                            Convert2Cartesian_coords_1, Convert2Cartesian_coords_2,&
                                            Convert2Cartesian_coords_3, Convert2Cartesian_coords_4, &
                                            Convert2Cartesian_vectors_1, Convert2Cartesian_vectors_2,&
                                            Convert2Cartesian_vectors_3, Convert2Cartesian_vectors_4
    !> \public Convert curvilinear to cartesian coordinates
    GENERIC, PUBLIC :: Convert2Curvilinear => Convert2Curvilinear_coords, Convert2Curvilinear_vectors, &
                                              Convert2Curvilinear_coords_1, Convert2Curvilinear_coords_2,&
                                              Convert2Curvilinear_coords_3, Convert2Curvilinear_coords_4,&
                                              Convert2Curvilinear_vectors_1, Convert2Curvilinear_vectors_2,&
                                              Convert2Curvilinear_vectors_3, Convert2Curvilinear_vectors_4
    GENERIC         :: SetScale => SetScale1, SetScale2, SetScale3
    GENERIC         :: GetScale => GetScale1, GetScale2
  END TYPE geometry_base

  ABSTRACT INTERFACE
    PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN) :: this
      REAL, DIMENSION(:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:),  INTENT(OUT) :: hx,hy,hz
    END SUBROUTINE
    PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)   :: this
      REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:,:),  INTENT(OUT) :: hx,hy,hz
    END SUBROUTINE
    PURE SUBROUTINE ScaleFactors_3(this,coords,hx,hy,hz)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)     :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:,:,:),  INTENT(OUT) :: hx,hy,hz
    END SUBROUTINE
    PURE SUBROUTINE ScaleFactors_4(this,coords,hx,hy,hz)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)       :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: hx,hy,hz
    END SUBROUTINE
    PURE SUBROUTINE Radius_1(this,coords,r)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN) :: this
      REAL, DIMENSION(:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:),  INTENT(OUT) :: r
    END SUBROUTINE
    PURE SUBROUTINE Radius_2(this,coords,r)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)   :: this
      REAL, DIMENSION(:,:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:,:),  INTENT(OUT) :: r
    END SUBROUTINE
    PURE SUBROUTINE Radius_3(this,coords,r)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)     :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:,:,:),  INTENT(OUT) :: r
    END SUBROUTINE
    PURE SUBROUTINE Radius_4(this,coords,r)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)       :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: coords
      REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: r
    END SUBROUTINE
    PURE SUBROUTINE PositionVector_1(this,coords,posvec)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:), INTENT(IN)  :: coords
      REAL, DIMENSION(:,:), INTENT(OUT) :: posvec
    END SUBROUTINE
    PURE SUBROUTINE PositionVector_2(this,coords,posvec)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)    :: this
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: coords
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: posvec
    END SUBROUTINE
    PURE SUBROUTINE PositionVector_3(this,coords,posvec)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)      :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: coords
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: posvec
    END SUBROUTINE
    PURE SUBROUTINE PositionVector_4(this,coords,posvec)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)        :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: coords
      REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: posvec
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:), INTENT(OUT) :: cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)    :: this
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_coords_3(this,curv,cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)      :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_coords_4(this,curv,cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)        :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:), INTENT(IN)  :: cart
      REAL, DIMENSION(:,:), INTENT(OUT) :: curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)    :: this
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: cart
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_coords_3(this,cart,curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)      :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: cart
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_coords_4(this,cart,curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)        :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: cart
      REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:), INTENT(IN)  :: v_curv
      REAL, DIMENSION(:,:), INTENT(OUT) :: v_cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_vectors_2(this,curv,v_curv,v_cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: v_curv
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: v_cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_vectors_3(this,curv,v_curv,v_cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: v_curv
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: v_cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Cartesian_vectors_4(this,curv,v_curv,v_cart)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: v_curv
      REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: v_cart
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_vectors_1(this,curv,v_cart,v_curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:), INTENT(IN)  :: v_cart
      REAL, DIMENSION(:,:), INTENT(OUT) :: v_curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_vectors_2(this,curv,v_cart,v_curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:), INTENT(IN)  :: v_cart
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: v_curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_vectors_3(this,curv,v_cart,v_curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: v_cart
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: v_curv
    END SUBROUTINE
    PURE SUBROUTINE Convert2Curvilinear_vectors_4(this,curv,v_cart,v_curv)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)  :: this
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: curv
      REAL, DIMENSION(:,:,:,:,:), INTENT(IN)  :: v_cart
      REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: v_curv
    END SUBROUTINE
    SUBROUTINE Finalize(this)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(INOUT) :: this
    END SUBROUTINE
  END INTERFACE

  !> \name Public Attributes
  !! #### geometries
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: CARTESIAN         = 1
!  INTEGER, PARAMETER :: SINHCARTESIAN     = 2
!  INTEGER, PARAMETER :: POLAR             = 20
!  INTEGER, PARAMETER :: LOGPOLAR          = 21
!  INTEGER, PARAMETER :: TANPOLAR          = 22
!  INTEGER, PARAMETER :: SINHPOLAR         = 23
!  INTEGER, PARAMETER :: SINHTANHPOLAR     = 24
!  INTEGER, PARAMETER :: POLYPOLAR         = 25
!  INTEGER, PARAMETER :: ELLIPTIC          = 27
!  INTEGER, PARAMETER :: BIPOLAR           = 29
  INTEGER, PARAMETER :: CYLINDRICAL       = 30
  INTEGER, PARAMETER :: LOGCYLINDRICAL    = 31
!  INTEGER, PARAMETER :: TANCYLINDRICAL    = 32
!  INTEGER, PARAMETER :: LNCOSHCYLINDRICAL = 33
  INTEGER, PARAMETER :: SPHERICAL         = 40
  INTEGER, PARAMETER :: LOGSPHERICAL      = 41
  INTEGER, PARAMETER :: SPHERICAL_PLANET  = 42
!  INTEGER, PARAMETER :: SINHSPHERICAL     = 43
!  INTEGER, PARAMETER :: OBLATE_SPHEROIDAL = 50
!  INTEGER, PARAMETER :: CHANNEL           = 60
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       geometry_base, &
       PI, &
       CARTESIAN, &
       CYLINDRICAL, LOGCYLINDRICAL, &
       SPHERICAL, LOGSPHERICAL, SPHERICAL_PLANET
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of generic geometry module
  SUBROUTINE InitGeometry(this,gnum,gname,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    TYPE(DICT_TYP),POINTER              :: config
    INTEGER                             :: gnum, gt !\todo{same variable}
    CHARACTER(LEN=*)                    :: gname
    REAL                                :: gs_def,gs_def2,gs_def3
    !------------------------------------------------------------------------!
    CHARACTER(LEN=8)                    :: gs_str
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: gnum,gname
    !------------------------------------------------------------------------!
    CALL this%logging_base%InitLogging(gnum,gname)

    CALL GetAttr(config, "geometry", gt)

    ! check if geometry parameters were given
    ! and set to defaults if not
    SELECT CASE(gnum)
    CASE(CARTESIAN,CYLINDRICAL,SPHERICAL,SPHERICAL_PLANET)
       ! do nothing (no parameters needed)
    CASE DEFAULT
       ! geometries with at least one parameter
       ! gparam defaults to 1 except for CHANNEL
!       IF (gt.EQ.CHANNEL) THEN
!          CALL GetAttr(config, "gparam", gs_def, 0.5)
!       ELSE
          CALL GetAttr(config, "gparam", gs_def, 1.0)
!       END IF
       ! geometries with three parameters
!       SELECT CASE(gt)
!       CASE(SINHPOLAR)
!          ! gp2,gp3 defaults to 0
!          CALL GetAttr(config, "gparam2", gs_def2, 0.0)
!          CALL GetAttr(config, "gparam3", gs_def3, 0.0)
!       CASE(SINHTANHPOLAR,LNCOSHCYLINDRICAL)
!          ! gp2,gp3 defaults to 1
!          CALL GetAttr(config, "gparam2", gs_def2, 1.0)
!          CALL GetAttr(config, "gparam3", gs_def3, 1.0)
!       END SELECT
    END SELECT

    ! print some information
    CALL this%Info(" GEOMETRY-> coordinates:       " // TRIM(this%GetName()))
    SELECT CASE(gt)
    CASE(LOGCYLINDRICAL,LOGSPHERICAL)
       WRITE (gs_str,'(ES8.1)') this%GetScale()
       CALL this%Info( "            geometry scale:   " // TRIM(gs_str))
!    CASE(LNCOSHCYLINDRICAL,SINHTANHPOLAR)
!       WRITE (gs_str,'(ES8.1)') this%GetScale()
!       CALL this%Info( "            geometry scale:    " // TRIM(gs_str))
!       WRITE (gs_str,'(ES8.1)') this%GetScale(2)
!       CALL this%Info("            geometry scale:    " // TRIM(gs_str))
!       WRITE (gs_str,'(ES8.1)') this%GetScale(3)
!       CALL this%Info("            geometry scale:    " // TRIM(gs_str))

    END SELECT
  END SUBROUTINE InitGeometry

  PURE FUNCTION GetScale1(this) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN) :: this
    REAL                             :: gp
    !------------------------------------------------------------------------!
    gp = this%geoparam(1)
  END FUNCTION GetScale1

  PURE FUNCTION GetScale2(this,i) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN) :: this
    INTEGER, INTENT(IN)              :: i
    REAL                             :: gp
    !------------------------------------------------------------------------!
    gp = this%geoparam(i)
  END FUNCTION GetScale2

  PURE SUBROUTINE SetScale1(this,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    REAL, INTENT(IN)                    :: gp
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
  END SUBROUTINE SetScale1

  PURE SUBROUTINE SetScale2(this,gp,gp2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    REAL, INTENT(IN)                    :: gp,gp2
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
    this%geoparam(2) = gp2
  END SUBROUTINE SetScale2

  PURE SUBROUTINE SetScale3(this,gp,gp2,gp3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    REAL, INTENT(IN)                    :: gp,gp2,gp3
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
    this%geoparam(2) = gp2
    this%geoparam(3) = gp3
  END SUBROUTINE SetScale3


  !> \public Compute scale factors
  !!
  !! Computes the scale factors of the given geometry.
  !! \f[
  !!   h_x = \left|\frac{\partial \vec{r}}{\partial x}\right|,\quad
  !!   h_y = \left|\frac{\partial \vec{r}}{\partial y}\right|,\quad
  !!   h_z = \left|\frac{\partial \vec{r}}{\partial z}\right|
  !! \f]
  !! ATTENTION: hx,hy,hz should actually be intent(out), but this breaks the
  !!            pointer assignments to hx%data1d,... in the mesh array
  PURE SUBROUTINE ScaleFactors_0(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)   :: this
    TYPE(marray_cellvector), INTENT(IN)    :: coords
    TYPE(marray_cellscalar), INTENT(INOUT) :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL this%ScaleFactors(coords%data2d(:,:),hx%data1d(:),hy%data1d(:),hz%data1d(:))
  END SUBROUTINE ScaleFactors_0

  !> \public Compute radial distances to the origin
  !!
  !! Computes the radial distances to the origin for the given
  !! coordinates depending on the geometry.
  !! ATTENTION: radius should actually be intent(out), but this breaks the
  !!            pointer assignment to radius%data1d in the mesh array
  PURE SUBROUTINE Radius_0(this,coords,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)     :: this
    TYPE(marray_cellvector), INTENT(IN)    :: coords
    TYPE(marray_cellscalar), INTENT(INOUT) :: radius
    !------------------------------------------------------------------------!
    CALL this%Radius(coords%data2d(:,:),radius%data1d(:))
  END SUBROUTINE Radius_0

  !> \public compute position vector components for all cell positions
  !!
  !! Computes the curvilinear position vector components with respect to
  !! the given geometry.
  !! ATTENTION: posvec should actually be intent(out), but this breaks the
  !!            pointer assignment to posvec%data1d in the mesh array
  PURE SUBROUTINE PositionVector_0(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN) :: this
    TYPE(marray_cellvector), INTENT(IN)    :: coords
    TYPE(marray_cellvector), INTENT(INOUT) :: posvec
    !------------------------------------------------------------------------!
    CALL this%PositionVector(coords%data2d(:,:),posvec%data2d(:,:))
  END SUBROUTINE PositionVector_0

  !> \public Convert curvilinear to cartesian coordinates
  !!
  !! ATTENTION: cart should actually be intent(out), but this breaks the
  !!            pointer assignment to cart%data1d in the mesh array
  PURE SUBROUTINE Convert2Cartesian_coords(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)       :: this  !> \param [in] this geometry class
    TYPE(marray_cellvector), INTENT(IN)    :: curv  !> \param [in] curv curvilinear coordinates
    TYPE(marray_cellvector), INTENT(INOUT) :: cart  !> \param [out] cart cartesian coordinates
    !------------------------------------------------------------------------!
    CALL this%Convert2Cartesian(curv%data2d(:,:),cart%data2d(:,:))
  END SUBROUTINE Convert2Cartesian_coords

  !> \public Convert cartesian to curvilinear coordinates
  !!
  !! ATTENTION: curv should actually be intent(out), but this breaks the
  !!            pointer assignment to curv%data1d in the mesh array
  PURE SUBROUTINE Convert2Curvilinear_coords(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)       :: this  !> \param [in] this geometry type
    TYPE(marray_cellvector), INTENT(IN)    :: cart  !> \param [in] cart cartesian coordinates
    TYPE(marray_cellvector), INTENT(INOUT) :: curv  !> \param [out] curv curvilinear coordinates
    !------------------------------------------------------------------------!
    CALL this%Convert2Curvilinear(cart%data2d(:,:),curv%data2d(:,:))
  END SUBROUTINE Convert2Curvilinear_coords

  !> \public Convert curvilinear vector components to cartesian vector components
  !!
  !! ATTENTION: v_cart should actually be intent(out), but this breaks the
  !!            pointer assignment to v_cart%data1d in the mesh array
  PURE SUBROUTINE Convert2Cartesian_vectors(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)       :: this   !> \param [in] this geometry class
    TYPE(marray_cellvector), INTENT(IN)    :: curv   !> \param [in] curv curvilinear coordinates
    TYPE(marray_cellvector), INTENT(IN)    :: v_curv !> \param [in] v_curv curvilinear vector comp.
    TYPE(marray_cellvector), INTENT(INOUT) :: v_cart !> \param [out] v_cart cartesian vector comp.
    !------------------------------------------------------------------------!
    CALL this%Convert2Cartesian(curv%data2d(:,:),v_curv%data2d(:,:),v_cart%data2d(:,:))
  END SUBROUTINE Convert2Cartesian_vectors


  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)       :: this   !> \param [in] this geometry class
    TYPE(marray_cellvector), INTENT(IN)    :: curv   !> \param [in] curv curvilinear coordinates
    TYPE(marray_cellvector), INTENT(IN)    :: v_cart !> \param [in] v_cart cartesian vector comp.
    TYPE(marray_cellvector), INTENT(INOUT) :: v_curv !> \param [out] v_curv curvilinear vector comp.
    !------------------------------------------------------------------------!
    CALL this%Convert2Curvilinear(curv%data2d(:,:),v_cart%data2d(:,:),v_curv%data2d(:,:))
  END SUBROUTINE Convert2Curvilinear_vectors

  !> \public sets the coordinate index of the azimuthal angle
  !!
  !! the default is 0, if there is no azimuthal angle, e.g. for cartesian coordinates
  PURE SUBROUTINE SetAzimuthIndex(this,idx)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    INTEGER, INTENT(IN)                 :: idx
    !------------------------------------------------------------------------!
    this%azimuthIndex = idx
  END SUBROUTINE SetAzimuthIndex

  !> \public returns the coordinate index of the azimuthal angle
  !!
  !! the default is 0, if there is no azimuthal angle, e.g. for cartesian coordinates
  PURE FUNCTION GetAzimuthIndex(this) RESULT(idx)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN) :: this
    INTEGER                          :: idx
    !------------------------------------------------------------------------!
    idx = this%azimuthIndex
  END FUNCTION GetAzimuthIndex

  !> \public Destructor of generic geometry module
  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
         CALL this%Error("CloseGeometry","not initialized")
    !CALL this%logging_base%Finalize()
  END SUBROUTINE Finalize_base

END MODULE geometry_base_mod
