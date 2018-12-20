!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_generic.f90                                              #
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
!! \author Manuel Jung
!!
!! \brief generic module for geometrical properties
!!
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_base_mod
  USE logging_base_mod
  USE array, ONLY : MArrayS_TYP, MArrayV_TYP !\todo{new type - not yet implemented as class}
  USE common_dict
!  USE geometry_cartesian, InitGeometry_common => InitGeometry, &
!       CloseGeometry_common => CloseGeometry
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
  CONTAINS

!    PRIVATE
    PROCEDURE :: ScaleFactors_all
    PROCEDURE (ScaleFactors_0), DEFERRED :: ScaleFactors_0
    PROCEDURE :: ScaleFactors_1
    PROCEDURE :: ScaleFactors_2
    PROCEDURE :: Radius_all
    PROCEDURE (Radius_0), DEFERRED :: Radius_0
    PROCEDURE :: Radius_1
    PROCEDURE :: Radius_2
    PROCEDURE :: PositionVector_all
    PROCEDURE (PositionVector_0), DEFERRED :: PositionVector_0
    PROCEDURE :: PositionVector_1
    PROCEDURE :: PositionVector_2
    PROCEDURE :: Convert2Cartesian_all
    PROCEDURE (Convert2Cartesian_coords_0), DEFERRED :: Convert2Cartesian_coords_0
    PROCEDURE :: Convert2Cartesian_coords_1
    PROCEDURE :: Convert2Cartesian_coords_2
    PROCEDURE (Convert2Cartesian_vectors_0), DEFERRED :: Convert2Cartesian_vectors_0
    PROCEDURE :: Convert2Cartesian_vectors_1
    PROCEDURE :: Convert2Cartesian_vectors_2
    PROCEDURE :: Convert2Cartesian_vectors_3
    PROCEDURE :: Convert2Curvilinear_all
    PROCEDURE (Convert2Curvilinear_coords_0), DEFERRED :: Convert2Curvilinear_coords_0
    PROCEDURE :: Convert2Curvilinear_coords_1
    PROCEDURE :: Convert2Curvilinear_coords_2
    PROCEDURE (Convert2Curvilinear_vectors_0), DEFERRED :: Convert2Curvilinear_vectors_0
    PROCEDURE :: Convert2Curvilinear_vectors_1
    PROCEDURE :: Convert2Curvilinear_vectors_2
    PROCEDURE :: Convert2Curvilinear_vectors_3
    PROCEDURE :: SetScale1
    PROCEDURE :: SetScale2
    PROCEDURE :: SetScale3
    PROCEDURE :: GetScale1
    PROCEDURE :: GetScale2

!    PUBLIC
    PROCEDURE, PUBLIC   :: InitGeometry
    PROCEDURE, PUBLIC   :: FinalizeGeometry

    GENERIC, PUBLIC     :: ScaleFactors => ScaleFactors_all, ScaleFactors_0, ScaleFactors_1, ScaleFactors_2
    !> \public Convert curvilinear to cartesian coordinates
    GENERIC, PUBLIC     :: Convert2Cartesian => Convert2Cartesian_all, &
                                                Convert2Cartesian_coords_0, Convert2Cartesian_coords_1,&
                                                Convert2Cartesian_coords_2, &
                                                Convert2Cartesian_vectors_0, Convert2Cartesian_vectors_1,&
                                                Convert2Cartesian_vectors_2, Convert2Cartesian_vectors_3
    !> \public Convert curvilinear to cartesian coordinates
    GENERIC, PUBLIC     :: Convert2Curvilinear => Convert2Curvilinear_all, &
                                                  Convert2Curvilinear_coords_0, Convert2Curvilinear_coords_1,&
                                                  Convert2Curvilinear_coords_2, &
                                                  Convert2Curvilinear_vectors_0, Convert2Curvilinear_vectors_1,&
                                                  Convert2Curvilinear_vectors_2, Convert2Curvilinear_vectors_3
    !> \public Compute radial distances to the origin
    GENERIC, PUBLIC     :: Radius => Radius_all, Radius_0, Radius_1, Radius_2
    !> \public Compute position vector components
    GENERIC, PUBLIC     :: PositionVector => PositionVector_all, PositionVector_0, PositionVector_1, PositionVector_2
    GENERIC             :: SetScale => SetScale1, SetScale2, SetScale3
    GENERIC             :: GetScale => GetScale1, GetScale2
  END TYPE geometry_base

  ABSTRACT INTERFACE
    PURE ELEMENTAL SUBROUTINE ScaleFactors_0(this,xi,eta,phi,hx,hy,hz)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)                         :: xi,eta,phi
      REAL, INTENT(OUT)                        :: hx,hy,hz
    END SUBROUTINE
    PURE ELEMENTAL SUBROUTINE Radius_0(this,xi,eta,phi,radius)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)                         :: xi,eta,phi
      REAL, INTENT(OUT)                        :: radius
    END SUBROUTINE
    PURE ELEMENTAL SUBROUTINE PositionVector_0(this,xi,eta,phi,x,y,z)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)                         :: xi,eta,phi
      REAL, INTENT(OUT)                        :: x,y,z
    END SUBROUTINE
    PURE ELEMENTAL SUBROUTINE Convert2Cartesian_coords_0(this,xi,eta,phi,x,y,z)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)  :: xi,eta,phi
      REAL, INTENT(OUT) :: x,y,z
    END SUBROUTINE
    PURE ELEMENTAL SUBROUTINE Convert2Cartesian_vectors_0(this,xi,eta,phi,vxi,veta,vphi,vx,vy,vz)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)                         :: xi,eta,phi,vxi,veta,vphi
      REAL, INTENT(OUT)                        :: vx,vy,vz
    END SUBROUTINE
    PURE ELEMENTAL SUBROUTINE Convert2Curvilinear_coords_0(this,x,y,z,xi,eta,phi)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)                         :: x,y,z
      REAL, INTENT(OUT)                        :: xi,eta,phi
    END SUBROUTINE
    PURE ELEMENTAL SUBROUTINE Convert2Curvilinear_vectors_0(this,xi,eta,phi,vx,vy,vz,vxi,veta,vphi)
      IMPORT geometry_base
      IMPLICIT NONE
      CLASS(geometry_base), INTENT(IN)         :: this
      REAL, INTENT(IN)                         :: xi,eta,phi,vx,vy,vz
      REAL, INTENT(OUT)                        :: vxi,veta,vphi
    END SUBROUTINE
  END INTERFACE

  !> \name Public Attributes
  !! #### geometries
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: CARTESIAN         = 1
  INTEGER, PARAMETER :: SINHCARTESIAN     = 2
  INTEGER, PARAMETER :: POLAR             = 20
  INTEGER, PARAMETER :: LOGPOLAR          = 21
  INTEGER, PARAMETER :: TANPOLAR          = 22
  INTEGER, PARAMETER :: SINHPOLAR         = 23
  INTEGER, PARAMETER :: SINHTANHPOLAR     = 24
  INTEGER, PARAMETER :: POLYPOLAR         = 25
  INTEGER, PARAMETER :: ELLIPTIC          = 27
  INTEGER, PARAMETER :: BIPOLAR           = 29
  INTEGER, PARAMETER :: CYLINDRICAL       = 30
  INTEGER, PARAMETER :: TANCYLINDRICAL    = 31
  INTEGER, PARAMETER :: LNCOSHCYLINDRICAL = 32
  INTEGER, PARAMETER :: SPHERICAL         = 40
  INTEGER, PARAMETER :: SINHSPHERICAL     = 41
  INTEGER, PARAMETER :: BIANGLESPHERICAL  = 42
  INTEGER, PARAMETER :: OBLATE_SPHEROIDAL = 50
  INTEGER, PARAMETER :: CHANNEL           = 60
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       geometry_base, &
       PI, &
       CARTESIAN, POLAR, LOGPOLAR, TANPOLAR, SINHPOLAR, SINHTANHPOLAR, BIPOLAR, &
       CYLINDRICAL, TANCYLINDRICAL, LNCOSHCYLINDRICAL, SPHERICAL, SINHSPHERICAL, &
       BIANGLESPHERICAL, OBLATE_SPHEROIDAL, CHANNEL, POLYPOLAR, ELLIPTIC, &
       SINHCARTESIAN
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of generic geometry module
  SUBROUTINE InitGeometry(this,gnum,gname,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT)    :: this
    TYPE(DICT_TYP),POINTER                 :: config
    INTEGER                                :: gnum, gt !\todo{same variable}
    CHARACTER(LEN=*)                       :: gname
    REAL                                   :: gs_def,gs_def2,gs_def3, dr
    !------------------------------------------------------------------------!
    CHARACTER(LEN=8)                       :: gs_str
    REAL                                   :: dz
    !------------------------------------------------------------------------!
    INTENT(IN)                             :: gnum,gname
    !------------------------------------------------------------------------!

    !\todo{Here initcommon is commentet!}
    CALL this%logging_base%InitLogging(gnum,gname)
    ! set geometry parameter to default value
    CALL this%SetScale(1.0,1.0,1.0) !\todo{what is this for?}

    CALL GetAttr(config, "geometry", gt)


    ! check if geometry parameters were given
    ! and set to defaults if not
    SELECT CASE(gnum)
    CASE(CARTESIAN,POLAR,CYLINDRICAL,SPHERICAL)
       ! do nothing (no parameters needed)
    CASE DEFAULT
       ! geometries with at least one parameter
       ! gparam defaults to 1 except for CHANNEL
       IF (gt.EQ.CHANNEL) THEN
          CALL GetAttr(config, "gparam", gs_def, 0.5)
       ELSE
          CALL GetAttr(config, "gparam", gs_def, 1.0)
       END IF
       ! geometries with three parameters
       SELECT CASE(gt)
       CASE(SINHPOLAR)
          ! gp2,gp3 defaults to 0
          CALL GetAttr(config, "gparam2", gs_def2, 0.0)
          CALL GetAttr(config, "gparam3", gs_def3, 0.0)
       CASE(SINHTANHPOLAR,LNCOSHCYLINDRICAL)
          ! gp2,gp3 defaults to 1
          CALL GetAttr(config, "gparam2", gs_def2, 1.0)
          CALL GetAttr(config, "gparam3", gs_def3, 1.0)
       END SELECT
    END SELECT


    !SELECT CASE(gnum)
    !CASE(CARTESIAN)
    !   CALL InitGeometry_cartesian(this,gt)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(SINHCARTESIAN)
    !   CALL InitGeometry_sinhcartesian(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(POLAR)
    !   CALL InitGeometry_polar(this,gt)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(LOGPOLAR)
    !   CALL InitGeometry_logpolar(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(TANPOLAR)
    !   CALL InitGeometry_tanpolar(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(SINHPOLAR)
    !   CALL InitGeometry_sinhpolar(this,gt,gs_def,gs_def2,gs_def3)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(SINHTANHPOLAR)
    !   CALL InitGeometry_sinhtanh(this,gt,gs_def,gs_def2,gs_def3)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(POLYPOLAR)
    !   CALL InitGeometry_polypolar(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(ELLIPTIC)
    !   CALL InitGeometry_elliptic(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(BIPOLAR)
    !   CALL InitGeometry_bipolar(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(CYLINDRICAL)
    !   CALL InitGeometry_cylindrical(this,gt)
    !   CALL GetAttr(config, "dz", dz, 2.0*PI)
    !CASE(TANCYLINDRICAL)
    !   CALL InitGeometry_tancyl(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 2.0*PI)
    !CASE(LNCOSHCYLINDRICAL)
    !   CALL InitGeometry_lncoshcyl(this,gt,gs_def,gs_def2,gs_def3)
    !   CALL GetAttr(config, "dz", dz, 2.0*PI)
    !CASE(SPHERICAL)
    !   CALL InitGeometry_spherical(this,gt)
    !   CALL GetAttr(config, "dz", dz, 2.0*PI)
    !CASE(SINHSPHERICAL)
    !   CALL InitGeometry_sinhspher(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 2.0*PI)
    !CASE(BIANGLESPHERICAL)
    !   CALL InitGeometry_bianglespher(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE(OBLATE_SPHEROIDAL)
    !   CALL InitGeometry_oblatespher(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 2.0*PI)
    !CASE(CHANNEL)
    !   CALL InitGeometry_channel(this,gt,gs_def)
    !   CALL GetAttr(config, "dz", dz, 1.0)
    !CASE DEFAULT
    !   CALL Error(this,"InitGeometry",  "Unknown geometry.")
    !END SELECT

    ! print some information
    CALL this%Info(" GEOMETRY-> coordinates:       " // TRIM(this%GetName()))
    SELECT CASE(gt)
    CASE(LOGPOLAR,TANPOLAR,SINHPOLAR,TANCYLINDRICAL,SINHSPHERICAL,&
         BIANGLESPHERICAL,OBLATE_SPHEROIDAL,CHANNEL,POLYPOLAR,ELLIPTIC,&
         SINHCARTESIAN)
       WRITE (gs_str,'(ES8.1)') this%GetScale()
       CALL this%Info( "            geometry scale:    " // TRIM(gs_str))
    CASE(LNCOSHCYLINDRICAL,SINHTANHPOLAR)
       WRITE (gs_str,'(ES8.1)') this%GetScale()
       CALL this%Info( "            geometry scale:    " // TRIM(gs_str))
       WRITE (gs_str,'(ES8.1)') this%GetScale(2)
       CALL this%Info("            geometry scale:    " // TRIM(gs_str))
       WRITE (gs_str,'(ES8.1)') this%GetScale(3)
       CALL this%Info("            geometry scale:    " // TRIM(gs_str))

    END SELECT
  END SUBROUTINE InitGeometry

  PURE FUNCTION GetScale1(this) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN) :: this
    REAL :: gp
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
    REAL, INTENT(IN) :: gp
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
  END SUBROUTINE SetScale1

  PURE SUBROUTINE SetScale2(this,gp,gp2)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: gp,gp2
    !------------------------------------------------------------------------!
    this%geoparam(1) = gp
    this%geoparam(2) = gp2
  END SUBROUTINE SetScale2

  PURE SUBROUTINE SetScale3(this,gp,gp2,gp3)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    REAL, INTENT(IN) :: gp,gp2,gp3
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
  PURE SUBROUTINE ScaleFactors_all(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)   :: this
    TYPE(MArrayV_TYP), INTENT(IN)      :: coords    !\TODO{ARRAY KLASSEN HIER NOCH EINTRAGEN (BISHER IM MESH)}
    TYPE(MArrayS_TYP), INTENT(OUT)     :: hx,hy,hz
    !------------------------------------------------------------------------!
    CALL this%ScaleFactors_1(coords%center ,hx%center ,hy%center ,hz%center)
    CALL this%ScaleFactors_1(coords%bcenter,hx%bcenter,hy%bcenter,hz%bcenter)
    CALL this%ScaleFactors_2(coords%faces  ,hx%faces  ,hy%faces  ,hz%faces)
    CALL this%ScaleFactors_2(coords%corners,hx%corners,hy%corners,hz%corners)
  END SUBROUTINE ScaleFactors_all

  PURE SUBROUTINE ScaleFactors_1(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)      :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
!CDIR IEXPAND    
       CALL this%ScaleFactors(coords(:,:,:,1),coords(:,:,:,2),coords(:,:,:,3), &
                                  hx(:,:,:)  ,    hy(:,:,:)  ,    hz(:,:,:))
  END SUBROUTINE ScaleFactors_1

  !> \public Compute scale factors
  !!
  !! Computes the scale factors of the given geometry.
  !! \f[
  !!   h_x = \left|\frac{\partial \vec{r}}{\partial x}\right|,\quad
  !!   h_y = \left|\frac{\partial \vec{r}}{\partial y}\right|,\quad
  !!   h_z = \left|\frac{\partial \vec{r}}{\partial z}\right|
  !! \f]
  PURE SUBROUTINE ScaleFactors_2(this,coords,hx,hy,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)        :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: hx,hy,hz
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%ScaleFactors(coords(:,:,:,:,1),coords(:,:,:,:,2),coords(:,:,:,:,3), &
                               hx(:,:,:,:)  ,    hy(:,:,:,:)  ,    hz(:,:,:,:))
  END SUBROUTINE ScaleFactors_2


  !> \public Compute radial distances to the origin
  !! for all cell positions
  PURE SUBROUTINE Radius_all(this,coords,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)     :: this
    TYPE(MArrayV_TYP), INTENT(IN)        :: coords
    TYPE(MArrayS_TYP), INTENT(OUT)       :: radius
    !------------------------------------------------------------------------!
    CALL this%Radius_1(coords%center ,radius%center)
    CALL this%Radius_1(coords%bcenter,radius%bcenter)
    CALL this%Radius_2(coords%faces  ,radius%faces)
    CALL this%Radius_2(coords%corners,radius%corners)
  END SUBROUTINE Radius_all

  !> \public Compute radial distances to the origin
  !!
  !! Computes the radial distances to the origin for the given
  !! coordinates depending on the geometry.
  PURE SUBROUTINE Radius_1(this,coords,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)      :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:)   :: radius
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Radius(coords(:,:,:,1),coords(:,:,:,2),coords(:,:,:,3),radius(:,:,:))
  END SUBROUTINE Radius_1

  !> \public Compute radial distances to the origin
  !!
  !! Computes the radial distances to the origin for the given
  !! coordinates depending on the geometry.
  PURE SUBROUTINE Radius_2(this,coords,radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)        :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:)   :: radius
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Radius(coords(:,:,:,:,1),coords(:,:,:,:,2),coords(:,:,:,:,3),radius(:,:,:,:))
  END SUBROUTINE Radius_2

  !> \public compute position vector components
  !! for all cell positions
  PURE SUBROUTINE PositionVector_all(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN) :: this
    TYPE(MArrayV_TYP), INTENT(IN)    :: coords
    TYPE(MArrayV_TYP), INTENT(OUT)   :: posvec
    !------------------------------------------------------------------------!
    CALL this%PositionVector_1(coords%center ,posvec%center)
    CALL this%PositionVector_1(coords%bcenter,posvec%bcenter)
    CALL this%PositionVector_2(coords%faces  ,posvec%faces)
    CALL this%PositionVector_2(coords%corners,posvec%corners)
  END SUBROUTINE PositionVector_all

  !> \public Compute position vector components
  !!
  !! Computes the curvilinear position vector components with respect to
  !! the given geometry.
  PURE SUBROUTINE PositionVector_1(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)      :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: posvec
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%PositionVector(coords(:,:,:,1),coords(:,:,:,2),coords(:,:,:,3), &
                             posvec(:,:,:,1),posvec(:,:,:,2),posvec(:,:,:,3)) 
  END SUBROUTINE PositionVector_1

  !> \public Compute position vector components
  !!
  !! Computes the curvilinear position vector components with respect to
  !! the given geometry.
  PURE SUBROUTINE PositionVector_2(this,coords,posvec)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)        :: this
    REAL, INTENT(IN),  DIMENSION(:,:,:,:,:) :: coords
    REAL, INTENT(OUT), DIMENSION(:,:,:,:,:) :: posvec
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%PositionVector(coords(:,:,:,:,1),coords(:,:,:,:,2),coords(:,:,:,:,3), &
                             posvec(:,:,:,:,1),posvec(:,:,:,:,2),posvec(:,:,:,:,3)) 
  END SUBROUTINE PositionVector_2

  !> \public Convert curvilinear to cartesian coordinates
  PURE SUBROUTINE Convert2Cartesian_all(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)     :: this  !> \param [in] this geometry class
    TYPE(MArrayV_TYP), INTENT(OUT)       :: curv  !> \param [in] curv curvilinear coordinates
    TYPE(MArrayV_TYP), INTENT(OUT)       :: cart  !> \param [out] cart cartesian coordinates
    !------------------------------------------------------------------------!
    CALL this%Convert2Cartesian_coords_1(curv%center ,cart%center)
    CALL this%Convert2Cartesian_coords_1(curv%bcenter,cart%bcenter)
    CALL this%Convert2Cartesian_coords_2(curv%faces  ,cart%faces)
    CALL this%Convert2Cartesian_coords_2(curv%corners,cart%corners)
  END SUBROUTINE Convert2Cartesian_all

  !> \public Convert curvilinear to cartesian coordinates
  PURE SUBROUTINE Convert2Cartesian_coords_1(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:)        :: curv, cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv
    INTENT(OUT)                     :: cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Cartesian(curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3), &
                                cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_1


  !> \public Convert curvilinear to cartesian coordinates
  PURE SUBROUTINE Convert2Cartesian_coords_2(this,curv,cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:)      :: curv, cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv
    INTENT(OUT)                     :: cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Cartesian(curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3), &
                                cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_coords_2


  !> \public Convert cartesian to curvilinear coordinates
  PURE SUBROUTINE Convert2Curvilinear_all(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)     :: this  !> \param [in] this geometry type
    TYPE(MArrayV_TYP), INTENT(IN)        :: cart  !> \param [in] cart cartesian coordinates
    TYPE(MArrayV_TYP), INTENT(OUT)       :: curv  !> \param [out] curv curvilinear coordinates
    !------------------------------------------------------------------------!
    CALL this%Convert2Curvilinear_coords_1(cart%center ,curv%center)
    CALL this%Convert2Curvilinear_coords_1(cart%bcenter,curv%bcenter)
    CALL this%Convert2Curvilinear_coords_2(cart%faces  ,curv%faces)
    CALL this%Convert2Curvilinear_coords_2(cart%corners,curv%corners)
  END SUBROUTINE Convert2Curvilinear_all


  !> \public Convert cartesian to curvilinear coordinates
  PURE SUBROUTINE Convert2Curvilinear_coords_1(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)     :: this  !> \param [in] this geometry type
    REAL, DIMENSION(:,:,:,:)             :: cart  !> \param [in] cart cartesian coordinates
    REAL, DIMENSION(:,:,:,:)             :: curv  !> \param [out] curv curvilinear coordinates
    !------------------------------------------------------------------------!
    INTENT(IN)                           :: cart
    INTENT(OUT)                          :: curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Curvilinear(cart(:,:,:,1),cart(:,:,:,2),cart(:,:,:,3), &
                                  curv(:,:,:,1),curv(:,:,:,2),curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_1


  !> \public Convert cartesian to curvilinear coordinates
  PURE SUBROUTINE Convert2Curvilinear_coords_2(this,cart,curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(IN)     :: this  !> \param [in] this geometry type
    REAL, DIMENSION(:,:,:,:,:)           :: cart, curv
    !------------------------------------------------------------------------!
    INTENT(IN)                           :: cart
    INTENT(OUT)                          :: curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Curvilinear(cart(:,:,:,:,1),cart(:,:,:,:,2),cart(:,:,:,:,3), &
                                  curv(:,:,:,:,1),curv(:,:,:,:,2),curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_coords_2

  !> \public Convert curvilinear vector components to cartesian vector components
  PURE SUBROUTINE Convert2Cartesian_vectors_1(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:)        :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv,v_curv
    INTENT(OUT)                     :: v_cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Cartesian(  curv(:,:,:,1),  curv(:,:,:,2),  curv(:,:,:,3), &
                                v_curv(:,:,:,1),v_curv(:,:,:,2),v_curv(:,:,:,3), &
                                v_cart(:,:,:,1),v_cart(:,:,:,2),v_cart(:,:,:,3))
  END SUBROUTINE Convert2Cartesian_vectors_1 

  !> \public Convert curvilinear vector components to cartesian vector components
  PURE SUBROUTINE Convert2Cartesian_vectors_2(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:)      :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv,v_curv
    INTENT(OUT)                     :: v_cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Cartesian(  curv(:,:,:,:,1),  curv(:,:,:,:,2),  curv(:,:,:,:,3), &
                                v_curv(:,:,:,:,1),v_curv(:,:,:,:,2),v_curv(:,:,:,:,3), &
                                v_cart(:,:,:,:,1),v_cart(:,:,:,:,2),v_cart(:,:,:,:,3))
  END SUBROUTINE Convert2Cartesian_vectors_2

  
  !> \public Convert curvilinear vector components to cartesian vector components
  PURE SUBROUTINE Convert2Cartesian_vectors_3(this,curv,v_curv,v_cart)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:)              :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv,v_curv
    INTENT(OUT)                     :: v_cart
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Cartesian(  curv(1),  curv(2),  curv(3), &
                                v_curv(1),v_curv(2),v_curv(3), &
                                v_cart(1),v_cart(2),v_cart(3))
  END SUBROUTINE Convert2Cartesian_vectors_3

  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors_1(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:)        :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv,v_cart
    INTENT(OUT)                     :: v_curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Curvilinear(  curv(:,:,:,1),  curv(:,:,:,2),  curv(:,:,:,3), &
                                  v_cart(:,:,:,1),v_cart(:,:,:,2),v_cart(:,:,:,3), &
                                  v_curv(:,:,:,1),v_curv(:,:,:,2),v_curv(:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_vectors_1

  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors_2(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:,:,:,:,:)      :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv,v_cart
    INTENT(OUT)                     :: v_curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Curvilinear(  curv(:,:,:,:,1),  curv(:,:,:,:,2),  curv(:,:,:,:,3), &
                                  v_cart(:,:,:,:,1),v_cart(:,:,:,:,2),v_cart(:,:,:,:,3), &
                                  v_curv(:,:,:,:,1),v_curv(:,:,:,:,2),v_curv(:,:,:,:,3))
  END SUBROUTINE Convert2Curvilinear_vectors_2

  !> \public Convert cartesian vector components to curvilinear vector components
  PURE SUBROUTINE Convert2Curvilinear_vectors_3(this,curv,v_cart,v_curv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base),INTENT(IN) :: this
    REAL, DIMENSION(:)              :: curv, v_curv, v_cart
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: curv,v_cart
    INTENT(OUT)                     :: v_curv
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    CALL this%Convert2Curvilinear(  curv(1),  curv(2),  curv(3),   &
                                  v_cart(1),v_cart(2),v_cart(3), &
                                  v_curv(1),v_curv(2),v_curv(3))
  END SUBROUTINE Convert2Curvilinear_vectors_3

  !> \public Destructor of generic geometry module
  SUBROUTINE FinalizeGeometry(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(geometry_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
         CALL this%Error("CloseGeometry","not initialized")
    !CALL this%logging_base%Finalize()
  END SUBROUTINE FinalizeGeometry

END MODULE geometry_base_mod
