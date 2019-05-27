!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2019                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
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
!> \defgroup mesh Mesh
!! \{
!! \brief Family of mesh modules
!! \}
!----------------------------------------------------------------------------!
!> \addtogroup mesh
!! - general parameters of mesh group as key-values
!! \key{type,INTEGER,spatial integration method}
!! \key{geometry,INTEGER,geometry of the mesh
!!      (see \ref geometry_base_mod for currently supported geometries)}
!! \key{shearingbox,INTEGER,enable/disable shearing box with shear along x-/y-coordinate
!!      \{1\,2\}\, 0 disables the shearing box,0}
!! \key{inum,INTEGER,x-resolution}
!! \key{jnum,INTEGER,y-resolution}
!! \key{xmin,REAL,minimum of x-coordinates}
!! \key{xmax,REAL,maximum of x-coordinates}
!! \key{ymin,REAL,minimum of y-coordinates}
!! \key{ymax,REAL,maximum of y-coordinates}
!! \key{gparam,REAL,1st geometry parameter}
!! \key{gparam2,REAL,2nd geometry parameter}
!! \key{gparam3,REAL,3rd geometry parameter}
!! \key{omega,REAL,angular speed of rotating frame of reference,0.0}
!! \key{Qshear,REAL,power law exponent of the rotation law \f$ Q=-\frac{d\ln\Omega}{d\ln r}\f$
!!      in the shearing box,1.5}
!! \key{rotcent,REAL,cartesian (x\,y)-coordiantes for center of rotation,(0\,0)}
!! - enable/disable output of certain mesh arrays
!! \key{output/bary,INTEGER,cell bary center in cartesian coordinates,1}
!! \key{output/bary_curv,INTEGER,cell bary center in curvilinear coordinates,1}
!! \key{output/corners,INTEGER,cell corners in cartesian coordinates,1}
!! \key{output/dl,INTEGER,line elements,0}
!! \key{output/bh,INTEGER,scale factors,0}
!! \key{output/volume,INTEGER,volume elements,0}
!! \key{output/dA,INTEGER,surface elements,0}
!! \key{output/radius,INTEGER,distance to the origin,0}
!! \key{output/position_vector,INTEGER,curvilinear position vector components,0}
!! \key{output/rotation,INTEGER,rotation angle of local orthonormal frames,0}

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic mesh module
!!
!! \ingroup mesh
!----------------------------------------------------------------------------!
MODULE mesh_base_mod
  USE logging_base_mod
  USE marray_base_mod
  USE marray_cellscalar_mod
  USE marray_cellvector_mod
  USE geometry_base_mod
  USE geometry_generic_mod
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
!--------------------------------------------------------------------------!
  PRIVATE
  !> \name Public Attributes
  !! #### mesh types
  INTEGER, PARAMETER :: MIDPOINT     = 1 !< use midpoint rule to approximate flux integrals
!  INTEGER, PARAMETER :: TRAPEZOIDAL  = 2 !< use trapezoidal rule to approximate flux integrals
  !> \}
  !! #### parameters depending on dimensionality
  INTEGER, PARAMETER :: NFACES(3)    = (/ 2, 4, 6 /)  !< number of faces
  INTEGER, PARAMETER :: NCORNERS(3)  = (/ 2, 4, 8 /)  !< number of corners
  !> \todo Check why this is not in boundary base?!
  INTEGER, PARAMETER :: &
     WEST   = 1, & !< named constant for western  boundary
     EAST   = 2, & !< named constant for eastern  boundary
     SOUTH  = 3, & !< named constant for southern boundary
     NORTH  = 4, & !< named constant for northern boundary
     BOTTOM = 5, & !< named constant for bottom   boundary
     TOP    = 6    !< named constant for top      boundary
  !> flags to check which vector components are enabled
  INTEGER, PARAMETER :: &
     VECTOR_X = INT(B'001'), &
     VECTOR_Y = INT(B'010'), &
     VECTOR_Z = INT(B'100')
  !> mesh data structure
  !PRIVATE
  TYPE,ABSTRACT, EXTENDS(logging_base) :: mesh_base
    !> \name Variables
    CLASS(geometry_base),ALLOCATABLE :: Geometry        !< geometrical properties
    TYPE(selection_base) :: without_ghost_zones !< for masking part of comp. domain
    INTEGER           :: GNUM              !< number of ghost cells
    INTEGER           :: GINUM,GJNUM,GKNUM !< number of ghost cells in any direction
    INTEGER           :: INUM,JNUM,KNUM    !< resolution
    INTEGER           :: IMIN,IMAX         !< minimal & maximal index in x-direction
    INTEGER           :: JMIN,JMAX         !< minimal & maximal index in y-direction
    INTEGER           :: KMIN,KMAX         !< minimal & maximal index in z-direction
    INTEGER           :: IGMIN,IGMAX       !< minimal & maximal index in x-direction with ghost cells
    INTEGER           :: JGMIN,JGMAX       !< minimal & maximal index in y-direction with ghost cells
    INTEGER           :: KGMIN,KGMAX       !< minimal & maximal index in z-direction with ghost cells
    INTEGER           :: NDIMS             !< dimensionality of the mesh: 1 (1D), 2 (2D), 3 (3D)
    INTEGER           :: NFACES            !< amount of faces, 2 (1D), 4 (2D), 6 (3D)
    INTEGER           :: NCORNERS          !< amount of corners, 2 (1D), 4 (2D), 8 (3D)
    INTEGER           :: Ip1,Ip2,Im1,Im2   !< access indices, which might become zero without x-dim
    INTEGER           :: Jp1,Jp2,Jm1,Jm2   !< access indices, which might become zero without y-dim
    INTEGER           :: Kp1,Kp2,Km1,Km2   !< access indices, which might become zero without z-dim
    REAL              :: xmin, xmax        !< spatial extent of computational domain in x-direction
    REAL              :: ymin, ymax        !< spatial extent of computational domain in y-direction
    REAL              :: zmin, zmax        !< spatial extent of computational domain in z-direction
    REAL              :: dx,dy,dz          !< curvilinear spatial differences
    REAL              :: invdx,invdy,invdz !< inverse of curvilinear spatial differences
    REAL              :: omega             !< speed of the rotating frame of ref.
    REAL              :: Q                 !< shearing parameter
    INTEGER           :: use_fargo         !< Fargo parameter (0 disabled, 1 enabled)
    INTEGER           :: shear_dir         !< Enables shearingbox with shear direction (0 disabled, 1 enabled)
    INTEGER           :: fargo             !< Fargo parameter (0 disabled, 1,2,3 enabled)
    INTEGER           :: ROTSYM            !< assume rotational symmetry (> 0 enabled, contains index of azimuthal angle)
    INTEGER           :: VECTOR_COMPONENTS !< enabled vector components
    REAL              :: rotcent(3)        !< center of the rotating frame of ref.
    !> \name
    !! #### cell coordinates
    TYPE(marray_cellvector) :: curv, &     !< curvilinear coordinates
                               cart        !< cartesian coordinates
    REAL, DIMENSION(:,:,:,:), POINTER :: &
                         center, &       !< geometrical centers
                         bcenter, &      !< bary centers
                         bccart          !< cartesian bary centers
    REAL, DIMENSION(:,:,:,:,:), POINTER :: &
                         fccart, &       !< cartesian face centered positions
                         ccart           !< cartesian corner positions
    !> \name
    !! #### line, area and volume elements
    TYPE(marray_base) :: dlx,dly,dlz, &        !< cell centered line elements
                         volume, &             !< cell volumes
                         dxdydV,dydzdV,dzdxdV  !< dx/volume and dy/volume
    REAL, DIMENSION(:,:,:,:), POINTER :: &
                         dAx,dAy,dAz, &            !< cell surface area elements
                         dAxdydz,dAydzdx,dAzdxdy  !< dAx/dydz, dAy/dxdz and dAz/dxdy
    !> \name
    !! #### scale factors and related quantities
    TYPE(marray_cellscalar) :: hx,hy,hz, &     !< scale factors
                               sqrtg, &        !< sqrt(det(metric))
                 cyxy,cyzy,cxzx,cxyx,czxz,czyz !< commutator coefficients
    !> \name
    !! #### radius and curvilinear position vector
    TYPE(marray_cellscalar) :: radius      !< real distance to coordinate origin
    TYPE(marray_cellvector) :: posvec      !< curvilinear position vector
    !> \name
    !! #### other geometrial quantities
    REAL, DIMENSION(:,:), POINTER :: &
                         rotation        !< rotation angle of local curvilinear orthonormal frames
    REAL, DIMENSION(:,:,:,:,:,:), POINTER :: &
                         weights         !< interpolation weights
#ifdef PARALLEL
    !> \name Variables in Parallel Mode
    INTEGER                    :: MAXINUM,MAXJNUM,MAXKNUM !< max. of local INUM,JNUM
    INTEGER                    :: MININUM,MINJNUM,MINKNUM !< min. of local INUM,JNUM
    INTEGER                    :: comm_cart               !< cartesian communicator
    INTEGER                    :: Icomm,Jcomm,Kcomm       !< communicators for cartesian rows and cols
    ! always setup a 3D cartesian process topology, even for 1D/2D simulations
    INTEGER, DIMENSION(NFACES(3)) :: comm_boundaries !< communicators for boundary processes
    INTEGER, DIMENSION(NFACES(3)) :: rank0_boundaries!< map rank0 -> world rank
    INTEGER, DIMENSION(NFACES(3)) :: neighbor        !< ranks of neighbor proc.
    INTEGER, DIMENSION(3)         :: dims            !< dimensions of cart comm
    INTEGER, DIMENSION(3)         :: mycoords        !< par. proc coordinates
#endif
  CONTAINS
    PROCEDURE :: RemapBounds_1
    PROCEDURE :: RemapBounds_2
    PROCEDURE :: RemapBounds_3
    PROCEDURE :: RemapBounds_4
    GENERIC   :: RemapBounds => RemapBounds_1, RemapBounds_2, RemapBounds_3, RemapBounds_4
    PROCEDURE :: InitMesh
    PROCEDURE :: Finalize_base
    PROCEDURE (Finalize), DEFERRED  :: Finalize
    PROCEDURE :: InternalPoint
    PROCEDURE (TensorDivergence3D),   DEFERRED :: TensorDivergence3D
    PROCEDURE (VectorDivergence3D),   DEFERRED :: VectorDivergence3D
    PROCEDURE (TensorDivergence2D_1),  DEFERRED :: TensorDivergence2D_1
    PROCEDURE (VectorDivergence2D_1),  DEFERRED :: VectorDivergence2D_1
    GENERIC  :: Divergence => TensorDivergence3D, VectorDivergence3D, &
                                     VectorDivergence2D_1, &! VectorDivergence2D_2, &
                                     TensorDivergence2D_1 !, TensorDivergence2D_2, &

  END TYPE mesh_base
  ABSTRACT INTERFACE
    PURE SUBROUTINE TensorDivergence3D(this,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz, &
                                       divTx,divTy,divTz)
      IMPORT mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(mesh_base), INTENT(IN) :: this
      REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
                    INTENT(IN) :: Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
      REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
                    INTENT(OUT) :: divTx,divTy,divTz
    END SUBROUTINE
    PURE SUBROUTINE VectorDivergence3D(this,vx,vy,vz,divv)
      IMPORT mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(mesh_base),INTENT(IN) :: this
      REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX) &
                        :: vx,vy,vz,divv
      !------------------------------------------------------------------------!
      INTENT(IN)        :: vx,vy,vz
      INTENT(OUT)       :: divv
    END SUBROUTINE
    PURE SUBROUTINE VectorDivergence2D_1(this,vx,vy,divv)
      IMPORT mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(mesh_base),INTENT(IN) :: this
      REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX) &
                        :: vx,vy,divv
      !------------------------------------------------------------------------!
      INTENT(IN)        :: vx,vy
      INTENT(OUT)       :: divv
    END SUBROUTINE
    PURE SUBROUTINE TensorDivergence2D_1(this,Txx,Txy,Tyx,Tyy,divTx,divTy)
      IMPORT mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(mesh_base),INTENT(IN)   :: this
      REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX) &
                        :: Txx,Txy,Tyx,Tyy,divTx,divTy
      !------------------------------------------------------------------------!
      INTENT(IN)        :: Txx,Txy,Tyx,Tyy
      INTENT(OUT)       :: divTx,divTy
    END SUBROUTINE
    SUBROUTINE Finalize(this)
      IMPORT mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(mesh_base), INTENT(INOUT) :: this
    END SUBROUTINE
  END INTERFACE
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       mesh_base, &
       ! constants
       PI, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       MIDPOINT, &
       CARTESIAN, &
       CYLINDRICAL, LOGCYLINDRICAL, &
       SPHERICAL, LOGSPHERICAL, SPHERICAL_PLANET, &
       VECTOR_X, VECTOR_Y, VECTOR_Z, &
       WEST, EAST, SOUTH, NORTH, BOTTOM, TOP
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public Constructor of generic mesh module
  !!
  !! This subroutine reads the necessary config data for setting up the mesh.
  !! It initializes the geometry and various mesh data arrays. Some of those
  !! are marked for output. For a detailed discription of the various geometries
  !! and how to setup those see \link geometry \endlink .
  SUBROUTINE InitMesh(this,config,IO,mtype,mname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base)        :: this   !< \param [in,out] this all mesh data
    TYPE(Dict_TYP),POINTER  :: config !< \param [in,out] config sub-dictionary
                                      !! with mesh configuration data
    TYPE(Dict_TYP),POINTER  :: IO     !< \param [in,out] IO mesh specific i/o data
    INTEGER                 :: mtype
    CHARACTER(LEN=32)       :: mname
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32)       :: xres,yres,zres,somega
    INTEGER                 :: use_fargo,fargo
    INTEGER                 :: i,j,k
    REAL                    :: mesh_dx,mesh_dy,mesh_dz
#ifdef PARALLEL
    INTEGER                 :: inum,jnum,knum,err
#endif
    REAL, DIMENSION(6,3)    :: cfaces
    REAL, DIMENSION(8,3)    :: ccorners
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    !------------------------------------------------------------------------!
    CALL this%logging_base%InitLogging(mtype,mname)

    CALL new_geometry(this%geometry, config)

    ! check whether shearing box simulation is enabled and along which
    ! coordinate direction the shear should be applied
    ! default: disable shearing box
    CALL GetAttr(config,"shearingbox",this%shear_dir,0)
    IF (this%shear_dir.LT.0.OR.this%shear_dir.GT.2) &
      CALL this%Error("mesh_base::InitMesh","direction of the shear should be one of {1,2} or 0")

    ! Here the angular velocity and the center of rotation are defined, if
    ! the mesh is in a rotating frame of reference. The fictious forces must
    ! be added through one of the following methods:
    ! 1. rotframe external source module
    ! 2. special physics module with angular momentum transport *IAMT, *IAMROT
    ! 3. special rhstype with angular momentum conservation
    ! Only the IAMROT physics module and rotframe source module allow
    ! for a different center of rotation than (/0, 0, 0/)

    ! angular velocity of the rotating reference frame
    CALL GetAttr(config, "omega", this%omega, 0.0)

    ! shearing force in shearingsheet (default: Keplerian rotation)
    CALL GetAttr(config, "Qshear", this%Q, 1.5)

    ! center of rotation
    this%rotcent = (/ 0., 0., 0./)
    CALL GetAttr(config, "rotation_center", this%rotcent, this%rotcent)



    ! total resolution
    ! IMPORTANT: The resolution is the key value in order to determine the
    !            used dimensions below
    CALL GetAttr(config, "inum", this%INUM)
    CALL GetAttr(config, "jnum", this%JNUM)
    CALL GetAttr(config, "knum", this%KNUM)

    ! do some sanity checks
    IF (ANY((/ this%INUM.LE.0, this%JNUM.LE.0, this%KNUM.LE.0 /))) &
      CALL this%Error("InitMesh","zero or negative resolution not allowed.")
    IF (ALL((/ this%INUM.EQ.1, this%JNUM.EQ.1, this%KNUM.EQ.1 /))) &
      CALL this%Error("InitMesh","at least one of inum,jnum,knum must be > 1.")

    ! coordinate domain
    CALL GetAttr(config, "xmin", this%xmin)
    CALL GetAttr(config, "xmax", this%xmax)
    CALL GetAttr(config, "ymin", this%ymin)
    CALL GetAttr(config, "ymax", this%ymax)
    CALL GetAttr(config, "zmin", this%zmin)
    CALL GetAttr(config, "zmax", this%zmax)

    !> \todo check if min/max coordinates are compatible with
    !!       the restrictions imposed by the geometry

    ! check if rotational symmetry is assumed, i.e., if the azimuthal
    ! direction consists of only one cell and a full 2*PI range
    ! this%ROTSYM is either set to the direction of the azimuthal angle {1,2,3}
    ! or 0 if rotational symmetry is disabled
    this%ROTSYM = 0     ! disable rotational symmetry by default
    SELECT CASE(this%geometry%GetAzimuthIndex())
    CASE(1)
      IF (this%INUM.EQ.1.AND.(ABS(this%xmax-this%xmin-2*PI).LE.EPSILON(this%xmin))) &
        this%ROTSYM = 1
    CASE(2)
      IF (this%JNUM.EQ.1.AND.(ABS(this%ymax-this%ymin-2*PI).LE.EPSILON(this%ymin))) &
        this%ROTSYM = 2
    CASE(3)
      IF (this%KNUM.EQ.1.AND.(ABS(this%zmax-this%zmin-2*PI).LE.EPSILON(this%zmin))) &
        this%ROTSYM = 3
    END SELECT

    ! These constants are used to access date from adjacent cells, i.e., with
    ! indices i,j,k +/- 1 and +/- 2. This should always be used instead of direct
    ! indexing, e.g. i-1 should become i+this%IM1. Some of the constants are
    ! set to 0 below, if there is only one cell in the particular direction.
    this%ip1 = 1
    this%ip2 = 2
    this%im1 = -1
    this%im2 = -2
    this%jp1 = 1
    this%jp2 = 2
    this%jm1 = -1
    this%jm2 = -2
    this%kp1 = 1
    this%kp2 = 2
    this%km1 = -1
    this%km2 = -2

    ! number of ghost cells; default is 2 at all boundaries
    ! it may be set to 0, if there is only one cell in
    ! a particular dimension (see below)
    this%GNUM  = 2
    this%GINUM = this%GNUM
    this%GJNUM = this%GNUM
    this%GKNUM = this%GNUM

    ! check dimensionality and suppress boundaries
    ! if there is only one cell in that direction
    this%NDIMS = 3 ! assume 3D simulation
    this%VECTOR_COMPONENTS = INT(B'111') ! enable all vector components
    IF (this%INUM.EQ.1) THEN
      IF (this%shear_dir.EQ.1) &
        CALL this%Error("mesh_base:Initmesh","shearing box enabled with direction 1, but INUM set to 1")
      this%NDIMS = this%NDIMS-1
      IF (this%ROTSYM.NE.1) this%VECTOR_COMPONENTS = IEOR(this%VECTOR_COMPONENTS,VECTOR_X)
      this%GINUM = 0
      this%ip1 = 0
      this%ip2 = 0
      this%im1 = 0
      this%im2 = 0
    END IF
    IF (this%JNUM.EQ.1) THEN
      IF (this%shear_dir.EQ.2) &
        CALL this%Error("mesh_base:Initmesh","shearing box enabled with direction 2, but JNUM set to 1")
      this%NDIMS = this%NDIMS-1
      IF (this%ROTSYM.NE.2) this%VECTOR_COMPONENTS = IEOR(this%VECTOR_COMPONENTS,VECTOR_Y)
      this%GJNUM = 0
      this%jp1 = 0
      this%jp2 = 0
      this%jm1 = 0
      this%jm2 = 0
    END IF
    IF (this%KNUM.EQ.1) THEN
      IF (this%shear_dir.EQ.3) &
        CALL this%Error("mesh_base:Initmesh","shearing box enabled with direction 3, but KNUM set to 1")
      this%NDIMS = this%NDIMS-1
      IF (this%ROTSYM.NE.3) this%VECTOR_COMPONENTS = IEOR(this%VECTOR_COMPONENTS,VECTOR_Z)
      this%GKNUM = 0
      this%kp1 = 0
      this%kp2 = 0
      this%km1 = 0
      this%km2 = 0
    END IF
    ! set number of faces and corners depending on dimensionality
    this%NFACES   = NFACES(this%NDIMS)
    this%NCORNERS = NCORNERS(this%NDIMS)

    ! set index ranges
    ! ATTENTION: all variables are local on MPI processes, hence this%IMIN on
    !            one process may differ from this%IMIN on another process
#ifdef PARALLEL
    ! compute the MPI partitioning, i.e. the decomposition of the computational domain,
    ! and setup various data structures MPI uses for communication
    CALL InitMesh_parallel(this, config)
    CALL MPI_Barrier(MPI_COMM_WORLD,err)
    ! determine the maximal amount of cells in 1st dimension
    INUM = this%IMAX - this%IMIN + 1
    CALL MPI_Allreduce(inum,this%MININUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(inum,this%MAXINUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    ! determine the maximal amount of cells in 2nd dimension
    JNUM = this%JMAX - this%JMIN + 1
    CALL MPI_Allreduce(jnum,this%MINJNUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(jnum,this%MAXJNUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    ! determine the maximal amount of cells in 3rd dimension
    KNUM = this%KMAX - this%KMIN + 1
    CALL MPI_Allreduce(knum,this%MINKNUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(knum,this%MAXKNUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
#else
    this%IMIN   = 1
    this%IMAX   = this%INUM
    this%JMIN   = 1
    this%JMAX   = this%JNUM
    this%KMIN   = 1
    this%KMAX   = this%KNUM
#endif

    ! index ranges with ghost cells
    this%IGMIN = this%IMIN - this%GINUM
    this%IGMAX = this%IMAX + this%GINUM
    this%JGMIN = this%JMIN - this%GJNUM
    this%JGMAX = this%JMAX + this%GJNUM
    this%KGMIN = this%KMIN - this%GKNUM
    this%KGMAX = this%KMAX + this%GKNUM

    ! initialize mesh arrays using indices with ghost cells
    CALL InitMeshProperties(this%IGMIN,this%IGMAX,this%JGMIN,this%JGMAX,this%KGMIN,this%KGMAX)

    ! create selection for the internal region
    this%without_ghost_zones = selection_base((/this%IMIN,this%IMAX,this%JMIN,this%JMAX,this%KMIN,this%KMAX/))

    ! coordinate differences in each direction
    this%dx = (this%xmax - this%xmin) / this%INUM
    mesh_dx = this%dx
    this%dy = (this%ymax - this%ymin) / this%JNUM
    mesh_dy = this%dy
    this%dz = (this%zmax - this%zmin) / this%KNUM
    mesh_dz = this%dz
    ! reset dx,dy,dz if direction is suppressed
    ! additionally set a mesh_dx,mesh_dy,mesh_dz to set corners, etc. 
    ! for output correctly (see below)
    IF (this%dx.LT.2*EPSILON(this%dx)) THEN
      IF (this%INUM.EQ.1) THEN
        this%dx = 1.0
        mesh_dx = 0.0
      ELSE
        CALL this%Error("mesh_base::InitMesh","INUM > 1 conflicts with zero x-extent")
      END IF
    END IF
    IF (this%dy.LT.2*EPSILON(this%dy)) THEN
      IF (this%JNUM.EQ.1) THEN
        this%dy = 1.0
        mesh_dy = 0.0
      ELSE
        CALL this%Error("mesh_base::InitMesh","JNUM > 1 conflicts with zero y-extent")
      END IF
    END IF
    IF (this%dz.LT.2*EPSILON(this%dz)) THEN
      IF (this%KNUM.EQ.1) THEN
        this%dz = 1.0
        mesh_dz = 0.0
      ELSE
        CALL this%Error("mesh_base::InitMesh","KNUM > 1 conflicts with zero z-extent")
      END IF
    END IF

    ! inverse coordinate differences
    this%invdx = 1./this%dx
    this%invdy = 1./this%dy
    this%invdz = 1./this%dz

    ! allocate memory for curvilinear positions
    this%curv = marray_cellvector()
    this%center  => this%curv%RemapBounds(this%curv%center)
    this%bcenter => this%curv%RemapBounds(this%curv%bcenter)

    ! create mesh arrays for scale factors
    this%hx = marray_cellscalar()
    this%hy = marray_cellscalar()
    this%hz = marray_cellscalar()

    ! create mesh array for square root of determinant of metric
    this%sqrtg = marray_cellscalar()

    ! create mesh arrays for commutator coefficients
    this%cxyx = marray_cellscalar()
    this%cxzx = marray_cellscalar()
    this%cyxy = marray_cellscalar()
    this%cyzy = marray_cellscalar()
    this%czxz = marray_cellscalar()
    this%czyz = marray_cellscalar()

    ! create mesh arrays for line and volume elements and related arrays
    this%volume = marray_base()
    this%dxdydV = marray_base()
    this%dydzdV = marray_base()
    this%dzdxdV = marray_base()
    this%dlx = marray_base()
    this%dly = marray_base()
    this%dlz = marray_base()

    ! nullify remaining mesh arrays
    NULLIFY(this%dAx,this%dAy,this%dAz, &
            this%dAxdydz,this%dAydzdx,this%dAzdxdy, &
            this%rotation,this%weights)

    ! translation vectors for cell faces and cell corners
    ! (with respect to geometrical cell center)
    cfaces(1,1) = -0.5*mesh_dx   ! western  x coordinate
    cfaces(1,2) =  0.0           ! western  y coordinate
    cfaces(1,3) =  0.0           ! western  z coordinate
    cfaces(2,1) =  0.5*mesh_dx   ! eastern  x coordinate
    cfaces(2,2) =  0.0           ! eastern  y coordinate
    cfaces(2,3) =  0.0           ! eastern  z coordinate
    cfaces(3,1) =  0.0           ! southern x coordinate
    cfaces(3,2) = -0.5*mesh_dy   ! southern y coordinate
    cfaces(3,3) =  0.0           ! southern z coordinate
    cfaces(4,1) =  0.0           ! northern x coordinate
    cfaces(4,2) =  0.5*mesh_dy   ! northern y coordinate
    cfaces(4,3) =  0.0           ! northern z coordinate
    cfaces(5,1) =  0.0           ! bottom   x coordinate
    cfaces(5,2) =  0.0           ! bottom   y coordinate
    cfaces(5,3) = -0.5*mesh_dz   ! bottom   z coordinate
    cfaces(6,1) =  0.0           ! top      x coordinate
    cfaces(6,2) =  0.0           ! top      y coordinate
    cfaces(6,3) =  0.5*mesh_dz   ! top      z coordinate

    ccorners(1,1) = cfaces(1,1) ! bottom-south-west
    ccorners(1,2) = cfaces(3,2)
    ccorners(1,3) = cfaces(5,3)
    ccorners(2,1) = cfaces(2,1) ! bottom-south-east
    ccorners(2,2) = cfaces(3,2)
    ccorners(2,3) = cfaces(5,3)
    ccorners(3,1) = cfaces(1,1) ! bottom-north-west
    ccorners(3,2) = cfaces(4,2)
    ccorners(3,3) = cfaces(5,3)
    ccorners(4,1) = cfaces(2,1) ! bottom-north-east
    ccorners(4,2) = cfaces(4,2)
    ccorners(4,3) = cfaces(5,3)
    ccorners(5,1) = cfaces(1,1) ! top-south-west
    ccorners(5,2) = cfaces(3,2)
    ccorners(5,3) = cfaces(6,3)
    ccorners(6,1) = cfaces(2,1) ! top-south-east
    ccorners(6,2) = cfaces(3,2)
    ccorners(6,3) = cfaces(6,3)
    ccorners(7,1) = cfaces(1,1) ! top-north-west
    ccorners(7,2) = cfaces(4,2)
    ccorners(7,3) = cfaces(6,3)
    ccorners(8,1) = cfaces(2,1) ! top-north-east
    ccorners(8,2) = cfaces(4,2)
    ccorners(8,3) = cfaces(6,3)
    ! calculate coordinates (geometric centers)
    DO k=this%KGMIN,this%KGMAX
       DO j=this%JGMIN,this%JGMAX
          DO i=this%IGMIN,this%IGMAX
             ! geometrical cell centers
             this%curv%center(i,j,k,1)    = this%xmin + (2*i-1)*0.5*mesh_dx    ! x coord
             this%curv%center(i,j,k,2)    = this%ymin + (2*j-1)*0.5*mesh_dy    ! y coord
             this%curv%center(i,j,k,3)    = this%zmin + (2*k-1)*0.5*mesh_dz    ! z coord
             ! cell face centered positions
             this%curv%faces(i,j,k,1,:)   = this%curv%center(i,j,k,:) + cfaces(1,:)   ! western
             this%curv%faces(i,j,k,2,:)   = this%curv%center(i,j,k,:) + cfaces(2,:)   ! eastern
             this%curv%faces(i,j,k,3,:)   = this%curv%center(i,j,k,:) + cfaces(3,:)   ! southern
             this%curv%faces(i,j,k,4,:)   = this%curv%center(i,j,k,:) + cfaces(4,:)   ! northern
             this%curv%faces(i,j,k,5,:)   = this%curv%center(i,j,k,:) + cfaces(5,:)   ! bottom
             this%curv%faces(i,j,k,6,:)   = this%curv%center(i,j,k,:) + cfaces(6,:)   ! top
             ! cell corner positions
             this%curv%corners(i,j,k,1,:) = this%curv%center(i,j,k,:) + ccorners(1,:) ! bottom-south-west
             this%curv%corners(i,j,k,2,:) = this%curv%center(i,j,k,:) + ccorners(2,:) ! bottom-south-east
             this%curv%corners(i,j,k,3,:) = this%curv%center(i,j,k,:) + ccorners(3,:) ! bottom-north-west
             this%curv%corners(i,j,k,4,:) = this%curv%center(i,j,k,:) + ccorners(4,:) ! bottom-north-east
             this%curv%corners(i,j,k,5,:) = this%curv%center(i,j,k,:) + ccorners(5,:) ! top-south-west
             this%curv%corners(i,j,k,6,:) = this%curv%center(i,j,k,:) + ccorners(6,:) ! top-south-east
             this%curv%corners(i,j,k,7,:) = this%curv%center(i,j,k,:) + ccorners(7,:) ! top-north-west
             this%curv%corners(i,j,k,8,:) = this%curv%center(i,j,k,:) + ccorners(8,:) ! top-north-east
          END DO
       END DO
    END DO

    ! set bary center curvilinear coordinates to geometric centers
    ! will be overwritten later
    this%curv%bcenter = this%curv%center

    use_fargo = 0
    IF(this%shear_dir.GT.0) use_fargo = 1
    CALL GetAttr(config, "use_fargo", this%use_fargo, use_fargo)
    ! enable fargo timestepping for polar geometries
    ! 1 = calculated mean background velocity w
    ! 2 = fixed user set background velocity w
    ! 3 = shearingsheet fixed background velocity w
    fargo = 0
    IF (this%use_fargo.NE.0) THEN
      IF(this%shear_dir.GT.0) fargo = 3
      CALL GetAttr(config, "fargo", this%FARGO, fargo)
    ELSE
      this%fargo = fargo
    END IF

    ! basic mesh initialization

    ! get geometrical scale factors for all cell positions
    ! bary center values are overwritten below
    CALL this%Geometry%ScaleFactors(this%curv,this%hx,this%hy,this%hz)

    ! get square root of determinant of the metric
    ! bary center values are overwritten below
    this%sqrtg = this%hx*(this%hy*this%hz)

    ! create mesh array for cartesian coordinates
    this%cart = marray_cellvector()
    this%bccart => this%cart%RemapBounds(this%cart%bcenter)
    this%fccart => this%cart%RemapBounds(this%cart%faces)
    this%ccart  => this%cart%RemapBounds(this%cart%corners)

    ! compute cartesian coordinates for all cell positions (center,bcenter,faces,corners)
    CALL this%Geometry%Convert2Cartesian(this%curv,this%cart)

    ! create mesh array for distance to the origin of the mesh
    this%radius = marray_cellscalar()
    CALL this%geometry%Radius(this%curv,this%radius)

    ! create mesh array for position vector
    this%posvec = marray_cellvector()
    CALL this%geometry%PositionVector(this%curv,this%posvec)

!    ! This is a quick hack for VTK and XDMF output on bipolar mesh.
!    ! It maps infinite coordinates on the boundary in cart%corners to finite values.
!    ! The bipolar mesh has a degeneracy at sigma=0(or 2*PI) and tau=0. If we
!    ! we transform to cartesian coordinates, these points would be mapped onto
!    ! a circle at infinity.
!    !> \todo NOT VERIFIED here is most probably something wrong in the 3D version
!    !! in most recent fosite 2d version this part is declared obsolete
!    IF (this%Geometry%GetType().EQ.BIPOLAR) THEN
!       ! get the maximal cartesian y coordinate at multiply by some appropriate number > 1
!       cart_max = 1.5*MAXVAL(ABS(this%cart%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,2)))
!       DO j=this%JMIN,this%JMAX+1
!          IF (this%IMIN.EQ.1) THEN          ! this is important for MPI
!             i = 1  ! boundary at tau_min
!             k = 1
!             ! check if tau changes sign
!             IF ((this%curv%corners(i,j,k,1,2).GT.0.0).AND.(this%curv%corners(i,j-1,k,1,2).LE.0.0)) THEN
!                ! map the cartesian corner to a value (x,y) = (0,cart_max)
!                IF (ABS(this%cart%corners(i,j-1,k,1,1)).LT.ABS(this%cart%corners(i,j,k,1,1))) THEN
!                    this%cart%corners(i,j,k,1,1)   = 0.0
!                    this%cart%corners(i,j,k,1,2)   = cart_max
!                    this%cart%corners(i,j-1,k,3,:) = this%cart%corners(i,j,k,1,:)
!                ELSE
!                    this%cart%corners(i,j-1,k,1,1) = 0.0
!                    this%cart%corners(i,j-1,k,1,2) = cart_max
!                    this%cart%corners(i,j-2,k,3,:) = this%cart%corners(i,j-1,k,1,:)
!                END IF
!             END IF
!          END IF
!          IF (this%IMAX.EQ.this%INUM) THEN  ! this is important for MPI
!             i = this%IMAX+1  ! boundary at tau_max
!             ! check if tau changes sign
!             IF ((this%curv%corners(i,j,k,1,2).GT.0.0).AND.(this%curv%corners(i,j-1,k,1,2).LE.0.0)) THEN
!                ! map the cartesian corner to a value (x,y) = (0,cart_max)
!                IF (ABS(this%cart%corners(i,j-1,k,1,1)).GT.ABS(this%cart%corners(i,j,k,1,1))) THEN
!                   this%cart%corners(i,j-1,k,1,1)   = 0.0
!                   this%cart%corners(i,j-1,k,1,2)   = -cart_max
!                   this%cart%corners(i-1,j-1,k,2,:) = this%cart%corners(i,j-1,k,1,:)
!                ELSE
!                   this%cart%corners(i,j,k,1,1)   = 0.0
!                   this%cart%corners(i,j,k,1,2)   = -cart_max
!                   this%cart%corners(i-1,j,k,2,:) = this%cart%corners(i,j,k,1,:)
!                END IF
!                EXIT
!             END IF
!          END IF
!       END DO
!    END IF

    ! print some information
    CALL this%Info(" MESH-----> quadrature rule:   " // TRIM(this%GetName()))
    WRITE (xres, '(I0)') this%INUM    ! this is just for better looking output
    WRITE (yres, '(I0)') this%JNUM
    WRITE (zres, '(I0)') this%KNUM
    CALL this%Info("            resolution:        " // TRIM(xres) // " x " // TRIM(yres) // " y "// TRIM(zres) // " z ")
    WRITE (xres, '(ES9.2,A,ES9.2)') this%xmin, " ..", this%xmax
    WRITE (yres, '(ES9.2,A,ES9.2)') this%ymin, " ..", this%ymax
    WRITE (zres, '(ES9.2,A,ES9.2)') this%zmin, " ..", this%zmax
    CALL this%Info("            computat. domain:  x=" // TRIM(xres))
    CALL this%Info("                               y=" // TRIM(yres))
    CALL this%Info("                               z=" // TRIM(zres))
    IF(this%omega.NE.0.) THEN
      WRITE (somega, '(ES9.2)') this%omega
      CALL this%Info("            ang. velocity:    " // TRIM(somega))
    END IF
    SELECT CASE(this%ROTSYM)
    CASE(0)
      somega = "disabled"
    CASE(1)
      somega = "x-coordinate"
    CASE(2)
      somega = "y-coordinate"
    CASE(3)
      somega = "z-coordinate"
    END SELECT
    CALL this%Info("            rot. symmetry:     " // TRIM(somega))
#ifdef PARALLEL
    WRITE (xres, '(I0)') this%dims(1)
    WRITE (yres, '(I0)') this%dims(2)
    WRITE (zres, '(I0)') this%dims(3)
    CALL this%Info("            MPI partition:     " // TRIM(xres) // " x " // TRIM(yres) // " y "// TRIM(zres) // " z ")
#endif

    ! \todo implement a correct check of the
    ! curvilinear range (e.g. xi > 0 in polar)


    CALL SetOutput(this,config,IO)

  END SUBROUTINE InitMesh

  !> \private Setup mesh fields for i/o
  SUBROUTINE SetOutput(this,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(INOUT) :: this   !< \param [in,out] this all mesh data
    TYPE(Dict_TYP),POINTER         :: config,IO
    !------------------------------------------------------------------------!
    INTEGER                        :: writefields
    !------------------------------------------------------------------------!
    !Set OutputDict
    CALL GetAttr(config, "output/corners", writefields, 1)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%ccart)) &
        !CALL SetAttr(IO,"corners",this%ccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX,:,:))
!                      Dict("name" / "coordinates:corners"))

    CALL GetAttr(config, "output/grid", writefields, 1)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%ccart)) THEN
        CALL SetAttr(IO,"grid_x",this%ccart(this%IMIN:this%IMAX+this%ip1, &
                                            this%JMIN:this%JMAX+this%jp1, &
                                            this%KMIN:this%KMAX+this%kp1,1,1))
        CALL SetAttr(IO,"grid_y",this%ccart(this%IMIN:this%IMAX+this%ip1, &
                                            this%JMIN:this%JMAX+this%jp1, &
                                            this%KMIN:this%KMAX+this%kp1,1,2))
        CALL SetAttr(IO,"grid_z",this%ccart(this%IMIN:this%IMAX+this%ip1, &
                                            this%JMIN:this%JMAX+this%jp1, &
                                            this%KMIN:this%KMAX+this%kp1,1,3))
    END IF

    CALL GetAttr(config, "output/bary_curv", writefields, 1)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%bcenter)) &
        CALL SetAttr(IO,"bary_curv",this%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:))
!                      Dict("name" / "coordinates:bary_curv"))
    CALL GetAttr(config, "output/bary", writefields, 1)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%bccart)) &
        CALL SetAttr(IO,"bary_centers",this%bccart(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:))
!                      Dict("name" / "coordinates:bary"))

    CALL GetAttr(config, "output/rotation", writefields, 0)
    IF(writefields.EQ.1) THEN
        CALL CalculateRotation(this)
        CALL SetAttr(IO,"rotation",this%rotation(this%IMIN:this%IMAX,this%JMIN:this%JMAX))
    END IF

    CALL GetAttr(config, "output/dl", writefields, 0)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%dlx%data3d)) &
           CALL SetAttr(IO,"dlx",this%dlx%data3d(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%dly%data3d)) &
           CALL SetAttr(IO,"dly",this%dly%data3d(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%dlz%data3d)) &
           CALL SetAttr(IO,"dlz",this%dlz%data3d(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
    END IF

    CALL GetAttr(config, "output/bh", writefields, 0)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%hx%bcenter)) &
           CALL SetAttr(IO,"bhx",this%hx%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%hy%bcenter)) &
           CALL SetAttr(IO,"bhy",this%hy%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%hz%bcenter)) &
           CALL SetAttr(IO,"bhz",this%hz%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
    END IF

    CALL GetAttr(config, "output/commutator", writefields, 0)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%cxyx%bcenter)) &
           CALL SetAttr(IO,"cxyx",this%cxyx%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%cxzx%bcenter)) &
           CALL SetAttr(IO,"cxzx",this%cxzx%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%cyxy%bcenter)) &
           CALL SetAttr(IO,"cyxy",this%cyxy%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%cyzy%bcenter)) &
           CALL SetAttr(IO,"cyzy",this%cyzy%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%czxz%bcenter)) &
           CALL SetAttr(IO,"czxz",this%czxz%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%czyz%bcenter)) &
           CALL SetAttr(IO,"czyz",this%czyz%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
    END IF

    CALL GetAttr(config, "output/volume", writefields, 0)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%volume%data3d)) &
        CALL SetAttr(IO,"volume",this%volume%data3d(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))

    CALL GetAttr(config, "output/dA", writefields, 0)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%dAx)) &
           CALL SetAttr(IO,"dAx",this%dAx(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:))
        IF (ASSOCIATED(this%dAy)) &
           CALL SetAttr(IO,"dAy",this%dAy(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:))
        IF (ASSOCIATED(this%dAz)) &
           CALL SetAttr(IO,"dAz",this%dAz(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:))
    END IF

    CALL GetAttr(config, "output/radius", writefields, 0)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%radius%bcenter)) &
        CALL SetAttr(IO,"radius",this%radius%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))

    CALL GetAttr(config, "output/position_vector", writefields, 0)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%posvec%bcenter)) &
        CALL SetAttr(IO,"position_vector",this%posvec%bcenter(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:))

  END SUBROUTINE SetOutput


  !> \private initialize array for rotation angle
  !!
  !! the rotation angle contains the local angle between the direction of the
  !! curve along which the second curvilinear coordinate varies (coordinate
  !! curve along the 1st curvilinear coordinate) and the cartesian x-direction
  !! (coordinate curve along the 1st cartesian coordinate)
  SUBROUTINE CalculateRotation(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base) :: this   !< \param [in,out] this all mesh data
    !------------------------------------------------------------------------!
    INTEGER          :: status
    REAL,DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX,3) :: cart,curv
    !------------------------------------------------------------------------!
        ALLOCATE(this%rotation(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
                 STAT=status)
        IF(status.NE.0) &
            CALL this%Error("CalculateRotation", "Couldn't allocate memory")
        ! initialize with vectors pointing in x-direction
        cart(:,:,:,1) = 1.
        cart(:,:,:,2) = 0.
        cart(:,:,:,3) = 0.
        CALL this%geometry%Convert2Curvilinear(this%bcenter, cart, curv)
        this%rotation(:,:) = ATAN2(curv(:,:,1,2),curv(:,:,1,1))  !!/todo{3D}
  END SUBROUTINE CalculateRotation

  !> \public Check if a given coordinate pair represents an internal point
  !!
  !! This function checks if a point given by its curvilinear coordinates (x,y)
  !! lies inside the computational domain or (if mask is given) inside a
  !! rectangular (with respect to the curvilinear mesh) region specified by an
  !! index mask.
  PURE FUNCTION InternalPoint(this,x,y,z,mask) RESULT(ip)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN)     :: this   !< \param [in,out] this all mesh data
    INTEGER, DIMENSION(6), OPTIONAL :: mask
    REAL                            :: x,y,z
    LOGICAL                         :: ip
    !------------------------------------------------------------------------!
    INTEGER                         :: imin,imax,jmin,jmax,kmin,kmax
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: x,y,z,mask
    !------------------------------------------------------------------------!
    ip = .FALSE.
    ! compare with the curvilinear extend (curvilinear domain is always rectangular)
    !
    ! the function behaves in 2 different ways: global and local
    !
    !   1. if mask is given, we do a local comparison, i.e.
    !      we check if the point lies inside the specified rectangular
    !      domain given by its cell indices: imin=mask(1),imax=mask(2),...
    !   2. otherwise we do a global check, i.e. with the whole computational
    !      domain
    IF (PRESENT(mask)) THEN
       ! restrict indices to the actual domain,
       ! which is different for each MPI process in parallel mode
       imin = MAX(this%IGMIN,mask(1))
       imax = MIN(this%IGMAX,mask(2))
       jmin = MAX(this%JGMIN,mask(3))
       jmax = MIN(this%JGMAX,mask(4))
       kmin = MAX(this%KGMIN,mask(5))
       kmax = MIN(this%KGMAX,mask(6))
       ! first check if the masked region is at least a subdomain of the actual domain
       IF (((imin.LE.this%IGMAX).AND.(imax.GE.imin)).AND. &
           ((jmin.LE.this%JGMAX).AND.(jmax.GE.jmin)).AND. &
           ((kmin.LE.this%KGMAX).AND.(kmax.GE.kmin))) THEN
          ! compare the curvilinear coordinates at the boundaries of the masked region
          ! with the transformed coordinates of the given point
          IF ((x.GE.this%curv%faces(imin,jmin,kmin,1,1).AND.x.LE.this%curv%faces(imax,jmax,kmax,2,1)).AND. &
              (y.GE.this%curv%faces(imin,jmin,kmin,1,2).AND.y.LE.this%curv%faces(imax,jmax,kmax,2,2)).AND. &
              (z.GE.this%curv%faces(imin,jmin,kmin,1,3).AND.z.LE.this%curv%faces(imax,jmax,kmax,2,3))) THEN
             ip = .TRUE.
          END IF
       END IF
    ELSE
       ! do the global check
       IF (((x.GE.this%xmin).AND.(x.LE.this%xmax)).AND. &
           ((y.GE.this%ymin).AND.(y.LE.this%ymax)).AND. &
           ((z.GE.this%zmin).AND.(z.LE.this%zmax))) THEN
          ip = .TRUE.
       END IF
    END IF
  END FUNCTION InternalPoint


#ifdef PARALLEL
  !> \public Initialize MPI (parallel mode only)
  !!
  !! Create a 3D cartesian topology of processes and subdivide
  !! the computational domain.
  !! ATTENTION: The dimensionality of the process topology is always 3D
  !!            even for 1D/2D simulations!
  SUBROUTINE InitMesh_parallel(this, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base), INTENT(INOUT)    :: this   !< \param [in,out] this all mesh data
    TYPE(Dict_TYP),POINTER             :: config
    !------------------------------------------------------------------------!
    CHARACTER(LEN=64)                  :: decomp_str
    LOGICAL, DIMENSION(3)              :: remain_dims, periods
    INTEGER                            :: num,rem
    INTEGER                            :: ierror
    INTEGER                            :: i,j,k
    INTEGER                            :: worldgroup,newgroup
    INTEGER, DIMENSION(1)              :: rank0in, rank0out
    INTEGER, DIMENSION(3)              :: coords,ncells,dims
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ranks
    !------------------------------------------------------------------------!
    ! 1. Check if the user has requests a particular decomposition.
    ! There are two ways a user may control automatic decomposition:
    ! (a) explicitly set the number of processes along a certain direction,
    !     e.g., decomposition = (/ 5,4,1 /) generates a cartesian communicator
    !     with 5 processes along the x-direction and 4 processes along the y-direction
    ! (b) if a number specified for certain dimension is negative let the algorithm
    !     find an optimal number of processes along that direction, e.g.,
    !     decomposition = (/ -1, 4, 1 /) generates a cartesian communicator
    !     with NumProcs / 4 processes along the x-direction and 4 processes
    !     along the y-direction with NumProcs beeing the total number of MPI
    !     processes given at the command line. If 4 is not a prime factor of
    !     NumProcs the program aborts with an error message.
    ! The default is automatic domain decomposition, i.e. finding the optimal
    ! distribution of processes along all directions considering the number of
    ! cells and the hardware vector length.
    ! This is done in CalculateDecomposition() (see below).

    ! defaults to -1, i.e., fully automatic decomposition
    dims(:)= -1
    CALL GetAttr(config, "decomposition", dims, dims(1:3))

    ! perform the decomposition on rank 0 and broadcast the result 
    ! to other processes (see below)
    IF (this%GetRank().EQ.0) THEN
      ! for more elaborate error output
      WRITE (decomp_str,'(3(A,I0),A)') "[ ", dims(1), ", ", dims(2), ", ",dims(3), " ]"

      ! set number of cells along each direction again
      ncells(1) = this%INUM
      ncells(2) = this%JNUM
      ncells(3) = this%KNUM

      ! perform some sanity checks
      IF (ALL(dims(:).GE.1).AND.PRODUCT(dims(:)).NE.this%GetNumProcs()) &
        CALL this%Error("InitMesh_parallel","total number of processes in domain decomposition " &
          // TRIM(decomp_str) // ACHAR(10) // REPEAT(' ',7) // &
          "does not match the number passed to mpirun")

      IF (ANY(dims(:).EQ.0)) &
        CALL this%Error("InitMesh_parallel","numbers in decomposition should not be 0")

      IF (MOD(this%GetNumProcs(),PRODUCT(dims(:),dims(:).GT.1)).NE.0) &
        CALL this%Error("InitMesh_parallel","numbers in decomposition " &
          // TRIM(decomp_str) // ACHAR(10) // REPEAT(' ',7) // &
          "are not devisors of the total number of processes passed to mpirun")

      ! balance number of processes if requested
      IF (ALL(dims(:).GT.0)) THEN
        ! (a) all dims user supplied -> do nothing
      ELSE IF (PRODUCT(dims(:)).GT.0) THEN
        ! (b) two dims < 0 one dim > 0 -> find optimal decomposition fixing one of the dims
        !     using the user supplied number
        k = MAXLOC(dims(:),1) ! get the one index k with dims(k) > 0
        ! suppress decomposition along the k-direction
        ncells(k) = 1
        i = dims(k) ! remember dims(k)
        IF (k.EQ.1) THEN
          ! switch entries 1 and 3 to make sure the first entry is not the one we are not decomposing
          dims(3) = this%GetNumProcs() / dims(k) ! reduced total number of processes
          dims(2) = 1
          dims(1) = 1
          CALL CalculateDecomposition(ncells(3),ncells(2),ncells(1),this%GKNUM, &
                                      dims(3),dims(2),dims(1))
        ELSE
          dims(1) = this%GetNumProcs() / dims(k) ! reduced total number of processes
          dims(2) = 1
          dims(3) = 1
          CALL CalculateDecomposition(ncells(1),ncells(2),ncells(3),this%GKNUM, &
                                      dims(1),dims(2),dims(3))
        END IF
        dims(k) = i
      ELSE IF (ALL(dims(:).LT.0)) THEN
        ! (c) all dims < 0 -> find optimal decomposition using all dimensions
        dims(1) = this%GetNumProcs()
        dims(2) = 1
        dims(3) = 1
        CALL CalculateDecomposition(ncells(1),ncells(2),ncells(3),this%GINUM, &
                                    dims(1),dims(2),dims(3))
      ELSE
        ! (d) one dim < 0 two dims > 0 -> fix the two dimensions and set the 3rd
        !     using the given number of processes
        k = MINLOC(dims(:),1) ! get the one index k with dims(k) < 0
        dims(k) = this%GetNumProcs() / PRODUCT(dims(:),dims(:).GT.1)
      END IF

      WRITE (decomp_str,'(3(A,I0),A)') "[ ", dims(1), ", ", dims(2), ", ",dims(3), " ]"

      ! set number of cells along each direction again
      ncells(1) = this%INUM
      ncells(2) = this%JNUM
      ncells(3) = this%KNUM
      DO i=1,3
        IF (dims(i).GT.ncells(i)) &
          CALL this%Error("InitMesh_parallel","number of processes exceeds number of cells " &
            // ACHAR(10) // REPEAT(' ',7) // "in dimension " // ACHAR(48+i) // ", check decomposition " &
            // TRIM(decomp_str))
      END DO

      IF (dims(3).LE.0) THEN
        CALL this%Error("InitMesh_parallel","automatic domain decomposition failed.")
      END IF

      this%dims(:) = dims(:)

    END IF
    
    CALL MPI_Bcast(this%dims,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    
! PRINT *,this%dims(:)
! CALL this%Error("InitMesh_parallel","debug breakpoint")

    ! 2. create the cartesian communicator
    ! IMPORTANT: disable reordering of nodes
    periods = .FALSE.
    CALL MPI_Cart_create(MPI_COMM_WORLD,3,this%dims,periods,.TRUE., &
         this%comm_cart,ierror)

    ! 3. inquire and save the own position
    this%mycoords(:) = 0
    CALL MPI_Cart_get(this%comm_cart,3,this%dims,periods,this%mycoords,ierror)

    ! 4. create communicators for every column and row of the cartesian topology
    remain_dims = (/ .TRUE., .FALSE., .FALSE. /)
    CALL MPI_Cart_Sub(this%comm_cart,remain_dims,this%Icomm,ierror)
    remain_dims = (/ .FALSE., .TRUE., .FALSE. /)
    CALL MPI_Cart_Sub(this%comm_cart,remain_dims,this%Jcomm,ierror)
    remain_dims = (/ .FALSE., .FALSE., .TRUE. /)
    CALL MPI_Cart_Sub(this%comm_cart,remain_dims,this%Kcomm,ierror)

    ! subdivide the mesh and set mesh indices
    ! x-direction
    rem = MOD(this%INUM,this%dims(1))            ! remainder
    num = (this%INUM-rem) / this%dims(1)         ! fraction
    IF (this%mycoords(1).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%IMIN = this%mycoords(1) * (num + 1) + 1
       this%IMAX = this%IMIN + num
    ELSE
       this%IMIN = rem * (num+1) + (this%mycoords(1) - rem) * num + 1
       this%IMAX = this%IMIN + num - 1
    END IF
    ! y-direction
    rem = MOD(this%JNUM,this%dims(2))            ! remainder
    num = (this%JNUM-rem) / this%dims(2)         ! fraction
    IF (this%mycoords(2).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%JMIN = this%mycoords(2) * (num + 1) + 1
       this%JMAX = this%JMIN + num
    ELSE
       this%JMIN = rem * (num+1) + (this%mycoords(2) - rem) * num + 1
       this%JMAX = this%JMIN + num - 1
    END IF
    ! z-direction
    rem = MOD(this%KNUM,this%dims(3))            ! remainder
    num = (this%KNUM-rem) / this%dims(3)         ! fraction
    IF (this%mycoords(3).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%KMIN = this%mycoords(3) * (num + 1) + 1
       this%KMAX = this%KMIN + num
    ELSE
       this%KMIN = rem * (num+1) + (this%mycoords(3) - rem) * num + 1
       this%KMAX = this%KMIN + num - 1
    END IF

    ! create communicators for all boundaries
    ALLOCATE(ranks(MAX(this%dims(1)*this%dims(2),this%dims(1)*this%dims(3),this%dims(2)*this%dims(3))))
    ranks = 0
    CALL MPI_Comm_group(MPI_COMM_WORLD,worldgroup,ierror)
    ! western boundary
    coords(1) = 0
    DO k=0,this%dims(3)-1
      DO j=0,this%dims(2)-1
        coords(2) = j
        coords(3) = k
        CALL MPI_Cart_rank(this%comm_cart,coords,i,ierror)
        ranks(j+1+k*this%dims(2))=i
      END DO
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(2)*this%dims(3),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(1),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(1) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! eastern boundary
    coords(1) = this%dims(1)-1
    DO k=0,this%dims(3)-1
      DO j=0,this%dims(2)-1
        coords(2) = j
        coords(3) = k
        CALL MPI_Cart_rank(this%comm_cart,coords,i,ierror)
        ranks(j+1+k*this%dims(2))=i
      END DO
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(2)*this%dims(3),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(2),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(2) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! southern boundary
    coords(2) = 0
    DO k=0,this%dims(3)-1
      DO i=0,this%dims(1)-1
        coords(1) = i
        coords(3) = k
        CALL MPI_Cart_rank(this%comm_cart,coords,j,ierror)
        ranks(i+1+k*this%dims(1))=j
      END DO
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1)*this%dims(3),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(3),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(3) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! northern boundary
    coords(2) = this%dims(2)-1
    DO k=0,this%dims(3)-1
      DO i=0,this%dims(1)-1
        coords(1) = i
        coords(3) = k
        CALL MPI_Cart_rank(this%comm_cart,coords,j,ierror)
        ranks(i+1+k*this%dims(1))=j
      END DO
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1)*this%dims(3),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(4),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(4) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! bottomer boundaries
    coords(3) = 0
    DO j=0,this%dims(2)-1
      DO i=0,this%dims(1)-1
        coords(1) = i
        coords(2) = j
        CALL MPI_Cart_rank(this%comm_cart,coords,k,ierror)
        ranks(i+1+j*this%dims(1))=k
      END DO
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1)*this%dims(2),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(5),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(5) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! topper boundaries
    coords(3) = this%dims(3)-1
    DO j=0,this%dims(2)-1
      DO i=0,this%dims(1)-1
        coords(1) = i
        coords(2) = j
        CALL MPI_Cart_rank(this%comm_cart,coords,k,ierror)
        ranks(i+1+j*this%dims(1))=k
      END DO
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1)*this%dims(2),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(6),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(6) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)

    CALL MPI_Group_free(worldgroup,ierror)
    DEALLOCATE(ranks)
  END SUBROUTINE InitMesh_parallel


 !> return the best partitioning of processes
 !!
 !! compute pi x pj x pk for a given mesh with resolution ni x nj x nk and
 !! number of ghost cells ginum in the 1st dimension (used for optimizing vector length)
 !! pi = pnum initially with total number of processes pnum<MAXNUM  (see module "factors")
 !! pj, pk = 1 initially
 !! pk returns 0, for erroneous input
 SUBROUTINE CalculateDecomposition(ni,nj,nk,ginum,pi,pj,pk)
   USE factors
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   INTEGER, INTENT(IN)    :: ni,nj,nk,ginum
   INTEGER, INTENT(INOUT) :: pi,pj,pk
   !------------------------------------------------------------------------!
   CHARACTER(LEN=8) :: svl_char
   INTEGER :: svl,err
   !------------------------------------------------------------------------!
   ! get the system vector length from preprocessor macro (string)
   svl_char = VECTOR_LENGTH
   READ (svl_char,'(I8)',IOSTAT=err) svl
   ! return immediatly for malformed input
   IF (ANY((/pi.LT.1,pi.GT.MAXNUM,pj.NE.1,pk.NE.1,err.NE.0, &
            MIN(ni,nj,nk).LT.1/))) THEN
      pk = 0
      RETURN
   END IF
   CALL Decompose(pi,pj,pk)

 CONTAINS

   !> searches for the best domain decomposition
   !!
   !! Accounts for the costs due to MPI communication (internal boundaries)
   !! and optimizes for a given vector length (VECTOR_LENGTH) on vector computers;
   !! parameters:
   !!   svl: system vector length
   !!   ni : number of grid cells in first dimension
   !!   nj : number of grid cells in second dimension
   !!   nk : number of grid cells in third dimension
   !!   pi : number of processes  in first dimension (=NumProcs first call)
   !!   pj : number of processes  in second dimension (=1 at first call)
   !!   pk : number of processes  in third dimension (=1 at first call)
   RECURSIVE SUBROUTINE Decompose(pi,pj,pk)
     IMPLICIT NONE
     !-------------------------------------------------------------------!
     INTEGER, INTENT(INOUT)  :: pi,pj,pk
     !-------------------------------------------------------------------!
     INTEGER :: pp,ptot
     INTEGER :: p1,p2,p3,pinew,pjnew,pknew,piold,pjold,pkold,pinew_,pjnew_,pknew_
     INTEGER :: pfmin,pfnew,pfold
     INTEGER :: bl,vl,blnew,vlnew,blnew_,vlnew_
     REAL    :: bl_gain,vl_gain
     !-------------------------------------------------------------------!
     ! return immediatly if the number of processes exceeds the number of
     ! grid cells, i.e. no subdivision possible
     IF (pj.GT.nj.OR.pk.GT.nk) RETURN
     ! measure the costs of the given configuration
     CALL GetCosts(ni,nj,nk,pi,pj,pk,bl,vl)
! PRINT '(3(A,I4),A,I7,A,I4)'," pi=",pi," pj=",pj," pk=",pk," boundary length=",bl," vector length=",vl
     ! save the configuration
     piold=pi
     pjold=pj
     pkold=pk
     p1 = pi
     p2 = pj
     p3 = pk
     ! compute the total number of processes
     ptot = pi*pj*pk
     pfmin = GetFactor(pj)      ! get smallest prime factor of pj
     pfold = 1
     pp = pi
     DO ! loop over all prime factors of pi which are larger than
        ! the smallest prime factor of pj
        IF (pp.LE.1) EXIT       ! if pj has been reduced to 1
        ! get smallest prime factor of pp, i.e. pj
        pfnew = GetFactor(pp)
        IF ((pfnew.GT.pfmin).AND.(pfmin.NE.1)) EXIT
        pp = pp/pfnew
        ! skip multiple prime factors
        IF (pfnew.NE.pfold) THEN
           ! create new configuration
           p1 = p1*pfold/pfnew
           p2 = p2/pfold*pfnew
           pfold = pfnew
           ! get the best configuration possible with p1 x p2 x p3 processes
           ! by reducing p1 => recursion
           pinew = p1
           pjnew = p2
           pknew = p3
           CALL Decompose(pinew,pjnew,pknew)
           CALL GetCosts(ni,nj,nk,pinew,pjnew,pknew,blnew,vlnew)
           ! do the same with p2 <-> p3 interchanged
           pinew_= p1
           pjnew_= p3
           pknew_= p2
           CALL Decompose(pinew_,pjnew_,pknew_)
           CALL GetCosts(ni,nj,nk,pinew_,pjnew_,pknew_,blnew_,vlnew_)
           ! check which of the 2 configurations is better
           bl_gain = blnew*(1.0/blnew_)             ! smaller is better
           vl_gain = vlnew_*(1.0/vlnew)             ! larger is better
!!$PRINT '(4I7)',pjold,pkold,bl,vl
!!$PRINT '(4I7)',pjnew,pknew,blnew,vlnew
!!$PRINT *,"--------------------------------------------"
!!$PRINT '(A,3F7.1)',"              ",bl_gain,vl_gain,bl_gain*vl_gain
           IF (vl_gain*bl_gain.GT.1.0) THEN
             ! 2nd configuration is better
             pinew = pinew_
             pjnew = pjnew_
             pknew = pknew_
             blnew = blnew_
             vlnew = vlnew_
           END IF
           ! check if the old configuration is better
           bl_gain = bl*(1.0/blnew)             ! smaller is better
           vl_gain = vlnew*(1.0/vl)             ! larger is better
           IF (vl_gain*bl_gain.GT.1.0) THEN
             ! new configuration is better
             piold = pinew
             pjold = pjnew
             pkold = pknew
             bl    = blnew
             vl    = vlnew
           END IF
        END IF
     END DO
     ! return optimized number of processes in the first dimension
     ! (for the initial configuration pj x pk)
     pj=pjold
     pk=pkold
     pi=ptot/pj/pk
   END SUBROUTINE Decompose

   ! computes the sum of the length of all internal boundaries (bl)
   ! and the mean vector length (vl)
   PURE SUBROUTINE GetCosts(n1,n2,n3,p1,p2,p3,bl,vl)
     IMPLICIT NONE
     !------------------------------------------------------------------------!
     INTEGER, INTENT(IN)  :: n1,n2,n3,p1,p2,p3
     INTEGER, INTENT(OUT) :: bl,vl
     !------------------------------------------------------------------------!
     INTEGER :: num
     !------------------------------------------------------------------------!
     ! length of internal boundaries
     bl = n3*n2*(p1-1) + n1*n3*(p2-1)  + n1*n2*(p3-1)
     ! estimate the average vector length if the system vector length is large than 1
     IF (svl.GT.1) THEN
        ! get the maximal number of cells along the first dimension
        ! including ghost cells on both ends
        num = n1 / p1  + 2*ginum
        IF (MOD(n1,p1).NE.0) num = num + 1 ! if the remainder is not zero add 1
        IF (num.LE.svl) THEN
          ! num fits into one vector
          vl = num
        ELSE
          ! we need more the one vector which may be filled completely or only partly
          ! return the closest integer to the average filling
          vl = num / ((num / svl) + 1)
        END IF
     ELSE
       vl = 1
     END IF
   END SUBROUTINE GetCosts

 END SUBROUTINE CalculateDecomposition
#endif
!---------------------END PARALLEL-------------------------------------------!


  !> \public remap lower bounds in the first 2 dimensions of rank 2 subarrays
  !!
  !! This is a short hack to obviate a restriction in the generation of
  !! subarray pointers. In Fortran 90/95 the indices of subarrays always
  !! start with a lower bound of 1, but Fosite requires that all mesh data
  !! arrays start with lower bounds of IGMIN and JGMIN, which are not equal to
  !! 1 in general.
  FUNCTION RemapBounds_1(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base)                :: this !< \param [in,out] this all mesh data
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:,this%KGMIN:), TARGET :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                      :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_1

  !> \public remap lower bounds in the first 2 dimensions of rank 3 subarrays
  FUNCTION RemapBounds_2(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base)                  :: this !< \param [in,out] this all mesh data
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:,this%KGMIN:,:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                        :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_2

  !> \public remap lower bounds in the first 2 dimensions of rank 4 subarrays
  FUNCTION RemapBounds_3(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base)                    :: this !< \param [in,out] this all mesh data
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:,this%KGMIN:,:,:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_3

  !> \public remap lower bounds in the first 2 dimensions of rank 5 subarrays
  FUNCTION RemapBounds_4(this,array) RESULT(ptr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base)                      :: this !< \param [in,out] this all mesh data
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:,this%KGMIN:,:,:,:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:,:,:), POINTER :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)                            :: array
    !------------------------------------------------------------------------!
    ptr => array
  END FUNCTION RemapBounds_4


  !> \public Destructor of mesh class
  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(INOUT) :: this   !< \param [in,out] this all mesh data
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER :: i,ierror
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("CloseMesh","not initialized")

    DO i=1,6
       IF (this%comm_boundaries(i).NE.MPI_COMM_NULL) &
          CALL MPI_Comm_free(this%comm_boundaries(i),ierror)
    END DO
#endif


    CALL this%curv%Destroy()
    CALL this%hx%Destroy()
    CALL this%hy%Destroy()
    CALL this%hz%Destroy()
    CALL this%sqrtg%Destroy()
    CALL this%cxyx%Destroy()
    CALL this%cxzx%Destroy()
    CALL this%cyxy%Destroy()
    CALL this%cyzy%Destroy()
    CALL this%czxz%Destroy()
    CALL this%czyz%Destroy()

    CALL this%volume%Destroy()
    CALL this%dxdydV%Destroy()
    CALL this%dydzdV%Destroy()
    CALL this%dzdxdV%Destroy()

    CALL this%dlx%Destroy()
    CALL this%dly%Destroy()
    CALL this%dlz%Destroy()

    CALL this%cart%Destroy()
    CALL this%radius%Destroy()
    CALL this%posvec%Destroy()

    IF (ASSOCIATED(this%rotation)) DEALLOCATE(this%rotation)

    CALL this%without_ghost_zones%Destroy()

    CALL CloseMeshProperties

    CALL this%Geometry%Finalize()
    DEALLOCATE(this%Geometry)
  END SUBROUTINE Finalize_base


END MODULE mesh_base_mod
