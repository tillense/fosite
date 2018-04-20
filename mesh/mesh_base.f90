!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
!> \addtogroup mesh
!! - general parameters of mesh group as key-values
!! \key{type,INTEGER,spatial integration method
!!      (see \link mesh_generic \endlink for currently supported mesh types)}
!! \key{geometry,INTEGER,geometry of the mesh
!!      (see \link geometry_generic \endlink for currently supported geometries)}
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
!! \key{rotcent,REAL,cartesian (x\,y)-coordiantes for center of rotation,(0\,0)
!! \key{dz,REAL,extent of 3rd dimension}
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
!! \brief base mesh module
!!
!! \ingroup mesh
!----------------------------------------------------------------------------!
MODULE mesh_base_mod
  USE logging_base_mod
  USE array
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
  INTEGER, PARAMETER :: TRAPEZOIDAL  = 2 !< use trapezoidal rule to approximate flux integrals
  !> \}
  !> \todo
  INTEGER, PARAMETER :: NDIMS  = 3       ! dimensions of cartesian topology
  INTEGER, PARAMETER :: NFACES = 6       ! amount of faces
  INTEGER, PARAMETER :: VECLEN = &       ! vector length ...
#if defined(NECSX8) || defined(NECSX9) || defined(NECSXACE)
  256                                    ! ... of NEC SX8/SX9/SX-ACE CPUs
#else
  1                                      ! ... of everthing else
#endif
  !> \todo Check why this is not in boundary base?!
  INTEGER, PARAMETER :: &
     WEST   = 1, & !< named constant for western  boundary
     EAST   = 2, & !< named constant for eastern  boundary
     SOUTH  = 3, & !< named constant for southern boundary
     NORTH  = 4, & !< named constant for northern boundary
     BOTTOM = 5, & !< named constant for bottom   boundary
     TOP    = 6    !< named constant for top      boundary
  !> index selection type
  TYPE Selection_TYP
     !> \name Variables
     INTEGER           :: IMIN,IMAX       !< selection in x-direction
     INTEGER           :: JMIN,JMAX       !< selection in y-direction
     INTEGER           :: KMIN,KMAX       !< selection in z-direction
     LOGICAL, POINTER  :: mask(:,:)       !< optional selection mask
  END TYPE Selection_TYP
  !> mesh data structure
  TYPE,ABSTRACT, EXTENDS(logging_base) :: mesh_base
    !PRIVATE
    !> \name Variables
    CLASS(geometry_base),ALLOCATABLE :: Geometry        !< geometrical properties
    INTEGER           :: GNUM              !< number of ghost cells
    INTEGER           :: GINUM,GJNUM,GKNUM !< number of ghost cells in any direction
    INTEGER           :: INUM,JNUM,KNUM    !< resolution
    INTEGER           :: IMIN,IMAX         !< minimal & maximal index in x-direction
    INTEGER           :: JMIN,JMAX         !< minimal & maximal index in y-direction
    INTEGER           :: KMIN,KMAX         !< minimal & maximal index in z-direction
    INTEGER           :: IGMIN,IGMAX       !< minimal & maximal index in x-direction with ghost cells
    INTEGER           :: JGMIN,JGMAX       !< minimal & maximal index in y-direction with ghost cells
    INTEGER           :: KGMIN,KGMAX       !< minimal & maximal index in z-direction with ghost cells
    INTEGER           :: NDIMS              !< amount of dimension, 1 (1D), 2 (2D), 3 (3D)
    INTEGER           :: NFACES            !< amount of faces, 2 (1D), 4 (2D), 6 (3D)
    INTEGER           :: NCORNERS          !< amount of faces, 2 (1D), 4 (2D), 8 (3D)
    INTEGER           :: Ip1,Ip2,Im1,Im2   !< access indices, which might become zero without x-dim
    INTEGER           :: Jp1,Jp2,Jm1,Jm2   !< access indices, which might become zero without y-dim
    INTEGER           :: Kp1,Kp2,Km1,Km2   !< access indices, which might become zero without z-dim
    REAL              :: xmin, xmax        !< spatial extent of computational domain in x-direction
    REAL              :: ymin, ymax        !< spatial extent of computational domain in y-direction
    REAL              :: zmin, zmax        !< spatial extent of computational domain in z-direction
    REAL              :: dx,dy,dz          !< curvilinear spatial differences
    REAL              :: invdx,invdy,invdz !< inverse of curvilinear spatial differences
    REAL              :: omega             !< speed of the rotating frame of ref.
    REAL              :: rotcent(2)        !< center of the rotating frame of ref.
    !> \name
    !! #### cell coordinates
    TYPE(MArrayV_TYP) :: curv, &         !< curvilinear coordinates
                         cart            !< cartesian coordinates
    REAL, DIMENSION(:,:,:,:), POINTER :: &
                         center, &       !< geometrical centers
                         bcenter, &      !< bary centers
                         bccart          !< cartesian bary centers
    REAL, DIMENSION(:,:,:,:,:), POINTER :: &
                         fccart, &       !< cartesian face centered positions
                         ccart           !< cartesian corner positions
    !> \name
    !! #### line, area and volume elements
    REAL, DIMENSION(:,:,:), POINTER :: &
                         dlx,dly,dlz, &        !< cell centered line elements
                         volume, &             !< cell volumes
                         dxdydV,dydzdV,dzdxdV  !< dx/volume and dy/volume
    REAL, DIMENSION(:,:,:,:), POINTER :: &
                         dAx,dAy,dAz, &            !< cell surface area elements
                         dAxdydz,dAydzdx,dAzdxdy  !< dAx/dydz, dAy/dxdz and dAz/dxdy
    !> \name
    !! #### scale factors and related quantities
    TYPE(MArrayS_TYP) :: hx,hy,hz, &     !< scale factors
                         sqrtg, &        !< sqrt(det(metric))
                         cyxy,cyzy,cxzx,cxyx,czxz,czyz !< commutator coefficients
    !> \name
    !! #### radius and curvilinear position vector
    TYPE(MArrayS_TYP) :: radius      !< real distance to coordinate origin
    TYPE(MArrayV_TYP) :: posvec      !< curvilinear position vector
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
    !todo: allocable, in initialisierung mit this%nfaces allokieren
    INTEGER, DIMENSION(NFACES) :: comm_boundaries !< communicators for boundary processes
    INTEGER, DIMENSION(NFACES) :: rank0_boundaries!< map rank0 -> world rank
    INTEGER, DIMENSION(NFACES) :: neighbor        !< ranks of neighbor proc.
    INTEGER, DIMENSION(NDIMS)  :: dims            !< dimensions of cart comm
    INTEGER, DIMENSION(NDIMS)  :: mycoords        !< par. proc coordinates
#endif
  CONTAINS
    PRIVATE
    PROCEDURE :: AllocateMesharrayS
    PROCEDURE :: AllocateMesharrayV
    PROCEDURE :: AllocateMesharrayT
    PROCEDURE :: DeallocateMesharrayS
    PROCEDURE :: DeallocateMesharrayV
    PROCEDURE :: DeallocateMesharrayT
    PROCEDURE :: RemapBounds_1
    PROCEDURE :: RemapBounds_2
    PROCEDURE :: RemapBounds_3
    PROCEDURE :: RemapBounds_4
    GENERIC, PUBLIC :: AllocateMesharray => AllocateMesharrayS, AllocateMesharrayV, AllocateMesharrayT
    GENERIC, PUBLIC :: DeallocateMesharray => DeallocateMesharrayS, DeallocateMesharrayV, DeallocateMesharrayT
    GENERIC, PUBLIC :: RemapBounds => RemapBounds_1, RemapBounds_2, RemapBounds_3, RemapBounds_4
    PROCEDURE, PUBLIC :: InitMesh
    PROCEDURE, PUBLIC :: Finalize
    PROCEDURE, PUBLIC :: InternalPoint
    PROCEDURE (TensorDivergence3D), PUBLIC,  DEFERRED :: TensorDivergence3D
    PROCEDURE (VectorDivergence3D), PUBLIC,  DEFERRED :: VectorDivergence3D
    PROCEDURE (TensorDivergence2D_1), PUBLIC,  DEFERRED :: TensorDivergence2D_1
    PROCEDURE (VectorDivergence2D_1), PUBLIC,  DEFERRED :: VectorDivergence2D_1
    GENERIC, PUBLIC :: DIVERGENCE => TensorDivergence3D, VectorDivergence3D, & 
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
  END INTERFACE
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       mesh_base, &
       Selection_TYP, &
       ! constants
       PI, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       MIDPOINT, TRAPEZOIDAL, &
       CARTESIAN, POLAR, LOGPOLAR, TANPOLAR, SINHPOLAR, SINHTANHPOLAR, BIPOLAR, &
       CYLINDRICAL, TANCYLINDRICAL, LNCOSHCYLINDRICAL, SPHERICAL, SINHSPHERICAL, &
       BIANGLESPHERICAL, OBLATE_SPHEROIDAL, CHANNEL, ELLIPTIC, SINHCARTESIAN, &
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
    INTEGER                 :: meshtype
    INTEGER                 :: i,j,k,err
    REAL                    :: mesh_dx,mesh_dy,mesh_dz
#ifdef PARALLEL
    INTEGER                 :: inum, jnum, knum
#endif
    REAL, DIMENSION(6,3)    :: cfaces
    REAL, DIMENSION(8,3)    :: ccorners
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    !------------------------------------------------------------------------!
    CALL this%logging_base%InitLogging(mtype,mname)

    CALL new_geometry(this%geometry, config)

    ! total resolution
    ! IMPORTANT: The resolution is the key value in order to determine the
    !            used dimensions below
    CALL GetAttr(config, "inum", this%inum)
    CALL GetAttr(config, "jnum", this%jnum)
    CALL GetAttr(config, "knum", this%knum)

    ! set access scalars depending of used dimension
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

    ! number of ghost rows/columns
    this%GNUM  = 2
    this%GINUM = this%GNUM
    this%GJNUM = this%GNUM
    this%GKNUM = this%GNUM

    ! reset values if dimensions do not exist
    IF (this%INUM.EQ.1) THEN
      this%GINUM = 0
      this%ip1 = 0
      this%ip2 = 0
      this%im1 = 0
      this%im2 = 0
    ELSE IF (this%JNUM.EQ.1) THEN
      this%GJNUM = 0
      this%jp1 = 0
      this%jp2 = 0
      this%jm1 = 0
      this%jm2 = 0
    ELSE IF (this%KNUM.EQ.1) THEN
      this%GKNUM = 0
      this%kp1 = 0
      this%kp2 = 0
      this%km1 = 0
      this%km2 = 0
    END IF

    ! set amount of faces
    IF  & ! case 1D
      ((((this%INUM.EQ.1).AND.(this%JNUM.EQ.1)).AND.(this%KNUM.GT.1)) .OR. &
       (((this%JNUM.EQ.1).AND.(this%KNUM.EQ.1)).AND.(this%INUM.GT.1)) .OR. &
       (((this%KNUM.EQ.1).AND.(this%INUM.EQ.1)).AND.(this%JNUM.GT.1))) THEN
      this%NFACES = 2
      this%NCORNERS = 2
      this%NDIMS = 1
    ELSE IF  & ! case 2D
      ((((this%INUM.EQ.1).AND.(this%JNUM.GT.1)).AND.(this%KNUM.GT.1)) .OR. &
       (((this%JNUM.EQ.1).AND.(this%KNUM.GT.1)).AND.(this%INUM.GT.1)) .OR. &
       (((this%KNUM.EQ.1).AND.(this%INUM.GT.1)).AND.(this%JNUM.GT.1))) THEN
      this%NFACES = 4
      this%NCORNERS = 4
      this%NDIMS = 2
    ELSE IF  & ! case 3D
      (((this%INUM.GT.1).AND.(this%JNUM.GT.1)).AND.(this%KNUM.GT.1)) THEN
      this%NFACES = 6
      this%NCORNERS = 6
      this%NDIMS = 3
    ELSE
      CALL this%Error("InitMesh","Cell numbering is not allowed.")
    END IF

    ! coordinate domain
    CALL GetAttr(config, "xmin", this%xmin)
    CALL GetAttr(config, "xmax", this%xmax)
    CALL GetAttr(config, "ymin", this%ymin)
    CALL GetAttr(config, "ymax", this%ymax)
    CALL GetAttr(config, "zmin", this%zmin)
    CALL GetAttr(config, "zmax", this%zmax)

    ! coordinate differences in each direction
    this%dx = (this%xmax - this%xmin) / this%inum
    mesh_dx = this%dx
    this%dy = (this%ymax - this%ymin) / this%jnum
    mesh_dy = this%dy
    this%dz = (this%zmax - this%zmin) / this%knum
    mesh_dz = this%dz
    ! reset dx,dy,dz if dimension does not exists,
    ! additionally set a mesh_dx,mesh_dy,mesh_dz to set corners, etc.
    ! for output correctly
    IF (this%inum.EQ.1 .AND. this%dx.EQ.0.0) THEN
      this%dx = 1.0
      mesh_dx = 0.0
    ELSE IF(this%jnum.EQ.1 .AND. this%dy.EQ.0.0) THEN
      this%dy = 1.0
      mesh_dy = 0.0
    ELSE IF(this%knum.EQ.1 .AND. this%dz.EQ.0.0) THEN
      this%dz = 1.0
      mesh_dz = 0.0
    END IF

    ! inverse coordinate differences
    this%invdx = 1./this%dx
    this%invdy = 1./this%dy
    this%invdz = 1./this%dz

    ! set index ranges
#ifdef PARALLEL
    CALL InitMesh_parallel(this, config)
    CALL MPI_Barrier(MPI_COMM_WORLD,err)
    inum = this%IMAX - this%IMIN + 1
    CALL MPI_Allreduce(inum,this%MININUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(inum,this%MAXINUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    jnum = this%JMAX - this%JMIN + 1
    CALL MPI_Allreduce(jnum,this%MINJNUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(jnum,this%MAXJNUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    knum = this%KMAX - this%KMIN + 1
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

    ! allocate memory for curvilinear positions
    CALL this%AllocateMesharray(this%curv)
    this%center  => this%RemapBounds(this%curv%center)
    this%bcenter => this%RemapBounds(this%curv%bcenter)

    ! allocate memory for scale factors
    CALL this%AllocateMesharray(this%hx)
    CALL this%AllocateMesharray(this%hy)
    CALL this%AllocateMesharray(this%hz)

    ! allocate memory for square root of determinant of metric
    CALL this%AllocateMesharray(This%sqrtg)

    ! allocate memory for commutator coefficients
    CALL this%AllocateMesharray(this%cxyx)
    CALL this%AllocateMesharray(this%cxzx)
    CALL this%AllocateMesharray(this%cyxy)
    CALL this%AllocateMesharray(this%cyzy)
    CALL this%AllocateMesharray(this%czxz)
    CALL this%AllocateMesharray(this%czyz)

    ! allocate memory for all pointers that are independent of fluxtype
    ALLOCATE( &
         this%volume(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
         this%dxdydV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
         this%dydzdV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
         this%dzdxdV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
         this%dlx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
         this%dly(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
         this%dlz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX), &
        STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("InitMesh_common","Unable to allocate memory!")
    END IF

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
    cfaces(5,3) =  -0.5*mesh_dz  ! bottom   z coordinate
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

    ! set bary center curvilinear coordinates to geometric centers
    ! will be overwritten later
    CALL GetAttr(config, "meshtype", meshtype)

    ! Here the angular velocity and the center of rotation are defined, if
    ! the mesh is in a rotating frame of reference. The fictious forces must
    ! be added through one of the following methods:
    ! 1. rotframe external source module
    ! 2. special physics module with angular momentum transport *IAMT, *IAMROT
    ! 3. special rhstype with angular momentum conservation
    ! Only the IAMROT physics module and rotframe source module allow
    ! for a different center of rotation than (/0, 0/)

    ! angular velocity of the rotating reference frame
    CALL GetAttr(config, "omega", this%omega, 0.0)

    ! center of rotation
    this%rotcent = (/ 0., 0./)
    CALL GetAttr(config, "rotation_center", this%rotcent, this%rotcent)

    ! basic mesh initialization

    ! get geometrical scale factors for all cell positions
    ! bary center values are overwritten below
    CALL this%geometry%ScaleFactors(this%curv,this%hx,this%hy,this%hz)

    ! get square root of determinant of the metric
    ! bary center values are overwritten below
    this%sqrtg = this%hx*this%hy*this%hz

    ! allocate memory for cartesian positions
    !> Should not exist anymore in 3D coordinates
!    IF (this%Geometry%GetType().EQ.BIANGLESPHERICAL) THEN
!      !> \todo: insert passage in mesh_common.f90 (ll. 264), problem: mixed dimensions in
!      !! bianglespherical for curvilinear and cartesian coordinates -> no generality
!      !! see also geometry_generic
!      CALL this%AllocateMesharray(this%cart,3)
!      CALL this%AllocateMesharray(this%posvec,3)
!    ELSE
      CALL this%AllocateMesharray(this%cart)
      CALL this%AllocateMesharray(this%posvec)
!    END IF
    this%bccart => this%RemapBounds(this%cart%bcenter)
    this%fccart => this%RemapBounds(this%cart%faces)
    this%ccart  => this%RemapBounds(this%cart%corners)

    ! compute cartesian coordinates for all cell positions (center,bcenter,faces,corners)
    CALL this%Geometry%Convert2Cartesian(this%curv,this%cart)

    ! allocate memory for radius
    CALL this%AllocateMesharray(this%radius)

    ! get radii and position vectors for all cell positions
    CALL this%geometry%Radius(this%curv,this%radius)
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
        !CALL SetAttr(IO,"corners",this%ccart(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX,:,:))
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
        IF (ASSOCIATED(this%dlx)) &
           CALL SetAttr(IO,"dlx",this%dlx(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%dly)) &
           CALL SetAttr(IO,"dly",this%dly(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
        IF (ASSOCIATED(this%dlz)) &
           CALL SetAttr(IO,"dlz",this%dlz(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))
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
    IF((writefields.EQ.1).AND.ASSOCIATED(this%volume)) &
        CALL SetAttr(IO,"volume",this%volume(this%IMIN:this%IMAX,this%JMIN:this%JMAX,this%KMIN:this%KMAX))

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
              (y.GE.this%curv%faces(imin,jmin,kmin,1,2).AND.x.LE.this%curv%faces(imax,jmax,kmax,2,2)).AND. &
              (z.GE.this%curv%faces(imin,jmin,kmin,1,3).AND.x.LE.this%curv%faces(imax,jmax,kmax,2,3))) THEN
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
  SUBROUTINE InitMesh_parallel(this, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base), INTENT(INOUT)    :: this   !< \param [in,out] this all mesh data
    TYPE(Dict_TYP),POINTER             :: config
    !------------------------------------------------------------------------!
    LOGICAL, DIMENSION(NDIMS) :: remain_dims, periods = .FALSE.
    INTEGER                            :: num,rem
    INTEGER                            :: ierror
    INTEGER                            :: i,j,k
    INTEGER                            :: worldgroup,newgroup
    INTEGER, DIMENSION(1)              :: rank0in, rank0out
    INTEGER, DIMENSION(NDIMS)          :: coords
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ranks
    !------------------------------------------------------------------------!

    ! create a cartesian topology of processes
    ! 1. balance number of processes per direction
    this%dims(1)=this%GetNumProcs()
    this%dims(2)=1
    this%dims(3)=1
    ! account for vector length of vector CPUs
    ! \todo Comment this in again when parallelization is working!
!    CALL CalculateDecomposition(this%INUM,this%JNUM,this%KNUM,this%GNUM, &
!         this%dims(1),this%dims(2),this%dims(3))
    IF (this%dims(2).LE.0) THEN
       CALL this%Error("InitMesh_parallel","Domain decomposition algorithm failed.")
    END IF

    ! Check if the user set the decomposition dims himself and override the
    ! automatic settings
    CALL GetAttr(config, "decomposition", this%dims, this%dims(:))

    ! If a dimension equals -1, replace it with the total number of processors
    ! This makes it easy to define pure annular ring decompositions in polar
    ! coordinates like as "decomposition" / (/ -1, 1, 1 /)
    WHERE(this%dims.EQ.-1) this%dims=this%GetNumProcs()

    IF((this%dims(1).LE.0).OR.(this%dims(2)).LE.0.OR.(this%dims(3).LE.0).OR.&
       (this%dims(1)*this%dims(2)*this%dims(3).NE.this%GetNumProcs())) THEN
        CALL this%Error("InitMesh_parallel","Invalid user-defined MPI domain "&
            //"decomposition with key='/mesh/decomposition'")
    END IF

    ! 2. create the cartesian communicator
    ! IMPORTANT: disable reordering of nodes
    CALL MPI_Cart_create(MPI_COMM_WORLD,this%NDIMS,this%dims,periods,.TRUE., &
         this%comm_cart,ierror)

    ! 3. inquire and save the own position
    CALL MPI_Cart_get(this%comm_cart,this%NDIMS,this%dims,periods,this%mycoords,ierror)

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


!  !> return the best partitioning of processes
!  !!
!  !! pi x pj x pk for a given mesh with resolution ni x nj x nk with
!  !! ghost number of ghost cells GINUM, GJNUM, GKNUM
!  !! pi : total number of processes (<MAXNUM  see module "factors")
!  !! pj, pk = 1 initially
!  !! pk return 0, for erroneous input
!  !! \todo NOT VERIFIED: Not running in 3D
!  SUBROUTINE CalculateDecomposition(ni,nj,nk,ginum,gjnum,gknum,pi,pj,pk)
!    USE factors
!    IMPLICIT NONE
!    !------------------------------------------------------------------------!
!    INTEGER, INTENT(IN)    :: ni,nj,nk,ginum,gjnum,gknum
!    INTEGER, INTENT(INOUT) :: pi,pj,pk
!    !------------------------------------------------------------------------!
!    ! return immediatly for malformed input
!    IF (((pi.LT.2).AND.(pj.AND.pk.NE.1)).OR.(pj.GT.MAXNUM) .OR.(nj*nk.LT.pj)) THEN
!       pk = 0
!       RETURN
!    END IF
!    Decompose(ni,nj,nk,ginum,gjnum,gknum,pi,pj,pk)
!
!  CONTAINS
!
!    !> searches for the best domain decomposition
!    !!
!    !! Accounts for the costs due to MPI communication (internal boundaries)
!    !! and optimizes for a given vector length (VECLEN) on vector computers;
!    !! parameters:
!    !!   ni : number of grid cells in first dimension
!    !!   nj : number of grid cells in second dimension
!    !!   nk : number of grid cells in third dimension
!    !!   pi : number of processes  in first dimension (=NumProcs first call)
!    !!   pj : number of processes  in second dimension (=1 at first call)
!    !!   pk : number of processes  in third dimension (=1 at first call)
!    !! \todo NOT VERIFIED: Not running for 3D
!    RECURSIVE SUBROUTINE Decompose(ni,nj,nk,pi,pj,pk,pires,pjres,pkres)
!      IMPLICIT NONE
!      !-------------------------------------------------------------------!
!      INTEGER, INTENT(IN)  :: ni,nj,nk,pi,pj,pk
!      INTEGER, INTENT(OUT) :: pires,pjres,pkres
!      !-------------------------------------------------------------------!
!      INTEGER :: pp,ptot
!      INTEGER :: p1,p2,p3,pinew,pjnew,pknew,piold,pjold,pkold
!      INTEGER :: pfmin,pfnew,pfold
!      INTEGER :: bl,vl,blnew,vlnew
!      REAL    :: bl_gain,vl_gain
!      !-------------------------------------------------------------------!
!      ! measure the costs of the given configuration
!      CALL GetCosts(ni,nj,nk,pi,pj,pk,bl,vl)
!!!$PRINT '(A,2(I7),I12,I4)'," costs: ",pj,pk,bl,vl
!      ! save the configuration
!      piold=pi
!      pjold=pj
!      pkold=pk
!      p1 = pi
!      p2 = pj
!      p3 = pk
!      !> todo from here nothing changed yet
!      ! compute the total number of processes
!      ptot = pi*pj*pk
!      pfmin = GetFactor(pk)      ! get smallest prime factor of pk
!      pfold = 1
!      pp = pi
!      DO ! loop over all prime factors of pj which are larger than
!         ! the smallest prime factor of pk
!         IF (pp.LE.1) EXIT       ! if pj has been reduced to 1
!         ! get smallest prime factor of pp, i.e. pj
!         pfnew = GetFactor(pp)
!         IF ((pfnew.GT.pfmin).AND.(pfmin.NE.1)) EXIT
!         pp = pp/pfnew
!         ! skip multiple prime factors
!         IF (pfnew.NE.pfold) THEN
!            ! create new configuration
!            p1 = p1*pfold/pfnew
!            p2 = p2/pfold*pfnew
!            pfold = pfnew
!            ! get the best configuration possible with p1 x p2 x p3 processes
!            ! by reducing p1 => recursion
!            pjnew = Decompose(ni,nj,nk,p1,p2,p3,pires,pjres,pkres)
!            pknew = ptot/pjnew  ! compute the second factor using the
!                                ! total number of processes
!            CALL GetCosts(ni,nj,nk,pjnew,pknew,blnew,vlnew)
!            bl_gain = bl*(1.0/blnew)             ! smaller is better
!            vl_gain = vlnew*(1.0/vl)             ! larger is better
!!!$PRINT '(4I7)',pjold,pkold,bl,vl
!!!$PRINT '(4I7)',pjnew,pknew,blnew,vlnew
!!!$PRINT *,"--------------------------------------------"
!!!$PRINT '(A,3F7.1)',"              ",bl_gain,vl_gain,bl_gain*vl_gain
!            ! compare new with old configuration
!            IF (vl_gain*bl_gain.GT.1.0) THEN
!               ! new configuration is better
!               piold = pinew
!               pjold = pjnew
!               pkold = pknew
!               bl = blnew
!               vl = vlnew
!            END IF
!         END IF
!      END DO
!      ! return optimized number of processes in the first dimension
!      ! (for the initial configuration pj x pk)
!      pjres=pjold
!    END FUNCTION Decompose
!
!    ! computes the sum of the length of all internal boundaries (bl)
!    ! and the maximal vector length (vl)
!    PURE SUBROUTINE GetCosts(n1,n2,n3,p1,p2,p3,bl,vl)
!      IMPLICIT NONE
!      !------------------------------------------------------------------------!
!      INTEGER, INTENT(IN)  :: n1,n2,n3,p1,p2,p3
!      INTEGER, INTENT(OUT) :: bl,vl
!      !------------------------------------------------------------------------!
!      INTEGER :: num,rem
!      !------------------------------------------------------------------------!
!      ! length of internal boundaries
!      !bl = n1*(p2-1) + n2*(p1-1)
!      bl = n3*n2*(p1-1) + n1*n3*(p2-1)  + n1*n2*(p3-1)
!      ! maximal possible vector length (first array dimension)
!      ! if n1 > VECLEN return the length of the remainder
!      rem = MOD(n1,p1)
!      num = (n1-rem) / p1 + MIN(rem,1) + 2*ginum  ! account for ghost cells
!      IF (num.GT.VECLEN) num=MOD(num,VECLEN)
!      vl  = VECLEN-MOD(ABS(num-VECLEN),VECLEN)
!    END SUBROUTINE GetCosts
!
!  END SUBROUTINE CalculateDecomposition
#endif
!---------------------END PARALLEL-------------------------------------------!

  !> \public allocates memory for scalar mesh array
  SUBROUTINE AllocateMesharrayS(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: this !< \param [in,out] this all mesh data
    TYPE(MArrayS_TYP)           :: ma   !< \param [out] ma scalar mesh array
    !------------------------------------------------------------------------!
    INTEGER                     :: err
    !------------------------------------------------------------------------!
    INTENT(OUT)                 :: ma
    !------------------------------------------------------------------------!
    ALLOCATE(ma%data(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX,(1+1+6+8)),STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("AllocateMesharrayS","Unable to allocate memory!")
    END IF
    ma%center  => this%RemapBounds(ma%data(:,:,:,1))
    ma%bcenter => this%RemapBounds(ma%data(:,:,:,2))
    ma%faces   => this%RemapBounds(ma%data(:,:,:,3:8))
    ma%corners => this%RemapBounds(ma%data(:,:,:,9:16))
  END SUBROUTINE AllocateMesharrayS

  !> \public allocates memory for vector mesh array
  SUBROUTINE AllocateMesharrayV(this,ma,n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: this !< \param [in,out] this all mesh data
    TYPE(MArrayV_TYP)           :: ma   !< \param [out] ma vector mesh array
    INTEGER,OPTIONAL            :: n    !< \param [in] n vector dimension
    !------------------------------------------------------------------------!
    INTEGER                     :: err,n_def
    !------------------------------------------------------------------------!
    INTENT(OUT)                 :: ma
    !------------------------------------------------------------------------!
    IF (PRESENT(n)) THEN
       n_def = n
    ELSE
       ! dimension defaults to 2
       n_def = 3
    END IF
    ALLOCATE(ma%data(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX,(1+1+6+8),n_def),STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("AllocateMesharrayV","Unable to allocate memory!")
    END IF
    ma%center  => this%RemapBounds(ma%data(:,:,:,1,:))
    ma%bcenter => this%RemapBounds(ma%data(:,:,:,2,:))
    ma%faces   => this%RemapBounds(ma%data(:,:,:,3:8,:))
    ma%corners => this%RemapBounds(ma%data(:,:,:,9:16,:))
  END SUBROUTINE AllocateMesharrayV


  !> \public allocates memory for tensor mesh array
  SUBROUTINE AllocateMesharrayT(this,ma,n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: this !< \param [in,out] this all mesh data
    TYPE(MArrayT_TYP)           :: ma   !< \param [out] ma tensor mesh array
    INTEGER,OPTIONAL            :: n    !< \param [in] n tensor dimension
    !------------------------------------------------------------------------!
    INTEGER                     :: err,n_def
    !------------------------------------------------------------------------!
    INTENT(OUT)                 :: ma
    !------------------------------------------------------------------------!
    IF (PRESENT(n)) THEN
       n_def = n
    ELSE
       ! dimension defaults to 2
       n_def = 3
    END IF
    ALLOCATE(ma%data(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,this%KGMIN:this%KGMAX,(1+1+6+8),n_def,n_def),STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("AllocateMesharrayT","Unable to allocate memory!")
    END IF
    ma%center  => this%RemapBounds(ma%data(:,:,:,1,:,:))
    ma%bcenter => this%RemapBounds(ma%data(:,:,:,2,:,:))
    ma%faces   => this%RemapBounds(ma%data(:,:,:,3:8,:,:))
    ma%corners => this%RemapBounds(ma%data(:,:,:,9:16,:,:))
  END SUBROUTINE AllocateMesharrayT


  !> \public deallocates memory of scalar mesh array
  SUBROUTINE DeallocateMesharrayS(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: this !< \param [in,out] this all mesh data
    TYPE(MArrayS_TYP)           :: ma
    !------------------------------------------------------------------------!
    INTEGER                     :: err
    !------------------------------------------------------------------------!
    INTENT(OUT)                 :: ma
    !------------------------------------------------------------------------!
    DEALLOCATE(ma%data,STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("DeallocateMesharrayS","Unable to deallocate memory!")
    END IF
    NULLIFY(ma%center,ma%bcenter,ma%faces,ma%corners)
  END SUBROUTINE DeallocateMesharrayS

  !> \public deallocates memory of vector mesh array
  SUBROUTINE DeallocateMesharrayV(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: this !< \param [in,out] this all mesh data
    TYPE(MArrayV_TYP)           :: ma
    !------------------------------------------------------------------------!
    INTEGER                     :: err
    !------------------------------------------------------------------------!
    INTENT(OUT)                 :: ma
    !------------------------------------------------------------------------!
    DEALLOCATE(ma%data,STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("DeallocateMesharrayV","Unable to deallocate memory!")
    END IF
    NULLIFY(ma%center,ma%bcenter,ma%faces,ma%corners)
  END SUBROUTINE DeallocateMesharrayV

  !> \public deallocates memory of tensor mesh array
  SUBROUTINE DeallocateMesharrayT(this,ma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(IN) :: this !< \param [in,out] this all mesh data
    TYPE(MArrayT_TYP)           :: ma
    !------------------------------------------------------------------------!
    INTEGER                     :: err
    !------------------------------------------------------------------------!
    INTENT(OUT)                 :: ma
    !------------------------------------------------------------------------!
    DEALLOCATE(ma%data,STAT=err)
    IF (err.NE.0) THEN
       CALL this%Error("DeallocateMesharrayT","Unable to deallocate memory!")
    END IF
    NULLIFY(ma%center,ma%bcenter,ma%faces,ma%corners)
  END SUBROUTINE DeallocateMesharrayT


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
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),INTENT(INOUT) :: this   !< \param [in,out] this all mesh data
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER :: i,ierror
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("CloseMesh","not initialized")

    DO i=1,4
       IF (this%comm_boundaries(i).NE.MPI_COMM_NULL) &
          CALL MPI_Comm_free(this%comm_boundaries(i),ierror)
    END DO
#endif
    CALL this%DeallocateMesharray(this%curv)
    CALL this%DeallocateMesharray(this%hx)
    CALL this%DeallocateMesharray(this%hy)
    CALL this%DeallocateMesharray(this%hz)
    CALL this%DeallocateMesharray(this%sqrtg)
    CALL this%DeallocateMesharray(this%cxyx)
    CALL this%DeallocateMesharray(this%cxzx)
    CALL this%DeallocateMesharray(this%cyxy)
    CALL this%DeallocateMesharray(this%cyzy)
    CALL this%DeallocateMesharray(this%czxz)
    CALL this%DeallocateMesharray(this%czyz)

    DEALLOCATE(this%volume,this%dxdydV,this%dydzdV,this%dzdxdV,this%dlx,this%dly,this%dlz)
    IF(ASSOCIATED(this%rotation)) &
        DEALLOCATE(this%rotation)


    CALL this%DeallocateMesharray(this%cart)
    CALL this%DeallocateMesharray(this%radius)
    CALL this%DeallocateMesharray(this%posvec)

    IF (ASSOCIATED(this%rotation)) DEALLOCATE(this%rotation)

    DEALLOCATE(this%Geometry)
  END SUBROUTINE Finalize


END MODULE mesh_base_mod
