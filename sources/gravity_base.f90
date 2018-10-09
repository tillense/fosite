!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_base.f03                                                  #
!#                                                                           #
!# Copyright (C) 2014-2018                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup gravity
!! - general parameters of gravity group as key-values
!! \key{gtype,INTEGER,Type of gravity source}
!! \key{energy,INTEGER,Add source terms to energy equation?}
!! \key{output/accel,INTEGER,enable(=1) output of acceleration}
!! \key{output/height,INTEGER,enable(=1) output of disc height}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!! \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief generic gravity terms module providing functionaly common to all
!! gravity terms
!----------------------------------------------------------------------------!
MODULE gravity_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE boundary_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, ABSTRACT, EXTENDS(logging_base) :: gravity_base
     !> \name Variables
     CLASS(logging_base), ALLOCATABLE :: gravitytype  !< type of gravity term
     CLASS(gravity_base), POINTER     :: next => null() !< next gravity in list
     REAL                             :: time         !< last update
     REAL, POINTER                    :: mass2        !< 2nd mass for binaries
     REAL                             :: excent       !< excentricity
     REAL                             :: semaaxis     !< semi major axis
     REAL                             :: period       !< period of binaries
     REAL                             :: eps1,eps2    !< softening parameter
     REAL                             :: alphag,alphah!< viscosity parameter
     LOGICAL                          :: CALCALPHA    !< true if alpha output is on
     !> time when the pointmass is fully switched on
     REAL                             :: switchon
     REAL                             :: switchon2    !< same for secondary in binary systems
     REAL                             :: scaling      !< scaling for switchon procedure
     !> angular velocity of rotating reference frame
     REAL                             :: omega_rot
     !> time of a orbital period at the inner and outer boundaries
     REAL, DIMENSION(2)               :: tau
     INTEGER                          :: outbound     !< outflow boundary
     REAL, DIMENSION(:,:,:,:), POINTER   :: accel        !< acceleration
     REAL, DIMENSION(:,:,:), POINTER     :: radius,radius3!< distance to origin
     REAL, DIMENSION(:,:,:), POINTER     :: r_sec        !<    and to secondary point mass
     REAL, DIMENSION(:,:,:,:), POINTER   :: posvec_sec   !<   secondary to all cell bary centers
     REAL, DIMENSION(:,:,:,:), POINTER   :: fr_sec
     REAL, DIMENSION(:,:,:,:,:), POINTER :: fposvec_sec
     REAL, DIMENSION(:,:), POINTER      :: joff, jrem   !< shifting indices (in SB)
     REAL, DIMENSION(:,:,:), POINTER    :: den_ip       !< interpolated density
     REAL, DIMENSION(:,:,:), POINTER    :: omega        !< angular velocity
     REAL, DIMENSION(:,:,:,:), POINTER  :: omega2       !< Omega Kepler squared
     REAL, DIMENSION(:,:,:,:), POINTER  :: gposvecr3    !< = GN*x/radius**3
     REAL, DIMENSION(:,:,:), POINTER    :: cellmass     !< rho*dV
     REAL, DIMENSION(:,:,:), POINTER    :: enclmass     !< enclosed mass
     REAL, DIMENSION(:,:,:), POINTER    :: rho_c        !< disk density in equatorial plane
     REAL, DIMENSION(:,:,:), POINTER    :: inv2PIr2     !< = 1/ (2*PI*r**2)
     REAL, DIMENSION(:,:,:), POINTER    :: tmp,tmp2,tmp3!< temp arrays
     REAL, DIMENSION(:,:,:), POINTER    :: dVpart       !< partial cell volumes (monopole approx.)
     LOGICAL                         :: addtoenergy
     LOGICAL                         :: update_disk_height !< enable/disable computation of disk scale height
     INTEGER                         :: timeid       !<    update of this id?

     !> \name todo this part is of the multigrid method
     !! which was not understandable in fosite 2D anymore
     !!#### Poisson module
!     TYPE(Common_TYP)                :: poissontype   !< type of source term
!     TYPE(Multipole_TYP)             :: multipole     !< multipole expansion
!     TYPE(Grid_TYP), POINTER         :: grid(:)       !< coarse grids
!     REAL                            :: MAXRESIDNORM  !< max error of residuum
!     !> type of relaxation method
!     INTEGER                         :: RELAXTYPE
!     INTEGER                         :: NGRID         !< number of grids
!     INTEGER                         :: NPRE,NPOST    !< pre and post smoothing
!     INTEGER                         :: NMAXCYCLE     !< max of iterations
!     INTEGER                         :: MINRES        !< min resolution
     REAL, DIMENSION(:,:,:), POINTER    :: phi           !< potential
     REAL, DIMENSION(:,:,:,:), POINTER  :: pot           !< general potential
     REAL, DIMENSION(:,:,:,:), POINTER  :: pot_prim      !< potential second component
     REAL, DIMENSION(:,:,:,:), POINTER  :: pot_sec       !< potential second component
     REAL, DIMENSION(:,:,:,:,:), POINTER :: mpot        !< multiple potentials
     REAL, DIMENSION(:,:,:), POINTER    :: rho_ext       !< disk density in equatorial plane
                                                      !< of the external potentials
     REAL, DIMENSION(:,:,:,:,:), POINTER :: mrho_ext    !< multiple rho_ext
     REAL, DIMENSION(:,:,:), POINTER  :: invheight2   !< 1/h**2
     REAL, DIMENSION(:,:,:), POINTER  :: height,h_ext !< disk scale height
     INTEGER                          :: n             !< number of potentials
     REAL,DIMENSION(:,:),POINTER      :: s0, sdelta    !< ramp fn: s0 + sdelta*t
     REAL,DIMENSION(:),POINTER        :: lastfac       !< last switchon factors
     INTEGER, DIMENSION(4)            :: Boundary      !< boundary condition
     LOGICAL                          :: DIRICHLET     !< true if min ONE bound.
                                                      !< boundary cond. is
#if defined(HAVE_FFTW) || defined(HAVE_FFTKEISAN)
    INTEGER, DIMENSION(:), POINTER    :: sizes
    INTEGER                           :: MNUM          !< number of modes
    REAL,DIMENSION(:), POINTER        :: kx            !< wave numbers for FFT (x)
    REAL,DIMENSION(:), POINTER        :: ky            !< wave numbers for FFT (y)
    REAL                              :: Lx, Ly
#endif
#if defined(HAVE_FFTW)
    !> \name
    !!#### spectral poisson solver
    !> plan for real to complex fourier transforms
    TYPE(C_PTR)                      :: plan_r2c
    !> plan for complex to real fourier transforms
    TYPE(C_PTR)                      :: plan_c2r
    TYPE(C_PTR)                      :: pFdensity, pFphi
    REAL(C_DOUBLE), DIMENSION(:,:), POINTER &
                                     :: Fdensity,Fphi,block
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:), POINTER &
                                     :: cFdensity, cFphi, cblock
    !> Important precalculated matrix - fourier transformed I
    REAL(C_DOUBLE), DIMENSION(:,:,:), POINTER &
                                     :: FI
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), POINTER &
                                     :: cFI
    TYPE(C_PTR)                      :: p_FI
    !> \name
    !!#### spectral poisson solver shearing box
    REAL(C_DOUBLE), POINTER          :: mass2D(:,:)    !< temporary variable
    COMPLEX(C_DOUBLE_COMPLEX), POINTER &
                                     :: Fmass2D(:,:)   !< temporary variable
    REAL, DIMENSION(:,:), POINTER    :: Fmass2D_real   !< temporary variable
    INTEGER(C_INTPTR_T)              :: local_joff
#elif defined(HAVE_FFTKEISAN)
    REAL, DIMENSION(:,:,:), POINTER  :: FI
    REAL, DIMENSION(:,:), POINTER    :: Fdensity,Fphi,block
    INTEGER                          :: local_joff
#endif
#if defined(HAVE_FFTW) && defined(PARALLEL)
    INTEGER(C_INTPTR_T)              :: C_INUM, C_JNUM
    INTEGER(C_INTPTR_T)              :: alloc_local, local_JNUM
    TYPE(C_PTR)                      :: mass2D_pointer
    TYPE(C_PTR)                      :: Fmass2D_pointer
#elif defined(HAVE_FFTKEISAN) && defined(PARALLEL)
    COMPLEX, DIMENSION(:,:), POINTER :: mpi_sendbuf
    COMPLEX, DIMENSION(:,:), POINTER :: mpi_recvbuf
#endif
#ifdef HAVE_FFTW_LEGACY
    !> plan for real to complex fourier transforms
    INTEGER*8                        :: plan_r2c
    !> plan for complex to real fourier transforms
    INTEGER*8                        :: plan_c2r
#endif
#if defined(HAVE_FFTW_LEGACY) || defined(HAVE_FFTKEISAN)
#if defined(HAVE_FFTKEISAN) && defined(PARALLEL)
    COMPLEX, DIMENSION(:,:), POINTER :: mass2D         !< temporary variable
#else
    REAL, DIMENSION(:,:), POINTER    :: mass2D         !< temporary variable
#endif
    COMPLEX, POINTER                 :: Fmass2D(:,:)   !< temporary variable
    REAL, DIMENSION(:,:), POINTER    :: Fmass2D_real   !< temporary variable
#endif
    INTEGER                          :: green
    REAL                             :: sigma, ecut
!    REAL, DIMENSION(:), POINTER      :: height
    INTEGER, POINTER                 :: mcut
    !> local IMAX, INUM
    INTEGER                          :: IMAX, INUM
    !> \name Variables in Parallel Mode
#ifdef PARALLEL
    INTEGER                          :: error         !< MPI error flag
    REAL,DIMENSION(:,:),POINTER      :: sbuf1,sbuf2,rbuf1,rbuf2
    !> displacment and length of domain
    INTEGER,DIMENSION(:),POINTER     :: displ, num
#endif
    INTEGER, DIMENSION(:),POINTER :: relaxcount

  CONTAINS

    PROCEDURE :: InitGravity
    PROCEDURE :: SetOutput
    PROCEDURE (InfoGravity),     DEFERRED :: InfoGravity
    PROCEDURE :: GravitySources
    PROCEDURE :: UpdateGravity
    PROCEDURE (UpdateGravity_single), DEFERRED :: UpdateGravity_single
    PROCEDURE :: CloseGravity
    PROCEDURE :: CalcDiskHeight
    PROCEDURE (CalcDiskHeight_single), DEFERRED :: CalcDiskHeight_single
    PROCEDURE :: GetGravityPointer

  END TYPE gravity_base
  ABSTRACT INTERFACE
    SUBROUTINE InfoGravity(this,Mesh)
      IMPORT gravity_base, mesh_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(IN)    :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
    END SUBROUTINE
    SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
      IMPORT gravity_base, mesh_base, physics_base, fluxes_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                           INTENT(IN)    :: pvar
      REAL,                INTENT(IN)    :: time,dt
    END SUBROUTINE
    SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
      IMPORT gravity_base, mesh_base, physics_base
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(gravity_base), INTENT(INOUT) :: this
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN)    :: Physics
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                           INTENT(IN)    :: pvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(IN)    :: bccsound
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(OUT)   :: h_ext
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                           INTENT(INOUT) :: height
    END SUBROUTINE
  END INTERFACE
  ! flags for source terms
  INTEGER, PARAMETER :: POINTMASS        = 1
  INTEGER, PARAMETER :: POINTMASS_BINARY = 2
  INTEGER, PARAMETER :: MONOPOL          = 3
  INTEGER, PARAMETER :: MULTIGRID        = 4
  INTEGER, PARAMETER :: SPECTRAL         = 5
  INTEGER, PARAMETER :: POTENTIAL        = 6
  INTEGER, PARAMETER :: SBOXSPECTRAL     = 7
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       gravity_base, &
       ! constants
       POINTMASS, POINTMASS_BINARY, MONOPOL, &
       MULTIGRID, SPECTRAL, POTENTIAL, SBOXSPECTRAL
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGravity(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base),  INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    INTEGER                             :: stype
    TYPE(Dict_TYP),POINTER              :: config,IO
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: dir,src,IOsrc
    INTEGER                :: gtype,err,i
    LOGICAL                :: external_potential = .FALSE.
    !------------------------------------------------------------------------!
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DIM), &
             this%pot(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4), &
             STAT=err)
    IF (err.NE.0) CALL this%Error("InitGravity", "Unable allocate memory!")

    NULLIFY(this%height,this%h_ext,this%invheight2)
    this%accel(:,:,:,:) = 0.

    ! Add source terms to energy equation?
    ! Set this to zero, if a potential is defined in physics_euler2Diamt
    CALL GetAttr(config, "energy", i, 1)
    IF(i.EQ.0) THEN
      this%addtoenergy = .FALSE.
    ELSE
      this%addtoenergy = .TRUE.
    END IF

    CALL this%SetOutput(Mesh,config,IO)

    ! enable update of disk scale height if requested
    CALL GetAttr(config, "update_disk_height", i, 0)
    IF (i.EQ.1) THEN
       IF (Physics%DIM.EQ.2) THEN
          this%update_disk_height = .TRUE.
          IF (.NOT.ASSOCIATED(this%height)) &
             ALLOCATE(this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      STAT=err)
          IF (err.EQ.0) &
             ALLOCATE(this%h_ext(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%invheight2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      STAT=err)
          IF (err.NE.0) CALL this%Error("InitGravity","Memory allocation failed!")
          IF (external_potential) &
             CALL this%Warning("InitGravity", &
             "calculation of disk height not fully supported for external potential")
       ELSE
          CALL this%Error("InitGravity", "DiskHeight is only supported in 2D")
       END IF
    ELSE
       this%update_disk_height = .FALSE.
    END IF

    ! reset start value for time variable
    this%time = -1.0
    this%timeid = 0


!    ! add new gravity term to beginning of
!    ! list of gravity terms
!    IF (.NOT.ASSOCIATED(this)) THEN
!       this => newgrav
!       NULLIFY(this%next)
!    ELSE
!       tmpgrav => this
!       this => newgrav
!       this%next => tmpgrav
!    END IF
!    !- END COMMON PART -!

    CALL this%Info(" GRAVITY--> gravity term:      " // this%GetName())
    CALL this%InfoGravity(Mesh)
  END SUBROUTINE InitGravity

  SUBROUTINE GravitySources(this,Mesh,Physics,Fluxes,time,dt,pvar,cvar,gterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), TARGET        :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    REAL                               :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN)    :: cvar,pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(OUT)   :: gterm
    !------------------------------------------------------------------------!
    ! update acceleration of all gravity sources
    CALL this%UpdateGravity(Mesh,Physics,Fluxes,time,dt,pvar)

    ! gravitational source terms
    CALL Physics%ExternalSources(Mesh,this%accel,pvar,cvar,gterm)

    ! Set src term in energy equation to zero, if it is handeled in the physics
    ! module
    IF((.NOT.this%addtoenergy).AND.(Physics%ENERGY.GT.0)) THEN
      gterm(:,:,:,Physics%ENERGY) = 0.
    END IF
  END SUBROUTINE GravitySources

  !> Loops over all gravity terms and yiels acceleration and optionally the potential
  SUBROUTINE UpdateGravity(this,Mesh,Physics,Fluxes,time,dt,pvar)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), TARGET        :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(INOUT) :: Physics
    CLASS(fluxes_base),  INTENT(IN)    :: Fluxes
    REAL                               :: time,dt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER       :: gravptr
    !------------------------------------------------------------------------!
    ! reset gterm
    this%pot(:,:,:,:) = 0.
    this%accel(:,:,:,:) = 0.

    ! go through all gravity terms in the list
    gravptr => this
    DO WHILE (ASSOCIATED(gravptr))
      ! call specific subroutine
!CDIR IEXPAND
      CALL gravptr%UpdateGravity_single(Mesh,Physics,Fluxes,pvar,time,dt)

      ! add to the sources
      this%accel(:,:,:,:) = this%accel(:,:,:,:) + gravptr%accel(:,:,:,:)
      IF(ASSOCIATED(gravptr%pot)) &
        this%pot(:,:,:,:) = this%pot(:,:,:,:) + gravptr%pot(:,:,:,:)

      ! next source term
      gravptr => gravptr%next
    END DO

    ! update disk scale height if requested
    IF (this%update_disk_height) CALL this%CalcDiskHeight(Mesh,Physics,pvar)

  END SUBROUTINE UpdateGravity

  SUBROUTINE CalcDiskHeight(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), TARGET, INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER :: grav_ptr,selfgrav_ptr => null()
    LOGICAL                      :: has_external_potential = .FALSE.
    !------------------------------------------------------------------------!
    ! reset inverse scale height^2
    this%invheight2(:,:,:) = 0.0
    ! go through all gravity terms in the list
    grav_ptr => this
    ! \todo Most probably here it is necessar to make a selection for all types
    DO WHILE(ASSOCIATED(grav_ptr))
      ! call specific subroutine
!CDIR IEXPAND
      CALL grav_ptr%CalcDiskHeight_single(Mesh,Physics,pvar,Physics%bccsound,this%h_ext,this%height)
      this%invheight2(:,:,:) = this%invheight2(:,:,:) + 1./this%h_ext(:,:,:)**2
      has_external_potential = .TRUE.
!     \todo add this routines to underlying sub-classes
!      CASE(MONOPOL,SPECTRAL,SBOXSPECTRAL)
!        ! remember self-gravity pointer
!        selfgrav_ptr => grav_ptr

      ! next gravity term
      grav_ptr => grav_ptr%next
    END DO

    ! self-gravity of the disk needs special treatment
    IF (ASSOCIATED(selfgrav_ptr)) THEN
      IF (has_external_potential) THEN
        ! compute the resultant height due to all external gravitational forces
        this%h_ext(:,:,:) = 1./SQRT(this%invheight2(:,:,:))
      END IF
      CALL selfgrav_ptr%CalcDiskHeight_single(Mesh,Physics,pvar,Physics%bccsound, &
                                  this%h_ext,this%height)
    ELSE
      ! non-selfgravitating disk
      this%height(:,:,:) = 1./SQRT(this%invheight2(:,:,:))
    END IF

  END SUBROUTINE CalcDiskHeight

  SUBROUTINE SetOutput(this,Mesh,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    TYPE(Dict_TYP),      POINTER       :: config,IO
    !------------------------------------------------------------------------!
    CHARACTER(LEN=1) :: xyz(3) = (/"x","y","z"/)
    INTEGER          :: valwrite,err,k
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "output/accel", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       DO k=1,SIZE(this%accel,3)
          CALL SetAttr(IO, (xyz(k) // "accel"),&
             this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,k))
       END DO
    END IF

    CALL GetAttr(config, "output/potential", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "potential",&
         this%pot(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1))

    CALL GetAttr(config, "output/height", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       CALL SetAttr(config, "update_disk_height", 1)
       IF (.NOT.ASSOCIATED(this%height)) THEN
          ALLOCATE(this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                   STAT=err)
          IF (err.NE.0) CALL this%Error("SetOutput", &
                                   "Memory allocation failed for this%height!")
       END IF
       CALL SetAttr(IO, "height", &
         this%height(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    END IF

  END SUBROUTINE SetOutput

  FUNCTION GetGravityPointer(list,stype) RESULT(gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base), TARGET  :: list
    INTEGER, INTENT(IN)          :: stype
    CLASS(gravity_base), POINTER :: gp
    !------------------------------------------------------------------------!
    gp => list
    DO
       IF (ASSOCIATED(gp).EQV..FALSE.) EXIT
!CDIR IEXPAND
       IF (gp%GetType().EQ.stype) RETURN
       gp => gp%next
    END DO
  END FUNCTION GetGravityPointer

  SUBROUTINE CloseGravity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_base) :: this
    !------------------------------------------------------------------------!
    CLASS(gravity_base), POINTER :: gravptr,ptemp
    !------------------------------------------------------------------------!
    ! release temporary/global storage
    DEALLOCATE(this%accel,this%pot)
    ! TODO Wait for Lars' solution
!    IF (this%update_disk_height) DEALLOCATE(this%height,this%h_ext,this%invheight2)
!    gravptr => this%glist
!    ! call deallocation procedures for all source terms
!    DO
!       IF (.NOT.ASSOCIATED(gravptr)) EXIT
!       IF (.NOT.Initialized(gravptr)) &
!            CALL Error(gravptr,"CloseGravity","not initialized")
!       ! call specific deconstructor
!!CDIR IEXPAND
!       ! deallocate source term structure
!       ptemp=>gravptr
!       gravptr=>gravptr%next
!       DEALLOCATE(ptemp)
!    END DO
  END SUBROUTINE CloseGravity

END MODULE gravity_base_mod
