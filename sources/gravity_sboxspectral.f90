!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_sboxspectral.f90                                          #
!#                                                                           #
!# Copyright (C) 2015-2019                                                   #
!# Jannes Klee <jklee@astrophysik.uni-kiel.de>                               #
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
!> \addtogroup gravity
!! - parameters of \link gravity_sboxspectral_mod gravity_sboxspectral \endlink
!!   as key-value pairs
!!  \key{order,INTEGER,order for the approximation of the gradient of the
!!       potential \f$ \in \\{2\,4\\}\f$ ,2}
!!  \key{print_plans,INTEGER,enable=1/disable=0 output of FFTW plans,0}
!!  \key{output/phi,INTEGER,(enable=1) output of gravitational potential}
!!  \key{output/accel_x,INTEGER,(enable=1) output of accel in x-direction}
!!  \key{output/accel_y,INTEGER,(enable=1) output of accel in y-direction}
!!  \key{output/mass2D,INTEGER,(enable=1) output of density that will be
!!                              fourier transformed}
!!  \key{output/Fmass2D,INTEGER,(enable=1) output of density in spectral
!!                              space}
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!! \author Tobias Illenseer
!!
!! \brief Poisson solver via spectral methods within the shearingbox.
!!
!! \extends gravity_common
!! \ingroup gravity
!!
!! Program and data initialization for a shearing-box simulation
!!
!! There are currently two different implementations. Please be aware that this module
!! requires a certain domain decomposition in pencils and cannot be chosen
!! arbitrarily.
!!
!! 1. FFTW on one core
!! 2. FFTW distributed parallel
!
!! \todo Remove unnecessary allocations like ip_den, etc. At the moment there
!!       are certainly possibilities to save memory, because for every task
!!       a new array is allocated.
!
!! \attention The parallel version with this module needs a pencil decomposition
!!            of Mesh%INUM X Mesh%JNUM/nprocs for each process, in order to use
!!            the FFT-method with period boundaries in one direction.
!!
!! References:
!!
!! \cite gammie2001 , \cite gammiecode, \cite frigo2005, \cite mathkeisan2018
!!
!! \extends gravity_spectral
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_sboxspectral_mod
  USE gravity_base_mod
  USE gravity_spectral_mod, ONLY : gravity_spectral
  USE boundary_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE logging_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE common_dict
  USE fftw
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
  REAL, PARAMETER              :: SQRTTWOPI &
    = 2.50662827463100050241576528481104525300698674
  TYPE C_DOUBLE_PTR2D_TYP
    REAL(C_DOUBLE), DIMENSION(:,:), POINTER :: ptr2D
  END TYPE C_DOUBLE_PTR2D_TYP
  TYPE C_DOUBLE_COMPLEX_PTR2D_TYP
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:), POINTER :: ptr2D
  END TYPE C_DOUBLE_COMPLEX_PTR2D_TYP

  TYPE, EXTENDS(gravity_spectral) :: gravity_sboxspectral
#ifdef HAVE_FFTW
    !> \name
    !!#### spectral poisson solver shearing box
    TYPE(C_DOUBLE_PTR2D_TYP), DIMENSION(:), ALLOCATABLE &
                                     :: density        !< z-layered density field
    TYPE(C_DOUBLE_COMPLEX_PTR2D_TYP), DIMENSION(:), ALLOCATABLE &
                                     :: FTdensity      !< z-layered fourier transformed density
    REAL(C_DOUBLE), DIMENSION(:,:), POINTER, CONTIGUOUS &
                                     :: mass2D, &      !< 2D mass density within a z-layer
                                        Kxy2           !< length of wave vectors squared
    COMPLEX(C_DOUBLE_COMPLEX), POINTER, CONTIGUOUS &
                                     :: Fmass2D(:,:),&   !< fourier transformed 2D mass density
                                        Fmass3D(:,:,:),& !< fourier transformed 3D mass density
                                        Fsum3D(:,:,:), & !< sum used for reconstruction
                                                         !<   of transformed potential (3D only)
                                        qk(:,:,:)        !< weight factors (3D only)
    REAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS &
                                     :: Fmass3D_real => null(), & !< just for output
                                        den_ip         !< interpolated 3D density
    INTEGER(C_INTPTR_T)              :: local_joff
    REAL,DIMENSION(:), POINTER       :: kx,ky          !< x/y wave numbers for FFT
    REAL                             :: shiftconst     !< constant for shift
    REAL                             :: maxKxy2        !< max
    REAL, DIMENSION(:), POINTER      :: joff, jrem     !< shifting indices (in SB)
    INTEGER                          :: order
#ifdef PARALLEL
    INTEGER(C_INTPTR_T)              :: C_INUM, C_JNUM
    INTEGER(C_INTPTR_T)              :: alloc_local, local_JNUM
    TYPE(C_PTR)                      :: mass2D_pointer
    TYPE(C_PTR)                      :: Fmass2D_pointer
#endif
#endif

  CONTAINS

    PROCEDURE :: InitGravity_sboxspectral
    PROCEDURE :: UpdateGravity_single
!     PROCEDURE :: InfoGravity
    PROCEDURE :: SetOutput
    PROCEDURE :: CalcDiskHeight_single
    PROCEDURE :: Finalize
#ifdef HAVE_FFTW
    PROCEDURE :: CalcPotential
    PROCEDURE :: FFT_Forward
    PROCEDURE :: FFT_Backward
#endif
    PROCEDURE :: FieldShift
  END TYPE

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       gravity_sboxspectral
  !--------------------------------------------------------------------------!

  CONTAINS

  !> \public Constructor of gravity sboxspectral module
  !!
  !! This subroutine reads the necessary config data for the the specral
  !! gravity solver within the shearingbox. It initializes the gravity type
  !! and various mesh data arrays. Some of those are marked for output.
  !! The possible combinations are
  !! 1. FFTW - serial, parallel
  SUBROUTINE InitGravity_sboxspectral(this,Mesh,Physics,config,IO)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    TYPE(Dict_TYP),      POINTER    :: config,IO
    !------------------------------------------------------------------------!
    INTEGER             :: gravity_number
    INTEGER             :: err, valwrite, i,j,k
#if defined(HAVE_FFTW) && defined(PARALLEL)
    INTEGER(C_INTPTR_T) :: C_INUM,C_JNUM
    INTEGER(C_INTPTR_T) :: alloc_local, local_JNUM, local_joff
#endif
    !------------------------------------------------------------------------!
    CALL this%InitGravity(Mesh,Physics,"shearingbox spectral solver",config,IO)

    !-------------- checks & warnings for initial conditions ----------------!
#if !defined(HAVE_FFTW)
    CALL this%Error("InitGravity_sboxspectral", &
         "No fftw package could be loaded.")
#else
    CALL GetAttr(config, "order", this%order, 2)
    SELECT CASE(this%order)
    CASE(2,4)
      ! ok, do nothing
    CASE DEFAULT
      CALL this%Error("InitGravity_sboxspectral", &
        "Order must be one of {2,4}.")
    END SELECT

#ifdef PARALLEL
    CALL fftw_mpi_init()
    C_INUM = Mesh%INUM
    C_JNUM = Mesh%JNUM
#endif

    ! rotational symmetry is not allowed
    IF (Mesh%ROTSYM.NE.0) &
      CALL this%Error("InitGravity_sboxspectral", &
         "Rotational symmetry not supported.")

    ! 1D simulations are not allowed
    IF (Mesh%NDIMS.EQ.1) &
      CALL this%Error("InitGravity_sboxspectral", &
         "Only 2D (shearingsheet) and 3D (shearingbox) simulations allowed with this module.")

    ! check physics
    SELECT TYPE (phys => Physics)
    CLASS IS(physics_eulerisotherm)
      ! do nothing
    CLASS DEFAULT
      CALL this%Error("InitGravity_sboxspectral", &
         "Physics modules with density necessary for this module.")
    END SELECT

    !------------------------------------------------------------------------!

#ifdef PARALLEL
    ! check the dimensions if fftw should be used in parallel in order to
    ! know how to allocate the local arrays
    IF (this%GetNumProcs().GT.1 .AND. Mesh%shear_dir.EQ.2) THEN
      CALL this%Error("InitGravity_sboxspectral", &
                 "Execution of the parallel Fourier solver is only possible with pencil "// &
                 "decomposition. The first dimension should be fully accessible by the "//&
                 "first core. The shear thereto needs to vary along the second " // &
                 "dimension. In order to achieve this set the boundaries accordingly. "// &
                 "Check InitBoundary where the shear direction is set.")
    END IF
    IF (Mesh%dims(1).GT.1 .AND. Mesh%shear_dir.EQ.1) THEN
      CALL this%Error("InitGravity_sboxspectral", &
                 "The first dimension needs to be fully accessible by each process. "// &
                 "Please change domain decomposition to pencil decompositon. " // &
                 "This needs to be done during initialization. In the Mesh dictionary " // &
                 "use the key 'decompositon' and the value (/ 1,-1, 1/) to force decomposition " // &
                 "along second direction")
    END IF
#endif

    ! check for vertical extent (3D)
    IF (Mesh%KNUM.GT.1) THEN
      IF (Mesh%zmin.NE.-Mesh%zmax) &
        CALL this%Error("InitGravity_sboxspectral","z_min != -zmax not allowed: " // &
                        ACHAR(13) // " computational domain must be symmetric " // &
                        "with respect to z=0 plane")
#ifdef PARALLEL
      CALL this%Error("InitGravity_sboxspectral","sboxspectral does currently not support parallelization in 3D")
#endif
    END IF

    ALLOCATE( &
             this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
#if !defined(PARALLEL)
             this%Fmass3D(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
             this%qk(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,0:Mesh%KNUM-1),&
             this%density(Mesh%KMIN:Mesh%KMAX), &
             this%FTdensity(Mesh%KMIN:Mesh%KMAX), &
#else
             ! initialization handled by FFTW (see below)
#endif
!!!! TODO: only allocate these for KNUM > 1
             this%Fsum3D(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
!!!! END TODO
             this%Kxy2(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX), &
             this%kx(Mesh%INUM), &
             this%ky(Mesh%JNUM), &
             this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
             STAT=err)
    IF (err.NE.0) &
        CALL this%Error("InitGravity_sboxspectral","Memory allocation failed.")

    this%local_joff = 0

    ! use special allocation pattern from fftw when using MPI in order to
    ! assure good alignment
#if defined(PARALLEL)
    this%alloc_local     = fftw_mpi_local_size_2d(C_JNUM, C_INUM, &
                              MPI_COMM_WORLD, this%local_JNUM, this%local_joff)
    this%mass2D_pointer  = fftw_alloc_real(2*this%alloc_local)
    this%Fmass2D_pointer = fftw_alloc_complex(this%alloc_local)
    call c_f_pointer(this%mass2D_pointer,this%mass2D, &
                     [2*(C_INUM/2+1),this%local_JNUM])
    call c_f_pointer(this%Fmass2D_pointer,this%Fmass2D, &
                     [C_INUM/2+1,this%local_JNUM])
#endif

    !------------------- create plans for fftw ------------------------------!
    ! Pay attention to the argument order of the dimension (JNUM and INUM    !
    ! are switched because of C -> row-major, Fortran -> column-major),      !
    ! BUT ONLY in modern Fortran UNLIKE the legacy version                   !
    ! ------------ plans are allocated in dictionary ------------------------!
    CALL this%Info("            initializing FFTW: " // &
#if  !defined(PARALLEL)
                   "serial mode")
    DO k=Mesh%KMIN,Mesh%KMAX
      this%density(k)%ptr2D => this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)
      this%FTdensity(k)%ptr2D => this%Fmass3D(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,k)
    END DO
!!! REMOVE this, if this%mass2D and this%Fmass2D has been replaced everywhere
    this%mass2D(1:,1:) => this%density(Mesh%KMIN)%ptr2D
    this%Fmass2D(1:,1:) => this%FTdensity(Mesh%KMIN)%ptr2D
    this%plan_r2c = fftw_plan_dft_r2c_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1,this%density(Mesh%KMIN)%ptr2D, &
                                         this%FTdensity(Mesh%KMIN)%ptr2D,FFTW_MEASURE)
    this%plan_c2r = fftw_plan_dft_c2r_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1, this%FTdensity(Mesh%KMIN)%ptr2D, &
                                         this%density(Mesh%KMIN)%ptr2D, FFTW_MEASURE)
    IF ((.NOT. C_ASSOCIATED(this%plan_r2c)) .OR. (.NOT. c_associated(this%plan_c2r))) THEN
       CALL this%Error("InitGravity_sboxspectral","FFT plan could not be created.")
    END IF
#elif defined(PARALLEL)
                   "parallel mode")
    this%plan_r2c = fftw_mpi_plan_dft_r2c_2d(C_JNUM,C_INUM, &
                                             this%mass2D,this%Fmass2D, &
                                             MPI_COMM_WORLD, FFTW_MEASURE)
    this%plan_c2r = fftw_mpi_plan_dft_c2r_2d(C_JNUM,C_INUM, &
                                             this%Fmass2D,this%mass2D, &
                                             MPI_COMM_WORLD, FFTW_MEASURE)
#endif
    CALL GetAttr(config, "print_plans", valwrite, 0)
    IF (valwrite.GT.0) THEN
      IF (this%GetRank().EQ.0) THEN
        CALL fftw_print_plan(this%plan_r2c)
        PRINT *,ACHAR(13)
        CALL fftw_print_plan(this%plan_c2r)
        PRINT *,ACHAR(13)
      END IF
    END IF

    !------------------------------------------------------------------------!
    ! initialize arrays
    this%phi(:,:,:) = 0.
    this%Fmass3D(:,:,:) = CMPLX(0.,0)
    this%Fsum3D(:,:,:) = CMPLX(0.,0)
    this%kx(:) = 0.
    this%ky(:) = 0.
    this%qk(:,:,:) = 0.
    ! nullify not used potential explicitely
    NULLIFY(this%pot)

    ! precompute wave numbers
    this%kx(:) = CSHIFT((/(i-(Mesh%INUM+1)/2,i=0,Mesh%INUM-1)/),(Mesh%INUM+1)/2) &
                   *2.*PI/(Mesh%XMAX-Mesh%XMIN)
    this%ky(:) = CSHIFT((/(j-(Mesh%JNUM+1)/2,j=0,Mesh%JNUM-1)/),(Mesh%JNUM+1)/2) &
                   *2.*PI/(Mesh%YMAX-Mesh%YMIN)

    ! cut-off value for wave numbers
    this%maxKxy2 = 0.5*MIN(PI/Mesh%dx,PI/Mesh%dy)**2

    ! precompute shift constant
    SELECT CASE (Mesh%shear_dir)
    CASE(1)
      this%shiftconst = Mesh%Q*Mesh%OMEGA*(Mesh%ymax-Mesh%ymin)/Mesh%dx
      ALLOCATE(&
             this%joff(Mesh%JMIN:Mesh%JMAX), &
             this%jrem(Mesh%JMIN:Mesh%JMAX), &
             STAT=err)
    CASE(2)
      this%shiftconst = Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)/Mesh%dy
      ALLOCATE(&
             this%joff(Mesh%IMIN:Mesh%IMAX), &
             this%jrem(Mesh%IMIN:Mesh%IMAX), &
             STAT=err)
    CASE DEFAULT
      CALL this%Error("InitGravity_sboxspectral", &
        "Shear direction must be one of WE_shear (x-direction) or SN_shear (y-direction).")
    END SELECT
    IF (err.NE.0) &
      CALL this%Error("InitGravity_sboxspectral","Memory allocation failed for joff & jrem.")
    this%joff(:) = 0.
    this%jrem(:) = 0.

#endif
  END SUBROUTINE InitGravity_sboxspectral

  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    TYPE(Dict_TYP),      POINTER       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER          :: valwrite,err
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    valwrite = 0
    CALL GetAttr(config, "output/potential", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
      IF (ASSOCIATED(this%phi)) &
        CALL SetAttr(IO, "potential", &
              this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    END IF
    valwrite = 0
    CALL GetAttr(config, "output/Fmass3D", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
      ALLOCATE(this%Fmass3D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
               STAT=err)
      IF (err.NE.0) &
        CALL this%Error("sboxspectral::SetOutput","Memory allocation failed for Fmass3D_real.")
      this%Fmass3D_real(:,:,:) = 0.0
      CALL SetAttr(IO, "Fmass3D_real", &
              this%Fmass3D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    END IF
#endif
  END SUBROUTINE SetOutput

  !> \public Calculates the acceleration from potential.
  !!
  !! Calculates
  !! \f[
  !!     \mathbf{a} = -\nabla \Phi.
  !! \f]
  !! Uses second order symmetric difference quotient.
  SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),            INTENT(IN) :: Mesh
    CLASS(physics_base),         INTENT(IN) :: Physics
    CLASS(fluxes_base),          INTENT(IN) :: Fluxes
    CLASS(marray_compound),   INTENT(INOUT) :: pvar
    REAL,                        INTENT(IN) :: time,dt
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    REAL    :: w1,w2
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! calc potential first
    CALL this%CalcPotential(Mesh,Physics,time,pvar)

    IF (this%order.EQ.2) THEN
!NEC$ ivdep
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ ivdep
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ ivdep
          DO i = Mesh%IMIN,Mesh%IMAX
            ! second order approximation
            this%accel%data4d(i,j,k,1) = -1.0*(this%phi(i+1,j,k)-this%phi(i-1,j,k))/ &
                                 (2*Mesh%dlx%data3d(i,j,k))
            this%accel%data4d(i,j,k,2) = -1.0*(this%phi(i,j+1,k)-this%phi(i,j-1,k))/ &
                                 (2*Mesh%dly%data3d(i,j,k))
         END DO
        END DO
      END DO
    ELSE ! this%order.EQ.4
      w1 = 3./48.
      w2 = 30./48.
!NEC$ ivdep
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ ivdep
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ ivdep
          DO i = Mesh%IMIN,Mesh%IMAX
            ! fourth order
            this%accel%data4d(i,j,k,1) = -1.0*(w1*this%phi(i-2,j,k)-w2*this%phi(i-1,j,k)+ &
                                        w2*this%phi(i+1,j,k)-w1*this%phi(i+2,j,k))/ &
                                       (Mesh%dlx%data3d(i,j,k))
            this%accel%data4d(i,j,k,2) = -1.0*(w1*this%phi(i,j-2,k)-w2*this%phi(i,j-1,k)+ &
                                        w2*this%phi(i,j+1,k)-w1*this%phi(i,j+2,k)) / &
                                       (Mesh%dly%data3d(i,j,k))
          END DO
        END DO
      END DO
    END IF
#endif
  END SUBROUTINE UpdateGravity_single

  !> \public Computes the potential with FFT method within a shearingsheet.
  !!
  !! Calculates Poisson's equation:
  !! \f[
  !!    \Delta \Phi = 4 \pi G \Sigma \delta (z),
  !! \f]
  !! where \f$ \delta (z) \f$ is the \f$ \delta \f$-distribution, thus for
  !! razor-thin disks.
  !!
  !! The implemenation is done by the following steps:
  !! 1. Shift surface density field to next periodic point
  !!  \f[
  !!      \Sigma(x,y) \longrightarrow \Sigma(x,y').
  !!  \f]
  !!    See \link gravity_sboxspectral.fieldshift \endlink for more informations.
  !! 2. FFT of shifted surface density field
  !!  \f[
  !!      \Sigma(x,y') \longrightarrow \Sigma(\mathbf{k})
  !!  \f]
  !! 3. Calculate
  !!  \f[
  !!      \Phi(\mathbf{k}) = -\frac{2\pi G}{|\mathbf{k}|}\Sigma(\mathbf{k})
  !!                e^{i\mathbf{k\cdot x - |kz|}}.
  !!  \f]
  !! 4. FFT\f$^{-1}\f$ of Potential
  !!  \f[
  !!      \Phi(\mathbf{k}) \longrightarrow \Phi(x,y')
  !!  \f]
  !! 5. Backshift of Potential
  !!  \f[
  !!      \Phi(x,y') \longrightarrow \Phi(x,y).
  !!  \f]
  !!
  !! The acceleration is eventually calculated in
  !! \link gravity_sboxspectral.updategravity_single \endlink.
#ifdef HAVE_FFTW
  SUBROUTINE CalcPotential(this,Mesh,Physics,time,pvar)
    USE physics_eulerisotherm_mod, ONLY : statevector_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    REAL,                INTENT(IN) :: time
    CLASS(marray_compound),   INTENT(INOUT) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k,kk
    REAL    :: joff2,jrem2,delt,time0
#ifdef PARALLEL
    INTEGER :: status(MPI_STATUS_SIZE),ierror
    REAL    :: mpi_buf(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) !TODO: ONLY 2D
#endif
    !------------------------------------------------------------------------!
    !---------------- fourier transformation of density ---------------------!
    ! calculate the shift of the indice at time t                            !
    ! shift density to pretend periodic behavior with interpolation          !

    IF (Mesh%OMEGA .EQ. 0) THEN
      time0 = 0.0
    ELSE IF (Mesh%shear_dir.EQ.2) THEN
      time0 = NINT(time*Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)/ &
          (Mesh%ymax-Mesh%ymin))*(Mesh%ymax-Mesh%ymin)/(Mesh%Q*Mesh%OMEGA &
          *(Mesh%xmax-Mesh%xmin))
    ELSE IF (Mesh%shear_dir.EQ.1) THEN
      time0 = NINT(time*Mesh%Q*Mesh%OMEGA*(Mesh%ymax-Mesh%ymin)/ &
          (Mesh%xmax-Mesh%xmin))*(Mesh%xmax-Mesh%xmin)/(Mesh%Q*Mesh%OMEGA &
          *(Mesh%ymax-Mesh%ymin))
    END IF
    delt = time - time0

    !----------------- shift field to periodic point ------------------------!
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_eulerisotherm)
      CALL this%FieldShift(Mesh,Physics,delt, &
                p%density%data3d(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    CLASS DEFAULT
      CALL this%Error("gravity_sboxspectral::CalcPotential","unsupported state vector")
    END SELECT


#if defined(PARALLEL)
    this%mass2D(1:Mesh%INUM,1:this%local_JNUM) = &
      RESHAPE(this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
              (/ Mesh%INUM,Mesh%JNUM/this%GetNumProcs() /))
#endif

    !-------------- fourier transform (forward) of shifted density ----------!
! for performance checks
#ifdef _FTRACE
CALL ftrace_region_begin("forward_fft")
#endif

    !> \todo Violation of best practices. It should take Fmass2D and mass2D as
    !!      dummy argument, but this is different for every combination
!!NEC$ IEXPAND
    CALL this%FFT_Forward(Mesh,Physics)

! performance checks
#ifdef _FTRACE
CALL ftrace_region_end("forward_fft")
#endif

    !------------- calculate wave number vector squared --------------------!
    IF (Mesh%shear_dir.EQ.2) THEN
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          this%Kxy2(i,j) = (this%kx(i) + Mesh%Q*Mesh%OMEGA*this%ky(j)*delt)**2 + this%ky(j)**2
        END DO
      END DO
    ELSE IF (Mesh%shear_dir.EQ.1) THEN
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          this%Kxy2(i,j) = this%kx(i)**2 + (this%ky(j)-Mesh%Q*mesh%OMEGA*this%kx(i)*delt)**2
        END DO
      END DO
    END IF

    !------------- calculate potential in fourier domain --------------------!
    IF (ABS(Mesh%zmax-Mesh%zmin).LT.TINY(Mesh%zmin)) THEN
      ! 2D razor thin disk
      k=Mesh%KMIN
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          IF ((this%Kxy2(i,j).LT.this%maxKxy2).AND.(this%Kxy2(i,j).GT. 0)) THEN
              this%Fmass3D(i,j-this%local_joff,k) = -2.*PI*Physics%Constants%GN &
                *this%Fmass3D(i,j-this%local_joff,k)/SQRT(this%Kxy2(i,j))
          ELSE
              this%Fmass3D(i,j-this%local_joff,k) = 0.0
          END IF
        END DO
      END DO
!       k=Mesh%KMIN
!       WHERE ((this%Kxy2(:,:).LT.this%maxKxy2).AND.(this%Kxy2(:,:).GT. 0))
!         this%Fmass3D(:,:,k) = -2.*PI*Physics%Constants%GN &
!           *this%Fmass3D(:,:,k)/SQRT(this%Kxy2(:,:))
!       ELSEWHERE
!         this%Fmass3D(:,:,k) = 0.0
!       END WHERE
    ELSE
      ! disk with vertical extent
      this%qk(:,:,0) = EXP(-0.5*Mesh%dz*SQRT(this%Kxy2(:,:)))
      ! compute weight factors for disks with more than one vertical cell
      IF (Mesh%KNUM.EQ.1) THEN
        ! 2D disk with one zone with finite vertical size
        k=Mesh%KMIN
!NEC$ IVDEP
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX/2+1
            IF ((this%Kxy2(i,j).LT.this%maxKxy2).AND.(this%Kxy2(i,j).GT. 0)) THEN
              this%Fmass3D(i,j-this%local_joff,k) = -4.*PI*Physics%Constants%GN &
                *(1.-this%qk(i,j,0))/this%Kxy2(i,j) * this%Fmass3D(i,j-this%local_joff,k)
            ELSE
              this%Fmass3D(i,j-this%local_joff,k) = 0.0
            END IF
          END DO
        END DO
      ELSE
        ! 3D disk with more than one vertical zone
        this%qk(:,:,1) = this%qk(:,:,0)*this%qk(:,:,0)
        DO k=2,Mesh%KNUM-1
          this%qk(:,:,k) = this%qk(:,:,k-1) * this%qk(:,:,1)
        END DO
        DO k=Mesh%KMIN,Mesh%KMAX
          this%Fsum3D(:,:,k) = 0.0
!!!! the following seems to be equivalent to a matrix vector multiplication
!!!! where the matrix is a skew symmetric band matrix
          DO kk=Mesh%KMIN,k-1
            this%Fsum3D(:,:,k) = this%Fsum3D(:,:,k) + this%qk(:,:,ABS(k-kk)) &
                                    * SIGN(1,k-kk) * this%Fmass3D(:,:,kk)
          END DO
          ! skip k=kk
          DO kk=k+1,Mesh%KMAX
            this%Fsum3D(:,:,k) = this%Fsum3D(:,:,k) + this%qk(:,:,ABS(k-kk)) &
                                    * SIGN(1,k-kk) * this%Fmass3D(:,:,kk)
          END DO
!NEC$ IVDEP
          DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i = Mesh%IMIN,Mesh%IMAX/2+1
              IF ((this%Kxy2(i,j).LT.this%maxKxy2).AND.(this%Kxy2(i,j).GT. 0)) THEN
                this%Fmass3D(i,j,k) = -4.*PI*Physics%Constants%GN &
                  * (1.-this%qk(i,j,0))/this%Kxy2(i,j) &
                  * (this%Fmass3D(i,j,k) + 0.5*(1.0+1./this%qk(i,j,0)) &
                    * this%Fsum3D(i,j,k))
              ELSE
                this%Fmass3D(i,j-this%local_joff,k) = 0.0
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF

    !----------- fourier transform (backward) of shifted density -----------!
#ifdef _FTRACE
CALL ftrace_region_begin("backward_fft")
#endif

    !> \todo Violation of best practices. It should take Fmass3D and den_ip as
    !!      dummy argument, but this is different for every combination.
    !! The potential is now in this%den_ip.
    CALL this%FFT_Backward(Mesh,Physics)

#ifdef _FTRACE
CALL ftrace_region_end("backward_fft")
#endif

    !------ calculate final potential with backshift and normalization ------!
!NEC$ IVDEP
    DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX
          this%phi(i,j,k) = this%mass2D(i,j-this%local_joff)/ &
               (Mesh%JNUM*Mesh%INUM)                        ! no norm. by FFTW !
        END DO
      END DO
    END DO
    !-------------------------- shift field back ----------------------------!
!!NEC$ IEXPAND
    CALL this%FieldShift(Mesh,Physics,-delt, &
            this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
            this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

!NEC$ IVDEP
    DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX
          this%phi(i,j,k) = this%den_ip(i,j,k)
        END DO
      END DO
    END DO

    !----- copy values to boundaries in order to calculate acceleration -----!
    IF(Mesh%shear_dir.EQ.2) THEN
!NEC$ SHORTLOOP
      DO j = 1,Mesh%GJNUM
        ! southern northern (periodic)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX)
      END DO
      joff2 = -this%shiftconst*delt
      jrem2 = joff2 - FLOOR(joff2)
!NEC$ SHORTLOOP
      DO i = 1,Mesh%GINUM
        ! western (periodic)
        this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
        ! western integer shift
        this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)  = &
              CSHIFT(this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
        ! western residual shift
        this%phi(Mesh%IMIN-i,Mesh%JMAX+1,:) = this%phi(Mesh%IMIN-i,Mesh%JMIN,:)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            this%phi(Mesh%IMIN-i,j,k) = &
                (1.0 - jrem2)*this%phi(Mesh%IMIN-i,j,k) + jrem2*this%phi(Mesh%IMIN-i,j+1,k)
          END DO
        END DO
      END DO
      joff2 = this%shiftconst*delt
      jrem2 = joff2 - FLOOR(joff2)
!NEC$ SHORTLOOP
      DO i = 1,Mesh%GINUM
        ! eastern (periodic)
        this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
        ! eastern integer shift
        this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)  = &
              CSHIFT(this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
        ! eastern residual shift
        this%phi(Mesh%IMAX+i,Mesh%JMAX+1,:) = this%phi(Mesh%IMAX+i,Mesh%JMIN,:)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            this%phi(Mesh%IMAX+i,j,k) = &
                (1.0 - jrem2)*this%phi(Mesh%IMAX+i,j,k) + jrem2*this%phi(Mesh%IMAX+i,j+1,k)
          END DO
        END DO
      END DO
! Only north-south direction has parallelization allowed
! Attention: The order of the copies plays a role. First the non-shifted direction needs to be
! copied, afterwards the shifted.
    ELSE IF(Mesh%shear_dir.EQ.1) THEN
!NEC$ SHORTLOOP
      DO i = 1,Mesh%GINUM
        ! western and eastern (always periodic)
        this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMAX-i+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
        this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN+i-1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
      END DO
#ifdef PARALLEL
      IF(Mesh%dims(2).GT.1) THEN
        mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) = &
          this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-Mesh%GJNUM+1:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
        CALL MPI_Sendrecv_replace(&
          mpi_buf,&
          2*(Mesh%IMAX-Mesh%IMIN+1), &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(NORTH), 53+NORTH, &
          Mesh%neighbor(SOUTH), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,Mesh%KMIN:Mesh%KMAX) = &
          mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX)

        mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) = &
          this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMIN+Mesh%GJNUM-1,Mesh%KMIN:Mesh%KMAX)
        CALL MPI_Sendrecv_replace(&
          mpi_buf,&
          2*(Mesh%IMAX-Mesh%IMIN+1), &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(SOUTH), 53+SOUTH, &
          Mesh%neighbor(NORTH), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX) = &
          mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX)
      ELSE
!NEC$ SHORTLOOP
        DO j = 1,Mesh%GJNUM
          ! southern northern (periodic in first step - further shift-treatment below)
          this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX)
          this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX)
        END DO
      END IF
#endif

#ifdef PARALLEL
      IF (Mesh%mycoords(2).EQ.0) THEN
#endif
!NEC$ SHORTLOOP
      DO j = 1,Mesh%GJNUM
        ! southern (shorn periodic) - residual and integer shift
#ifndef PARALLEL
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX)
#endif
        joff2 = this%shiftconst*delt
        jrem2 = joff2 - FLOOR(joff2)
        !------- integral shift ---------------------------------------------!
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX)  = &
          CSHIFT(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
        !------- residual shift ---------------------------------------------!
        this%phi(Mesh%IMAX+1,Mesh%JMIN-j,:) = this%phi(Mesh%IMIN,Mesh%JMIN-j,:)
        DO k = Mesh%KMIN,Mesh%KMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMIN-j,k) = &
              (1.0 - jrem2)*this%phi(i,Mesh%JMIN-j,k) + jrem2*this%phi(i+1,Mesh%JMIN-j,k)
          END DO
        END DO
      END DO
#ifdef PARALLEL
      END IF
#endif

#ifdef PARALLEL
      IF (Mesh%mycoords(2).EQ.Mesh%dims(2)-1) THEN
#endif
!NEC$ SHORTLOOP
      DO j = 1,Mesh%GJNUM
#ifndef PARALLEL
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX)
#endif
        ! northern (shorn periodic) - residual and integer shift
        joff2 = -this%shiftconst*delt
        jrem2 = joff2 - FLOOR(joff2)
        !------- integral shift ---------------------------------------------!
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX)  = &
          CSHIFT(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
        !------- residual shift ---------------------------------------------!
        this%phi(Mesh%IMAX+1,Mesh%JMAX+j,:) = this%phi(Mesh%IMIN,Mesh%JMAX+j,:)
        DO k = Mesh%KMIN,Mesh%KMAX
          DO i = Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMAX+j,k) = &
              (1.0 - jrem2)*this%phi(i,Mesh%JMAX+j,k) + jrem2*this%phi(i+1,Mesh%JMAX+j,k)
          END DO
        END DO
      END DO
#ifdef PARALLEL
      END IF
#endif
    END IF
  END SUBROUTINE CalcPotential


  !> \private Calculates the FFT forward
  SUBROUTINE FFT_Forward(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT):: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    !------------------------------------------------------------------------!

#ifdef _FTRACE
CALL ftrace_region_begin("forward FFT")
#endif

#if defined(PARALLEL)
    ! 3D is currently not supported in parallel
    k=Mesh%KMIN
    CALL fftw_mpi_execute_dft_r2c(this%plan_r2c,this%density(k)%ptr2D, &
                                  this%FTdensity(k)%ptr2D)
#else
    DO k=Mesh%KMIN,Mesh%KMAX
      CALL fftw_execute_dft_r2c(this%plan_r2c,this%density(k)%ptr2D, &
                                this%FTdensity(k)%ptr2D)
    END DO
#endif

    IF (ASSOCIATED(this%Fmass3D_real)) THEN
!NEC$ IVDEP
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX/2+1
            this%Fmass3D_real(2*i-1,j,k) = REAL(REAL(this%Fmass3D(i,j-this%local_joff,k)))
            this%Fmass3D_real(2*i,j,k)   = REAL(AIMAG(this%Fmass3D(i,j-this%local_joff,k)))
          END DO
        END DO
      END DO
    END IF

#ifdef _FTRACE
CALL ftrace_region_end("foward FFT")
#endif

  END SUBROUTINE


  !> \private Calculates the FFT backward transformation
  SUBROUTINE FFT_Backward(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    !------------------------------------------------------------------------!
    INTEGER           :: k
    !------------------------------------------------------------------------!
#if defined(PARALLEL)
    ! 3D is currently not supported in parallel
    k=Mesh%KMIN
    CALL fftw_mpi_execute_dft_c2r(this%plan_c2r,this%FTdensity(k)%ptr2D, &
                                  this%density(k)%ptr2D)
#else
    DO k=Mesh%KMIN,Mesh%KMAX
      CALL fftw_execute_dft_c2r(this%plan_c2r,this%FTdensity(k)%ptr2D, &
                                this%density(k)%ptr2D)
    END DO
#endif
  END SUBROUTINE
#endif


  !> \public Compute disk pressure scale height for geometrically thin
  !! self-gravitating shearingsheet.
  !!
  !! The algorithm solves the equations
  !! \f{eqnarray}{
  !!     1/h^2 &=& 1/h_{sg}^2 + 1/h_{ext}^2 & (1)\\
  !!     \left.\frac{\partial^2\Phi_{sg}}{\partial z^2}\right|_{z=0} =
  !!           (c_s/h_{sg})^2 &=& 4\pi G \rho_c & (2) \\
  !!     \rho_c &=& 1/\sqrt{2\pi} \Sigma / h & (3) \\
  !!     1/h_{ext}^2 &=& \Omega_{K}^2 / c_{s}^2
  !! \f}
  !! for the disk scale height \f$ h \f$ where \f$ h_{sg} \f$ is the self-
  !! gravitating pressure scale height due to the material in the disk and
  !! \f$ h_{ext} \f$ the scale height caused by an external (point mass)
  !! potential. \f$ \Sigma \f$ is the surface density, \f$ c_s \f$ is the
  !! speed of sound and \f$ \rho_c \f$ is the density in the equatorial plane.
  !! \f$ \Phi_{sg} \f$ is the gravitational potential in the shearing sheet.
  !! Replacing \f$ \rho_c \f$ in (2) by (3) and inserting the result in (1)
  !! together with (4) yields a quadratic equation  with the two solutions
  !! \f[
  !!     h = -p \pm \sqrt{q+p^2}
  !!       = p (-1 \pm \sqrt{1 + q/p^2})
  !! \f]
  !! where \f$ p = \sqrt{2\pi} G \Sigma / \Omega_{K}^2 \f$ and
  !! \f$ q = c_{s}^2 / \Omega_{K}^2 \f$.
  !! The sign of the root has to be positive, because in the non self-
  !! gravitating limit \f$ p=0, q=1 \f$ and \f$ h \f$ needs to be
  !! \f$ c_{s}/\Omega_{K} \f$ not \f$ -c_{s}/\Omega_{K} \f$.
  SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),            INTENT(IN)    :: Mesh
    CLASS(physics_base),         INTENT(IN)    :: Physics
    CLASS(marray_compound),      INTENT(INOUT) :: pvar
    TYPE(marray_base),           INTENT(INOUT) :: bccsound,h_ext,height
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    REAL              :: cs2,p,q
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! pure self-gravitating shearing sheet with external point mass potential
!NEC$ COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ IVDEP
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
           cs2 = bccsound%data3d(i,j,k)*bccsound%data3d(i,j,k)
           p = -SQRTTWOPI*Physics%Constants%GN*pvar%data4d(i,j,k,Physics%DENSITY) &
                  /Mesh%OMEGA**2.
           q = cs2/Mesh%OMEGA**2.
           ! return the new disk height
           height%data3d(i,j,k) = p+SQRT(q+p*p)
        END DO
      END DO
    END DO
#endif
  END SUBROUTINE CalcDiskHeight_single

  !> \public Shifts the whole field to the next periodic point.
  !!
  !! Implementation is done by
  !! \f[
  !!     y' = y - q \Omega x (t-t_p),
  !! \f]
  !! with \f$ t_p = \text{NINT}(q\Omega t) / (q\Omega) \f$. In order to map
  !! the continuous shift at the discrete field linear interpolation is used,
  !! and assumes periodic behaviour along the y-direction.
  SUBROUTINE FieldShift(this,Mesh,Physics,delt,field,shifted_field)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)  :: Mesh
    CLASS(physics_base), INTENT(IN)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)  :: field
#if defined(PARALLEL)
    REAL, DIMENSION(1:Mesh%INUM,1:this%local_JNUM,Mesh%KMIN:Mesh%KMAX), &
                         INTENT(OUT) :: shifted_field
#else
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
                         INTENT(OUT) :: shifted_field
#endif
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,local_joff
    REAL              :: delt
    !------------------------------------------------------------------------!
    local_joff = this%local_joff

    IF (Mesh%shear_dir.EQ.1) THEN
!NEC$ IVDEP
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j = Mesh%JMIN,Mesh%JMAX
          this%joff(j)   = Mesh%Q*Mesh%OMEGA*Mesh%bcenter(Mesh%IMIN,j,k,2)* &
                      delt/Mesh%dx
          this%jrem(j)  = this%joff(j) - FLOOR(this%joff(j))
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX
            shifted_field(i,j-local_joff,k) = &
             (1.0-this%jrem(j))*field(1+MODULO(i-1+FLOOR(this%joff(j)), &
              Mesh%IMAX-Mesh%IMIN+1),j,k) + this%jrem(j)*field(1+MODULO(i+ &
              FLOOR(this%joff(j)),Mesh%IMAX-Mesh%IMIN+1),j,k)
          END DO
        END DO
      END DO
    ELSE ! must be Mesh%shear_dir.EQ.2, because otherwise initialization would raise an error
!NEC$ IVDEP
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX
          this%joff(i)   = -Mesh%Q*Mesh%OMEGA*Mesh%bcenter(i,Mesh%JMIN,k,1)* &
                      delt/Mesh%dy
          this%jrem(i)  = this%joff(i) - FLOOR(this%joff(i))
        END DO
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX
            shifted_field(i,j,k) = &
             (1.0-this%jrem(i))*field(i,1+MODULO(j-1+FLOOR(this%joff(i)), &
              Mesh%JMAX-Mesh%JMIN+1),k) + this%jrem(i)*field(i,1+MODULO(j+ &
              FLOOR(this%joff(i)),Mesh%JMAX-Mesh%JMIN+1),k)
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE FieldShift

  !> \public Closes the gravity term of the shearingsheet spectral solver.
  SUBROUTINE Finalize(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
   !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    CALL fftw_destroy_plan(this%plan_r2c)
    CALL fftw_destroy_plan(this%plan_c2r)
#if defined(PARALLEL)
    CALL fftw_free(this%mass2D_pointer)
    CALL fftw_free(this%Fmass2D_pointer)
#else
    DEALLOCATE(this%Fmass3D,this%Fsum3D,this%density,this%FTdensity)
#endif
    ! Free memory
    IF (ASSOCIATED(this%Fmass3D_real)) DEALLOCATE(this%Fmass3D_real)
    DEALLOCATE(&
               this%phi, &
               this%kx, this%ky, this%Kxy2, &
               this%qk, &
               this%joff, &
               this%jrem, &
               this%den_ip &
               )
#endif
    END SUBROUTINE Finalize

END MODULE gravity_sboxspectral_mod
