!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_sboxspectral.f90                                          #
!#                                                                           #
!# Copyright (C) 2015-2019                                                   #
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
!> \addtogroup gravity
!! - parameters of \link gravity_sboxspectral_mod gravity_sboxspectral \endlink as key-values
!!  \key{output/phi,INTEGER,(enable=1) output of gravitational potential}
!!  \key{output/accel_x,INTEGER,(enable=1) output of accel in x-direction}
!!  \key{output/accel_y,INTEGER,(enable=1) output of accel in y-direction}
!!  \key{output/mass2D,INTEGER,(enable=1) output of density that will be
!!                              fourier transformed}
!!  \key{output/Fmass2D,INTEGER,(enable=1) output of density in spectral
!!                              space}
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!!
!! \brief Poisson solver via spectral methods within the shearingsheet.
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

  TYPE, EXTENDS(gravity_spectral) :: gravity_sboxspectral
#ifdef HAVE_FFTW
    !> \name
    !!#### spectral poisson solver shearing box
    REAL(C_DOUBLE), POINTER          :: mass2D(:,:)    !< temporary variable
    COMPLEX(C_DOUBLE_COMPLEX), POINTER &
                                     :: Fmass2D(:,:)   !< temporary variable
    REAL, DIMENSION(:,:,:), POINTER  :: Fmass2D_real   !< temporary variable
    INTEGER(C_INTPTR_T)              :: local_joff
    REAL,DIMENSION(:), POINTER       :: kx             !< wave numbers for FFT (x)
    REAL,DIMENSION(:), POINTER       :: ky             !< wave numbers for FFT (y)
    REAL                             :: shiftconst     !< constant for shift
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
    PROCEDURE :: InfoGravity
    PROCEDURE :: SetOutput
    PROCEDURE :: CalcDiskHeight_single
    PROCEDURE :: Finalize
#ifdef HAVE_FFTW
    PROCEDURE :: CalcPotential
    PROCEDURE :: FFT_Forward
    PROCEDURE :: FFT_Backward
    PROCEDURE :: FieldShift
#endif
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
    INTEGER             :: err, valwrite, i
#if defined(HAVE_FFTW) && defined(PARALLEL)
    INTEGER(C_INTPTR_T) :: C_INUM,C_JNUM
    INTEGER(C_INTPTR_T) :: alloc_local, local_JNUM, local_joff
    INTEGER             :: nprocs
#endif
    !------------------------------------------------------------------------!
    CALL this%InitGravity(Mesh,Physics,"shearingbox spectral solver",config,IO)

    !-------------- checks & warnings for initial conditions ----------------!
#if !defined(HAVE_FFTW)
    CALL this%Error("InitGravity_sboxspectral", &
         "No fftw package could be loaded.")
#else
    CALL GetAttr(config, "order", this%order, 2)

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

    ! check the dimensions if fftw should be used in parallel in order to
    ! know how to allocate the local arrays
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)
    IF (nprocs.GT.1 .AND. Mesh%shear_dir.EQ.2) THEN
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
    ALLOCATE( &
             this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
#if !defined(PARALLEL)
             this%mass2D(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
             this%Fmass2D(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX), &
#else
             ! initialization handled by FFTW (see below)
#endif
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
#if  !defined(PARALLEL)
    this%plan_r2c = fftw_plan_dft_r2c_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1,this%mass2D, &
                                         this%Fmass2D,FFTW_MEASURE)
    this%plan_c2r = fftw_plan_dft_c2r_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1, this%Fmass2D, &
                                         this%mass2D, FFTW_MEASURE)
    IF ((.NOT. C_ASSOCIATED(this%plan_r2c)) .OR. (.NOT. c_associated(this%plan_c2r))) THEN
       CALL this%Error("InitGravity_sboxspectral","FFT plan could not be created.")
    END IF
#elif defined(PARALLEL)
    this%plan_r2c = fftw_mpi_plan_dft_r2c_2d(C_JNUM,C_INUM, &
                                             this%mass2D,this%Fmass2D, &
                                             MPI_COMM_WORLD, FFTW_MEASURE)
    this%plan_c2r = fftw_mpi_plan_dft_c2r_2d(C_JNUM,C_INUM, &
                                             this%Fmass2D,this%mass2D, &
                                             MPI_COMM_WORLD, FFTW_MEASURE)
#endif
    !------------------------------------------------------------------------!
    ! set potential and acceleration to zero
    this%phi(:,:,:) = 0.
    this%mass2D(:,:) = 0.
    this%Fmass2D(:,:) = CMPLX(0.,0)
    this%kx(:) = 0.
    this%ky(:) = 0.
    ! nullify not used potential explicitely
    NULLIFY(this%pot)

    ! precompute wavenumbers
    this%kx = cshift((/(i-(Mesh%INUM+1)/2,i=0,Mesh%INUM-1)/), &
                +(Mesh%INUM+1)/2)*2.*PI/(Mesh%XMAX-Mesh%XMIN)
    this%ky = cshift((/(i-(Mesh%JNUM+1)/2,i=0,Mesh%JNUM-1)/), &
                +(Mesh%JNUM+1)/2)*2.*PI/(Mesh%YMAX-Mesh%YMIN)

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
      CALL this%Error("InitGravity_sboxspectral","Memory allocation failed.")
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
    valwrite = 0
    CALL GetAttr(config, "output/potential", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
      IF (ASSOCIATED(this%phi)) &
        CALL SetAttr(IO, "potential", &
              this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    END IF
#ifdef HAVE_FFTW
    valwrite = 0
    CALL GetAttr(config, "output/Fmass2D", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
      ALLOCATE(this%Fmass2D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
        STAT=err)
      IF (err.NE.0) &
         CALL this%Error("sboxspectral::SetOutput","Memory allocation failed.")
      this%Fmass2D_real(:,:,:) = 0.0
      CALL SetAttr(IO, "Fmass2D_real", &
              this%Fmass2D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
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
    ELSE IF (this%order.EQ.4) THEN
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
    ELSE
      CALL this%Error("InitGravity_sboxspectral", &
        "The chosen order approximation does not exist (only 2 and 3).")
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
    INTEGER :: i,j,k
    REAL    :: K2max,K2,KO
    REAL    :: joff2,jrem2,delt,time0
    INTEGER :: nprocs
#ifdef PARALLEL
    INTEGER :: status(MPI_STATUS_SIZE),ierror,ier
    REAL    :: mpi_buf(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) !TODO: ONLY 2D
#endif
    !------------------------------------------------------------------------!
    !---------------- fourier transformation of density ---------------------!
    ! calculate the shift of the indice at time t                            !
    ! shift density to pretend periodic behavior with interpolation          !
    nprocs=1
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ier)
#endif

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
      RESHAPE(this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), (/ Mesh%INUM,Mesh%JNUM/nprocs /))
#else
    this%mass2D(1:Mesh%INUM,1:Mesh%JNUM/nprocs) = &
      RESHAPE(this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), (/ Mesh%INUM,Mesh%JNUM/nprocs /))
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

    !------------- calculate potential in fourier domain --------------------!
    IF (Mesh%dx .GT. Mesh%dy) THEN
      K2max = 0.5*PI**2/(Mesh%dx**2)
    ELSE
      K2max = 0.5*PI**2/(Mesh%dy**2)
    END IF

    IF (Mesh%shear_dir.EQ.2) THEN
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          K2 = (this%kx(i) + Mesh%Q*Mesh%OMEGA*this%ky(j)*delt)**2 + this%ky(j)**2
          KO = SQRT(K2)
          IF ((K2 .LT. K2max) .AND. (K2 .GT. 0)) THEN
            this%Fmass2D(i,j-this%local_joff) = -2.*PI*Physics%Constants%GN* &
                                              this%Fmass2D(i,j-this%local_joff)/KO
          ELSE
            this%Fmass2D(i,j-this%local_joff) = 0.0
          END IF
        END DO
      END DO
    ELSE IF (Mesh%shear_dir.EQ.1) THEN
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          K2 = this%kx(i)**2 + (this%ky(j)-Mesh%Q*mesh%OMEGA*this%kx(i)*delt)**2
           KO = SQRT(K2)
          IF ((K2 .LT. K2max) .AND. (K2 .GT. 0)) THEN
            this%Fmass2D(i,j-this%local_joff) = -2.*PI*Physics%Constants%GN* &
                                              this%Fmass2D(i,j-this%local_joff)/KO
          ELSE
            this%Fmass2D(i,j-this%local_joff) = 0.0
          END IF
        END DO
      END DO
    END IF

    !----------- fourier transform (backward) of shifted density -----------!
#ifdef _FTRACE
CALL ftrace_region_begin("backward_fft")
#endif

    !> \todo Violation of best practices. It should take Fmass2D and mass2D as
    !!      dummy argument, but this is different for every combination.
!!NEC$ EXPAND
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
#ifdef PARALLEL
    INTEGER :: nprocs,status(MPI_STATUS_SIZE)
    INTEGER :: ier
#endif
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ier)
#endif

#ifdef _FTRACE
CALL ftrace_region_begin("forward FFT")
#endif

#if !defined(PARALLEL)
    CALL fftw_execute_dft_r2c(this%plan_r2c, this%mass2D, this%Fmass2D)
#else
    CALL fftw_mpi_execute_dft_r2c(this%plan_r2c,this%mass2D, this%Fmass2D)
#endif

    IF (ASSOCIATED(this%Fmass2D_real)) THEN
!NEC$ IVDEP
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX/2+1
            this%Fmass2D_real(2*i-1,j,k) = REAL(REAL(this%Fmass2D(i,j-this%local_joff)))
            this%Fmass2D_real(2*i,j,k)   = REAL(AIMAG(this%Fmass2D(i,j-this%local_joff)))
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
#ifdef PARALLEL
    INTEGER           :: nprocs,status(MPI_STATUS_SIZE)
    INTEGER           :: ier
#endif
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ier)
#endif

#if defined(PARALLEL)
    CALL fftw_mpi_execute_dft_c2r(this%plan_c2r,this%Fmass2D, this%mass2D)
#else
    CALL fftw_execute_dft_c2r(this%plan_c2r, this%Fmass2D, this%mass2D)
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

  !> Prints out information
  SUBROUTINE InfoGravity(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(IN) :: this
    CLASS(mesh_base),            INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    CALL this%Info( "            FFT-Package:       FFTW " // &
#if defined(PARALLEL)
                                                    "- parallel mode")
#else
                                                    "- serial mode")
#endif
    CALL this%Info("            .. done initializing")
#endif
  END SUBROUTINE InfoGravity

  !> \public Shifts the whole field to the next periodic point.
  !!
  !! Implementation is done by
  !! \f[
  !!     y' = y - q \Omega x (t-t_p),
  !! \f]
  !! with \f$ t_p = \text{NINT}(q\Omega t) / (q\Omega) \f$. In order to map
  !! the continuous shift at the discrete field linear interpolation is used,
  !! and assumes periodic behaviour along the y-direction.
#ifdef HAVE_FFTW
  SUBROUTINE FieldShift(this,Mesh,Physics,delt,field,mass2D)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)  :: Mesh
    CLASS(physics_base), INTENT(IN)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)  :: field
#if defined(PARALLEL)
    REAL, DIMENSION(1:Mesh%INUM,1:this%local_JNUM), &
                         INTENT(OUT) :: mass2D
#else
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
                         INTENT(OUT) :: mass2D
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
            mass2D(i,j-local_joff) = &
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
            mass2D(i,j) = &
             (1.0-this%jrem(i))*field(i,1+MODULO(j-1+FLOOR(this%joff(i)), &
              Mesh%JMAX-Mesh%JMIN+1),k) + this%jrem(i)*field(i,1+MODULO(j+ &
              FLOOR(this%joff(i)),Mesh%JMAX-Mesh%JMIN+1),k)
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE FieldShift
#endif
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
    DEALLOCATE(this%mass2D,this%Fmass2D)
#endif
    ! Free memory
    IF (ASSOCIATED(this%Fmass2D_real)) DEALLOCATE(this%Fmass2D)
    DEALLOCATE(&
               this%phi, &
               this%kx, this%ky, &
               this%joff, &
               this%jrem, &
               this%den_ip &
               )
#endif
    END SUBROUTINE Finalize

END MODULE gravity_sboxspectral_mod
