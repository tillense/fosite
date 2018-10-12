!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_sboxspectral.f03                                          #
!#                                                                           #
!# Copyright (C) 2015-2018                                                   #
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
!! - parameters of \link gravity_sboxspectral \endlink as key-values
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
!! There are currently two different implementations, which all need
!! different datatypes and array boundaries and optional plannings, etc.
!! This is why there there are lots of preprocessing statements. Generally
!! the following settings are possible. Please be aware that this module
!! requires a certain domain decomposition in pencils and cannot be chosen
!! arbitrarily.
!!
!! 1. FFTW on one core
!! 2. FFTW distributed parallel
!
!! \todo Remove unnecessary allocations like ip_den, etc. At the moment there
!!       are certainly possibilities to save memory, because for every task
!!       a new array is allocated.
!! \todo General cleanup: Many different array boundaries from 1:Mesh%INUM or
!!       Mesh%IMIN:Mesh%IMAX. This should be somehow coherently. Although
!!       Mesh%IMIN:Mesh%IMAX is used by Fosite, using Mesh%INUM.
!!
!! \note \b Concerning the parallel mode on SX Ace: \b Performance measurements
!!       showed that there is a large performance gain, if the individual
!!       strides are \b smaller \b than \f$ 256 \f$ cells. Simulations where
!!       performed on a \f$ 512 \times 512 \f$ grid and \f$ 1, 2, 4, 8 \f$ and
!!       \f$ 16 \f$ nodes (one core per node). Using \f$ 1 \f$ or \f$ 2 \f$
!!       nodes was magnitudes slower (!) than using \f$ 4 \f$ nodes. For even
!!       more nodes there was a typical performance gain.
!!
!! \attention The parallel version with this module needs a pencil decomposition
!!            of Mesh%INUM X Mesh%JNUM/nprocs for each process, in order to use
!!            the FFT-method with period boundaries in one direction.
!!
!! References:
!!
!! \cite gammie2001 , \cite gammiecode, \cite frigo2005, \cite mathkeisan2018
!----------------------------------------------------------------------------!
MODULE gravity_sboxspectral_mod
  USE gravity_base_mod
  USE boundary_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE logging_base_mod
  USE common_dict
#if defined(HAVE_FFTW)
  USE fftw
#endif
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
  CHARACTER(LEN=32), PARAMETER :: solver_name  = "sboxspectral"
  REAL, PARAMETER              :: SQRTTWOPI &
    = 2.50662827463100050241576528481104525300698674

  TYPE, EXTENDS(gravity_base) :: gravity_sboxspectral
    CHARACTER(LEN=32) :: gravity_name = "shearingbox spectral solver"
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
    REAL, DIMENSION(:,:,:), POINTER  :: Fmass2D_real   !< temporary variable
    INTEGER(C_INTPTR_T)              :: local_joff
    INTEGER, DIMENSION(:), POINTER   :: sizes
    INTEGER                          :: MNUM          !< number of modes
    REAL,DIMENSION(:), POINTER       :: kx            !< wave numbers for FFT (x)
    REAL,DIMENSION(:), POINTER       :: ky            !< wave numbers for FFT (y)
    REAL                             :: Lx, Ly
    REAL, DIMENSION(:), POINTER      :: joff, jrem   !< shifting indices (in SB)
#endif
#if defined(HAVE_FFTW) && defined(PARALLEL)
    INTEGER(C_INTPTR_T)              :: C_INUM, C_JNUM
    INTEGER(C_INTPTR_T)              :: alloc_local, local_JNUM
    TYPE(C_PTR)                      :: mass2D_pointer
    TYPE(C_PTR)                      :: Fmass2D_pointer
#endif


  CONTAINS

    PROCEDURE :: InitGravity_sboxspectral
    PROCEDURE :: UpdateGravity_single
    PROCEDURE :: InfoGravity
    PROCEDURE :: CalcPotential
    PROCEDURE :: FFT_Forward
    PROCEDURE :: FFT_Backward
    PROCEDURE :: CalcDiskHeight_single
    PROCEDURE :: FieldShift
    PROCEDURE :: FinalizeGravity_sboxspectral
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
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    TYPE(Dict_TYP),      POINTER    :: config,IO
    !------------------------------------------------------------------------!
    INTEGER             :: gravity_number
    INTEGER             :: err, valwrite, i
    CHARACTER(LEN=32)   :: info_str
#if defined(HAVE_FFTW) && defined(PARALLEL)
    INTEGER(C_INTPTR_T) :: C_INUM,C_JNUM
    INTEGER(C_INTPTR_T) :: alloc_local, local_JNUM, local_joff
    INTEGER             :: nprocs
#endif
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", gravity_number)
    CALL this%InitLogging(gravity_number,this%gravity_name)

    !>\todo{implement a check for equidistant mesh}
#if defined(HAVE_FFTW) && defined(PARALLEL)
    CALL fftw_mpi_init()
    C_INUM = Mesh%INUM
    C_JNUM = Mesh%JNUM
#endif

    !-------------- checks & warnings for initial conditions ----------------!
#if !defined(HAVE_FFTW) && !defined(HAVE_FFTW_LEGACY) && !defined(HAVE_FFTKEISAN)
    CALL this%Error("InitGravity_sboxspectral", &
         "No fft package could be loaded.")
#endif

#if defined(HAVE_FFTW)
    ! Check even number of cells
    IF(.NOT.(MOD(Mesh%JMAX-Mesh%JMIN+1,2)==0)) THEN
      CALL this%Error("InitGravity_sboxspectral", &
                 "The spectral poisson solver needs an even number of "// &
                 "cells in the phi direction due to the discrete cosinus "// &
                 "transform.")
    END IF
    !------------------------------------------------------------------------!

    ! check the dimensions if fftw should be used in parallel in order to
    ! know how to allocate the local arrays
#if (defined(HAVE_FFTW) && defined(PARALLEL))
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)
    IF (nprocs.GT.1 .AND. Mesh%WE_shear) THEN
      CALL this%Error("InitGravity_sboxspectral", &
                 "Execution of the parallel Fourier solver is only possible with pencil "// &
                 "decomposition. The first dimension should be fully accessible by the "//&
                 "first core. The shear thereto needs to vary along the second " // &
                 "dimension. In order to achieve this set the boundaries accordingly. "// &
                 "Check InitBoundary where the shear direction is set.")
    END IF
    IF (Mesh%dims(1).GT.1 .AND. Mesh%SN_shear) THEN
      CALL this%Error("InitGravity_sboxspectral", &
                 "The first dimension needs to be fully accessible by each process. "// &
                 "Please change domain decomposition to pencil decompositon by. " // &
                 "This needs to be done during initialization. In the Mesh dictionary " // &
                 "use the key 'decompositon' and the value (/ 1,-1/) to force decomposition " // &
                 "along second direction")
    END IF
#endif
    ALLOCATE( &
             this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KMAX), &
             this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DIM), &
#if defined(HAVE_FFTW) && !defined(PARALLEL)
             this%mass2D(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), & !TODO TODO TODO This needs to be 2D, make difference for 3D
             this%Fmass2D(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX), &
#elif defined(HAVE_FFTW) && defined(PARALLEL)
             ! initialization handled by FFTW (see below)
#endif
             this%Fmass2D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
             this%kx(Mesh%INUM), &
             this%ky(Mesh%JNUM), &
             this%joff(Mesh%IMIN:Mesh%IMAX), &
             this%jrem(Mesh%IMIN:Mesh%IMAX), &
             this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX), &
             STAT=err)
    IF (err.NE.0) &
        CALL this%Error("InitGravity_sboxspectral","Memory allocation failed.")

    this%local_joff = 0

    ! use special allocation pattern from fftw when using MPI in order to
    ! assure good alignment
#if defined(HAVE_FFTW) && defined(PARALLEL)
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
#if defined(HAVE_FFTW) && !defined(PARALLEL)
    this%plan_r2c = fftw_plan_dft_r2c_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1,this%mass2D, &
                                         this%Fmass2D,FFTW_MEASURE)
    this%plan_c2r = fftw_plan_dft_c2r_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1, this%Fmass2D, &
                                         this%mass2D, FFTW_MEASURE)
#elif defined(HAVE_FFTW) && defined(PARALLEL)
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
    this%accel(:,:,:,:) = 0.
    this%mass2D(:,:) = 0.
    this%Fmass2D(:,:) = CMPLX(0.,0)
    this%Fmass2D_real(:,:,:) = 0.
    this%joff(:) = 0.
    this%jrem(:) = 0.
    this%kx(:) = 0.
    this%ky(:) = 0.
    ! nullify not used potential explicitely
    NULLIFY(this%pot)

    ! precompute wavenumbers
    this%kx = cshift((/(i-(Mesh%INUM+1)/2,i=0,Mesh%INUM-1)/), &
                +(Mesh%INUM+1)/2)*2.*PI/(Mesh%XMAX-Mesh%XMIN)
    this%ky = cshift((/(i-(Mesh%JNUM+1)/2,i=0,Mesh%JNUM-1)/), &
                +(Mesh%JNUM+1)/2)*2.*PI/(Mesh%YMAX-Mesh%YMIN)

    !------------------------------- output ---------------------------------!
    valwrite = 0
    CALL GetAttr(config, "output/phi", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "phi", &
              this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
    valwrite = 0
    CALL GetAttr(config, "output/accel_x", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "accel_x", &
              this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,1))
    valwrite = 0
    CALL GetAttr(config, "output/accel_y", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "accel_y", &
              this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,2))
    valwrite = 0
    CALL GetAttr(config, "output/Fmass2D", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "Fmass2D_real", &
              this%Fmass2D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    CALL this%InitGravity(Mesh,Physics,config,IO)
#endif
  END SUBROUTINE InitGravity_sboxspectral

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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                                 INTENT(IN) :: pvar
    REAL,                        INTENT(IN) :: time,dt
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    REAL    :: w1,w2
    !------------------------------------------------------------------------!
    ! calc potential first
    CALL this%CalcPotential(Mesh,Physics,time,pvar)

    w1 = 3./48.
    w2 = 30./48.

    !\todo{Here a more robust approximation needs to be searched}
    ! Maybe Physics%VNUM is greater than 2 => set all to zero
    this%accel(:,:,:,:) = 0.
    DO k = Mesh%KMIN,Mesh%KMAX
      DO j = Mesh%JMIN,Mesh%JMAX
        DO i = Mesh%IMIN,Mesh%IMAX
          ! second order approximation
          this%accel(i,j,k,1) = -1.0*(this%phi(i+1,j,k)-this%phi(i-1,j,k))/ &
                               (2*Mesh%dlx(i,j,k))
          this%accel(i,j,k,2) = -1.0*(this%phi(i,j+1,k)-this%phi(i,j-1,k))/ &
                               (2*Mesh%dly(i,j,k))
          ! fourth order
!          this%accel(i,j,1) = -1.0*(w1*this%phi(i-2,j)-w2*this%phi(i-1,j)+ &
!                                    w2*this%phi(i+1,j)-w1*this%phi(i+2,j))/ &
!                                   (Mesh%dlx(i,j))
!          this%accel(i,j,2) = -1.0*(w1*this%phi(i,j-2)-w2*this%phi(i,j-1)+ &
!                                    w2*this%phi(i,j+1)-w1*this%phi(i,j+2)) / &
!                                   (Mesh%dly(i,j))
        END DO
      END DO
    END DO
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
  !!    See \link gravity_sboxspectral::fieldshift \endlink for more informations.
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
  !! The acceleration is eventually calculated in \link
  !! gravity_sboxspectral::updategravity_sboxspectral \endlink.
  SUBROUTINE CalcPotential(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    REAL,                INTENT(IN) :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN) :: pvar
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    REAL    :: kxp,kyp,K2max,K2,KO,kx,kern,ky
    REAL    :: joff2,jrem2,delt,time0,phi_old,phi_tmp,joff3
    INTEGER :: nx,ny,nprocs
    INTEGER :: iopt,ier
#ifdef PARALLEL
    INTEGER :: status(MPI_STATUS_SIZE), ierror
    REAL    :: mpi_buf(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM) !TODO: ONLY 2D
#endif
    !------------------------------------------------------------------------!
#if defined(HAVE_FFTW)
    !---------------- fourier transformation of density ---------------------!
    ! calculate the shift of the indice at time t                            !
    ! shift density to pretend periodic behavior with interpolation          !
    nprocs=1
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ier)
#endif

    IF (Mesh%OMEGA .EQ. 0) THEN
      time0 = 0.0
    ELSE IF (Mesh%WE_shear) THEN
      time0 = NINT(time*Mesh%Q*Mesh%OMEGA*(Mesh%xmax-Mesh%xmin)/ &
          (Mesh%ymax-Mesh%ymin))*(Mesh%ymax-Mesh%ymin)/(Mesh%Q*Mesh%OMEGA &
          *(Mesh%xmax-Mesh%xmin))
    ELSE IF (Mesh%SN_shear) THEN
      time0 = NINT(time*Mesh%Q*Mesh%OMEGA*(Mesh%ymax-Mesh%ymin)/ &
          (Mesh%xmax-Mesh%xmin))*(Mesh%xmax-Mesh%xmin)/(Mesh%Q*Mesh%OMEGA &
          *(Mesh%ymax-Mesh%ymin))
    END IF
    delt = time - time0

    !----------------- shift field to periodic point ------------------------!
!NEC$ IEXPAND
    CALL this%FieldShift(Mesh,Physics,delt, &
              pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DENSITY), &
              this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

! TODO POTENTIAL ERROR BECAUSE OF RESHAPE FUNCTIONS
#if defined(HAVE_FFTW) && defined(PARALLEL)
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
    !!      Do this pretty when everything else works.
!NEC$ IEXPAND
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

    DO j = Mesh%JMIN,Mesh%JMAX
#if defined(HAVE_FFTKEISAN) && defined(PARALLEL)
      DO i = Mesh%IMIN,Mesh%IMAX
#else
      DO i = Mesh%IMIN,Mesh%IMAX/2+1
#endif
        IF (Mesh%WE_shear) THEN
          K2 = (this%kx(i) + Mesh%Q*Mesh%OMEGA*this%ky(j)*delt)**2 + this%ky(j)**2
        ELSE IF (Mesh%SN_shear) THEN
          K2 = this%kx(i)**2 + (this%ky(j)-Mesh%Q*mesh%OMEGA*this%kx(i)*delt)**2
        END IF

        KO = SQRT(K2)
        IF ((K2 .LT. K2max) .AND. (K2 .GT. 0)) THEN
          this%Fmass2D(i,j-this%local_joff) = -2.*PI*Physics%Constants%GN* &
                                            this%Fmass2D(i,j-this%local_joff)/KO
        ELSE
          this%Fmass2D(i,j-this%local_joff) = 0.0
        END IF
      END DO
    END DO

    !----------- fourier transform (backward) of shifted density -----------!
#ifdef _FTRACE
CALL ftrace_region_begin("backward_fft")
#endif

    !> \todo Violation of best practices. It should take Fmass2D and mass2D as
    !!      dummy argument, but this is different for every combination.
    !!      Do this pretty when everything else works.
!NEC$ EXPAND
    CALL this%FFT_Backward(Mesh,Physics)

#ifdef _FTRACE
CALL ftrace_region_end("backward_fft")
#endif

    !------ calculate final potential with backshift and normalization ------!
    DO k = Mesh%KMIN,Mesh%KMAX
      DO j = Mesh%JMIN,Mesh%JMAX
        DO i = Mesh%IMIN,Mesh%IMAX
#if defined(HAVE_FFTW)
          this%phi(i,j,k) = this%mass2D(i,j-this%local_joff)/ &
               (Mesh%JNUM*Mesh%INUM)                        ! no norm. by FFTW !
#endif
        END DO
      END DO
    END DO
    !-------------------------- shift field back ----------------------------!
!NEC$ IEXPAND
    CALL this%FieldShift(Mesh,Physics,-delt, &
            this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
            this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    DO k = Mesh%KMIN,Mesh%KMAX
      DO j = Mesh%JMIN,Mesh%JMAX
        DO i = Mesh%IMIN,Mesh%IMAX
          this%phi(i,j,k) = this%den_ip(i,j,k)
        END DO
      END DO
    END DO

    !----- copy values to boundaries in order to calculate acceleration -----!

    IF(Mesh%WE_shear) THEN
!NEC$ NODEP
      DO j = 1,Mesh%GNUM
        ! southern northern (periodic)
        this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX-j+1,Mesh%KMIN:Mesh%KMAX)
        this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX) = this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN+j-1,Mesh%KMIN:Mesh%KMAX)
      END DO

      DO i = 1,Mesh%GNUM
        ! western (shorn periodic) - residual and integer shift
        joff2 = -Mesh%Q*Mesh%OMEGA*(Mesh%XMAX-Mesh%XMIN)*delt/Mesh%dy
        jrem2 = joff2 - FLOOR(joff2)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            this%phi(Mesh%IMIN-i,j,k) = &
                (1.0 - jrem2)*this%phi(Mesh%IMAX-i+1,j,k) + jrem2*this%phi(Mesh%IMAX-i+1,j+1,k)
          END DO
        END DO
        this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)  = &
              CSHIFT(this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
        ! eastern (shorn periodic) - residual and integer shift
        joff2 = Mesh%Q*Mesh%OMEGA*(Mesh%XMAX-Mesh%XMIN)*delt/Mesh%dy
        jrem2 = joff2 - FLOOR(joff2)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            this%phi(Mesh%IMAX+i,j,k) = &
                (1.0 - jrem2)*this%phi(Mesh%IMIN+i-1,j,k) + jrem2*this%phi(Mesh%IMIN+i-1,j+1,k)
          END DO
        END DO
        this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)  = &
              CSHIFT(this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
      END DO
! Only north-south direction has parallelization allowed
! \TODO ONLY 2D (Reshapes necessary below, when mpi_buf is in use)
! Approach below:
! Either with or without parallel mode apply peridic boundaries to all cells. In a second step shift the
! cells at the boundaries where it is necessary.
    ELSE IF(Mesh%SN_shear) THEN
#ifdef PARALLEL
      IF(Mesh%dims(2).GT.1) THEN
        mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GNUM) = &
        RESHAPE(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),(/ Mesh%IMAX-Mesh%IMIN+1,Mesh%GNUM /))
        CALL MPI_Sendrecv_replace(&
          mpi_buf,&
          2*(Mesh%IMAX-Mesh%IMIN+1), &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(NORTH), 53+NORTH, &
          Mesh%neighbor(SOUTH), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,Mesh%KMIN:Mesh%KMAX) = &
          RESHAPE(mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GNUM), (/ Mesh%IMAX-Mesh%IMAX+1, Mesh%GJNUM, Mesh%KMAX-Mesh%KMIN+1 /))

        mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GNUM) = &
        RESHAPE(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN+Mesh%GNUM-1:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX),(/ Mesh%IMAX-Mesh%IMIN+1,Mesh%GNUM /))
        CALL MPI_Sendrecv_replace(&
          mpi_buf,&
          2*(Mesh%IMAX-Mesh%IMIN+1), &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(SOUTH), 53+SOUTH, &
          Mesh%neighbor(NORTH), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX) = &
          RESHAPE(mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GNUM), (/ Mesh%IMAX-Mesh%IMAX+1, Mesh%GJNUM, Mesh%KMAX-Mesh%KMIN+1 /))
      ELSE
        DO j = 1,Mesh%GNUM
          ! southern northern (periodic in first step - further shift-treatment below)
          this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN-j,Mesh%KGMIN:Mesh%KGMAX) = this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX-j+1,Mesh%KGMIN:Mesh%KGMAX)
          this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+j,Mesh%KGMIN:Mesh%KGMAX) = this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN+j-1,Mesh%KGMIN:Mesh%KGMAX)
        END DO
      END IF
#else
      DO j = 1,Mesh%GNUM
        ! southern northern (periodic in first step - further shift-treatment below)
        this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN-j,Mesh%KGMIN:Mesh%KGMAX) = this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX-j+1,Mesh%KGMIN:Mesh%KGMAX)
        this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+j,Mesh%KGMIN:Mesh%KGMAX) = this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN+j-1,Mesh%KGMIN:Mesh%KGMAX)
      END DO
#endif
      DO i = 1,Mesh%GNUM
        ! western and eastern (always periodic)
        this%phi(Mesh%IMIN-i,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) = this%phi(Mesh%IMAX-i+1,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX)
        this%phi(Mesh%IMAX+i,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) = this%phi(Mesh%IMIN+i-1,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX)
      END DO

#ifdef PARALLEL
      IF (Mesh%mycoords(2).EQ.0) THEN
#endif
      DO j = 1,Mesh%GNUM
        ! southern (shorn periodic) - residual and integer shift
        joff2 = Mesh%Q*Mesh%OMEGA*(Mesh%YMAX-Mesh%YMIN)*delt/Mesh%dx
        jrem2 = joff2 - FLOOR(joff2)
        !------- residual shift ---------------------------------------------!
        DO k = Mesh%KMIN,Mesh%KMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMIN-j,k) = &
              (1.0 - jrem2)*this%phi(i,Mesh%JMIN-j,k) + jrem2*this%phi(i+1,Mesh%JMIN-j,k)
          END DO
        END DO
        !------- integral shift ---------------------------------------------!
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX)  = &
          CSHIFT(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
      END DO
#ifdef PARALLEL
      END IF
#endif

#ifdef PARALLEL
      IF (Mesh%mycoords(2).EQ.Mesh%dims(2)-1) THEN
#endif
      DO j = 1,Mesh%GNUM
        ! northern (shorn periodic) - residual and integer shift
        joff2 = -Mesh%Q*Mesh%OMEGA*(Mesh%YMAX-Mesh%YMIN)*delt/Mesh%dx
        jrem2 = joff2 - FLOOR(joff2)
        !------- residual shift ---------------------------------------------!
        DO k = Mesh%KMIN,Mesh%KMAX
          DO i = Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMAX+j,k) = &
              (1.0 - jrem2)*this%phi(i,Mesh%JMAX+j,k) + jrem2*this%phi(i+1,Mesh%JMAX+j,k)
          END DO
        END DO
        !------- integral shift ---------------------------------------------!
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX)  = &
          CSHIFT(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX),FLOOR(joff2))
      END DO
#ifdef PARALLEL
      END IF
#endif
    END IF
#endif
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
    INTEGER :: ier,iopt
#ifdef PARALLEL
    INTEGER :: nprocs
    INTEGER :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ier)
#endif

#ifdef _FTRACE
CALL ftrace_region_begin("forward FFT")
#endif

#if defined(HAVE_FFTW) && !defined(PARALLEL)
    CALL fftw_execute_dft_r2c(this%plan_r2c, this%mass2D, this%Fmass2D)
    ! turn complex array to real one for output
    DO k = Mesh%KMIN,Mesh%KMAX
      DO j = Mesh%JMIN,Mesh%JMAX
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          this%Fmass2D_real(2*i-1,j,k) = REAL(REAL(this%Fmass2D(i,j)))
          this%Fmass2D_real(2*i,j,k)   = REAL(AIMAG(this%Fmass2D(i,j)))
        END DO
      END DO
    END DO
#elif defined(HAVE_FFTW) && defined(PARALLEL)
    CALL fftw_mpi_execute_dft_r2c(this%plan_r2c,this%mass2D, this%Fmass2D)
    DO k = Mesh%KMIN,Mesh%KMAX
      DO j = Mesh%JMIN,Mesh%JMAX
        DO i = Mesh%IMIN,Mesh%IMAX/2+1
          this%Fmass2D_real(2*i-1,j,k) = REAL(REAL(this%Fmass2D(i,j-this%local_joff)))
          this%Fmass2D_real(2*i,j,k)   = REAL(AIMAG(this%Fmass2D(i,j-this%local_joff)))
        END DO
      END DO
    END DO
#endif

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
    INTEGER           :: i,j
    INTEGER           :: ier,iopt
#ifdef PARALLEL
    INTEGER           :: nprocs
    INTEGER           :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ier)
#endif

#if defined(HAVE_FFTW) && !defined(PARALLEL)
    CALL fftw_execute_dft_c2r(this%plan_c2r, this%Fmass2D, this%mass2D)
#elif defined(HAVE_FFTW) && defined(PARALLEL)
    CALL fftw_mpi_execute_dft_c2r(this%plan_c2r,this%Fmass2D, this%mass2D)
#endif
  END SUBROUTINE



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
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX), &
                         INTENT(IN)    :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(OUT)   :: h_ext
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX), &
                         INTENT(INOUT) :: height
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    REAL              :: cs2,p,q
    !------------------------------------------------------------------------!
    ! pure self-gravitating shearing sheet with external point mass potential
!NEC$ COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        DO i=Mesh%IGMIN,Mesh%IGMAX
           cs2 = bccsound(i,j,k)*bccsound(i,j,k)
           p = -SQRTTWOPI*Physics%Constants%GN*pvar(i,j,k,Physics%DENSITY) &
                  /Mesh%OMEGA**2.
           q = cs2/Mesh%OMEGA**2.
           ! return the new disk height
           height(i,j,k) = p+SQRT(q+p*p)
        END DO
      END DO
    END DO

  END SUBROUTINE CalcDiskHeight_single

  !> Prints out information
  SUBROUTINE InfoGravity(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(IN) :: this
    CLASS(mesh_base),            INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
#if defined(HAVE_FFTW)
    CALL this%Info( "            FFT-Package:       FFTW " // &
#if defined(PARALLEL)
                                                    "- parallel mode")
#elif !defined(PARALLEL)
                                                    "- serial mode")
#endif
#endif
    CALL this%Info("            .. done initializing")
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
  SUBROUTINE FieldShift(this,Mesh,Physics,delt,field,mass2D)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)  :: Mesh
    CLASS(physics_base), INTENT(IN)  :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)  :: field
#if defined(HAVE_FFTW) && defined(PARALLEL)
    REAL, DIMENSION(1:Mesh%INUM,1:this%local_JNUM), &
                         INTENT(OUT) :: mass2D
#else
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
                         INTENT(OUT) :: mass2D
#endif
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,local_joff
    REAL              :: cs2,p,q,delt
    REAL              :: mass2D_old, mass2D_tmp
    !------------------------------------------------------------------------!
    local_joff = this%local_joff

    IF (Mesh%SN_shear) THEN
      DO k = Mesh%KMIN,Mesh%KMAX
        DO j = Mesh%JMIN,Mesh%JMAX
          this%joff(j)   = Mesh%Q*Mesh%OMEGA*Mesh%bcenter(Mesh%IMIN,j,k,2)* &
                      delt/Mesh%dx
          this%jrem(j)  = this%joff(j) - FLOOR(this%joff(j))
          DO i = Mesh%IMIN,Mesh%IMAX
            mass2D(i,j-local_joff) = &
             (1.0-this%jrem(j))*field(1+MODULO(i-1+FLOOR(this%joff(j)), &
              Mesh%IMAX-Mesh%IMIN+1),j,k) + this%jrem(j)*field(1+MODULO(i+ &
              FLOOR(this%joff(j)),Mesh%IMAX-Mesh%IMIN+1),j,k)
          END DO
        END DO
      END DO
    ELSE IF (Mesh%WE_shear) THEN
      DO k = Mesh%KMIN,Mesh%KMAX
        DO i = Mesh%IMIN,Mesh%IMAX
          this%joff(i)   = -Mesh%Q*Mesh%OMEGA*Mesh%bcenter(i,Mesh%JMIN,k,1)* &
                      delt/Mesh%dy
          this%jrem(i)  = this%joff(i) - FLOOR(this%joff(i))
        END DO
!NEC$ COLLAPSE
        DO j = Mesh%JMIN,Mesh%JMAX
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

  !> \public Closes the gravity term of the shearingsheet spectral solver.
  SUBROUTINE FinalizeGravity_sboxspectral(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   CLASS(gravity_sboxspectral) :: this
   !------------------------------------------------------------------------!
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY) || defined(HAVE_FFTKEISAN)
    ! Destroy plans
#ifdef HAVE_FFTW
    CALL fftw_destroy_plan(this%plan_r2c)
    CALL fftw_destroy_plan(this%plan_c2r)
#endif
#if defined(HAVE_FFTW) && defined(PARALLEL)
!    fftw_free(this%mass2D,this%Fmass2D)
#endif
    ! Free memory
    DEALLOCATE(&
               this%phi, &
               this%accel, &
               this%mass2D, &
               this%Fmass2D, &
               this%Fmass2D_real, &
               this%joff, &
               this%jrem, &
               this%den_ip &
               )
#endif
    END SUBROUTINE FinalizeGravity_sboxspectral

END MODULE gravity_sboxspectral_mod
