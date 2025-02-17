!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_sboxspectral.f90                                          #
!#                                                                           #
!# Copyright (C) 2015-2024                                                   #
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
!!  \key{output/potential,INTEGER,(enable=1) output gravitational potential,0}
!!  \key{output/Fmass3D,INTEGER,(enable=1) output z-slices of 2D (x-y) fourier
!!                              transformed density,0}
!----------------------------------------------------------------------------!
!> \author Jannes Klee
!! \author Tobias Illenseer
!!
!! \brief Poisson solver using spectral methods within the shearingbox.
!! \parblock
!! There are currently two different implementations. Please be aware that this module
!! requires a certain domain decomposition in pencils and cannot be chosen
!! arbitrarily.
!!
!! 1. FFTW on one core
!! 2. FFTW distributed parallel
!! \attention The parallel execution of this spectral solver requires pencil
!!            decomposition of the computation domain, i.e., splitting
!!            is only allowed in 2nd and/or 3rd dimension. Check the key
!!            "decomposition" when setting the mesh properties (see \ref mesh ).
!! \endparblock
!! #### References
!! \cite gammie2001 , \cite gammiecode, \cite frigo2005
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
  TYPE Field_TYP
    REAL(C_DOUBLE), DIMENSION(:,:), POINTER :: original
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:), POINTER :: transformed
  END TYPE Field_TYP

  TYPE, EXTENDS(gravity_base) :: gravity_sboxspectral
#ifdef HAVE_FFTW
    TYPE(C_PTR)           :: plan_r2c                  !> plan for real to complex FT
    TYPE(C_PTR)           :: plan_c2r                  !> plan for complex to real FT
    TYPE(Field_TYP), DIMENSION(:), ALLOCATABLE &       !< z-layered density/potential field
                          :: field                     !<   and their fourier transformed counterparts
    REAL(C_DOUBLE), DIMENSION(:,:), POINTER, CONTIGUOUS &
                          :: Kxy2 => null()            !< length of wave vectors squared
    COMPLEX(C_DOUBLE_COMPLEX), POINTER, CONTIGUOUS &
                          :: Fmass3D(:,:,:) => null(),&!< fourier transformed 3D mass density
                             Fsum3D(:,:,:) => null(), &!< sum used for reconstruction
                                                       !<   of transformed potential (3D only)
                             qk(:,:,:) => null()       !< weight factors (3D only)
    REAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS &
                          :: Fmass3D_real => null(), & !< just for output
                             den_ip => null(), &       !< interpolated 3D density
                             phi => null()             !< potential
    INTEGER(C_INTPTR_T)              :: local_joff
    REAL,DIMENSION(:), POINTER &
                          :: kx => null(),ky => null() !< x/y wave numbers for FFT
    REAL                             :: shiftconst     !< constant for shift
    REAL                             :: maxKxy2        !< max
    REAL, DIMENSION(:), POINTER &
                          :: joff=>null(),jrem=>null() !< shifting indices (in SB)
    INTEGER                          :: order
#ifdef PARALLEL
    REAL, DIMENSION(:,:,:), POINTER  :: mpi_buf => null()
    INTEGER(C_INTPTR_T)              :: C_INUM, C_JNUM
    INTEGER(C_INTPTR_T)              :: alloc_local, local_JNUM
    TYPE(C_PTR)                      :: fftw_real_pointer, &
                                        fftw_complex_pointer
#endif
#endif

  CONTAINS

    PROCEDURE :: InitGravity_sboxspectral
    PROCEDURE :: UpdateGravity_single
!     PROCEDURE :: InfoGravity
    PROCEDURE :: SetOutput
    PROCEDURE :: CalcDiskHeight_single
    FINAL :: Finalize
#ifdef HAVE_FFTW
    PROCEDURE :: CalcPotential
    PROCEDURE :: FFT_Forward
    PROCEDURE :: FFT_Backward
    PROCEDURE :: SetBoundaryData
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
    INTEGER             :: err, valwrite, i,j,k
#if defined(HAVE_FFTW) && defined(PARALLEL)
    INTEGER(C_INTPTR_T) :: C_INUM,C_JNUM
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
    CASE(2)
      ! ok, do nothing
    CASE(4)
      IF (Mesh%KNUM.GT.1) &
        CALL this%Error("InitGravity_sboxspectral", &
          "order 4 is currently not supported in 3D")
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
      CALL this%Warning("InitGravity_sboxspectral","sboxspectral 3D is experimental")
#endif
    END IF

    ALLOCATE( &
             this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
             this%qk(Mesh%IMIN/2+1:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,0:Mesh%KNUM+1),&
             this%field(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
#ifndef PARALLEL
             ! allocate this field and set pointer references into z-slices (see below)
             ! in parallel mode initialization handled by FFTW (see below)
             this%Fmass3D(Mesh%IMIN/2+1:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
#endif
!!!! TODO: only allocate these for KNUM > 1
             this%Fsum3D(Mesh%IMIN/2+1:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
!!!! END TODO
             this%Kxy2(Mesh%IMIN/2+1:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX), &
             this%kx(Mesh%INUM), &
             this%ky(Mesh%JNUM), &
             this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
             STAT=err)
    IF (err.NE.0) &
        CALL this%Error("InitGravity_sboxspectral","Memory allocation failed.")

#ifdef HAVE_FFTW
    this%local_joff = 0
#endif

    ! use special allocation pattern from fftw when using MPI in order to
    ! assure good alignment
#if defined(PARALLEL)
    ALLOCATE(this%mpi_buf(Mesh%IMIN:Mesh%IMAX,Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX), &
             STAT=err)
    IF (err.NE.0) &
        CALL this%Error("InitGravity_sboxspectral","Memory allocation failed for mpi_buf.")

    this%alloc_local     = fftw_mpi_local_size_2d(C_JNUM, C_INUM, &
                              MPI_COMM_WORLD, this%local_JNUM, this%local_joff)
    DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
      this%fftw_real_pointer    = fftw_alloc_real(2*this%alloc_local)
      this%fftw_complex_pointer = fftw_alloc_complex(this%alloc_local)
      CALL c_f_pointer(this%fftw_real_pointer,this%field(k)%original, &
                      [2*(C_INUM/2+1),this%local_JNUM])
      CALL c_f_pointer(this%fftw_complex_pointer,this%field(k)%transformed, &
                      [C_INUM/2+1,this%local_JNUM])
    END DO
#endif

    !------------------- create plans for fftw ------------------------------!
    ! Pay attention to the argument order of the dimension (JNUM and INUM    !
    ! are switched because of C -> row-major, Fortran -> column-major),      !
    ! BUT ONLY in modern Fortran UNLIKE the legacy version                   !
    ! ------------ plans are allocated in dictionary ------------------------!
    CALL this%Info("            initializing FFTW: " // &
#if  !defined(PARALLEL)
                   "serial mode")
    DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
      this%field(k)%original => this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k)
      this%field(k)%transformed => this%Fmass3D(Mesh%IMIN/2+1:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX,k)
    END DO
    this%plan_r2c = fftw_plan_dft_r2c_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1,this%field(Mesh%KMIN)%original, &
                                         this%field(Mesh%KMIN)%transformed,FFTW_MEASURE)
    this%plan_c2r = fftw_plan_dft_c2r_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1, this%field(Mesh%KMIN)%transformed, &
                                         this%field(Mesh%KMIN)%original, FFTW_MEASURE)
    IF ((.NOT. C_ASSOCIATED(this%plan_r2c)) .OR. (.NOT. c_associated(this%plan_c2r))) THEN
       CALL this%Error("InitGravity_sboxspectral","FFT plan could not be created.")
    END IF
#elif defined(PARALLEL)
                   "parallel mode")
    this%plan_r2c = fftw_mpi_plan_dft_r2c_2d(C_JNUM,C_INUM,this%field(Mesh%KMIN)%original, &
                                             this%field(Mesh%KMIN)%transformed, &
                                             MPI_COMM_WORLD, FFTW_MEASURE)
    this%plan_c2r = fftw_mpi_plan_dft_c2r_2d(C_JNUM,C_INUM, this%field(Mesh%KMIN)%transformed, &
                                             this%field(Mesh%KMIN)%original, &
                                             MPI_COMM_WORLD, FFTW_MEASURE)
#endif
    CALL GetAttr(config, "print_plans", valwrite, 0)
    IF (valwrite.GT.0) THEN
#ifndef NECSXAURORA
      IF (this%GetRank().EQ.0) THEN
        CALL fftw_print_plan(this%plan_r2c)
        PRINT *,ACHAR(13)
        CALL fftw_print_plan(this%plan_c2r)
        PRINT *,ACHAR(13)
      END IF
#else
      CALL this%Warning("InitGravity_sboxspectral", "fftw_print_plan currently not supported on SX-Aurora",0)
#endif
    END IF

    !------------------------------------------------------------------------!
    ! initialize arrays
    this%phi(:,:,:) = 0.
    DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
      this%field(k)%transformed(:,:) = CMPLX(0.,0.)
      ! in serial mode this also zeros this%Fmass3D
      this%Fsum3D(:,:,k) = CMPLX(0.,0)
      this%qk(:,:,k) = 0.
    END DO
    this%kx(:) = 0.
    this%ky(:) = 0.

    ! precompute wave numbers
    this%kx(:) = CSHIFT((/(i-(Mesh%INUM+1)/2,i=0,Mesh%INUM-1)/),(Mesh%INUM+1)/2) &
                   *2.*PI/(Mesh%XMAX-Mesh%XMIN)
    this%ky(:) = CSHIFT((/(j-(Mesh%JNUM+1)/2,j=0,Mesh%JNUM-1)/),(Mesh%JNUM+1)/2) &
                   *2.*PI/(Mesh%YMAX-Mesh%YMIN)

    ! cut-off value for wave numbers
    this%maxKxy2 = 0.5*MIN(PI/Mesh%dx,PI/Mesh%dy)**2

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

    ! register data fields for output
    CALL this%SetOutput(Mesh,Physics,config,IO)
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
            IF (Physics%VDIM.EQ.3) &
              this%accel%data4d(i,j,k,3) = -1.0*(this%phi(i,j,k+1)-this%phi(i,j,k-1))/ &
                                 (2*Mesh%dlz%data3d(i,j,k))
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
            IF (Physics%VDIM.EQ.3) &
              this%accel%data4d(i,j,k,3) = -1.0*(w1*this%phi(i,j,k-2)-w2*this%phi(i,j,k-1)+ &
                                        w2*this%phi(i,j,k+1)-w1*this%phi(i,j,k+2)) / &
                                       (Mesh%dlz%data3d(i,j,k))
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
    REAL    :: delt,time0
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
                this%field(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1))
    CLASS DEFAULT
      CALL this%Error("gravity_sboxspectral::CalcPotential","unsupported state vector")
    END SELECT

    !-------------- fourier transform (forward) of shifted density ----------!
! for performance checks
#ifdef _FTRACE
CALL ftrace_region_begin("forward_fft")
#endif

    !> \todo Violation of best practices. It should take Fmass2D and mass2D as
    !!      dummy argument, but this is different for every combination
!!NEC$ IEXPAND
    CALL this%FFT_Forward(Mesh,Physics,this%field(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1))

! performance checks
#ifdef _FTRACE
CALL ftrace_region_end("forward_fft")
#endif

    !------------- calculate wave number vector squared --------------------!
    IF (Mesh%shear_dir.EQ.2) THEN
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
          this%Kxy2(i,j) = (this%kx(i) + Mesh%Q*Mesh%OMEGA*this%ky(j)*delt)**2 + this%ky(j)**2
        END DO
      END DO
    ELSE IF (Mesh%shear_dir.EQ.1) THEN
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
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
        DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
          IF ((this%Kxy2(i,j).LT.this%maxKxy2).AND.(this%Kxy2(i,j).GT. 0)) THEN
              this%field(k)%transformed(i,j-this%local_joff) = -2.*PI*Physics%Constants%GN &
                *this%field(k)%transformed(i,j-this%local_joff)/SQRT(this%Kxy2(i,j))
          ELSE
              this%field(k)%transformed(i,j-this%local_joff) = 0.0
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
          DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
            IF ((this%Kxy2(i,j).LT.this%maxKxy2).AND.(this%Kxy2(i,j).GT. 0)) THEN
              this%field(k)%transformed(i,j-this%local_joff) = -4.*PI*Physics%Constants%GN &
                *(1.-this%qk(i,j,0))/this%Kxy2(i,j) * this%field(k)%transformed(i,j-this%local_joff)
            ELSE
              this%field(k)%transformed(i,j-this%local_joff) = 0.0
            END IF
          END DO
        END DO
      ELSE
        ! 3D disk with more than one vertical zone
        this%qk(:,:,1) = this%qk(:,:,0)*this%qk(:,:,0)
        DO k=2,Mesh%KNUM+1
          this%qk(:,:,k) = this%qk(:,:,k-1) * this%qk(:,:,1)
        END DO
        DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
          this%Fsum3D(:,:,k) = 0.0
          ! the following operations are equivalent to a matrix vector multiplication
          ! where the matrix is a symmetric band matrix with zeros on the main diagonal
          DO kk=Mesh%KMIN-Mesh%KP1,k-1
!NEC$ IVDEP
            DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
              DO i = Mesh%IMIN,Mesh%IMAX/2+1
                this%Fsum3D(i,j,k) = this%Fsum3D(i,j,k) + this%qk(i,j,ABS(k-kk)) &
                                   * this%field(kk)%transformed(i,j-this%local_joff)
              END DO
            END DO
          END DO
          ! skip k=kk
!NEC$ IVDEP
          DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
              DO kk=k+1,Mesh%KMAX+Mesh%KP1
                this%Fsum3D(i,j,k) = this%Fsum3D(i,j,k) + this%qk(i,j,ABS(k-kk)) &
                                   * this%field(kk)%transformed(i,j-this%local_joff)
              END DO
            END DO
          END DO
        END DO
        DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
!NEC$ IVDEP
          DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
              IF ((this%Kxy2(i,j).LT.this%maxKxy2).AND.(this%Kxy2(i,j).GT. 0)) THEN
                this%field(k)%transformed(i,j-this%local_joff) = -4.*PI*Physics%Constants%GN &
                  * (1.-this%qk(i,j,0))/this%Kxy2(i,j) &
                  * (this%field(k)%transformed(i,j-this%local_joff) + 0.5*(1.0+1./this%qk(i,j,0)) &
                    * this%Fsum3D(i,j,k))
              ELSE
                this%field(k)%transformed(i,j-this%local_joff) = 0.0
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF

    !----------- fourier transform (backward) of shifted potential  ---------!
#ifdef _FTRACE
CALL ftrace_region_begin("backward_fft")
#endif

    ! take the fourier transformed potential this%field(:)%transformed,
    ! apply the inverse transforms and store the result in this%field(:)%original
    CALL this%FFT_Backward(Mesh,Physics,this%field(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1))

#ifdef _FTRACE
CALL ftrace_region_end("backward_fft")
#endif

    !------ calculate final potential with backshift and normalization ------!
!NEC$ IVDEP
    DO k = Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX
          this%phi(i,j,k) = this%field(k)%original(i,j-this%local_joff)/ &
               (Mesh%JNUM*Mesh%INUM)                        ! no norm. by FFTW !
        END DO
      END DO
    END DO
    !-------------------------- shift field back ----------------------------!
!!NEC$ IEXPAND
    CALL this%FieldShift(Mesh,Physics,-delt, &
            this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
            this%field(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1))

!NEC$ IVDEP
    DO k = Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
!NEC$ IVDEP
      DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX
          this%phi(i,j,k) = this%field(k)%original(i,j-this%local_joff)
        END DO
      END DO
    END DO

    CALL this%SetBoundaryData(Mesh,Physics,delt)
  END SUBROUTINE CalcPotential

  !> \public Set ghost cell data for the potential
  !!
  !! Parallelization in x-y-plane is only supported for shear along 1st dimension,
  !! i.e., in west-east direction, so we need MPI communication only in this case.
  SUBROUTINE SetBoundaryData(this,Mesh,Physics,delt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT):: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    REAL,                INTENT(IN) :: delt
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    REAL    :: joff,jrem
#ifdef PARALLEL
    INTEGER :: status(MPI_STATUS_SIZE),ierror
#endif
    !------------------------------------------------------------------------!
    IF(Mesh%shear_dir.EQ.2) THEN
      ! southern / northern boundary: periodic
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ SHORTLOOP
        DO j = 1,Mesh%GJNUM
          DO i = Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMIN-j,k) = this%phi(i,Mesh%JMAX-j+1,k)
            this%phi(i,Mesh%JMAX+j,k) = this%phi(i,Mesh%JMIN+j-1,k)
          END DO
        END DO
      END DO

      ! western / eastern boundary: shifted periodic
      joff = -this%shiftconst*delt
      jrem = joff - FLOOR(joff)
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ SHORTLOOP
        DO i = 1,Mesh%GINUM
          DO j = Mesh%JMIN,Mesh%JMAX
            ! copy eastern data
            this%phi(Mesh%IMIN-i,j,k) = this%phi(Mesh%IMAX-i+1,j,k)
          END DO
          ! apply integer shift
          this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,k) = &
            CSHIFT(this%phi(Mesh%IMIN-i,Mesh%JMIN:Mesh%JMAX,k),FLOOR(joff))
          this%phi(Mesh%IMIN-i,Mesh%JMAX+1,k) = this%phi(Mesh%IMIN-i,Mesh%JMIN,k)
          ! apply residual shift (interpolation)
          DO j = Mesh%JMIN,Mesh%JMAX
            this%phi(Mesh%IMIN-i,j,k) = &
                (1.0 - jrem)*this%phi(Mesh%IMIN-i,j,k) + jrem*this%phi(Mesh%IMIN-i,j+1,k)
          END DO
        END DO
      END DO

      joff = this%shiftconst*delt
      jrem = joff - FLOOR(joff)
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ SHORTLOOP
        DO i = 1,Mesh%GINUM
          DO j = Mesh%JMIN,Mesh%JMAX
            ! copy western data
            this%phi(Mesh%IMAX+i,j,k) = this%phi(Mesh%IMIN+i-1,j,k)
          END DO
          ! apply integer shift
          this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,k) = &
            CSHIFT(this%phi(Mesh%IMAX+i,Mesh%JMIN:Mesh%JMAX,k),FLOOR(joff))
          this%phi(Mesh%IMAX+i,Mesh%JMAX+1,:) = this%phi(Mesh%IMAX+i,Mesh%JMIN,:)
          ! apply residual shift (interpolation)
          DO j = Mesh%JMIN,Mesh%JMAX
            this%phi(Mesh%IMAX+i,j,k) = &
                (1.0 - jrem)*this%phi(Mesh%IMAX+i,j,k) + jrem*this%phi(Mesh%IMAX+i,j+1,k)
          END DO
        END DO
      END DO
    ELSE IF(Mesh%shear_dir.EQ.1) THEN
      ! western / eastern boundary: periodic (no MPI communication)
      DO k = Mesh%KMIN,Mesh%KMAX
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ SHORTLOOP
          DO i = 1,Mesh%GJNUM
            this%phi(Mesh%IMIN-i,j,k) = this%phi(Mesh%IMAX-i+1,j,k)
            this%phi(Mesh%IMAX+i,j,k) = this%phi(Mesh%IMIN+i-1,j,k)
          END DO
        END DO
      END DO
#ifdef PARALLEL
      ! ATTENTION: send order is important for MPI:
      !  1. send non-shifted direction
      !  2. send shifted direction
      IF(Mesh%dims(2).GT.1) THEN
        this%mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) = &
          this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-Mesh%GJNUM+1:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX)
        CALL MPI_Sendrecv_replace(&
          this%mpi_buf,&
          SIZE(this%mpi_buf), &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(NORTH), 53+NORTH, &
          Mesh%neighbor(SOUTH), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)

        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,Mesh%KMIN:Mesh%KMAX) = &
          this%mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX)

        this%mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX) = &
          this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMIN+Mesh%GJNUM-1,Mesh%KMIN:Mesh%KMAX)

        CALL MPI_Sendrecv_replace(&
          this%mpi_buf,&
          SIZE(this%mpi_buf), &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(SOUTH), 53+SOUTH, &
          Mesh%neighbor(NORTH), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX) = &
          this%mpi_buf(Mesh%IMIN:Mesh%IMAX,1:Mesh%GJNUM,Mesh%KMIN:Mesh%KMAX)
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
        joff = this%shiftconst*delt
        jrem = joff - FLOOR(joff)
        !------- integral shift ---------------------------------------------!
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX)  = &
          CSHIFT(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN-j,Mesh%KMIN:Mesh%KMAX),FLOOR(joff))
        !------- residual shift ---------------------------------------------!
        this%phi(Mesh%IMAX+1,Mesh%JMIN-j,:) = this%phi(Mesh%IMIN,Mesh%JMIN-j,:)
        DO k = Mesh%KMIN,Mesh%KMAX
          DO i=Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMIN-j,k) = &
              (1.0 - jrem)*this%phi(i,Mesh%JMIN-j,k) + jrem*this%phi(i+1,Mesh%JMIN-j,k)
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
        joff = -this%shiftconst*delt
        jrem = joff - FLOOR(joff)
        !------- integral shift ---------------------------------------------!
        this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX)  = &
          CSHIFT(this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+j,Mesh%KMIN:Mesh%KMAX),FLOOR(joff))
        !------- residual shift ---------------------------------------------!
        this%phi(Mesh%IMAX+1,Mesh%JMAX+j,:) = this%phi(Mesh%IMIN,Mesh%JMAX+j,:)
        DO k = Mesh%KMIN,Mesh%KMAX
          DO i = Mesh%IMIN,Mesh%IMAX
            this%phi(i,Mesh%JMAX+j,k) = &
              (1.0 - jrem)*this%phi(i,Mesh%JMAX+j,k) + jrem*this%phi(i+1,Mesh%JMAX+j,k)
          END DO
        END DO
      END DO
#ifdef PARALLEL
      END IF
#endif
    END IF
  END SUBROUTINE SetBoundaryData

  !> \private Calculates the FFT forward
  SUBROUTINE FFT_Forward(this,Mesh,Physics,field)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT):: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    TYPE(Field_TYP), DIMENSION(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
                      INTENT(INOUT) :: field
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    !------------------------------------------------------------------------!

#ifdef _FTRACE
CALL ftrace_region_begin("forward FFT")
#endif

    DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
#if defined(PARALLEL)
      CALL fftw_mpi_execute_dft_r2c(this%plan_r2c,field(k)%original, &
                                  this%field(k)%transformed)
#else
      CALL fftw_execute_dft_r2c(this%plan_r2c,this%field(k)%original, &
                                this%field(k)%transformed)
#endif
    END DO

    IF (ASSOCIATED(this%Fmass3D_real)) THEN
!NEC$ IVDEP
      DO k = Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i = Mesh%IMIN/2+1,Mesh%IMAX/2+1
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
  SUBROUTINE FFT_Backward(this,Mesh,Physics,field)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    TYPE(Field_TYP), DIMENSION(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
                      INTENT(INOUT) :: field
    !------------------------------------------------------------------------!
    INTEGER           :: k
    !------------------------------------------------------------------------!
    DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
#if defined(PARALLEL)
      CALL fftw_mpi_execute_dft_c2r(this%plan_c2r,field(k)%transformed, &
                                    field(k)%original)
#else
      CALL fftw_execute_dft_c2r(this%plan_c2r,field(k)%transformed, &
                                field(k)%original)
#endif
    END DO
  END SUBROUTINE
#endif


  !> \public Compute disk pressure scale height for geometrically thin
  !! self-gravitating shearingsheet.
  !!
  !! For razor-thin 2D disks the algorithm solves the equations
  !! \f{eqnarray}{
  !!     1/h^2 &=& 1/h_{sg}^2 + 1/h_{ext}^2 & (1)\\
  !!     \left.\frac{\partial^2\Phi_{sg}}{\partial z^2}\right|_{z=0} =
  !!           (c_s/h_{sg})^2 &=& 4\pi G \rho_c & (2) \\
  !!     \rho_c &=& 1/\sqrt{2\pi} \Sigma / h & (3) \\
  !!     1/h_{ext}^2 &=& \Omega_{K}^2 / c_s^2
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
  !! \f$ q = c_s^2 / \Omega_{K}^2 \f$.
  !! The sign of the root has to be positive, because in the non-selfgravitating
  !! limit \f$ p=0, q=1 \f$ and \f$ h \f$ needs to be
  !! \f$ c_s/\Omega_{K} \f$ not \f$ -c_s/\Omega_{K} \f$.
  !!
  !! In the 3D case the midplane density \f$ \rho_c \f$ can be obtained
  !! directly from the data. Thus there is no need for equation (3) and
  !! the scale height is simply given by
  !! \f[
  !!     h = \frac{c_s}{\sqrt{4\pi G \rho_c + \Omega_{K}^2}}.
  !! \f]
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
    REAL              :: cs2,p,q,rho_c
    !------------------------------------------------------------------------!
    ! pure self-gravitating shearing sheet with external point mass potential
    IF (ABS(Mesh%zmax-Mesh%zmin).LT.TINY(Mesh%zmin)) THEN
      ! 2D razor thin disk
      k = Mesh%KMIN
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
    ELSE
      k = Mesh%KNUM / 2
#ifdef PARALLEL
      IF (k.LT.Mesh%KMIN.OR.k.GT.Mesh%KMAX) &
        CALL this%Error("gravity_sboxspectral::CalcDiskHeight_single","not yet parallelized")
#endif
!NEC$ IVDEP
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! check for even number grid cells in z-direction
          IF (MOD(Mesh%KNUM,2).EQ.0) THEN
            ! average density around midplane
            rho_c = 0.5*(pvar%data4d(i,j,k,Physics%DENSITY)+pvar%data4d(i,j,k+1,Physics%DENSITY))
          ELSE
            ! take midplane density
            rho_c = pvar%data4d(i,j,k+1,Physics%DENSITY)
          END IF
          ! return the new disk height
          height%data3d(i,j,Mesh%KGMIN:Mesh%KGMAX) = bccsound%data3d(i,j,k) &
            / SQRT(4*PI*Physics%Constants%GN * rho_c + Mesh%OMEGA**2)
        END DO
      END DO
    END IF
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
#ifdef HAVE_FFTW
  SUBROUTINE FieldShift(this,Mesh,Physics,delt,field,shifted_field)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_sboxspectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)  :: Mesh
    CLASS(physics_base), INTENT(IN)  :: Physics
    REAL,                INTENT(IN)  :: delt
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)  :: field
    TYPE(Field_TYP), DIMENSION(Mesh%KMIN-Mesh%KP1:Mesh%KMAX+Mesh%KP1), &
                         INTENT(OUT) :: shifted_field
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,local_joff
    !------------------------------------------------------------------------!
    local_joff = this%local_joff

    IF (Mesh%shear_dir.EQ.1) THEN
!NEC$ IVDEP
      DO k = Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
!NEC$ IVDEP
        DO j = Mesh%JMIN,Mesh%JMAX
          this%joff(j)   = Mesh%Q*Mesh%OMEGA*Mesh%bcenter(Mesh%IMIN,j,k,2)* &
                      delt/Mesh%dx
          this%jrem(j)  = this%joff(j) - FLOOR(this%joff(j))
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX
            shifted_field(k)%original(i,j-local_joff) = &
             (1.0-this%jrem(j))*field(1+MODULO(i-1+FLOOR(this%joff(j)), &
              Mesh%IMAX-Mesh%IMIN+1),j,k) + this%jrem(j)*field(1+MODULO(i+ &
              FLOOR(this%joff(j)),Mesh%IMAX-Mesh%IMIN+1),j,k)
          END DO
        END DO
      END DO
    ELSE ! must be Mesh%shear_dir.EQ.2, because otherwise initialization would raise an error
!NEC$ IVDEP
      DO k = Mesh%KMIN-Mesh%KP1,Mesh%KMAX+Mesh%KP1
!NEC$ IVDEP
        DO i = Mesh%IMIN,Mesh%IMAX
          this%joff(i)   = -Mesh%Q*Mesh%OMEGA*Mesh%bcenter(i,Mesh%JMIN,k,1)* &
                      delt/Mesh%dy
          this%jrem(i)  = this%joff(i) - FLOOR(this%joff(i))
        END DO
        DO j = Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i = Mesh%IMIN,Mesh%IMAX
            shifted_field(k)%original(i,j-local_joff) = &
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
   TYPE(gravity_sboxspectral), INTENT(INOUT) :: this
   !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! Destroy plans
    CALL fftw_destroy_plan(this%plan_r2c)
    CALL fftw_destroy_plan(this%plan_c2r)
#if defined(PARALLEL)
    CALL fftw_free(this%fftw_real_pointer)
    CALL fftw_free(this%fftw_complex_pointer)
    DEALLOCATE(this%mpi_buf)
#else
    DEALLOCATE(this%Fmass3D,this%Fsum3D,this%field)
#endif
    ! Free all temporary memory of FFTW
    CALL fftw_cleanup()

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
    CALL this%Finalize_base()
   END SUBROUTINE Finalize

END MODULE gravity_sboxspectral_mod
