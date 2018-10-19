!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_spectral.f03                                              #
!#                                                                           #
!# Copyright (C) 2011-2018                                                   #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!! - parameters of \link gravity_spectral \endlink as key-values
!! \key{green,INTEGER, type of Green-function,1}
!! \key{sigma,REAL,standard deviation,0.05}
!! \key{output/potential,INTEGER,enable(=1) output of grav. potential}
!----------------------------------------------------------------------------!
!> \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief poisson solver via spectral methods and direct integration
!!
!! \attention This is a 2D solver and works only in flat, polar geometries. This
!!            is why it should be paid attention to the dimension of the
!!            arrays. They where only made allocated 3D when necessary.
!!
!! \cite chan2006
!! \cite li2009
!----------------------------------------------------------------------------!
MODULE gravity_spectral_mod
  USE gravity_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE logging_base_mod
  USE functions
  USE common_dict
  USE fftw
#ifdef PARALLEL
  USE mpi
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: solver_name  = "spectral"
  REAL, PARAMETER              :: SQRTTWOPI &
    = 2.50662827463100050241576528481104525300698674

  TYPE, EXTENDS(gravity_base) :: gravity_spectral
#if defined(HAVE_FFTW)
    !> \name
    !!#### spectral poisson solver
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
    REAL, DIMENSION(:,:), POINTER    :: tmp2D
    REAL, DIMENSION(:,:), POINTER    :: phi2D
    REAL, DIMENSION(:), POINTER      :: height1D
    INTEGER                          :: green
    REAL                             :: sigma, ecut
    INTEGER, DIMENSION(:), POINTER   :: sizes
    INTEGER, POINTER                 :: mcut
    INTEGER                          :: INUM,KNUM
    INTEGER                          :: MNUM          !< number of modes
#if defined(PARALLEL)
    !> displacment and length of domain
    INTEGER,DIMENSION(:),POINTER     :: displ, num
    INTEGER                          :: mpi_error
#endif
#endif
  CONTAINS
    PROCEDURE :: InitGravity_spectral
    PROCEDURE :: UpdateGravity_single
    PROCEDURE :: InfoGravity
    PROCEDURE :: CalcDiskHeight_single
    PROCEDURE :: Finalize
#ifdef HAVE_FFTW
    PROCEDURE :: CalcPotential
    PROCEDURE :: PrecomputeI
    PROCEDURE :: CalcMcut
#endif
  END TYPE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       gravity_spectral
  !--------------------------------------------------------------------------!
  CONTAINS

  SUBROUTINE InitGravity_spectral(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_spectral), INTENT(INOUT) :: this
    CLASS(mesh_base),            INTENT(IN) :: Mesh
    CLASS(physics_base),         INTENT(IN) :: Physics
    TYPE(Dict_TYP),              POINTER    :: config,IO
    !------------------------------------------------------------------------!
    INTEGER           :: err, valwrite, i, solver
    CHARACTER(LEN=32) :: info_str
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", solver)
    CALL this%InitLogging(solver,solver_name)
#if !defined(HAVE_FFTW)
    CALL this%Error("InitGravity_spectral", &
         "Mandatory requirement fft has not been enabled. "//&
         "Please add --with-fftw=$FFTWDIR or similar to configure call.")
#endif
#if defined(HAVE_FFTW)
#ifdef PARALLEL
    ! Check domains a annular rings
    IF(.NOT.(Mesh%dims(2).EQ.1)) &
      CALL this%Error("InitGravity_spectral", &
                 "Only domains shaped as annular rings are allowed (N x 1 decompositions).")
#endif

!    IF(.NOT.(GetType(Boundary(NORTH)).EQ.PERIODIC .AND. &
!             GetType(Boundary(SOUTH)).EQ.PERIODIC)) THEN
!      CALL Error(this,"InitGravity_spectral", &
!                 "The boundary conditions in north and south direction have " // &
!                 "to be periodic! Don't forget: This kind of " // &
!                 "self-gravitation only works for polar-like coordinate " // &
!                 "systems.")
!    END IF


    IF(.NOT.(MOD(Mesh%JNUM,2)==0)) THEN
      CALL this%Error("InitGravity_spectral", &
                 "The spectral poisson solver needs an even number of cells " // &
                 "in the phi direction due to the discrete cosinus transform.")
    END IF


! If this is the innermost domain include 1 ghost cell,
! if this is the outermost domain include also 1 ghost cell.
#ifdef PARALLEL
    IF(Mesh%IMAX.EQ.Mesh%INUM) THEN
      this%IMAX = Mesh%IMAX+1
    ELSE
      this%IMAX = Mesh%IMAX
    END IF
#else
    !always innermost and outermost domain
    this%IMAX = Mesh%IMAX+1
#endif

    this%INUM = this%IMAX - Mesh%IMIN + 1

    ! number of fourier modes
    this%MNUM = Mesh%JNUM/2+1

    ALLOCATE(this%phi2D(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DIM), &
             this%tmp2D(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%pot(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4), &
             this%height1D(1:Mesh%INUM), &
             this%sizes(1),&
             this%mcut, &
#ifdef PARALLEL
             this%sbuf1(Mesh%GNUM,Mesh%JNUM), &
             this%rbuf1(Mesh%GNUM,Mesh%JNUM), &
             this%sbuf2(Mesh%GNUM,Mesh%JNUM), &
             this%rbuf2(Mesh%GNUM,Mesh%JNUM), &
             this%displ(0:this%GetNumProcs()-1),&
             this%num(0:this%GetNumProcs()-1),&
#endif
             STAT=err)

    CALL GetAttr(config, "green", this%green, 1)

    CALL GetAttr(config, "sigma", this%sigma, 0.05)
    this%height1D(:) = this%sigma

    ! mode cut off.
    !  0: disabled
    ! >0: constant cutoff number
    !     (can be overwritten by ecut>0.)
    CALL GetAttr(config, "mcut", this%mcut, this%MNUM)
    this%mcut = MIN(this%mcut,this%MNUM)

    ! mode cutoff relative error boundary
    ! <=0: no automatic cutoff calculation
    ! >0 : automatic setting of mcut
    CALL GetAttr(config, "ecut", this%ecut, 0.)

    CALL this%Info(" POISSON--> spectral method")
    CALL this%Info(" POISSON--> Initializing")

    WRITE (info_str, '(I8)') this%green
    CALL this%Info(" POISSON--> green-fn type:     " // TRIM(info_str))
    WRITE (info_str, '(ES8.2)') this%sigma
    CALL this%Info(" POISSON--> sigma:             " // TRIM(info_str))


    !\todo only 2D variables
    this%p_FI = fftw_alloc_real(INT(2*this%MNUM * (this%INUM)*(Mesh%INUM), C_SIZE_T))
    CALL C_F_POINTER(this%p_FI, this%FI, &
                     [2*this%MNUM,this%INUM,Mesh%INUM])
    CALL C_F_POINTER(this%p_FI, this%cFI, &
                     [this%MNUM,this%INUM,Mesh%INUM])


    this%pFdensity = fftw_alloc_real(INT(2*this%MNUM*Mesh%INUM, C_SIZE_T))
    CALL C_F_POINTER(this%pFdensity, this%Fdensity, &
                     [2*this%MNUM,Mesh%INUM])
    CALL C_F_POINTER(this%pFdensity, this%cFdensity, &
                     [Mesh%JNUM/2+1,Mesh%INUM])

    this%pFphi = fftw_alloc_real(INT(2*this%MNUM*this%INUM,C_SIZE_T))
    CALL C_F_POINTER(this%pFphi, this%Fphi, &
                     [2*this%MNUM, this%INUM])
    CALL C_F_POINTER(this%pFphi, this%cFphi, &
                     [Mesh%JNUM/2+1, this%INUM])

    IF (err.NE.0) &
        CALL this%Error("InitGravity_spectral","Memory allocation failed.")

    this%block => this%Fdensity(:,Mesh%IMIN:Mesh%IMAX)
    this%cblock => this%cFdensity(:,Mesh%IMIN:Mesh%IMAX)

    this%phi2D(:,:) = 0.
    this%accel(:,:,:,:) = 0.
    this%pot(:,:,:,:) = 0.
    this%tmp2D(:,:) = 0.

    ! Create plans for fftw

    ! Use FFTW_MEASURE for calculating the fastest plan, but this
    ! costs some extra seconds
    this%sizes(1) = Mesh%JNUM
    this%plan_r2c = fftw_plan_many_dft_r2c(1, this%sizes, (Mesh%IMAX-Mesh%IMIN+1),&
                                           this%block, this%sizes, &
                                           1, 2*this%MNUM, &
                                           this%cblock, this%sizes, &
                                           1, this%MNUM, &
                                           FFTW_MEASURE)

    this%sizes(1) = Mesh%JNUM
    this%plan_c2r = fftw_plan_many_dft_c2r(1, this%sizes, this%INUM, &
                                           this%cFphi, this%sizes, &
                                           1, this%MNUM, &
                                           this%Fphi, this%sizes, &
                                           1, 2*this%MNUM, &
                                           FFTW_MEASURE)

    CALL GetAttr(config, "output/potential", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "potential",this%phi2D(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))

    IF(this%ecut.GT.0.) &
      CALL SetAttr(IO, "mcut",Ref(this%mcut))

#ifdef PARALLEL
    this%displ((this%GetRank())) = (Mesh%IMIN-1)*2*this%mcut
    this%num((this%GetRank())) = (Mesh%IMAX-Mesh%IMIN+1)*2*this%mcut
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       this%displ, 1, MPI_INTEGER, MPI_COMM_WORLD, this%mpi_error)
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       this%num, 1, MPI_INTEGER, MPI_COMM_WORLD, this%mpi_error)
#endif

    CALL this%PrecomputeI(Mesh, Physics)

    CALL this%Info(" POISSON--> .. done initializing")

    CALL this%InitGravity(Mesh,Physics,config,IO)
#endif

  END SUBROUTINE InitGravity_spectral


#ifdef HAVE_FFTW
  FUNCTION CalcMcut(this,Mesh) RESULT(res)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_spectral)      :: this
    CLASS(mesh_base), INTENT(IN) :: Mesh
    INTEGER                      :: res
    !------------------------------------------------------------------------!
    INTEGER           :: m,i,ier
    REAL              :: cumsum(this%MNUM)
    INTEGER           :: mcut(Mesh%IMIN:Mesh%IMAX)
    !------------------------------------------------------------------------!
    DO i=Mesh%IMIN,Mesh%IMAX
      cumsum(1) = this%Fdensity(2*1-1,i)**2 + this%Fdensity(2*1,i)**2
      DO m=2,this%MNUM
        cumsum(m) = cumsum(m-1) + this%Fdensity(2*m-1,i)**2 + this%Fdensity(2*m,i)**2
      END DO
      DO m=1,this%MNUM
        IF(ABS(cumsum(m)/cumsum(this%MNUM)-1.).LT.this%ecut) THEN
          mcut(i) = m
          EXIT
        END IF
      END DO
    END DO
    res = MAXVAL(mcut)

#ifdef PARALLEL
    CALL MPI_Allreduce(MPI_IN_PLACE,res,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
#endif
  END FUNCTION CalcMcut
#endif

  !> Prints out information
  SUBROUTINE InfoGravity(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_spectral), INTENT(IN) :: this
    CLASS(mesh_base),            INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
  END SUBROUTINE InfoGravity


#if defined(HAVE_FFTW)
  SUBROUTINE CalcPotential(this,Mesh,Physics,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_spectral), INTENT(INOUT):: this
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                         INTENT(IN) :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k,m,i0,ier,oldmcut
#ifdef PARALLEL
    INTEGER           :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
    ! Fourier transform the density with respect to Phi
    DO k=Mesh%KMIN, Mesh%KMAX
      DO i=Mesh%IMIN, Mesh%IMAX
        this%Fdensity(1:Mesh%JNUM,i) = pvar(i,Mesh%JMIN:Mesh%JMAX,k,Physics%DENSITY)
      END DO
    END DO

    CALL fftw_execute_dft_r2c(this%plan_r2c, this%block, this%cblock)

    IF(this%ecut.GT.0.) THEN
      oldmcut = this%mcut

      this%mcut = this%CalcMcut(Mesh)

#ifdef PARALLEL
      this%displ = this%displ / oldmcut * this%mcut
      this%num = this%num / oldmcut * this%mcut
#endif
    END IF

#ifdef PARALLEL
    ! Distribute Fdensity to all processes
      CALL MPI_AllGatherV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&
                          this%Fdensity(1:2*this%mcut,1:Mesh%INUM), &
                          this%num, this%displ, DEFAULT_MPI_REAL, &
                          MPI_COMM_WORLD, this%mpi_error)
#endif

    ! Integrate in radial direction by numerical quadrature
    ! Fphi_m(t,r) = int_rmin^rmax Fdensity_m(t,r') * I_m(r,r') dr'
    ! dr' is already multiplicated into FI.

    ! Direct summation and multiplication with complex numbers:
    !DO m=1, this%MNUM
    !  ! integrate
    !  DO i=1, this%INUM
    !    this%cFPhi(m,i) = SUM(this%cFdensity(m,:) * this%cFI(m,i,:))
    !  END DO
    !END DO

    ! or with real numbers, real and imaginary part seperate:
    this%FPhi(:,:) = 0.
    DO i0=1,Mesh%INUM
      DO i=1,this%INUM
        DO m=1,this%mcut
           ! real part
           this%FPhi(2*m-1,i) = this%FPhi(2*m-1,i) &
             + this%Fdensity(2*m-1,i0) * this%FI(2*m-1,i,i0)
           ! imag part
           this%FPhi(2*m,i) = this%FPhi(2*m,i) &
             + this%Fdensity(2*m,i0) * this%FI(2*m-1,i,i0)
        END DO
      END DO
    END DO

    ! Inverse fourier transform Fphi with respect to Phi
    CALL fftw_execute_dft_c2r(this%plan_c2r, this%cFphi, this%Fphi)

    DO i=1,this%INUM
      this%phi2D(i+Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX) = this%FPhi(1:Mesh%JNUM,i)
    END DO

    ! calculate boundary data (only one ghost cell is enough)
#ifdef PARALLEL
    ! send boundary data to western and receive from eastern neighbor
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         this%sbuf1 = this%phi2D(Mesh%IMIN:Mesh%IMIN+Mesh%GNUM-1,Mesh%JMIN:Mesh%JMAX)
    CALL MPI_Sendrecv(this%sbuf1,Mesh%GNUM*Mesh%JNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,this%rbuf1, &
         Mesh%GNUM*Mesh%JNUM,DEFAULT_MPI_REAL,Mesh%neighbor(EAST), &
         MPI_ANY_TAG,Mesh%comm_cart,status,this%mpi_error)
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         this%phi2D(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX) = this%rbuf1
    ! send boundary data to western and receive from eastern neighbor
    IF (Mesh%neighbor(EAST).NE.MPI_PROC_NULL) &
         this%sbuf2 = this%phi2D(Mesh%IMAX-Mesh%GNUM+1:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
    CALL MPI_Sendrecv(this%sbuf2,Mesh%GNUM*Mesh%JNUM, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,this%rbuf2, &
         Mesh%GNUM*Mesh%JNUM,DEFAULT_MPI_REAL,Mesh%neighbor(WEST), &
         MPI_ANY_TAG,Mesh%comm_cart,status,this%mpi_error)
    IF (Mesh%neighbor(WEST).NE.MPI_PROC_NULL) &
         this%phi2D(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX) = this%rbuf2
#endif
    this%phi2D(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JMIN-1) &
      = this%phi2D(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX)
    this%phi2D(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMAX+1:Mesh%JGMAX) &
      = this%phi2D(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN:Mesh%JMIN+Mesh%GNUM-1)

  END SUBROUTINE CalcPotential
#endif

  !> Green's function
  !!
  !! This is integral, which has to be calculated for a general \f$ Z(r,z) \f$
  !! \f[
  !!    G = \int_{-\infty}^{\infty} -\frac{Z(r_1,z_1)}{\sqrt{r^2+r'^2 - 2\,r\,r'\cos{(\phi)} +
  !!    \epsilon^2 + z'^2} dz_1.
  !! \f]
  !! For a vertical Gaussian density distribution \f$ Z(r,z) \f$
  !! \f[
  !!    Z(r,z) = (2\pi(H(r))^2)^{-0.5}\exp{(-\frac{z^2}{2(H(r))^2})}
  !! \f]
  !! this integral can be evaluated analyticly, which results in
  !! \f[
  !!    G(r,r',\phi) = - \frac{\exp{(R^2/4)} * K_0(R^2/4)}{\sqrt{2*\pi} * H(r')},
  !! \f]
  !! with \f$ R^2 = r^2 + r'^2 - 2\,r\,r'\cos{(\phi)} + \epsilon^2)/(H(r'))^2 \f$
  !! and \f$ K_0 \f$ the modified Bessel function of the second kind.
#ifdef HAVE_FFTW
  ELEMENTAL FUNCTION GreenFunction(dr2, green, sigma) RESULT(G)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL,    INTENT(IN) :: dr2
    INTEGER, INTENT(IN) :: green
    REAL,    INTENT(IN) :: sigma
    REAL                :: G
    !------------------------------------------------------------------------!
    SELECT CASE(green)
      CASE(1)  ! Razor sharp disc
        G = -1.0/SQRT(dr2)
      CASE(2)  ! gaussian spheres
        G = -1.0 *  Bessel_K0e(dr2/(4.0*sigma*sigma)) &
            / (SQRTTWOPI*sigma)
      CASE(3)  ! Used for the orbiting cylinder example
        G = LOG(dr2)
      CASE DEFAULT
        ! should never happen
        G = 0.0
    END SELECT
  END FUNCTION GreenFunction
#endif

  !> Precomputes the fourier transform
  !!
  !! Precompute the fourier transform with respect to \f$ \Phi-\Phi' \f$ of
  !! \f[
  !!    I(r,r',\Phi-\Phi') = 2 * \pi * r' * G(r,r',\Phi-\Phi'),
  !! \f}
  !! where G is the softened Green's function.
  !!
  !! In the case of a gaussian density distribution in the z direction,
  !! we have:
  !! \f[
  !!    Z(r,z) = (2*\pi*H(r))^{0.5} * \exp{( -z^2 / (2*(H(r))^2) )},
  !! \f]
  !! and therefore as Green's function:
  !! \f[
  !!    G(r,r',Phi-Phi') = -\left(exp(R^2/4)*K_0(R^2/4))/(\sqrt{2*pi}*H(r')\right),
  !! \f]
  !! with \f$ R^2 = \left(r^2 + r'^2 - 2\,r\,r'\,\cos{(\Phi-\Phi')} + \epsilon^2) /
  !! (H(r')\right)^2 \f$
  !!
  !! (epsilon is a small softening parameter)
#ifdef HAVE_FFTW
  SUBROUTINE PrecomputeI(this, Mesh, Physics)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   CLASS(gravity_spectral)             :: this
   CLASS(mesh_base),        INTENT(IN) :: Mesh
   CLASS(physics_base)                 :: Physics
   !------------------------------------------------------------------------!
   TYPE(C_PTR)        :: plan_r2c
   INTEGER            :: i0, i1, i, j, k, m
   REAL               :: r0, notzero
   CHARACTER(LEN=128) :: str
   REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX) &
                      :: phi
   REAL, DIMENSION(1:Mesh%INUM+1) :: r,hx,hy
   REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX,Mesh%IMIN:this%IMAX) :: dr2
   TYPE(C_PTR)        :: p
   COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), POINTER &
                      :: cFI
   INTEGER            :: rank, howmany, istride, idist, ostride, odist
   INTEGER, DIMENSION(1) &
                      :: n
#ifdef PARALLEL
   INTEGER, DIMENSION(0:(this%GetNumProcs())-1) :: num,displ
#endif
   !------------------------------------------------------------------------!
    n(1)    = Mesh%JNUM
    rank    = 1
    howmany = this%INUM * Mesh%INUM
    istride = 1
    idist = 2*this%MNUM
    ostride = istride
    odist   = this%MNUM
    plan_r2c = fftw_plan_many_dft_r2c(rank, n, howmany, &
                                  this%FI, n, &
                                  istride, idist, &
                                  this%cFI, n, &
                                  ostride, odist, &
                                  FFTW_MEASURE &
                                  )

    phi = Mesh%curv%center(Mesh%IMIN,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN,2)&
        - Mesh%curv%center(Mesh%IMIN,Mesh%JMIN,Mesh%KMIN,2)

    r(Mesh%IMIN:this%IMAX) = Mesh%radius%bcenter(Mesh%IMIN:this%IMAX,Mesh%JMIN,Mesh%KMIN)
    hx(Mesh%IMIN:this%IMAX) = Mesh%hx%bcenter(Mesh%IMIN:this%IMAX,Mesh%JMIN,Mesh%KMIN)
    hy(Mesh%IMIN:this%IMAX) = Mesh%hy%bcenter(Mesh%IMIN:this%IMAX,Mesh%JMIN,Mesh%KMIN)
#ifdef PARALLEL
    displ(this%GetRank()) = (Mesh%IMIN-1)
    num(this%GetRank()) = this%INUM
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       displ, 1, MPI_INTEGER, MPI_COMM_WORLD, this%mpi_error)
    CALL MPI_AllGather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       num, 1, MPI_INTEGER, MPI_COMM_WORLD, this%mpi_error)
    CALL MPI_AllGatherV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        r, num, displ, &
                        DEFAULT_MPI_REAL, MPI_COMM_WORLD, this%mpi_error)
    CALL MPI_AllGatherV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        hx, num, displ, &
                        DEFAULT_MPI_REAL, MPI_COMM_WORLD, this%mpi_error)
    CALL MPI_AllGatherV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        hy, num, displ, &
                        DEFAULT_MPI_REAL, MPI_COMM_WORLD, this%mpi_error)
#endif

    DO k=Mesh%KMIN,Mesh%KMAX
      DO i1=1,Mesh%INUM
        DO i=Mesh%IMIN,this%IMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            r0 = Mesh%radius%faces(i,j,k,1)
            dr2(j,i) = r(i1)**2 + r0**2 - 2.*r(i1)*r0*COS(phi(j))
          END DO
        END DO
!CDIR IEXPAND
        this%FI(1:Mesh%JNUM,1:this%INUM,i1) &
          = 2.0 * PI * hx(i1) * hy(i1) * Mesh%dx &
          * Physics%Constants%GN &
          * GreenFunction(dr2(:,:),this%green,this%height1D(i1)) &
          / (Mesh%JNUM)**2
      END DO
    END DO

    CALL fftw_execute_dft_r2c(plan_r2c, this%FI, this%cFI)

    ! Free plan_r2c
    CALL fftw_destroy_plan(plan_r2c)

    ! I has an even symmetry around 0. Normally a cosinus tranform (r2r)
    ! would be sufficent. But this is not available on all platforms. There-
    ! fore do a standard r2c dft and ditch the imaginary part of the result,
    ! which is approx zero numerically (and exactly zero analytically).
    !this%cFI = REAL(this%cFI)
    notzero = SUM(ABS(this%FI(2:2*this%MNUM:2,:,:)))/(this%MNUM*Mesh%JNUM*this%INUM)
#ifdef PARALLEL
    CALL MPI_Allreduce(MPI_IN_PLACE, notzero, 1, DEFAULT_MPI_REAL, &
      MPI_MAX, &
      MPI_COMM_WORLD, this%mpi_error)
#endif
    IF(notzero.GT.1.) THEN
      WRITE(str,'(A,ES12.4,A)')"Imag(FI) should be zero, but ",notzero,&
        " > 1."
      CALL this%Warning("ProcomputeI_spectral",TRIM(str))
    END IF

    this%FI(2:2*this%MNUM:2,:,:) = 0.
  END SUBROUTINE PrecomputeI
#endif

  !> \attention This routine only works in 2D
  SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_spectral), INTENT(INOUT) :: this
    CLASS(mesh_base),            INTENT(IN) :: Mesh
    CLASS(physics_base),         INTENT(IN) :: Physics
    CLASS(fluxes_base),          INTENT(IN) :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                                 INTENT(IN) :: pvar
    REAL,                        INTENT(IN) :: time,dt
    !------------------------------------------------------------------------!
    INTEGER              :: i, j, k
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! calc potential first
    CALL this%CalcPotential(Mesh,Physics,pvar)

    ! compute gravitational acceleration using the gradient of the potential
    ! ATTENTION: ghost cell values and any other component of the gravitational
    !            acceleration are set to zero in InitGravity_spectral
    ! \todo This difference quotient is still 2D and only the third component is
    !       added. It will not work in 3D.
    DO k = Mesh%KMIN,Mesh%KMAX
      DO j = Mesh%JMIN-1,Mesh%JMAX+1
        DO i = Mesh%IMIN-1,Mesh%IMAX
          this%accel(i,j,k,1) = -1.0*(this%phi2D(i+1,j)-this%phi2D(i,j))/Mesh%dlx(i,j,k)
          this%accel(i,j,k,2) = -1.0*(this%phi2D(i+1,j+1)+this%phi2D(i,j+1) &
                           -this%phi2D(i+1,j-1)-this%phi2D(i,j-1))&
                            /(4.0*Mesh%dly(i,j,k))
          this%pot(i,j,k,1) = 0.5*(this%phi2D(i,j)+this%phi2D(i+1,j))
          this%pot(i,j,k,2) = this%phi2D(i+1,j)
          this%pot(i,j,k,3) = 0.25*(this%phi2D(i,j)+this%phi2D(i+1,j)+this%phi2D(i,j+1)+this%phi2D(i+1,j+1))
        END DO
      END DO
    END DO
#endif
  END SUBROUTINE UpdateGravity_single


  !> \public compute disk pressure scale height for geometrically thin self-gravitating disks
  !!
  !! The algorithm solves the equations
  !! \f{eqnarray}{
  !!     1/h^2 &=& 1/h_{sg}^2 + 1/h_{ext}^2 & (1)\\
  !!     \left.\frac{\partial^2\Phi_{sg}}{\partial z^2}\right|_{z=0} = (c_s/h_{sg})^2
  !!           &=& 4\pi G \rho_c - \left.\Delta_{\xi\eta} \Phi_{sg}\right|_{z=0} & (2) \\
  !!     \rho_c &=& 1/\sqrt{2\pi} \Sigma / h & (3)
  !! \f}
  !! for the disk scale height \f$ h \f$ where \f$ h_{sg} \f$ is the self-gravitating
  !! pressure scale height due to the material in the disk and \f$ h_{ext} \f$ the scale
  !! height caused by an external (probably point mass) potential.
  !! \f$ \Sigma \f$ is the surface density, \f$ c_s \f$ is the speed of sound and
  !! \f$ \rho_c \f$ is the density in the equatorial plane.
  !! \f$ \Delta_{\xi\eta} \f$ denotes the 2D Laplacian in curvilinear coordinates
  !! and \f$ \Phi_{sg} \f$ is the gravitational potential of the disk.
  !! One can rewrite the last term in equation (2) replacing the Laplacian of the
  !! potential with the divergence of the 2D gravitational acceleration according to
  !! \f[
  !!     \left.\Delta_{\xi\eta} \Phi_{sg}\right|_{z=0}
  !!  = -\nabla_{\xi\eta}\cdot\vec{g}_{sg}.
  !! \f]
  !! Replacing \f$ \rho_c \f$ in (2) by (3) and inserting the result in (1) yields
  !! a quadratic equation in \f$ h_{ext}/h \f$ with the two solutions
  !! \f[
  !!     h_{ext}/h = p \pm \sqrt{q+p^2}
  !! \f]
  !! where
  !! \f[
  !!     p = \sqrt{2\pi} G \Sigma h_{ext}/c_s^2
  !! \f]
  !! and
  !! \f[
  !!     q = 1 + \left(h_{ext}/c_s\right)^2 + \nabla_{\xi\eta}\cdot\vec{g}_{sg}.
  !! \f]
  !! The sign of the root must be positive, because in the non self-gravitating limit
  !! \f$ p=0, q=1 \f$ and \f$ h_{ext}/h \f$ must become \f$ 1 \f$ not \f$ -1 \f$ .
  !!
  !! \attention This routine only works in 2D
  SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_spectral), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX,Physics%VNUM), &
                         INTENT(IN)    :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX), &
                         INTENT(IN)    :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KMIN:Mesh%KMAX), &
                         INTENT(INOUT) :: h_ext,height
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: phi
    INTEGER           :: i,j,k
    REAL              :: cs2,p,q
    !------------------------------------------------------------------------!
#ifdef HAVE_FFTW
    ! compute the laplacian of the gravitational potential in the disks
    ! equatorial plane; gravity_generic updates the gravitational acceleration
    DO k=Mesh%KMIN-1,Mesh%KMAX
      DO j=Mesh%JMIN-1,Mesh%JMAX
        DO i=Mesh%IMIN-1,Mesh%IMAX
          phi(i,j,1) = Mesh%hy%faces(i,j,k,2)*Mesh%hz%faces(i,j,k,2)/Mesh%hx%faces(i,j,k,2)&
           * this%pot(i,j,k,2)
          phi(i,j,2) = Mesh%hx%faces(i,j,k,4)*Mesh%hz%faces(i,j,k,4)/Mesh%hy%faces(i,j,k,4)&
           * this%pot(i,j,k,3)
        END DO
      END DO
    END DO
    ! div(-grad(phi)))
    DO k=Mesh%KMIN,Mesh%KMAX
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          this%tmp2D(i,j) = - 1./Mesh%sqrtg%bcenter(i,j,k) &
            * (  (phi(i,j,1) - phi(i-1,j,1)) / Mesh%dlx(i,j,k) &
               + (phi(i,j,2) - phi(i,j-1,2)) / Mesh%dly(i,j,k))
        END DO
      END DO
    END DO

!CDIR COLLAPSE
    ! disk height with external potential
!CDIR COLLAPSE
!    DO k=Mesh%KGMIN,Mesh%KGMAX
!      DO j=Mesh%JGMIN,Mesh%JGMAX
!        DO i=Mesh%IGMIN,Mesh%IGMAX
!          cs2 = bccsound(i,j,k)**2
!          p = SQRTTWOPI*Physics%Constants%GN*pvar(i,j,k,Physics%DENSITY)*h_ext(i,j,k)/cs2
!          q = 1. + this%tmp2D(i,j)*h_ext(i,j,k)**2/cs2
!          ! return the new disk height
!          height(i,j,k) = h_ext(i,j,k) / (p+SQRT(q+p*p))
!        END DO
!      END DO
!    END DO
       ! pure self-gravitating disk without external potential
!CDIR COLLAPSE
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            cs2 = bccsound(i,j,k)**2
            p = SQRTTWOPI*Physics%Constants%GN*pvar(i,j,k,Physics%DENSITY)*Mesh%radius%bcenter(i,j,k)/cs2
            q = this%tmp2D(i,j)*h_ext(i,j,k)**2/cs2
            ! return the new disk height
            height(i,j,k) = Mesh%radius%bcenter(i,j,k) / (p+SQRT(q+p*p))
          END DO
        END DO
      END DO
#endif
  END SUBROUTINE CalcDiskHeight_single


  SUBROUTINE Finalize(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   CLASS(gravity_spectral), INTENT(INOUT) :: this
   !------------------------------------------------------------------------!
#if defined(HAVE_FFTW)
    ! Destroy plans
    CALL fftw_destroy_plan(this%plan_r2c)
    CALL fftw_destroy_plan(this%plan_c2r)

    ! Free memomry
    DEALLOCATE(&
#ifdef PARALLEL
               this%sbuf1, &
               this%rbuf1, &
               this%sbuf2, &
               this%rbuf2, &
               this%displ,this%num,&
#endif
               this%phi2D, &
               this%pot, &
               this%accel, &
               this%tmp2D, &
               this%sizes, &
               this%mcut &
               )

    CALL fftw_free(this%p_FI)
    CALL fftw_free(this%pFdensity)
    CALL fftw_free(this%pFphi)
#endif
    END SUBROUTINE Finalize

END MODULE gravity_spectral_mod
