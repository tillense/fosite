!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_sboxspectral.f90                                          #
!#                                                                           #
!# Copyright (C) 2015                                                        #
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
!! References:
!!
!! \cite gammie2001 , \cite gammiecode
!----------------------------------------------------------------------------!
MODULE gravity_sboxspectral
  USE gravity_common
  USE mesh_generic
  USE physics_generic
  USE boundary_generic
  USE functions
  USE common_dict
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY)
  USE fftw
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER &
                    :: solver_name  = "sboxspectral"
  REAL, PARAMETER              :: SQRTTWOPI &
    = 2.50662827463100050241576528481104525300698674
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       Grid_TYP, &
       ! constants
       ! methods
       InitGravity_sboxspectral, &
       UpdateGravity_sboxspectral, &
       CalcPotential_sboxspectral, &
       CalcDiskHeight_sboxspectral, &
       CloseGravity_sboxspectral, &
       FieldShift
  !--------------------------------------------------------------------------!
  CONTAINS
  !> \public Constructor of gravity sboxspectral module
  !!
  !! This subroutine reads the necessary config data for the the specral
  !! gravity solver within the shearingbox. It initializes the gravity type
  !! and various mesh data arrays. Some of those are marked for output.
  !! There are three different spectral solvers used:
  !! 1. FFTW (normally used)
  !! 2. FFTW-legacy (if old compilers are used)
  !! 3. FFT-Keisan (optimized for NEC vector machines)
  SUBROUTINE InitGravity_sboxspectral(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP),POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Dict_TYP), POINTER :: config,IO
    INTEGER           :: solver
    !------------------------------------------------------------------------!
    INTEGER           :: err, valwrite
    CHARACTER(LEN=32) :: info_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    !>\todo{implement a boundary check (only periodic or sbox boundaries)}
    !>\todo{implement a check for equidistant mesh}

    CALL GetAttr(config, "gtype", solver)
    CALL InitGravity(this,solver,solver_name)

    !-------------- checks & warnings for initial conditions ----------------!
#if !defined(HAVE_FFTW) && !defined(HAVE_FFTW_LEGACY) && !defined(HAVE_FFTKEISAN)
    CALL Error(this,"InitGravity_sboxspectral", &
         "No fft package could be loaded. "// &
         "Try --with-fftw=$FFTWDIR or --with-fftkeisan=$FFTKEISANDIR or "// &
         "similar to configure call.")
#endif
#if !defined(HAVE_ISO_C_BINDING) && defined(HAVE_FFTW)
    CALL Error(this,"InitGravity_sboxspectral", &
         "No ISO_C_BINDINGs are available for this compiler, but they are "// &
         "a Mandatory requirement for fftw with modern fortran.")
#endif

#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY) || defined(HAVE_FFTKEISAN)
    ! Check even number of cells
    IF(.NOT.(MOD(Mesh%JMAX-Mesh%JMIN+1,2)==0)) THEN
      CALL Error(this,"InitGravity_sboxspectral", &
                 "The spectral poisson solver needs an even number of "// &
                 "cells in the phi direction due to the discrete cosinus "// &
                 "transform.")
    END IF
    !------------------------------------------------------------------------!

    ALLOCATE( &
             this%phi(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
             this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY)
             this%mass2D(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &     ! input density
#elif defined(HAVE_FFTKEISAN)
             this%mass2D(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX), &
#endif
             this%Fmass2D(Mesh%IMIN:Mesh%IMAX/2+1,Mesh%JMIN:Mesh%JMAX), &
             this%Fmass2D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX), &
             this%joff(Mesh%IMIN:Mesh%IMAX), &
             this%jrem(Mesh%IMIN:Mesh%IMAX), &
             this%den_ip(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX), &
             STAT=err)
    IF (err.NE.0) &
        CALL Error(this,"InitGravity_sboxspectral","Memory allocation failed.")

    !------------------- create plans for fftw ------------------------------!
    ! Pay attention to the argument order of the dimension (JNUM and INUM    !
    ! are switched because of C -> row-major, Fortran -> column-major),      !
    ! BUT ONLY in modern Fortran UNLIKE the legacy version                   !
    ! ------------ plans are allocated in dictionary ------------------------!
#if defined(HAVE_FFTW)
    this%plan_r2c = fftw_plan_dft_r2c_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1,this%mass2D, &
                                         this%Fmass2D,FFTW_MEASURE)
    this%plan_c2r = fftw_plan_dft_c2r_2d(Mesh%JMAX-Mesh%JMIN+1, &
                                         Mesh%IMAX-Mesh%IMIN+1, this%Fmass2D, &
                                         this%mass2D, FFTW_MEASURE)
#elif defined(HAVE_FFTW_LEGACY)
    CALL dfftw_plan_dft_r2c_2d(this%plan_r2c,Mesh%IMAX-Mesh%IMIN+1, &
                                         Mesh%JMAX-Mesh%JMIN+1, this%mass2D, &
                                         this%Fmass2D, FFTW_MEASURE)
    CALL dfftw_plan_dft_c2r_2d(this%plan_c2r,Mesh%IMAX-Mesh%IMIN+1,
                                         Mesh%JMAX-Mesh%JMIN+1, this%Fmass2D, &
                                         this%mass2D, FFTW_MEASURE)
#endif
    !------------------------------------------------------------------------!
    ! set potential and acceleration to zero
    this%phi(:,:) = 0.
    this%accel(:,:,:) = 0.
    this%mass2D(:,:) = 0.
    this%Fmass2D(:,:) = 0.
    this%Fmass2D_real(:,:) = 0.
    this%joff(:) = 0.
    this%jrem(:) = 0.
    ! nullify not used potential explicitely
    NULLIFY(this%pot)

    !------------------------------- output ---------------------------------!
    ! fargo = 0 (disabled), = 3 (enabled)
!    CALL GetAttr(config, "fargo", this%fargo)

    valwrite = 0
    CALL GetAttr(config, "output/phi", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "phi", &
              this%phi(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
    valwrite = 0
    CALL GetAttr(config, "output/accel_x", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "accel_x", &
              this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1))
    valwrite = 0
    CALL GetAttr(config, "output/accel_y", valwrite, 0)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "accel_y", &
              this%accel(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,2))
    valwrite = 0
    CALL GetAttr(config, "output/mass2D", valwrite, 1)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "mass2D", &
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY)
              this%mass2D(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
#elif defined(HAVE_FFTKEISAN)
              this%mass2D(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX))
#endif
    valwrite = 0
    CALL GetAttr(config, "output/Fmass2D", valwrite, 1)
    IF (valwrite .EQ. 1) &
      CALL SetAttr(IO, "Fmass2D_real", &
              this%Fmass2D_real(Mesh%IMIN:Mesh%IMAX+2,Mesh%JMIN:Mesh%JMAX))

    CALL Info(this, "            .. done initializing")
#endif


  END SUBROUTINE InitGravity_sboxspectral


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
  SUBROUTINE CalcPotential_sboxspectral(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum)&
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: kxp,kyp,Ksqmax,Kmax,Ksq,K,kx,kern
    REAL              :: joff2,jrem2,delt,time0
    INTEGER           :: ip,jp
    INTEGER           :: iopt,ier
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar
    !------------------------------------------------------------------------!
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY) || defined(HAVE_FFTKEISAN)
    !---------------- fourier transformation of density ---------------------!
    ! calculate the shift of the indice at time t                            !
    ! shift density to pretend periodic behavior with interpolation          !

    IF (Mesh%omega .EQ. 0) THEN
      time0 = 0.0       ! avoid division by 0
    ELSE
      time0 = NINT(time*Mesh%Q*Mesh%omega*(Mesh%xmax-Mesh%xmin)/ &
          (Mesh%ymax-Mesh%ymin))*(Mesh%ymax-Mesh%ymin)/(Mesh%Q*Mesh%omega &
          *(Mesh%xmax-Mesh%xmin))
    END IF
    delt = time - time0

    !----------------- shift field to periodic point ------------------------!
    CALL FieldShift(this,Mesh,Physics,delt,pvar(:,:,Physics%DENSITY), &
                    this%mass2D(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))

    !-------------- fourier transform of shifted density --------------------!
! for performance checks
#ifdef _FTRACE
CALL ftrace_region_begin("forward_fft")
#endif

#ifdef HAVE_FFTW
    CALL fftw_execute_dft_r2c(this%plan_r2c, this%mass2D, this%Fmass2D)
#elif HAVE_FFTW_LEGACY
    CALL dfftw_execute_dft_r2c(this%plan_r2c, this%mass2D, this%Fmass2D)
#elif HAVE_FFTKEISAN
    iopt=0    ! forward transformation
    CALL DRC2FT(this%mass2D,Mesh%IMAX-Mesh%IMIN+1,Mesh%JMAX-Mesh%JMIN+1, &
                Mesh%IMAX-Mesh%IMIN+3,iopt,ier)
    IF ( IER .NE. 0 ) THEN
      CALL Error(this, "IER", " FROM DRC2FT")
      STOP
    END IF
    ! move real data in a fortran complex array
    DO i = Mesh%IMIN,Mesh%IMAX/2+1
      DO j = Mesh%JMIN,Mesh%JMAX
        this%Fmass2D(i,j) = CMPLX(this%mass2D(2*i-1,j),this%mass2D(2*i,j))
        this%Fmass2D_real(2*i-1,j) = this%mass2D(2*i-1,j)
        this%Fmass2D_real(2*i,j) = this%mass2D(2*i,j)
      END DO
    END DO
#endif

! performance checks
#ifdef _FTRACE
CALL ftrace_region_end("forward_fft")
#endif

#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY)
    ! turn complex array to real one for output
    DO j = Mesh%JMIN,Mesh%JMAX
      DO i = Mesh%IMIN,Mesh%IMAX/2+1
        this%Fmass2D_real(2*i-1,j) = REAL(REAL(this%Fmass2D(i,j)))
        this%Fmass2D_real(2*i,j)   = REAL(AIMAG(this%Fmass2D(i,j)))
      END DO
    END DO
#endif

    !------------- calculate potential in fourier domain --------------------!
    IF (Mesh%dx .GT. Mesh%dy) THEN
      Ksqmax = 0.5*PI**2/(Mesh%dx**2)
    ! fourier transform of shifted density for different libraries
    ELSE
      Ksqmax = 0.5*PI**2/(Mesh%dy**2)
    END IF

    DO i = Mesh%IMIN-1,Mesh%IMAX/2
      DO j = Mesh%JMIN-1,Mesh%JMAX-1
        ip = MODULO(i+(Mesh%IMAX-Mesh%IMIN+1)/2, &
                    (Mesh%IMAX-Mesh%IMIN+1)) - (Mesh%IMAX-Mesh%IMIN+1)/2
        kxp = 2.*PI*ip/(Mesh%xmax-Mesh%xmin)
        jp = MODULO(j+(Mesh%JMAX-Mesh%JMIN+1)/2,(Mesh%JMAX-Mesh%JMIN+1)) &
                    - (Mesh%JMAX-Mesh%JMIN+1)/2
        kyp = 2.*PI*jp/(Mesh%ymax-Mesh%ymin)

        kx = kxp + Mesh%Q*Mesh%omega*kyp*delt
        Ksq = kx*kx+kyp*kyp

        K = SQRT(Ksq)
        IF ((Ksq .LT. Ksqmax) .AND. (Ksq .GT. 0)) THEN
          kern = -2.*PI*Physics%Constants%GN*(1./K)
          this%Fmass2D(i+1,j+1) = kern * this%Fmass2D(i+1,j+1)
        ELSE
          this%Fmass2D(i+1,j+1) = 0.0
        END IF
      END DO
    END DO

    !----------- fourier transform of shifted density -----------------------!
#ifdef _FTRACE
CALL ftrace_region_begin("backward_fft")
#endif

#ifdef HAVE_FFTW
    CALL fftw_execute_dft_c2r(this%plan_c2r, this%Fmass2D, this%mass2D)
#elif HAVE_FFTW_LEGACY
    CALL dfftw_execute_dft_c2r(this%plan_c2r, this%Fmass2D, this%mass2D)
#elif HAVE_FFTKEISAN
    ! move complex data in a fortran real array
    DO i = Mesh%IMIN,Mesh%IMAX/2+1
      DO j = Mesh%JMIN,Mesh%JMAX
        this%mass2D(2*i-1,j) = REAL(this%Fmass2D(i,j))
        this%mass2D(2*i,j)   = IMAG(this%Fmass2D(i,j))
      END DO
    END DO
    iopt = -1   ! backward transformation
    CALL DRC2FT(this%mass2D,(Mesh%IMAX-Mesh%IMIN+1),(Mesh%JMAX-Mesh%JMIN+1), &
                            (Mesh%IMAX-Mesh%IMIN+3),iopt,ier)
    IF ( IER .NE. 0 ) THEN
      CALL Error(this, "IER", " FROM DRC2FT")
      STOP
    END IF
#endif

#ifdef _FTRACE
CALL ftrace_region_end("backward_fft")
#endif

    !------ calculate final potential with backshift and normalization ------!
    DO j = Mesh%JMIN,Mesh%JMAX
      DO i = Mesh%IMIN,Mesh%IMAX
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY)
        this%phi(i,j) = this%mass2D(i,j)/ &
             ((Mesh%JMAX-Mesh%JMIN+1)*(Mesh%IMAX-Mesh%IMIN+1)) ! no norm. by FFTW !
#elif defined(HAVE_FFTKEISAN)
        this%phi(i,j) = this%mass2D(i,j)                       ! norm. by KEISAN  !
#endif
      END DO
    END DO
    !-------------------------- shift field back ----------------------------!
    CALL FieldShift(this,Mesh,Physics,-delt,this%phi,this%den_ip)

    DO j = Mesh%JMIN,Mesh%JMAX
      DO i = Mesh%IMIN,Mesh%IMAX
        this%phi(i,j) = this%den_ip(i,j)
      END DO
    END DO

    !----- copy values to boundaries in order to calculate acceleration -----!
!CDIR NODEP
    DO j = Mesh%JGMIN, Mesh%JGMAX
      DO i = 1,Mesh%GNUM
        ! northern and southern (periodic)
        this%phi(j,Mesh%JMAX+i) = this%phi(j,Mesh%JMIN+i-1)
        this%phi(j,Mesh%JMIN-i) = this%phi(j,Mesh%JMAX-i+1)
        ! western and eastern (shorn periodic)
        joff2 = -Mesh%Q*Mesh%omega*(Mesh%xmax-Mesh%xmin)*delt/Mesh%dy
        jrem2 = joff2 - FLOOR(joff2)
        this%phi(Mesh%IMIN-i,j) = &
          (1.0-jrem2)*this%phi(Mesh%IMAX-i+1,1+MODULO(j-1+FLOOR(joff2), &
          (Mesh%JMAX-Mesh%JMIN+1))) + jrem2*this%phi(Mesh%IMAX-i+1, &
          1+MODULO(j+FLOOR(joff2), (Mesh%JMAX-Mesh%JMIN+1)))
        joff2 = -joff2
        jrem2 = joff2 - FLOOR(joff2)
        this%phi(Mesh%IMAX+i,j) = &
          (1.0-jrem2)*this%phi(Mesh%IMIN+i-1,1+MODULO(j-1+FLOOR(joff2), &
          (Mesh%JMAX-Mesh%JMIN+1))) + jrem2*this%phi(Mesh%IMIN+i-1, &
          1+MODULO(j+FLOOR(joff2), (Mesh%JMAX-Mesh%JMIN+1)))
      END DO
    END DO

!    ! move potential to borders, but pay attention! Now the borders are deleted
!    DO i = Mesh%IMIN-1,Mesh%IMAX+1
!      DO j = Mesh%JMIN-1,Mesh%JMAX+1
!        this%phi(i,j) = 0.5*(this%phi(i,j)+this%phi(i+1,j))
!      END DO
!    END DO

#endif
  END SUBROUTINE CalcPotential_sboxspectral

  !> \public Calculates the acceleration from potential.
  !!
  !! Calculates
  !! \f[
  !!     \mathbf{a} = -\nabla \Phi.
  !! \f]
  !! Uses second order symmetric difference quotient.
  SUBROUTINE UpdateGravity_sboxspectral(this,Mesh,Physics,time,pvar)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   TYPE(Mesh_TYP)       :: Mesh
   TYPE(Physics_TYP)    :: Physics
   REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                        :: pvar
   REAL                 :: time
   !------------------------------------------------------------------------!
   INTEGER              :: i, j
   REAL                 :: w1,w2
   !------------------------------------------------------------------------!
   INTENT(IN)           :: Mesh,Physics,pvar
   !------------------------------------------------------------------------!
   ! calc potential first
   CALL CalcPotential_sboxspectral(this,Mesh,Physics,time,pvar)

   w1 = 3./48.
   w2 = 30./48.

   !\todo{Here a more robust approximation needs to be searched}
   ! Maybe Physics%VNUM is greater than 2 => set all to zero
   this%accel(:,:,:) = 0.
   DO j = Mesh%JMIN,Mesh%JMAX
     DO i = Mesh%IMIN,Mesh%IMAX
       ! second order approximation
       this%accel(i,j,1) = -1.0*(this%phi(i+1,j)-this%phi(i-1,j))/ &
                            (2*Mesh%dlx(i,j))
       this%accel(i,j,2) = -1.0*(this%phi(i,j+1)-this%phi(i,j-1))/ &
                            (2*Mesh%dly(i,j))
!       ! fourth order
!       this%accel(i,j,1) = -1.0*(w1*this%phi(i-2,j)-w2*this%phi(i-1,j)+ &
!                                 w2*this%phi(i+1,j)-w1*this%phi(i+2,j))/ &
!                                (Mesh%dlx(i,j))
!       this%accel(i,j,2) = -1.0*(w1*this%phi(i,j-2)-w2*this%phi(i,j-1)+ &
!                                 w2*this%phi(i,j+1)-w1*this%phi(i,j+2)) / &
!                                (Mesh%dly(i,j))
     END DO
   END DO
  END SUBROUTINE UpdateGravity_sboxspectral

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
  SUBROUTINE CalcDiskHeight_sboxspectral(this,Mesh,Physics,pvar,cs,height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: cs,height
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: cs2,p,q
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,cs
    INTENT(OUT)       :: height
    !------------------------------------------------------------------------!
    ! pure self-gravitating shearing sheet with external point mass potential
!CDIR COLLAPSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          cs2 = cs(i,j)*cs(i,j)
          p = -SQRTTWOPI*Physics%Constants%GN*pvar(i,j,Physics%DENSITY) &
                 /Mesh%omega**2.
          q = cs2/Mesh%omega**2.
          ! return the new disk height
          height(i,j) = p+SQRT(q+p*p)
       END DO
    END DO

  END SUBROUTINE CalcDiskHeight_sboxspectral

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
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: field
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) :: mass2D
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: cs2,p,q
    REAL              :: delt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,delt,field
    INTENT(OUT)       :: mass2D
    !------------------------------------------------------------------------!

    DO i = Mesh%IMIN,Mesh%IMAX
      this%joff(i)   = -Mesh%Q*Mesh%omega*Mesh%bcenter(i,Mesh%JMIN,1)* &
                  delt/Mesh%dy
      this%jrem(i)  = this%joff(i) - FLOOR(this%joff(i))
    END DO

    DO j = Mesh%JMIN,Mesh%JMAX
      DO i = Mesh%IMIN,Mesh%IMAX
        mass2D(i,j) = &
         (1.0-this%jrem(i))*field(i,1+MODULO(j-1+FLOOR(this%joff(i)), &
          Mesh%JMAX-Mesh%JMIN+1)) + this%jrem(i)*field(i,1+MODULO(j+ &
          FLOOR(this%joff(i)),Mesh%JMAX-Mesh%JMIN+1))
      END DO
    END DO
  END SUBROUTINE FieldShift

  !> \public Closes the gravity term of the shearingsheet spectral solver.
  SUBROUTINE CloseGravity_sboxspectral(this)
   IMPLICIT NONE
   !------------------------------------------------------------------------!
   TYPE(Gravity_TYP), POINTER :: this
   !------------------------------------------------------------------------!
#if defined(HAVE_FFTW) || defined(HAVE_FFTW_LEGACY) || defined(HAVE_FFTKEISAN)
    ! Destroy plans
#ifdef HAVE_FFTW
    CALL fftw_destroy_plan(this%plan_r2c)
    CALL fftw_destroy_plan(this%plan_c2r)
#endif
#ifdef HAVE_FFTW_LEGACY
    CALL dfftw_destroy_plan(this%plan_r2c)
    CALL dfftw_destroy_plan(this%plan_c2r)
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
    END SUBROUTINE CloseGravity_sboxspectral

END MODULE gravity_sboxspectral
