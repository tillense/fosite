!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: agndisk.f90                                                       #
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

!----------------------------------------------------------------------------!
!> Program and data initialization for 3D rotating flow
!!
!! \example agndisk.f90
!!
!! initial condition:
!!
!! hydrostatic equilibrium for slightly non-keplerian flow deviation from
!! Keplerian:
!! \f{eqnarray*}{
!!       v_{\varphi}(s,z) &=& \frac{v_{\mathrm{kep}}(s,z)}{\lambda (z)} \\
!!       \lambda(z)       &=& 1 + \frac{A_0}{1+(z/Z_0)^2} \\
!!       P_\mathrm{c}(s)  &=& P_{\mathrm{0}}\exp{(-s/S_0)},
!! \f}
!! with \f$ P_\mathrm{c}(s) \f$ the pressure in equatorial plane and
!! constant  \f$ A_0 = 0.1 \f$, \f$ Z_0 = 300 \, L_{\mathrm{scale}} \f$,
!! \f$ S_0 = 150 L_{\mathrm{scale}} \f$.
!----------------------------------------------------------------------------!
PROGRAM AGNDISK
  USE fosite_mod
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
  ! general constants
  REAL, PARAMETER :: GN      = 6.6742D-11           ! Newtons grav. constant !
  REAL, PARAMETER :: CC      = 2.99792458D+8        ! speed of light         !
  REAL, PARAMETER :: MSUN    = 1.989D+30            ! solar mass [kg]        !
  REAL, PARAMETER :: AU      = 1.49597870691E+11    ! astronomical unit [m]  !
  REAL, PARAMETER :: PARSEC  = 0.5/PI*1.296E+6 * AU ! parsec [m]             !
  REAL, PARAMETER :: YEAR    = 3.15576E+7           ! Julian year [sec]      !
  REAL, PARAMETER :: ETA_ACC = 0.1                  ! accretion efficiency   !
  ! simulation parameters
  REAL, PARAMETER :: TSIM    = 1.0E+4*YEAR          ! simulation time        !
  REAL, PARAMETER :: M_BH    = 1.0E+6 * MSUN        ! central point mass     !
  REAL, PARAMETER :: M_DISK  = 1.0E+0 * M_BH        ! accretion disk mass    !
  REAL, PARAMETER :: RS = 2 * GN * M_BH / CC**2     ! Schwarzschild radius   !
  REAL, PARAMETER :: LSCALE  = 1.0E+3 * RS          ! typical length scale   !
  REAL, PARAMETER :: ETASGS  = 1.0E-6
  REAL, PARAMETER :: gamma   = 1.4
  ! computational domain
  REAL, PARAMETER :: RMIN    = 5.0D+0 * LSCALE
  REAL, PARAMETER :: RMAX    = 3.0D+2 * LSCALE
  REAL, PARAMETER :: ZMIN    = -1.5D+2 * LSCALE
  REAL, PARAMETER :: ZMAX    = 1.5D+2 * LSCALE
  REAL, PARAMETER :: ZTAN    = 1.0E+2 * LSCALE      ! scale for tancyl. coord.!
  ! initial condition:
  !!   hydrostatic equilibrium for slightly non-keplerian flow
  !!   deviation from Keplerian: v_phi(s,z) = v_kep(s,z) / lambda(z)
  !!                             lambda(z)  = 1 + A0/(1+(z/Z0)**2)
  !!   pressure in equatorial plane: Pc(s) = P0 * EXP(-s/S0)
  REAL, PARAMETER :: Z0      = 3.0E+02 * LSCALE     ! vertical scale         !
  REAL, PARAMETER :: A0      = 1.0E-01              ! velocity deviation     !
  REAL, PARAMETER :: S0      = 1.5E+02 * LSCALE     ! radial scale           !
  REAL, PARAMETER :: B0      = -1.0 ! <1 !          ! slope den. power law   !
  REAL, PARAMETER :: RHO_INF = 1.0E-15              ! den. at infinity       !
  REAL, PARAMETER :: P_INF   = 1.4E-10              ! pressure at infinity   !
  INTEGER, PARAMETER :: PHYS = EULER                ! transport model        !
  ! resolution
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
  INTEGER, PARAMETER :: RRES = 100                  ! res. in r-direction    !
  INTEGER, PARAMETER :: PHIRES = 1                  ! res. in phi-direction  !
  INTEGER, PARAMETER :: ZRES = 101                  ! res. in z-direction    !
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 1000                  ! number of outputs      !
  CHARACTER(LEN=256), PARAMETER :: ODIR &! output directory
                                 = "./"
  CHARACTER(LEN=256), PARAMETER &                   ! output data file name  !
                     :: OFNAME = 'agndisk'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!

  ALLOCATE(Sim)

  CALL Sim%InitFosite()

  CALL MakeConfig(Sim%config)

  CALL Sim%Setup()

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc,Sim%Sources)

  CALL Sim%Run()

  CALL Sim%Finalize()

CONTAINS

   SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               timedisc,fluxes,stemp, grav
    !------------------------------------------------------------------------!
   ! mesh settings
   SELECT CASE(MGEO)
   CASE(CYLINDRICAL)
    mesh => Dict( &
            "meshtype"       / MIDPOINT, &
            "geometry"       / MGEO, &
            "inum"           / RRES, &
            "jnum"           / PHIRES, &
            "knum"           / ZRES, &
            "xmin"           / RMIN, &
            "xmax"           / RMAX, &
            "ymin"           / (0.0), &
            "ymax"           / (2.0*PI), &
            "zmin"           / ZMIN, &
            "zmax"           / ZMAX)

    ! boundary conditions
    boundary => Dict( &
!         "western"  / REFLECTING, &
!         "eastern"  / REFLECTING, &
         "western"  / CUSTOM, &
!         "eastern"  / CUSTOM, &
!         "western"  / NO_GRADIENTS, &
!         "eastern"  / NO_GRADIENTS, &
!         "western"  / FARFIELD, &
         "eastern"  / FARFIELD, &
         "bottomer"  / FARFIELD, &
         "topper"    / FARFIELD, &
!         "southern" / AXIS, &
         "southern" / PERIODIC, &
!         "southern" / REFLECTING, &
         "northern" / PERIODIC)!, &
!         "northern" / REFLECTING)
!         "bottomer" / NO_GRADIENTS, &
!         "northern" / CUSTOM)
!         "southern" / NO_GRADIENTS, &
!         "topper" / NO_GRADIENTS)
!         "southern" / FARFIELD, &
!         "northern" / FARFIELD)
!         "northern" / REFLECTING)

    CASE(SPHERICAL)
      mesh => Dict( &
            "meshtype"       / MIDPOINT, &
            "geometry"       / MGEO, &
            "inum"           / RRES, &
            "jnum"           / ZRES, &
            "knum"           / PHIRES, &
            "xmin"           / RMIN, &
            "xmax"           / RMAX, &
            "ymin"           / (0.0), &
            "ymax"           / (PI), &
            "zmin"           / (0.0), &
            "zmax"           / (2.0*PI))

    ! boundary conditions
    boundary => Dict( &
!         "western"  / REFLECTING, &
!         "eastern"  / REFLECTING, &
         "western"  / CUSTOM, &
!         "eastern"  / CUSTOM, &
!         "western"  / NO_GRADIENTS, &
!         "eastern"  / NO_GRADIENTS, &
!         "western"  / FARFIELD, &
         "eastern"  / FARFIELD, &
         "bottomer"  / PERIODIC, &
         "topper"  / PERIODIC, &
!         "southern" / AXIS, &
         "southern" / AXIS, &
!         "southern" / REFLECTING, &
         "northern" / AXIS)!, &
!         "northern" / REFLECTING)
!         "bottomer" / NO_GRADIENTS, &
!         "northern" / CUSTOM)
!         "southern" / NO_GRADIENTS, &
!         "topper" / NO_GRADIENTS)
!         "southern" / FARFIELD, &
!         "northern" / FARFIELD)
!         "northern" / REFLECTING)


    END SELECT


    ! physics settings
    physics => Dict( &
            "problem"        / PHYS, &
            "gamma"          / gamma, &          ! ratio of specific heats        !
            "mu"             / 6.02E-4, &      ! mean molecular weight          !
            "dpmax"          / 1.0E+3)         ! for advanced time step control !

    ! numerical scheme for flux calculation
    fluxes => Dict(&
            "fluxtype"       / KT, &
            "order"          / LINEAR, &
            "variables"      / PRIMITIVE, &   ! vars. to use for reconstruction!
            "limiter"        / MINMOD, &    ! one of: minmod, monocent,...   !
            "theta"          / 1.2)           ! optional parameter for limiter !

    ! source term due to a point mass
    grav => Dict("stype"    / GRAVITY, &
           "output/accel"    / 1, &
           "pmass/gtype"  / POINTMASS, &    ! grav. accel. of a point mass     !
           "pmass/potential" / NEWTON, &
           "pmass/outbound"  / 1, &
!           "pmass/acclimit"  / (2.2E-9/ETA_ACC/YEAR), &
           "pmass/mass" / M_BH)           ! mass of the accreting object[kg] !


    ! account for energy losses due to radiative cooling
!!$    stemp => Dict(, &
!!$         "stype" / COOLING,, &
!!$         "cvis"  / 0.9)                     ! CFL number for cooling time    !
!!$    CALL SetAttr(sources, "cooling", stemp)


!     ! viscosity source term
!    stemp => Dict("stype"   / VISCOSITY, &
!         "cvis"     / 0.6, &
!!         "vismodel" / ALPHA, &
!!         "dynconst" / 0.1, &
!         "vismodel" / BETA, &
!         "dynconst" / 0.001, &
!           "output/stress" / 1, &
!           "output/kinvis" / 1, &
!           "output/dynvis" / 1)
!!$    CALL SetAttr(sources, "viscosity", stemp)


    ! time discretization settings
    timedisc => Dict(&
         "method"   / SSPRK, &
         "order"    / 5, &
         "cfl"      / 0.4, &
         "stoptime" / TSIM, &
         "dtlimit"  / 1.0E-3, &
         "maxiter"  / 1000000000, &
         "output/density" / 1, &
         "output/xvelocity" / 1, &
         "output/yvelocity" / 1, &
         "output/zvelocity" / 1, &
         "output/pressure" / 1, &
         "output/fluxes"  / 1, &
         "output/rhs" / 1, &
         "output/sgspressure" / 1, &
         "output/geometrical_sources" / 1, &
         "output/external_sources" / 1)

    sources => Dict( &
        "grav"            / grav)

!    ! initialize log input/output
!    logfile => Dict(, &
!         "fileformat" / BINARY, &
!         +"filename"   / (TRIM(ODIR) // TRIM(OFNAME) // "log"), &
!         +"dtwall"     / 3600, &
!         +"filecycles" / 1)

    ! initialize data input/output
    datafile => Dict( &
          "fileformat" / VTK, &
        !  "fileformat" / XDMF, &
          "unit" / 5555, &
          "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
          "stoptime"   / TSIM, &
          "count"      / ONUM, &
          "filecycles" / (ONUM+1))

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)

  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Fluxes,Timedisc,Sources)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Mesh_base),     INTENT(IN)    :: Mesh
    CLASS(Physics_base),  INTENT(INOUT) :: Physics
    CLASS(Fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(Timedisc_base), INTENT(INOUT) :: Timedisc
    CLASS(Sources_base),  POINTER       :: Sources
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp
    CLASS(sources_gravity), POINTER :: gp => null()
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX)   :: velo_periodic
    REAL              :: s,z,tau,lambda,tau1,tau2,dPs,dPz,s1,s2,z1,z2,P0,rho0,cs_inf
    REAL              :: mdisk,am
#ifdef PARALLEL
    INTEGER           :: ierror
    REAL              :: mdisk_local,am_local
#endif
    INTEGER           :: i,j,k,n
    CHARACTER(LEN=20) :: mdisk_str,am_str
    !------------------------------------------------------------------------!

    sp => Sources
    DO
      IF (.NOT.ASSOCIATED(sp)) EXIT 
      SELECT TYPE(sp)
      CLASS IS(sources_gravity)
        gp => sp
        EXIT
      END SELECT
      sp => sp%next
    END DO
    IF (.NOT.ASSOCIATED(sp)) CALL Physics%Error("AGNDISK::InitData","no gravity term initialized")

    ! initial condition (computational domain)
    ! the total disk mass is proportional to the central pressure;
    ! to determine the necessary central pressure for a given
    ! disk mass, the initial condition is computed twice
    Timedisc%pvar%data2d(:,Physics%XVELOCITY) = 0.0
    Timedisc%pvar%data2d(:,Physics%ZVELOCITY) = 0.0
    P0 = 1.0E+10
    rho0 = 1.0E-09
    cs_inf = SQRT(gamma*P_INF/RHO_INF)
    DO n=1,2
       DO k=Mesh%KGMIN,Mesh%KGMAX
         DO j=Mesh%JGMIN,Mesh%JGMAX
            DO i=Mesh%IGMIN,Mesh%IGMAX
             s = abs(Mesh%bccart(i,j,k,1))
             z = Mesh%bccart(i,j,k,3)
!if (s .eq. 0.0 .or. z .eq. 0.0) print *,s,z,i,j

             tau = func_tau(s,z)
             lambda = func_lambda(z)
             Timedisc%pvar%data4d(i,j,k,Physics%PRESSURE) = Pc(tau,rho0,cs_inf,gamma)
             z1 = Mesh%cart%faces(i,j,k,5,3) ! western cell faces z coordinate
             z2 = Mesh%cart%faces(i,j,k,6,3) ! eastern cell faces z coordinate
             tau1 = func_tau(s,z1)
             tau2 = func_tau(s,z2)
             dPz = (Pc(tau2,rho0,cs_inf,gamma)-Pc(tau1,rho0,cs_inf,gamma))/Mesh%dlz%data3d(i,j,k)
             IF (ABS(z).GT.0.0) THEN
!!$                Timedisc%pvar(i,j,Physics%DENSITY) = rho0 * rho(s,z,tau,lambda)
                Timedisc%pvar%data4d(i,j,k,Physics%DENSITY)   = dPz / gz(s,z)
             ELSE
                Timedisc%pvar%data4d(i,j,k,Physics%DENSITY)   = rho0 * (s/S0)**B0
             END IF
             s1 = abs(Mesh%cart%faces(i,j,k,1,1)) ! southern cell faces s coordinate
             s2 = Mesh%cart%faces(i,j,k,2,1) ! northern cell faces s coordinate
             tau1 = func_tau(s1,z)
             tau2 = func_tau(s2,z)
if (tau1 .eq. 0.0 .or. tau2 .eq. 0.0) print *,tau1,tau2,i,j
             dPs = (Pc(tau2,rho0,cs_inf,gamma)-Pc(tau1,rho0,cs_inf,gamma))/Mesh%dlx%data3d(i,j,k)
!             Timedisc%pvar%data4d(i,j,k,Physics%YVELOCITY) = SQRT(MAX(0.0, &
!                  s*(dPs/Timedisc%pvar%data4d(i,j,k,Physics%DENSITY) - gs(s,z))))
             Timedisc%pvar%data4d(i,j,k,Physics%YVELOCITY) = SQRT(MAX(0.0,-s*gs(s,z) / lambda))
          END DO
        END DO
       END DO
      ! 2. initialize gravitational acceleration
!      CALL gp%UpdateGravity(Mesh,Physics,Fluxes,Timedisc%pvar,0.0,0.0)
!       Timedisc%pvar%data4d(:,:,:,2:2+Physics%VDIM) = &
!                Timedisc%GetCentrifugalVelocity(Mesh,Physics,Fluxes,Sources,(/0.,0.,1./),gp%accel%data4d)
       ! compute disk mass
       mdisk = SUM(Timedisc%pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY) &
            * Mesh%volume%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
       ! compute angular momentum
       am = SUM(Mesh%hy%bcenter(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX) &
            * Timedisc%pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%YVELOCITY) &
            * Timedisc%pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,Physics%DENSITY) &
            * Mesh%volume%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

#ifdef PARALLEL
       mdisk_local = mdisk
       CALL MPI_AllReduce(mdisk_local,mdisk,1,DEFAULT_MPI_REAL,MPI_SUM, &
            Mesh%comm_cart,ierror)
       am_local = am
       CALL MPI_AllReduce(am_local,am,1,DEFAULT_MPI_REAL,MPI_SUM, &
            Mesh%comm_cart,ierror)
#endif
!!$       ! adjust central pressure
!!$       P0 = M_DISK / mdisk * P0
       ! adjust density at (S0,Z=0)
       rho0 = M_DISK / mdisk * rho0
    END DO

    ! add velocity perturbations
!!$    CALL RANDOM_SEED
!!$    CALL RANDOM_NUMBER(dv)
!!$    Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
!!$         + (dv(:,:,1)-0.5)*2.0E+2
!!$    Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
!!$         + (dv(:,:,2)-0.5)*2.0E+2


!      ! eastern
       SELECT TYPE(beast => Timedisc%Boundary%Boundary(EAST)%p)
       CLASS IS (boundary_farfield)
         DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            DO i=1,Mesh%GINUM
               beast%data(i,j,k,:) = Timedisc%pvar%data4d(Mesh%IMAX+i,j,k,:)
            END DO
         END DO
        END DO
      END SELECT

      ! southern
      SELECT TYPE(bwest => Timedisc%Boundary%Boundary(WEST)%p)
      CLASS IS(boundary_custom)
        CALL bwest%SetCustomBoundaries(Mesh,Physics, &
          (/CUSTOM_NOGRAD,CUSTOM_OUTFLOW,CUSTOM_KEPLER,CUSTOM_NOGRAD,CUSTOM_NOGRAD/))
      END SELECT

      ! bottom
      SELECT TYPE (bbottom => Timedisc%Boundary%Boundary(BOTTOM)%p)
      CLASS IS(boundary_farfield)
        DO k=1,Mesh%GKNUM
         DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
               bbottom%data(i,j,k,:) = Timedisc%pvar%data4d(i,j,Mesh%KMIN-k,:)
            END DO
         END DO
        END DO
       END SELECT

      ! Top
      SELECT TYPE (btop => Timedisc%Boundary%Boundary(Top)%p)
      CLASS IS(boundary_farfield)
        DO k=1,Mesh%GKNUM
         DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
               btop%data(i,j,k,:) = Timedisc%pvar%data4d(i,j,Mesh%KMAX+k,:)
            END DO
         END DO
        END DO
      END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)

    CALL Timedisc%Boundary%CenterBoundary(Mesh,Physics,0.0,Timedisc%pvar, &
         Timedisc%cvar)

    ! print some information
    WRITE (mdisk_str, '(ES8.2)') mdisk/MSUN
    WRITE (am_str, '(ES8.2)') am
    CALL Mesh%Info(" DATA-----> initial condition: non Keplerian flow")
    CALL Mesh%Info("            disk mass:         " // TRIM(mdisk_str) // " M_sun")
    CALL Mesh%Info("            angular momentum:  " // TRIM(am_str) // " kg/m^2/s")

  END SUBROUTINE InitData

    ELEMENTAL FUNCTION vkepler(sc,zc) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: sc,zc
      REAL :: res
      res = SQRT(GN*M_BH/sc) * (1.0D0+(zc/sc)**2)**(-0.75)
    END FUNCTION vkepler

    ELEMENTAL FUNCTION gs(sc,zc) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: sc,zc
      REAL :: res
      res = -(GN*M_BH/sc**2) / SQRT(1.0D0+(zc/sc)**2)**3
    END FUNCTION gs

    ELEMENTAL FUNCTION gz(sc,zc) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: sc,zc
      REAL :: res
      res = -(GN*M_BH/sc**2) * (zc/sc) / SQRT(1.0D0+(zc/sc)**2)**3
    END FUNCTION gz

    ELEMENTAL FUNCTION func_lambda(zc) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: zc
      REAL :: res
      res = 1.0D0 + A0/(1.0D0+(zc/Z0)**2)
    END FUNCTION func_lambda

    ELEMENTAL FUNCTION func_tau(sc,zc) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: sc,zc
      REAL :: res
      res = sc*SQRT(1.0D0 + (zc/sc)**2*(1.0D0 + (1.0D0+0.5D0*(zc/Z0)**2)/A0))
    END FUNCTION func_tau

    ELEMENTAL FUNCTION Pc(tau,rho_s0,cs_inf,gamma) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: tau,rho_s0,cs_inf,gamma
      REAL :: res
!      res = P0 * EXP(-tau/S0)
      res = P_INF * (1.0 + A0/(1.0+A0) * 1.0/(1.0-B0) * 0.5 * gamma * &
           (CC/cs_inf)**2 * rho_s0/RHO_INF * (RS/S0) * (tau/S0)**(B0-1.0))
    END FUNCTION Pc

!!$    ELEMENTAL FUNCTION dPc(tau) RESULT(res)
!!$      IMPLICIT NONE
!!$      REAL,INTENT(IN) :: tau
!!$      REAL :: res
!!$      res = -1.0D0/S0 * Pc(tau)
!!$    END FUNCTION dPc

!    ELEMENTAL FUNCTION rho(sc,zc,tau,lambda) RESULT(res)
!      IMPLICIT NONE
!      REAL,INTENT(IN) :: sc,zc,tau,lambda
!      REAL :: res
!      res = A0/(1.0+A0) * lambda/(lambda-1.0) * sqrt(1.0+(zc/sc)**2)**3 &
!           * (tau/S0)**B0 * (s/tau)**3
!    END FUNCTION rho


END PROGRAM AGNDISK
