!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: ttauridisk.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006-2009                                                   #
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
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE sources_generic
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
  ! general constants
  REAL, PARAMETER :: GN      = 6.6742D-11     ! Newtons grav. constant       !
  REAL, PARAMETER :: CC      = 2.99792458D+8  ! speed of light               !
  REAL, PARAMETER :: MSUN    = 1.989D+30      ! solar mass [kg]              !
  REAL, PARAMETER :: AU      = 1.49597870691E+11 ! astronomical unit [m]     !
  REAL, PARAMETER :: PARSEC  = 0.5/PI*1.296E+6 * AU ! parsec [m]             !
  ! simulation parameters
  REAL, PARAMETER :: M_BH    = 0.5E+0 * MSUN  ! mass of central point mass   !
  REAL, PARAMETER :: M_DISK  = 1.0E-1 * M_BH  ! accretion disk mass          !
  REAL, PARAMETER :: RS = 2 * GN * M_BH / CC**2 ! Schwarzschild radius       !
  REAL, PARAMETER :: LSCALE  = AU             ! typical length scale         !
  REAL, PARAMETER :: ETASGS = 1.0E-6
  ! computational domain
  REAL, PARAMETER :: RMIN = 2.0D+0 * LSCALE
  REAL, PARAMETER :: RMAX = 3.0D+2 * LSCALE
  REAL, PARAMETER :: ZMIN = -1.3D+2 * LSCALE
  REAL, PARAMETER :: ZMAX = 1.3D+2 * LSCALE
  REAL, PARAMETER :: ZTAN = 1.0E+2 * LSCALE   ! scale for tancyl. coordinates!
  ! initial condition:
  !   hydrostatic equilibrium for slightly non-keplerian flow
  !   deviation from Keplerian: v_phi(s,z) = v_kep(s,z) / lambda(z)
  !                             lambda(z)  = 1 + A0/(1+(z/Z0)**2)
  !   pressure in equatorial plane: Pc(s) = P0 * EXP(-s/S0)
  REAL, PARAMETER :: Z0      = 3.0E+02 * LSCALE ! vertical scale             !
  REAL, PARAMETER :: A0      = 1.0E-01        ! velocity deviation           !
  REAL, PARAMETER :: S0      = 1.5E+02 * LSCALE ! radial scale               !
  REAL, PARAMETER :: B0      = -1.0 ! <1 !    ! slope for density power law  !
  REAL, PARAMETER :: RHO_INF = 1.0E-15        ! density at infinity          !
  REAL, PARAMETER :: P_INF   = 1.4E-10        ! pressure at infinity         !
  !INTEGER, PARAMETER :: PROBLEM = EULER3D_ROTSYMSGS
  INTEGER, PARAMETER :: PROBLEM = EULER3D_ROTSYM
  ! resolution
  INTEGER, PARAMETER :: RRES = 300
  INTEGER, PARAMETER :: ZRES = 252
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)  :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Fluxes, Sim%Timedisc)

  CALL RunFosite(Sim)
  
  CALL CloseFosite(Sim)

CONTAINS

   SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               timedisc,fluxes,stemp
    REAL              :: test_stoptime 
    INTEGER           :: onum
    CHARACTER(LEN=64) :: ofname,lfname
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    
    test_stoptime = 1.0E+10
    onum = 100
    ofname = "nonkepdisk_SGS3"
!    lfname = "nonkepdisklog"
!    test_stoptime = 5.0E+13
!    onum = 500
!    ofname = "/gfs/work-sh1/supas214/ttauridisk_252x300-sgs-ms5e-1-md5e-2"
!    lfname = "/gfs/work-sh1/supas214/ttauridisk_252x300-sgs-ms5e-1-md5e-2log"

   ! mesh settings
    mesh => Dict(&
             "meshtype" / MIDPOINT, &
             "geometry" / CYLINDRICAL, &
             "inum" / ZRES, &        ! resolution in x and            !
             "jnum" / RRES, &        !   y direction                  !             
             "xmin" / ZMIN, &
             "xmax" / ZMAX, &
             "ymin" / RMIN, &
             "ymax" / RMAX)
!!$         +"geometry" / TANCYLINDRICAL, &
!!$             +"inum" / ZRES, &        ! resolution in x and            !
!!$             +"jnum" / RRES, &        !   y direction                  !             
!!$             +"xmin" / ATAN(ZMIN/ZTAN), &
!!$             +"xmax" / ATAN(ZMAX/ZTAN), &
!!$             +"ymin" / RMIN, &
!!$             +"ymax" / RMAX, &
!!$           +"gparam" / ZTAN)
!!$    OPEN(UNIT=1,FILE="volume_508x252-cyl.bin",form='UNFORMATTED',, &
!!$          ACTION='WRITE')
!!$    WRITE(UNIT=1) Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
!!$    CLOSE(UNIT=1)
!!$    STOP


    ! physics settings
    physics => Dict(&
         "problem" / PROBLEM, &
         "gamma"   / 1.4, &          ! ratio of specific heats        !
         "mu"      / 6.02E-4, &      ! mean molecular weight          !
         "dpmax"   / 1.0E+3)         ! for advanced time step control !

    ! numerical scheme for flux calculation
    fluxes => Dict(&
         "fluxtype"  / KT, &
         "order"     / LINEAR, &
         "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
         "limiter"   / MINMOD, &    ! one of: minmod, monocent,...   !
         "theta"     / 1.2)           ! optional parameter for limiter !

     ! boundary conditions
    boundary => Dict( &
!!$         "western"  / REFLECTING, &
!!$         +"eastern"  / REFLECTING, &
!!$         +"western"  / CUSTOM, &
!!$         +"eastern"  / CUSTOM, &
!!$         +"western"  / NO_GRADIENTS, &
!!$         +"eastern"  / NO_GRADIENTS, &
         "western"  / FARFIELD, &
         "eastern"  / FARFIELD, &
!!$         +"southern" / AXIS, &
!!$         +"southern" / REFLECTING, &
!!$         +"northern" / REFLECTING)
         "southern" / CUSTOM, &
!!$         +"northern" / CUSTOM)
!!$         +"southern" / NO_GRADIENTS, &
!!$         +"northern" / NO_GRADIENTS)
!!$         +"southern" / FARFIELD, &
         "northern" / FARFIELD)
!!$         +"northern" / REFLECTING)

    ! source term due to a point mass
    stemp => Dict("stype" / GRAVITY, &
           "pmass/gtype"  / POINTMASS, &            ! grav. accel. of a point mass     !
           "pmass/mass" / M_BH)               ! mass of the accreting object[kg] !

    NULLIFY(sources)
    CALL SetAttr(sources, "grav", stemp)


    ! account for energy losses due to radiative cooling
!!$    stemp => Dict(, &
!!$         "stype" / COOLING,, &
!!$         "cvis"  / 0.9)                     ! CFL number for cooling time    !
!!$    CALL SetAttr(sources, "cooling", stemp)

     stemp => Dict("stype" / SGS, &
                   "output/stress"    / 0, &
                   "output/diffusion" / 0, &
                   "output/rhoeps"    / 0, &
                   "output/sigma"     / 0 )

    ! viscosity source term
    stemp => Dict("stype"   / VISCOSITY, &
         "cvis"     / 0.6, &
!!$         +"vismodel" / ALPHA, &
!!$         +"dynconst" / 0.1)
         "vismodel" / BETA, &
         "dynconst" / 0.001, &
           "output/stress" / 1, &
!           "output/kinvis" / 1, &
           "output/dynvis" / 1)

    IF (PROBLEM.EQ.EULER3D_ROTSYMSGS) THEN
      CALL SetAttr(sources, "sgs", stemp)
    ELSE
      CALL SetAttr(sources, "vis", stemp)
    END IF


    ! time discretization settings
    timedisc => Dict(&
         "method"   / MODIFIED_EULER, &
         "order"    / 3, &
         "cfl"      / 0.4, &
         "stoptime" / test_stoptime, &
         "dtlimit"  / 1.0E-3, &
         "maxiter"  / 1000000000, &
         "output/external_sources" / 1, &
         "output/geometrical_sources" / 1, &
         "output/fluxes" / 0)

!    ! initialize log input/output
!    logfile => Dict(, &
!         "fileformat" / BINARY, &
!         +"filename"   / TRIM(lfname), &
!         +"dtwall"     / 3600, &
!         +"filecycles" / 1)

    ! initialize data input/output
    datafile => Dict("fileformat" / XDMF, &
          "unit" / 5555, &
          "filename"   / TRIM(ofname), &
          "stoptime"   / test_stoptime, &
          "count"      / onum, &
          "filecycles" / (onum+1))

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)

  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Fluxes,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    REAL              :: s,z,tau,lambda,tau1,tau2,dPs,dPz,s1,s2,z1,z2,P0,rho0,cs_inf
    REAL              :: mdisk,am
#ifdef PARALLEL
    INTEGER           :: ierror
    REAL              :: mdisk_local,am_local
#endif
    INTEGER           :: i,j,n
    CHARACTER(LEN=20) :: mdisk_str,am_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition (computational domain)
    ! the total disk mass is proportional to the central pressure;
    ! to determine the necessary central pressure for a given
    ! disk mass, the initial condition is computed twice
    P0 = 1.0E+10
    rho0 = 1.0E-09
    cs_inf = SQRT(Physics%gamma*P_INF/RHO_INF)
    DO n=1,2
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             s = Mesh%bccart(i,j,1)
             z = Mesh%bccart(i,j,2)
if (s .eq. 0.0 .or. z .eq. 0.0) print *,s,z,i,j

             tau = func_tau(s,z)
             lambda = func_lambda(z)
             Timedisc%pvar(i,j,Physics%XVELOCITY) = 0.0
             Timedisc%pvar(i,j,Physics%YVELOCITY) = 0.0
             Timedisc%pvar(i,j,Physics%PRESSURE) = Pc(tau,rho0,cs_inf)
             z1 = Mesh%cart%faces(i,j,1,2) ! western cell faces z coordinate
             z2 = Mesh%cart%faces(i,j,2,2) ! eastern cell faces z coordinate
             tau1 = func_tau(s,z1)
             tau2 = func_tau(s,z2)
             dPz = (Pc(tau2,rho0,cs_inf)-Pc(tau1,rho0,cs_inf))/Mesh%dlx(i,j)
             IF (ABS(z).GT.0.0) THEN
!!$                Timedisc%pvar(i,j,Physics%DENSITY) = rho0 * rho(s,z,tau,lambda)
                Timedisc%pvar(i,j,Physics%DENSITY)   = dPz / gz(s,z)
             ELSE
                Timedisc%pvar(i,j,Physics%DENSITY)   = rho0 * (s/S0)**B0
             END IF
             s1 = Mesh%cart%faces(i,j,3,1) ! southern cell faces s coordinate
             s2 = Mesh%cart%faces(i,j,4,1) ! northern cell faces s coordinate
             tau1 = func_tau(s1,z)
             tau2 = func_tau(s2,z)
             dPs = (Pc(tau2,rho0,cs_inf)-Pc(tau1,rho0,cs_inf))/Mesh%dly(i,j)
             Timedisc%pvar(i,j,Physics%ZVELOCITY) = SQRT(MAX(0.0, &
                  s*(dPs/Timedisc%pvar(i,j,Physics%DENSITY) - gs(s,z))))
!!$             Timedisc%pvar(i,j,Physics%ZVELOCITY) = SQRT(MAX(0.0,-s*gs(s,z) / lambda))
          END DO
       END DO
       ! compute disk mass
       mdisk = SUM(Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
            * Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
       ! compute angular momentum
       am = SUM(Mesh%hz%bcenter(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) &
            * Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%ZVELOCITY) &
            * Timedisc%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Physics%DENSITY) &
            * Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))

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

    ! specific angular momentum
    IF (GetType(Physics).EQ.EULER3D_ROTAMT) &
         Timedisc%pvar(:,:,Physics%ZVELOCITY) = Timedisc%pvar(:,:,Physics%ZVELOCITY) &
         * Mesh%hz%bcenter(:,:)

    IF (GetType(Physics).EQ.EULER3D_ROTSYMSGS) &
       Timedisc%pvar(:,:,Physics%SGSPRESSURE) = ETASGS*Timedisc%pvar(:,:,Physics%PRESSURE)

   IF (GetType(Timedisc%Boundary(WEST)).EQ.FARFIELD) THEN
         DO j=Mesh%JMIN,Mesh%JMAX
            DO i=1,Mesh%GNUM
               Timedisc%Boundary(WEST)%data(i,j,:) = Timedisc%pvar(Mesh%IMIN-i,j,:)
!!$               Timedisc%Boundary(WEST)%data(i,j,Physics%DENSITY)   = RHO_INF
!!$               Timedisc%Boundary(WEST)%data(i,j,Physics%XVELOCITY) = 0.0
!!$               Timedisc%Boundary(WEST)%data(i,j,Physics%YVELOCITY) = 0.0
!!$               Timedisc%Boundary(WEST)%data(i,j,Physics%ZVELOCITY) = 0.0
!!$               Timedisc%Boundary(WEST)%data(i,j,Physics%PRESSURE)  = P_INF
            END DO
         END DO
      END IF
      ! eastern
      IF (GetType(Timedisc%Boundary(EAST)).EQ.FARFIELD) THEN
         DO j=Mesh%JMIN,Mesh%JMAX
            DO i=1,Mesh%GNUM
               Timedisc%Boundary(EAST)%data(i,j,:) = Timedisc%pvar(Mesh%IMAX+i,j,:)
!!$               Timedisc%Boundary(EAST)%data(i,j,Physics%DENSITY)   = RHO_INF
!!$               Timedisc%Boundary(EAST)%data(i,j,Physics%XVELOCITY) = 0.0
!!$               Timedisc%Boundary(EAST)%data(i,j,Physics%YVELOCITY) = 0.0
!!$               Timedisc%Boundary(EAST)%data(i,j,Physics%ZVELOCITY) = 0.0
!!$               Timedisc%Boundary(EAST)%data(i,j,Physics%PRESSURE)  = P_INF
            END DO
         END DO
      END IF
      ! southern
      IF (GetType(Timedisc%Boundary(SOUTH)).EQ.CUSTOM) THEN
         DO j=1,Mesh%GNUM
            DO i=Mesh%IMIN,Mesh%IMAX
               Timedisc%Boundary(SOUTH)%data(i,j,:) = Timedisc%pvar(i,Mesh%JMIN-j,:)
            END DO
         END DO
         WHERE (ABS(Mesh%bccart(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,2)).GE. 2*LSCALE)
            ! axis boundary
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%DENSITY)   = CUSTOM_REFLECT
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%XVELOCITY) = CUSTOM_REFLECT
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%YVELOCITY) = CUSTOM_REFLNEG
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%ZVELOCITY) = CUSTOM_REFLNEG
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%PRESSURE)  = CUSTOM_REFLECT
         ELSEWHERE
            ! subsonic outflow (no gradients)
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%DENSITY)   = CUSTOM_NOGRAD
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%XVELOCITY) = CUSTOM_NOGRAD
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%YVELOCITY) = CUSTOM_NOGRAD
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%ZVELOCITY) = CUSTOM_NOGRAD
            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%PRESSURE)  = CUSTOM_NOGRAD
         END WHERE
!     IF (GetType(Physics).EQ.EULER3D_ROTSYMSGS) THEN
!         WHERE (ABS(Mesh%bccart(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN,2)).GE. 2*LSCALE)
!            ! axis boundary
!            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%SGSPRESSURE)  = CUSTOM_REFLECT
!         ELSEWHERE
!            ! subsonic outflow (no gradients)
!            Timedisc%Boundary(SOUTH)%cbtype(:,Physics%SGSPRESSURE)  = CUSTOM_NOGRAD
!         END WHERE
!     END IF

      END IF
      ! northern
      IF (GetType(Timedisc%Boundary(NORTH)).EQ.FARFIELD) THEN
         DO j=1,Mesh%GNUM
            DO i=Mesh%IMIN,Mesh%IMAX
               Timedisc%Boundary(NORTH)%data(i,j,:) = Timedisc%pvar(i,Mesh%JMAX+j,:)
!!$               Timedisc%Boundary(NORTH)%data(i,j,Physics%DENSITY)   = RHO_INF
!!$               Timedisc%Boundary(NORTH)%data(i,j,Physics%XVELOCITY) = 0.0
!!$               Timedisc%Boundary(NORTH)%data(i,j,Physics%YVELOCITY) = 0.0
!!$               Timedisc%Boundary(NORTH)%data(i,j,Physics%ZVELOCITY) = 0.0
!!$               Timedisc%Boundary(NORTH)%data(i,j,Physics%PRESSURE) &
!!$                    = Timedisc%pvar(i,Mesh%JMAX+j,Physics%PRESSURE)
!!$               Timedisc%Boundary(NORTH)%data(i,j,Physics%PRESSURE)  = P_INF
            END DO
         END DO
    END IF

    CALL Convert2Conservative(Physics,Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN, &
             Mesh%JMAX,Timedisc%pvar,Timedisc%cvar)
    CALL CenterBoundary(Timedisc%boundary,Mesh,Fluxes,Physics,0.0,Timedisc%pvar, &
         Timedisc%cvar)

    ! print some information
    WRITE (mdisk_str, '(ES8.2)') mdisk/MSUN
    WRITE (am_str, '(ES8.2)') am
    CALL Info(Mesh, " DATA-----> initial condition: non Keplerian flow")
    CALL Info(Mesh, "            disk mass:         " // TRIM(mdisk_str) // " M_sun")
    CALL Info(Mesh, "            angular momentum:  " // TRIM(am_str) // " kg/m^2/s")

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

    ELEMENTAL FUNCTION Pc(tau,rho_s0,cs_inf) RESULT(res)
      IMPLICIT NONE
      REAL,INTENT(IN) :: tau,rho_s0,cs_inf
      REAL :: res
!      res = P0 * EXP(-tau/S0)
      res = P_INF * (1.0 + A0/(1.0+A0) * 1.0/(1.0-B0) * 0.5 * Sim%Physics%gamma * &
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

 
END PROGRAM Init
