!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: kepdisk.f90                                                       #
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
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER :: GN      = 6.6742D-11     ! Newtons grav. constant       !
  REAL, PARAMETER :: CC      = 2.99792458D+8  ! speed of light               !
  REAL, PARAMETER :: MSUN    = 1.989D+30      ! solar mass [kg]              !
  REAL, PARAMETER :: AU      = 1.49597870691E+11 ! astronomical unit [m]     !
  REAL, PARAMETER :: PARSEC  = 0.5/PI*1.296E+6 * AU ! parsec [m]             !
  ! simulation parameters
  REAL, PARAMETER :: M_BH    = 1.0E+6 * MSUN  ! mass of central point mass   !
  REAL, PARAMETER :: M_DISK  = 1.0E-1 * M_BH  ! accretion disk mass          !
  REAL, PARAMETER :: RS = 2 * GN * M_BH / CC**2 ! Schwarzschild radius       !
  REAL, PARAMETER :: LSCALE  = RS             ! typical length scale         !
  ! computational domain
  REAL, PARAMETER :: RMIN = 3.0D+1 * LSCALE
  REAL, PARAMETER :: RMAX = 1.0D+3 * LSCALE
  REAL, PARAMETER :: ZMIN = -1.0D+2 * LSCALE
  REAL, PARAMETER :: ZMAX = 1.0D+2 * LSCALE
  REAL, PARAMETER :: ZTAN = 1.0E+2 * LSCALE   ! scale for tancyl. coordinates!
  ! initial condition: Keplerian motion
  ! 
  REAL, PARAMETER :: Z0      = 1.0E+01 * LSCALE ! vertical scale             !
  REAL, PARAMETER :: RHO_INF = 1.0E-15        ! density at infinity          !
  REAL, PARAMETER :: P_INF   = 1.4E-10        ! pressure at infinity         !
  ! resolution
  INTEGER, PARAMETER :: RRES = 100
  INTEGER, PARAMETER :: ZRES = 20
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
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
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               grav, timedisc, fluxes, cool, vis
    REAL              :: x1,x2,y1,y2
    CHARACTER(LEN=128):: ofname,lfname
    INTEGER           :: onum
    REAL              :: test_stoptime
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
!!$         geometry = TANCYLINDRICAL, &
!!$             inum = ZRES, &        ! resolution in x and            !
!!$             jnum = RRES, &        !   y direction                  !             
!!$             xmin = ATAN(ZMIN/ZTAN), &
!!$             xmax = ATAN(ZMAX/ZTAN), &
!!$             ymin = RMIN, &
!!$             ymax = RMAX, &
!!$           gparam = ZTAN)
!!$    OPEN(UNIT=1,FILE="volume_252x300-cyl.bin",form='UNFORMATTED', &
!!$          ACTION='WRITE')
!!$    WRITE(UNIT=1) Mesh%volume(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX)
!!$    CLOSE(UNIT=1)
!!$    STOP
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CYLINDRICAL, &
           "inum"     / ZRES, &
           "jnum"     / RRES, &
           "xmin"     / ZMIN, &
           "xmax"     / ZMAX, &
           "ymin"     / RMIN, &
           "ymax"     / RMAX)
    test_stoptime = 1.0E+6
    onum = 10
    ofname = "kepdisk"
    lfname = "kepdisklog"
!!$    test_stoptime = 1.0E+8
!!$    onum = 1000
!!$    ofname => Dict("/data2/tillense/kepdisk/kepdisk_100x500-md1e5"
!!$    lfname => Dict("/data2/tillense/kepdisk/kepdisk_100x500-md1e5log"

    ! physics settings
    physics => Dict("problem" / EULER3D_ROTSYM, &
              "gamma"   / (5./3.), &        ! ratio of specific heats        !
              "mu"      / 6.02E-4, &        ! mean molecular weight          !
              "dpmax"   / 1.0E+3)          ! for advanced time step control !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
!!$           "variables" / CONSERVATIVE, &  ! vars. to use for reconstruction!
             "variables" / PRIMITIVE, &      ! vars. to use for reconstruction!
             "limiter"   / MINMOD, &         ! one of: minmod, monocent,...   !
!!$           "limiter"   / MONOCENT, &      ! one of: minmod, monocent,...   !
!!$           "limiter"   / OSPRE, &         ! one of: minmod, monocent,...   !
             "theta"     / 1.2)             ! optional parameter for limiter !

    ! boundary conditions
    boundary => Dict( &
!!$         western  = REFLECTING,, &
!!$         eastern  = REFLECTING,, &
         "western"  / NOSLIP, &
           "eastern"  / NOSLIP, &
!!$         western  = CUSTOM,, &
!!$         eastern  = CUSTOM,, &
!!$         western  = NO_GRADIENTS,, &
!!$         eastern  = NO_GRADIENTS,, &
!!$         western  = FARFIELD,, &
!!$         eastern  = FARFIELD,, &
!!$         southern = AXIS,, &
!!$         southern = REFLECTING,, &
!!$         northern = REFLECTING)
           "southern" / NOSLIP, &
           "northern" / NOSLIP &
!!$         southern = CUSTOM,, &
!!$         northern = CUSTOM)
!!$         southern = NO_GRADIENTS,, &
!!$         northern = NO_GRADIENTS)
!!$         southern = FARFIELD,, &
!!$         northern = FARFIELD)
    )
    ! source term due to a point mass
    grav => Dict("stype" / GRAVITY, &
            "pmass/gtype" / POINTMASS, &             ! grav. accel. of a point mass     !
            "pmass/mass"  / M_BH)                   ! mass of the accreting object[kg] !

    ! account for energy losses due to radiative cooling
    cool => Dict("stype" / COOLING, &
           "cvis"  / 0.9)                    ! CFL number for cooling time    !

    ! viscosity source term
    vis => Dict("stype"    / VISCOSITY, &
          "cvis"     / 0.6, &
          "vismodel" / ALPHA, &
          "dynconst" / 0.1)
!!$          "vismodel" / BETA, &
!!$          "dynconst" / 0.001

    sources => Dict(&
           "grav" / grav)
!           "vis" / vis, &
!           "cooling" / cool
         

    ! time discretization settings
    timedisc => Dict("method"   / MODIFIED_EULER, &
               "order"    / 3, &
               "cfl"      / 0.4, &
               "stoptime" / test_stoptime, &
               "dtlimit"  / 1.0E-3, &
               "maxiter"  / 1000000000)

    ! initialize log input/output
!     CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,, &
!          fileformat = BINARY,, &
!          filename   = TRIM(lfname),, &
!          dtwall     = 3600,, &
!          filecycles = 1)

    ! initialize data input/output
    datafile => Dict("fileformat" / VTK, &
               "filename"   / TRIM(ofname), &
!               "stoptime"   / Timedisc%stoptime, &
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

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
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
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Sources_TYP), POINTER :: sp
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: wecart,eacart,socart,nocart
    REAL              :: s,z,P0,dPz,z1,z2
    REAL              :: mdisk,am
#ifdef PARALLEL
    INTEGER           :: ierror
    REAL              :: mdisk_local,am_local
#endif
    INTEGER           :: i,j,n
    CHARACTER(LEN=20) :: mdisk_str,am_str
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition (computational domain)
    ! the total disk mass is proportional to the central pressure;
    ! to determine the necessary central pressure for a given
    ! disk mass, the initial condition is computed twice
    P0 = 1.0E-06
    DO n=1,2
       DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
             s = Mesh%cart%bcenter(i,j,1)
             z = Mesh%cart%bcenter(i,j,2)
             Timedisc%pvar(i,j,Physics%XVELOCITY) = 0.0
             Timedisc%pvar(i,j,Physics%YVELOCITY) = 0.0
             Timedisc%pvar(i,j,Physics%PRESSURE)  = func_P(P0,z)
             z1 = Mesh%cart%faces(i,j,1,2) ! western cell faces z coordinate
             z2 = Mesh%cart%faces(i,j,2,2) ! eastern cell faces z coordinate
             dPz = (func_P(P0,z2)-func_P(P0,z1))/Mesh%dlx(i,j)
             IF (ABS(z).GT.0.0) THEN
                Timedisc%pvar(i,j,Physics%DENSITY)   = dPz / gz(s,z)
             ELSE
                Timedisc%pvar(i,j,Physics%DENSITY)   = func_rho(P0,s,z)
             END IF
             Timedisc%pvar(i,j,Physics%ZVELOCITY) = vkepler(s,z)
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
       ! adjust central pressure
       P0 = M_DISK / mdisk * P0
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

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    ! western
    IF (GetType(Timedisc%Boundary(WEST)).EQ.NOSLIP) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             Timedisc%Boundary(WEST)%data(i,j,:) = Timedisc%pvar(Mesh%IMIN-i,j,:)
          END DO
       END DO
       Timedisc%Boundary(WEST)%data(:,:,Physics%PRESSURE) = -1.0
    END IF
    ! eastern
    IF (GetType(Timedisc%Boundary(EAST)).EQ.NOSLIP) THEN
       DO j=Mesh%JMIN,Mesh%JMAX
          DO i=1,Mesh%GNUM
             Timedisc%Boundary(EAST)%data(i,j,:) = Timedisc%pvar(Mesh%IMAX+i,j,:)
          END DO
       END DO
       Timedisc%Boundary(EAST)%data(:,:,Physics%PRESSURE) = -1.0
    END IF
    ! southern
    IF (GetType(Timedisc%Boundary(SOUTH)).EQ.NOSLIP) THEN
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             Timedisc%Boundary(SOUTH)%data(i,j,:) = Timedisc%pvar(i,Mesh%JMIN-j,:)
          END DO
       END DO
       Timedisc%Boundary(SOUTH)%data(:,:,Physics%PRESSURE) = -1.0
    END IF
    ! northern
    IF (GetType(Timedisc%Boundary(NORTH)).EQ.NOSLIP) THEN
       DO j=1,Mesh%GNUM
          DO i=Mesh%IMIN,Mesh%IMAX
             Timedisc%Boundary(NORTH)%data(i,j,:) = Timedisc%pvar(i,Mesh%JMAX+j,:)
          END DO
       END DO
       Timedisc%Boundary(NORTH)%data(:,:,Physics%PRESSURE) = -1.0
    END IF

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
  
  ELEMENTAL FUNCTION func_P(Pc,z) RESULT(res)
    IMPLICIT NONE
    REAL,INTENT(IN) :: Pc,z
    REAL :: res
    res = Pc * EXP(-0.5*(z/Z0)**2)
  END FUNCTION func_P
  
  ELEMENTAL FUNCTION func_rho(Pc,s,z) RESULT(res)
    IMPLICIT NONE
    REAL,INTENT(IN) :: Pc,s,z
    REAL :: res
    res = (s/Z0)**2 * func_P(Pc,z)/vkepler(s,z)**2
  END FUNCTION func_rho
  
END PROGRAM Init
