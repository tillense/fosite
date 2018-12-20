!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: linearwave.f90                                                    #
!#                                                                           #
!# Copyright (C) 2014                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> Linear Wave
!! \author Manuel Jung
!!
!! References:
!! [1] http://www.astro.virginia.edu/VITA/ATHENA/linear_waves.html
!!
!! Warning: Running all of the convergence tests in this configuration
!!          is taking about 7 hours on 16 cores.
!----------------------------------------------------------------------------!
PROGRAM acousticwave
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
  USE solutions
#include "tap.h"
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
  REAL, PARAMETER    :: LX      = 2.236068
  REAL, PARAMETER    :: LY      = 1.118034
  REAL, PARAMETER    :: ALPHA   = 1.10714871779409050301706546017853
                               !=  ATAN2(LX,LY)
  INTEGER, PARAMETER :: PERIODS = 4
  REAL, PARAMETER    :: GAMMA   = 5./3.
  INTEGER, DIMENSION(6,2), PARAMETER &
                     :: resolutions = RESHAPE((/100, 200, 400, 800, 1600, 3200, &
                                                32,  64, 128, 256,  512,  1024/),&
                                                (/ 6,2 /))
  ! mesh settings
  ! output parameters
  CHARACTER(LEN=256), PARAMETER &        ! output file name
                     :: OFNAME = 'linearwave'
  CHARACTER(LEN=256), PARAMETER &
                     :: resultname = 'result.dat'
  INTEGER, PARAMETER :: ONUM = 1
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  INTEGER            :: res, c, maxwave, dims, fluxtype
  REAL               :: L1, CSISO, lambda
  INTEGER            :: XRES, YRES, WAVE
  LOGICAL            :: existing
  !--------------------------------------------------------------------------!

  INQUIRE(file=resultname,exist=existing)
  IF(existing) THEN
    print *,"Result file is already existing. Delete it first."
    STOP
  END IF

  DO dims=1,2
  !DO dims=2,2
    DO c=0,1
      CSISO = REAL(c) ! c=0: adibatic, c=1: isothermal
      maxwave = 4-c ! Only 3 waves for isothermal simulations
      DO fluxtype=1,3
      !DO fluxtype=3,3
        DO WAVE=1,maxwave
          DO res=1,SIZE(resolutions(:,dims))
            XRES = resolutions(res,dims)
            IF(dims.EQ.2) THEN
              YRES = XRES / 2
              lambda = LX * COS(alpha)
            ELSE
              YRES = 1
              lambda = LX
            END IF

            CALL InitFosite(Sim)

            CALL MakeConfig(Sim, Sim%config)

            CALL SetupFosite(Sim)

            ! set initial condition
            CALL Run(Sim, Sim%Mesh, Sim%Physics, Sim%Timedisc, L1)

            IF(GetRank(Sim).EQ.0) THEN
              OPEN(unit=4, file=resultname, form='formatted', &
                action='write',position='append')
              WRITE(4, '(ES12.4,2X,I1,2X,I1,2X,I4,2X,I4,2X,ES12.4)') &
                CSISO, fluxtype, WAVE, XRES, YRES, L1
              CLOSE(UNIT=4)
            END IF
            !print *,"L1 error:", L1
          END DO
        END DO
      END DO
    END DO
  END DO

  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               timedisc, fluxes
    !------------------------------------------------------------------------!
    CHARACTER(LEN=256) :: csname, wavename, xresname, yresname, fname, fluxname
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / CARTESIAN, &
               "inum" / XRES, &             ! resolution in x and            !
               "jnum" / YRES, &             !   y direction                  !
               "xmin" / 0.0, &
               "xmax" / LX, &
               "ymin" / 0.0, &
               "ymax" / LY)

    IF(CSISO.LE.0.) THEN
      physics => Dict(&
         "problem" / EULER2D, &
         "gamma"   / GAMMA)                   ! ratio of specific heats        !
    ELSE
      physics => Dict(&
         "problem" / EULER2D_ISOTHERM, &
         "cs"      / CSISO)
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict(&
             "fluxtype"  / fluxtype, &
             "order"     / LINEAR, &
             "variables" / PRIMITIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.4, &
           "stoptime" / (PERIODS * lambda / 1.), & ! CS = 1.
           "dtlimit"  / 1.0E-10, &
           "maxiter"  / 100000000)

    ! boundary conditions
    boundary => Dict( &
           "western"  / PERIODIC, &
           "eastern"  / PERIODIC, &
           "southern" / PERIODIC, &
           "northern" / PERIODIC)

    IF(c.EQ.0) THEN
      csname = "adibatic"
    ELSE
      csname = "isotherm"
    END IF
    SELECT CASE(fluxtype)
    CASE(KT)
      fluxname = "KT"
    CASE(HLL)
      fluxname = "HLL"
    CASE(HLLC)
      fluxname = "HLLC"
    CASE(EXACT)
      fluxname = "EXACT"
    END SELECT
    WRITE(wavename, '(I1)') WAVE
    WRITE(xresname, '(I4.4)') XRES
    WRITE(yresname, '(I4.4)') YRES
    fname = TRIM(OFNAME) // '_' // &
            TRIM(csname) // '_' // &
            TRIM(fluxname) // '_' // &
            TRIM(wavename) // '_' // &
            TRIM(xresname) // '_' // &
            TRIM(yresname)

    ! initialize data input/output
    datafile => Dict( &
           "fileformat" / XDMF, &
           "filename"   / (TRIM(fname)), &
           "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig

  SUBROUTINE Run(this,Mesh,Physics,Timedisc,L1)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fosite_TYP)  :: this
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    REAL              :: L1
    !------------------------------------------------------------------------!
    REAL              :: A, rho, u, v, p, cs, H, E
    REAL, DIMENSION(Physics%VNUM,Physics%VNUM) :: R
    REAL, DIMENSION(Physics%VNUM) :: background, sumdU
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                      :: cvar0, dU
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: x
    INTEGER           :: i,j,k,ierr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc,this
    INTENT(OUT)       :: L1
    !------------------------------------------------------------------------!

    A   = 1.E-6
    rho = 1.0
    IF((WAVE.NE.1).AND.(WAVE.NE.maxwave)) THEN
      ! Contant wave
      u = 1.0
    ELSE
      u = 0.0
    END IF
    v   = 0.0
    p   = 1./GAMMA

    IF(YRES.GT.1) THEN
      x = Mesh%bccart(:,:,1) * COS(ALPHA) + Mesh%bccart(:,:,2) * SIN(ALPHA)
    ELSE
      x = Mesh%bccart(:,:,1)
    END IF

    IF(Physics%PRESSURE.GT.0) THEN
      cs = SQRT(GAMMA*p/rho)
      E = p/(GAMMA-1.) + 0.5*rho*(u*u+v*v)
      H = (E + p)/rho
      R(1,:) = (/ 1., u - cs,  v,      H - cs*u /)
      R(2,:) = (/ 1.,      u,  v, 0.5*(u*u+v*v) /)
      R(3,:) = (/ 0.,     0., 1.,            v  /)
      R(4,:) = (/ 1., u + cs,  v,      H + cs*u /)

      background = (/ rho, u, v, p /)
    ELSE
      cs = CSISO

      R(1,:) = (/ 1., u - cs,  v /)
      R(2,:) = (/ 0.,     0., 1. /)
      R(3,:) = (/ 1., u + cs,  v /)

      background = (/ rho, u, v /)
    END IF

    DO k=1,Physics%VNUM
      Timedisc%pvar(:,:,k) = background(k)
    END DO

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)

    DO j=Mesh%JGMIN,Mesh%JGMAX
      DO i=Mesh%IGMIN,Mesh%IGMAX
        Timedisc%cvar(i,j,:) = Timedisc%cvar(i,j,:) &
          + A * R(WAVE,:) * SIN(2.*PI/LAMBDA*x(i,j))
      END DO
    END DO

    IF(YRES.GT.1) THEN
      cvar0 = Timedisc%cvar
      ! Rotate wave by alpha
      Timedisc%cvar(:,:,Physics%XMOMENTUM) &
        = cvar0(:,:,Physics%XMOMENTUM) * COS(ALPHA) &
          - cvar0(:,:,Physics%YMOMENTUM) * SIN(ALPHA)
      Timedisc%cvar(:,:,Physics%YMOMENTUM) &
        = cvar0(:,:,Physics%XMOMENTUM) * SIN(ALPHA) &
          + cvar0(:,:,Physics%YMOMENTUM) * COS(ALPHA)
    END IF

    CALL Convert2Primitive(Physics,Mesh,Timedisc%cvar,Timedisc%pvar)
    CALL Info(Mesh, " DATA-----> initial condition: Linear Wave convergence test")

    cvar0 = Timedisc%cvar

    CALL RunFosite(Sim)

    dU = ABS(Timedisc%cvar - cvar0)

    ! Calculate L1 error
    sumdU = 0.
    DO k=1,Physics%VNUM
      sumdU(k) = SUM(dU(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,k))
#ifdef PARALLEL
      CALL MPI_Allreduce(MPI_IN_PLACE,sumdU(k),1,DEFAULT_MPI_REAL,MPI_SUM,Mesh%comm_cart,ierr)
#endif
    END DO

    L1 = SQRT(SUM(sumdU**2)) / (XRES*YRES)


  END SUBROUTINE

END PROGRAM acousticwave
