!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: KHI3D.f90                                                         #
!#                                                                           #
!# Copyright (C) 2006-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jubin Lirawi     <jlirawi@astrophysik.uni-kiel.de>                        #
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
!> Kelvin-Helmholtz instability
!! \author Tobias Illenseer
!! \author Jubin Lirawi
!! \author Lars Boesch

!! References:
!! \cite thomson1871
!! \cite helmholtz1868
!----------------------------------------------------------------------------!
PROGRAM KHI
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 30.0       ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  REAL, PARAMETER    :: RE    = 1.0E+4     ! Reynolds number (HUGE(RE) disables viscosity)
  ! initial condition
  REAL, PARAMETER    :: RHO0 = 1.0         ! outer region: density
  REAL, PARAMETER    :: V0   =-0.5         !   velocity
  REAL, PARAMETER    :: P0   = 2.5         !   pressure
  REAL, PARAMETER    :: RHO1 = 2.0         ! inner region: density
  REAL, PARAMETER    :: V1   = 0.5         !   x-velocity
  REAL, PARAMETER    :: P1   = P0          !   pressure
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
  INTEGER, PARAMETER :: RES  = 8          ! resolution
  REAL, PARAMETER    :: XYZLEN= 1.0        ! spatial extend
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 1           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'KHI'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!
  INTEGER            :: i,j,k,n,clock
  REAL               :: sigma_1,sigma_2,sigma_3,sigma_4,sigma_5,sigma_6
  INTEGER, DIMENSION(:), ALLOCATABLE    :: seed
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: &
         pvar_xy,pvar_yx,pvar_xz,pvar_zx,pvar_yz,pvar_zy, pvar_temp, dv, dv_temp
  !--------------------------------------------------------------------------!
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
#ifdef PARALLEL
  CALL Sim%Warning("KHI3d::main","in parallel mode the symmetry test currently works only for xy-direction")
  TAP_PLAN(1)
#else
  TAP_PLAN(6)
#endif
  ! allocate memory for random seed variable
  CALL RANDOM_SEED(SIZE=n)
  ALLOCATE(seed(n))
  ! initialize the seed with a mix from current time and mpi rank
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + (Sim%GetRank()+1) * (/(i-1, i=1,n)/)
  CALL RANDOM_SEED(PUT=seed)

  ! flow parallel to x-y-plane
  ! 1. along the x-direction
  CALL MakeConfig(Sim, Sim%config)
#ifdef PARALLEL
  ! decompose computational domain along z-direction only
  CALL SetAttr(Sim%config, "/mesh/decomposition", (/1,1,-1/))
#endif
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_xy"))
  CALL Sim%Setup()

  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'xy')
  CALL Sim%Run()
  ! allocate memory and store the results
  ALLOCATE(pvar_xy,SOURCE=Sim%Timedisc%pvar%data4d,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_xy")
  CALL Sim%Finalize(mpifinalize_=.FALSE.)
  DEALLOCATE(Sim)

  ! 2. along the y-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
#ifdef PARALLEL
  ! decompose computational domain along z-direction only
  CALL SetAttr(Sim%config, "/mesh/decomposition", (/1,1,-1/))
#endif
  ! decompose computational domain along z-direction only
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_yx"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'yx')
  CALL Sim%Run()
  ! allocate memory and store the results
  ALLOCATE(pvar_yx,SOURCE=Sim%Timedisc%pvar%data4d,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_yx")
#ifndef PARALLEL
  CALL Sim%Finalize(mpifinalize_=.FALSE.)
  DEALLOCATE(Sim)

  ! flow parallel to x-z-plane
  ! 1. along the x-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_xz"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'xz')
  CALL Sim%Run()
  ! allocate memory and store the results
  ALLOCATE(pvar_xz,SOURCE=Sim%Timedisc%pvar%data4d,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_xz")
  CALL Sim%Finalize(mpifinalize_=.FALSE.)
  DEALLOCATE(Sim)

  ! 2. along the z-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_zx"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'zx')
  CALL Sim%Run()
  ! allocate memory and store the results
  ALLOCATE(pvar_zx,SOURCE=Sim%Timedisc%pvar%data4d,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_zx")
  CALL Sim%Finalize(mpifinalize_=.FALSE.)
  DEALLOCATE(Sim)

  ! flow parallel to y-z-plane
  ! 1. along the z-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_zy"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'zy')
  CALL Sim%Run()
  ! allocate memory and store the results
  ALLOCATE(pvar_zy,SOURCE=Sim%Timedisc%pvar%data4d,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_zy")
  CALL Sim%Finalize(mpifinalize_=.FALSE.)
  DEALLOCATE(Sim)

  ! 1. along the y-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_yz"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'yz')
  CALL Sim%Run()
  ! allocate memory and store the results
  ALLOCATE(pvar_yz,SOURCE=Sim%Timedisc%pvar%data4d,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_yz")
#endif

  ALLOCATE(pvar_temp,MOLD=pvar_yx,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_temp")
  DO k=Sim%Mesh%KGMIN,Sim%Mesh%KGMAX
      pvar_temp(:,:,k,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%DENSITY))
      pvar_temp(:,:,k,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%YVELOCITY))
      pvar_temp(:,:,k,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%XVELOCITY))
      pvar_temp(:,:,k,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%ZVELOCITY))
      pvar_temp(:,:,k,Sim%Physics%PRESSURE)  = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%PRESSURE))
  END DO
  sigma_1 = SQRT(SUM((pvar_temp(:,:,:,:) - pvar_yx(:,:,:,:))**2)/SIZE(pvar_temp))
  DEALLOCATE(pvar_temp)

#ifndef PARALLEL
  ALLOCATE(pvar_temp,MOLD=pvar_zx,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_temp")
  DO j=Sim%Mesh%JGMIN, Sim%Mesh%JGMAX
    pvar_temp(:,j,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%DENSITY))
    pvar_temp(:,j,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%ZVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%YVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%XVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%ENERGY))
  END DO
  sigma_2 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_zx(:,:,:,:))**2)/SIZE(pvar_temp))
  DEALLOCATE(pvar_temp)

  ALLOCATE(pvar_temp,MOLD=pvar_yz,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_temp")
  DO i=Sim%Mesh%IGMIN,Sim%Mesh%IGMAX
    pvar_temp(i,:,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%DENSITY))
    pvar_temp(i,:,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%XVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%ZVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%YVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%ENERGY))
  END DO
  sigma_3 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_yz(:,:,:,:))**2)/SIZE(pvar_temp))
  DEALLOCATE(pvar_temp)

  ALLOCATE(pvar_temp,MOLD=pvar_zy,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_temp")
  DO k=Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
    pvar_temp(:,:,k,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%DENSITY))
    pvar_temp(:,:,k,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%YVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%XVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%ZVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%ENERGY))
  END DO
  sigma_4 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_zy(:,:,:,:))**2)/SIZE(pvar_temp))
  DEALLOCATE(pvar_temp)

  ALLOCATE(pvar_temp,MOLD=pvar_yz,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_temp")
  DO j=Sim%Mesh%JGMIN, Sim%Mesh%JGMAX
    pvar_temp(:,j,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%DENSITY))
    pvar_temp(:,j,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%ZVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%YVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%XVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%ENERGY))
  END DO
  sigma_5 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_yz(:,:,:,:))**2)/SIZE(pvar_temp))
  DEALLOCATE(pvar_temp)

  ALLOCATE(pvar_temp,MOLD=pvar_xz,STAT=Sim%err)
  IF (Sim%err.NE.0) CALL Sim%Error("KHI3D:main","memory allocation failed for pvar_temp")
  DO i=Sim%Mesh%IGMIN, Sim%Mesh%IGMAX
    pvar_temp(i,:,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%DENSITY))
    pvar_temp(i,:,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%XVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%ZVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%YVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%ENERGY))
  END DO
  sigma_6 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_xz(:,:,:,:))**2)/SIZE(pvar_temp))
  DEALLOCATE(pvar_xz,pvar_zx,pvar_yz,pvar_zy,pvar_temp)
#endif

  CALL Sim%Finalize()
  DEALLOCATE(Sim,pvar_xy,pvar_yx,dv,dv_temp,seed)

  TAP_CHECK_SMALL(sigma_1,4*EPSILON(sigma_1),"xy-yx symmetry test")
#ifndef PARALLEL
  TAP_CHECK_SMALL(sigma_2,4*EPSILON(sigma_2),"xz-zx symmetry test")
  TAP_CHECK_SMALL(sigma_3,4*EPSILON(sigma_3),"zy-yz symmetry test")
  TAP_CHECK_SMALL(sigma_4,4*EPSILON(sigma_4),"zx-zy symmetry test")
  TAP_CHECK_SMALL(sigma_5,4*EPSILON(sigma_5),"yx-yz symmetry test")
  TAP_CHECK_SMALL(sigma_6,4*EPSILON(sigma_6),"xy-xz symmetry test")
#endif
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)           :: Sim
    TYPE(Dict_TYP),POINTER  :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               sources, timedisc, fluxes
    REAL                    :: dynvis
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype"            /      MIDPOINT, &
           "geometry"            /          MGEO, &
               "inum"            /           RES, & ! resolution in x,       !
               "jnum"            /           RES, & !   y and                !
               "knum"            /           RES, & !   z direction          !
               "xmin"            / (-0.5*XYZLEN), &
               "xmax"            /  (0.5*XYZLEN), &
               "ymin"            / (-0.5*XYZLEN), &
               "ymax"            /  (0.5*XYZLEN), &
               "zmin"            / (-0.5*XYZLEN), &
               "zmax"            /  (0.5*XYZLEN), &
               "output/dl"       /             0, &
               "output/bh"       /             0, &
               "output/rotation" /             0  &
    )

    ! physics settings
    physics => Dict( &
              "problem" /         EULER, &
              "gamma"   /         GAMMA  &         ! ratio of specific heats !
    )

    ! flux calculation and reconstruction method
    fluxes => Dict( &
             "fluxtype"  /           KT, &
             "order"     /       LINEAR, &
             "variables" / CONSERVATIVE, & ! vars. to use for reconstruction !
             "limiter"   /     MONOCENT, & ! one of: minmod, monocent,...    !
             "theta"     /          1.2  & ! optional parameter for limiter  !
             )

    ! boundary conditions
    boundary => Dict( &
             "western"  / PERIODIC, &
             "eastern"  / PERIODIC, &
             "southern" / PERIODIC, &
             "northern" / PERIODIC, &
             "bottomer" / PERIODIC, &
             "topper"   / PERIODIC  &
    )


    NULLIFY(sources)
    ! viscosity source term
    ! compute dynamic viscosity constant using typical scales and Reynolds number
    dynvis = ABS(RHO0 * XYZLEN * (V0-V1) / RE)
    IF (dynvis.GT.TINY(1.0)) THEN
       sources => Dict( &
          "vis/stype"          /       VISCOSITY, &
          "vis/vismodel"       /       MOLECULAR, &
          "vis/dynconst"       /          dynvis, &
          "vis/bulkconst"      / (-2./3.*dynvis), &
          "vis/output/dynvis"  /               0, &
          "vis/output/stress"  /               0, &
          "vis/output/kinvis"  /               0, &
          "vis/output/bulkvis" /               0  &
       )
    END IF

    ! time discretization settings
    timedisc => Dict( &
               "method"           / MODIFIED_EULER, &
               "order"            /              3, &
               "cfl"              /            0.4, &
               "stoptime"         /           TSIM, &
               "dtlimit"          /        1.0E-10, &
               "maxiter"          /         100000 &
!               "output/energy"  /              1, &
!               "output/xmomentum" /              1, &
!               "output/ymomentum" /              1, &
!               "output/zmomentum" /              1,  &
!               "output/rhs"       /              1 &
    )

    ! initialize data input/output
    datafile => Dict(&
       "fileformat" /                          VTK, &
!        "fileformat" /                         XDMF, &
!        "fileformat" /    GNUPLOT, "filecycles" / 0, &
        "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
        "count"      /                         ONUM  &
    )

    config => Dict( &
             "mesh"     /     mesh, &
             "physics"  /  physics, &
             "boundary" / boundary, &
             "fluxes"   /   fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  /  logfile, &
             "datafile" /  datafile &
    )

    IF (ASSOCIATED(sources)) &
       CALL SetAttr(config, "sources", sources)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Mesh,Physics,Timedisc,flow_dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)  :: Physics
    CLASS(mesh_base)     :: Mesh
    CLASS(timedisc_base) :: Timedisc
    CHARACTER(2)         :: flow_dir
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER              :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics,flow_dir
    INTENT(INOUT)        :: Timedisc
    !------------------------------------------------------------------------!
    SELECT TYPE(pvar=>Timedisc%pvar)
    CLASS IS(statevector_euler)
      IF (.NOT.ALLOCATED(dv)) THEN
         ALLOCATE(dv,MOLD=pvar%velocity%data4d,STAT=Sim%err)
         ! compute velocity perturbations
         CALL RANDOM_NUMBER(dv)
      END IF
      IF (.NOT.ALLOCATED(dv_temp)) ALLOCATE(dv_temp,MOLD=pvar%velocity%data4d,STAT=Sim%err)

      ! initial condition
      SELECT CASE(flow_dir)
      CASE('xy') ! flow along x-direction, density & velocity gradient in y-direction
        ! y and z-velocity vanish everywhere
        pvar%velocity%data2d(:,2) = 0.0
        pvar%velocity%data2d(:,3) = 0.0
        WHERE ((Mesh%bcenter(:,:,:,2).LT.(Mesh%ymin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,2).GT.(Mesh%ymin+0.75*XYZLEN)))
            pvar%density%data3d(:,:,:) = RHO0
            pvar%velocity%data4d(:,:,:,1) = V0
            pvar%pressure%data3d(:,:,:) = P0
        ELSEWHERE
            pvar%density%data3d(:,:,:) = RHO1
            pvar%velocity%data4d(:,:,:,1) = V1
            pvar%pressure%data3d(:,:,:) = P1
        END WHERE
          ! add perturbation to the velocity field
          pvar%velocity%data4d(:,:,:,:) = pvar%velocity%data4d(:,:,:,:) + (dv(:,:,:,:)-0.5)*0.2
      CASE('yx') ! flow along y-direction, density & velocity gradient in x-direction
        ! x and z-velocity vanish everywhere
        pvar%velocity%data2d(:,1) = 0.0
        pvar%velocity%data2d(:,3) = 0.0
        WHERE ((Mesh%bcenter(:,:,:,1).LT.(Mesh%xmin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,1).GT.(Mesh%xmin+0.75*XYZLEN)))
            pvar%density%data3d(:,:,:) = RHO0
            pvar%velocity%data4d(:,:,:,2) = V0
            pvar%pressure%data3d(:,:,:) = P0
        ELSEWHERE
            pvar%density%data3d(:,:,:) = RHO1
            pvar%velocity%data4d(:,:,:,2) = V1
            pvar%pressure%data3d(:,:,:) = P1
        END WHERE
        ! add perturbation to the velocity field
        DO k = Mesh%KGMIN, Mesh%KGMAX
           pvar%velocity%data4d(:,:,k,1) = pvar%velocity%data4d(:,:,k,1) + TRANSPOSE(dv(:,:,k,2)-0.5)*0.2
           pvar%velocity%data4d(:,:,k,2) = pvar%velocity%data4d(:,:,k,2) + TRANSPOSE(dv(:,:,k,1)-0.5)*0.2
           pvar%velocity%data4d(:,:,k,3) = pvar%velocity%data4d(:,:,k,3) + TRANSPOSE(dv(:,:,k,3)-0.5)*0.2
        END DO
      CASE('xz') ! flow along x-direction, density & velocity gradient in z-direction
        ! y and z-velocity vanish everywhere
        pvar%velocity%data2d(:,2) = 0.0
        pvar%velocity%data2d(:,3) = 0.0
        WHERE ((Mesh%bcenter(:,:,:,3).LT.(Mesh%zmin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,3).GT.(Mesh%zmin+0.75*XYZLEN)))
            pvar%density%data3d(:,:,:) = RHO0
            pvar%velocity%data4d(:,:,:,1) = V0
            pvar%pressure%data3d(:,:,:) = P0
        ELSEWHERE
            pvar%density%data3d(:,:,:) = RHO1
            pvar%velocity%data4d(:,:,:,1) = V1
            pvar%pressure%data3d(:,:,:) = P1
        END WHERE
        ! add perturbation to the velocity field
        DO i = Mesh%IGMIN, Mesh%IGMAX
           pvar%velocity%data4d(i,:,:,1) = pvar%velocity%data4d(i,:,:,1) + TRANSPOSE(dv(i,:,:,1)-0.5)*0.2
           pvar%velocity%data4d(i,:,:,2) = pvar%velocity%data4d(i,:,:,2) + TRANSPOSE(dv(i,:,:,3)-0.5)*0.2
           pvar%velocity%data4d(i,:,:,3) = pvar%velocity%data4d(i,:,:,3) + TRANSPOSE(dv(i,:,:,2)-0.5)*0.2
        END DO
      CASE('zx') ! flow along z-direction, density & velocity gradient in x-direction
        ! x and y-velocity vanish everywhere
        pvar%velocity%data2d(:,1) = 0.0
        pvar%velocity%data2d(:,2) = 0.0
        WHERE ((Mesh%bcenter(:,:,:,1).LT.(Mesh%xmin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,1).GT.(Mesh%xmin+0.75*XYZLEN)))
            pvar%density%data3d(:,:,:) = RHO0
            pvar%velocity%data4d(:,:,:,3) = V0
            pvar%pressure%data3d(:,:,:) = P0
        ELSEWHERE
            pvar%density%data3d(:,:,:) = RHO1
            pvar%velocity%data4d(:,:,:,3) = V1
            pvar%pressure%data3d(:,:,:) = P1
        END WHERE
        ! transpose velocity perturbations
        DO k = Mesh%KGMIN,Mesh%KGMAX
          dv_temp(:,:,k,1) = TRANSPOSE(dv(:,:,k,2))
          dv_temp(:,:,k,2) = TRANSPOSE(dv(:,:,k,1))
          dv_temp(:,:,k,3) = TRANSPOSE(dv(:,:,k,3))
        END DO
        ! add perturbation to the velocity field
        DO i = Mesh%IGMIN, Mesh%IGMAX
           pvar%velocity%data4d(i,:,:,1) = pvar%velocity%data4d(i,:,:,1) + TRANSPOSE(dv_temp(i,:,:,1)-0.5)*0.2
           pvar%velocity%data4d(i,:,:,2) = pvar%velocity%data4d(i,:,:,2) + TRANSPOSE(dv_temp(i,:,:,3)-0.5)*0.2
           pvar%velocity%data4d(i,:,:,3) = pvar%velocity%data4d(i,:,:,3) + TRANSPOSE(dv_temp(i,:,:,2)-0.5)*0.2
        END DO
      CASE('zy') ! flow along z-direction, density & velocity gradient in y-direction
        ! x and y-velocity vanish everywhere
        pvar%velocity%data2d(:,1) = 0.0
        pvar%velocity%data2d(:,2) = 0.0
        WHERE ((Mesh%bcenter(:,:,:,2).LT.(Mesh%ymin+0.25*XYZLEN)).OR. &
                (Mesh%bcenter(:,:,:,2).GT.(Mesh%ymin+0.75*XYZLEN)))
            pvar%density%data3d(:,:,:) = RHO0
            pvar%velocity%data4d(:,:,:,3) = V0
            pvar%pressure%data3d(:,:,:) = P0
        ELSEWHERE
            pvar%density%data3d(:,:,:) = RHO1
            pvar%velocity%data4d(:,:,:,3) = V1
            pvar%pressure%data3d(:,:,:) = P1
        END WHERE 
        ! add perturbation to the velocity field
        DO j = Mesh%JGMIN, Mesh%JGMAX
           pvar%velocity%data4d(:,j,:,1) = pvar%velocity%data4d(:,j,:,1) + TRANSPOSE(dv(:,j,:,3)-0.5)*0.2
           pvar%velocity%data4d(:,j,:,2) = pvar%velocity%data4d(:,j,:,2) + TRANSPOSE(dv(:,j,:,2)-0.5)*0.2
           pvar%velocity%data4d(:,j,:,3) = pvar%velocity%data4d(:,j,:,3) + TRANSPOSE(dv(:,j,:,1)-0.5)*0.2
        END DO
      CASE('yz') ! flow along y-direction, density & velocity gradient in z-direction
        ! x and z-velocity vanish everywhere
        pvar%velocity%data2d(:,1) = 0.0
        pvar%velocity%data2d(:,3) = 0.0
        WHERE ((Mesh%bcenter(:,:,:,3).LT.(Mesh%zmin+0.25*XYZLEN)).OR. &
                (Mesh%bcenter(:,:,:,3).GT.(Mesh%zmin+0.75*XYZLEN)))
          pvar%density%data3d(:,:,:) = RHO0
          pvar%velocity%data4d(:,:,:,2) = V0
          pvar%pressure%data3d(:,:,:) = P0
        ELSEWHERE
          pvar%density%data3d(:,:,:) = RHO1
          pvar%velocity%data4d(:,:,:,2) = V1
          pvar%pressure%data3d(:,:,:) = P1
        END WHERE
        ! transpose velocity perturbations
        DO i = Mesh%IGMIN, Mesh%IGMAX
          dv_temp(i,:,:,1) = TRANSPOSE(dv(i,:,:,1))
          dv_temp(i,:,:,2) = TRANSPOSE(dv(i,:,:,3))
          dv_temp(i,:,:,3) = TRANSPOSE(dv(i,:,:,2))
        END DO
        ! add perturbation to the velocity field
        DO k = Mesh%KGMIN, Mesh%KGMAX
           pvar%velocity%data4d(:,:,k,1) = pvar%velocity%data4d(:,:,k,1) + TRANSPOSE(dv_temp(:,:,k,2)-0.5)*0.2
           pvar%velocity%data4d(:,:,k,2) = pvar%velocity%data4d(:,:,k,2) + TRANSPOSE(dv_temp(:,:,k,1)-0.5)*0.2
           pvar%velocity%data4d(:,:,k,3) = pvar%velocity%data4d(:,:,k,3) + TRANSPOSE(dv_temp(:,:,k,3)-0.5)*0.2
        END DO
     CASE DEFAULT
        CALL Mesh%Error("KHI INIT","Unknown flow direction")
      END SELECT
    END SELECT

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability " // flow_dir //"-direction")

  END SUBROUTINE InitData

END PROGRAM KHI
