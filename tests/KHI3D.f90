!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: KHI.f90                                                           #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
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
  INTEGER, PARAMETER :: RES  = 10          ! resolution
  REAL, PARAMETER    :: XYZLEN= 1.0        ! spatial extend
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'KHI'
  !--------------------------------------------------------------------------!
!  TYPE(fosite_TYP)   :: Sim
  INTEGER            :: n
  INTEGER            :: i,j,k
  REAL               :: sigma_1,sigma_2,sigma_3,sigma_4,sigma_5,sigma_6
  INTEGER, DIMENSION(:), ALLOCATABLE    :: seed
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: &
         pvar_xy,pvar_xz,pvar_yx,pvar_yz,pvar_zx,pvar_zy, pvar_temp
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!
  TAP_PLAN(6)

  ! allocate memory for random seed variable
  CALL RANDOM_SEED(size=n)
  ALLOCATE(seed(n))

  ! test flow along x-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_xy"))
  CALL Sim%Setup()
  ! allocate memory to store the result of the first run
  ALLOCATE(pvar_xy(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_xz(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_yx(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_yz(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_zx(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_zy(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_temp(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'xy',.FALSE.)
  CALL Sim%Run()
  pvar_xy(:,:,:,:) = SIM%Timedisc%pvar%data4d(:,:,:,:)
  CALL Sim%Finalize(.FALSE.)
  DEALLOCATE(SIM)

  ! test flow along y-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_yx"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'yx',.TRUE.)
  CALL Sim%Run()
  pvar_yx(:,:,:,:) = Sim%Timedisc%pvar%data4d(:,:,:,:)
  CALL Sim%Finalize(.FALSE.)
  DEALLOCATE(Sim)

  ! test flow along z-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_xz"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'xz',.TRUE.)
  CALL Sim%Run()
  pvar_xz(:,:,:,:) = Sim%Timedisc%pvar%data4d(:,:,:,:)
  CALL Sim%Finalize(.FALSE.)
  DEALLOCATE(Sim)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_zx"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'zx',.TRUE.)
  CALL Sim%Run()
  pvar_zx(:,:,:,:) = Sim%Timedisc%pvar%data4d(:,:,:,:)
  CALL Sim%Finalize(.FALSE.)
  DEALLOCATE(Sim)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_zy"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'zy',.TRUE.)
  CALL Sim%Run()
  pvar_zy(:,:,:,:) = Sim%Timedisc%pvar%data4d(:,:,:,:)
  CALL Sim%Finalize(.FALSE.)
  DEALLOCATE(Sim)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_yz"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,'yz',.TRUE.)
  CALL Sim%Run()
  pvar_yz(:,:,:,:) = Sim%Timedisc%pvar%data4d(:,:,:,:)

   DO k=Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
    pvar_temp(:,:,k,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%DENSITY))
    pvar_temp(:,:,k,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%YVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%XVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%ZVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_xy(:,:,k,Sim%Physics%ENERGY))
  END DO
  sigma_1 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_yx(:,:,:,:))**2)/SIZE(pvar_temp))

  DO j=Sim%Mesh%JGMIN, Sim%Mesh%JGMAX
    pvar_temp(:,j,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%DENSITY))
    pvar_temp(:,j,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%ZVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%YVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%XVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_xz(:,j,:,Sim%Physics%ENERGY))
  END DO
  sigma_2 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_zx(:,:,:,:))**2)/SIZE(pvar_temp))

 DO i=Sim%Mesh%IGMIN, Sim%Mesh%IGMAX
    pvar_temp(i,:,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%DENSITY))
    pvar_temp(i,:,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%XVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%ZVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%YVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_zy(i,:,:,Sim%Physics%ENERGY))
  END DO
  sigma_3 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_yz(:,:,:,:))**2)/SIZE(pvar_temp))

 DO k=Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
    pvar_temp(:,:,k,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%DENSITY))
    pvar_temp(:,:,k,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%YVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%XVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%ZVELOCITY))
    pvar_temp(:,:,k,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_zx(:,:,k,Sim%Physics%ENERGY))
  END DO
  sigma_4 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_zy(:,:,:,:))**2)/SIZE(pvar_temp))

  DO j=Sim%Mesh%JGMIN, Sim%Mesh%JGMAX
    pvar_temp(:,j,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%DENSITY))
    pvar_temp(:,j,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%ZVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%YVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%XVELOCITY))
    pvar_temp(:,j,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_yx(:,j,:,Sim%Physics%ENERGY))
  END DO
  sigma_5 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_yz(:,:,:,:))**2)/SIZE(pvar_temp))

 DO i=Sim%Mesh%IGMIN, Sim%Mesh%IGMAX
    pvar_temp(i,:,:,Sim%Physics%DENSITY)   = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%DENSITY))
    pvar_temp(i,:,:,Sim%Physics%XVELOCITY) = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%XVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%YVELOCITY) = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%ZVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ZVELOCITY) = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%YVELOCITY))
    pvar_temp(i,:,:,Sim%Physics%ENERGY)    = TRANSPOSE(pvar_xy(i,:,:,Sim%Physics%ENERGY))
  END DO
  sigma_6 = SQRT(SUM((pvar_temp(:,:,:,:) -  pvar_xz(:,:,:,:))**2)/SIZE(pvar_temp))

  CALL Sim%Finalize(.TRUE.)
  DEALLOCATE(Sim,pvar_xy,pvar_xz,pvar_yx,pvar_yz,pvar_zx,pvar_zy,pvar_temp,seed)
  TAP_CHECK_SMALL(sigma_1,3*EPSILON(sigma_1),"xy-yx symmetry test")
  TAP_CHECK_SMALL(sigma_2,3*EPSILON(sigma_2),"xz-zx symmetry test")
  TAP_CHECK_SMALL(sigma_3,3*EPSILON(sigma_3),"zy-yz symmetry test")
  TAP_CHECK_SMALL(sigma_4,3*EPSILON(sigma_4),"zx-zy symmetry test")
  TAP_CHECK_SMALL(sigma_5,3*EPSILON(sigma_5),"yx-yz symmetry test")
  TAP_CHECK_SMALL(sigma_6,3*EPSILON(sigma_6),"xy-xz symmetry test")


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
!        "fileformat" /                          HDF, &
        "fileformat" /                          VTK, &
!        "fileformat" /                         XDMF, &
!        "fileformat" /    GNUPLOT, "filecycles" / 0, &
!        "fileformat" /                       BINARY, &
!        "fileformat" /                       NETCDF, &
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




  SUBROUTINE InitData(Mesh,Physics,Timedisc,flow_dir,reuse_random_seed)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)  :: Physics
    CLASS(mesh_base)     :: Mesh
    CLASS(timedisc_base) :: Timedisc
    LOGICAL, OPTIONAL    :: reuse_random_seed
    CHARACTER(2), OPTIONAL  :: flow_dir
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER              :: i,j,k
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) :: dv,dv_temp
    INTEGER              :: clock
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics,flow_dir,reuse_random_seed
    INTENT(INOUT)        :: Timedisc
    !------------------------------------------------------------------------!
    ! Seed the random number generator
    IF (.NOT.(PRESENT(reuse_random_seed).AND.reuse_random_seed)) THEN
       ! initialize with a mix from current time and mpi rank
       CALL SYSTEM_CLOCK(COUNT=clock)
       seed = clock + (Timedisc%getrank()+1) * (/(i-1, i=1,n)/)
    END IF
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(dv)
    ! initial condition
    IF (PRESENT(flow_dir).AND.flow_dir.EQ.'yx') THEN
       ! flow along yx-direction:
       print *,"yx-Direction"
       ! x and z-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
       Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.

       WHERE ((Mesh%bcenter(:,:,:,1).LT.(Mesh%xmin+0.25*XYZLEN)).OR. &
            (Mesh%bcenter(:,:,:,1).GT.(Mesh%xmin+0.75*XYZLEN)))
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = V0
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
       ELSEWHERE
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = V1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
       END WHERE
        DO k = Mesh%KGMIN, Mesh%KGMAX
         ! add perturbation to the velocity field
         Timedisc%pvar%data4d(:,:,k,Physics%XVELOCITY) = Timedisc%pvar%data4d(:,:,k,Physics%XVELOCITY) &
              + TRANSPOSE(dv(:,:,k,2)-0.5)*0.2
         Timedisc%pvar%data4d(:,:,k,Physics%YVELOCITY) = Timedisc%pvar%data4d(:,:,k,Physics%YVELOCITY) &
              + TRANSPOSE(dv(:,:,k,1)-0.5)*0.2
         Timedisc%pvar%data4d(:,:,k,Physics%ZVELOCITY) = Timedisc%pvar%data4d(:,:,k,Physics%ZVELOCITY) &
              + TRANSPOSE(dv(:,:,k,3)-0.5)*0.2
       END DO
    ! initial condition
    ELSEIF (PRESENT(flow_dir).AND.flow_dir.EQ.'zx') THEN
       ! flow along zx-direction:
       print *,"zx-Direction"
       ! x and y-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
       Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.

      WHERE ((Mesh%bcenter(:,:,:,1).LT.(Mesh%xmin+0.25*XYZLEN)).OR. &
            (Mesh%bcenter(:,:,:,1).GT.(Mesh%xmin+0.75*XYZLEN)))
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = V0
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
       ELSEWHERE
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = V1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
       END WHERE
       DO k = Mesh%KGMIN,Mesh%KGMAX
        dv_temp(:,:,k,1) = TRANSPOSE(dv(:,:,k,2))
        dv_temp(:,:,k,2) = TRANSPOSE(dv(:,:,k,1))
        dv_temp(:,:,k,3) = TRANSPOSE(dv(:,:,k,3))
       END DO
       DO i= Mesh%IGMIN, Mesh%IGMAX
        ! add perturbation to the velocity field
        Timedisc%pvar%data4d(i,:,:,Physics%XVELOCITY) = Timedisc%pvar%data4d(i,:,:,Physics%XVELOCITY) &
             + TRANSPOSE(dv_temp(i,:,:,1)-0.5)*0.2
        Timedisc%pvar%data4d(i,:,:,Physics%YVELOCITY) = Timedisc%pvar%data4d(i,:,:,Physics%YVELOCITY) &
             + TRANSPOSE(dv_temp(i,:,:,3)-0.5)*0.2
        Timedisc%pvar%data4d(i,:,:,Physics%ZVELOCITY) = Timedisc%pvar%data4d(i,:,:,Physics%ZVELOCITY) &
             + TRANSPOSE(dv_temp(i,:,:,2)-0.5)*0.2
      END DO
    ELSEIF (PRESENT(flow_dir).AND.flow_dir.EQ.'xy') THEN
       ! flow along xy-direction:
       print *,"xy-Direction"
       ! y and z-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
       Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
        WHERE ((Mesh%bcenter(:,:,:,2).LT.(Mesh%ymin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,2).GT.(Mesh%ymin+0.75*XYZLEN)))
            Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
            Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = V0
            Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
          ELSEWHERE
            Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
            Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = V1
            Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
          END WHERE
        ! add perturbation to the velocity field
        Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) &
             + (dv(:,:,:,1)-0.5)*0.2
        Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) &
             + (dv(:,:,:,2)-0.5)*0.2
        Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) &
             + (dv(:,:,:,3)-0.5)*0.2

     ELSEIF (PRESENT(flow_dir).AND.flow_dir.EQ.'xz') THEN
       ! flow along xz-direction:
       print *,"xz-Direction"
       ! y and z-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
       Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
       WHERE ((Mesh%bcenter(:,:,:,3).LT.(Mesh%zmin+0.25*XYZLEN)).OR. &
            (Mesh%bcenter(:,:,:,3).GT.(Mesh%zmin+0.75*XYZLEN)))     
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = V0
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
       ELSEWHERE
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = V1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
       END WHERE
       DO i = Mesh%IGMIN, Mesh%IGMAX
        ! add perturbation to the velocity field
        Timedisc%pvar%data4d(i,:,:,Physics%XVELOCITY) = Timedisc%pvar%data4d(i,:,:,Physics%XVELOCITY) &
             + TRANSPOSE(dv(i,:,:,1)-0.5)*0.2
        Timedisc%pvar%data4d(i,:,:,Physics%YVELOCITY) = Timedisc%pvar%data4d(i,:,:,Physics%YVELOCITY) &
             + TRANSPOSE(dv(i,:,:,3)-0.5)*0.2
        Timedisc%pvar%data4d(i,:,:,Physics%ZVELOCITY) = Timedisc%pvar%data4d(i,:,:,Physics%ZVELOCITY) &
             + TRANSPOSE(dv(i,:,:,2)-0.5)*0.2
      END DO
   ELSEIF (PRESENT(flow_dir).AND.flow_dir.EQ.'zy') THEN
       ! flow along zy-direction:
       print *,"zy-Direction"
       ! x and y-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
       Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
       WHERE ((Mesh%bcenter(:,:,:,2).LT.(Mesh%ymin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,2).GT.(Mesh%ymin+0.75*XYZLEN)))     
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = V0
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
        ELSEWHERE
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = V1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
        END WHERE 
        DO j = Mesh%JGMIN, Mesh%JGMAX
        ! add perturbation to the velocity field
        Timedisc%pvar%data4d(:,j,:,Physics%XVELOCITY) = Timedisc%pvar%data4d(:,j,:,Physics%XVELOCITY) &
             + TRANSPOSE(dv(:,j,:,3)-0.5)*0.2
        Timedisc%pvar%data4d(:,j,:,Physics%YVELOCITY) = Timedisc%pvar%data4d(:,j,:,Physics%YVELOCITY) &
             + TRANSPOSE(dv(:,j,:,2)-0.5)*0.2
        Timedisc%pvar%data4d(:,j,:,Physics%ZVELOCITY) = Timedisc%pvar%data4d(:,j,:,Physics%ZVELOCITY) &
             + TRANSPOSE(dv(:,j,:,1)-0.5)*0.2
       END DO
   ELSEIF (PRESENT(flow_dir).AND.flow_dir.EQ.'yz') THEN
       ! flow along yz-direction:
       print *,"yz-Direction"
       ! x and z-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
       Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) = 0.
       WHERE ((Mesh%bcenter(:,:,:,3).LT.(Mesh%zmin+0.25*XYZLEN)).OR. &
              (Mesh%bcenter(:,:,:,3).GT.(Mesh%zmin+0.75*XYZLEN)))     
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = V0
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
        ELSEWHERE
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = V1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
        END WHERE
        DO i = Mesh%IGMIN, Mesh%IGMAX
          dv_temp(i,:,:,1) = TRANSPOSE(dv(i,:,:,1))
          dv_temp(i,:,:,2) = TRANSPOSE(dv(i,:,:,3))
          dv_temp(i,:,:,3) = TRANSPOSE(dv(i,:,:,2))
        END DO
        DO k = Mesh%KGMIN, Mesh%KGMAX
        ! add perturbation to the velocity field
        Timedisc%pvar%data4d(:,:,k,Physics%XVELOCITY) = Timedisc%pvar%data4d(:,:,k,Physics%XVELOCITY) &
             + TRANSPOSE(dv_temp(:,:,k,2)-0.5)*0.2
        Timedisc%pvar%data4d(:,:,k,Physics%YVELOCITY) = Timedisc%pvar%data4d(:,:,k,Physics%YVELOCITY) &
             + TRANSPOSE(dv_temp(:,:,k,1)-0.5)*0.2
        Timedisc%pvar%data4d(:,:,k,Physics%ZVELOCITY) = Timedisc%pvar%data4d(:,:,k,Physics%ZVELOCITY) &
             + TRANSPOSE(dv_temp(:,:,k,3)-0.5)*0.2
      END DO
    ELSE 
      CALL Mesh%Error("KHI INIT","Unknown flow direction")
    END IF

    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability")

  END SUBROUTINE InitData

END PROGRAM KHI
