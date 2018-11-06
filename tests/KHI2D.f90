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
!!
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
                     :: OFNAME = 'KHI2D2D'
  !--------------------------------------------------------------------------!
!  TYPE(fosite_TYP)   :: Sim
  INTEGER            :: n
  INTEGER            :: i
  REAL               :: sigma,sigma_dens,sigma_xvel,sigma_yvel,sigma_pres
  INTEGER, DIMENSION(:), ALLOCATABLE    :: seed
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: pvar,pvar_1,pvar_2
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  !--------------------------------------------------------------------------!
  TAP_PLAN(1)

  ! allocate memory for random seed variable
  CALL RANDOM_SEED(size=n)
  ALLOCATE(seed(n))

  ! test flow along x-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(SIM,Sim%config)
  CALL Sim%Setup()
  ! allocate memory to store the result of the first run
  ALLOCATE(pvar(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_1(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
  ALLOCATE(pvar_2(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Mesh%KGMIN:Sim%Mesh%KGMAX,Sim%Physics%VNUM))
 
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,.FALSE.,.FALSE.)
  DO i = Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
   ! transpose result and store for comparison with 2nd run
    pvar_1(:,:,i,Sim%Physics%DENSITY)   = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%DENSITY))
    pvar_1(:,:,i,Sim%Physics%XVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%YVELOCITY))
    pvar_1(:,:,i,Sim%Physics%YVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%XVELOCITY))
   ! pvar(:,:,i,Sim%Physics%ZVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%ZVELOCITY))
    pvar_1(:,:,i,Sim%Physics%PRESSURE)  = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%PRESSURE))
  END DO
  CALL Sim%Run()
  DO i = Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
   ! transpose result and store for comparison with 2nd run
    pvar(:,:,i,Sim%Physics%DENSITY)   = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%DENSITY))
    pvar(:,:,i,Sim%Physics%XVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%YVELOCITY))
    pvar(:,:,i,Sim%Physics%YVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%XVELOCITY))
   ! pvar(:,:,i,Sim%Physics%ZVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%ZVELOCITY))
    pvar(:,:,i,Sim%Physics%PRESSURE)  = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%PRESSURE))
  END DO
  CALL Sim%Finalize()
  DEALLOCATE(SIM)
  ALLOCATE(SIM)
  ! test flow along y-direction
  CALL SIM%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_rotate"))
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,.TRUE.,.TRUE.)
  DO i = Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
   ! transpose result and store for comparison with 2nd run
    pvar_2(:,:,i,Sim%Physics%DENSITY)   = (Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%DENSITY))
    pvar_2(:,:,i,Sim%Physics%XVELOCITY) = (Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%XVELOCITY))
    pvar_2(:,:,i,Sim%Physics%YVELOCITY) = (Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%YVELOCITY))
   ! pvar(:,:,i,Sim%Physics%ZVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%ZVELOCITY))
    pvar_2(:,:,i,Sim%Physics%PRESSURE)  = (Sim%Timedisc%pvar%data4d(:,:,i,Sim%Physics%PRESSURE))
  END DO
  CALL Sim%Run()
  ! compare results
  sigma = SQRT(SUM((Sim%Timedisc%pvar%data4d(:,:,:,:)-pvar(:,:,:,:))**2)/SIZE(pvar))
  sigma_dens = SQRT(SUM((pvar_1(:,:,:,1)-pvar_2(:,:,:,1))**2))
  sigma_xvel = SQRT(SUM((pvar_1(:,:,:,2)-pvar_2(:,:,:,2))**2))
  sigma_yvel = SQRT(SUM((pvar_1(:,:,:,3)-pvar_2(:,:,:,3))**2))
  sigma_pres = SQRT(SUM((pvar_1(:,:,:,4)-pvar_2(:,:,:,4))**2))
  
  CALL Sim%Finalize()
  DEALLOCATE(pvar,seed,Sim)
  TAP_CHECK_SMALL(sigma,TINY(sigma),"x-y symmetry test")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite),INTENT(INOUT) :: Sim
    TYPE(Dict_TYP),POINTER  :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               sources, timedisc, fluxes, vis
    REAL                    :: dynvis
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype"            /      MIDPOINT, &
           "geometry"            /          MGEO, &
               "inum"            /           RES, & ! resolution in x,       !
               "jnum"            /           RES, & !   y and                !
               "knum"            /             1, & !   z direction          !
               "xmin"            / (-0.5*XYZLEN), &
               "xmax"            /  (0.5*XYZLEN), &
               "ymin"            / (-0.5*XYZLEN), &
               "ymax"            /  (0.5*XYZLEN), &
               "zmin"            /        (-0.0), &
               "zmax"            /         (0.0), &
               "output/dl"       /             0, &
               "output/bh"       /             0, &
               "output/rotation" /             0  &
    )

    ! physics settings
    physics => Dict( &
              "problem" /       EULER2D, &
              "gamma"   /         GAMMA  &         ! ratio of specific heats !
    )

    ! flux calculation and reconstruction method
    fluxes => Dict( &
             "fluxtype"  /           KT, &
             "order"     /       LINEAR, &
             "variables" / CONSERVATIVE, & ! vars. to use for reconstruction !
             "limiter"   /     MONOCENT, & ! one of: minmod, monocent,...    !
             "theta"     /          1.2, & ! optional parameter for limiter  !
             "output/pfluxes" /       1  &
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
!    dynvis = ABS(RHO0 * XYZLEN * (V0-V1) / RE)
!    IF (dynvis.GT.TINY(1.0)) THEN
!       sources => Dict( &
!          "vis/stype"          /       VISCOSITY, &
!          "vis/vismodel"       /       MOLECULAR, &
!          "vis/dynconst"       /          dynvis, &
!          "vis/bulkconst"      / (-2./3.*dynvis), &
!          "vis/output/dynvis"  /               0, &
!          "vis/output/stress"  /               0, &
!          "vis/output/kinvis"  /               0, &
!          "vis/output/bulkvis" /               0  &
!       )
!    END IF

    ! time discretization settings
    timedisc => Dict( &
               "method"           / MODIFIED_EULER, &
               "order"            /              3, &
               "cfl"              /            0.4, &
               "stoptime"         /           TSIM, &
               "dtlimit"          /        1.0E-10, &
               "maxiter"          /         100000, &
               "output/pressure"  /              1, &
               "output/density"   /              1, &
               "output/xvelocity" /              1, &
               "output/yvelocity" /              1, &
               "output/zvelocity" /              1  &
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

  SUBROUTINE InitData(Mesh,Physics,Timedisc,rotate90deg,reuse_random_seed)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)  :: Physics
    CLASS(mesh_base)     :: Mesh
    CLASS(timedisc_base) :: Timedisc
    LOGICAL, OPTIONAL    :: rotate90deg,reuse_random_seed
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER              :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) :: dv
    INTEGER              :: clock
    REAL, DIMENSION(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX) :: test
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics,rotate90deg,reuse_random_seed
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
    IF (PRESENT(rotate90deg).AND.rotate90deg) THEN
       ! flow along y-direction:
       ! x and z-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = 0.
!       Timedisc%pvar(:,:,:,Physics%ZVELOCITY) = 0.
       WHERE ((Mesh%bcenter(:,:,:,1).LT.(Mesh%xmin+0.25*XYZLEN)).OR. &
            (Mesh%bcenter(:,:,:,1).GT.(SIM%Mesh%xmin+0.75*XYZLEN)))
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = V0
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P0
       ELSEWHERE
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = V1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = P1
       END WHERE
       DO i = Sim%Mesh%KGMIN, Sim%Mesh%KGMAX
         ! add perturbation to the velocity field
         Timedisc%pvar%data4d(:,:,i,Physics%XVELOCITY) = Timedisc%pvar%data4d(:,:,i,Physics%XVELOCITY) &
              + (TRANSPOSE(dv(:,:,i,2))-0.5)*0.02
         Timedisc%pvar%data4d(:,:,i,Physics%YVELOCITY) = Timedisc%pvar%data4d(:,:,i,Physics%YVELOCITY) &
              + (TRANSPOSE(dv(:,:,i,1))-0.5)*0.02
!         Timedisc%pvar(:,:,i,Physics%ZVELOCITY) = Timedisc%pvar%data4d(:,:,i,Physics%ZVELOCITY) &
!              + (TRANSPOSE(dv(:,:,i,3))-0.5)*0.02
       END DO
    ELSE
       ! flow along x-direction:
       ! y and z-velocity vanish everywhere
       Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = 0.
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
            + (dv(:,:,:,1)-0.5)*0.02
       Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) &
            + (dv(:,:,:,2)-0.5)*0.02
!       Timedisc%pvar(:,:,:,Physics%ZVELOCITY) = Timedisc%pvar%data4d(:,:,:,Physics%ZVELOCITY) &
!            + (dv(:,:,:,3)-0.5)*0.02
    END IF
    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar%data4d,Timedisc%cvar%data4d)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability")
  END SUBROUTINE InitData

END PROGRAM KHI
