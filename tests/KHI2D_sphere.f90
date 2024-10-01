!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: KHI2D_sphere.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jubin Lirawi      <jlirawi@astrophysik.uni-kiel.de>                       #
!# Lars Boesch       <lboesch@astrophysik.uni-kiel.de>                       #
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
  ! simulation parameters                       !                            !
  REAL, PARAMETER    :: TSIM  = 10.0            ! simulation time            !
  REAL, PARAMETER    :: GAMMA = 1.4             ! ratio of specific heats    !
  REAL, PARAMETER    :: RE    = 1.0E+4          ! Reynolds number (HUGE(RE) disables viscosity)
  REAL, PARAMETER    :: OMEGA   = 0.1           ! angular speed of rotational frame
  ! initial condition                           !                            !
  REAL, PARAMETER    :: RHO0 = 1.0              ! outer region: density      !
  REAL, PARAMETER    :: V0   = 0.0              !   velocity                 !
  REAL, PARAMETER    :: P0   = 2.5              !   pressure                 !
  REAL, PARAMETER    :: RHO1 = 2.0              ! inner region: density      !
  REAL, PARAMETER    :: V1   = 1.0              !   velocity                 !
  REAL, PARAMETER    :: P1   = P0               !   pressure                 !
  ! mesh settings                               !                            !
  INTEGER, PARAMETER :: MGEO = SPHERICAL_PLANET ! geometry of the mesh       !
  INTEGER, PARAMETER :: XRES   = 20             ! x-resolution               !
  INTEGER, PARAMETER :: YRES   = 40             ! y-resolution               !
  REAL, PARAMETER :: GPAR   = 1.0               ! radius of sphere           !
  ! output file parameter                       !                            !
  INTEGER, PARAMETER :: ONUM = 10               ! number of output data sets !
  CHARACTER(LEN=256), PARAMETER &               ! output data dir            !
                     :: ODIR = './'             !                            !
  CHARACTER(LEN=256), PARAMETER &               ! output data file name      !
                     :: OFNAME = 'KHI2D_sphere' !                            !
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  LOGICAL            :: ok
  !--------------------------------------------------------------------------!
  TAP_PLAN(1)

  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim,Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  CALL Sim%Run()
  ok = .NOT.Sim%aborted
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

  TAP_CHECK(ok, "stoptime reached")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite),INTENT(INOUT) :: Sim
    TYPE(Dict_TYP),POINTER  :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes, vis, rotframe
    REAL                    :: dynvis
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype"   /       MIDPOINT, &
           "geometry"   /           MGEO, &
           "fargo/method" / 0, &
           "decomposition" /  (/ -1, 1, 1/), &
               "inum"   /           XRES, & ! resolution in x,          !
               "jnum"   /           YRES, & !               y and       !
               "knum"   /              1, & !               z direction !
               "xmin"   /          (0.1), &
               "xmax"   /       (PI-0.1), &
               "ymin"   /          (0.0), &
               "ymax"   /       (2.0*PI), &
               "zmin"   /         GPAR, &
               "zmax"   /         GPAR, &
!                "zmin"   /           GPAR, &
!                "zmax"   /    (GPAR + .1), &
               "gparam" /           GPAR  &
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
!              "variables" / PRIMITIVE, & ! vars. to use for reconstruction !
             "limiter"   /     MONOCENT, & ! one of: minmod, monocent,...    !
!              "output/pfluxes" /       1, &
!              "output/pstates" /       1, &
!              "output/cstates" /       1, &
             "theta"     /          1.2 &  ! optional parameter for limiter  !
             )

    ! boundary conditions
    boundary => Dict( &
               "western"  / REFLECTING, &
               "eastern"  / REFLECTING, &
               "southern" / PERIODIC, &
               "northern" / PERIODIC, &
               "bottomer" / REFLECTING, &
               "topper"   / REFLECTING  &
    )

    NULLIFY(vis,rotframe)
    ! viscosity source term
    ! compute dynamic viscosity constant using typical scales and Reynolds number
    dynvis = ABS(RHO0 * 2*PI * (V0-V1) / RE)
    IF (RE.LT.HUGE(RE)) THEN
        vis => Dict( &
          "stype"          /       VISCOSITY, &
          "vismodel"       /       MOLECULAR, &
          "dynconst"       /          dynvis, &
          "bulkconst"      / (-2./3.*dynvis), &
          "output/dynvis"  /               0, &
          "output/stress"  /               0, &
          "output/kinvis"  /               0, &
          "output/bulkvis" /               0  &
        )
    END IF
    ! activate inertial forces due to rotating frame if OMEGA > 0
    IF (OMEGA.GT.TINY(OMEGA)) THEN
       rotframe => Dict("stype" / ROTATING_FRAME, &
                        "disable_centaccel" / 1)
       CALL SetAttr(mesh,"omega",OMEGA)
    END IF

    ! time discretization settings
    timedisc => Dict( &
               "method"           / MODIFIED_EULER, &
!               "method"           / SSPRK, &
               "order"            /              3, &
               "cfl"              /            0.4, &
               "stoptime"         /           TSIM, &
               "tol_rel"          /         0.0095, &
               "dtlimit"          /        1.0E-15, &
               "maxiter"          /      100000000, &
               "output/pressure"  /              1, &
               "output/density"   /              1, &
               "output/xvelocity" /              1, &
               "output/yvelocity" /              1  &
    )

    ! initialize data input/output
    datafile => Dict(&
!         "fileformat" /                          VTK, &
        "fileformat" /                         XDMF, &
        "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
        "count"      /                         ONUM  &
    )

    config => Dict( &
             "mesh"     /     mesh, &
             "physics"  /  physics, &
             "boundary" / boundary, &
             "fluxes"   /   fluxes, &
             "timedisc" / timedisc, &
             "datafile" /  datafile &
    )
    IF (ASSOCIATED(vis)) &
       CALL SetAttr(config, "sources/vis", vis)
    IF (ASSOCIATED(rotframe)) &
       CALL SetAttr(config, "sources/rotframe", rotframe)
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    USE sources_base_mod, ONLY : sources_base
    USE sources_rotframe_mod, ONLY : sources_rotframe
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)  :: Physics
    CLASS(mesh_base)     :: Mesh
    CLASS(timedisc_base) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CLASS(sources_base), POINTER :: sp => null()
    INTEGER              :: k,n,clock
    INTEGER, DIMENSION(:), ALLOCATABLE    :: seed
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics
    INTENT(INOUT)        :: Timedisc
    !------------------------------------------------------------------------!
    ! Seed the random number generator
    ! initialize with a mix from current time and mpi rank
    CALL SYSTEM_CLOCK(COUNT=clock)
    CALL RANDOM_SEED(SIZE = n) ! determine the minimum size of the seed
    seed = clock + (Timedisc%getrank()+1) * (/(k-1, k=1,n)/)
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(dv)

    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_euler) ! non-isothermal HD
      ! flow along x-direction:
      WHERE ((Mesh%bcenter(:,:,:,1).LT.(Mesh%ymin+0.33*PI)).OR. &
        (Mesh%bcenter(:,:,:,1).GT.(SIM%Mesh%ymin+0.66*PI)))
        pvar%density%data3d(:,:,:)    = RHO0
        pvar%velocity%data4d(:,:,:,2) = V0
        pvar%velocity%data4d(:,:,:,1) = 0.0
        pvar%pressure%data3d(:,:,:)   = P0
      ELSEWHERE
        pvar%density%data3d(:,:,:)    = RHO1
        pvar%velocity%data4d(:,:,:,2) = V1
        pvar%velocity%data4d(:,:,:,1) = 0.0
        pvar%pressure%data3d(:,:,:)   = P1
      END WHERE
      sp => Sim%Sources%GetSourcesPointer(ROTATING_FRAME)
      IF (ASSOCIATED(sp)) THEN
        SELECT TYPE(sp)
        CLASS IS(sources_rotframe)
          CALL sp%Convert2RotatingFrame(Mesh,Physics,pvar)
        END SELECT
      END IF

      SELECT CASE(Mesh%fargo%GetType())
      CASE(1,2)
        ! 2D simulation -> 2nd velocity is the azimuthal velocity
        Timedisc%w(:,:) = pvar%velocity%data4d(:,Mesh%JMIN,:,2)
      END SELECT

      ! add velocity perturbations
      pvar%velocity%data4d(:,:,:,1) = pvar%velocity%data4d(:,:,:,1) &
            + (dv(:,:,:,1)-0.5)*0.2
      pvar%velocity%data4d(:,:,:,2) = pvar%velocity%data4d(:,:,:,2) &
            + (dv(:,:,:,2)-0.5)*0.2
    CLASS DEFAULT
      CALL Physics%Error("KHI2D_sphere::InitData","only non-isothermal HD supported")
    END SELECT

    ! initial condition
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "spherical shear flow")
  END SUBROUTINE InitData

END PROGRAM KHI
