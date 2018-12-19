!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: KHI2D.f90                                                         #
!#                                                                           #
!# Copyright (C) 2006-2018                                                   #
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
  REAL, PARAMETER    :: TSIM  = 10.0       ! simulation time
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
  INTEGER, PARAMETER :: RES  = 50          ! resolution
  REAL, PARAMETER    :: XYZLEN= 1.0        ! spatial extend
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'KHI2D'
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE   :: Sim
  CLASS(marray_compound), POINTER :: pvar,pvar_init
  INTEGER            :: n
  REAL               :: sigma
  INTEGER, DIMENSION(:), ALLOCATABLE    :: seed
  !--------------------------------------------------------------------------!
  TAP_PLAN(1)

  ! allocate memory for dynamic types and arrays
  CALL RANDOM_SEED(size=n)
  ALLOCATE(Sim,pvar,seed(n))

  ! simulate KHI with initial x-velocity along x-direction
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim,Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  ! store transposed initial data
  CALL Sim%Physics%new_statevector(pvar_init,PRIMITIVE)
  CALL TransposeData(Sim%Mesh,Sim%Timedisc%pvar,pvar_init,"xy")
  CALL Sim%Run()
  ! store transposed result of the first run
  CALL Sim%Physics%new_statevector(pvar,PRIMITIVE)
  CALL TransposeData(Sim%Mesh,Sim%Timedisc%pvar,pvar,"xy")
  ! finish the simulation
  CALL Sim%Finalize()
  DEALLOCATE(Sim)
  
  ! simulate KHI with initial y-velocity along y-direction
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_rotate"))
  CALL Sim%Setup()
  Sim%Timedisc%pvar%data1d(:) = pvar_init%data1d(:)
  CALL Sim%Physics%Convert2Conservative(Sim%Timedisc%pvar,Sim%Timedisc%cvar)
  CALL Sim%Run()
  ! compare results
  sigma = SQRT(SUM((Sim%Timedisc%pvar%data4d(:,:,:,:)-pvar%data4d(:,:,:,:))**2)/SIZE(pvar%data4d))

  CALL pvar%Destroy()
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
               "knum"            /             1, & !   z direction          !
               "xmin"            / (-0.5*XYZLEN), &
               "xmax"            /  (0.5*XYZLEN), &
               "ymin"            / (-0.5*XYZLEN), &
               "ymax"            /  (0.5*XYZLEN), &
!                "zmin"            / (-0.5*XYZLEN), &
!                "zmax"            /  (0.5*XYZLEN) &
               "zmin"            /        (-0.0), &
               "zmax"            /         (0.0) &
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
               "maxiter"          /         100000, &
               "output/pressure"  /              1, &
               "output/density"   /              1, &
               "output/xvelocity" /              1, &
               "output/yvelocity" /              1  &
    )

    ! initialize data input/output
    datafile => Dict(&
        "fileformat" /                          VTK, &
!        "fileformat" /                         XDMF, &
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
    IF (ASSOCIATED(sources)) &
       CALL SetAttr(config, "sources", sources)
  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base)  :: Physics
    CLASS(mesh_base)     :: Mesh
    CLASS(timedisc_base) :: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER              :: k,clock
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)           :: Mesh,Physics
    INTENT(INOUT)        :: Timedisc
    !------------------------------------------------------------------------!
    ! Seed the random number generator
    ! initialize with a mix from current time and mpi rank
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + (Timedisc%getrank()+1) * (/(k-1, k=1,n)/)
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(dv)
    
    ! initial condition
    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_euler) ! non-isothermal HD
      ! flow along x-direction:
      WHERE ((Mesh%bcenter(:,:,:,2).LT.(Mesh%ymin+0.25*XYZLEN)).OR. &
        (Mesh%bcenter(:,:,:,2).GT.(SIM%Mesh%ymin+0.75*XYZLEN)))
        pvar%density%data3d(:,:,:)    = RHO0
        pvar%velocity%data4d(:,:,:,1) = V0
        pvar%velocity%data4d(:,:,:,2) = 0.0 
        pvar%pressure%data3d(:,:,:)   = P0
      ELSEWHERE
        pvar%density%data3d(:,:,:)    = RHO1
        pvar%velocity%data4d(:,:,:,1) = V1
        pvar%velocity%data4d(:,:,:,2) = 0.0 
        pvar%pressure%data3d(:,:,:)   = P1
      END WHERE
      ! add perturbations
      DO k=Mesh%KGMIN,Mesh%KGMAX
          pvar%velocity%data4d(:,:,k,1) = pvar%velocity%data4d(:,:,k,1) &
                + (dv(:,:,k,1)-0.5)*0.02
          pvar%velocity%data4d(:,:,k,2) = pvar%velocity%data4d(:,:,k,2) &
                + (dv(:,:,k,2)-0.5)*0.02
      END DO
    CLASS DEFAULT
      CALL Physics%Error("KHI2D::InitData","only non-isothermal HD supported")
    END SELECT

    ! initial condition
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability")
  END SUBROUTINE InitData

  SUBROUTINE TransposeData(Mesh,pvar_in,pvar_out,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base), INTENT(IN)            :: Mesh
    CLASS(marray_compound), INTENT(INOUT)   :: pvar_in,pvar_out
    CHARACTER(LEN=2), INTENT(IN)            :: dir
    !------------------------------------------------------------------------!
    INTEGER :: k
    !------------------------------------------------------------------------!
    SELECT TYPE(pin => pvar_in)
    TYPE IS(statevector_euler) ! non-isothermal HD
      SELECT TYPE(pout => pvar_out)
      TYPE IS(statevector_euler) ! non-isothermal HD
        SELECT CASE(dir)
        CASE("xy")
          ! transpose x-y directions
          DO k=Mesh%KGMIN,Mesh%KGMAX
            ! simply transpose the 1st and 2nd indices ...
            pout%density%data3d(:,:,k) = TRANSPOSE(pin%density%data3d(:,:,k))
            pout%pressure%data3d(:,:,k) = TRANSPOSE(pin%pressure%data3d(:,:,k))
            ! ... and exchange x- and y-velocities
            pout%velocity%data4d(:,:,k,1) = TRANSPOSE(pin%velocity%data4d(:,:,k,2))
            pout%velocity%data4d(:,:,k,2) = TRANSPOSE(pin%velocity%data4d(:,:,k,1))
          END DO
        CASE DEFAULT
          CALL Mesh%Error("KHI2D::TransposeData","directions must be one of 'xy','xz' or 'yz'")
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE TransposeData
  
END PROGRAM KHI
