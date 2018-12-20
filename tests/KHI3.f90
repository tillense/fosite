!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: KHI.f90                                                           #
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
!> Kevin-Helmholtz instability
!! \author Tobias Illenseer
!!
!! References:
!! -# Lord Kelvin (William Thomson) (1871). "Hydrokinetic solutions and
!!     observations" Philosophical Magazine 42: 362-377
!! -# Hermann von Helmholtz (1868). "On the discontinuous movements of fluids"
!!     Monthly Reports of the Royal Prussian Academy of Philosophy in Berlin 23: 215-228
!----------------------------------------------------------------------------!
PROGRAM KHI
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE sources_generic
  USE fileio_generic
  USE timedisc_generic
  USE common_dict
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
  REAL, PARAMETER    :: XYLEN= 1.0         ! spatial extend
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'KHI'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  INTEGER            :: n
  REAL               :: sigma
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: pvar
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ! allocate memory for random seed variable
  CALL RANDOM_SEED(size=n)
  ALLOCATE(seed(n))

  ! test flow along x-direction
  CALL InitFosite(Sim)
  CALL MakeConfig(Sim%config)
  CALL SetupFosite(Sim)
  ! allocate memory to store the result of the first run
  ALLOCATE(pvar(Sim%Mesh%IGMIN:Sim%Mesh%IGMAX,Sim%Mesh%JGMIN:Sim%Mesh%JGMAX,Sim%Physics%VNUM))
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,.FALSE.,.FALSE.)
  CALL RunFosite(Sim)
  ! transpose result and store for comparison with 2nd run
  pvar(:,:,Sim%Physics%DENSITY) = TRANSPOSE(Sim%Timedisc%pvar(:,:,Sim%Physics%DENSITY))
  pvar(:,:,Sim%Physics%XVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar(:,:,Sim%Physics%YVELOCITY))
  pvar(:,:,Sim%Physics%YVELOCITY) = TRANSPOSE(Sim%Timedisc%pvar(:,:,Sim%Physics%XVELOCITY))
  pvar(:,:,Sim%Physics%PRESSURE) = TRANSPOSE(Sim%Timedisc%pvar(:,:,Sim%Physics%PRESSURE))

  ! test flow along y-direction
  CALL InitFosite(Sim)
  CALL MakeConfig(Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_rotate"))
  CALL SetupFosite(Sim)
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc,.TRUE.,.TRUE.)
  CALL RunFosite(Sim)
  ! compare results
  sigma = SQRT(SUM((Sim%Timedisc%pvar(:,:,:)-pvar(:,:,:))**2)/SIZE(pvar))
  CALL CloseFosite(Sim)

  DEALLOCATE(pvar,seed)

  TAP_CHECK_SMALL(sigma,TINY(sigma),"x-y symmetry test")
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               sources, timedisc, fluxes, vis
    REAL     :: dynvis
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
               "inum" / RES, &                 ! resolution in x and         !
               "jnum" / RES, &                 !   y direction               !             
               "xmin" / (-0.5*XYLEN), &
               "xmax" / (0.5*XYLEN), &
               "ymin" / (-0.5*XYLEN), &
               "ymax" / (0.5*XYLEN), &
               "output/dl" / 0, &
               "output/bh" / 0, &
               "output/rotation" / 0 &
    )

    ! physics settings
    physics => Dict( &
              "problem" / EULER2D, &
              "gamma"   / GAMMA &            ! ratio of specific heats        !
    )

    ! flux calculation and reconstruction method
    fluxes => Dict( &
             "fluxtype"  / KT, &
             "order"     / LINEAR, &
             "variables" / CONSERVATIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &       ! one of: minmod, monocent,...   !
             "theta"     / 1.2 &             ! optional parameter for limiter !
    )

    ! boundary conditions
    boundary => Dict( &
               "western"  / PERIODIC, &
               "eastern"  / PERIODIC, &
               "southern" / PERIODIC, &
               "northern" / PERIODIC &
    )

    NULLIFY(sources)
    ! viscosity source term
    ! compute dynamic viscosity constant using typical scales and Reynolds number
    dynvis = ABS(RHO0 * XYLEN * (V0-V1) / RE)
    IF (dynvis.GT.TINY(1.0)) THEN
       sources => Dict( &
          "vis/stype"     / VISCOSITY, &
          "vis/vismodel"  / MOLECULAR, &
          "vis/dynconst"  / dynvis, &
          "vis/bulkconst" / (-2./3.*dynvis), &
          "vis/output/dynvis" / 0, &
          "vis/output/stress" / 0, &
          "vis/output/kinvis" / 0, &
          "vis/output/bulkvis" / 0 &
       )
    END IF

    ! time discretization settings
    timedisc => Dict( &
               "method"   / MODIFIED_EULER, &
               "order"    / 3, &
               "cfl"      / 0.4, &
               "stoptime" / TSIM, &
               "dtlimit"  / 1.0E-10, &
               "maxiter"  / 100000, &
               "output/pressure" / 1, &
               "output/density" / 1, &
               "output/xvelocity" / 1, &
               "output/yvelocity" / 1 &
    )

    ! initialize data input/output
    datafile => Dict(&
!        "fileformat" / HDF, &
!        "fileformat" / VTK, &
        "fileformat" / GNUPLOT, "filecycles" / 0, &
!        "fileformat" / BINARY, & 
!        "fileformat" / NETCDF, &
        "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
        "count"      / ONUM &
    )

    config => Dict( &
             "mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile &
    )

    IF (ASSOCIATED(sources)) &
       CALL SetAttr(config, "sources", sources)

  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc,rotate90deg,reuse_random_seed)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    LOGICAL, OPTIONAL :: rotate90deg,reuse_random_seed
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    INTEGER           :: clock
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,rotate90deg,reuse_random_seed
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! Seed the random number generator
    IF (.NOT.(PRESENT(reuse_random_seed).AND.reuse_random_seed)) THEN
       ! initialize with a mix from current time and mpi rank
       CALL SYSTEM_CLOCK(COUNT=clock)
       seed = clock + (GetRank(Timedisc)+1) * (/(i-1, i=1,n)/)
    END IF
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(dv)

    ! initial condition
    IF (PRESENT(rotate90deg).AND.rotate90deg) THEN
       ! flow along y-direction:
       ! x-velocity vanishes everywhere
       Timedisc%pvar(:,:,Physics%XVELOCITY) = 0.
       WHERE ((Mesh%bcenter(:,:,1).LT.(Mesh%xmin+0.25*XYLEN)).OR. &
            (Mesh%bcenter(:,:,1).GT.(Mesh%xmin+0.75*XYLEN)))
          Timedisc%pvar(:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar(:,:,Physics%YVELOCITY) = V0
          Timedisc%pvar(:,:,Physics%PRESSURE) = P0
       ELSEWHERE
          Timedisc%pvar(:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar(:,:,Physics%YVELOCITY) = V1
          Timedisc%pvar(:,:,Physics%PRESSURE) = P1
       END WHERE
       ! add perturbation to velocity field
       Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
            + (TRANSPOSE(dv(:,:,2))-0.5)*0.02
       Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
            + (TRANSPOSE(dv(:,:,1))-0.5)*0.02
    ELSE   
       ! flow along x-direction:
       ! y-velocity vanishes everywhere
       Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.
       WHERE ((Mesh%bcenter(:,:,2).LT.(Mesh%ymin+0.25*XYLEN)).OR. &
            (Mesh%bcenter(:,:,2).GT.(Mesh%ymin+0.75*XYLEN)))
          Timedisc%pvar(:,:,Physics%DENSITY) = RHO0
          Timedisc%pvar(:,:,Physics%XVELOCITY) = V0
          Timedisc%pvar(:,:,Physics%PRESSURE) = P0
       ELSEWHERE
          Timedisc%pvar(:,:,Physics%DENSITY) = RHO1
          Timedisc%pvar(:,:,Physics%XVELOCITY) = V1
          Timedisc%pvar(:,:,Physics%PRESSURE) = P1
       END WHERE
       ! add perturbation to velocity field
       Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
            + (dv(:,:,1)-0.5)*0.02
       Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
            + (dv(:,:,2)-0.5)*0.02
    END IF

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability")

  END SUBROUTINE InitData
  
END PROGRAM KHI
