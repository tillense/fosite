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
!> Program and data initialization for Kelvin-Helmholtz instability
!! References: 
!! [1] Lord Kelvin (William Thomson) (1871). "Hydrokinetic solutions and 
!!     observations" Philosophical Magazine 42: 362-377 
!! [2] Hermann von Helmholtz (1868). "On the discontinuous movements of fluids"
!!     Monthly Reports of the Royal Prussian Academy of Philosophy in Berlin 23: 215-228
!----------------------------------------------------------------------------!
PROGRAM Init
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
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  = 10.0        ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  REAL, PARAMETER    :: ETA   = 0.0E-4        ! dynamic viscosity (0.0 disables)           
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
  INTEGER, PARAMETER :: XRES = 80         ! resolution
  INTEGER, PARAMETER :: YRES = 40
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 10         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'KHI2'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)
 
  CALL MakeConfig(Sim%config)

!  CALL PrintDict(Sim%config)

  CALL SetupFosite(Sim)

  ! set initial condition
  IF (Sim%Timedisc%time .EQ. 0.0) THEN
    CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  ELSE
    print *, "No InitData, time>0"
  END IF
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

CONTAINS

  SUBROUTINE MakeConfig(config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4),sgbc
    DOUBLE PRECISION, PARAMETER :: GN = 6.6742E-11     ! [m^3/kg/s^2]
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               sources, timedisc, fluxes, vis
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
               "inum" / XRES, &                 ! resolution in x and         !
               "jnum" / YRES, &                 !   y direction               !             
               "xmin" / (-0.5), &
               "xmax" / 0.5, &
               "ymin" / (-0.5), &
               "ymax" / 0.5 &
    )

    ! physics settings
    physics => Dict( &
              "problem" / EULER2D, &
              "gamma"   / GAMMA, &           ! ratio of specific heats        !
              "dpmax"   / 1.0 &              ! for advanced time step control !
    )

    ! flux calculation and reconstruction method
    fluxes => Dict( &
             "fluxtype"  / KT, &
             "order"     / LINEAR, &
             "variables" / CONSERVATIVE, &   ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &       ! one of: minmod, monocent,...   !
             "theta"     / 1.2 &            ! optional parameter for limiter !
    )

    ! boundary conditions
    boundary => Dict( &
               "western"  / PERIODIC, &
               "eastern"  / PERIODIC, &
               "southern" / PERIODIC, &
               "northern" / PERIODIC &
    )

    ! viscosity source term
    vis => Dict( &
          "stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / ETA, &
          "bulkconst" / (-2./3.*ETA), &
          "output/dynvis" / 0, &
          "output/stress" / 1, &
          "output/kinvis" / 0, &
          "output/bulkvis" / 0 &
    )

    sources => Dict("." / 0)
    IF (ETA.GT.TINY(ETA)) THEN
        sources => Dict("vis" / vis)
    END IF

    ! time discretization settings
    timedisc => Dict( &
               "method"   / MODIFIED_EULER, &
               "order"    / 3, &
               "cfl"      / 0.4, &
               "stoptime" / TSIM, &
               "dtlimit"  / 1.0E-4, &
               "maxiter"  / 100000, &
               "output/pressure" / 1, &
               "output/density" / 1, &
               "output/xvelocity" / 1, &
               "output/yvelocity" / 1 &
    )

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY,, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log',, &
!!$         filecycles = 1)

    ! initialize data input/output
    !datafile => Dict("fileformat" / VTK, &
    datafile => Dict(&
!        "fileformat" / HDF, &
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
             "sources"  / sources, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile &
    )

  END SUBROUTINE MakeConfig

  SUBROUTINE InitData(Mesh,Physics,Timedisc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Physics_TYP) :: Physics
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Timedisc_TYP):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    REAL              :: rho0, rho1, v0, v1, P0, P1
    REAL              :: ylen
    INTEGER           :: n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    
    ! outer regions
    rho0 = 1.0
    v0   = -0.5
    P0   = 2.5

    ! inner region
    rho1 = 2.0
    v1   = 0.5
    P1   = P0

    ! y-velocity vanishes everywhere
    Timedisc%pvar(:,:,Physics%YVELOCITY) = 0.

    ! extent along the the y-direction
    ylen = ABS(Mesh%ymax-Mesh%ymin)

    WHERE ((Mesh%bcenter(:,:,2).LT.(Mesh%ymin+0.25*ylen)).OR. &
         (Mesh%bcenter(:,:,2).GT.(Mesh%ymin+0.75*ylen)))
       Timedisc%pvar(:,:,Physics%DENSITY) = rho0
       Timedisc%pvar(:,:,Physics%XVELOCITY) = v0
       Timedisc%pvar(:,:,Physics%PRESSURE) = P0
    ELSEWHERE
       Timedisc%pvar(:,:,Physics%DENSITY) = rho1
       Timedisc%pvar(:,:,Physics%XVELOCITY) = v1
       Timedisc%pvar(:,:,Physics%PRESSURE) = P1
    END WHERE
       
    ! add velocity perturbations

    
    ! Seed the random number generator with a mix from current time and mpi rank
    n = 4       ! has been choosen arbitrary
    CALL RANDOM_SEED(size=n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + (GetRank(Timedisc)+1) * (/(i-1, i=1,n)/)
    CALL RANDOM_SEED(PUT=seed)
    DEALLOCATE(seed)

    CALL RANDOM_NUMBER(dv)
    Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
         + (dv(:,:,1)-0.5)*0.02
    Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
         + (dv(:,:,2)-0.5)*0.02

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "Kelvin-Helmholtz instability")

  END SUBROUTINE InitData

END PROGRAM Init
