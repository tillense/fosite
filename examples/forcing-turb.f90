!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: forcing-turb.f90                                                  #
!#                                                                           #
!# Copyright (C) 2012                                                        #
!# Björn Sperling <sperling@astrophysik.uni-kiel.de>                         #
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
!> Program and data initialization for forcing turbulence
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
#ifdef HAVE_FFTW
  USE fftw
#endif
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM  =  1.0E+0       ! simulation time
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  REAL, PARAMETER    :: ETA   = 0.0E-3     ! dynamic viscosity (0.0 disables)
  REAL, PARAMETER    :: P0    = 5.0
  REAL, PARAMETER    :: XI    = 1.0E-9     ! ratio of Psgs/P           
  REAL, PARAMETER    :: L    = 1.0E-1         ! length scale
  REAL, PARAMETER    :: V    = 4.0E+1      ! characteristic velocity
  REAL, PARAMETER    :: CHI  = 0.75         ! 
  REAL, PARAMETER    :: Lx   = 1.0
  REAL, PARAMETER    :: Ly   = Lx
  REAL, PARAMETER    :: fSTOP = TSIM * 0.1
  INTEGER, PARAMETER :: PROBLEM = EULER2D_SGS
!  INTEGER, PARAMETER :: PROBLEM = EULER2D
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
  INTEGER, PARAMETER :: XRES = 50 !512         ! resolution
  INTEGER, PARAMETER :: YRES = 50 !512
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 200          ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = ''
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = "spectest" 
!'forcing2D_L1.0E-1_CHI0.75_128x128_novis_SGS' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)
 
  CALL MakeConfig(Sim%config)

  CALL SetupFosite(Sim)

!  CALL PrintDict(Sim%IO)

!  CALL PrintDict(Sim%config)

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
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, &
                               sources, timedisc, fluxes, stemp
    !------------------------------------------------------------------------!
    ! mesh settings
    mesh => Dict( &
           "meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
               "inum" / XRES, &                 ! resolution in x and         !
               "jnum" / YRES, &                 !   y direction               !             
               "xmin" / 0.0, &
               "xmax" / Lx, &
               "ymin" / 0.0, &
               "ymax" / Ly &      
    )

    ! physics settings
    physics => Dict( &
              "problem" / PROBLEM, &
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


    ! sgs source term
    stemp => Dict( &
            "stype"     / SGS, &
            "output/stress" / 1, &
            "output/dynvis" / 1)

    NULLIFY(sources)
    IF (PROBLEM.EQ.EULER2D_SGS) &
      CALL SetAttr(sources, "sgs", stemp)

    ! viscosity source term
    stemp => Dict( &
          "stype"     / VISCOSITY, &
          "vismodel"  / MOLECULAR, &
          "dynconst"  / ETA, &
          "bulkconst" / (-2./3.*ETA), &
          "output/dynvis" / 0, &
          "output/stress" / 0, &
          "output/kinvis" / 0, &
          "output/bulkvis" / 0 &
    )

    IF (ETA.GT.TINY(ETA)) &
        CALL SetAttr(sources, "vis", stemp)

    stemp => Dict("stype" / FORCING, &
                  "L" / L, &
                  "V" / V, &
                  "CHI" / chi, &
                  "stoptime" / fSTOP, &
                  "output/accel" / 1, &
                  "output/f" / 1, &
                  "output/fk" / 1)
    CALL SetAttr(sources, "forcing", stemp)

    ! time discretization settings
    timedisc => Dict( &
               "method"   / MODIFIED_EULER, &
               "order"    / 3, &
               "cfl"      / 0.4, &
               "stoptime" / TSIM, &
               "dtlimit"  / 1.0E-9, &
               "maxiter"  / 10000000, &
               "output/pressure" / 1, &
               "output/density" / 1, &
               "output/xvelocity" / 1, &
               "output/yvelocity" / 1, &
               "output/xmomentum" / 1, &
               "output/ymomentum" / 1)


    ! initialize data input/output
    datafile => Dict(&
!        "fileformat" / HDF,   "filecycles" / 0, &
         "fileformat" / VTK, &
!        "fileformat" / GNUPLOT,   "filecycles" / 0, &
!        "fileformat" / BINARY, & 
!	"fileformat" / NETCDF, &
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
    REAL              :: rho0, v0
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    
    ! region
    rho0 = 1.0
    v0   = 0.0

    Timedisc%pvar(:,:,Physics%DENSITY) = rho0
    Timedisc%pvar(:,:,Physics%XVELOCITY) = v0
    Timedisc%pvar(:,:,Physics%YVELOCITY) = v0
    Timedisc%pvar(:,:,Physics%PRESSURE) = P0
    IF (GetType(Physics) .EQ. Euler2D_SGS .OR. GetType(Physics) .EQ. Euler3D_ROTSYMSGS ) &
        Timedisc%pvar(:,:,Physics%SGSPRESSURE) = XI*Timedisc%pvar(:,:,Physics%PRESSURE)

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh, " DATA-----> initial condition: " // &
         "forcing turbulence")

  END SUBROUTINE InitData

END PROGRAM Init



