!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: vortex2d.f90                                                      #
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
  USE sources_rotframe, ONLY : convert2RotatingFrame_rotframe
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! general constants
  REAL, PARAMETER    :: GN      = 6.6742D-11        ! Newtons grav. constant       !
  REAL, PARAMETER    :: CC      = 2.99792458D+8     ! speed of light               !
  REAL, PARAMETER    :: RG      = 8.31447           ! molar gas constant           !
  REAL, PARAMETER    :: AU      = 1.49597870691E+11 ! astronomical unit [m]        !
  REAL, PARAMETER    :: YEAR    = 3.15576E+7        ! Julian year [sec]            !
  REAL, PARAMETER    :: DAY     = 8.6400E+4         ! Day [sec]                    !
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 0.2*DAY          ! simulation stop time         !
  INTEGER, PARAMETER :: XRES    = 60               ! theta-resolution             !
  INTEGER, PARAMETER :: YRES    = 60               ! phi-resolution               !
  ! initial condition
  ! geometrical parameter
  REAL, PARAMETER    :: GPAR    = 6.371E+6          ! geometry scaling parameter   !
  REAL, PARAMETER    :: THETA0  = 0.41              ! axis-plane-angle             !
  REAL, PARAMETER    :: PHI0    = 0.0               !                              !
  REAL, PARAMETER    :: OMEGA   = -2.99E-7          ! ang. sp. of rot. frame[rad/s]!
  REAL, PARAMETER    :: OMEGASUN= -6.229E-7         ! frequency(!) of day-night[/s]!
  REAL, PARAMETER    :: RD      = 1*AU              ! distance between star and    !
                                                    ! planet in AU                 !
  REAL, PARAMETER    :: PLANET_YEAR = YEAR          ! tropical year of the planet  !
  ! gaseous parameter
  REAL, PARAMETER    :: GAMMA   = 1.4               ! ratio of specific heats      !
  REAL, PARAMETER    :: P0      = 1.014E+5          ! ambient surf. press.[N/m**3] !
  REAL, PARAMETER    :: MU      = 2.897E-2          ! molar mass of atmosphere     !
  REAL, PARAMETER    :: c_p     = 1.005E+3*MU       ! spec. heat cap. in J/(mol*K) !
  REAL, PARAMETER    :: T_0     = 287.76            ! med. init. temp. [K]         !
  REAL, PARAMETER    :: RHO0    = P0*MU/(T_0*RG)    ! ambient surf. dens.[kg/m**3] !
  ! other atmospherical parameter
  REAL, PARAMETER    :: ALBEDO  = 0.3               ! albedo of the Planet         !
  REAL, PARAMETER    :: INTENSITY= 1.4E+3           ! intesity of the sun at 1 AU  !
  REAL, PARAMETER    :: HEIGHT  = 2.0E+4            ! height of the atmosphere     !
  ! other parameters                                ! opacity: dt/dp=-1/g*kappa    !
  REAL, PARAMETER    :: GACC    = 9.81              ! grav. accel. [m/s**2]        ! 
  ! output parameters
  INTEGER, PARAMETER :: ONUM    = 10                ! number of output data sets  !
  CHARACTER(LEN=256), PARAMETER &                   ! output data dir              !
                     :: ODIR    = './'
  CHARACTER(LEN=256), PARAMETER &                   ! output data file name        !
                     :: OFNAME  = 'planet2d' 
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!



 !   open(11, file="energy.dat",status="new",action="write")


  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

  

!    close(unit=11)
  
CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fosite_TYP)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(4)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes, heating, sources,&
                               cooling, rotframe
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    !mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
         "geometry" / BIANGLESPHERICAL, &
         "omega"    / OMEGA, &
         "inum"     / XRES, &
         "jnum"     / YRES, &
         "xmin"     / (0.001), &
         "xmax"     / (PI-0.001), &
         "ymin"     / (0.0), &
         "ymax"     / (2.0*PI), &
         "gparam"   / GPAR, &
         "dz"       / HEIGHT, &
         "output/rotation" / 0, &
         "output/volume"   / 1, &
         "output/dAz"      / 1)

    ! boundary conditions
    boundary => Dict( &
         "western" / REFLECTING, &
         "eastern" / REFLECTING, &
         "southern" / PERIODIC, &
         "northern" / PERIODIC)

    ! physics settings
    physics => Dict( &
         "problem" / EULER2D, &
         "gamma"   / GAMMA, &                       ! ratio of specific heats        !
         "dpmax"   / 1.0)                           ! for advanced time step control !

    ! flux calculation and reconstruction method
    fluxes => Dict( &
         "fluxtype"  / KT, &
         "order"     / LINEAR, &
         "variables" / CONSERVATIVE, &              ! vars. to use for reconstruction!
!          "variables" / PRIMITIVE, &               ! vars. to use for reconstruction!
         "limiter"   / MONOCENT, &                  ! one of: minmod, monocent,...   !
         "theta"     / 1.2)                         ! optional parameter for limiter !

    rotframe => Dict( &
         "stype" / ROTATING_FRAME, &
         "gparam"   / GPAR, &
         "x"     / 0.0, &
         "y"     / 0.0)

    cooling => Dict( &
         "stype" / PLANET_COOLING, &
         "output/Qcool" / 1, &
         "output/T_s" / 1, &
         "output/RHO_s" / 1, &   
         "distance" / RD, &
         "output/P_s" / 1, &   
         "albedo" / ALBEDO, &
         "mu" / MU, &
         "intensity" / INTENSITY, &
         "dz" / HEIGHT, &
         "T_0" / T_0, &
         "cvis" / 0.1, &
         "c_p" / c_p, &
         "gacc" / GACC, &
         "gamma" / GAMMA) 
          
    heating => Dict( &
         "stype" / PLANET_HEATING, &
         "output/Qstar" / 1, &                   ! Qstar is the heating term 
         "distance"/ RD, &
         "year" / PLANET_YEAR, &
         "theta0" / THETA0, &
         "phi0" / PHI0, &
         "omegasun"/ OMEGASUN, &
!          "R_planet"/ GPAR, &
         "albedo"/ ALBEDO, &
         "mu" / MU, &
         "intensity"/ INTENSITY,&
         "dz"      / HEIGHT, &
!          "a_eff"/5.35e-7, &
         "cvis"  / 0.1,&
         "c_p"  / c_p,&
         "gacc" / GACC, &
         "gamma" / GAMMA)

    sources => Dict( &
         "rotframe" / rotframe, &
         "cooling"  / cooling, &
         "heating"  / heating)

    ! time discretization settings
    timedisc => Dict( &
         "method"   / MODIFIED_EULER, &
         "order"    / 3, &
         "cfl"      / 0.4, &
         "stoptime" / TSIM, &
         "dtlimit"  / 1.0E-15, &
         "maxiter"  / 1000000)

    datafile => Dict(&
!          "fileformat" / BINARY, "filecycles" / 0, &
          "fileformat" / GNUPLOT, "filecycles" / 0, &
!          "fileformat" / HDF, &
!          "fileformat" / XDMF, &
         "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
         "count"      / ONUM)

    config => Dict( &
         "mesh" / mesh, &
         "physics"  / physics, &
         "boundary" / boundary, &
         "fluxes"   / fluxes, &
         "timedisc" / timedisc, &
         "sources"  / sources, &
!          "logfile"  / logfile, &
         "datafile" / datafile)
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
    REAL              :: theta1,phi1,theta2,phi2,vtheta,vphi,xlen,SIGMA
    REAL              :: R0,P1,PHISTART,THETASTART
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: dv
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! some starting constants
    PHISTART   = 0.0
    THETASTART = PI/3.
    P1 = P0-40.0E+2
    R0 = 0.4*GPAR
    SIGMA   = 1.0E+5 

    ! different starting condition
    Timedisc%pvar(:,:,Physics%DENSITY)   = P0/GACC
    Timedisc%pvar(:,:,Physics%PRESSURE)  = Physics%Constants%RG*&
      Timedisc%pvar(:,:,Physics%DENSITY)*T_0/(MU*(1.+(GAMMA-1.)/GAMMA)) 
    ! 1. random velocity distribution
!    CALL RANDOM_NUMBER (dv)
!    Timedisc%pvar(:,:,Physics%XVELOCITY) = Timedisc%pvar(:,:,Physics%XVELOCITY) &
!         + (dv(:,:,1)-0.5)*2.
!    Timedisc%pvar(:,:,Physics%YVELOCITY) = Timedisc%pvar(:,:,Physics%YVELOCITY) &
!         + (dv(:,:,2)-0.5)*2.

     ! 2. test case for rotframe 
     !    (velocity strips along latitude)
     xlen = ABS(Mesh%xmax-Mesh%xmin)
     WHERE ((Mesh%bcenter(:,:,1).GT.(Mesh%xmin+0.37*xlen)).AND. &
          (Mesh%bcenter(:,:,1).LT.(Mesh%xmin+0.4*xlen)))
        Timedisc%pvar(:,:,Physics%YVELOCITY) = -80.0
     ELSEWHERE ((Mesh%bcenter(:,:,1).LT.(Mesh%xmin+(PI-0.37*xlen))).AND. &
          (Mesh%bcenter(:,:,1).GT.(Mesh%xmin+(PI-0.4*xlen))))
        Timedisc%pvar(:,:,Physics%YVELOCITY) = +80.0
     ELSEWHERE ((Mesh%bcenter(:,:,1).LT.(Mesh%xmin+(PI-0.2*xlen))).AND. &
          (Mesh%bcenter(:,:,1).GT.(Mesh%xmin+(PI-0.23*xlen))))
        Timedisc%pvar(:,:,Physics%YVELOCITY) = -80.0
     ELSEWHERE ((Mesh%bcenter(:,:,1).GT.(Mesh%xmin+0.2*xlen)).AND. &
          (Mesh%bcenter(:,:,1).LT.(Mesh%xmin+0.23*xlen)))
        Timedisc%pvar(:,:,Physics%YVELOCITY) = +80.0
     END WHERE

!    ! 3. low density area
!     WHERE ((Mesh%bcenter(:,:,1).LE.(THETASTART+R0/GPAR)).AND.&
!           (Mesh%bcenter(:,:,2).GE.(PHISTART-R0/GPAR/SIN(THETASTART))).AND.&
!           (Mesh%bcenter(:,:,1).GE.(THETASTART-R0/GPAR)).AND.& 
!           (Mesh%bcenter(:,:,2).LE.(PHISTART+R0/GPAR/SIN(THETASTART)))) 
!        ! behind the shock front
!        Timedisc%pvar(:,:,Physics%PRESSURE)  = P1
!     ELSEWHERE
!       Timedisc%pvar(:,:,Physics%PRESSURE)  = SIGMA*&
!               Timedisc%pvar(:,:,Physics%DENSITY)**GAMMA
!       
!     END WHERE
    
    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 2D planetary")

  END SUBROUTINE InitData

END PROGRAM Init
