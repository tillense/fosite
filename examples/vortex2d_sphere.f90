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

!----------------------------------------------------------------------------!
!> Program and data initialization for 2D isentropic vortex
!! References:
!! [1] Yee, H. C. et al.: Low-dissipative high-order shock-capturing methods
!!     using characteristic-based filters, J. Comput. Phys. 150 (1999), 199-238
!!     DOI: 10.1006/jcph.1998.6177
!!
!! [2] Ramachandran D. Nair, Christiane Jablonowski:
!!     "Moving Vortices on the Sphere: A Test Case for Horizontal Advection Problems"
!!     www-personal.umich.edu/~cjablono/MWR_Nair_Jablonowski_2007.pdf
!----------------------------------------------------------------------------!
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
  ! simulation parameters
  REAL, PARAMETER    :: YEAR    = 3.15576E+7        ! Julian year [sec]            !
  REAL, PARAMETER    :: TSIM    = 2.0        ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4        ! ratio of specific heats
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHOINF  = 1.0        ! ambient density
  REAL, PARAMETER    :: PINF    = 1.0        ! ambient pressure
  REAL, PARAMETER    :: VSTR    = 1.5        ! nondimensional vortex strength
  REAL, PARAMETER    :: LAMBDA  = 0.12       ! vortex strength between 0-1
  REAL, PARAMETER    :: UINF    = 0.0        ! cartesian components of constant
  REAL, PARAMETER    :: VINF    = 0.0        ! global velocity field
  REAL, PARAMETER    :: X0      = 0.0        ! vortex position (cart. coords.)
  REAL, PARAMETER    :: Y0      = 0.0      
  REAL, PARAMETER    :: THETA0  = PI/2.      ! vortex position (spheric. coords.)
  REAL, PARAMETER    :: PHI0    = PI/2. 
  REAL, PARAMETER    :: OMEGA   = 0.0        ! angular speed of rotational frame
                                             ! around [X0,Y0]
  ! mesh settings
!   INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry
!$$  INTEGER, PARAMETER :: MGEO = POLAR    
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
!!$  INTEGER, PARAMETER :: MGEO = BIPOLAR
 INTEGER, PARAMETER :: MGEO = BIANGLESPHERICAL

  INTEGER, PARAMETER :: XRES    = 100        ! x-resolution
  INTEGER, PARAMETER :: YRES    = 100        ! y-resolution
  REAL, PARAMETER    :: RMIN    = 1.0E-2     ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX    = 5.0        ! outer radius
  REAL, PARAMETER    :: GPAR    = 1.0        ! geometry scaling parameter 
  REAL, PARAMETER    :: R_0     = GPAR/5.    ! important for extention of the vortex
                                             ! on a sphere
  REAL, PARAMETER    :: SHEIGHT = 1.0        ! scale height
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10            ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &            ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &            ! output data file name
                     :: OFNAME = 'vortex2d_sphere'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc)
  
  CALL RunFosite(Sim)

  CALL CloseFosite(Sim)

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
                               timedisc, fluxes, sources, rotframe
    REAL              :: x1,x2,y1,y2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 =-RMAX*PI/2.
       x2 = RMAX*PI/2.
       y1 =-RMAX*PI 
       y2 = RMAX*PI
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(SINHPOLAR)
       x1 = RMIN/GPAR
       x1 = LOG(x1+SQRT(1.0+x1*x1))  ! = ASINH(RMIN/GPAR))
       x2 = RMAX/GPAR
       x2 = LOG(x2+SQRT(1.0+x2*x2))  ! = ASINH(RMAX/GPAR))
       y1 = 0.0 
       y2 = 2.0*PI       
    CASE(BIPOLAR)
       x1 = 0.0
       x2 = 2.0*PI
       y1 = GPAR/RMIN
       y1 = -LOG(y1+SQRT(1.0+y1*y1))  ! = -ASINH(GPAR/RMIN)
       y2 = -0.5*y1                     ! should be larger than |y1|
    CASE(BIANGLESPHERICAL)
!       x1 = 0.001
       x1 = 0.3
!       x2 = PI/2.
       x2 = PI-0.3
       y1 = 0.0
       y2 = PI
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for 2D isentropic vortex")
    END SELECT

    !mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "gparam"   / GPAR, &
           "dz"       / SHEIGHT, &
           "output/rotation" / 0)


    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       bc(WEST)  = PERIODIC
       bc(EAST)  = PERIODIC
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(POLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(LOGPOLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(TANPOLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(SINHPOLAR)
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE(BIPOLAR)
       bc(WEST)  = PERIODIC
       bc(EAST)  = PERIODIC
       bc(SOUTH) = REFLECTING
       bc(NORTH) = REFLECTING
    CASE(BIANGLESPHERICAL)
       bc(WEST)  = REFLECTING
       bc(EAST)  = REFLECTING
       bc(SOUTH) = PERIODIC
       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Error(Sim%Physics,"InitProgram","mesh geometry not supported for 2D isentropic vortex")
    END SELECT

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH))

    ! physics settings
    physics => Dict("problem" / EULER2D, &
              "gamma"   / GAMMA, &                 ! ratio of specific heats        !
              "dpmax"   / 1.0)                    ! for advanced time step control !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
              "fluxtype"  / KT, &
!             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
              "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
              "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
              "theta"     / 1.2)          ! optional parameter for limiter !

!    rotframe => Dict("stype" / ROTATING_FRAME, &
!               "omega" / OMEGA, &
!               "gparam"   / GPAR, &
!               "x"     / X0, &
!               "y"     / Y0)

!    sources => Dict("rotframe" / rotframe)

    ! time discretization settings
    timedisc => Dict(&
          "method"   / MODIFIED_EULER, &
          "order"    / 3, &
          "cfl"      / 0.4, &
          "stoptime" / TSIM, &
          "dtlimit"  / 1.0E-14, &
          "maxiter"  / 100000000)

    ! initialize log input/output
!!$    CALL InitFileIO(Logfile,Mesh,Physics,Timedisc,&
!!$         fileformat = BINARY,, &
!!$         filename   = TRIM(ODIR) // TRIM(OFNAME) // 'log',, &
!!$         filecycles = 1)

    ! initialize data input/output
!    datafile => Dict("fileformat" / VTK, &
!    datafile => Dict("fileformat" / BINARY, &
     datafile => Dict(&
!          "fileformat" / GNUPLOT, "filecycles" / 0, &
!           "fileformat" / HDF, &
           "fileformat" / XDMF, &
          "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
          "count"      / ONUM)

    config => Dict(&
          "mesh"       / mesh, &
          "physics"    / physics, &
          "boundary"   / boundary, &
          "fluxes"     / fluxes, &
          "timedisc"   / timedisc, &
!           "sources"   / sources, &
!           "logfile"   / logfile, &
          "datafile"   / datafile)
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: vxy
    REAL, DIMENSION(3,3,3) :: rotx,rotz
    REAL              :: theta1,phi1,theta2,phi2,x,y,r2
    REAL              :: du,dv,vtheta2,vphi2,vphi,vtheta
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! initial condition
!    open(11, file="theta2.dat",status="new",action="write")
    IF (MGEO.EQ.BIANGLESPHERICAL) THEN
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          ! standard curves
          theta1 = mesh%bcenter(i,j,1)
          phi1   = mesh%bcenter(i,j,2)

          ! transformed curves 
          theta2 = ACOS(COS(theta1)*COS(THETA0) + &   
                   SIN(theta1)*SIN(THETA0)*COS(phi1-PHI0))  
          phi2   = MODULO(ATAN2(SIN(theta1)*SIN(phi1-PHI0),& 
                   (-COS(theta1)*SIN(THETA0)+ &
                   SIN(theta1)*COS(THETA0)*COS(phi1-PHI0))),2*PI)
                   
!          theta2 = ACOS(COS(theta1)*COS(THETA0) - &
!                   SIN(theta1)*SIN(THETA0)*COS(phi1))
!          phi2   = MODULO(PHI0 + ATAN2(SIN(theta1)*SIN(phi1), &
!                   SIN(theta1)*COS(phi1)*COS(theta1)- &
!                   COS(theta1)*SIN(THETA0)),2*PI)

          ! velocities for pole-vortex
          IF ((theta2) .LT. (PI/3.)) THEN
          vtheta2 = 0.0
          vphi2   = (GPAR/R_0)*SQRT(GAMMA/(GAMMA-1.)*LAMBDA*2.*theta2&
                    *TAN(theta2))*EXP(0.5*(1.-((GPAR*theta2)/R_0)**(2.))) 
!          vphi2 = VSTR*GPAR*((theta2)*TAN(theta2))**(1./2.)*&
!          EXP(-(GPAR**(2.)*theta2**(2.)/2.))
          ELSE
          vtheta2 = 0.0
          vphi2   = 0.0
          END IF       

          ! transformed velocities for arbitrary vortex
!           vtheta   = vphi2*SIN(THETA0)*SIN(phi2)
!           vphi     = vphi2*(COS(THETA0)*SIN(theta2)+ &
!                      SIN(THETA0)*COS(theta2)*COS(phi2))

!          vtheta   = -vphi2*SIN(THETA0)*SIN(phi1-PHI0)
          vtheta   = -vphi2*SIN(THETA0)*SIN(phi2)/SIN(theta1)
!          vphi     = vphi2*(SIN(theta1)*COS(THETA0)-&
!                   COS(theta1)*SIN(THETA0)*COS(phi1-PHI0)) 
          vphi     = (vphi2/SIN(theta2))&!(SIN(phi2)/(SIN(theta1)*SIN(phi1-PHI0)))&
                     *(SIN(theta1)*COS(THETA0)-&
                     COS(theta1)*SIN(THETA0)*COS(phi1-PHI0)) 

          ! initial conditions
          Timedisc%pvar(i,j,Physics%XVELOCITY) = vtheta
          Timedisc%pvar(i,j,Physics%YVELOCITY) = vphi + VINF*SIN(theta1)
!          Timedisc%pvar(i,j,Physics%DENSITY)   = RHOINF * (1.0 - &
!                        (GAMMA-1.0)/(GAMMA*2.0)*VSTR**2.* &
!                        EXP((GPAR**(2.0)*theta2**(2.0))) )**(1./(GAMMA-1.))
          Timedisc%pvar(i,j,Physics%DENSITY)   = RHOINF*(1.0-LAMBDA*&
                         EXP(1.0-((GPAR*theta2)/R_0)**(2.)))**(1./(GAMMA-1.))
          Timedisc%pvar(i,j,Physics%PRESSURE)  = PINF &
                        *(Timedisc%pvar(i,j,Physics%DENSITY)/RHOINF)**GAMMA
               

          ! output to get vortex time
          IF ((Timedisc%pvar(i,j,Physics%DENSITY).GE.(0.9-0.02)*RHOINF).AND.&
             (Timedisc%pvar(i,j,Physics%DENSITY).LE.(0.9+0.02)*RHOINF)) THEN
          !Print *, theta2,vphi2, SQRT(vphi**2.+vtheta**2.),"T=",2.*PI*GPAR*theta2 &
          !                                           / SQRT(vphi**2.+vtheta**2.)
          ELSE
          ! do nothing
          END IF

       END DO
    END DO
    ELSE
    DO j=Mesh%JGMIN,Mesh%JGMAX
       DO i=Mesh%IGMIN,Mesh%IGMAX
          x  = Mesh%bccart(i,j,1) - X0
          y  = Mesh%bccart(i,j,2) - Y0
          r2 = x*x + y*y
          du = -0.5*VSTR/PI*y*EXP(0.5*(1.-r2))
          dv = -du*x/y
          ! cartesian velocity components
          vxy(i,j,1) = UINF + du
          vxy(i,j,2) = VINF + dv
          ! density
          ! ATTENTION: there's a factor of 1/PI missing in the density
          ! formula  eq. (3.3) in [1]
          Timedisc%pvar(i,j,Physics%DENSITY) = RHOINF * (1.0 - &
               (GAMMA-1.0)*VSTR**2/(8*PI**2*GAMMA) * EXP(1.-r2) )**(1./(GAMMA-1.))
          IF ((Timedisc%pvar(i,j,Physics%DENSITY).GE.(0.9-0.2)*RHOINF).AND.&
             (Timedisc%pvar(i,j,Physics%DENSITY).LE.(0.9+0.2)*RHOINF)) THEN
          Print *, x, y, r2, (vxy(i,j,1)**2.+vxy(i,j,2)**2.)**(1./2.)
          ELSE
          ! do nothing
          END IF
          ! pressure
          Timedisc%pvar(i,j,Physics%PRESSURE) = PINF &
               * (Timedisc%pvar(i,j,Physics%DENSITY)/RHOINF)**GAMMA
       END DO
    END DO
    ! transform velocities
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,vxy(:,:,:),&
         Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY))
    END IF
!    write(11,*)
!    write(11,*)
!    close(unit=11)

    !CALL Convert2RotatingFrame_rotframe(&
    !    GetSourcesPointer(Physics%sources, ROTATING_FRAME),&
    !    Mesh,&
    !    Physics,&
    !    Timedisc%pvar)

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: 2D isentropic vortex")

  END SUBROUTINE InitData
END PROGRAM Init
