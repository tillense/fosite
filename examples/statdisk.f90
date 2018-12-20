!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: stationary_disk.f90                                               #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Marc Junker      <maj@astrophysik.uni-kiel.de>                            #
!# Bjoern Sperling  <sperling@astrophysik.uni-kiel.de>                       #
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
!> stationary (local) isothermal disk in centrifugal balance
!----------------------------------------------------------------------------!
PROGRAM Init
  USE fosite
  USE physics_generic
  USE fluxes_generic
  USE mesh_generic
  USE reconstruction_generic
  USE boundary_generic
  USE fileio_generic
  USE sources_generic
  USE gravity_generic
  USE timedisc_generic
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  ! the equations are solved in non-dimensional units; 
  ! we use geometrical units and set the gravitational constant to unity
  ! simulation parameters
  REAL, PARAMETER :: TSIM  = 1.0E+02       ! simulation time [dyn. time]
  ! initial conditions (power law with respect to radial distance)
  REAL, PARAMETER :: C0    = -0.0          ! power for sound speed squared
  REAL, PARAMETER :: D0    = -0.0          ! power for surface density
  REAL, PARAMETER :: CSISO = 1.0E-2        ! speed of sound @ r=1
  REAL, PARAMETER :: OMEGA = 2.            ! rotating reference frame velocity
  ! mesh settings
  INTEGER, PARAMETER :: MGEO     = POLAR      ! geometry
!!$  INTEGER, PARAMETER :: MGEO     = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO     = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO     = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 1000        ! x-resolution
  INTEGER, PARAMETER :: YRES = 1           ! y-resolution (1D!)
  REAL, PARAMETER    :: RMIN = 0.2         ! min radius of comp. domain
  REAL, PARAMETER    :: RMAX = 1.0         ! max radius of comp. domain
  REAL, PARAMETER    :: GPAR = 0.3         ! geometry scaling parameter
  ! physics settings
  LOGICAL, PARAMETER :: WITH_IAMT = .TRUE. ! with inertial ang. momemtum transport
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 100         ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'statdisk'
  !--------------------------------------------------------------------------!
  TYPE(fosite_TYP)   :: Sim
  !--------------------------------------------------------------------------!

  CALL InitFosite(Sim)

  CALL MakeConfig(Sim, Sim%config)

!  CALL PrintDict(config)

  CALL SetupFosite(Sim)

  ! set initial condition
  CALL InitData(Sim%Timedisc, Sim%Mesh, Sim%Physics, Sim%Fluxes)
  
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
    INTEGER           :: bc(4),sgbc,phys
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, logfile, sources, &
                               timedisc, fluxes, grav, vis
    REAL              :: x1,x2,om
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
    CASE(LOGPOLAR)
       x1 = LOG(RMIN/GPAR)
       x2 = LOG(RMAX/GPAR)
    CASE(TANPOLAR)
       x1 = ATAN(RMIN/GPAR)
       x2 = ATAN(RMAX/GPAR)
    CASE(SINHPOLAR)
       x1 = LOG(RMIN/GPAR+SQRT(1.0+(RMIN/GPAR)**2))  ! = ASINH(RIN))
       x2 = LOG(RMAX/GPAR+SQRT(1.0+(RMAX/GPAR)**2))
    END SELECT
    

    mesh => Dict("meshtype" / MIDPOINT, &
!!$    mesh => Dict("meshtype" / TRAPEZOIDAL, &
           "geometry" / MGEO, &
           "omega"    / om, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / (-PI), &
           "ymax"     / ( PI), &
           "gparam"   / GPAR)

    ! physics settings
    IF (WITH_IAMT) THEN
      phys = EULER2D_ISOIAMT
      om = OMEGA
    ELSE
      phys = EULER2D_ISOTHERM
      om = 0.
    END IF
    IF (ABS(C0).LT.EPSILON(C0)) THEN
      physics => Dict("problem" / phys, &
           "units"   / GEOMETRICAL, &       ! Newton constant = 1 !
           "cs"      / CSISO)               ! isothermal sound speed
    ELSE
      physics => Dict("problem" / phys, &
           "units"   / GEOMETRICAL, &       ! Newton constant = 1 !
           "output/bccsound" / 1)
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / PRIMITIVE, &   ! vars. to use for reconstruction!
             "limiter"   / VANLEER, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    boundary => Dict( &
           "western"  / CUSTOM, &
           "eastern"  / CUSTOM, &
           "southern" / PERIODIC, &
           "northern" / PERIODIC)

    ! source term due to a point mass
    grav => Dict("stype" / GRAVITY, &
            "pmass/gtype"     / POINTMASS, &                ! grav. accel. of a point mass    !
            "pmass/potential" / NEWTON, &                   ! type of gravitational potential !
            "pmass/mass"      / 1.)

    ! viscosity source term
    vis => Dict("stype"     / VISCOSITY, &
          "vismodel"  / BETA, &
          "dynconst"  / (1.0E-5), &
          "bulkconst" / 0.0, &
          "cvis"      / 0.5)


    sources => Dict("grav"       / grav) !, &
!              "vis"         / vis), &

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.5, &
           "stoptime" / TSIM, &
           "tol_rel" / 1.0E-6, &
           "tol_abs" / (/ 0.0, 1.0E-12, 1.0E-12 /), &
           "dtlimit"  / (1.0E-8*TSIM), &
           "maxiter"  / 2000000000, &
           "output/fluxes" / 1, &
           "output/external_sources" / 1, &
           "output/geometrical_sources" / 1)


    ! initialize data input/output
    datafile => Dict("fileformat" / GNUPLOT, "filecycles" / 0, "decimals" / 10, &
               "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
               "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "sources"  / sources, &
             "timedisc" / timedisc, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig


  SUBROUTINE InitData(Timedisc,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: Timedisc
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    ! Local variable declaration
    CHARACTER(LEN=32) :: info_str
    INTEGER           :: i,j,k,ig,dir
    REAL              :: r
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,4) :: fcsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM) :: accel
    REAL, DIMENSION(:,:,:), POINTER :: radius
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh
    INTENT(INOUT)     :: Timedisc,Physics,Fluxes
    !------------------------------------------------------------------------!

    ! set radius to either face or corners values
    SELECT CASE(GetType(Mesh))
    CASE(MIDPOINT)
      radius => RemapBounds(Mesh,Mesh%radius%faces)
    CASE(TRAPEZOIDAL)
      radius => RemapBounds(Mesh,Mesh%radius%corners)
    CASE DEFAULT
      CALL Error(Mesh,"InitData","Mesh not supported in this example.")
    END SELECT
    ! compute local sound speed using power law
    DO i=Mesh%IGMIN,Mesh%IGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
        ! speed of sound: bary center values
        bccsound(i,j) = CSISO*Mesh%radius%bcenter(i,j)**(0.5*C0)
        ! speed of sound: face or corner values (see above)
        fcsound(i,j,:) = CSISO*radius(i,j,:)**(0.5*C0)
        ! set y-velocity to local speed of sound for later use
        Timedisc%pvar(i,j,Physics%YVELOCITY) = bccsound(i,j)
      END DO
    END DO
    CALL SetSoundSpeeds(Physics,Mesh,bccsound)
    CALL SetSoundSpeeds(Physics,Mesh,fcsound)

    Timedisc%cvar(:,:,Physics%DENSITY) = Mesh%radius%bcenter**D0
    Timedisc%cvar(:,:,Physics%XMOMENTUM:Physics%YMOMENTUM) = 0.
    CALL Convert2Primitive(Physics,Mesh,Timedisc%cvar,Timedisc%pvar)

    Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = &
      Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) &
      + GetCentrifugalVelocity(Timedisc,Mesh,Physics,Fluxes,(/0.,0.,1./))

    ! boundary conditions
    DO dir=WEST,EAST
       SELECT CASE(GetType(Timedisc%Boundary(dir)))
       CASE(CUSTOM)
         Timedisc%boundary(dir)%cbtype(:,Physics%DENSITY) = CUSTOM_REFLECT
         Timedisc%boundary(dir)%cbtype(:,Physics%XVELOCITY) = CUSTOM_REFLNEG
         Timedisc%boundary(dir)%cbtype(:,Physics%YVELOCITY) = CUSTOM_KEPLER
       END SELECT
    END DO

    CALL Convert2Conservative(Physics,Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Info(Mesh," DATA-----> initial condition: " // "centrifugal balance")
  END SUBROUTINE InitData

END PROGRAM Init
