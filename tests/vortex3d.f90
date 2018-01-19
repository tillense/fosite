!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: vortex2d.f90                                                      #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> 2D isentropic vortex
!! \author Tobias Illenseer
!!
!! References:
!! [1] Yee, H. C. et al.: Low-dissipative high-order shock-capturing methods
!!     using characteristic-based filters, J. Comput. Phys. 150 (1999), 199-238
!!     DOI: 10.1006/jcph.1998.6177
!----------------------------------------------------------------------------!
PROGRAM vortex3d
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: TSIM    = 30.0     ! simulation stop time
  REAL, PARAMETER    :: GAMMA   = 1.4      ! ratio of specific heats
  REAL, PARAMETER    :: CSISO   = 1.0      ! if .ne. 0.0 -> isothermal simulation
  ! initial condition (dimensionless units)
  REAL, PARAMETER    :: RHOINF  = 1.       ! ambient density
  REAL, PARAMETER    :: PINF    = 1.       ! ambient pressure
  REAL, PARAMETER    :: VSTR    = 5.0      ! nondimensional vortex strength
  REAL, PARAMETER    :: UINF    = 0.0      ! cartesian components of constant
  REAL, PARAMETER    :: VINF    = 0.0      !   global velocity field
  REAL, PARAMETER    :: WINF    = 0.0
  REAL, PARAMETER    :: X0      = 0.0      ! vortex position (cart. coords.)
  REAL, PARAMETER    :: Y0      = 0.0
  REAL, PARAMETER    :: Z0      = 0.0
  REAL, PARAMETER    :: R0      = 1.0      ! size of vortex
  REAL, PARAMETER    :: OMEGA   = 0.0      ! angular speed of rotational frame
                                           ! around [X0,Y0]
  ! mesh settings
!!$  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry
!!$  INTEGER, PARAMETER :: MGEO = POLAR
  INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
!!$  INTEGER, PARAMETER :: MGEO = BIPOLAR
!!$  INTEGER, PARAMETER :: MGEO = ELLIPTIC
  INTEGER, PARAMETER :: XRES = 40         ! x-resolution
  INTEGER, PARAMETER :: YRES = 40         ! y-resolution
  INTEGER, PARAMETER :: ZRES = 40           ! y-resolution
  REAL, PARAMETER    :: RMIN = 1.0E-2      ! inner radius for polar grids
  REAL, PARAMETER    :: RMAX = 5.0         ! outer radius
  REAL, PARAMETER    :: GPAR = 1.0         ! geometry scaling parameter     !
  ! physics settings
!!$  LOGICAL, PARAMETER :: WITH_IAR = .TRUE.  ! use EULER2D_IAMROT
  LOGICAL, PARAMETER :: WITH_IAR = .FALSE.
  REAL, PARAMETER    :: PTB_AMP = 0.0E-02  ! amplitude of velocity perturbations
                                           !   set to zero to disable this
  ! output parameters
  INTEGER, PARAMETER :: ONUM = 10           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'vortex3d'
  !--------------------------------------------------------------------------!
  REAL               :: sigma
  CLASS(fosite),ALLOCATABLE :: Sim
  !--------------------------------------------------------------------------!

  TAP_PLAN(1)

  ALLOCATE(Sim)

  CALL Sim%InitFosite()

  CALL MakeConfig(Sim, Sim%config)

  CALL Sim%Setup()

  sigma = Run(Sim%Mesh, Sim%Physics, Sim%Timedisc)

  TAP_CHECK_SMALL(sigma,3.8E-3,"PP")

  DEALLOCATE(Sim)
  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config)
    USE functions, ONLY : ASINH
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(fosite)  :: Sim
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes, sources, rotframe
    REAL              :: x1,x2,y1,y2,z1,z2
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       x1 =-RMAX
       x2 = RMAX
       y1 =-RMAX
       y2 = RMAX
       z1 =-RMAX
       z2 = RMAX
    CASE(POLAR)
       x1 = RMIN
       x2 = RMAX
       y1 = 0.0
       y2 = 2.0*PI
    CASE(CYLINDRICAL)
       x1 = RMIN
       x2 = RMAX
       y1 =  0.0
       y2 =  2.0*PI
       z1 =  0.0
       z2 = 10.0
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
       x1 = ASINH(RMIN/GPAR)
       x2 = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = 2.0*PI
    CASE(BIPOLAR)
       x1 = 0.0
       x2 = 2.0*PI
       y1 = -ASINH(GPAR/RMIN)
       y2 = -y1
    CASE(ELLIPTIC)
       x1 = RMIN
       x2 = ASINH(RMAX/GPAR)
       y1 = 0.0
       y2 = 2.0*PI
    CASE DEFAULT
       CALL Sim%Error("InitProgram","mesh geometry not supported for 3D isentropic vortex")
    END SELECT

    !mesh settings
    mesh => Dict( &
              "meshtype" / MIDPOINT, &
              "geometry" / MGEO, &
              "omega"    / OMEGA, &
              "inum"     / XRES, &
              "jnum"     / YRES, &
              "knum"     / ZRES, &
              "xmin"     / x1, &
              "xmax"     / x2, &
              "ymin"     / y1, &
              "ymax"     / y2, &
              "zmin"     / z1, &
              "zmax"     / z2, &
              "gparam"   / GPAR)


    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       bc(WEST)   = NO_GRADIENTS
       bc(EAST)   = NO_GRADIENTS
       bc(SOUTH)  = NO_GRADIENTS
       bc(NORTH)  = NO_GRADIENTS
       bc(BOTTOM) = NO_GRADIENTS
       bc(TOP)    = NO_GRADIENTS
    CASE(CYLINDRICAL)
       bc(WEST)   = NO_GRADIENTS
       bc(EAST)   = NO_GRADIENTS
       bc(SOUTH)  = NO_GRADIENTS
       bc(NORTH)  = NO_GRADIENTS
       bc(BOTTOM) = NO_GRADIENTS
       bc(TOP)    = NO_GRADIENTS
    CASE DEFAULT
       CALL Sim%Error("InitProgram","mesh geometry not supported for 3D isentropic vortex")
    END SELECT

    ! boundary conditions
    boundary => Dict( &
                "western"   / bc(WEST), &
                "eastern"   / bc(EAST), &
                "southern"  / bc(SOUTH), &
                "northern"  / bc(NORTH), &
                "bottomer"  / bc(BOTTOM), &
                "topper"    / bc(TOP))

    ! physics settings
    IF (CSISO.GT.TINY(CSISO)) THEN
       physics => Dict("problem" / EULER3D_ISOTH, &
                 "cs"      / CSISO)                      ! isothermal sound speed  !
    ELSE
       IF (WITH_IAR) THEN
          ! REMARK: the optimal softening parameter depends on mesh geometry, limiter and
          ! possibly other settings; modify this starting with the default of 1.0, if the
          ! results show odd behaviour near the center of rotation; larger values increase
          ! softening; 0.5 give reasonable results for PP limiter on cartesian mesh
          physics => Dict("problem"   / EULER2D_IAMT, &
                     "centrot_x" / X0,"centrot_y" / Y0, & ! center of rotation      !
                     "softening" / 0.5, &                 ! softening parameter     !
                     "gamma"     / GAMMA)                 ! ratio of specific heats !
       ELSE
          physics => Dict("problem"   / EULER2D, &
                    "gamma"     / GAMMA)
       END IF
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict( &
              "fluxtype"  / KT, &
!              "order"     / CONSTANT, &
              "order"     / LINEAR, &
              "variables" / PRIMITIVE, &
              "limiter"   / VANLEER, &                      ! PP limiter gives better results than
              "theta"     / 1.0E-20, &                 ! should be < EPS for PP limiter
              "output/slopes" / 0)

    ! activate inertial forces due to rotating frame if OMEGA > 0
    NULLIFY(sources)
!    IF (OMEGA.GT.TINY(OMEGA)) THEN
!       rotframe => Dict("stype" / ROTATING_FRAME, &
!               "x"     / X0, &
!               "y"     / Y0)
!       sources => Dict("rotframe" / rotframe)
!    END IF

    ! time discretization settings
    timedisc => Dict(&
         "method"   / MODIFIED_EULER, &
         "order"    / 3, &
         "cfl"      / 0.4, &
         "stoptime" / TSIM, &
         "dtlimit"  / 1.0E-4, &
         "maxiter"  / 1000000, &
!          "output/fluxes" / 1, &
         "output/iangularmomentum" / 1, &
         "output/rothalpy" / 1)

    ! initialize data input/output
     datafile => Dict(&
          "fileformat" / VTK, &
!          "fileformat" / BINARY, &
          "filename"   / (TRIM(ODIR) // TRIM(OFNAME)), &
          "count"      / ONUM)
!    NULLIFY(datafile)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
             "datafile" / datafile)

    ! add sources terms
    IF (ASSOCIATED(sources)) &
        CALL SetAttr(config, "sources", sources)


  END SUBROUTINE MakeConfig


  FUNCTION Run(Mesh,Physics,Timedisc) RESULT(sigma)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(physics_base),INTENT(IN) :: Physics
    CLASS(mesh_base),INTENT(IN)    :: Mesh
    CLASS(timedisc_base),INTENT(INOUT):: Timedisc
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: i,j,k
    INTEGER           :: dir,ig
    INTEGER           :: n,clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: radius
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) :: dist_rot
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                      :: posvec,ephi,v0
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM) &
                      :: pvar0
    REAL              :: csinf,domega,sigma
    !------------------------------------------------------------------------!
    IF (ABS(X0).LE.TINY(X0).AND.ABS(Y0).LE.TINY(Y0).AND.ABS(Z0).LE.TINY(Z0)) THEN
       ! no shift of point mass set radius and posvec to Mesh defaults
       radius(:,:,:) = Mesh%radius%bcenter(:,:,:)
       posvec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:)
       dist_rot(:,:,:) = SQRT(posvec(:,:,:,1)**2 + posvec(:,:,:,2)**2)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       posvec(:,:,:,1) = X0
       posvec(:,:,:,2) = Y0
       posvec(:,:,:,3) = Z0
       CALL Mesh%geometry%Convert2Curvilinear(Mesh%bcenter,posvec,posvec)
       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing
       ! from the point mass to the bary center of any cell on the mesh
       posvec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - posvec(:,:,:,:)
       ! compute its absolute value
       radius(:,:,:) = SQRT(posvec(:,:,:,1)**2+posvec(:,:,:,2)**2+posvec(:,:,:,3)**2)
       dist_rot(:,:,:) = SQRT(posvec(:,:,:,1)**2+posvec(:,:,:,2)**2)
   END IF

    ! curvilinear components of azimuthal unit vector
    ! (maybe with respect to shifted origin)
    ! from ephi = ez x er = ez x posvec/radius = ez x (rxi*exi + reta*eeta)/r
    !             = rxi/r*(ez x exi) + reta/r*(ez x eeta) = rxi/r*eeta - reta/r*exi
    ! because (ez,exi,eeta) is right handed orthonormal set of basis vectors
    ephi(:,:,:,1) = -posvec(:,:,:,2)/radius(:,:,:)
    ephi(:,:,:,2) =  posvec(:,:,:,1)/radius(:,:,:)
    ephi(:,:,:,3) = 0.0

    csinf = SQRT(GAMMA*PINF/RHOINF) ! sound speed at infinity (isentropic vortex)
    ! initial condition depends on physics;
    ! could be either isothermal or isentropic
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
         DO i=Mesh%IGMIN,Mesh%IGMAX
            ! local angular velocity of the vortex
            domega = 0.5*VSTR/PI*EXP(0.5*(1.-(dist_rot(i,j,k)/R0)**2))
            SELECT CASE(Physics%GetType())
            CASE(EULER2D_ISOTHERM,EULER3D_ISOTH,EULER2D_ISOIAMT)
               ! density
               Timedisc%pvar(i,j,k,Physics%DENSITY) = RHOINF * EXP(-0.5*(R0*domega/CSISO)**2)
            CASE(EULER2D,EULER2D_IAMROT,EULER2D_IAMT)
               ! density
               ! ATTENTION: there's a factor of 1/PI missing in the density
               ! formula  eq. (3.3) in [1]
               Timedisc%pvar(i,j,k,Physics%DENSITY) = RHOINF * (1.0 - &
                    0.5*(GAMMA-1.0)*(R0*domega/csinf)**2  )**(1./(GAMMA-1.))
               ! pressure
               Timedisc%pvar(i,j,k,Physics%PRESSURE) = PINF &
                    * (Timedisc%pvar(i,j,k,Physics%DENSITY)/RHOINF)**GAMMA
            CASE DEFAULT
               CALL Physics%Error("InitData","Physics must be either EULER2D or EULER2D_ISOTHERM")
            END SELECT
            Timedisc%pvar(i,j,k,Physics%XVELOCITY:Physics%ZVELOCITY) = &
                    domega*dist_rot(i,j,k)*ephi(i,j,k,1:3)
         END DO
      END DO
    END DO

    ! compute curvilinear components of constant background velocity field
    ! and add to the vortex velocity field
    IF (ABS(UINF).GT.TINY(UINF).OR.ABS(VINF).GT.TINY(VINF).AND.ABS(WINF).LE.TINY(WINF)) THEN
       v0(:,:,:,1) = UINF
       v0(:,:,:,2) = VINF
       v0(:,:,:,3) = WINF
       CALL Mesh%geometry%Convert2Curvilinear(Mesh%bcenter,v0,v0)
       Timedisc%pvar(:,:,:,Physics%XVELOCITY:Physics%ZVELOCITY) = &
           Timedisc%pvar(:,:,:,Physics%XVELOCITY:Physics%ZVELOCITY) + v0(:,:,:,1:3)
    END IF

!    IF (ASSOCIATED(Sources)) &
!      CALL Convert2RotatingFrame_rotframe(&
!          GetSourcesPointer(Physics%sources, ROTATING_FRAME),&
!          Mesh,&
!          Physics,&
!          Timedisc%pvar)

    ! boundary conditions
!    DO dir=WEST,EAST
!       DO j=Mesh%JMIN,Mesh%JMAX
!          DO ig=1,Mesh%GNUM
!             SELECT CASE(dir)
!             CASE(WEST)
!                i = Mesh%IMIN-ig
!             CASE(EAST)
!                i = Mesh%IMAX+ig
!             END SELECT
!             SELECT CASE(Timedisc%Boundary(dir)%p%GetType())
!             CASE(NOSLIP)
!                Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) = &
!                    Timedisc%pvar(i,j,Physics%YVELOCITY)
!             CASE(CUSTOM)
!                Timedisc%boundary(dir)%cbtype(j,Physics%DENSITY) = CUSTOM_REFLECT
!                Timedisc%boundary(dir)%cbtype(j,Physics%XVELOCITY) = CUSTOM_REFLNEG
!                Timedisc%boundary(dir)%cbtype(j,Physics%YVELOCITY) = CUSTOM_FIXED
!                Timedisc%Boundary(dir)%data(ig,j,Physics%YVELOCITY) = &
!                    Timedisc%pvar(i,j,Physics%YVELOCITY)
!             END SELECT
!          END DO
!       END DO
!    END DO
!    DO dir=SOUTH,NORTH
!       DO i=Mesh%IMIN,Mesh%IMAX
!          DO ig=1,Mesh%GNUM
!             SELECT CASE(dir)
!             CASE(SOUTH)
!                j = Mesh%JMIN-ig
!             CASE(NORTH)
!                j = Mesh%JMAX+ig
!             END SELECT
!             SELECT CASE(GetType(Timedisc%Boundary(dir)))
!             CASE(NOSLIP)
!                Timedisc%Boundary(dir)%data(i,ig,Physics%XVELOCITY) = &
!                    Timedisc%pvar(i,j,Physics%XVELOCITY)
!             CASE(CUSTOM)
!                Timedisc%boundary(dir)%cbtype(i,Physics%DENSITY) = CUSTOM_REFLECT
!                Timedisc%boundary(dir)%cbtype(i,Physics%XVELOCITY) = CUSTOM_FIXED
!                Timedisc%boundary(dir)%cbtype(i,Physics%YVELOCITY) = CUSTOM_REFLNEG
!                Timedisc%Boundary(dir)%data(i,ig,Physics%XVELOCITY) = &
!                    Timedisc%pvar(i,j,Physics%XVELOCITY)
!             END SELECT
!          END DO
!       END DO
!    END DO

    ! add velocity perturbations if requested
    !IF (PTB_AMP.GT.TINY(PTB_AMP)) THEN
    !   ! Seed the random number generator with a mix from current time and mpi rank
    !   n = 4       ! has been choosen arbitrary
    !   CALL RANDOM_SEED(size=n)
    !   ALLOCATE(seed(n))
    !   CALL SYSTEM_CLOCK(COUNT=clock)
    !   seed = clock + (Timedisc%GetRank()+1) * (/(i-1, i=1,n)/)
    !   CALL RANDOM_SEED(PUT=seed)
    !   DEALLOCATE(seed)

    !   CALL RANDOM_NUMBER(v0)
    !   Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY) = (v0(:,:,1:2)-0.5)*PTB_AMP &
    !        + Timedisc%pvar(:,:,Physics%XVELOCITY:Physics%YVELOCITY)
    !END IF

    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar,Timedisc%cvar)
    CALL Mesh%Info(" DATA-----> initial condition: 3D vortex")

    pvar0 = Timedisc%pvar
    CALL Sim%Run()
    sigma = SQRT(SUM((Timedisc%pvar(:,:,:,:)-pvar0(:,:,:,:))**2)/SIZE(pvar0))
  END Function Run
END PROGRAM vortex3d
