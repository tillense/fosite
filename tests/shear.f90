!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: shear.f03                                                         #
!#                                                                           #
!# Copyright (C) 2008-2018                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Jannes Klee      <jklee@astrophysik.uni-kiel.de>                          #
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
!> \test Test shear flow
!! \author Jannes Klee
!----------------------------------------------------------------------------!
PROGRAM shear
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameters
  REAL, PARAMETER    :: GN         = 1.0            ! grav. constant [GEOM]  !
  INTEGER, PARAMETER :: UNITS      = GEOMETRICAL
  ! make a shearingsheet simulation and denote the direction of shearing
  ! (applies boundaries, ficitious forces and fargo automatically)
  INTEGER, PARAMETER :: SHEARSHEET_DIRECTION = 2
  ! simulation parameter
  REAL, PARAMETER    :: CSISO    = &
                                     0.0      ! non-isothermal simulation
!                                     1.0      ! isothermal simulation
                                               !   with CSISO as sound speed
  REAL, PARAMETER    :: OMEGA      = 1.0
  REAL, PARAMETER    :: SIGMA0     = 1.0
  REAL, PARAMETER    :: TSIM       = 10./OMEGA
  REAL, PARAMETER    :: GAMMA      = 2.0
  REAL, PARAMETER    :: Q          = 1.5            ! shearing parameter     !
  ! mesh settings
  INTEGER, PARAMETER :: MGEO       = CARTESIAN
  INTEGER, PARAMETER :: XRES       = 128
  INTEGER, PARAMETER :: YRES       = 128
  INTEGER, PARAMETER :: ZRES       = 1
  REAL               :: DOMAINX    = 320.0
  REAL               :: DOMAINY    = 320.0
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM       = 10
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sbox"
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  CLASS(marray_compound), POINTER :: pvar,pvar_init
  REAL               :: sigma
  !--------------------------------------------------------------------------!
  TAP_PLAN(1)

  ! with west-east shear
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
  ! store transposed initial data
  CALL Sim%Physics%new_statevector(pvar_init,PRIMITIVE)
  CALL RotateData(Sim%Mesh,Sim%Timedisc%pvar,pvar_init,"xy")
  CALL Sim%Run()
  ! store transposed result of the first run
  CALL Sim%Physics%new_statevector(pvar,PRIMITIVE)
  CALL RotateData(Sim%Mesh,Sim%Timedisc%pvar,pvar,"xy")
  ! finish the simulation
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

  ! simulate south-north shear
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_rotate"))
  CALL SetAttr(Sim%config, "mesh/shear_dir", 1)
  CALL Sim%Setup()
  Sim%Timedisc%pvar%data1d(:) = pvar_init%data1d(:)
  CALL Sim%Physics%Convert2Conservative(Sim%Timedisc%pvar,Sim%Timedisc%cvar)
  CALL Sim%Run()

  ! compare results
  sigma = SQRT(SUM((Sim%Timedisc%pvar%data4d(:,:,:,:)-pvar%data4d(:,:,:,:))**2)/ &
                    SIZE(pvar%data4d(:,:,:,:)))

  CALL pvar%Destroy()
  CALL Sim%Finalize()
  DEALLOCATE(pvar,Sim)

  TAP_CHECK_SMALL(sigma,1e-13,"Shear x-y symmetry test")
  TAP_DONE

  CONTAINS

  !> Set all configurations and safe it in a dictionary
  SUBROUTINE MakeConfig(Sim, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(fosite)            :: Sim
    TYPE(Dict_TYP), POINTER  :: config
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP), POINTER  :: mesh, physics, datafile, &
                                sources, timedisc, fluxes, shear
    !------------------------------------------------------------------------!
    INTENT(INOUT)            :: Sim
    REAL                     :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
    !------------------------------------------------------------------------!
    DOMAINX    = DOMAINX*GN*SIGMA0/(OMEGA*OMEGA)
    DOMAINY    = DOMAINY*GN*SIGMA0/(OMEGA*OMEGA)
    XMIN       = -0.5*DOMAINX
    XMAX       = +0.5*DOMAINX
    YMIN       = -0.5*DOMAINY
    YMAX       = +0.5*DOMAINY
    ZMIN       = 0.0
    ZMAX       = 0.0

    ! mesh settings
    mesh =>     Dict(&
                "meshtype"    / MIDPOINT, &
                "geometry"    / MGEO, &
                "shear_dir"   / SHEARSHEET_DIRECTION, &
                "inum"        / XRES, &
                "jnum"        / YRES, &
                "knum"        / ZRES, &
                "xmin"        / XMIN, &
                "xmax"        / XMAX, &
                "ymin"        / YMIN, &
                "ymax"        / YMAX, &
                "zmin"        / ZMIN, &
                "zmax"        / ZMAX, &
                "omega"       / OMEGA, &
                "output/rotation" / 0, &
                "output/volume"   / 0, &
                "output/bh"   / 0, &
                "output/dl"   / 0  &
                )

    ! physics settings
    IF (CSISO.GT.TINY(CSISO)) THEN
      physics => Dict("problem" / EULER_ISOTHERM, &
                        "units" / UNITS, &
                           "cs" / CSISO)             ! isothermal speed of sound
    ELSE
      physics => Dict("problem" / EULER, &
                        "units" / UNITS, &
                        "gamma" / GAMMA)             ! ratio of specific heats
    END IF

    ! flux calculation and reconstruction method
    fluxes => Dict( &
              "fluxtype"      / KT, &
              "order"         / LINEAR, &
              "variables"     / PRIMITIVE, &
              "limiter"       / VANLEER &
              )

    ! fluxes settings
    fluxes =>   Dict(&
                "order"       / LINEAR, &
                "fluxtype"    / KT, &
                "variables"   / PRIMITIVE, &
                "limiter"     / VANLEER &
                )

    ! shearing box fictious forces
    shear => Dict( &
                "stype"           / SHEARBOX, &
                "output/accel_x"  / 0, &
                "output/accel_y"  / 0 &
                )

    ! collect sources in dictionary
    sources => Dict("shear"   / shear)

    ! time discretization settings
    timedisc => Dict( &
              "method"        / MODIFIED_EULER, &
              "order"         / 3, &
              "cfl"           / 0.4, &
              "stoptime"      / TSIM, &
              "dtlimit"       / 1.0E-40, &
              "tol_rel"       / 0.1, &
              "output/external_sources" / 1, &
              "maxiter"       / 100000000)

    ! initialize data input/output
    datafile => Dict( &
              "fileformat"    / VTK, &
              "filename"      / (TRIM(ODIR) // TRIM(OFNAME)), &
              "count"         / ONUM)

    ! collect all above dicts in the configuration dict
    config => Dict( &
              "mesh"          / mesh, &
              "physics"       / physics, &
              "fluxes"        / fluxes, &
              "sources"       / sources, &
              "timedisc"      / timedisc, &
              "datafile"      / datafile)
  END SUBROUTINE MakeConfig


  !> Set initial conditions
  SUBROUTINE InitData(Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),    INTENT(IN) :: Mesh
    CLASS(physics_base), INTENT(IN) :: Physics
    CLASS(marray_compound), POINTER, INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    ! Test for shearing and boundary module
    ! initial condition
    SELECT TYPE(p => pvar)
    TYPE IS(statevector_eulerisotherm) ! isothermal HD
      WHERE ((ABS(Mesh%bcenter(:,:,:,1)).LT.0.15*DOMAINX.AND. &
              ABS(Mesh%bcenter(:,:,:,2)).LT.0.15*DOMAINY))
        p%density%data3d(:,:,:) = SIGMA0
      ELSEWHERE
        p%density%data3d(:,:,:) = SIGMA0*1e-2
      END WHERE
      p%velocity%data2d(:,1) = 0.0
      p%velocity%data4d(:,:,:,2) = -Mesh%Q*Mesh%bcenter(:,:,:,1)*Mesh%Omega
    TYPE IS(statevector_euler) ! non-isothermal HD
      WHERE ((ABS(Mesh%bcenter(:,:,:,1)).LT.0.15*DOMAINX.AND. &
              ABS(Mesh%bcenter(:,:,:,2)).LT.0.15*DOMAINY))
        p%density%data3d(:,:,:) = SIGMA0
      ELSEWHERE
        p%density%data3d(:,:,:) = SIGMA0*1e-2
      END WHERE
      p%velocity%data2d(:,1) = 0.0
      p%velocity%data4d(:,:,:,2) = -Mesh%Q*Mesh%bcenter(:,:,:,1)*Mesh%Omega
      p%pressure%data1d(:) = 1e-1
    CLASS DEFAULT
      CALL Physics%Error("shear::InitData","only non-isothermal HD supported")
    END SELECT

    CALL Physics%Convert2Conservative(pvar,cvar)
    CALL Mesh%Info(" DATA-----> initial condition: " // &
         "Shearing patch")
  END SUBROUTINE InitData

  ! only for symmetric matrix
  SUBROUTINE RotateData(Mesh,pvar_in,pvar_out,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base), INTENT(IN)            :: Mesh
    CLASS(marray_compound), INTENT(INOUT)   :: pvar_in,pvar_out
    CHARACTER(LEN=2), INTENT(IN)            :: dir
    !------------------------------------------------------------------------!
    INTEGER :: i,j,k
    !------------------------------------------------------------------------!
    SELECT TYPE(pin => pvar_in)
    TYPE IS(statevector_euler) ! non-isothermal HD
      SELECT TYPE(pout => pvar_out)
      TYPE IS(statevector_euler) ! non-isothermal HD
        SELECT CASE(dir)
        CASE("xy")
          ! rotate at middle of the field, because x' = -y in shearingsheet.
          DO k = Mesh%KGMIN,Mesh%KGMAX
            DO j = Mesh%JGMIN,Mesh%JGMAX
              DO i = Mesh%IGMIN,Mesh%IGMAX
                pout%density%data3d(i,j,k) = pin%density%data3d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k)
                pout%pressure%data3d(i,j,k) = pin%pressure%data3d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k)
                pout%velocity%data4d(i,j,k,1) = pin%velocity%data4d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k,2)
                pout%velocity%data4d(i,j,k,2) = pin%velocity%data4d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k,1)
              END DO
            END DO
          END DO
        CASE DEFAULT
          CALL Mesh%Error("shear::RotateData","directions must be one of 'xy','xz' or 'yz'")
        END SELECT
      END SELECT
    END SELECT
  END SUBROUTINE RotateData

END PROGRAM shear
