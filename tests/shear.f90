!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: shear.f90                                                         #
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
  INTEGER, PARAMETER :: RES_XY     = 128            ! resolution in x/y-direction
  INTEGER, PARAMETER :: RES_Z      = 1              ! resolution in z-direction
  REAL, PARAMETER    :: BOX_SIZE   = 320*GN*SIGMA0/(OMEGA*OMEGA)
  REAL, PARAMETER    :: BOX_HEIGHT = 0.0
  ! number of output time steps
  INTEGER, PARAMETER :: ONUM       = 10
  ! output directory and output name
  CHARACTER(LEN=256), PARAMETER :: ODIR   = "./"
  CHARACTER(LEN=256), PARAMETER :: OFNAME = "sbox"
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE :: Sim
  CLASS(marray_compound), POINTER :: pvar=>null(),pvar_init=>null()
  REAL               :: sigma
  LOGICAL            :: ok
  !--------------------------------------------------------------------------!
  TAP_PLAN(1)

  ! with south-north shear
  ALLOCATE(Sim)
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL Sim%Setup()
  CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc%pvar, Sim%Timedisc%cvar)
#ifndef PARALLEL
  ! store transposed initial data
  CALL Sim%Physics%new_statevector(pvar_init,PRIMITIVE)
  CALL RotateData(Sim%Mesh,Sim%Timedisc%pvar,pvar_init,"xy")
#endif
  CALL Sim%Run()
  ok = .NOT.Sim%aborted
#ifndef PARALLEL
  ! store transposed result of the first run
  CALL Sim%Physics%new_statevector(pvar,PRIMITIVE)
  CALL RotateData(Sim%Mesh,Sim%Timedisc%pvar,pvar,"xy")

  ! simulate west-east shear
  CALL Sim%InitFosite()
  CALL MakeConfig(Sim, Sim%config)
  CALL SetAttr(Sim%config, "/datafile/filename", (TRIM(ODIR) // TRIM(OFNAME) // "_rotate"))
  CALL SetAttr(Sim%config, "mesh/shearingbox", 2)
  CALL Sim%Setup()
  Sim%Timedisc%pvar%data1d(:) = pvar_init%data1d(:)
  SELECT TYPE(phys => Sim%Physics)
  CLASS IS(physics_eulerisotherm)
    CALL phys%Convert2Conservative(Sim%Timedisc%pvar,Sim%Timedisc%cvar)
  CLASS DEFAULT
    CALL phys%Error("shear","only (non-)isothermal HD supported")
  END SELECT
  CALL Sim%Run()

  ! compare results
  sigma = SQRT(SUM((Sim%Timedisc%pvar%data4d(:,:,:,:)-pvar%data4d(:,:,:,:))**2)/ &
                    SIZE(pvar%data4d(:,:,:,:)))
  DEALLOCATE(pvar,pvar_init)
#endif
  CALL Sim%Finalize()
  DEALLOCATE(Sim)

#ifndef PARALLEL
  ! x-y symmetry test does not work in parallel, because west-east shearing is not implemented
  TAP_CHECK_SMALL(sigma,1e-13,"Shear x-y symmetry test")
#else
  TAP_CHECK(ok,"stoptime reached")
#endif
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
    XMIN       = -0.5*BOX_SIZE
    XMAX       = +0.5*BOX_SIZE
    YMIN       = -0.5*BOX_SIZE
    YMAX       = +0.5*BOX_SIZE
    ZMIN       = -0.5*BOX_HEIGHT
    ZMAX       = +0.5*BOX_HEIGHT

    ! mesh settings
    mesh =>     Dict(&
                "meshtype"    / MIDPOINT, &
                "geometry"    / MGEO, &
                "shearingbox" / 1, &
                "decomposition" / (/1,-1,1/), & ! decompose along y-direction
                "inum"        / RES_XY, &
                "jnum"        / RES_XY, &
                "knum"        / RES_Z, &
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
    CLASS IS(statevector_eulerisotherm) ! isothermal HD
      WHERE ((ABS(Mesh%bcenter(:,:,:,1)).LT.0.15*BOX_SIZE.AND. &
              ABS(Mesh%bcenter(:,:,:,2)).LT.0.15*BOX_SIZE))
        p%density%data3d(:,:,:) = SIGMA0
      ELSEWHERE
        p%density%data3d(:,:,:) = SIGMA0*1e-2
      END WHERE
      SELECT CASE(Mesh%shear_dir)
      CASE(1)
        p%velocity%data2d(:,2) = 0.0
        p%velocity%data4d(:,:,:,1) = Mesh%Q*Mesh%bcenter(:,:,:,2)*Mesh%Omega
      CASE(2)
        p%velocity%data2d(:,1) = 0.0
        p%velocity%data4d(:,:,:,2) = -Mesh%Q*Mesh%bcenter(:,:,:,1)*Mesh%Omega
      CASE DEFAULT
        CALL Physics%Error("shear::InitData","shearing disabled or unsupported direction")
      END SELECT
    CLASS DEFAULT
      CALL Physics%Error("shear::InitData","only (non-)isothermal HD supported")
    END SELECT

    SELECT TYPE(p => pvar)
    CLASS IS(statevector_euler) ! non-isothermal HD
      p%pressure%data1d(:) = 1e-1
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
    SELECT CASE(dir)
    CASE("xy")
      SELECT TYPE(pin => pvar_in)
      CLASS IS(statevector_eulerisotherm) ! (non-)isothermal HD
        SELECT TYPE(pout => pvar_out)
        CLASS IS(statevector_eulerisotherm) ! (non-)isothermal HD
          ! rotate at middle of the field, because x' = -y in shearingsheet.
          FORALL(i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX,k=Mesh%KGMIN:Mesh%KGMAX)
            pout%density%data3d(i,j,k) = pin%density%data3d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k)
            pout%velocity%data4d(i,j,k,1) = pin%velocity%data4d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k,2)
            pout%velocity%data4d(i,j,k,2) = -pin%velocity%data4d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k,1)
          END FORALL
        END SELECT
      END SELECT
      SELECT TYPE(pin => pvar_in)
      CLASS IS(statevector_euler) ! non-isothermal HD
        SELECT TYPE(pout => pvar_out)
        CLASS IS(statevector_euler) ! non-isothermal HD
          ! rotate at middle of the field, because x' = -y in shearingsheet.
          FORALL(i=Mesh%IGMIN:Mesh%IGMAX,j=Mesh%JGMIN:Mesh%JGMAX,k=Mesh%KGMIN:Mesh%KGMAX)
            pout%pressure%data3d(i,j,k) = pin%pressure%data3d(Mesh%IGMAX-Mesh%IGMIN-j-2,i,k)
          END FORALL
        END SELECT
      END SELECT
    CASE DEFAULT
      CALL Mesh%Error("shear::RotateData","only 'xy'-plane is currently supported for rotating data")
    END SELECT
  END SUBROUTINE RotateData

END PROGRAM shear
