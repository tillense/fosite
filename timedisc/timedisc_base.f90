!r#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_base.f90                                                 #
!#                                                                           #
!# Copyright (C) 2007-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Manuel Jung      <mjung@astrophysik.uni-kiel.de>                          #
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
!> \addtogroup timedisc
!! \brief Base module for time discretization
!! - general parameters of timedisc group as key-values
!! \key{method,INTEGER,time integration method}
!! \key{rhstype,INTEGER,enable(1)/disabel(0) special right hand side treatment,0}
!! \key{stoptime,REAL,physical stop time of simulation}
!! \key{cfl,REAL,CFL number,0.4}
!! \key{dtlimit,REAL,time step minimum,EPSILON(dtlimit)*stoptime}
!! \key{dtmax,REAL,time step maximum in units of [CFL timestep]
!!   (Used in Dumka), 5}
!! \key{maxiter,INTEGER,maximum iterations,0}
!! \key{tol_rel,REAL,relative tolerance for adaptive step size control,0.01}
!! \key{tol_abs,REAL\, DIMENSION(Physics%VNUM), absolute tolerance for adaptive
!!   step size control, (/0.001\,0.001\,../)}
!! \key{beta,REAL,time step friction parameter for PI-Controller,0.08}
!! \key{always_update_bccsound,INTEGER,enable(1)/disable(0) update of bccsound at each
!!      substep,1}
!! \key{checkdata,INTEGER,bitmask to verify data validity/integrity,CHECK_ALL}
!! \key{pmin,REAL,pressure minimum to check if CHECK_PMIN has been enabled,TINY(REAL)}
!! \key{tmin,REAL,temperature minimum to check if CHECK_TMIN has been enabled,TINY(REAL)}
!! \key{rhomin,REAL,density minimum to check if CHECK_RHOMIN has been enabled,TINY(REAL)}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_base_mod
  USE logging_base_mod
  USE boundary_generic_mod
  USE mesh_base_mod
  USE marray_compound_mod
  USE marray_base_mod
  USE physics_base_mod
  USE physics_generic_mod
  USE sources_base_mod
  USE sources_generic_mod
  USE sources_gravity_mod
  USE fluxes_base_mod
  USE reconstruction_base_mod
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
!----------------------------------------------------------------------------!
PRIVATE
  TYPE, ABSTRACT, EXTENDS (logging_base) ::  timedisc_base
     !> \name Variables
     CLASS(boundary_generic), ALLOCATABLE :: Boundary  !< one for each boundary
     CLASS(marray_compound), POINTER  &
                      :: pvar => null(),cvar => null(), &     !< prim/cons state vectors
                         ptmp => null(),ctmp => null(), &
                         cold => null(), &
                         src  => null(),geo_src => null(), &  !< source terms
                         rhs => null(), &                     !< ODE right hand side
                         cerr => null(),cerr_max => null(), & !< error control & output
                         solution => null()                   !< analytical solution
     INTEGER          :: order                         !< time order
     REAL             :: cfl                           !< Courant number
     REAL             :: dt                            !< actual time step
     REAL             :: dtold                         !< last time step
     REAL             :: dtmin                         !< min dt of act. calc
     INTEGER          :: dtcause,dtmincause            !< cause of dt/dtmin
     REAL             :: dtmax                         !< max dt of act. calc
     REAL,POINTER     :: dtmean, dtstddev              !< mean and stddev of timestep
     INTEGER          :: dtaccept                      !< no of accepted ts
     REAL,POINTER     :: time                          !< actual time
     REAL             :: stoptime                      !< end of simulation
     REAL             :: dtlimit                       !< lower limit for dt
     REAL             :: tol_rel                       !< rel. error tolerance
     REAL             :: maxerrold                     !< old maxerr
     LOGICAL          :: break                         !< stop fosite run?
     LOGICAL          :: always_update_bccsound        !< controls update of bccsound
                      !! \details .TRUE. if bccsound should be updated at each substep
     INTEGER          :: maxiter                       !< maximal iterations
     INTEGER          :: n_adj                         !< num. of adjustments
     INTEGER          :: m                             !< cal steps in emb. RK
     INTEGER          :: degree                        !< polynom degree index in dumka
     INTEGER          :: rhstype                       !< rhs type ID
     INTEGER          :: checkdatabm                   !< Check Data Bitmask
     REAL             :: pmin,rhomin,tmin              !< minimum allowed pressure,
                                                       !< density,temperature
     TYPE(marray_base) :: dq, flux, dql                !< fargo fields
     !> old error and used in the dumka method
     REAL                              :: ERR_N, H_N
     REAL, DIMENSION(:), POINTER       :: tol_abs          !< abs. error tolerance
     REAL                              :: beta             !< time step friction

     !> multistep vars
     REAL, DIMENSION(:,:,:,:), POINTER :: phi,oldphi_s,&
                                          newphi_s
     REAL, DIMENSION(:), POINTER       :: gamma
     INTEGER                           :: pc               !< = 1 predictor-corrector
     !> numerical fluxes divided by dy or dx
     REAL, DIMENSION(:,:,:,:), POINTER :: xfluxdydz=>null(),yfluxdzdx=>null(),zfluxdxdy=>null()
     REAL, DIMENSION(:,:,:,:), POINTER :: amax=>null()    !< max. wave speeds
     REAL, DIMENSION(:,:), POINTER     :: bflux=>null()   !< boundary fluxes for output
     LOGICAL                           :: write_error     !< enable err writing
     INTEGER, DIMENSION(:,:), POINTER  :: shift=>null()   !< fargo annulus shift
     REAL, DIMENSION(:,:), POINTER     :: buf=>null()     !< fargo MPI buffer
     REAL, DIMENSION(:,:), POINTER     :: w=>null()       !< fargo background velocity
     REAL, DIMENSION(:,:), POINTER     :: wfactor=>null() !< fargo geometrical factor
     REAL, DIMENSION(:,:), POINTER     :: delxy =>null()  !< fargo residual shift

  CONTAINS

    PROCEDURE           :: InitTimedisc
    PROCEDURE           :: SetOutput
    PROCEDURE           :: IntegrationStep
    PROCEDURE           :: ComputeRHS
    PROCEDURE           :: AdjustTimestep
    PROCEDURE           :: CalcTimestep
    PROCEDURE           :: AcceptSolution
    PROCEDURE           :: RejectSolution
    PROCEDURE           :: CheckData
    PROCEDURE           :: ComputeError
    PROCEDURE           :: GetCentrifugalVelocity
    PROCEDURE           :: GetOrder
    PROCEDURE           :: GetCFL
    PROCEDURE           :: FargoAdvectionX
    PROCEDURE           :: FargoAdvectionY
    PROCEDURE           :: FargoAdvectionZ
    PROCEDURE           :: CalcBackgroundVelocity
    PROCEDURE (SolveODE), DEFERRED :: SolveODE
    PROCEDURE           :: Finalize_base
    PROCEDURE (Finalize), DEFERRED :: Finalize
  END TYPE timedisc_base
!----------------------------------------------------------------------------!
  ABSTRACT INTERFACE
    SUBROUTINE SolveODE(this,Mesh,Physics,Sources,Fluxes,time,dt,err)
      IMPORT timedisc_base, mesh_base, physics_base, fluxes_base, sources_base, sources_list
      IMPLICIT NONE
      CLASS(timedisc_base), INTENT(INOUT) :: this
      CLASS(mesh_base),     INTENT(IN)    :: Mesh
      CLASS(physics_base),  INTENT(INOUT) :: Physics
      CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
      CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
      REAL,                 INTENT(IN)    :: time
      REAL,                 INTENT(INOUT) :: dt, err
    END SUBROUTINE
  SUBROUTINE Finalize(this)
    IMPORT timedisc_base
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base) :: this
  END SUBROUTINE

  END INTERFACE

  INTEGER, PARAMETER :: MODIFIED_EULER   = 1
  INTEGER, PARAMETER :: RK_FEHLBERG      = 2
  INTEGER, PARAMETER :: CASH_KARP        = 3
  INTEGER, PARAMETER :: DORMAND_PRINCE   = 4
  INTEGER, PARAMETER :: SSPRK            = 5
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: DTCAUSE_CFL      =  0 ! smallest ts due to cfl cond. !
  INTEGER, PARAMETER :: DTCAUSE_ERRADJ   = -1 ! smallest ts due to err adj.  !
  INTEGER, PARAMETER :: DTCAUSE_SMALLERR = -2 ! smallest ts due to err
  INTEGER, PARAMETER :: DTCAUSE_FILEIO   = -4 !< smallest ts due to fileio
  ! CheckData_modeuler Bit Mask constants
  INTEGER, PARAMETER :: CHECK_NOTHING    = INT(b'0')
  INTEGER, PARAMETER :: CHECK_ALL        = NOT(CHECK_NOTHING)
  INTEGER, PARAMETER :: CHECK_CSOUND     = INT(b'1')
  INTEGER, PARAMETER :: CHECK_PMIN       = INT(b'10')
  INTEGER, PARAMETER :: CHECK_RHOMIN     = INT(b'100')
  INTEGER, PARAMETER :: CHECK_INVALID    = INT(b'1000')
  INTEGER, PARAMETER :: CHECK_TMIN       = INT(b'10000')
  !--------------------------------------------------------------------------!
  CHARACTER(LEN=40), PARAMETER, DIMENSION(3) :: FARGO_METHOD = (/ &
              "dynamic background velocity             ", &
              "user supplied fixed background velocity ", &
              "shearing box                            " /)
  PUBLIC :: &
       ! types
       timedisc_base, &
       ! constants
       MODIFIED_EULER, RK_FEHLBERG, CASH_KARP, DORMAND_PRINCE, &
       SSPRK, &
       DTCAUSE_CFL,DTCAUSE_ERRADJ,DTCAUSE_SMALLERR, DTCAUSE_FILEIO, &
       CHECK_ALL, CHECK_NOTHING, CHECK_CSOUND, CHECK_RHOMIN, CHECK_PMIN, &
       CHECK_INVALID, CHECK_TMIN
  !--------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE InitTimedisc(this,Mesh,Physics,config,IO,ttype,tname)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE physics_euler_mod, ONLY : physics_euler
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    USE geometry_logcylindrical_mod, ONLY: geometry_logcylindrical
    USE geometry_spherical_mod, ONLY: geometry_spherical
    USE geometry_cartesian_mod, ONLY: geometry_cartesian
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(INOUT) :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    TYPE(Dict_TYP),       POINTER       :: config,IO
    INTEGER,              INTENT(IN)    :: ttype
    CHARACTER(LEN=32),    INTENT(IN)    :: tname
    !------------------------------------------------------------------------!
    INTEGER              :: err, d, i,j,k
    CHARACTER(LEN=32)    :: order_str,cfl_str,stoptime_str,dtmax_str,beta_str
    CHARACTER(LEN=32)    :: info_str,shear_direction
    INTEGER              :: method
    !------------------------------------------------------------------------!
    CALL this%InitLogging(ttype,tname)

    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("InitTimedisc","physics and/or mesh module uninitialized")

    ! allocate memory for data structures needed in all timedisc modules
    ALLOCATE( &
      this%xfluxdydz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
      this%yfluxdzdx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
      this%zfluxdxdy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
      this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),                 &
      this%tol_abs(Physics%VNUM),        &
      this%dtmean, this%dtstddev,        &
      this%time,                         &
      STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("InitTimedisc", "Unable to allocate memory.")
    END IF

    ! initialize state vectors
    CALL Physics%new_statevector(this%pvar,PRIMITIVE)
    CALL Physics%new_statevector(this%ptmp,PRIMITIVE)
    CALL Physics%new_statevector(this%cvar,CONSERVATIVE)
    CALL Physics%new_statevector(this%ctmp,CONSERVATIVE)
    CALL Physics%new_statevector(this%cold,CONSERVATIVE)
    CALL Physics%new_statevector(this%geo_src,CONSERVATIVE)
    CALL Physics%new_statevector(this%src,CONSERVATIVE)
    CALL Physics%new_statevector(this%rhs,CONSERVATIVE)

    ! initialize all variables
    this%pvar%data1d(:)    = 0.
    this%ptmp%data1d(:)    = 0.
    this%ctmp%data1d(:)    = 0.
    this%cvar%data1d(:)    = 0.
    this%cold%data1d(:)    = 0.
    this%geo_src%data1d(:) = 0.
    this%src%data1d(:)     = 0.
    this%rhs%data1d(:)     = 0.
    this%xfluxdydz = 0.
    this%yfluxdzdx = 0.
    this%zfluxdxdy = 0.
    this%amax      = 0.
    this%tol_abs   = 0.
    this%dtmean    = 0.
    this%dtstddev  = 0.
    this%break     = .FALSE.
    this%n_adj     = 0
    this%maxerrold = 0.
    this%dtaccept  = 0

    CALL GetAttr(config, "starttime", this%time, 0.0)
    CALL GetAttr(config, "method", method)
    CALL GetAttr(config, "stoptime", this%stoptime)
    this%dt    = this%stoptime
    this%dtold = this%dt
    this%dtmin = this%dt

    ! set default values
    ! CFL number
    CALL GetAttr(config, "cfl", this%cfl, 0.4)

    ! time step minimum
    CALL GetAttr(config, "dtlimit", this%dtlimit, EPSILON(this%dtlimit)*this%stoptime)

    ! time step maximum in units of [CFL timestep] (Used in Dumka)
    CALL GetAttr(config, "dtmax", this%dtmax, 5.0)

    ! maximum iterations, maxiter <= 0 means no iteration limit
    CALL GetAttr(config, "maxiter", this%maxiter, 0)

    ! relative tolerance for adaptive step size control
    CALL GetAttr(config, "tol_rel", this%tol_rel, 0.01)

    ! absolute tolerance for adaptive step size control
    this%tol_abs(:) = 0.001
    CALL GetAttr(config, "tol_abs", this%tol_abs, this%tol_abs)

    ! rhs type. 0 = default, 1 = conserve angular momentum
    CALL GetAttr(config, "rhstype", this%rhstype, 0)
!    IF (this%rhstype.NE.0.AND.Physics%VDIM.NE.2) THEN
!      CALL this%Error("InitTimedisc", "Alternative rhstype only works in 2D.")
!    END IF

    ! time step friction parameter for PI-Controller
    CALL GetAttr(config, "beta", this%beta, 0.08)

    ! enable/disable update of bccsound at each substep
    ! Usually bccsound is only needed at the beginning of each integration step
    ! to determine the next time step. Hence it is not necessary to update
    ! bccsound at substeps when using higher order integration schemes. Since
    ! non-isothermal physics involve evaluation of the SQRT function to compute
    ! the speed of sound it might speed up the simulation a bit if one disables
    ! the update at substeps.
    ! Diabling the update is not possible if bccsound is needed somewhere else,
    ! e.g. some source terms depend on bccsound.
    CALL GetAttr(config, "always_update_bccsound", d, 1)
    IF (d.EQ.0) THEN
       this%always_update_bccsound = .FALSE.
    ELSE
       this%always_update_bccsound = .TRUE.
    END IF

    ! check data bit mask
    ! \todo{expected argument list...}
    CALL GetAttr(config, "checkdata", this%checkdatabm, CHECK_ALL)

    ! pressure minimum to check if CHECK_PMIN is active
    CALL GetAttr(config, "pmin", this%pmin, TINY(this%pmin))

    ! temperature minimum to check if CHECK_TMIN is active
    CALL GetAttr(config, "tmin", this%tmin, TINY(this%tmin))

    ! density minimum to check if CHECK_RHOMIN is active
    CALL GetAttr(config, "rhomin", this%rhomin, TINY(this%rhomin))

    CALL GetAttr(config, "output/"//TRIM(Physics%pvarname(1)), d, 1)

    CALL this%SetOutput(Mesh,Physics,config,IO)

    ! prepare data arrays for fargo transport
    ! check physics first; fargo is currently supported for physics_euler(isotherm) and derived types
    SELECT TYPE(phys => Physics)
    CLASS IS(physics_eulerisotherm)
      ! check for dimensionality of velocity vector
      IF (Mesh%fargo%GetType().GT.0) THEN
        IF (Physics%VDIM.LT.2) &
          CALL this%Error("timedisc_base::InitTimedisc","FARGO transport with 1D velocity vector not supported!")
      END IF
      ! check for fargo transport direction and allocate some arrays
      ! if direction is not in {1,2,3} skip this (fargo disabled)
      ! see mesh_base::InitMesh for general initialization of fargo transport
      SELECT CASE(Mesh%fargo%GetDirection())
      CASE(1) ! fargo transport along x-direction
        ALLOCATE( &
                  this%w(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                  this%delxy(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                  this%shift(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
#ifdef PARALLEL
                  this%buf(Physics%VNUM+Physics%PNUM,1:Mesh%MININUM), & !!! NOT CHANGED, because not clear were used !
#endif
                  STAT = err)
      CASE(2) ! fargo transport along y-direction
        ALLOCATE( &
                  this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                  this%delxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                  this%shift(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
#ifdef PARALLEL
                  this%buf(Physics%VNUM+Physics%PNUM,1:Mesh%MINJNUM), & !!! NOT CHANGED, because not clear were used !
#endif
                  STAT = err)
      CASE(3) ! fargo transport along z-direction
        ALLOCATE( &
                  this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                  this%delxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
                  this%shift(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
#ifdef PARALLEL
                  this%buf(Physics%VNUM+Physics%PNUM,1:Mesh%MINJNUM), & !!! NOT CHANGED, because not clear where used !
#endif
                  STAT = err)
      CASE DEFAULT
        err = 0
      END SELECT
      ! additional array used for computation of background velocity
      IF (err.EQ.0.AND.Mesh%fargo%GetType().EQ.1) ALLOCATE(this%wfactor,MOLD=this%w,STAT=err)
      IF (err.NE.0) &
          CALL this%Error("timedisc_base::InitTimedisc", "Unable to allocate memory for fargo advection.")

      ! initialize background velocity field w
      SELECT CASE(Mesh%fargo%GetType())
      CASE(0) ! fargo disabled -> do nothing
      CASE(1,2) ! fargo transport with fixed user supplied or dynamical velocity
        CALL this%Warning("timedisc_base::InitTimedisc","FARGO transport is experimental, use with care!")
                ! fargo advection type 1: w is computed in each time step (see CalcBackgroundVelocity below)
                ! fargo advection type 2: w is provided by the user, e.g. in InitData
        this%w(:,:) = 0.0
        this%shift(:,:) = 0.0
        this%dq = marray_base()
        this%dq%data1d(:) = 0.0
        this%dql = marray_base()
        this%dql%data1d(:) = 0.0
        this%flux = marray_base()
        this%flux%data1d(:) = 0.0
        IF (ASSOCIATED(this%wfactor)) THEN
          ! special geometrical factor for computation of dynamic background velocity;
          ! this factor is the ratio of the shortest width (in fargo transport direction)
          ! of a cell over all faces (perpendicular to the fargo transport direction) and
          ! the extend of the cell (in fargo transport direction) through the middle of the
          ! cell (see CalcBackgroundVelocity)
          SELECT CASE(Mesh%fargo%GetDirection())
          CASE(1) ! fargo transport along x-direction
            DO k=Mesh%KGMIN,Mesh%KGMAX
              DO j=Mesh%JGMIN,Mesh%JGMAX
                this%wfactor(j,k) = MIN(MINVAL(Mesh%hx%faces(Mesh%IMIN,j,k,SOUTH:NORTH)),&
                                        MINVAL(Mesh%hx%faces(Mesh%IMIN,j,k,BOTTOM:TOP))) &
                                     / Mesh%hx%bcenter(Mesh%IMIN,j,k)
              END DO
            END DO
          CASE(2) ! fargo transport along y-direction
            DO k=Mesh%KGMIN,Mesh%KGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%wfactor(i,k) = MIN(MINVAL(Mesh%hy%faces(i,Mesh%JMIN,k,WEST:EAST)),&
                                        MINVAL(Mesh%hy%faces(i,Mesh%JMIN,k,BOTTOM:TOP))) &
                                     / Mesh%hy%bcenter(i,Mesh%JMIN,k)
              END DO
            END DO
          CASE(3) ! fargo transport along z-direction
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%wfactor(i,j) = MIN(MINVAL(Mesh%hz%faces(i,j,Mesh%KMIN,WEST:EAST)),&
                                        MINVAL(Mesh%hz%faces(i,j,Mesh%KMIN,SOUTH:NORTH))) &
                                     / Mesh%hz%bcenter(i,j,Mesh%KMIN)
              END DO
            END DO
          END SELECT
        END IF
      CASE(3) ! fixed background velocity in shearing box
        IF(Mesh%shear_dir.EQ.2) THEN
          this%w(:,:) = -Mesh%Q*Mesh%omega*Mesh%bcenter(:,Mesh%JMIN,:,1) !-Q*Omega*x
          this%shift(:,:) = 0.0
          shear_direction = "west<->east"
        ELSE IF(Mesh%shear_dir.EQ.1) THEN
          this%w(:,:) = Mesh%Q*Mesh%omega*Mesh%bcenter(Mesh%IMIN,:,:,2) !Q*Omega*y
          this%shift(:,:) = 0.0
          shear_direction = "south<->north"
        END IF
        this%dq = marray_base()
        this%dq%data1d(:) = 0.0
        this%dql = marray_base()
        this%dql%data1d(:) = 0.0
        this%flux = marray_base()
        this%flux%data1d(:) = 0.0
      CASE DEFAULT
        CALL this%Error("timedisc_base::InitTimedisc","unknown fargo advection scheme")
      END SELECT
    CLASS DEFAULT
      CALL this%Error("timedisc_base::InitTimedisc","selected physics doesn't support fargo transport")
    END SELECT

    ! print some information
    WRITE (order_str, '(I0)') this%GetOrder()
    WRITE (cfl_str, '(F4.2)') this%GetCFL()
    WRITE (stoptime_str, '(ES12.4)') this%stoptime
    WRITE (dtmax_str, '(ES12.4)') this%dtmax
    WRITE (beta_str, '(ES12.4)') this%beta
    CALL this%Info(" TIMEDISC-> ODE solver:        " //TRIM(this%GetName()))
    CALL this%Info("            order:             " //TRIM(order_str))
    CALL this%Info("            CFL number:        " //TRIM(cfl_str))
    CALL this%Info("            dtmax:             " //TRIM(dtmax_str))
    CALL this%Info("            stoptime:          " //TRIM(stoptime_str))
    CALL this%Info("            beta:              " //TRIM(beta_str))
    ! adaptive step size control
    IF (this%tol_rel.LT.1.0.AND.this%order.GT.1) THEN
      ! create state vector to store the error
      CALL Physics%new_statevector(this%cerr,CONSERVATIVE)
      this%cerr%data1d(:)    = 0.
      ! create selection for the internal region
      WRITE (info_str,'(ES9.1)') this%tol_rel*100
      CALL this%Info("            step size control: enabled")
      CALL this%Info("            rel. precision:    "//TRIM(info_str)//" %")
    ELSE
       WRITE (info_str,'(A)') "disabled"
    END IF
    CALL this%Info("            fargo:             " //TRIM(Mesh%fargo%GetName()))
    IF(Mesh%fargo%GetType().GT.0) &
      CALL this%Info("            transport dir.:    " //TRIM(Mesh%fargo%GetDirectionName()))
    SELECT CASE(this%rhstype)
    CASE(0)
       ! special rhs disabled, print nothing
    CASE(1)
       IF (.NOT.ASSOCIATED(this%w)) THEN
          SELECT CASE(Mesh%Geometry%GetAzimuthIndex())
          CASE(2)
            ALLOCATE(this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX),STAT = err)
          CASE(3)
            ALLOCATE(this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),STAT = err)
          END SELECT
          IF (err.NE.0) THEN
             CALL this%Error("InitTimedisc", "Unable to allocate memory for special RHS time stepping.")
          END IF
          this%w(:,:) = 0.0
       END IF
       CALL this%Info("            special rhs:       enabled")
    CASE DEFAULT
       CALL this%Error("InitTimedisc","unknown rhstype")
    END SELECT
  END SUBROUTINE InitTimedisc


  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    TYPE(Dict_TYP),       POINTER       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER                             :: valwrite,i
    CHARACTER(LEN=128)                  :: key,pvar_key
    LOGICAL                             :: writeSolution
    !------------------------------------------------------------------------!
    valwrite = 0
    CALL GetAttr(config, "output/error", valwrite, 0)
    IF((valwrite.EQ.1).AND.this%tol_rel.LE.1.0) THEN
      this%write_error = .TRUE.
      CALL Physics%new_statevector(this%cerr_max,CONSERVATIVE)
      this%cerr_max%data1d(:) = 0.
    ELSE
      this%write_error = .FALSE.
    END IF

    valwrite = 0
    CALL GetAttr(config, "output/solution", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      writeSolution = .TRUE.
      CALL Physics%new_statevector(this%solution,PRIMITIVE)
    ELSE
      NULLIFY(this%solution)
      writeSolution = .FALSE.
    END IF

    CALL GetAttr(config, "output/time", valwrite, 1)
    IF(valwrite.EQ.1) &
      CALL SetAttr(IO, "time", Ref(this%time))

    DO i=1, Physics%VNUM
      !prim
      key = TRIM(Physics%pvarname(i))
      pvar_key = key
      valwrite = 0
      CALL GetAttr(config, "output/" // TRIM(key), valwrite, 1)
      ! second argument is important if pvarname is used twice
      IF (valwrite.EQ.1) THEN
        CALL SetAttr(IO, TRIM(key), &
                     this%pvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      IF(writeSolution) THEN
        CALL SetAttr(IO, TRIM(key)//"_ref", &
          this%solution%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      !cons
      key = TRIM(Physics%cvarname(i))
      IF(key.EQ.pvar_key) key = TRIM(key) // "_cvar"
      valwrite = 0
      CALL GetAttr(config, "output/" // TRIM(key), valwrite, 0)
      ! second argument is important if pvarname is used twice
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key), &
                      this%cvar%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      key = TRIM(Physics%cvarname(i))
      IF(this%write_error) THEN
        CALL SetAttr(IO, TRIM(key) // "_error", &
                     this%cerr_max%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write geometrical sources
      CALL GetAttr(config, "output/" // "geometrical_sources", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_geo_src", &
                        this%geo_src%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write external sources
      CALL GetAttr(config, "output/" // "external_sources", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_src", &
                        this%src%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write right hand side
      CALL GetAttr(config, "output/" // "rhs", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_rhs", &
                        this%rhs%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write fluxes
      ! ATTENTION: this are the numerical fluxes devided by dy or dx respectively
      CALL GetAttr(config, "output/" // "fluxes", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_xfluxdydz", &
                        this%xfluxdydz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
           CALL SetAttr(IO, TRIM(key)//"_yfluxdzdx", &
                        this%yfluxdzdx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
           CALL SetAttr(IO, TRIM(key)//"_zfluxdxdy", &
                        this%zfluxdxdy(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

    END DO
    ! write bfluxes
    CALL GetAttr(config, "output/" // "bflux", valwrite, 0)
    IF(valwrite.EQ.1) THEN
        ALLOCATE(this%bflux(Physics%VNUM,6))
        CALL SetAttr(IO, "bflux", this%bflux)
    ELSE
        NULLIFY(this%bflux)
    END IF
  END SUBROUTINE SetOutput


  SUBROUTINE IntegrationStep(this,Mesh,Physics,Sources,Fluxes,iter,config,IO)
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(INOUT) :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    INTEGER                             :: iter
    TYPE(Dict_TYP),       POINTER       :: config,IO
    !------------------------------------------------------------------------!
    REAL                                :: err,dtold,dt,time
    !------------------------------------------------------------------------!
    INTENT(INOUT)                       :: iter
    !------------------------------------------------------------------------!
    ! determine the background velocity if fargo advection type 1 is enabled
    IF (Mesh%fargo%GetType().EQ.1) &
       CALL this%CalcBackgroundVelocity(Mesh,Physics,this%pvar,this%cvar)
    ! transform to comoving frame if fargo is enabled
    SELECT CASE(Mesh%fargo%GetDirection())
    CASE(1)
      CALL Physics%SubtractBackgroundVelocityX(Mesh,this%w,this%pvar,this%cvar)
    CASE(2)
      CALL Physics%SubtractBackgroundVelocityY(Mesh,this%w,this%pvar,this%cvar)
    CASE(3)
      CALL Physics%SubtractBackgroundVelocityZ(Mesh,this%w,this%pvar,this%cvar)
    END SELECT

    time = this%time
    dt   = this%dt
    ! save dtmin if it is not determined by the fileio, i.e. time
    ! between two successive write operations
    IF (dt.LT.this%dtmin .AND. this%dtcause .NE. DTCAUSE_FILEIO) THEN
      this%dtmin = dt
      this%dtmincause = this%dtcause
    END IF

    ! loop until new solution for time+dt is accepted
    timestep: DO WHILE (.NOT.this%break.AND.time+dt.LE.this%time+this%dt)
      dtold = dt

      CALL this%SolveODE(Mesh,Physics,Sources,Fluxes,time,dt,err)

      ! check truncation error and restart if necessary,
      ! but first check data integrity of err
      IF (err.NE.err) THEN
        ! err = NaN
        CALL this%Warning("timedisc_base::IntegrationStep","max error returned by SolveODE is NaN")
        this%break = .TRUE.
      ELSE IF (err.GT.HUGE(err)) THEN
        ! err = Inf
        CALL this%Warning("timedisc_base::IntegrationStep","max error returned by SolveODE is Inf")
        this%break = .TRUE.
      ELSE IF (err.GE.1.0) THEN
        ! err >= 1.0
        CALL this%RejectSolution(Mesh,Physics,Sources,Fluxes,time,dt)
      ELSE
        ! err < 1.0
        CALL this%AcceptSolution(Mesh,Physics,Sources,Fluxes,time,dtold,iter)      
      END IF

      ! abort integration if lower time step limit reached
      IF (dt.LT.this%dtlimit) this%break = .TRUE.

    END DO timestep

    IF (.NOT.this%break) THEN
      ! Save true advanced time (for fargo linear advection)
      this%dt = time - this%time

      this%time  = time
      this%dtold = dt

      ! perform the fargo advection step if enabled
      IF (Mesh%fargo%GetType().GT.0) THEN
        ! ATTENTION: Boundary conditions are applied to data with subtracted background
        !            velocity here which may yield errornous results for the real boundaries
        !            but not for the periodic boundaries which is essential for the
        !            fargo advection step.
        !            Since we call ComputeRHS after the fargo step, boundary conditions
        !            are applied again and this time with added background velocity,
        !            at least if fargo method is set to 1 or 2 (see ComputeRHS).
        CALL this%boundary%CenterBoundary(Mesh,Physics,this%time,this%pvar,this%cvar)
        SELECT CASE(Mesh%fargo%GetDirection())
        CASE(1)
          CALL this%FargoAdvectionX(Fluxes,Mesh,Physics,Sources)
        CASE(2)
          CALL this%FargoAdvectionY(Fluxes,Mesh,Physics,Sources)
        CASE(3)
          CALL this%FargoAdvectionZ(Fluxes,Mesh,Physics,Sources)
        END SELECT
        ! convert conservative to primitive variables
        CALL Physics%Convert2Primitive(this%cvar,this%pvar)
        ! Calculate RHS after the Advection Step
        CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,this%time,0.,this%pvar,this%cvar,&
                        this%checkdatabm,this%rhs)
      ELSE
        CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,this%time,dtold,this%pvar,this%cvar,&
                        this%checkdatabm,this%rhs)
      END IF
      this%cold = this%cvar
    END IF
  END SUBROUTINE IntegrationStep


  !> \public adjust the time step
  !!
  !! This function implements adaptive step size control based on an
  !! error estimate.
  !!
  !! see E. Hairer, Solving Ordinary Differential Equ. II, 2ed, Springer (2.43c)
  FUNCTION AdjustTimestep(this,maxerr,dtold) RESULT(dtnew)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base) :: this
    REAL                 :: maxerr,dtold,dtnew
    !------------------------------------------------------------------------!
    REAL                 :: dttmp
    !------------------------------------------------------------------------!
    INTENT(IN)           :: maxerr,dtold
    INTENT(INOUT)        :: this
    !------------------------------------------------------------------------!
    ! adaptive step size control disabled
    IF (this%tol_rel.GE.1.0) THEN
      dtnew = HUGE(dtnew)
      RETURN
    END IF

    ! time step estimate with maxerr determined from the comparison of numerical
    ! results of two explicit schemes with different order where GetOrder(this)
    ! returns the order of the higher order scheme (P-Controller)

    IF(maxerr.GT.0.) THEN
      dttmp = 0.9*dtold*EXP(-LOG(maxerr)/GetOrder(this))
    ELSE
      ! ths limit for maxerr->0 is dttmp -> oo. But this cut off to
      ! 4*dtold later in this function.
      dttmp = 4.*dtold
    END IF

    ! this controls the increase/decrease of time steps based on the
    ! old value of maxerr from the previous time step (PI-Controller)
    IF (this%maxerrold.GT.0.) THEN
       dtnew = dttmp * EXP(-this%beta*LOG(this%maxerrold/maxerr))
    ELSE
       dtnew = dttmp
    END IF

    ! adjust the time step based on the error estimate
    IF (maxerr.LT.1.0) THEN
       ! increase time step
       dtnew = MIN(dtnew,4.*dtold)   ! enlarge at most by a factor of 4
       ! it is possible: maxerr < 1 =>  0.9*dt < dtnew => maybe dtnew < dt
       IF (dtnew .LT. dtold) this%dtcause = DTCAUSE_SMALLERR
    ELSE
       ! If maxerrold >> maxerr dtnew can become larger than dtold even
       ! if the timestep is rejected. If this is the case fall back to
       ! the simple P-Controller.
       IF (dtnew.GT.dtold) dtnew = dttmp
       ! decrease time step
       dtnew = MAX(dtnew,0.25*dtold)    ! reduce at most by a factor of 1/4
    END IF

    ! store maxerr for next time step
    this%maxerrold = maxerr

  END FUNCTION AdjustTimestep


  !> \public Determines the CFL time step and time step limits due to source terms
  REAL FUNCTION CalcTimestep(this,Mesh,Physics,Sources,Fluxes,time,dtcause) RESULT(dt)
    USE physics_euler_mod, ONLY : physics_euler, statevector_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL,                 INTENT(IN)    :: time
    INTEGER,              INTENT(INOUT) :: dtcause
    !------------------------------------------------------------------------!
    REAL                                :: invdt
    REAL                                :: dt_cfl, dt_src
    REAL                                :: invdt_x,invdt_y,invdt_z
    !------------------------------------------------------------------------!
    ! CFL condition:
    ! maximal wave speeds in each direction
    IF (.NOT.this%always_update_bccsound) THEN
      SELECT TYPE(phys => Physics)
      CLASS IS(physics_euler)
        SELECT TYPE(pvar => this%pvar)
        CLASS IS(statevector_euler)
          CALL phys%UpdateSoundSpeed(pvar)
        END SELECT
      END SELECT
    END IF

    SELECT TYPE(phys => Physics)
      TYPE IS(physics_eulerisotherm)
        SELECT TYPE(pvar => this%pvar)
        CLASS IS(statevector_eulerisotherm)
          CALL phys%CalculateWaveSpeeds(Mesh,pvar,Fluxes%minwav,Fluxes%maxwav)
        END SELECT
      TYPE IS(physics_euler)
        SELECT TYPE(pvar => this%pvar)
        CLASS IS(statevector_euler)
          CALL phys%CalculateWaveSpeeds(Mesh,pvar,Fluxes%minwav,Fluxes%maxwav)
        END SELECT
    END SELECT

    ! compute maximum of inverse time for CFL condition
    IF ((Mesh%JNUM.EQ.1).AND.(Mesh%KNUM.EQ.1)) THEN
       ! 1D, only x-direction
       invdt = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dlx%data1d(:))
    ELSE IF ((Mesh%INUM.EQ.1).AND.(Mesh%KNUM.EQ.1)) THEN
       ! 1D, only y-direction
       invdt = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dly%data1d(:))
    ELSE IF ((Mesh%INUM.EQ.1).AND.(Mesh%JNUM.EQ.1)) THEN
       ! 1D, only z-direction
       invdt = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dlz%data1d(:))
    ELSE IF ((Mesh%INUM.GT.1).AND.(Mesh%JNUM.GT.1).AND.(Mesh%KNUM.EQ.1)) THEN
       ! 2D, x-y-plane
       invdt = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dlx%data1d(:) &
                    + MAX(Fluxes%maxwav%data2d(:,2),-Fluxes%minwav%data2d(:,2)) / Mesh%dly%data1d(:))
    ELSE IF ((Mesh%INUM.GT.1).AND.(Mesh%KNUM.GT.1).AND.(Mesh%JNUM.EQ.1)) THEN
       ! 2D, x-z-plane
       invdt = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dlx%data1d(:) &
                    + MAX(Fluxes%maxwav%data2d(:,2),-Fluxes%minwav%data2d(:,2)) / Mesh%dlz%data1d(:))
    ELSE IF ((Mesh%JNUM.GT.1).AND.(Mesh%KNUM.GT.1).AND.(Mesh%INUM.EQ.1)) THEN
       ! 2D, y-z-plane
       invdt = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dly%data1d(:) &
                    + MAX(Fluxes%maxwav%data2d(:,2),-Fluxes%minwav%data2d(:,2)) / Mesh%dlz%data1d(:))
    ELSE
       ! full 3D
       !TODO: Achtung: Hier wurde fuer eine bessere Symmetrie fuer jede Richtung ein eigenes invdt
       ! berechnet. Dies koennte jedoch einen Verlust an Stabilitaet bewirken. Hier muesste mal eine
       ! Stabilitaetsanalyse gemacht werden
       invdt_x = MAXVAL(MAX(Fluxes%maxwav%data2d(:,1),-Fluxes%minwav%data2d(:,1)) / Mesh%dlx%data1d(:))
       invdt_y = MAXVAL(MAX(Fluxes%maxwav%data2d(:,2),-Fluxes%minwav%data2d(:,2)) / Mesh%dly%data1d(:))
       invdt_z = MAXVAL(MAX(Fluxes%maxwav%data2d(:,3),-Fluxes%minwav%data2d(:,3)) / Mesh%dlz%data1d(:))
       invdt = MAX(invdt_y,invdt_z,invdt_x)
     !  invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,3),-Fluxes%minwav(:,:,:,3)) / Mesh%dlz(:,:,:) &
     !               + MAX(Fluxes%maxwav(:,:,:,2),-Fluxes%minwav(:,:,:,2)) / Mesh%dly(:,:,:) &
     !               + MAX(Fluxes%maxwav(:,:,:,1),-Fluxes%minwav(:,:,:,1)) / Mesh%dlx(:,:,:))
    END IF

    ! largest time step due to CFL condition
    dt_cfl = this%cfl / invdt
    ! time step limitation due to cfl -> dtcause = 0
    dtcause = 0

    ! check for sources
    IF (ALLOCATED(Sources)) THEN
      ! initialize this to be sure dt_src > 0
      dt_src = dt_cfl
      CALL Sources%CalcTimestep(Mesh,Physics,Fluxes,this%pvar,this%cvar,time,dt_src,dtcause)
      dt = MIN(dt_cfl,dt_src)
    ELSE
      dt = dt_cfl
    END IF
  END FUNCTION CalcTimestep

  SUBROUTINE AcceptSolution(this,Mesh,Physics,Sources,Fluxes,time,dt,iter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL                                :: time,dt
    INTEGER                             :: iter
    !------------------------------------------------------------------------!
    REAL                                :: dtmeanold
    REAL, DIMENSION(:), POINTER         :: bflux,bfold
    !------------------------------------------------------------------------!
    INTENT(INOUT)                       :: time,dt
    !------------------------------------------------------------------------!
    time=time+dt
    this%dtaccept = this%dtaccept + 1
    dtmeanold = this%dtmean
    this%dtmean = this%dtmean + (dt - this%dtmean)/this%dtaccept
    this%dtstddev = this%dtstddev + (dt - dtmeanold)*(dt-this%dtmean)
    IF (time.LT.this%time+this%dt) THEN
       CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,time,dt,this%pvar,this%cvar,&
                            this%checkdatabm,this%rhs)
       this%cold = this%cvar
    END IF
    IF (Mesh%INUM.GT.1) THEN
      ! collapse the 4D arrays to improve vector performance;
      ! maybe we can switch to the old style assignments if
      ! automatic collapsing of NEC nfort compiler works
      bflux(1:SIZE(Fluxes%bxflux)) => Fluxes%bxflux
      bfold(1:SIZE(Fluxes%bxfold)) => Fluxes%bxfold
      bfold(:) = bflux(:)
!       Fluxes%bxfold(:,:,:,:) = Fluxes%bxflux(:,:,:,:)
    END IF
    IF (Mesh%JNUM.GT.1) THEN
      bflux(1:SIZE(Fluxes%byflux)) => Fluxes%byflux
      bfold(1:SIZE(Fluxes%byfold)) => Fluxes%byfold
      bfold(:) = bflux(:)
!       Fluxes%byfold(:,:,:,:) = Fluxes%byflux(:,:,:,:)
    END IF
    IF (Mesh%KNUM.GT.1) THEN
      bflux(1:SIZE(Fluxes%bzflux)) => Fluxes%bzflux
      bfold(1:SIZE(Fluxes%bzfold)) => Fluxes%bzfold
      bfold(:) = bflux(:)
!       Fluxes%bzfold(:,:,:,:) = Fluxes%bzflux(:,:,:,:)
    END IF
    iter = iter + 1
  END SUBROUTINE AcceptSolution

  SUBROUTINE RejectSolution(this,Mesh,Physics,Sources,Fluxes,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL                                :: time,dt
    !------------------------------------------------------------------------!
    INTEGER                             :: i,j,k,l,m
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: time
    INTENT(INOUT)                       :: dt
    !------------------------------------------------------------------------!
    this%cvar = this%cold
    ! This data has already been checked before in AcceptSolution
    CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,time,dt,this%pvar, &
      this%cvar,CHECK_NOTHING,this%rhs)
!NEC$ SHORTLOOP
    DO m=1,Physics%VNUM
!NEC$ UNROLL(2)
       DO l=1,2
          DO k=Mesh%KGMIN,Mesh%KGMAX
             DO j=Mesh%JGMIN,Mesh%JGMAX
                Fluxes%bxflux(j,k,l,m) = Fluxes%bxfold(j,k,l,m)
             END DO
          END DO
          DO i=Mesh%IGMIN,Mesh%IGMAX
             DO k=Mesh%KGMIN,Mesh%KGMAX
                Fluxes%byflux(k,i,l,m) = Fluxes%byfold(k,i,l,m)
             END DO
          END DO
          DO j=Mesh%JGMIN,Mesh%JGMAX
             DO i=Mesh%IGMIN,Mesh%IGMAX
                Fluxes%bzflux(i,j,l,m) = Fluxes%bzfold(i,j,l,m)
             END DO
          END DO
       END DO
    END DO
    ! count adjustments for information
    this%n_adj = this%n_adj + 1
    this%dtcause = DTCAUSE_ERRADJ
    ! only save dtmin if the reason is not the fileio (automatically satisfied)
    IF (dt.LT.this%dtmin) THEN
       this%dtmin = dt
       this%dtmincause = this%dtcause
    END IF
  END SUBROUTINE RejectSolution


  FUNCTION ComputeError(this,Mesh,Physics,cvar_high,cvar_low) RESULT(maxerr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    CLASS(marray_compound),INTENT(INOUT):: cvar_high,cvar_low
    REAL               :: maxerr
    !------------------------------------------------------------------------!
    INTEGER            :: i,l
#ifdef PARALLEL
    INTEGER            :: ierror
#endif
    REAL               :: rel_err(Physics%VNUM)
    !------------------------------------------------------------------------!
!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM
      DO CONCURRENT (i=1:SIZE(this%cerr%data2d,DIM=1))
        ! compute the local error (including ghost zones)
        this%cerr%data2d(i,l) = ABS(cvar_high%data2d(i,l)-cvar_low%data2d(i,l)) &
                         / (this%tol_rel*ABS(cvar_high%data2d(i,l)) + this%tol_abs(l))
      END DO
      ! determine the global maximum on the whole grid except for ghost zones
      rel_err(l) = MAXVAL(this%cerr%data2d(:,l),MASK=Mesh%without_ghost_zones%mask1d(:))
      ! store the maximum between two output time steps
      IF (this%write_error) THEN
        DO CONCURRENT (i=1:SIZE(this%cerr_max%data2d,DIM=1))
          this%cerr_max%data2d(i,l) = MAX(this%cerr_max%data2d(i,l),this%cerr%data2d(i,l))
        END DO
      END IF
    END DO

    ! compute the maximum of all variables
    maxerr = MAXVAL(rel_err(:))

#ifdef PARALLEL
    ! find the global maximum of all processes
    CALL MPI_Allreduce(MPI_IN_PLACE,maxerr,1,DEFAULT_MPI_REAL,MPI_MAX,&
                       Mesh%comm_cart,ierror)
#endif

  END FUNCTION ComputeError



  !> \public compute the RHS of the spatially discretized PDE
  !!
  !! This subroutine updates all data on the right hand side (RHS) of the
  !! ordinary differential equation obtained via spatial discretization
  !! of the system of partial differential equations (PDE). It expects that
  !! the conservative variables (cvar) contain valid data inside the
  !! computational domain. On the basis of this data the following
  !! update steps are performed in the given order:
  !!
  !! 1. update of ghost cell data (this implies an update of the
  !!    primitive variables (pvar) on the whole grid including ghost cells)
  !! 2. update of intercell numerical fluxes (this implies an update of
  !!    all face states and additional data necessary to compute the fluxes)
  !! 3. update geometrical source terms
  !! 4. update external source terms (this implies an update of all
  !!    auxiliary data arrays used for sources terms)
  SUBROUTINE ComputeRHS(this,Mesh,Physics,Sources,Fluxes,time,dt,pvar,cvar,checkdatabm,rhs)
    USE physics_euler_mod, ONLY : physics_euler,statevector_euler
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE geometry_spherical_mod, ONLY: geometry_spherical
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL,                 INTENT(IN)    :: time, dt
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar,rhs
    INTEGER,              INTENT(IN)    :: checkdatabm
    !------------------------------------------------------------------------!
    INTEGER                             :: i,j,k,l
    LOGICAL                             :: have_potential
    REAL                                :: t
    REAL                                :: phi,wp
    REAL, DIMENSION(:,:,:,:), POINTER   :: pot
    CLASS(sources_base),      POINTER   :: sp => NULL()
    CLASS(sources_gravity),   POINTER   :: grav => NULL()
    !------------------------------------------------------------------------!
    t = time
    SELECT CASE(Mesh%fargo%GetType())
    CASE(1,2)
        ! transform to real velocity, i.e. v_residual + w_background,
        ! before setting the boundary conditions
        SELECT CASE(Mesh%fargo%GetDirection())
        CASE(1)
          CALL Physics%AddBackgroundVelocityX(Mesh,this%w,pvar,cvar)
        CASE(2)
          CALL Physics%AddBackgroundVelocityY(Mesh,this%w,pvar,cvar)
        CASE(3)
          CALL Physics%AddBackgroundVelocityZ(Mesh,this%w,pvar,cvar)
        END SELECT
    CASE(3)
        ! boundary conditions are set for residual velocity in shearing box simulations
        ! usually it's not necessary to subtract the background velocity here, but
        ! in case someone adds it before, we subtract it here; the  subroutine checks,
        ! if the velocities have already been transformed
        SELECT CASE(Mesh%fargo%GetDirection())
        CASE(1)
          CALL Physics%SubtractBackgroundVelocityX(Mesh,this%w,pvar,cvar)
        CASE(2)
          CALL Physics%SubtractBackgroundVelocityY(Mesh,this%w,pvar,cvar)
        CASE(3)
          CALL Physics%SubtractBackgroundVelocityZ(Mesh,this%w,pvar,cvar)
        END SELECT
        ! ATTENTION: the time must be the initial time of the whole time step
        !            not the time of a substep if fargo type 3 (shearing box)
        !            is enabled
        t = this%time
    CASE DEFAULT
        ! fargo disabled (do nothing)
    END SELECT

    ! set boundary values and convert conservative to primitive variables
    ! ATTENTION: if fargo type 3 is enabled, background velocity is subtracted
    !            afterwards, otherwise not
    CALL this%boundary%CenterBoundary(Mesh,Physics,t,pvar,cvar)

    IF(IAND(checkdatabm,CHECK_TMIN).NE.CHECK_NOTHING.AND.&
      this%tmin.GT.1.E-10) THEN
      SELECT TYPE(p => pvar)
      CLASS IS(statevector_euler)
        SELECT TYPE(c => cvar)
        CLASS IS(statevector_euler)
          ! If temperature is below TMIN limit pressure to density*RG/MU*TMIN.
          p%pressure%data1d(:) = MAX(p%pressure%data1d(:), &
            p%density%data1d(:)*Physics%Constants%RG/Physics%MU*this%TMIN)
          CALL Physics%Convert2Conservative(p,c)
        END SELECT
      END SELECT
    END IF

    IF(IAND(checkdatabm,CHECK_PMIN).NE.CHECK_NOTHING.AND.&
       Physics%PRESSURE.GT.0) THEN
      ! Check if the pressure is below pmin. If it is, increase the pressure
      ! to reach pmin
      SELECT TYPE(p => pvar)
      CLASS IS(statevector_euler)
        SELECT TYPE(c => cvar)
        CLASS IS(statevector_euler)
          ! If temperature is below PMIN limit pressure to PMIN
          p%pressure%data1d(:) &
            = MAX(p%pressure%data1d(:),this%pmin)
          CALL Physics%Convert2Conservative(p,c)
        END SELECT
      END SELECT
    END IF

    ! update the speed of sound (non-isotherml physics only)
    IF (this%always_update_bccsound) THEN
      SELECT TYPE(phys => Physics)
      CLASS IS(physics_euler)
        SELECT TYPE(pvar => this%pvar)
        CLASS IS(statevector_euler)
          CALL phys%UpdateSoundSpeed(pvar)
        END SELECT
      END SELECT
    END IF

    ! check for illegal data
    IF(checkdatabm.NE.CHECK_NOTHING) &
      CALL this%CheckData(Mesh,Physics,Fluxes,pvar,cvar,checkdatabm)

    ! get geometrical sources
    CALL Physics%GeometricalSources(Mesh,pvar,cvar,this%geo_src)

    ! get source terms due to external forces if present
    IF (ALLOCATED(Sources)) &
       CALL Sources%ExternalSources(Mesh,Physics,Fluxes,Sources, &
            t,dt,pvar,cvar,this%src)

    ! if fargo advection is enabled additional source terms occur;
    ! furthermore computation of numerical fluxes should always be
    ! carried out with residual velocity
    SELECT CASE(Mesh%fargo%GetType())
    CASE(1,2)
        ! add fargo source terms to geometrical source terms and
        ! subtract background velocity
        ! zero geometry sources for cartesian simulations first
        IF (Mesh%geometry%GetType().EQ.CARTESIAN) this%geo_src%data1d(:) = 0.0
        SELECT CASE(Mesh%fargo%GetDirection())
        CASE(1)
          CALL Physics%AddFargoSourcesX(Mesh,this%w,pvar,cvar,this%geo_src)
          CALL Physics%SubtractBackgroundVelocityX(Mesh,this%w,pvar,cvar)
        CASE(2)
          CALL Physics%AddFargoSourcesY(Mesh,this%w,pvar,cvar,this%geo_src)
          CALL Physics%SubtractBackgroundVelocityY(Mesh,this%w,pvar,cvar)
        CASE(3)
          CALL Physics%AddFargoSourcesZ(Mesh,this%w,pvar,cvar,this%geo_src)
          CALL Physics%SubtractBackgroundVelocityZ(Mesh,this%w,pvar,cvar)
        END SELECT
    CASE(3)
        ! background velocity field has already been subtracted (do nothing);
        ! fargo specific source terms are handled in the shearing box source
        ! term module
    CASE DEFAULT
        ! fargo disabled (do nothing)
    END SELECT

    ! get the numerical fluxes
    CALL Fluxes%CalculateFluxes(Mesh,Physics,pvar,cvar,this%xfluxdydz, &
                                this%yfluxdzdx,this%zfluxdxdy)

    ! compute the right hand side for boundary flux computation;
    ! this is probably wrong for the special rhs (rhstype=1, see above)
!NEC$ SHORTLOOP
    rhsbdflux: DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ IVDEP
      DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO j=Mesh%JMIN,Mesh%JMAX
          ! western and eastern boundary fluxes
          rhs%data4d(Mesh%IMIN-Mesh%Ip1,j,k,l) = Mesh%dy*Mesh%dz * this%xfluxdydz(Mesh%IMIN-Mesh%Ip1,j,k,l)
          rhs%data4d(Mesh%IMAX+Mesh%Ip1,j,k,l) = -Mesh%dy*Mesh%dz * this%xfluxdydz(Mesh%IMAX,j,k,l)
        END DO
      END DO
!NEC$ IVDEP
      DO k=Mesh%KMIN,Mesh%KMAX
!NEC$ IVDEP
        DO i=Mesh%IMIN,Mesh%IMAX
          ! southern and northern boundary fluxes
          rhs%data4d(i,Mesh%JMIN-Mesh%Jp1,k,l) = Mesh%dz*Mesh%dx * this%yfluxdzdx(i,Mesh%JMIN-Mesh%Jp1,k,l)
          rhs%data4d(i,Mesh%JMAX+Mesh%Jp1,k,l) = -Mesh%dz*Mesh%dx * this%yfluxdzdx(i,Mesh%JMAX,k,l)
        END DO
      END DO
!NEC$ IVDEP
      DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
        DO i=Mesh%IMIN,Mesh%IMAX
          ! bottomer and upper boundary fluxes
          rhs%data4d(i,j,Mesh%KMIN-Mesh%Kp1,l) = Mesh%dx*Mesh%dy * this%zfluxdxdy(i,j,Mesh%KMIN-Mesh%Kp1,l)
          rhs%data4d(i,j,Mesh%KMAX+Mesh%Kp1,l) = -Mesh%dx*Mesh%dy * this%zfluxdxdy(i,j,Mesh%KMAX,l)
        END DO
      END DO
    END DO rhsbdflux

    SELECT CASE(this%rhstype)
    CASE(0)
!NEC$ SHORTLOOP
      DO l=1,Physics%VNUM+Physics%PNUM
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i=Mesh%IMIN,Mesh%IMAX
              ! update right hand side of ODE
              rhs%data4d(i,j,k,l) = &
                    Mesh%dydzdV%data3d(i,j,k)*(this%xfluxdydz(i,j,k,l) - this%xfluxdydz(i-Mesh%Ip1,j,k,l)) &
                  + Mesh%dzdxdV%data3d(i,j,k)*(this%yfluxdzdx(i,j,k,l) - this%yfluxdzdx(i,j-Mesh%Jp1,k,l)) &
                  + Mesh%dxdydV%data3d(i,j,k)*(this%zfluxdxdy(i,j,k,l) - this%zfluxdxdy(i,j,k-Mesh%Kp1,l)) &
                  - this%geo_src%data4d(i,j,k,l) - this%src%data4d(i,j,k,l)
            END DO
          END DO
        END DO
      END DO

    CASE(1)
      !> \todo Very hacky implementation of pluto style angular momentum
      !! conservation. A better implementation is needed, but we would need to
      !! restucture the timedisc module
      IF (ALLOCATED(Sources)) &
        sp => Sources%GetSourcesPointer(GRAVITY)

      NULLIFY(pot)
      have_potential = .FALSE.
      IF(ASSOCIATED(grav)) THEN
        IF(.NOT.grav%addtoenergy) THEN
          pot => Mesh%RemapBounds(grav%pot%data4d(:,:,:,:))
          have_potential=.TRUE.
        END IF
      END IF

      phi = 0.

    SELECT CASE(Mesh%geometry%GetAzimuthIndex())
    CASE(2)  ! Cylindrical
      DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX
        DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX
            wp = 0.5*(this%w(i,k)+this%w(i+1,k)) + Mesh%hy%faces(i+1,j,k,1)*Mesh%OMEGA
            IF(Physics%PRESSURE.GT.0) THEN
              this%xfluxdydz(i,j,k,Physics%ENERGY) &
                = this%xfluxdydz(i,j,k,Physics%ENERGY) &
                  + wp * (0.5 * wp * this%xfluxdydz(i,j,k,Physics%DENSITY) &
                  + this%xfluxdydz(i,j,k,Physics%YMOMENTUM))
              IF (have_potential) THEN
                this%xfluxdydz(i,j,k,Physics%ENERGY) = this%xfluxdydz(i,j,k,Physics%ENERGY) &
                  + pot(i,j,k,2)*this%xfluxdydz(i,j,k,Physics%DENSITY)
              END IF
            END IF

            this%xfluxdydz(i,j,k,Physics%YMOMENTUM) &
              = (this%xfluxdydz(i,j,k,Physics%YMOMENTUM) &
                + wp * this%xfluxdydz(i,j,k,Physics%DENSITY)) * Mesh%hy%faces(i+Mesh%IP1,j,k,1)

            wp = this%w(i,k) + Mesh%hy%bcenter(i,j,k)*Mesh%OMEGA

            IF(Physics%PRESSURE.GT.0) THEN
              this%yfluxdzdx(i,j,k,Physics%ENERGY) &
                = this%yfluxdzdx(i,j,k,Physics%ENERGY) &
                 + wp * ( 0.5 * wp * this%yfluxdzdx(i,j,k,Physics%DENSITY) &
                 + this%yfluxdzdx(i,j,k,Physics%YMOMENTUM))
              IF (have_potential) THEN
                this%yfluxdzdx(i,j,k,Physics%ENERGY) = this%yfluxdzdx(i,j,k,Physics%ENERGY) &
                  + pot(i,j,k,3)*this%yfluxdzdx(i,j,k,Physics%DENSITY)
              END IF
            END IF

            this%yfluxdzdx(i,j,k,Physics%YMOMENTUM) &
              = this%yfluxdzdx(i,j,k,Physics%YMOMENTUM) &
                + wp * this%yfluxdzdx(i,j,k,Physics%DENSITY)

            wp = this%w(i,k) + Mesh%hy%bcenter(i,j,k)*Mesh%OMEGA

            IF(Physics%PRESSURE.GT.0) THEN
              this%zfluxdxdy(i,j,k,Physics%ENERGY) &
                = this%zfluxdxdy(i,j,k,Physics%ENERGY) &
                 + wp * ( 0.5 * wp * this%zfluxdxdy(i,j,k,Physics%DENSITY) &
                 + this%zfluxdxdy(i,j,k,Physics%YMOMENTUM))
              IF (have_potential) THEN
                this%zfluxdxdy(i,j,k,Physics%ENERGY) = this%zfluxdxdy(i,j,k,Physics%ENERGY) &
                  + pot(i,j,k,4)*this%zfluxdxdy(i,j,k,Physics%DENSITY)
              END IF
            END IF

            this%zfluxdxdy(i,j,k,Physics%YMOMENTUM) &
              = (this%zfluxdxdy(i,j,k,Physics%YMOMENTUM) &
                + wp * this%zfluxdxdy(i,j,k,Physics%DENSITY)) * Mesh%hy%faces(i,j,k+Mesh%KP1,1)

          END DO
        END DO
      END DO

      ! set isothermal sound speeds
      SELECT TYPE (phys => Physics)
      CLASS IS (physics_eulerisotherm)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i=Mesh%IMIN,Mesh%IMAX
              wp = pvar%data4d(i,j,k,Physics%YVELOCITY) + this%w(i,k) + Mesh%hy%bcenter(i,j,k)*Mesh%OMEGA

              this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                 = pvar%data4d(i,j,k,Physics%DENSITY) * wp * (wp*Mesh%cyxy%bcenter(i,j,k) - &
                                                            pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxyx%bcenter(i,j,k))

              !cyxy and cyzy are eleminated due to the modified divergence operator
              this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                = cvar%data4d(i,j,k,Physics%XMOMENTUM) * ( pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxyx%center(i,j,k))

              IF(Physics%ZMOMENTUM.GT.0) THEN
                this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) = this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                  + cvar%data4d(i,j,k,Physics%ZVELOCITY) * (pvar%data4d(i,j,k,Physics%ZVELOCITY)*Mesh%czxz%bcenter(i,j,k) - &
                                                            pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxzx%bcenter(i,j,k))
                this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) = this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                  + cvar%data4d(i,j,k,Physics%ZMOMENTUM) * ( pvar%data4d(i,j,k,Physics%ZVELOCITY) * Mesh%czyz%center(i,j,k))
                this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) = cvar%data4d(i,j,k,Physics%XMOMENTUM) &
                  * ( pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxzx%center(i,j,k) &
                    - pvar%data4d(i,j,k,Physics%ZVELOCITY) *  Mesh%czxz%center(i,j,k) ) &
                  + pvar%data4d(i,j,k,Physics%DENSITY)*wp* ( wp*Mesh%cyzy%center(i,j,k) &
                  - pvar%data4d(i,j,k,Physics%ZVELOCITY) * Mesh%czyz%center(i,j,k) )
              END IF

              IF(Physics%PRESSURE.GT.0) THEN
                this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                   +pvar%data4d(i,j,k,Physics%PRESSURE) &
                      *( Mesh%cyxy%center(i,j,k) + Mesh%czxz%center(i,j,k) )
                this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%PRESSURE) &
                      *( Mesh%cxyx%center(i,j,k) + Mesh%czyz%center(i,j,k) )
              IF(Physics%ZMOMENTUM.GT.0) &
                this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%PRESSURE) &
                      *( Mesh%cxzx%center(i,j,k) + Mesh%cyzy%center(i,j,k) )

                this%geo_src%data4d(i,j,k,Physics%ENERGY) = 0.
              ELSE

                this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                 = this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                   +pvar%data4d(i,j,k,Physics%DENSITY)*phys%bccsound%data3d(i,j,k)**2 &
                      * ( Mesh%cyxy%center(i,j,k) + Mesh%czxz%center(i,j,k) )
                this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%DENSITY)*phys%bccsound%data3d(i,j,k)**2 &
                      *( Mesh%cxyx%center(i,j,k) + Mesh%czyz%center(i,j,k) )
              IF(Physics%ZMOMENTUM.GT.0) &
                this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%DENSITY)*phys%bccsound%data3d(i,j,k)**2 &
                      *( Mesh%cxzx%center(i,j,k) + Mesh%cyzy%center(i,j,k) )

              END IF
            END DO
          END DO
        END DO
      CLASS DEFAULT
        CALL this%Error("timedisc_base::ComputeRHS","physics not supported with rhstype=1")
      END SELECT


!NEC$ NOVECTOR
      DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ OUTERLOOP_UNROLL(8)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i=Mesh%IMIN,Mesh%IMAX
              ! update right hand side of ODE
              IF(l.EQ.Physics%YMOMENTUM) THEN
                rhs%data4d(i,j,k,l) = Mesh%dydzdV%data3d(i,j,k)*(this%xfluxdydz(i,j,k,l) - this%xfluxdydz(i-Mesh%IP1,j,k,l)) &
                                          / Mesh%hy%center(i,j,k) &
                                    + Mesh%dzdxdV%data3d(i,j,k)*(this%yfluxdzdx(i,j,k,l) - this%yfluxdzdx(i,j-Mesh%JP1,k,l)) &
                                    + Mesh%dxdydV%data3d(i,j,k)*(this%zfluxdxdy(i,j,k,l) - this%zfluxdxdy(i,j,k-Mesh%KP1,l)) &
                                          / Mesh%hy%center(i,j,k)
              ELSE
                rhs%data4d(i,j,k,l) = Mesh%dydzdV%data3d(i,j,k)*(this%xfluxdydz(i,j,k,l) - this%xfluxdydz(i-Mesh%IP1,j,k,l)) &
                                    + Mesh%dzdxdV%data3d(i,j,k)*(this%yfluxdzdx(i,j,k,l) - this%yfluxdzdx(i,j-Mesh%JP1,k,l)) &
                                    + Mesh%dxdydV%data3d(i,j,k)*(this%zfluxdxdy(i,j,k,l) - this%zfluxdxdy(i,j,k-Mesh%KP1,l))
              END IF
            END DO
          END DO
        END DO
      END DO

      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN,Mesh%IMAX
            wp = this%w(i,k) + Mesh%radius%bcenter(i,j,k) * Mesh%OMEGA
            rhs%data4d(i,j,k,Physics%YMOMENTUM) = rhs%data4d(i,j,k,Physics%YMOMENTUM) &
              - wp * rhs%data4d(i,j,k,Physics%DENSITY)
            IF(Physics%PRESSURE.GT.0) THEN
              rhs%data4d(i,j,k,Physics%ENERGY) = rhs%data4d(i,j,k,Physics%ENERGY) &
                - wp*( rhs%data4d(i,j,k,Physics%YMOMENTUM) &
                       + 0.5 * wp* rhs%data4d(i,j,k,Physics%DENSITY))
              IF (have_potential) THEN
                rhs%data4d(i,j,k,Physics%ENERGY) = &
                rhs%data4d(i,j,k,Physics%ENERGY) - pot(i,j,k,1) * rhs%data4d(i,j,k,Physics%DENSITY)
              END IF
            END IF
          END DO
        END DO
      END DO

    CASE(3)   !spherical, tancylindrical
      DO k=Mesh%KMIN-Mesh%KP1,Mesh%KMAX
        DO j=Mesh%JMIN-Mesh%JP1,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN-Mesh%IP1,Mesh%IMAX
            wp = 0.5*(this%w(i,j)+this%w(i+1,j)) + Mesh%hz%faces(i+1,j,k,1)*Mesh%OMEGA
            IF(Physics%PRESSURE.GT.0) THEN
              this%xfluxdydz(i,j,k,Physics%ENERGY) &
                = this%xfluxdydz(i,j,k,Physics%ENERGY) &
                  + wp * (0.5 * wp * this%xfluxdydz(i,j,k,Physics%DENSITY) &
                  + this%xfluxdydz(i,j,k,Physics%ZMOMENTUM))
              IF (have_potential) THEN
                this%xfluxdydz(i,j,k,Physics%ENERGY) = this%xfluxdydz(i,j,k,Physics%ENERGY) &
                  + pot(i,j,k,2)*this%xfluxdydz(i,j,k,Physics%DENSITY)
              END IF
            END IF

            this%xfluxdydz(i,j,k,Physics%ZMOMENTUM) &
              = (this%xfluxdydz(i,j,k,Physics%ZMOMENTUM) &
                + wp * this%xfluxdydz(i,j,k,Physics%DENSITY)) * Mesh%hz%faces(i+1,j,k,1)

            wp = this%w(i,j) + Mesh%hz%bcenter(i,j,k)*Mesh%OMEGA

            IF(Physics%PRESSURE.GT.0) THEN
              this%yfluxdzdx(i,j,k,Physics%ENERGY) &
                = this%yfluxdzdx(i,j,k,Physics%ENERGY) &
                 + wp * ( 0.5 * wp * this%yfluxdzdx(i,j,k,Physics%DENSITY) &
                 + this%yfluxdzdx(i,j,k,Physics%ZMOMENTUM))
              IF (have_potential) THEN
                this%yfluxdzdx(i,j,k,Physics%ENERGY) = this%yfluxdzdx(i,j,k,Physics%ENERGY) &
                  + pot(i,j,k,3)*this%yfluxdzdx(i,j,k,Physics%DENSITY)
              END IF
            END IF

            this%yfluxdzdx(i,j,k,Physics%ZMOMENTUM) &
              = (this%yfluxdzdx(i,j,k,Physics%ZMOMENTUM) &
                + wp * this%yfluxdzdx(i,j,k,Physics%DENSITY)) * Mesh%hz%faces(i,j+1,k,1)

            wp = this%w(i,j) + Mesh%hz%bcenter(i,j,k)*Mesh%OMEGA

            IF(Physics%PRESSURE.GT.0) THEN
              this%zfluxdxdy(i,j,k,Physics%ENERGY) &
                = this%zfluxdxdy(i,j,k,Physics%ENERGY) &
                 + wp * ( 0.5 * wp * this%zfluxdxdy(i,j,k,Physics%DENSITY) &
                 + this%zfluxdxdy(i,j,k,Physics%ZMOMENTUM))
              IF (have_potential) THEN
                this%zfluxdxdy(i,j,k,Physics%ENERGY) = this%zfluxdxdy(i,j,k,Physics%ENERGY) &
                  + pot(i,j,k,4)*this%zfluxdxdy(i,j,k,Physics%DENSITY)
              END IF
            END IF

            this%zfluxdxdy(i,j,k,Physics%ZMOMENTUM) &
              = this%zfluxdxdy(i,j,k,Physics%ZMOMENTUM) &
                + wp * this%zfluxdxdy(i,j,k,Physics%DENSITY)

          END DO
        END DO
      END DO

      ! set isothermal sound speeds
      SELECT TYPE (phys => Physics)
      CLASS IS (physics_eulerisotherm)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i=Mesh%IMIN,Mesh%IMAX
              wp = pvar%data4d(i,j,k,Physics%ZVELOCITY) + this%w(i,j) + Mesh%hz%bcenter(i,j,k)*Mesh%OMEGA

              this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                 = cvar%data4d(i,j,k,Physics%YMOMENTUM) * (pvar%data4d(i,j,k,Physics%YVELOCITY)*Mesh%cyxy%bcenter(i,j,k) - &
                                                            pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxyx%bcenter(i,j,k))  &
                  + pvar%data4d(i,j,k,Physics%DENSITY)*wp * (wp*Mesh%czxz%bcenter(i,j,k) - &
                                                            pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxzx%bcenter(i,j,k))

              this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                = cvar%data4d(i,j,k,Physics%XMOMENTUM) &
                  * ( pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxyx%center(i,j,k) &
                    - pvar%data4d(i,j,k,Physics%YVELOCITY) * Mesh%cyxy%center(i,j,k) ) &
                + pvar%data4d(i,j,k,Physics%DENSITY) * wp &
                  * ( wp * Mesh%czyz%center(i,j,k) - pvar%data4d(i,j,k,Physics%YVELOCITY) * Mesh%cyzy%center(i,j,k) )

              ! czxz and czyz are eleminated due to the modified divergence operator
              IF(Physics%ZMOMENTUM.GT.0) &
              this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                = cvar%data4d(i,j,k,Physics%XMOMENTUM) * ( pvar%data4d(i,j,k,Physics%XVELOCITY) * Mesh%cxzx%center(i,j,k) ) &
                + cvar%data4d(i,j,k,Physics%YMOMENTUM) * ( pvar%data4d(i,j,k,Physics%YVELOCITY) * Mesh%cyzy%center(i,j,k) )

              IF(Physics%PRESSURE.GT.0) THEN
                this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                   +pvar%data4d(i,j,k,Physics%PRESSURE) &
                      *( Mesh%cyxy%center(i,j,k) + Mesh%czxz%center(i,j,k) )
                this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%PRESSURE) &
                      *( Mesh%cxyx%center(i,j,k) + Mesh%czyz%center(i,j,k) )
              IF(Physics%ZMOMENTUM.GT.0) &
                this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%PRESSURE) &
                      *( Mesh%cxzx%center(i,j,k) + Mesh%cyzy%center(i,j,k) )

                this%geo_src%data4d(i,j,k,Physics%ENERGY) = 0.
              ELSE

                this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                 = this%geo_src%data4d(i,j,k,Physics%XMOMENTUM) &
                   +pvar%data4d(i,j,k,Physics%DENSITY)*phys%bccsound%data3d(i,j,k)**2 &
                      * ( Mesh%cyxy%center(i,j,k) + Mesh%czxz%center(i,j,k) )
                this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%YMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%DENSITY)*phys%bccsound%data3d(i,j,k)**2 &
                      *( Mesh%cxyx%center(i,j,k) + Mesh%czyz%center(i,j,k) )
              IF(Physics%ZMOMENTUM.GT.0) &
                this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                  = this%geo_src%data4d(i,j,k,Physics%ZMOMENTUM) &
                    + pvar%data4d(i,j,k,Physics%DENSITY)*phys%bccsound%data3d(i,j,k)**2 &
                      *( Mesh%cxzx%center(i,j,k) + Mesh%cyzy%center(i,j,k) )

              END IF
            END DO
          END DO
        END DO

      CLASS DEFAULT
        ! abort
        CALL this%Error("timedisc_base::ComputeRHS","physics not supported with rhstype=1")
      END SELECT


!NEC$ NOVECTOR
      DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ OUTERLOOP_UNROLL(8)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i=Mesh%IMIN,Mesh%IMAX
              ! update right hand side of ODE
              IF(l.EQ.Physics%ZMOMENTUM) THEN
                rhs%data4d(i,j,k,l) = Mesh%dydzdV%data3d(i,j,k)*(this%xfluxdydz(i,j,k,l) - this%xfluxdydz(i-Mesh%IP1,j,k,l)) &
                                        / Mesh%hz%center(i,j,k) &
                                    + Mesh%dzdxdV%data3d(i,j,k)*(this%yfluxdzdx(i,j,k,l) - this%yfluxdzdx(i,j-Mesh%JP1,k,l)) &
                                        / Mesh%hz%center(i,j,k) &
                                    + Mesh%dxdydV%data3d(i,j,k)*(this%zfluxdxdy(i,j,k,l) - this%zfluxdxdy(i,j,k-Mesh%KP1,l))
              ELSE
                rhs%data4d(i,j,k,l) = Mesh%dydzdV%data3d(i,j,k)*(this%xfluxdydz(i,j,k,l) - this%xfluxdydz(i-Mesh%IP1,j,k,l)) &
                                    + Mesh%dzdxdV%data3d(i,j,k)*(this%yfluxdzdx(i,j,k,l) - this%yfluxdzdx(i,j-Mesh%JP1,k,l)) &
                                    + Mesh%dxdydV%data3d(i,j,k)*(this%zfluxdxdy(i,j,k,l) - this%zfluxdxdy(i,j,k-Mesh%KP1,l))
              END IF
            END DO
          END DO
        END DO
      END DO

      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
          DO i=Mesh%IMIN,Mesh%IMAX
            wp = this%w(i,j) + Mesh%radius%bcenter(i,j,k) * Mesh%OMEGA
            rhs%data4d(i,j,k,Physics%ZMOMENTUM) = rhs%data4d(i,j,k,Physics%ZMOMENTUM) &
              - wp * rhs%data4d(i,j,k,Physics%DENSITY)
            IF(Physics%PRESSURE.GT.0) THEN
              rhs%data4d(i,j,k,Physics%ENERGY) = rhs%data4d(i,j,k,Physics%ENERGY) &
                - wp*( rhs%data4d(i,j,k,Physics%ZMOMENTUM) &
                       + 0.5 * wp* rhs%data4d(i,j,k,Physics%DENSITY))
              IF (have_potential) THEN
                rhs%data4d(i,j,k,Physics%ENERGY) = &
                rhs%data4d(i,j,k,Physics%ENERGY) - pot(i,j,k,1) * rhs%data4d(i,j,k,Physics%DENSITY)
              END IF
            END IF
          END DO
        END DO
      END DO

    CASE DEFAULT
      CALL this%Error("timedisc_base::ComputeRHS","rhstype 1 is not supported for this geometry")
    END SELECT



!NEC$ NOVECTOR
      DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ OUTERLOOP_UNROLL(8)
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
!NEC$ IVDEP
            DO i=Mesh%IMIN,Mesh%IMAX
              ! update right hand side of ODE
              rhs%data4d(i,j,k,l) = rhs%data4d(i,j,k,l) - this%geo_src%data4d(i,j,k,l) - this%src%data4d(i,j,k,l)
            END DO
          END DO
        END DO
      END DO
    CASE DEFAULT
      CALL this%Error("timedisc_base::ComputeRHS","only rhstype 0 or 1 supported")
    END SELECT
  END SUBROUTINE ComputeRHS

  !> \public compute the RHS of the spatially discretized PDE
  !!
  SUBROUTINE CheckData(this,Mesh,Physics,Fluxes,pvar,cvar,checkdatabm)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(fluxes_base),   INTENT(IN)    :: Fluxes
    CLASS(marray_compound),INTENT(INOUT):: pvar,cvar
    INTEGER,              INTENT(IN)    :: checkdatabm
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER                             :: err
#endif
    REAL                                :: val
    !------------------------------------------------------------------------!
    IF(IAND(checkdatabm,CHECK_CSOUND).NE.CHECK_NOTHING) THEN
      ! check for speed of sound
      SELECT TYPE(phys => Physics)
      CLASS IS(physics_eulerisotherm)
        val = MINVAL(phys%bccsound%data1d(:))
        IF(val.LE.0.) THEN
          ! warn now and stop after file output
          CALL this%Warning("CheckData","Illegal speed of sound value less than 0.")
          this%break = .TRUE.
        END IF
      CLASS DEFAULT
        CALL this%Warning("CheckData","check speed of sound selected, but bccsound not defined")
      END SELECT
    END IF
    IF((IAND(checkdatabm,CHECK_PMIN).NE.CHECK_NOTHING)) THEN
      ! check for non-isothermal physics with pressure defined
      SELECT TYPE(p => pvar)
      CLASS IS(statevector_euler)
        val = MINVAL(p%pressure%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX, &
                          Mesh%KMIN:Mesh%KMAX))
        IF(val.LT.this%pmin) THEN
          ! warn now and stop after file output
          CALL this%Warning("CheckData","Pressure below allowed pmin value.")
          this%break = .TRUE.
        END IF
      END SELECT
    END IF

    IF(IAND(checkdatabm,CHECK_RHOMIN).NE.CHECK_NOTHING) THEN
      ! check for physics with density defined
      SELECT TYPE(p => pvar)
      CLASS IS(statevector_eulerisotherm)
        val = MINVAL(p%density%data3d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX, &
                          Mesh%KMIN:Mesh%KMAX))
        IF(val.LT.this%rhomin) THEN
            ! warn now and stop after file output
            CALL this%Warning("CheckData","Density below allowed rhomin value.")
            this%break = .TRUE.
        END IF
      END SELECT
    END IF

    IF(IAND(checkdatabm,CHECK_INVALID).NE.CHECK_NOTHING) THEN
      IF(ANY((cvar%data1d.NE.cvar%data1d).OR.(pvar%data1d.NE.pvar%data1d))) THEN
        ! warn now and stop after file output
        CALL this%Warning("CheckData","Found NaN in pvar or cvar.")
        this%break = .TRUE.
      END IF
      IF(ANY((cvar%data1d.GT.HUGE(cvar%data1d)).OR.pvar%data1d.GT.HUGE(pvar%data1d))) THEN
        ! warn now and stop after file output
        CALL this%Warning("CheckData","Found Infinity in pvar or cvar.")
        this%break = .TRUE.
      END IF
    END IF

#ifdef PARALLEL
  CALL MPI_AllReduce(MPI_IN_PLACE, this%break, 1, MPI_LOGICAL, MPI_LOR, &
                             Mesh%comm_cart, err)
#endif
  END SUBROUTINE CheckData

  PURE FUNCTION GetOrder(this) RESULT(odr)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(IN) :: this
    INTEGER                          :: odr
    !------------------------------------------------------------------------!
    odr = this%order
  END FUNCTION GetOrder

  PURE FUNCTION GetCFL(this) RESULT(cfl)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(IN) :: this
    REAL                             :: cfl
    !------------------------------------------------------------------------!
    cfl = this%CFL
  END FUNCTION GetCFL


  !> \public Calculates the linear transport step in Fargo Scheme along x-axis \cite mignone2012 .
  !!
  !! The linear transport calculation is done by operator splitting. This
  !! yields an equation
  !! \f[
  !!    \frac{\partial q}{\partial t} + \mathbf{w}\cdot \nabla q = 0,
  !! \f]
  !! which can be descritized as
  !! \f[
  !!    q_j^{n+1}=q_j^{n}+ \frac{1}{\Delta y} \left(
  !!            \int_{y_{j+\frac{1}{2}}-w\Delta t}^{y_{j+\frac{1}{2}}}
  !!            q_j^{n}(\xi)\mathrm{d}\xi -
  !!            \int_{y_{j-\frac{1}{2}}-w\Delta t}^{y_{j-\frac{1}{2}}}
  !!            q_j^{n}(\xi)\mathrm{d}\xi \right).
  !! \f]
  !!
  !! FargoAdvectionX does the advection step along the x-axis with
  !! \f$ \mathbf{w} = w \mathbf{e_x} \f$.
  SUBROUTINE FargoAdvectionX(this,Fluxes,Mesh,Physics,Sources)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k,l
#ifdef PARALLEL
    CHARACTER(LEN=80)    :: str
    INTEGER              :: status(MPI_STATUS_SIZE)
    INTEGER              :: ierror
    REAL                 :: mpi_buf(2*Mesh%GNUM)
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX) :: buf
#endif
    !------------------------------------------------------------------------!
    ! determine step size of integer shift and length of remaining transport step
    ! first compute the whole step
    this%delxy(:,:)  = this%w(:,:) * this%dt / Mesh%dlx%data3d(Mesh%IMIN,:,:)

#ifdef PARALLEL
    ! make sure all MPI processes use the same step if domain is decomposed
    ! along the x-direction (can be different due to round-off errors)
    IF (Mesh%dims(1).GT.1) THEN
      CALL MPI_Allreduce(MPI_IN_PLACE,this%delxy,(Mesh%JGMAX-Mesh%JGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1), &
                         DEFAULT_MPI_REAL,MPI_MIN,Mesh%Icomm,ierror)
    END IF
#endif

    ! then subdivide into integer shift and remaining linear advection step
!NEC$ COLLAPSE
    DO k = Mesh%KGMIN,Mesh%KGMAX
!NEC$ IVDEP
      DO j = Mesh%JGMIN,Mesh%JGMAX
        this%shift(j,k) = NINT(this%delxy(j,k))
        this%delxy(j,k)  = this%delxy(j,k)-DBLE(this%shift(j,k))
      END DO
    END DO

!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM+Physics%PNUM
      this%dq%data3d(Mesh%IGMIN,:,:) = this%cvar%data4d(Mesh%IGMIN+1,:,:,l) - this%cvar%data4d(Mesh%IGMIN,:,:,l)
!NEC$ IVDEP
      DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
        this%dq%data3d(i,:,:) = this%cvar%data4d(i+1,:,:,l)-this%cvar%data4d(i,:,:,l)
        ! apply minmod limiter
        WHERE (SIGN(1.0,this%dq%data3d(i-1,:,:))*SIGN(1.0,this%dq%data3d(i,:,:)).GT.0.)
          this%dql%data3d(i,:,:) = SIGN(MIN(ABS(this%dq%data3d(i-1,:,:)),ABS(this%dq%data3d(i,:,:))),this%dq%data3d(i-1,:,:))
        ELSEWHERE
          this%dql%data3d(i,:,:) = 0.
        END WHERE
      END DO
!NEC$ IVDEP
      DO i=Mesh%IMIN-1,Mesh%IMAX
        WHERE(this%delxy(:,:).GT.0.)
          this%flux%data3d(i,:,:) = this%cvar%data4d(i,:,:,l) + .5 * this%dql%data3d(i,:,:) * (1. - this%delxy(:,:))
        ELSEWHERE
          this%flux%data3d(i,:,:) = this%cvar%data4d(i+1,:,:,l) - .5*this%dql%data3d(i+1,:,:)*(1. + this%delxy(:,:))
        END WHERE
        this%cvar%data4d(i,:,:,l) = this%cvar%data4d(i,:,:,l) - this%delxy(:,:)*(this%flux%data3d(i,:,:) - this%flux%data3d(i-1,:,:))
      END DO
    END DO

!#ifdef PARALLEL
!    ! We only need to do something, if we (also) are dealing with domain decomposition in
!    ! the second (phi) direction
!    IF(PAR_DIMS.GT.1) THEN
!        DO i=GMIN1,GMAX1
!          IF(this%shift(i,k).GT.0) THEN
!            DO l=1,Physics%VNUM+Physics%PNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:this%shift(i)) = this%cvar%data4d(MAX2-this%shift(i)+1:MAX2,i,k)
!              ELSE
!                this%buf(k,1:this%shift(i)) = this%cvar%data4d(i,MAX2-this%shift(i)+1:MAX2,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(EAST), i+Mesh%GNUM, &
!              Mesh%neighbor(WEST), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%cvar%data4d(MAX2-this%shift(i)+1:MAX2,i,k) = this%buf(k,1:this%shift(i))
!              ELSE
!                this%cvar%data4d(i,MAX2-this%shift(i)+1:MAX2,k) = this%buf(k,1:this%shift(i))
!              END IF
!            END DO
!          ELSE IF(this%shift(i).LT.0) THEN
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:-this%shift(i)) = this%cvar%data4d(MIN2:MIN2-this%shift(i)-1,i,k)
!              ELSE
!                this%buf(k,1:-this%shift(i)) = this%cvar%data4d(i,MIN2:MIN2-this%shift(i)-1,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              -this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(EAST), i+Mesh%GNUM, &
!              Mesh%neighbor(WEST), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF (Mesh%SN_shear) THEN
!                this%cvar%data4d(MIN2:MIN2-this%shift(i)-1,i,k) = this%buf(k,1:-this%shift(i))
!              ELSE
!                this%cvar%data4d(i,MIN2:MIN2-this%shift(i)-1,k) = this%buf(k,1:-this%shift(i))
!              END IF
!            END DO
!          END IF
!        END DO
!      END DO
!    END IF
!#endif

    ! Integer shift along x- or y-direction
!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ IVDEP
      DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ IVDEP
        DO j=Mesh%JGMIN,Mesh%JGMAX
          this%cvar%data4d(Mesh%IMIN:Mesh%IMAX,j,k,l) &
            = CSHIFT(this%cvar%data4d(Mesh%IMIN:Mesh%IMAX,j,k,l),-this%shift(j,k),1)
        END DO
      END DO
    END DO

  END SUBROUTINE FargoAdvectionX


  !> \public Calculates the linear transport step in Fargo Scheme along y-axis \cite mignone2012 .
  !!
  !! The linear transport calculation is done by operator splitting. This
  !! yields an equation
  !! \f[
  !!    \frac{\partial q}{\partial t} + \mathbf{w}\cdot \nabla q = 0,
  !! \f]
  !! which can be descritized as
  !! \f[
  !!    q_j^{n+1}=q_j^{n}+ \frac{1}{\Delta y} \left(
  !!            \int_{y_{j+\frac{1}{2}}-w\Delta t}^{y_{j+\frac{1}{2}}}
  !!            q_j^{n}(\xi)\mathrm{d}\xi -
  !!            \int_{y_{j-\frac{1}{2}}-w\Delta t}^{y_{j-\frac{1}{2}}}
  !!            q_j^{n}(\xi)\mathrm{d}\xi \right).
  !! \f]
  !!
  !! FargoAdvectionY does the advection step along the x-axis with
  !! \f$ \mathbf{w} = w \mathbf{e_y} \f$.
  SUBROUTINE FargoAdvectionY(this,Fluxes,Mesh,Physics,Sources)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k,l
#ifdef PARALLEL
    CHARACTER(LEN=80)    :: str
    INTEGER              :: status(MPI_STATUS_SIZE)
    INTEGER              :: ierror
    REAL                 :: mpi_buf(2*Mesh%GNUM)
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX) :: buf
#endif
    !------------------------------------------------------------------------!
    ! determine step size of integer shift and length of remaining transport step
    ! first compute the whole step
    this%delxy(:,:)  = this%w(:,:) * this%dt / Mesh%dly%data3d(:,Mesh%JMIN,:)

#ifdef PARALLEL
    ! make sure all MPI processes use the same step if domain is decomposed
    ! along the y-direction (can be different due to round-off errors)
    IF (Mesh%dims(2).GT.1) THEN
      CALL MPI_Allreduce(MPI_IN_PLACE,this%delxy,(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%KGMAX-Mesh%KGMIN+1), &
                         DEFAULT_MPI_REAL,MPI_MIN,Mesh%Jcomm,ierror)
    END IF
#endif

    ! then subdivide into integer shift and remaining linear advection step
!NEC$ COLLAPSE
    DO k = Mesh%KGMIN,Mesh%KGMAX
!NEC$ IVDEP
      DO i = Mesh%IGMIN,Mesh%IGMAX
        this%shift(i,k) = NINT(this%delxy(i,k))
        this%delxy(i,k)  = this%delxy(i,k)-DBLE(this%shift(i,k))
      END DO
    END DO

!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM+Physics%PNUM
      this%dq%data3d(:,Mesh%JGMIN,:) = this%cvar%data4d(:,Mesh%JGMIN+1,:,l) - this%cvar%data4d(:,Mesh%JGMIN,:,l)
!NEC$ IVDEP
      DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
        this%dq%data3d(:,j,:) = this%cvar%data4d(:,j+1,:,l)-this%cvar%data4d(:,j,:,l)
        ! apply minmod limiter
        WHERE (SIGN(1.0,this%dq%data3d(:,j-1,:))*SIGN(1.0,this%dq%data3d(:,j,:)).GT.0.)
          this%dql%data3d(:,j,:) = SIGN(MIN(ABS(this%dq%data3d(:,j-1,:)),ABS(this%dq%data3d(:,j,:))),this%dq%data3d(:,j-1,:))
        ELSEWHERE
          this%dql%data3d(:,j,:) = 0.
        END WHERE
      END DO
!NEC$ IVDEP
      DO j=Mesh%JMIN-1,Mesh%JMAX
        WHERE(this%delxy(:,:).GT.0.)
          this%flux%data3d(:,j,:) = this%cvar%data4d(:,j,:,l) + .5 * this%dql%data3d(:,j,:) * (1. - this%delxy(:,:))
        ELSEWHERE
          this%flux%data3d(:,j,:) = this%cvar%data4d(:,j+1,:,l) - .5*this%dql%data3d(:,j+1,:)*(1. + this%delxy(:,:))
        END WHERE
        this%cvar%data4d(:,j,:,l) = this%cvar%data4d(:,j,:,l) - this%delxy(:,:)*(this%flux%data3d(:,j,:) - this%flux%data3d(:,j-1,:))
      END DO
    END DO


!#ifdef PARALLEL
!    ! We only need to do something, if we (also) are dealing with domain decomposition in
!    ! the second (phi) direction
!    IF(PAR_DIMS.GT.1) THEN
!        DO i=GMIN1,GMAX1
!          IF(this%shift(i,k).GT.0) THEN
!            DO l=1,Physics%VNUM+Physics%PNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:this%shift(i)) = this%cvar%data4d(MAX2-this%shift(i)+1:MAX2,i,k)
!              ELSE
!                this%buf(k,1:this%shift(i)) = this%cvar%data4d(i,MAX2-this%shift(i)+1:MAX2,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(EAST), i+Mesh%GNUM, &
!              Mesh%neighbor(WEST), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%cvar%data4d(MAX2-this%shift(i)+1:MAX2,i,k) = this%buf(k,1:this%shift(i))
!              ELSE
!                this%cvar%data4d(i,MAX2-this%shift(i)+1:MAX2,k) = this%buf(k,1:this%shift(i))
!              END IF
!            END DO
!          ELSE IF(this%shift(i).LT.0) THEN
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:-this%shift(i)) = this%cvar%data4d(MIN2:MIN2-this%shift(i)-1,i,k)
!              ELSE
!                this%buf(k,1:-this%shift(i)) = this%cvar%data4d(i,MIN2:MIN2-this%shift(i)-1,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              -this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(EAST), i+Mesh%GNUM, &
!              Mesh%neighbor(WEST), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF (Mesh%SN_shear) THEN
!                this%cvar%data4d(MIN2:MIN2-this%shift(i)-1,i,k) = this%buf(k,1:-this%shift(i))
!              ELSE
!                this%cvar%data4d(i,MIN2:MIN2-this%shift(i)-1,k) = this%buf(k,1:-this%shift(i))
!              END IF
!            END DO
!          END IF
!        END DO
!      END DO
!    END IF
!#endif

    ! Integer shift along x- or y-direction
!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ IVDEP
      DO k=Mesh%KGMIN,Mesh%KGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          this%cvar%data4d(i,Mesh%JMIN:Mesh%JMAX,k,l) &
            = CSHIFT(this%cvar%data4d(i,Mesh%JMIN:Mesh%JMAX,k,l),-this%shift(i,k),1)
        END DO
      END DO
    END DO

  END SUBROUTINE FargoAdvectionY

  !> \public Calculates the linear transport step in Fargo Scheme along z-axis \cite mignone2012 .
  !!
  !! The linear transport calculation is done by operator splitting. This
  !! yields an equation
  !! \f[
  !!    \frac{\partial q}{\partial t} + \mathbf{w}\cdot \nabla q = 0,
  !! \f]
  !! which can be descritized as
  !! \f[
  !!    q_j^{n+1}=q_j^{n}+ \frac{1}{\Delta y} \left(
  !!            \int_{y_{j+\frac{1}{2}}-w\Delta t}^{y_{j+\frac{1}{2}}}
  !!            q_j^{n}(\xi)\mathrm{d}\xi -
  !!            \int_{y_{j-\frac{1}{2}}-w\Delta t}^{y_{j-\frac{1}{2}}}
  !!            q_j^{n}(\xi)\mathrm{d}\xi \right).
  !! \f]
  !!
  !! FargoAdvectionY does the advection step along the x-axis with
  !! \f$ \mathbf{w} = w \mathbf{e_y} \f$.
  SUBROUTINE FargoAdvectionZ(this,Fluxes,Mesh,Physics,Sources)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    !------------------------------------------------------------------------!
    INTEGER              :: i,j,k,l
#ifdef PARALLEL
    CHARACTER(LEN=80)    :: str
    INTEGER              :: status(MPI_STATUS_SIZE)
    INTEGER              :: ierror
    REAL                 :: mpi_buf(2*Mesh%GNUM)
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX) :: buf
#endif
    !------------------------------------------------------------------------!
    ! determine step size of integer shift and length of remaining transport step
    ! first compute the whole step
    this%delxy(:,:)  = this%w(:,:) * this%dt / Mesh%dlz%data3d(:,:,Mesh%KMIN)

#ifdef PARALLEL
    ! make sure all MPI processes use the same step if domain is decomposed
    ! along the z-direction (can be different due to round-off errors)
    IF (Mesh%dims(3).GT.1) THEN
      CALL MPI_Allreduce(MPI_IN_PLACE,this%delxy,(Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1), &
                         DEFAULT_MPI_REAL,MPI_MIN,Mesh%Kcomm,ierror)
    END IF
#endif

    ! then subdivide into integer shift and remaining linear advection step
!NEC$ COLLAPSE
    DO j = Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
      DO i = Mesh%IGMIN,Mesh%IGMAX
        this%shift(i,j) = NINT(this%delxy(i,j))
        this%delxy(i,j)  = this%delxy(i,j)-DBLE(this%shift(i,j))
      END DO
    END DO
! PRINT *,MINVAL(this%w),MAXVAL(this%w),MINVAL(this%shift),MAXVAL(this%shift),MINVAL(this%delxy),MAXVAL(this%delxy)

!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM+Physics%PNUM
      this%dq%data3d(:,:,Mesh%KGMIN) = this%cvar%data4d(:,:,Mesh%KGMIN+1,l) - this%cvar%data4d(:,:,Mesh%KGMIN,l)
!NEC$ IVDEP
      DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
        this%dq%data3d(:,:,k) = this%cvar%data4d(:,:,k+1,l)-this%cvar%data4d(:,:,k,l)
        ! apply minmod limiter
        WHERE (this%dq%data3d(:,:,k-1)*this%dq%data3d(:,:,k).LE.0.)
          ! different left and right sign -> no gradient
          this%dql%data3d(:,:,k) = 0.
        ELSEWHERE
          ! both positive/negative -> take the minium of the absolute value
          this%dql%data3d(:,:,k) = SIGN(MIN(ABS(this%dq%data3d(:,:,k-1)),ABS(this%dq%data3d(:,:,k))),this%dq%data3d(:,:,k))
        END WHERE
!         WHERE (SIGN(1.0,this%dq%data3d(:,:,k-1))*SIGN(1.0,this%dq%data3d(:,:,k)).GT.0.)
!           this%dql%data3d(:,:,k) = SIGN(MIN(ABS(this%dq%data3d(:,:,k-1)),ABS(this%dq%data3d(:,:,k))),this%dq%data3d(:,:,k))
!         ELSEWHERE
!           this%dql%data3d(:,:,k) = 0.
!         END WHERE
      END DO
! PRINT '(2(ES14.5,3(I4)))',MINVAL(this%dql%data3d),MINLOC(this%dql%data3d),MAXVAL(this%dql%data3d),MAXLOC(this%dql%data3d)
!NEC$ IVDEP
      DO k=Mesh%KMIN-1,Mesh%KMAX
!         IF (k.EQ.Mesh%KMAX) THEN
!           this%flux%data3d(:,:,k) = this%flux%data3d(:,:,Mesh%KMIN-1)
!         ELSE
          WHERE(this%delxy(:,:).GT.0.)
            this%flux%data3d(:,:,k) = this%cvar%data4d(:,:,k,l) + .5 * this%dql%data3d(:,:,k) * (1. - this%delxy(:,:))
          ELSEWHERE
            this%flux%data3d(:,:,k) = this%cvar%data4d(:,:,k+1,l) - .5*this%dql%data3d(:,:,k+1)*(1. + this%delxy(:,:))
          END WHERE
!         END IF
        IF (k.LT.Mesh%KMIN) CONTINUE
        this%cvar%data4d(:,:,k,l) = this%cvar%data4d(:,:,k,l) - this%delxy(:,:)*(this%flux%data3d(:,:,k) - this%flux%data3d(:,:,k-1))
! IF (l.GT.1) STOP
! PRINT *,k,this%flux%data3d(Mesh%IMIN,Mesh%JMIN,k),this%dq%data3d(Mesh%IMIN,Mesh%JMIN,k),this%dql%data3d(Mesh%IMIN,Mesh%JMIN,k)
      END DO
    END DO
!#ifdef PARALLEL
!    ! We only need to do something, if we (also) are dealing with domain decomposition in
!    ! the second (phi) direction
!    IF(PAR_DIMS.GT.1) THEN
!        DO i=GMIN1,GMAX1
!          IF(this%shift(i,k).GT.0) THEN
!            DO l=1,Physics%VNUM+Physics%PNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:this%shift(i)) = this%cvar%data4d(MAX2-this%shift(i)+1:MAX2,i,k)
!              ELSE
!                this%buf(k,1:this%shift(i)) = this%cvar%data4d(i,MAX2-this%shift(i)+1:MAX2,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(EAST), i+Mesh%GNUM, &
!              Mesh%neighbor(WEST), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%cvar%data4d(MAX2-this%shift(i)+1:MAX2,i,k) = this%buf(k,1:this%shift(i))
!              ELSE
!                this%cvar%data4d(i,MAX2-this%shift(i)+1:MAX2,k) = this%buf(k,1:this%shift(i))
!              END IF
!            END DO
!          ELSE IF(this%shift(i).LT.0) THEN
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:-this%shift(i)) = this%cvar%data4d(MIN2:MIN2-this%shift(i)-1,i,k)
!              ELSE
!                this%buf(k,1:-this%shift(i)) = this%cvar%data4d(i,MIN2:MIN2-this%shift(i)-1,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              -this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(EAST), i+Mesh%GNUM, &
!              Mesh%neighbor(WEST), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF (Mesh%SN_shear) THEN
!                this%cvar%data4d(MIN2:MIN2-this%shift(i)-1,i,k) = this%buf(k,1:-this%shift(i))
!              ELSE
!                this%cvar%data4d(i,MIN2:MIN2-this%shift(i)-1,k) = this%buf(k,1:-this%shift(i))
!              END IF
!            END DO
!          END IF
!        END DO
!      END DO
!    END IF
!#endif

    ! Integer shift along x- or y-direction
!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM+Physics%PNUM
!NEC$ IVDEP
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          this%cvar%data4d(i,j,Mesh%KMIN:Mesh%KMAX,l) &
            = CSHIFT(this%cvar%data4d(i,j,Mesh%KMIN:Mesh%KMAX,l),-this%shift(i,j),1)
        END DO
      END DO
    END DO

  END SUBROUTINE FargoAdvectionZ


  !> \public Calculates new background velocity for fargo advection
  !!
  !! \attention Only works when velocity is shifted in second direction.
  SUBROUTINE CalcBackgroundVelocity(this,Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(marray_compound), INTENT(INOUT) :: pvar,cvar
    !------------------------------------------------------------------------!
    REAL              :: wi
    INTEGER           :: i,j,k,v_idx
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    !------------------------------------------------------------------------!
    ! currently only eulerisotherm and derived physics are supported
    ! if fargo is enabled (see InitTimedisc)
    SELECT TYPE(p => pvar)
    CLASS IS(statevector_eulerisotherm)
      v_idx = Mesh%fargo%GetDirection() ! index of fargo transport velocity corresponds with transport direction
      SELECT CASE(Mesh%fargo%GetDirection())
      CASE(1)
        ! transform back to real, i.e. not co-moving, quantities
        CALL Physics%AddBackgroundVelocityX(Mesh,this%w,pvar,cvar)
        DO k=Mesh%KGMIN,Mesh%KGMAX
          DO j=Mesh%JGMIN,Mesh%JGMAX
            ! some up all xvelocities along the x-direction
            wi = SUM(p%velocity%data4d(Mesh%IMIN:Mesh%IMAX,j,k,v_idx)*this%wfactor(j,k))
#ifdef PARALLEL
            ! extend the sum over all partitions
            IF(Mesh%dims(1).GT.1) THEN
              CALL MPI_AllReduce(MPI_IN_PLACE, wi, 1, DEFAULT_MPI_REAL, MPI_SUM, &
                                Mesh%Icomm, ierror)
            END IF
#endif
            ! set new background velocity to the arithmetic mean of the
            ! xvelocity field along the x-direction
            this%w(j,k) = wi / Mesh%INUM
          END DO
        END DO
      CASE(2)
        ! transform back to real, i.e. not co-moving, quantities
        CALL Physics%AddBackgroundVelocityY(Mesh,this%w,pvar,cvar)
        ! if x-velocity is suppressed, i.e. zero bit not set -> first velocity component is y-velocity
        IF (.NOT.BTEST(Mesh%VECTOR_COMPONENTS,0)) v_idx = 1
        DO k=Mesh%KGMIN,Mesh%KGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! some up all yvelocities along the y-direction
            wi = SUM(p%velocity%data4d(i,Mesh%JMIN:Mesh%JMAX,k,v_idx)*this%wfactor(i,k))
#ifdef PARALLEL
            ! extend the sum over all partitions
            IF(Mesh%dims(2).GT.1) THEN
              CALL MPI_AllReduce(MPI_IN_PLACE, wi, 1, DEFAULT_MPI_REAL, MPI_SUM, &
                                  Mesh%Jcomm, ierror)
            END IF
#endif
            this%w(i,k) = wi / Mesh%JNUM
          END DO
        END DO
      CASE(3)
        ! transform back to real, i.e. not co-moving, quantities
        CALL Physics%AddBackgroundVelocityZ(Mesh,this%w,pvar,cvar)
        ! last velocity component should be the z-velocity
        v_idx = Physics%VDIM ! could be < 3
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            ! some up all zvelocities along the z-direction
            wi = SUM(p%velocity%data4d(i,j,Mesh%KMIN:Mesh%KMAX,v_idx)*this%wfactor(i,j))
#ifdef PARALLEL
            ! extend the sum over all partitions
            IF(Mesh%dims(3).GT.1) THEN
              CALL MPI_AllReduce(MPI_IN_PLACE, wi, 1, DEFAULT_MPI_REAL, MPI_SUM, &
                                  Mesh%Kcomm, ierror)
            END IF
#endif
            ! set new background velocity to the arithmetic mean of the
            ! yvelocity field along the y-direction
            this%w(i,j) = wi / Mesh%KNUM
          END DO
        END DO
      END SELECT
    CLASS DEFAULT
      CALL this%Error("timedisc_base::CalcBackgroundVelocity","physics currently not supported with fargo transport")
    END SELECT
  END SUBROUTINE CalcBackgroundVelocity


  !> Compute velocity leading to a centrifugal acceleration with
  !! respect to some given axis of rotation which balances any
  !! local radial acceleration. If the acceleration is not given
  !! explicitly determine it from a complete evaluation of
  !! of the right hand side of the transport problem.
  !!
  !! The algorithm follows these steps:
  !!
  !! 1. Get curvilinear components of the vectors pointing
  !!    from the origin/center of rotation into each cell.
  !! 2. Get curvilinear components of the vectors pointing
  !!    from the axis of rotation perpendicular to that axis
  !!    into each cell, i. e.
  !!      \f[ \vec{R} = \vec{r} - (\hat{e}_\Omega\cdot\vec{r})\hat{e}_\Omega \f].
  !!    and the components of the azimuthal unit vector in terms
  !!    of the local orthnormal basis
  !!      \f[ \hat{e}_\varphi = \hat{e}_\Omega \times \vec{r} / r \f]
  !! 3. Determine the acceleration to balance.
  !! 4. Compute absolute value of azimuthal velocity using the balance law
  !!    \f[ \Omega^2 R\hat{e}_R = -\vec{a} \Leftrightarrow
  !!        v_\varphi = \sqrt{\Omega^2 R^2} = \sqrt{-R \hat{e}_R\cdot\vec{a}} \f]
  !! 5. Compute components of \f$ v_\varphi \hat{e}_\varphi \f$ with respect
  !!    to the local orthonormal basis.
  !!
  !! \todo move code depending on physics into appropriate subroutines
  !! in physics modules
  FUNCTION GetCentrifugalVelocity(this,Mesh,Physics,Fluxes,Sources,&
                       dir_omega_,accel_,centrot) RESULT(velocity)
    USE physics_euler_mod, ONLY: physics_euler
    USE physics_eulerisotherm_mod, ONLY: physics_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(sources_list), ALLOCATABLE, INTENT(INOUT) :: Sources
    REAL, OPTIONAL,       INTENT(IN)    :: dir_omega_(3)
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VDIM), &
                OPTIONAL, INTENT(IN)    :: accel_
    REAL, OPTIONAL,       INTENT(IN)    :: centrot(3)
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VDIM) &
                                        :: velocity
    !------------------------------------------------------------------------!
    CLASS(marray_base), ALLOCATABLE &
                       :: accel,dist_axis_projected,ephi_projected,eomega,tmpvec, &
                          posvec,tmp
    REAL               :: omega2,dir_omega(3)
    INTEGER            :: k,l,err
    !------------------------------------------------------------------------!
    ALLOCATE(accel,dist_axis_projected,ephi_projected,posvec,eomega,tmpvec,tmp, &
             STAT=err)
    IF (err.NE.0) CALL this%Error("timedisc_base::GetCentrifugalVelocity", &
                       "Unable to allocate memory.")
    ! initialize marrays
    tmp    = marray_base()
    tmpvec = marray_base(3)
    eomega = marray_base(3)
    posvec = marray_base(3)
    dist_axis_projected = marray_base(Physics%VDIM)
    ephi_projected      = marray_base(Physics%VDIM)
    accel               = marray_base(Physics%VDIM)

    ! do some sanity checks on the input data
    IF (PRESENT(dir_omega_)) THEN
      dir_omega(:) = dir_omega_(:)
      omega2 = SUM(dir_omega(:)*dir_omega(:))
      ! omega must not be the zero vector
      IF (.NOT. (omega2.GT.0.0)) &
        CALL this%Error("timedisc_base::GetCentrifugalVelocity", &
           "omega must not be the zero vector")
      ! normalize dir_omega
      dir_omega(:) = dir_omega(:) / SQRT(omega2)
    ELSE
      ! default is rotation in a plane perpendicular to ê_z
      dir_omega(1:3) = (/ 0.0, 0.0, 1.0 /)
    END IF

    ! compute (local) curvilinear vector components of the unit
    ! vectors pointing in the direction of the rotational axis
    tmpvec%data2d(:,1) = dir_omega(1)
    tmpvec%data2d(:,2) = dir_omega(2)
    tmpvec%data2d(:,3) = dir_omega(3)
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,tmpvec%data4d,eomega%data4d)

    ! 1. Get the curvilinear components of all vectors pointing
    ! from the center of rotation into each cell using cartesian
    ! coordinates.
    ! The cartesian coordinates are the components of the position
    ! vector with respect to the cartesian standard orthonormal
    ! basis {ê_x, ê_y, ê_z}.
    IF (PRESENT(centrot)) THEN
      IF ((Mesh%rotsym.GT.0.0).AND.ANY(centrot(:).NE.0.0)) &
        CALL this%Error("timedisc_base::GetCentrifugalVelocity", &
           "if rotational symmetry is enabled center of rotation must be the origin")
      ! Translate the position vector to the center of rotation.
      DO k=1,3
        tmpvec%data4d(:,:,:,k) = Mesh%bccart(:,:,:,k) - centrot(k)
      END DO
      ! compute curvilinear components of translated position vectors
      CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,tmpvec%data4d,posvec%data4d)
    ELSE
      ! default center of rotation is the origin of the mesh
      ! -> use position vectors with respect to the origin
      posvec%data4d = Mesh%posvec%bcenter
    END IF

    ! 2. Compute vectors pointing from axis of rotation into each cell
    !    and the azimuthal unit vector ephi
    tmp%data1d(:) = SUM(posvec%data2d(:,:)*eomega%data2d(:,:),DIM=2)
    posvec%data2d(:,1) = posvec%data2d(:,1) - tmp%data1d(:)*eomega%data2d(:,1)
    posvec%data2d(:,2) = posvec%data2d(:,2) - tmp%data1d(:)*eomega%data2d(:,2)
    posvec%data2d(:,3) = posvec%data2d(:,3) - tmp%data1d(:)*eomega%data2d(:,3)

    ! distance to axis
    tmp%data1d(:) = SQRT(SUM(posvec%data2d(:,:)*posvec%data2d(:,:),DIM=2))
    ! use tmpvec for temporary storage of ephi
    CALL cross_product(eomega%data2d(:,1),eomega%data2d(:,2),eomega%data2d(:,3), &
      posvec%data2d(:,1)/tmp%data1d(:),posvec%data2d(:,2)/tmp%data1d(:), &
      posvec%data2d(:,3)/tmp%data1d(:),tmpvec%data2d(:,1),tmpvec%data2d(:,2),tmpvec%data2d(:,3))

    ! project it onto the mesh
    SELECT TYPE(phys => Physics)
    CLASS IS(physics_eulerisotherm)
      l = 1
      DO k=1,3
        ! check if vector component is used
        IF (BTEST(Mesh%VECTOR_COMPONENTS,k-1)) THEN
          dist_axis_projected%data2d(:,l) = posvec%data2d(:,k)
          ephi_projected%data2d(:,l) = tmpvec%data2d(:,k)
          l = l + 1
        END IF
      END DO
    CLASS DEFAULT
      CALL this%Error("timedisc_base::GetCentrifugalVelocity", &
        "physics not supported")
    END SELECT

    ! 3. determine the acceleration
    IF(PRESENT(accel_)) THEN
      accel%data4d(:,:,:,:) = accel_(:,:,:,:)
    ELSE
      ! Works only for centrot = (0,0,0)
      IF(PRESENT(centrot)) THEN
        IF(ANY(centrot(:).NE.0.0)) &
          CALL this%Error("timedisc_base::GetCentrifugalVelocity", &
            "You are not allowed to define centrot without accel.")
      END IF
      ! This may not work for physics with unusual conservative variables.
      ! We assume:  conservative momentum = density * velocity
      ! sanity check if someone derives new physics from the euler modules
      ! which does not have this property
      SELECT TYPE(phys => Physics)
      TYPE IS(physics_eulerisotherm)
        ! supported -> do nothing
      TYPE IS(physics_euler)
        ! supported -> do nothing
      CLASS DEFAULT
        CALL this%Error("timedisc_base::GetCentrifugalVelocity", &
          "physics not supported")
      END SELECT
      SELECT TYPE(phys => Physics)
      CLASS IS(physics_eulerisotherm)
        ! evaluate the right hand side for given (preliminary) initial data
        CALL this%ComputeRHS(Mesh,phys,Sources,Fluxes,this%time,0.0, &
                             this%pvar,this%cvar,this%checkdatabm,this%rhs)
        SELECT TYPE(pvar => this%pvar)
        CLASS IS(statevector_eulerisotherm)
          SELECT TYPE(rhs => this%rhs)
          CLASS IS(statevector_eulerisotherm)
            DO k=1,phys%VDIM
              accel%data2d(:,k) = -rhs%momentum%data2d(:,k) / pvar%density%data1d(:)
            END DO
          END SELECT
        END SELECT
      END SELECT
    END IF

    ! 4. compute absolute value of azimuthal velocity
    tmp%data1d(:) = SQRT(MAX(0.0,-SUM(accel%data2d(:,:)*dist_axis_projected%data2d(:,:),DIM=2)))

    ! 5. Compute components of the azimuthal velocity vector
    DO k=1,Physics%VDIM
      velocity(:,:,:,k) = tmp%data3d(:,:,:) * ephi_projected%data4d(:,:,:,k)
    END DO

    DEALLOCATE(accel,dist_axis_projected,ephi_projected,posvec,eomega,tmpvec,tmp)

  CONTAINS

    ELEMENTAL SUBROUTINE cross_product(ax,ay,az,bx,by,bz,aXbx,aXby,aXbz)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      REAL, INTENT(IN)  :: ax,ay,az,bx,by,bz
      REAL, INTENT(OUT) :: aXbx,aXby,aXbz
      !------------------------------------------------------------------------!
      aXbx = ay*bz - az*by
      aXby = az*bx - ax*bz
      aXbz = ax*by - ay*bx
    END SUBROUTINE cross_product

  END FUNCTION GetCentrifugalVelocity

  SUBROUTINE Finalize_base(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base)             :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
        CALL this%Error("CloseTimedisc","not initialized")
    ! call boundary destructor
    CALL this%Boundary%Finalize()

    DEALLOCATE( &
      this%pvar,this%cvar,this%ptmp,this%ctmp,this%cold, &
      this%geo_src,this%src,this%rhs, &
      this%xfluxdydz,this%yfluxdzdx,this%zfluxdxdy,this%amax,this%tol_abs,&
      this%dtmean,this%dtstddev,this%time)

    IF (ASSOCIATED(this%cerr)) DEALLOCATE(this%cerr)
    IF (ASSOCIATED(this%cerr_max)) DEALLOCATE(this%cerr_max)
    IF (ASSOCIATED(this%w)) DEALLOCATE(this%w)
    IF (ASSOCIATED(this%wfactor)) DEALLOCATE(this%wfactor)
    IF (ASSOCIATED(this%delxy))DEALLOCATE(this%delxy)
    IF (ASSOCIATED(this%shift))DEALLOCATE(this%shift)
#ifdef PARALLEL
    IF(ASSOCIATED(this%buf))  DEALLOCATE(this%buf)
#endif
    IF(ASSOCIATED(this%bflux)) DEALLOCATE(this%bflux)
    IF(ASSOCIATED(this%solution)) DEALLOCATE(this%solution)
  END SUBROUTINE Finalize_base


END MODULE timedisc_base_mod
