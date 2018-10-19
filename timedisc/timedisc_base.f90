!r#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: timedisc_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
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
!! - general parameters of timedisc group as key-values
!! \key{method,INTEGER,time integration method}
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
!! \key{always_update_bccsound,INTEGER,enable (=1) update of bccsound at each
!!      substep (0 disables),1}

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief generic subroutines for time discretization
!!
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_base_mod
  USE logging_base_mod
  USE boundary_generic_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE sources_base_mod
  USE fluxes_base_mod
  USE reconstruction_base_mod
  USE field_base_mod
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
     !> old error and used in the dumka method
     REAL                              :: ERR_N, H_N
     REAL, DIMENSION(:), POINTER       :: tol_abs          !< abs. error tolerance
     REAL                              :: beta             !< time step friction
     REAL, DIMENSION(:,:,:,:), POINTER :: pvar, cvar       !< prim/cons vars
     REAL, DIMENSION(:,:,:,:), POINTER :: solution=>Null() !< analytical solution
     REAL, DIMENSION(:,:,:,:), POINTER :: cold             !< old prim/cons vars
     REAL, DIMENSION(:,:,:,:), POINTER :: ptmp,ctmp        !< temporary cvars

     !> multistep vars
     REAL, DIMENSION(:,:,:,:), POINTER :: phi,oldphi_s,&
                                          newphi_s
     REAL, DIMENSION(:), POINTER       :: gamma
     INTEGER                           :: pc               !< = 1 predictor-corrector
     REAL, DIMENSION(:,:,:,:), POINTER :: src,geo_src      !< source terms
     REAL, DIMENSION(:,:,:,:), POINTER :: rhs              !< ODE right hand side
     !> numerical fluxes divided by dy or dx
     REAL, DIMENSION(:,:,:,:), POINTER :: xfluxdydz,yfluxdzdx,zfluxdxdy
     REAL, DIMENSION(:,:,:,:), POINTER :: amax            !< max. wave speeds
     REAL, DIMENSION(:,:), POINTER     :: bflux           !< boundary fluxes for output
     REAL, DIMENSION(:,:,:,:), POINTER :: errorval        !< max. wave speeds
     LOGICAL                           :: write_error     !< enable err writing
     INTEGER, DIMENSION(:,:), POINTER  :: shift=>null()   !< fargo annulus shift
     REAL, DIMENSION(:,:), POINTER     :: buf=>null()     !< fargo MPI buffer
     REAL, DIMENSION(:,:), POINTER     :: w=>null()       !< fargo background velocity
     REAL, DIMENSION(:,:), POINTER     :: delxy =>null()  !< fargo residual shift
     REAL, DIMENSION(:,:,:,:), POINTER :: fargo_src =>null() !< fargo source terms
     INTEGER                           :: fargo_order     !< used in FargoAdvection


!     TYPE(elem_TYP),POINTER            :: ts=>Null()       !< timesteps (multistep)
!  TYPE elem_TYP
!    TYPE(elem_TYP), POINTER            :: next,prev
!    REAL                               :: t                !< time
!  END TYPE elem_TYP
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
    PROCEDURE           :: FargoAdvection
    PROCEDURE           :: CalcBackgroundVelocity
    PROCEDURE (SolveODE), DEFERRED :: SolveODE
    PROCEDURE           :: Finalize_base
    PROCEDURE (Finalize), DEFERRED :: Finalize
  END TYPE timedisc_base
!----------------------------------------------------------------------------!
  ABSTRACT INTERFACE
    SUBROUTINE SolveODE(this,Mesh,Physics,Sources,Fluxes,time,dt,err)
      IMPORT timedisc_base, mesh_base, physics_base, fluxes_base, sources_base
      IMPLICIT NONE
      CLASS(timedisc_base), INTENT(INOUT) :: this
      CLASS(mesh_base),     INTENT(IN)    :: Mesh
      CLASS(physics_base),  INTENT(INOUT) :: Physics
      CLASS(sources_base),  POINTER       :: Sources
      CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
      REAL,                 INTENT(IN)    :: time
      REAL,                 INTENT(INOUT) :: dt, err
    END SUBROUTINE
!    REAL FUNCTION CalcTimestep(this,Mesh,Physics,Fluxes,time,dtcause) RESULT(dt)
!      IMPORT timedisc_base, mesh_base, physics_base, fluxes_base
!      IMPLICIT NONE
!      CLASS(timedisc_base), INTENT(INOUT) :: this
!      CLASS(mesh_base),     INTENT(IN)    :: Mesh
!      CLASS(physics_base),  INTENT(INOUT) :: Physics
!      CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
!      REAL,                 INTENT(IN)    :: time
!      INTEGER,              INTENT(INOUT) :: dtcause
!    END FUNCTION
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
!  INTEGER, PARAMETER :: MULTISTEP        = 5
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
       MODIFIED_EULER, RK_FEHLBERG, CASH_KARP, DORMAND_PRINCE, &! MULTISTEP, &
       SSPRK, &
       DTCAUSE_CFL,DTCAUSE_ERRADJ,DTCAUSE_SMALLERR, DTCAUSE_FILEIO, &
       CHECK_ALL, CHECK_NOTHING, CHECK_CSOUND, CHECK_RHOMIN, CHECK_PMIN, &
       CHECK_INVALID, CHECK_TMIN
  !--------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE InitTimedisc(this,Mesh,Physics,config,IO,ttype,tname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(INOUT) :: Mesh
    CLASS(physics_base),  INTENT(IN)    :: Physics
    TYPE(Dict_TYP),       POINTER       :: config,IO
    INTEGER                             :: ttype
    CHARACTER(LEN=32)                   :: tname
    !------------------------------------------------------------------------!
    INTEGER              :: err, d
    CHARACTER(LEN=32)    :: order_str,cfl_str,stoptime_str,dtmax_str,beta_str
    CHARACTER(LEN=32)    :: info_str,shear_direction
    INTEGER              :: method
    !------------------------------------------------------------------------!
    INTENT(IN)           :: ttype, tname
    !------------------------------------------------------------------------!
    CALL this%InitLogging(ttype,tname)

    IF (.NOT.Physics%Initialized().OR..NOT.Mesh%Initialized()) &
         CALL this%Error("InitTimedisc","physics and/or mesh module uninitialized")
    ! For other methods the data is stored in this%coeff and rhs points to it. Therefore an allocation is
    ! not necessary or rather leads to memory leaks.
    IF(this%GetType().EQ.MODIFIED_EULER) THEN
      ALLOCATE(this%rhs(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM))
      this%rhs  = 0.
    END IF

      ! allocate memory for data structures needed in all timedisc modules
    ALLOCATE( &
      this%pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),      &
      this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),      &
      this%cold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),      &
      this%ptmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),      &
      this%ctmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),      &
      this%geo_src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),   &
      this%src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),       &
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

    ! initialize all variables
    this%pvar      = 0.
    this%cvar      = 0.
    this%cold      = 0.
    this%ptmp      = 0.
    this%ctmp      = 0.
    this%geo_src   = 0.
    this%src       = 0.
    this%xfluxdydz = 0.
    this%yfluxdzdx = 0.
    this%zfluxdxdy = 0.
    this%amax      = 0.
    this%tol_abs   = 0.
    this%dtmean    = 0.
    this%dtstddev  = 0.
    this%time      = 0.
    this%break     = .FALSE.
    this%n_adj     = 0
    this%maxerrold = 0.
    this%dtaccept  = 0

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
    CALL GetAttr(config, "checkdata", this%checkdatabm, CHECK_INVALID)

    ! pressure minimum to check if CHECK_PMIN is active
    CALL GetAttr(config, "pmin", this%pmin, TINY(this%pmin))

    ! temperature minimum to check if CHECK_TMIN is active
    CALL GetAttr(config, "tmin", this%tmin, TINY(this%tmin))

    ! density minimum to check if CHECK_RHOMIN is active
    CALL GetAttr(config, "rhomin", this%rhomin, TINY(this%rhomin))

    CALL GetAttr(config, "output/"//TRIM(Physics%pvarname(1)), d, 1)

    CALL this%SetOutput(Mesh,Physics,config,IO)

    ! handles the order of the fargo step
    CALL GetAttr(config, "fargo_order", this%fargo_order, 2)

    ! check if fargo can be used
    SELECT CASE (Mesh%FARGO)
    CASE(0) ! fargo disabled
       ! do nothing
    CASE(1,2,3) ! fargo enabled
       ! check physics
       SELECT CASE(Physics%GetType())
       CASE(EULER2D_ISOIAMT,EULER2D_IAMT,EULER2D_ISOIAMROT,EULER2D_ISOTHERM,EULER2D)
          ! check geometry
          SELECT CASE(Mesh%Geometry%GetType())
          CASE(POLAR,TANPOLAR,LOGPOLAR,SINHPOLAR)
             ! allocate data arrays used for fargo
             ALLOCATE( &
                      this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%delxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%shift(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%fargo_src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM+Physics%PNUM), &
#ifdef PARALLEL
                      this%buf(Physics%VNUM+Physics%PNUM,1:Mesh%MINJNUM), & !!! NOT CHANGED, because not clear were used !
#endif
                      STAT = err)
             IF (err.NE.0) THEN
                CALL this%Error("InitTimedisc", "Unable to allocate memory for fargo advection.")
             END IF
             this%fargo_src(:,:,:,:) = 0.0
          CASE(CARTESIAN) ! in cartesian fargo shift can be chosen in either x- or y-direction
             IF(Mesh%WE_shear) THEN
                ALLOCATE( &
                      this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%delxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%shift(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
#ifdef PARALLEL
                      this%buf(Physics%VNUM+Physics%PNUM,1:Mesh%MINJNUM), &  !!! SEE ABOVE
#endif
                      STAT = err)
                IF (err.NE.0) THEN
                   CALL this%Error("InitTimedisc", "Unable to allocate memory for fargo advection.")
                END IF
             ELSE IF(Mesh%SN_shear) THEN
                ALLOCATE( &
                      this%w(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%delxy(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                      this%shift(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
#ifdef PARALLEL
                      this%buf(Physics%VNUM+Physics%PNUM,1:Mesh%MININUM), & !!! SEE ABOVE
#endif
                      STAT = err)
                IF (err.NE.0) THEN
                   CALL this%Error("InitTimedisc", "Unable to allocate memory for fargo advection.")
                END IF
             END IF
          CASE DEFAULT
             ! geometry not supported -> disable fargo
             Mesh%FARGO = 0
             CALL this%Warning("InitTimedisc", &
                 "fargo has been disabled, because the geometry is not supported.")
          END SELECT
       CASE DEFAULT
          ! geometry not supported -> disable fargo
          Mesh%FARGO = 0
          CALL this%Warning("InitTimedisc","fargo has been disabled, because the physics is not supported.")
       END SELECT
       ! initialize background velocity field w
       SELECT CASE(Mesh%FARGO)
       CASE(1,2) ! set to 0;
                 ! fargo advection type 1: w is computed in each time step (see FargoCalcVelocity)
                 ! fargo advection type 2: w is provided by the user, e.g. in InitData
          this%w(:,:) = 0.0
          this%shift(:,:) = 0.0
       CASE(3) ! fixed background velocity in shearing box
         IF(Mesh%WE_shear) THEN
           this%w(:,:) = -Mesh%Q*Mesh%omega*Mesh%bcenter(:,Mesh%JMIN,:,1) !-Q*Omega*x
           this%shift(:,:) = 0.0
           shear_direction = "west<->east"
         ELSE IF (Mesh%SN_shear) THEN
           this%w(:,:) = Mesh%Q*Mesh%omega*Mesh%bcenter(Mesh%IMIN,:,:,2) !Q*Omega*y
           this%shift(:,:) = 0.0
           shear_direction = "south<->north"
         END IF
       END SELECT
    CASE DEFAULT
       CALL this%Error("InitTimedisc","unknown fargo advection scheme")
    END SELECT


    ! print some information
    WRITE (order_str, '(I0)') this%GetOrder()
    WRITE (cfl_str, '(F4.2)') this%GetCFL()
    WRITE (stoptime_str, '(ES10.4)') this%stoptime
    WRITE (dtmax_str, '(ES10.4)') this%dtmax
    WRITE (beta_str, '(ES10.4)') this%beta
    CALL this%Info(" TIMEDISC-> ODE solver:        " //TRIM(this%GetName()))
    CALL this%Info("            order:             " //TRIM(order_str))
    CALL this%Info("            CFL number:        " //TRIM(cfl_str))
    CALL this%Info("            dtmax:             " //TRIM(dtmax_str))
    CALL this%Info("            stoptime:          " //TRIM(stoptime_str))
    CALL this%Info("            beta:              " //TRIM(beta_str))
    ! adaptive step size control
    IF (this%tol_rel.LT.1.0) THEN
       WRITE (info_str,'(ES7.1)') this%tol_rel*100
       CALL this%Info("            step size control: enabled")
       CALL this%Info("            rel. precision:    "//TRIM(info_str)//" %")
    ELSE
       WRITE (info_str,'(A)') "disabled"
    END IF
    IF (Mesh%FARGO.NE.0) &
       CALL this%Info("            fargo:             " //TRIM(FARGO_METHOD(Mesh%FARGO)))
    IF(Mesh%FARGO.EQ.3) &
       CALL this%Info("            shear-direction:   " //TRIM(shear_direction))
    SELECT CASE(this%rhstype)
    CASE(0)
       ! special rhs disabled, print nothing
    CASE(1)
       IF (.NOT.ASSOCIATED(this%w)) THEN
          ALLOCATE(this%w(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX),STAT = err)
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
    INTEGER                             :: valwrite,i,err
    CHARACTER(LEN=128)                  :: key,pvar_key
    LOGICAL                             :: writeSolution
    !------------------------------------------------------------------------!
    valwrite = 0
    CALL GetAttr(config, "output/error", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      this%write_error = .TRUE.
      ALLOCATE( &
        this%errorval(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
        STAT = err)
      IF (err.NE.0) &
        CALL this%Error("SetOutput_timedisc", "Unable to allocate memory.")
      this%errorval = 0.
    ELSE
      this%write_error = .FALSE.
    END IF

    valwrite = 0
    CALL GetAttr(config, "output/solution", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      writeSolution = .TRUE.
      ALLOCATE( &
        this%solution(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
        STAT = err)
      IF (err.NE.0) &
        CALL this%Error("SetOutput_timedisc", "Unable to allocate memory.")
    ELSE
      writeSolution = .False.
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
                     this%pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      IF(writeSolution) THEN
        CALL SetAttr(IO, TRIM(key)//"_solution", &
          this%solution(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      !cons
      key = TRIM(Physics%cvarname(i))
      IF(key.EQ.pvar_key) key = TRIM(key) // "_cvar"
      valwrite = 0
      CALL GetAttr(config, "output/" // TRIM(key), valwrite, 0)
      ! second argument is important if pvarname is used twice
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key), &
                      this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      key = TRIM(Physics%cvarname(i))
      IF(this%write_error) THEN
        CALL SetAttr(IO, "error_" // TRIM(key), &
                     this%errorval(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write geometrical sources
      CALL GetAttr(config, "output/" // "geometrical_sources", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_geo_src", &
                        this%geo_src(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write external sources
      CALL GetAttr(config, "output/" // "external_sources", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_src", &
                        this%src(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write right hand side
      CALL GetAttr(config, "output/" // "rhs", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_rhs", &
                        this%rhs(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END IF

      ! write fluxes
      ! ATTENTION: this are the numerical fluxes devided by dy or dx respectively
      CALL GetAttr(config, "output/" // "fluxes", valwrite, 0)
      IF (valwrite.EQ.1) THEN
           CALL SetAttr(IO, TRIM(key)//"_xfluxdy", &
                        this%xfluxdydz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
           CALL SetAttr(IO, TRIM(key)//"_yfluxdx", &
                        this%yfluxdzdx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
           CALL SetAttr(IO, TRIM(key)//"_zfluxdx", &
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
    CLASS(sources_base),  POINTER       :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    INTEGER                             :: iter
    TYPE(Dict_TYP),       POINTER       :: config,IO
    !------------------------------------------------------------------------!
    REAL                                :: err,dtold,dt,time
    !------------------------------------------------------------------------!
    INTENT(INOUT)                       :: iter
    !------------------------------------------------------------------------!
    ! transform to selenoidal velocities if fargo is enabled
    IF (Mesh%FARGO.GT.0) THEN
       IF (Mesh%SN_shear.AND.Mesh%FARGO.EQ.3) THEN
          CALL Physics%SubtractBackgroundVelocityX(Mesh,this%w,this%pvar,this%cvar)
       ELSE
          CALL Physics%SubtractBackgroundVelocityY(Mesh,this%w,this%pvar,this%cvar)
       END IF
    END IF


    time = this%time
    dt   = this%dt
    IF (dt.LT.this%dtmin .AND. this%dtcause .NE. DTCAUSE_FILEIO) THEN
      ! only save dtmin if the reasion is not the fileio
      this%dtmin = dt
      this%dtmincause = this%dtcause
    END IF
    timestep: DO WHILE (time+dt.LE.this%time+this%dt)
      dtold = dt


      CALL this%SolveODE(Mesh,Physics,Sources,Fluxes,time,dt,err)
      ! check truncation error and restart if necessary
      IF (err.LT.1.0) THEN
         CALL this%AcceptSolution(Mesh,Physics,Sources,Fluxes,time,dtold,iter)
      ELSE
         CALL this%RejectSolution(Mesh,Physics,Sources,Fluxes,time,dt)
      END IF

      IF (dt.LT.this%dtlimit) &
        this%break = .TRUE.
      ! Break if dt.LT.this%dtlimit or CheckData failed
      IF(this%break) THEN
        ! Do not attempt to fargo shift anymore
        Mesh%FARGO = 0
        EXIT timestep
      END IF
    END DO timestep

    ! Save true advanced time (for fargo linear advection)
    this%dt = time - this%time

    this%time  = time
    this%dtold = dt

     ! perform the fargo advection step if enabled
    IF (Mesh%FARGO.GT.0) THEN
      CALL this%FargoAdvection(Fluxes,Mesh,Physics,Sources)
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
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_base),  POINTER       :: Sources
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
    IF (.NOT.this%always_update_bccsound) &
       CALL Physics%UpdateSoundSpeed(Mesh,this%pvar)
    CALL Physics%CalculateWaveSpeeds(Mesh,this%pvar,Fluxes%minwav,Fluxes%maxwav)


    ! compute maximum of inverse time for CFL condition
    IF ((Mesh%JNUM.EQ.1).AND.(Mesh%KNUM.EQ.1)) THEN
       ! 1D, only x-direction
       invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,1),-Fluxes%minwav(:,:,:,1)) / Mesh%dlx(:,:,:))
    ELSE IF ((Mesh%INUM.EQ.1).AND.(Mesh%KNUM.EQ.1)) THEN
       ! 1D, only y-direction
       invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,2),-Fluxes%minwav(:,:,:,2)) / Mesh%dly(:,:,:))
    ELSE IF ((Mesh%INUM.EQ.1).AND.(Mesh%JNUM.EQ.1)) THEN
       ! 1D, only z-direction
       invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,3),-Fluxes%minwav(:,:,:,3)) / Mesh%dlz(:,:,:))
    ELSE IF ((Mesh%INUM.GT.1).AND.(Mesh%JNUM.GT.1).AND.(Mesh%KNUM.EQ.1)) THEN
       ! 2D, x-y-plane
       invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,1),-Fluxes%minwav(:,:,:,1)) / Mesh%dlx(:,:,:) &
                    + MAX(Fluxes%maxwav(:,:,:,2),-Fluxes%minwav(:,:,:,2)) / Mesh%dly(:,:,:))
    ELSE IF ((Mesh%INUM.GT.1).AND.(Mesh%KNUM.GT.1).AND.(Mesh%JNUM.EQ.1)) THEN
       ! 2D, x-z-plane
       invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,1),-Fluxes%minwav(:,:,:,1)) / Mesh%dlx(:,:,:) &
                    + MAX(Fluxes%maxwav(:,:,:,3),-Fluxes%minwav(:,:,:,3)) / Mesh%dlz(:,:,:))
    ELSE IF ((Mesh%JNUM.GT.1).AND.(Mesh%KNUM.GT.1).AND.(Mesh%INUM.EQ.1)) THEN
       ! 2D, y-z-plane
       invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,2),-Fluxes%minwav(:,:,:,2)) / Mesh%dly(:,:,:) &
                    + MAX(Fluxes%maxwav(:,:,:,3),-Fluxes%minwav(:,:,:,3)) / Mesh%dlz(:,:,:))
    ELSE
       ! full 3D
       !TODO: Achtung: Hier wurde fuer eine bessere Symmetrie fuer jede Richtung ein eigenes invdt
       ! berechnet. Dies koennte jedoch einen Verlust an Stabilitaet bewirken. Hier muesste mal eine
       ! Stabilitaetsanalyse gemacht werden
       invdt_x = MAXVAL(MAX(Fluxes%maxwav(:,:,:,1),-Fluxes%minwav(:,:,:,1)) / Mesh%dlx(:,:,:))
       invdt_y = MAXVAL(MAX(Fluxes%maxwav(:,:,:,2),-Fluxes%minwav(:,:,:,2)) / Mesh%dly(:,:,:))
       invdt_z = MAXVAL(MAX(Fluxes%maxwav(:,:,:,3),-Fluxes%minwav(:,:,:,3)) / Mesh%dlz(:,:,:))
       invdt = MAX(invdt_y,invdt_z,invdt_x)
     !  invdt = MAXVAL(MAX(Fluxes%maxwav(:,:,:,3),-Fluxes%minwav(:,:,:,3)) / Mesh%dlz(:,:,:) &
     !               + MAX(Fluxes%maxwav(:,:,:,2),-Fluxes%minwav(:,:,:,2)) / Mesh%dly(:,:,:) &
     !               + MAX(Fluxes%maxwav(:,:,:,1),-Fluxes%minwav(:,:,:,1)) / Mesh%dlx(:,:,:))
    END IF

    ! largest time step due to CFL condition
    dt_cfl = this%cfl / invdt
    ! due to cfl = 0
    dtcause = 0

    ! initialize this to be sure dt_src > 0
    dt_src = dt_cfl
    CALL Sources%CalcTimestep(Mesh,Physics,Fluxes,time,this%pvar,this%cvar,dt_src,dtcause)

    dt = MIN(dt_cfl,dt_src)
  END FUNCTION CalcTimestep

  SUBROUTINE AcceptSolution(this,Mesh,Physics,Sources,Fluxes,time,dt,iter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_base),  POINTER       :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL                                :: time,dt
    INTEGER                             :: iter
    !------------------------------------------------------------------------!
    REAL                                :: dtmeanold
    !------------------------------------------------------------------------!
    INTENT(INOUT)                       :: time,dt
    !------------------------------------------------------------------------!
    time=time+dt
    this%dtaccept = this%dtaccept + 1
    dtmeanold = this%dtmean
    this%dtmean = this%dtmean + (dt - this%dtmean)/this%dtaccept
    this%dtstddev = this%dtstddev + (dt - dtmeanold)*(dt-this%dtmean)
    CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,time,dt,this%pvar,&
      this%cvar,this%checkdatabm,this%rhs)
    this%cold(:,:,:,:) = this%cvar(:,:,:,:)
    Fluxes%bxfold(:,:,:,:) = Fluxes%bxflux(:,:,:,:)
    Fluxes%byfold(:,:,:,:) = Fluxes%byflux(:,:,:,:)
    Fluxes%bzfold(:,:,:,:) = Fluxes%bzflux(:,:,:,:)
    iter = iter + 1
  END SUBROUTINE AcceptSolution

  SUBROUTINE RejectSolution(this,Mesh,Physics,Sources,Fluxes,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_base),  POINTER       :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL                                :: time,dt
    !------------------------------------------------------------------------!
    INTENT(IN)                          :: time
    INTENT(INOUT)                       :: dt
    !------------------------------------------------------------------------!
    this%cvar(:,:,:,:) = this%cold(:,:,:,:)
    ! This data has already been checked at AcceptSolution
    CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,time,dt,this%pvar, &
      this%cvar,CHECK_NOTHING,this%rhs)
    Fluxes%bxflux(:,:,:,:) = Fluxes%bxfold(:,:,:,:)
    Fluxes%byflux(:,:,:,:) = Fluxes%byfold(:,:,:,:)
    Fluxes%bzflux(:,:,:,:) = Fluxes%bzfold(:,:,:,:)
    ! count adjustments for information
    this%n_adj = this%n_adj + 1
    this%dtcause = DTCAUSE_ERRADJ
    ! only save dtmin if the reasion is not the fileio (automatically satisfied)
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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX, Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX, Physics%VNUM)          &
                       :: cvar_high,cvar_low
    REAL               :: maxerr
    !------------------------------------------------------------------------!
    INTEGER            :: i,j,k,l
#ifdef PARALLEL
    INTEGER            :: ierror
#endif
    REAL               :: rel_err(Physics%VNUM),err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: cvar_high,cvar_low
    !------------------------------------------------------------------------!
    ! check for error output
    IF (this%write_error) THEN
       DO l=1,Physics%VNUM
          rel_err(l) = 0.
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
               DO i=Mesh%IMIN,Mesh%IMAX
                  ! estimate the error for all variables
                  err = ABS(cvar_high(i,j,k,l)-cvar_low(i,j,k,l)) &
                        / (this%tol_rel*ABS(cvar_high(i,j,k,l)) + this%tol_abs(l))
                  ! determine the global maximum on the whole grid for each variable
                  rel_err(l) = MAX(rel_err(l),err)
                  ! store the maximum error for output
                  this%errorval(i,j,k,l) = MAX(this%errorval(i,j,k,l),err)
               END DO
            END DO
          END DO
       END DO
    ELSE
       DO l=1,Physics%VNUM
          rel_err(l) = 0.
          DO k=Mesh%KMIN,Mesh%KMAX
            DO j=Mesh%JMIN,Mesh%JMAX
               DO i=Mesh%IMIN,Mesh%IMAX
                  ! estimate the error for all variables
                  err = ABS(cvar_high(i,j,k,l)-cvar_low(i,j,k,l)) &
                        / (this%tol_rel*ABS(cvar_high(i,j,k,l)) + this%tol_abs(l))
                  ! determine the global maximum on the whole grid for each variable
                  rel_err(l) = MAX(rel_err(l),err)
               END DO
            END DO
          END DO
       END DO

    END IF

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
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_base),  POINTER       :: Sources
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    REAL,                 INTENT(IN)    :: time, dt
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM+Physics%PNUM), &
                          INTENT(INOUT) :: cvar
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM+Physics%PNUM), &
                          INTENT(OUT)   :: pvar,rhs
    INTEGER,              INTENT(IN)    :: checkdatabm
    !------------------------------------------------------------------------!
    INTEGER                             :: i,j,k,l
    REAL                                :: t
    !------------------------------------------------------------------------!
    t = time
    SELECT CASE(Mesh%FARGO)
    CASE(1,2)
        ! transform to real velocity, i.e. v_residual + w_background,
        ! before setting the boundary conditions
        CALL Physics%AddBackgroundVelocityY(Mesh,this%w,pvar,cvar)
    CASE(3)
        ! boundary conditions are set for residual velocity in shearing box simulations
        ! usually it's not necessary to subtract the background velocity here, but
        ! in case someone adds it before, we subtract it here; the  subroutine checks,
        ! if the velocities have already been transformed
        IF (Mesh%SN_shear) THEN
          CALL Physics%SubtractBackgroundVelocityX(Mesh,this%w,pvar,cvar)
        ELSE IF (Mesh%WE_shear) THEN
          CALL Physics%SubtractBackgroundVelocityY(Mesh,this%w,pvar,cvar)
        END IF
        ! ATTENTION: the time must be the initial time of the whole time step
        !            not the time of a substep
        t = this%time
    CASE DEFAULT
        ! fargo disabled (do nothing)
    END SELECT

    ! set boundary values and convert conservative to primitive variables
    CALL this%boundary%CenterBoundary(Mesh,Physics,t,pvar,cvar)

    IF(IAND(checkdatabm,CHECK_TMIN).NE.CHECK_NOTHING.AND.&
      this%tmin.GT.1.E-10.AND.&
      Physics%PRESSURE.GT.0) THEN
      ! Check if the temperature is below tmin. If it is, increase the pressure
      ! to reach tmin
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
          DO i=Mesh%IGMIN,Mesh%IGMAX
            pvar(i,j,k,Physics%PRESSURE) &
              = MAX(pvar(i,j,k,Physics%PRESSURE), &
                    pvar(i,j,k,Physics%DENSITY)*Physics%Constants%RG/Physics%MU*this%TMIN)
          END DO
        END DO
      END DO
      CALL Physics%Convert2Conservative(Mesh,pvar,cvar)
    END IF

    ! update the speed of sound
    IF (this%always_update_bccsound) CALL Physics%UpdateSoundSpeed(Mesh,pvar)

    ! check for illegal data
    IF(checkdatabm.NE.CHECK_NOTHING) &
      CALL this%CheckData(Mesh,Physics,Fluxes,pvar,cvar,checkdatabm)

    ! get geometrical sources
    CALL Physics%GeometricalSources(Mesh,pvar,cvar,this%geo_src)

    ! get source terms due to external forces if present
    IF (ASSOCIATED(Sources)) &
       CALL Sources%ExternalSources(Mesh,Fluxes,Physics, &
            time,dt,pvar,cvar,this%src)

    ! if fargo advection is enabled additional source terms occur;
    ! furthermore computation of numerical fluxes should always be
    ! carried out with residual velocity
    SELECT CASE(Mesh%FARGO)
    CASE(1,2)
        ! compute fargo source terms
        CALL Physics%FargoSources(Mesh,this%w,pvar,cvar,this%fargo_src)
        ! add them to the geometrical source terms
        this%geo_src(:,:,:,:) = this%geo_src(:,:,:,:) + this%fargo_src(:,:,:,:)
        ! subtract background velocity
        CALL Physics%SubtractBackgroundVelocityY(Mesh,this%w,pvar,cvar)
    CASE(3)
        ! background velocity field has already been subtracted (do nothing);
        ! fargo specific source terms are handled in the shearing box source
        ! term module
    CASE DEFAULT
        ! fargo disabled (do nothing)
    END SELECT

    ! get the numerical fluxes
    CALL Fluxes%CalculateFluxes(Mesh,Physics,pvar,cvar,this%xfluxdydz,this%yfluxdzdx,this%zfluxdxdy)

    ! compute the right hand side for boundary flux computation;
    ! this is probably wrong for the special rhs (rhstype=1, see above)
    DO l=1,Physics%VNUM+Physics%PNUM
      DO k=Mesh%KMIN,Mesh%KMAX
        DO j=Mesh%JMIN,Mesh%JMAX
          ! western and eastern boundary fluxes
          rhs(Mesh%IMIN-Mesh%Ip1,j,k,l) = Mesh%dy*Mesh%dz * this%xfluxdydz(Mesh%IMIN-Mesh%Ip1,j,k,l)
          rhs(Mesh%IMAX+Mesh%Ip1,j,k,l) = -Mesh%dy*Mesh%dz * this%xfluxdydz(Mesh%IMAX,j,k,l)
        END DO
      END DO
      DO k=Mesh%KMIN,Mesh%KMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          ! southern and northern boundary fluxes
          rhs(i,Mesh%JMIN-Mesh%Jp1,k,l) = Mesh%dz*Mesh%dx * this%yfluxdzdx(i,Mesh%JMIN-Mesh%Jp1,k,l)
          rhs(i,Mesh%JMAX+Mesh%Jp1,k,l) = -Mesh%dz*Mesh%dx * this%yfluxdzdx(i,Mesh%JMAX,k,l)
        END DO
      END DO
      DO j=Mesh%JMIN,Mesh%JMAX
        DO i=Mesh%IMIN,Mesh%IMAX
          ! bottomer and upper boundary fluxes
          rhs(i,j,Mesh%KMIN-Mesh%Kp1,l) = Mesh%dx*Mesh%dy * this%zfluxdxdy(i,j,Mesh%KMIN-Mesh%Kp1,l)
          rhs(i,j,Mesh%KMAX+Mesh%Kp1,l) = -Mesh%dx*Mesh%dy * this%zfluxdxdy(i,j,Mesh%KMAX,l)
        END DO
      END DO
    END DO

!    !fixme: Very hacky implementation of pluto style angular momentum
!    !conservation. A better implementation is needed, but we would need to
!    !restucture the timedisc module
!    !also: use this with EULER2DISOTHM. EULER2D version is untested, may need
!    !small fixes.
    SELECT CASE(this%rhstype)
    CASE(0)
!NEC$ UNROLL(8)
      DO l=1,Physics%VNUM ! to be deleted later
        DO k=Mesh%KMIN,Mesh%KMAX
          DO j=Mesh%JMIN,Mesh%JMAX
            DO i=Mesh%IMIN,Mesh%IMAX
              ! update right hand side of ODE
              rhs(i,j,k,l) = &
                    Mesh%dydzdV(i,j,k)*(this%xfluxdydz(i,j,k,l) - this%xfluxdydz(i-Mesh%Ip1,j,k,l)) &
                  + Mesh%dzdxdV(i,j,k)*(this%yfluxdzdx(i,j,k,l) - this%yfluxdzdx(i,j-Mesh%Jp1,k,l)) &
                  + Mesh%dxdydV(i,j,k)*(this%zfluxdxdy(i,j,k,l) - this%zfluxdxdy(i,j,k-Mesh%Kp1,l)) &
                  - this%geo_src(i,j,k,l) - this%src(i,j,k,l)
            END DO
          END DO
        END DO
      END DO

   CASE(1)
!    grav => GetSourcesPointer(Physics%Sources, GRAVITY)
!    NULLIFY(pot)
!    IF(ASSOCIATED(grav)) THEN
!      IF(.NOT.grav%addtoenergy) THEN
!        pot => RemapBounds(Mesh,grav%pot(:,:,:))
!      END IF
!    END IF
!    phi = 0.
!    DO j=Mesh%JMIN-1,Mesh%JMAX
!      DO i=Mesh%IMIN-1,Mesh%IMAX
!        wp = 0.5*(this%w(i)+this%w(i+1)) + Mesh%hy%faces(i+1,j,1)*Mesh%omega
!        IF(ASSOCIATED(pot)) &
!          phi = pot(i,j,2)
!        IF(Physics%PRESSURE.GT.0) &
!          this%xfluxdy(i,j,Physics%ENERGY) &
!            = this%xfluxdy(i,j,Physics%ENERGY) &
!              + wp * (0.5 * wp * this%xfluxdy(i,j,Physics%DENSITY) &
!                      + this%xfluxdy(i,j,Physics%YMOMENTUM)) &
!              + phi*this%xfluxdy(i,j,Physics%DENSITY)
!
!        this%xfluxdy(i,j,Physics%YMOMENTUM) &
!          = (this%xfluxdy(i,j,Physics%YMOMENTUM) &
!            + wp * this%xfluxdy(i,j,Physics%DENSITY)) * Mesh%hy%faces(i+1,j,1)
!
!        wp = this%w(i) + Mesh%radius%bcenter(i,j)*Mesh%omega
!        IF(ASSOCIATED(pot)) &
!          phi = pot(i,j,3)
!        IF(Physics%PRESSURE.GT.0) &
!          this%yfluxdx(i,j,Physics%ENERGY) &
!            = this%yfluxdx(i,j,Physics%ENERGY) &
!             + wp * ( 0.5 * wp * this%yfluxdx(i,j,Physics%DENSITY) &
!                       + this%yfluxdx(i,j,Physics%YMOMENTUM)) &
!             + phi*this%yfluxdx(i,j,Physics%DENSITY)
!!
!        this%yfluxdx(i,j,Physics%YMOMENTUM) &
!          = this%yfluxdx(i,j,Physics%YMOMENTUM) &
!            + wp * this%yfluxdx(i,j,Physics%DENSITY) !* Mesh%hy%bcenter(i,j)
!      END DO
!    END DO
!
!    DO j=Mesh%JMIN,Mesh%JMAX
!     DO i=Mesh%IMIN,Mesh%IMAX
!        wp = pvar(i,j,Physics%YVELOCITY) + this%w(i) + Mesh%radius%bcenter(i,j)*Mesh%omega
!        this%geo_src(i,j,Physics%XMOMENTUM) &
!          = pvar(i,j,Physics%DENSITY) * wp**2 * Mesh%cyxy%bcenter(i,j) &
!            + pvar(i,j,Physics%DENSITY) * pvar(i,j,Physics%XVELOCITY) * wp * Mesh%cxyx%bcenter(i,j)
!
!        this%geo_src(i,j,Physics%YMOMENTUM) &
!          = cvar(i,j,Physics%XMOMENTUM) &
!            * ( pvar(i,j,Physics%XVELOCITY) * Mesh%cxyx%center(i,j) )
!
!        IF(Physics%PRESSURE.GT.0) THEN
!          this%geo_src(i,j,Physics%XMOMENTUM) &
!            = this%geo_src(i,j,Physics%XMOMENTUM) &
!             +pvar(i,j,Physics%PRESSURE) &
!                *( Mesh%cyxy%center(i,j) + Mesh%czxz%center(i,j) )
!          this%geo_src(i,j,Physics%YMOMENTUM) &
!            = this%geo_src(i,j,Physics%YMOMENTUM) &
!              + pvar(i,j,Physics%PRESSURE) &
!                *( Mesh%cxyx%center(i,j) + Mesh%czyz%center(i,j) )
!!          this%geo_src(i,j,Physics%XMOMENTUM) &
!!            = this%geo_src(i,j,Physics%XMOMENTUM) &
!!              - 0.5*(pvar(i+1,j,Physics%PRESSURE) - pvar(i-1,j,Physics%PRESSURE))/Mesh%dlx(i,j)
!!          this%geo_src(i,j,Physics%YMOMENTUM) &
!!            = this%geo_src(i,j,Physics%YMOMENTUM) &
!!              - 0.5*(pvar(i,j+1,Physics%PRESSURE) - pvar(i,j-1,Physics%PRESSURE))/Mesh%dly(i,j)
!          this%geo_src(i,j,Physics%ENERGY) = 0.
!        ELSE
!          this%geo_src(i,j,Physics%XMOMENTUM) &
!           = this%geo_src(i,j,Physics%XMOMENTUM) &
!             +pvar(i,j,Physics%DENSITY)*Physics%bccsound(i,j)**2 &
!                * ( Mesh%cyxy%center(i,j) + Mesh%czxz%center(i,j) )
!          this%geo_src(i,j,Physics%YMOMENTUM) &
!            = this%geo_src(i,j,Physics%YMOMENTUM) &
!              + pvar(i,j,Physics%DENSITY)*Physics%bccsound(i,j)**2 &
!                *( Mesh%cxyx%center(i,j) + Mesh%czyz%center(i,j) )
!        END IF
!      END DO
!    END DO
!
!!NEC$ NOVECTOR
!    DO k=1,Physics%VNUM
!!NEC$ OUTERLOOP_UNROLL(8)
!      DO j=Mesh%JMIN,Mesh%JMAX
!!NEC$ IVDEP
!        DO i=Mesh%IMIN,Mesh%IMAX
!          ! update right hand side of ODE
!          IF(k.EQ.Physics%YMOMENTUM) THEN
!            rhs(i,j,k) = Mesh%dydV(i,j)*(this%xfluxdy(i,j,k) - this%xfluxdy(i-1,j,k)) &
!                         / Mesh%hy%center(i,j)
!          ELSE
!            rhs(i,j,k) = Mesh%dydV(i,j)*(this%xfluxdy(i,j,k) - this%xfluxdy(i-1,j,k))
!          END IF
!          ! time step update
!          rhs(i,j,k) = rhs(i,j,k) &
!                       + Mesh%dxdV(i,j)*(this%yfluxdx(i,j,k) - this%yfluxdx(i,j-1,k))
!        END DO
!      END DO
!    END DO
!
!    DO j=Mesh%JMIN,Mesh%JMAX
!      DO i=Mesh%IMIN,Mesh%IMAX
!        IF(ASSOCIATED(pot)) &
!          phi = pot(i,j,1)
!        wp = this%w(i) + Mesh%radius%bcenter(i,j) * Mesh%omega
!        rhs(i,j,Physics%YMOMENTUM) = rhs(i,j,Physics%YMOMENTUM) &
!          - wp * rhs(i,j,Physics%DENSITY)
!        IF(Physics%PRESSURE.GT.0) &
!          rhs(i,j,Physics%ENERGY) = rhs(i,j,Physics%ENERGY) &
!            - wp*( rhs(i,j,Physics%YMOMENTUM) &
!                   + 0.5 * wp* rhs(i,j,Physics%DENSITY)) &
!            - phi * rhs(i,j,Physics%DENSITY)
!      END DO
!    END DO
!
!!NEC$ NOVECTOR
!    DO k=1,Physics%VNUM
!!NEC$ OUTERLOOP_UNROLL(8)
!      DO j=Mesh%JMIN,Mesh%JMAX
!!NEC$ IVDEP
!        DO i=Mesh%IMIN,Mesh%IMAX
!          ! update right hand side of ODE
!          rhs(i,j,k) = rhs(i,j,k) - this%geo_src(i,j,k) - this%src(i,j,k)
!        END DO
!      END DO
!    END DO
    END SELECT
  END SUBROUTINE ComputeRHS

  !> \public compute the RHS of the spatially discretized PDE
  !!
  SUBROUTINE CheckData(this,Mesh,Physics,Fluxes,pvar,cvar,checkdatabm)
    USE physics_euler2dit_mod, ONLY : physics_euler2dit
    USE physics_euler2d_mod, ONLY : physics_euler2d
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(fluxes_base),   INTENT(IN)    :: Fluxes
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                   Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                          INTENT(INOUT) :: pvar,cvar
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
      CLASS IS(physics_euler2dit)
        val = MINVAL(phys%bccsound)
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
      SELECT TYPE(phys => Physics)
      CLASS IS(physics_euler2d)
        val = MINVAL(pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX, &
                          Mesh%KMIN:Mesh%KMAX,phys%PRESSURE))
        IF(val.LT.this%pmin) THEN
          ! warn now and stop after file output
          CALL this%Warning("CheckData","Pressure below allowed pmin value.")
          this%break = .TRUE.
        END IF
      CLASS DEFAULT
        CALL this%Warning("CheckData","check pressure selected for isothermal physics")
      END SELECT
    END IF

    IF(IAND(checkdatabm,CHECK_RHOMIN).NE.CHECK_NOTHING) THEN
      ! check for physics with density defined
      SELECT TYPE(phys => Physics)
      CLASS IS(physics_euler2dit)
        val = MINVAL(pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX, &
                          Mesh%KMIN:Mesh%KMAX,phys%DENSITY))
        IF(val.LT.this%rhomin) THEN
            ! warn now and stop after file output
            CALL this%Warning("CheckData","Density below allowed rhomin value.")
            this%break = .TRUE.
        END IF
      CLASS DEFAULT
        CALL this%Warning("CheckData","check density selected, but density not defined physics")
      END SELECT
    END IF

    IF(IAND(checkdatabm,CHECK_INVALID).NE.CHECK_NOTHING) THEN
      IF(ANY((cvar.NE.cvar).OR.(pvar.NE.pvar))) THEN
        ! warn now and stop after file output
        CALL this%Warning("CheckData","Found NaN in pvar or cvar.")
        this%break = .TRUE.
      END IF
      IF(ANY((cvar.GT.HUGE(cvar)).OR.pvar.GT.HUGE(pvar))) THEN
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


  !> \public Calculates the linear transport step in Fargo Scheme \cite mignone2012 .
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
  !! The FargoAdvection step is done in strides along the x- or y-direction
  !! depending on the setup
  SUBROUTINE FargoAdvection(this,Fluxes,Mesh,Physics,Sources)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(sources_base),  POINTER       :: Sources
    !------------------------------------------------------------------------!
    REAL                 :: delxy
    INTEGER              :: i,j,k,l,sm,order
    REAL, DIMENSION(Mesh%JMIN:Mesh%JMAX) :: buf
    CHARACTER(LEN=80)    :: str
    REAL                 :: wi,wold
    INTEGER              :: GMIN1,GMIN2,GMAX1,GMAX2,MIN1,MIN2,MAX1,MAX2
#ifdef PARALLEL
    INTEGER              :: status(MPI_STATUS_SIZE)
    INTEGER              :: ierror
    INTEGER              :: PAR_DIMS,LINECOMM1,DIRECTION1,DIRECTION2,MINNUM1
    REAL                 :: mpi_buf(2*Mesh%GNUM)
#endif
    !------------------------------------------------------------------------!
    order = this%fargo_order
    IF(order.EQ.0) &
      CALL Fluxes%CalculateFaceData(Mesh,Physics,this%pvar,this%cvar)

    ! depending on the shear direction set universal values which are then
    ! used in the following
    IF (Mesh%SN_shear) THEN
      GMIN1     = Mesh%JGMIN
      GMAX1     = Mesh%JGMAX
      GMIN2     = Mesh%IGMIN
      GMAX2     = Mesh%IGMAX
      MIN1      = Mesh%JMIN
      MAX1      = Mesh%JMAX
      MIN2      = Mesh%IMIN
      MAX2      = Mesh%IMAX
#ifdef PARALLEL
      PAR_DIMS  = Mesh%dims(1)
      LINECOMM1 = Mesh%Icomm
      MINNUM1   = Mesh%MININUM
      DIRECTION1 = EAST
      DIRECTION2 = WEST
#endif
    ELSE
      GMIN1     = Mesh%IGMIN
      GMAX1     = Mesh%IGMAX
      GMIN2     = Mesh%JGMIN
      GMAX2     = Mesh%JGMAX
      MIN1      = Mesh%IMIN
      MAX1      = Mesh%IMAX
      MIN2      = Mesh%JMIN
      MAX2      = Mesh%JMAX
#ifdef PARALLEL
      PAR_DIMS   = Mesh%dims(2)
      LINECOMM1  = Mesh%Jcomm
      MINNUM1    = Mesh%MINJNUM
      DIRECTION1 = NORTH
      DIRECTION2 = SOUTH
#endif
    END IF

    ! determine step size of integer shift and length of remaining transport step
    ! first compute the whole step
    IF (Mesh%SN_shear) THEN
      this%delxy(:,:)  = this%w(:,:) * this%dt / Mesh%dlx(Mesh%IMIN,:,:)
    ELSE
      this%delxy(:,:)  = this%w(:,:) * this%dt / Mesh%dly(:,Mesh%JMIN,:)
    END IF

#ifdef PARALLEL
    ! make sure all MPI processes use the same step if domain is decomposed
    ! along the x- or y-direction (can be different due to round-off errors)
    IF (PAR_DIMS.GT.1) THEN
      CALL MPI_Allreduce(MPI_IN_PLACE,this%delxy,(GMAX1-GMIN1+1)*(Mesh%KGMAX-Mesh%KGMIN+1), &
                         DEFAULT_MPI_REAL,MPI_MIN,LINECOMM1,ierror)
    END IF
#endif

    ! then subdivide into integer shift and remaning linear advection step
    this%shift(:,:) = NINT(this%delxy(:,:))
    this%delxy(:,:)  = this%delxy(:,:)-DBLE(this%shift(:,:))

    ! advect with residual velocity
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO i=GMIN1,GMAX1
#ifdef PARALLEL
        IF(ABS(this%shift(i,k)).GT.MINNUM1) THEN
          WRITE(str,"(A,I4,A,I4,A)") "Shift of ",this%shift(i,k),&
            " at stride=", i, " too long."
          CALL this%Error("FargoAdvection",TRIM(str))
        END IF
#endif
        IF(abs(this%delxy(i,k)).GT.1.) THEN
          CALL this%Error("FargoAdvection","delxy bigger one should not be " &
          // "possible.")
        END IF

        ! residual shift along x- or y-direction
        IF(Mesh%SN_shear) THEN
          DO l=1,Physics%VNUM+Physics%PNUM
            CALL ResidualLineShift(Mesh,Fluxes,order,GMIN2,GMAX2,MIN2,MAX2, &
                                   this%delxy(i,k),this%cvar(:,i,k,l))
          END DO
        ELSE
          DO l=1,Physics%VNUM+Physics%PNUM
            CALL ResidualLineShift(Mesh,Fluxes,order,GMIN2,GMAX2,MIN2,MAX2, &
                                   this%delxy(i,k),this%cvar(i,:,k,l))
          END DO
        END IF
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
!                this%buf(k,1:this%shift(i)) = this%cvar(MAX2-this%shift(i)+1:MAX2,i,k)
!              ELSE
!                this%buf(k,1:this%shift(i)) = this%cvar(i,MAX2-this%shift(i)+1:MAX2,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(DIRECTION1), i+Mesh%GNUM, &
!              Mesh%neighbor(DIRECTION2), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%cvar(MAX2-this%shift(i)+1:MAX2,i,k) = this%buf(k,1:this%shift(i))
!              ELSE
!                this%cvar(i,MAX2-this%shift(i)+1:MAX2,k) = this%buf(k,1:this%shift(i))
!              END IF
!            END DO
!          ELSE IF(this%shift(i).LT.0) THEN
!            DO k=1,Physics%VNUM
!              IF(Mesh%SN_shear) THEN
!                this%buf(k,1:-this%shift(i)) = this%cvar(MIN2:MIN2-this%shift(i)-1,i,k)
!              ELSE
!                this%buf(k,1:-this%shift(i)) = this%cvar(i,MIN2:MIN2-this%shift(i)-1,k)
!              END IF
!            END DO
!            CALL MPI_Sendrecv_replace(&
!              this%buf,&
!              -this%shift(i)*Physics%VNUM, &
!              DEFAULT_MPI_REAL, &
!              Mesh%neighbor(DIRECTION2), i+Mesh%GNUM, &
!              Mesh%neighbor(DIRECTION1), i+Mesh%GNUM, &
!              Mesh%comm_cart, status, ierror)
!            DO k=1,Physics%VNUM
!              IF (Mesh%SN_shear) THEN
!                this%cvar(MIN2:MIN2-this%shift(i)-1,i,k) = this%buf(k,1:-this%shift(i))
!              ELSE
!                this%cvar(i,MIN2:MIN2-this%shift(i)-1,k) = this%buf(k,1:-this%shift(i))
!              END IF
!            END DO
!          END IF
!        END DO
!      END DO
!    END IF
!#endif

    ! Integer shift along x- or y-direction
    IF (Mesh%SN_shear) THEN
      DO l=1,Physics%VNUM+Physics%PNUM
        this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,l) &
          = CSHIFT(this%cvar(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,l), &
          -this%shift(Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),1)
      END DO
    ELSE
      DO l=1,Physics%VNUM+Physics%PNUM
        this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,l) &
          = CSHIFT(this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KGMIN:Mesh%KGMAX,l), &
          -this%shift(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX),2)
      END DO
    END IF

    ! Calculate RHS after the Advection Step
    CALL Physics%Convert2Primitive(Mesh,this%cvar,this%pvar)
    CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,this%time,0.,this%pvar,this%cvar,&
                    this%checkdatabm,this%rhs)

    this%cold(:,:,:,:) = this%cvar(:,:,:,:)

  END SUBROUTINE FargoAdvection


  SUBROUTINE ResidualLineShift(Mesh,Fluxes,order,GMIN2,GMAX2,MIN2,MAX2,shift,shifted_array)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(mesh_base),   INTENT(IN) :: Mesh
    CLASS(fluxes_base), INTENT(IN) :: Fluxes
    INTEGER,            INTENT(IN) :: order
    REAL,               INTENT(IN) :: shift
    INTEGER,            INTENT(IN) :: GMIN2,GMAX2,MIN2,MAX2
    REAL, DIMENSION(GMIN2:GMAX2), INTENT(INOUT) :: shifted_array
    !------------------------------------------------------------------------!
    REAL,DIMENSION(GMIN2:GMAX2) :: flux,dq,dql,qp,qm,q,dflux,qc
    REAL                 :: dqc, d2q, dqp, dqm, a, arg1, arg2
    INTEGER              :: i,j, ierror
#ifdef PARALLEL
    REAL                 :: mpi_buf(2*Mesh%GNUM)
    INTEGER              :: DIRECTION1, DIRECTION2, PAR_DIMS
    INTEGER              :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
    q(GMIN2:GMAX2) = shifted_array(GMIN2:GMAX2)
    dq(GMIN2) = q(GMIN2+1) - q(GMIN2)
    DO j=GMIN2+1,GMAX2-1
      dq(j) = q(j+1)-q(j)
      arg1 = dq(j-1)
      arg2 = dq(j)
      ! apply minmod limiter
      IF (SIGN(1.0,arg1)*SIGN(1.0,arg2).GT.0) THEN
        dql(j) = SIGN(MIN(ABS(arg1),ABS(arg2)),arg1)
      ELSE
        dql(j) = 0.
      END IF
      END DO
#ifdef PARALLEL
      IF (Mesh%SN_shear) THEN
        DIRECTION1 = EAST
        DIRECTION2 = WEST
        PAR_DIMS   = Mesh%dims(1)
      ELSE
        DIRECTION1 = NORTH
        DIRECTION2 = SOUTH
        PAR_DIMS   = Mesh%dims(2)
      END IF
      IF(PAR_DIMS.GT.1) THEN
        mpi_buf(1:Mesh%GNUM) = q(MAX2-Mesh%GNUM+1:MAX2)
        mpi_buf(Mesh%GNUM+1:2*Mesh%GNUM) = dq(MAX2-Mesh%GNUM+1:MAX2)
        CALL MPI_Sendrecv_replace(&
          mpi_buf,&
          4, &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(DIRECTION1), 51+DIRECTION1, &
          Mesh%neighbor(DIRECTION2), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        q(GMIN2:MIN2-1) = mpi_buf(1:Mesh%GNUM)
        dq(GMIN2:MIN2-1) = mpi_buf(Mesh%GNUM+1:2*Mesh%GNUM)

        mpi_buf(1:Mesh%GNUM) = q(MIN2:MIN2+Mesh%GNUM-1)
        mpi_buf(Mesh%GNUM+1:2*Mesh%GNUM) = dq(MIN2:MIN2+Mesh%GNUM-1)
        CALL MPI_Sendrecv_replace(&
          mpi_buf,&
          4, &
          DEFAULT_MPI_REAL, &
          Mesh%neighbor(DIRECTION2), 51+DIRECTION2, &
          Mesh%neighbor(DIRECTION1), MPI_ANY_TAG, &
          Mesh%comm_cart, status, ierror)
        q(MAX2+1:GMAX2) = mpi_buf(1:Mesh%GNUM)
        dq(MAX2+1:GMAX2) = mpi_buf(Mesh%GNUM+1:2*Mesh%GNUM)
      END IF
#endif
      SELECT CASE(order)
      ! first order
!      CASE(0)
!        IF(shift.GT.0.) THEN
!          DO j=MIN2-1,MAX2
!            flux(j) = Fluxes%prim(i,j,4,k)*(1. - shift) + shift * q(j) !order difference
!          END DO
!        ELSE
!          DO j=MIN2-1,MAX2
!            flux(j) = Fluxes%prim(i,j+1,3,k)*(1.+shift) - shift * q(j+1) !order difference
!          END DO
!        END IF
      ! second order
      CASE(2)
        IF(shift.GT.0.) THEN
          DO j=MIN2-1,MAX2
            flux(j) = q(j) + .5 * dql(j) * (1. - shift)
          END DO
        ELSE
          DO j=MIN2-1,MAX2
            flux(j) = q(j+1) - .5*dql(j+1)*(1. + shift)
          END DO
        END IF
      ! third order
      CASE(3)
        DO j=MIN2-1,MAX2+1
          dqp =  0.5*dq(j)   - (dql(j+1) - dql(j))/6.
          dqm = -0.5*dq(j-1) + (dql(j-1) - dql(j))/6.
          IF(dqp*dqm.GT.0.) THEN
            dqp = 0.
            dqm = 0.
          ELSE
            IF(ABS(dqp).GE.2.*ABS(dqm)) &
              dqp = -2. * dqm
            IF(ABS(dqm).GE.2.*ABS(dqp)) &
              dqm = -2.*dqp
          END IF
          qp(j) = q(j) + dqp
          qm(j) = q(j) + dqm
        END DO
        IF(shift.GT.0.) THEN
          DO j=MIN2-1,MAX2
            dqc = qp(j) - qm(j)
            d2q = qp(j)- 2.*q(j) + qm(j)
            flux(j) = qp(j) - .5 * shift * (dqc + d2q*(3. - 2.*shift))
          END DO
        ELSE
          DO j=MIN2-1,MAX2
            dqc =qp(j+1) - qm(j+1)
            d2q = qp(j+1) - 2.*q(j+1) + qm(j+1)
            flux(j) = qm(j+1) - 0.5*shift*(dqc - d2q*(3. + 2.*shift))
          END DO
        END IF
      END SELECT

      DO j=MIN2,MAX2
        shifted_array(j) = q(j) - shift*(flux(j) - flux(j-1))
      END DO

  END SUBROUTINE ResidualLineShift


  !> \public Calculates new background velocity for fargo advection
  !!
  !! \attention Only works when velocity is shifted in second direction.
  SUBROUTINE CalcBackgroundVelocity(this,Mesh,Physics,pvar,cvar,w)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM+Physics%PNUM), &
                          INTENT(INOUT) :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                          INTENT(OUT)   :: w
    !------------------------------------------------------------------------!
    REAL              :: wi
    INTEGER           :: i,k
#ifdef PARALLEL
    INTEGER           :: ierror
#endif
    !------------------------------------------------------------------------!
    ! make sure we are using true, i.e. non-selenoidal, quantities
    CALL Physics%AddBackgroundVelocityY(Mesh,w,pvar,cvar)
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO i=Mesh%IGMIN,Mesh%IGMAX
        ! some up all yvelocities along the y-direction
        wi = SUM(this%pvar(i,Mesh%JMIN:Mesh%JMAX,k,Physics%YVELOCITY))
#ifdef PARALLEL
        ! extend the sum over all partitions
        IF(Mesh%dims(2).GT.1) THEN
           CALL MPI_AllReduce(MPI_IN_PLACE, wi, 1, DEFAULT_MPI_REAL, MPI_SUM, &
                              Mesh%Jcomm, ierror)
        END IF
#endif
        ! set new background velocity to the arithmetic mean of the
        ! yvelocity field along the y-direction
        w(i,k) = wi / Mesh%JNUM
      END DO
    END DO
  END SUBROUTINE CalcBackgroundVelocity


  FUNCTION GetCentrifugalVelocity(this,Mesh,Physics,Fluxes,Sources,&
                       dir_omega_,accel_,centrot) RESULT(velo)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(timedisc_base), INTENT(INOUT) :: this
    CLASS(mesh_base),     INTENT(IN)    :: Mesh
    CLASS(physics_base),  INTENT(INOUT) :: Physics
    CLASS(fluxes_base),   INTENT(INOUT) :: Fluxes
    CLASS(sources_base),  POINTER       :: Sources
    REAL, DIMENSION(3),   INTENT(IN)    :: dir_omega_
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), OPTIONAL, &
                          INTENT(IN)    :: accel_
    REAL, DIMENSION(3), OPTIONAL, &
                          INTENT(IN)    :: centrot
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DIM) &
                                        :: velo
    !------------------------------------------------------------------------!
    REAL               :: time
    REAL, DIMENSION(3) :: dir_omega
    REAL               :: omega2
    INTEGER            :: k,i,j
    REAL               :: rotoemga
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3) &
                       :: bccart, bcposvec,  accel
    REAL,DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX) &
                       :: tmp
    !------------------------------------------------------------------------!

    IF(PRESENT(accel_)) THEN
      accel = accel_
    ELSE
      ! Works only for centrot = 0
      IF(PRESENT(centrot)) &
        CALL this%Error("GetCentrifugalVelocity","You are not allowed to "&
          //"define centrot without accel.")
      ! This may not work for physics with unusual conservative variables.
      ! We assume
      !   conservative momentum = density * velocity
      ! This is not true for the second component of physics with
      ! angular momentum transport, but that component should be zero.
      SELECT CASE(Physics%GetType())
      CASE(EULER2D,EULER2D_ISOTHERM,EULER3D_ROTSYM,EULER3D_ROTAMT,&
           EULER3D_ROTSYMSGS,EULER2D_SGS,EULER3D_ROTAMTSGS,EULER2D_ISOIAMT, &
           EULER2D_IAMT,EULER2D_ISOIAMROT)
        ! do nothing
      CASE DEFAULT
        CALL this%Error("GetCentrifugalVelocity","It is unknown, if the "&
          //"selected physics module works with this routine.")
      END SELECT
      CALL this%ComputeRHS(Mesh,Physics,Sources,Fluxes,this%time,0.0,this%pvar,this%cvar,this%checkdatabm,this%rhs)
      ! HERE DEPENDEND ON Physics
      DO k=Physics%XMOMENTUM,Physics%XMOMENTUM+Physics%DIM-1
        accel(:,:,:,k-Physics%XMOMENTUM+1) = -1. * this%rhs(:,:,:,k) &
                                           / this%pvar(:,:,:,Physics%DENSITY)
      END DO
    END IF

    dir_omega = dir_omega_

    ! be sure: |dir_omega| == 1
    omega2 = SUM(dir_omega(:)*dir_omega(:))
    ! omega must not be the zero vector
    IF (omega2 .EQ. 0.0) &
        CALL this%Error("GetCentrifugalVelocity", &
           "omega must not be the zero vector")
    ! norm must be one
    IF (omega2 .NE. 1.0) dir_omega(:) = dir_omega(:) / SQRT(omega2)

    IF ((Physics%DIM .EQ. 2) .AND. &
        ((dir_omega(1) .NE. 0.0) .OR. (dir_omega(2) .NE. 0.0))) &
        CALL this%Error("GetCentrifugalVelocity", &
           "the direction of omega should be (0,0,+-1) in case of two dimensions")


    IF (present(centrot)) THEN
      ! translate the position vector to the center of rotation
      bccart(:,:,:,1) = Mesh%bccart(:,:,:,1) - centrot(1)
      bccart(:,:,:,2) = Mesh%bccart(:,:,:,2) - centrot(2)

      ! compute curvilinear components of translated position vectors
      CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,bccart,bcposvec)
    ELSE
      bcposvec = Mesh%posvec%bcenter
    END IF


    ! compute distance to axis of rotation (It is automatically fulfilled in 2D.)
    IF (Physics%DIM .GT. 2) THEN
      tmp(:,:,:) =  bcposvec(:,:,:,1)*dir_omega(1) &
                  + bcposvec(:,:,:,2)*dir_omega(2) &
                  + bcposvec(:,:,:,3)*dir_omega(3)
      bcposvec(:,:,:,1) = bcposvec(:,:,:,1) - tmp(:,:,:)*dir_omega(1)
      bcposvec(:,:,:,2) = bcposvec(:,:,:,2) - tmp(:,:,:)*dir_omega(2)
      bcposvec(:,:,:,3) = bcposvec(:,:,:,3) - tmp(:,:,:)*dir_omega(3)
    END IF

    ! compute omega = SQRT(-dot(g,r)/|r|**2)
    tmp(:,:,:) = SQRT(MAX(0.0,-SUM(accel(:,:,:,1:3)*bcposvec(:,:,:,:),DIM=4))&
                    / SUM(bcposvec(:,:,:,:)*bcposvec(:,:,:,:),DIM=4))

    ! v / |omega| = dir_omega x r
    velo(:,:,:,:) = CROSS_PRODUCT(Mesh,Physics,dir_omega,bcposvec)
    ! v = |omega| * dir_omega x r
    DO k=1,Physics%DIM
      velo(:,:,:,k) = tmp(:,:,:)*velo(:,:,:,k)
    END DO

  CONTAINS
    FUNCTION CROSS_PRODUCT(Mesh,Physics,a,b) RESULT(cp)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      CLASS(mesh_base),    INTENT(IN)    :: Mesh
      CLASS(physics_base), INTENT(IN) :: Physics
      REAL, DIMENSION(3),  INTENT(IN) :: a
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
                           INTENT(IN) :: b
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DIM) &
                                      :: cp
      !------------------------------------------------------------------------!
      IF (Physics%DIM .EQ. 3) THEN
        cp(:,:,:,1) = a(2)*b(:,:,:,3) - a(3)*b(:,:,:,2)
        cp(:,:,:,2) = a(3)*b(:,:,:,1) - a(1)*b(:,:,:,3)
        cp(:,:,:,3) = a(1)*b(:,:,:,2) - a(2)*b(:,:,:,1)
      ELSE
        ! => a(1) = a(2) = 0 and a(3) = 1 or -1
        cp(:,:,:,1) = - a(3)*b(:,:,:,2)
        cp(:,:,:,2) = a(3)*b(:,:,:,1)
      END IF
     END FUNCTION CROSS_PRODUCT
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
      this%pvar,this%cvar,this%cold,this%ptmp,this%ctmp, &
      this%geo_src,this%src, &
      this%xfluxdydz,this%yfluxdzdx,this%zfluxdxdy,this%amax,this%tol_abs,&
      this%dtmean,this%dtstddev,this%time)

    IF (ASSOCIATED(this%w)) DEALLOCATE(this%w)
    IF (ASSOCIATED(this%delxy))DEALLOCATE(this%delxy)
    IF (ASSOCIATED(this%shift))DEALLOCATE(this%shift)
#ifdef PARALLEL
    IF(ASSOCIATED(this%buf))  DEALLOCATE(this%buf)
#endif
    IF(ASSOCIATED(this%bflux)) DEALLOCATE(this%bflux)
    IF(this%write_error) DEALLOCATE(this%errorval)
    IF(ASSOCIATED(this%solution)) DEALLOCATE(this%solution)
  END SUBROUTINE Finalize_base


END MODULE timedisc_base_mod
