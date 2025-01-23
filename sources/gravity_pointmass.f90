!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_pointmass.f90                                             #
!#                                                                           #
!# Copyright (C) 2007-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \addtogroup gravity
!! - parameters of \link gravity_pointmass_mod gravity_pointmass \endlink as key-values
!! \key{potential,INTEGER,type of the potential}
!! \key{mass,REAL,mass of the point mass, 1.0}
!! \key{x,REAL,cartesian x-position of the point mass,0.0}
!! \key{y,REAL,cartesian y-position of the point mass,0.0}
!! \key{z,REAL,cartesian z-position of the point mass,0.0}
!! \key{softening,REAL,Softening (e.g. for planets inside the computational domain),0.0}
!! \key{switchon,REAL,soft switch on,-1.0}
!! \key{acclimit,REAL,max accretion rate / mass [1/unit time] (if negative => disabled),-1.0}
!! \key{outbound,INTEGER,enable mass accretion by setting "outbound" to one of the four boundaries,depends on mesh}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Jannes Klee
!!
!! \brief source terms module for gravitational acceleration due to a point
!! mass at the center of the coordinate system
!!
!! \extends gravity_base
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_pointmass_mod
  USE gravity_base_mod
  USE boundary_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE logging_base_mod
  USE marray_base_mod
  USE marray_compound_mod
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
  PRIVATE
  CHARACTER(LEN=16), PARAMETER :: potential_name(2) = (/ &
                                  "Newton          ", &
                                  "Paczynski-Wiita " /)
  INTEGER, PARAMETER :: NEWTON = 1
  INTEGER, PARAMETER :: WIITA  = 2

  TYPE, EXTENDS(gravity_base) :: gravity_pointmass
    CLASS(logging_base), ALLOCATABLE    :: potential    !< newton or wiita
    INTEGER                             :: outbound     !< outflow boundary
    REAL, POINTER                       :: mass         !< mass of pointmass
    REAL, POINTER                       :: accrate      !< true (limited) accretion rate
    REAL, POINTER                       :: massloss     !< mass loss due to acc. limit
    REAL                                :: mdot         !< mass flux at the boundary
    REAL                                :: acclimit     !< accretion limit
    REAL                                :: switchon     !< duration of soft switch-on phase
    REAL, DIMENSION(:,:),       POINTER :: pos          !< 3D cart. positions
    REAL, DIMENSION(:),         POINTER :: r0           !< center of mass
    REAL, DIMENSION(:,:,:),     POINTER :: r_prim       !< distance to primary point mass
    REAL, DIMENSION(:,:,:,:),   POINTER :: fr_prim
    REAL, DIMENSION(:,:,:,:),   POINTER :: posvec_prim  !< pos. vectors from primary
    REAL, DIMENSION(:,:,:,:),   POINTER :: posvec_prim_tmp !< tmp. pos. vectors from primary
    REAL, DIMENSION(:,:,:,:,:), POINTER :: fposvec_prim !< face pos.
    TYPE(marray_base), ALLOCATABLE      :: omega, &     !< angular velocity
                                           omega2       !< Omega Kepler squared
  CONTAINS
    PROCEDURE :: InitGravity_pointmass
    PROCEDURE :: SetOutput
    PROCEDURE :: CalcPotential
    PROCEDURE :: UpdateGravity_single
    PROCEDURE :: InfoGravity
    PROCEDURE :: CalcDiskHeight_single
    FINAL     :: Finalize
  END TYPE

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       gravity_pointmass, &
       ! constants
       NEWTON, WIITA, &
       ! elemental functions
       GetDiskHeight
  !--------------------------------------------------------------------------!

CONTAINS

  !> \public constructor for point mass module
  !!
  !! Initializes several data arrays including the local Keplerian angular
  !! velocity
  !! \f[
  !!     \Omega = \sqrt{\frac{GM}{r(r-a)^2}}
  !! \f]
  !! with \f$ a=0 \f$ for Newtonian gravity and \f$ a=R_s \f$ set to the
  !! Schwarzschild radius for pseudo-Newtonian Paczynski-Wiita potential.
  !! The gravitational acceleration is given by
  !! \f[
  !!     \vec{g} = -\frac{d\Phi}{dr} \hat{e}_r
  !!             = -\frac{GM}{(r-a)^2} \vec{r}/r
  !!             = -\Omega^2 \vec{r}
  !! \f]
  SUBROUTINE InitGravity_pointmass(this,Mesh,Physics,config,IO)
    USE geometry_logcylindrical_mod, ONLY: geometry_logcylindrical
    USE geometry_spherical_mod, ONLY: geometry_spherical
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    TYPE(Dict_TYP),           POINTER       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER :: potential_number, valwrite, gtype
    INTEGER :: err,i,j,k,l
    REAL    :: r_schwarzschild=0.0,eps
    !------------------------------------------------------------------------!
    CALL this%InitGravity(Mesh,Physics,"central point mass",config,IO)
    !\todo last dimensions are absolutely not clear if right.. What are the standing for?
    ALLOCATE(this%potential, this%pot, this%omega, this%omega2, &
             this%r_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),&
             this%fr_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
             this%posvec_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
             this%posvec_prim_tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
             this%fposvec_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3,3),&
               ! last two entries (EAST,NORTH,TOP)x(dim1,dim2,dim3)
             this%mass, this%accrate, this%massloss, this%pos(1,3), &
         STAT = err)
    IF (err.NE.0) CALL this%Error("InitGravity_pointmass", "Memory allocation failed!")

    ! type of the potential
    CALL GetAttr(config, "potential", potential_number, NEWTON)
    CALL this%potential%InitLogging(potential_number,potential_name(potential_number))

    ! init mesh array for potential and angular velocity
    this%pot    = marray_base(4)
    this%omega  = marray_base()
    this%omega2 = marray_base()
    this%pot%data1d(:)    = 0.0
    this%omega%data1d(:)  = 0.0
    this%omega2%data1d(:) = 0.0

    ! set (initial) mass
    CALL GetAttr(config, "mass", this%mass, 1.0)
    CALL SetAttr(IO, "mass", Ref(this%mass))

    ! cartesian position of the point mass
    CALL GetAttr(config, "x", this%pos(1,1), 0.0)
    CALL GetAttr(config, "y", this%pos(1,2), 0.0)
    CALL GetAttr(config, "z", this%pos(1,3), 0.0)

    ! Softening (e.g. for planets inside the computational domain)
    CALL GetAttr(config, "softening", eps, 0.0)

    ! soft switch on
    CALL GetAttr(config, "switchon", this%switchon, -1.0)

    ! accretion limit (e.g. eddington limit) divided by central mass
    ! the units are therefore (e.g. SI): kg/s / central_mass = 1 / s
    ! -1.0 == off
    CALL GetAttr(config, "acclimit", this%acclimit, -1.0)

    CALL SetAttr(IO, "accrate", Ref(this%accrate))
    IF (this%acclimit.GT.0.0) &
       CALL SetAttr(IO, "massloss", Ref(this%massloss))

    SELECT CASE(potential_number)
    CASE(NEWTON)
       r_schwarzschild = 0.0
    CASE(WIITA)
       ! Schwarzschild radius
       r_schwarzschild = 2.*Physics%constants%GN * this%mass / Physics%constants%C**2
    CASE DEFAULT
       CALL this%Error("InitGravity_pointmass", "potential must be either NEWTON or WIITA")
    END SELECT

    ! reset mass flux and massloss and register for output
    this%mdot = 0.0
    this%massloss = 0.0
    this%accrate = 0.0

    ! define position vector from the central object to all cell bary centers
    ! and the corresponding distances
    IF (ALL(ABS(this%pos(1,:)).LE.TINY(this%pos))) THEN
       ! no shift: point mass is located at the origin of the mesh
       ! set position vector and inverse radius using appropriate mesh arrays
       this%posvec_prim(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:)
       this%r_prim(:,:,:) = Mesh%radius%bcenter(:,:,:)
       this%fr_prim(:,:,:,1) = Mesh%radius%faces(:,:,:,EAST)
       this%fr_prim(:,:,:,2) = Mesh%radius%faces(:,:,:,NORTH)
       this%fr_prim(:,:,:,3) = Mesh%radius%faces(:,:,:,TOP)
    ELSE
       ! shifted point mass position:
       ! compute curvilinear components of shift vector
       ! PLEASE NOTE: We need the shifted curvilinear components at the bary_centers and
       ! faces in EASTERN and NORTHERN direction. This is the case because we want to
       ! calculate the potential at these face positions further below.
       this%posvec_prim_tmp(:,:,:,1) = this%pos(1,1)
       this%posvec_prim_tmp(:,:,:,2) = this%pos(1,2)
       this%posvec_prim_tmp(:,:,:,3) = this%pos(1,3)

       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,EAST,:), &
            this%posvec_prim_tmp(:,:,:,:),this%fposvec_prim(:,:,:,1,:))
       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,NORTH,:), &
            this%posvec_prim_tmp(:,:,:,:),this%fposvec_prim(:,:,:,2,:))
       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,TOP,:), &
            this%posvec_prim_tmp(:,:,:,:),this%fposvec_prim(:,:,:,3,:))

       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,this%posvec_prim_tmp,this%posvec_prim)

       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing
       ! from the point mass to the bary center of any cell on the mesh
       this%posvec_prim(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - this%posvec_prim(:,:,:,:)
       this%fposvec_prim(:,:,:,1,:) = Mesh%posvec%faces(:,:,:,EAST,:) - this%fposvec_prim(:,:,:,1,:)  ! EAST
       this%fposvec_prim(:,:,:,2,:) = Mesh%posvec%faces(:,:,:,NORTH,:) - this%fposvec_prim(:,:,:,2,:) ! NORTH
       this%fposvec_prim(:,:,:,3,:) = Mesh%posvec%faces(:,:,:,TOP,:) - this%fposvec_prim(:,:,:,3,:)   ! TOP

       this%r_prim(:,:,:) = SQRT(this%posvec_prim(:,:,:,1)**2 &
         + this%posvec_prim(:,:,:,2)**2 + this%posvec_prim(:,:,:,3)**2)
       this%fr_prim(:,:,:,1) = SQRT(this%fposvec_prim(:,:,:,1,1)**2 &
         + this%fposvec_prim(:,:,:,1,2)**2 + this%fposvec_prim(:,:,:,1,3)**2) ! shifted EAST-faces
       this%fr_prim(:,:,:,2) = SQRT(this%fposvec_prim(:,:,:,2,1)**2 &
         + this%fposvec_prim(:,:,:,2,2)**2 + this%fposvec_prim(:,:,:,2,3)**2) ! shifted NORTH-faces
       this%fr_prim(:,:,:,3) = SQRT(this%fposvec_prim(:,:,:,3,1)**2 &
         +this%fposvec_prim(:,:,:,3,2)**2 + this%fposvec_prim(:,:,:,3,3)**2) ! shifted TOP-faces
    END IF

    ! Initialize gravitational acceleration and Keplerian angular velocity.
    ! Loop over ghost cells as well to ensure that all values are well defined
    ! and remain finite; this is of particular importance for the calculation
    ! of the disk scale height (see below).
!NEC$ COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! compute local Keplerian angular velocity
          ! Since omega usually becomes infinite in the limit r -> r_prim (or a, if a > 0 )
          ! (unless softening is enabled, i. e. eps > 0) we avoid devision by zero by first
          ! limiting the inverse of omega to some very small value.
          this%omega2%data3d(i,j,k) = 1./ MAX(10*TINY(eps), &
             (( this%r_prim(i,j,k) * (this%r_prim(i,j,k)-r_schwarzschild)**2) + eps**3) / (Physics%Constants%GN*this%mass) )
          this%omega%data3d(i,j,k) = SQRT(this%omega2%data3d(i,j,k))
          ! curvilinear components of the gravitational acceleration;
          ! account for dimensionality of the acceleration vector (could be < 3)
!NEC$ SHORTLOOP
          DO l=1,Physics%VDIM
            this%accel%data4d(i,j,k,l) = -this%omega2%data3d(i,j,k) * this%posvec_prim(i,j,k,Physics%VIDX(l))
          END DO
        END DO
      END DO
    END DO

    !! \todo implement computation for pseudo-Newton Paczynski-Wiita potential
    CALL this%CalcPotential(Mesh,Physics,this%mass,this%r_prim,this%fr_prim,this%pot%data4d)

    CALL this%SetOutput(Mesh,Physics,config,IO)

    ! enable mass accretion by setting "outbound" to one of the four boundaries
    ! of the computational domain (depends on mesh geometry)
    SELECT TYPE(geo=>Mesh%geometry)
    CLASS IS(geometry_cylindrical)
       CALL GetAttr(config, "outbound", this%outbound, WEST)
    CLASS IS(geometry_spherical)
       CALL GetAttr(config, "outbound", this%outbound, WEST)
    CLASS DEFAULT
       CALL GetAttr(config, "outbound", this%outbound, 0)
       CALL this%WARNING("GravitySources_pointmass","geometry does not support accretion")
    END SELECT
    CALL this%InfoGravity()
  END SUBROUTINE InitGravity_pointmass

  SUBROUTINE CalcPotential(this,Mesh,Physics,mass,dist,dist_faces,pot)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass)         :: this
    CLASS(mesh_base),    INTENT(IN)  :: Mesh
    CLASS(physics_base), INTENT(IN)  :: Physics
    REAL                             :: mass
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                         INTENT(IN)  :: dist
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
                         INTENT(IN)  :: dist_faces
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4), &
                         INTENT(OUT) :: pot
    !------------------------------------------------------------------------!

    !! \todo implement computation for pseudo-Newton Paczynski-Wiita potential
    pot(:,:,:,1) = - Physics%constants%GN * mass / dist(:,:,:)
    pot(:,:,:,2) = - Physics%constants%GN * mass / dist_faces(:,:,:,1) ! direction EAST !Mesh%radius%faces(:,:,:,EAST)
    pot(:,:,:,3) = - Physics%constants%GN * mass / dist_faces(:,:,:,2) ! direction NORTH !Mesh%radius%faces(:,:,:,NORTH)
    pot(:,:,:,4) = - Physics%constants%GN * mass / dist_faces(:,:,:,3) ! direction TOP !Mesh%radius%faces(:,:,:,TOP)
  END SUBROUTINE

  SUBROUTINE InfoGravity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(IN) :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: param_str
    !------------------------------------------------------------------------!
    WRITE (param_str,'(ES10.2)') this%mass
    CALL this%Info("            potential:         " // &
         TRIM(this%potential%GetName()))
    CALL this%Info("            initial mass:      " // TRIM(param_str))
    WRITE (param_str, '(3(ES10.2))') this%pos(1,1:3)
    CALL this%Info("            cart. position:   " // TRIM(param_str))
    IF (this%outbound.GT.0) THEN
      param_str = "enabled"
    ELSE
      param_str = "disabled"
    END IF
    CALL this%Info("            accretion:         " // TRIM(param_str))
    IF (this%outbound.GT.0.AND.this%acclimit.GT.0.0) THEN
      WRITE (param_str,'(ES10.2)') this%acclimit
      CALL this%Info("            accretion limit:   " // TRIM(param_str))
    END IF
    IF (this%switchon.GT.0.0) THEN
       WRITE (param_str,'(ES10.2)') this%switchon
       CALL this%Info("            switchon time:     " // TRIM(param_str))
    END IF
  END SUBROUTINE InfoGravity

  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
    CLASS(mesh_base),    INTENT(IN)    :: Mesh
    CLASS(physics_base), INTENT(IN)    :: Physics
    TYPE(Dict_TYP),      POINTER       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER          :: valwrite
    !------------------------------------------------------------------------!
    ! output cartesian positions of the point mass
    valwrite = 0
    CALL GetAttr(config, "output/position", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "position", this%pos)
    valwrite = 0
    CALL GetAttr(config, "output/potential", valwrite, 0)
    IF ((valwrite .EQ. 1).AND. ALLOCATED(this%pot)) THEN
      IF (ASSOCIATED(this%pot%data4d)) &
        CALL SetAttr(IO, "potential", &
          this%pot%data4d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,4))
    END IF
  END SUBROUTINE SetOutput

  !> \public updates the gravitational acceleration of the pointmass module
  !!
  !! The acceleration may not change at all if accretion and/or the switchon
  !! procedure is disabled. In that case this subroutine does nothing, because
  !! appropriate values have already been set in \link initgravity_pointmass \endlink .
  SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    CLASS(fluxes_base),       INTENT(IN)    :: Fluxes
    CLASS(marray_compound),   INTENT(INOUT) :: pvar
    REAL,                     INTENT(IN)    :: time,dt
    !------------------------------------------------------------------------!
    INTEGER                       :: i,j,k,l
    REAL, DIMENSION(Physics%VNUM) :: bflux
    REAL                          :: scaling,oldmass,massfac,sqrmassfac,dmass,dmasslim
    !------------------------------------------------------------------------!
    ! update accel and omega only in case of accretion
    IF (this%outbound.NE.0) THEN
       ! get boundary flux
#ifdef PARALLEL
       bflux(:)  = Fluxes%GetBoundaryFlux(Mesh,Physics,this%outbound,MPI_COMM_WORLD)
#else
       bflux(:)  = Fluxes%GetBoundaryFlux(Mesh,Physics,this%outbound)
#endif
       ! store old mass
       oldmass = this%mass
       ! compute the mass flux difference
       ! REMARK #1: bflux is the accumulated flux across the boundary
       ! REMARK #2: mdot is NOT the accretion rate, but the
       !            mass flux from the previous evalution (see below)
       dmass = (bflux(Physics%DENSITY) - this%mdot)
       IF(this%acclimit.GE.0.) THEN
         dmasslim  = MIN(dmass,this%acclimit*dt*this%mass)
         this%massloss = this%massloss + dmass - dmasslim
         dmass = dmasslim
       END IF
       ! compute new mass and current true, i. e. possibly limited, accretion rate
       this%mass = this%mass + dmass
       IF (dt.GT.0.0) this%accrate = dmass/dt
       ! scale gravitational acceleration and Keplerian angular velocity
       ! with newmass/oldmass to account for accretion
       massfac = this%mass/oldmass
       sqrmassfac = SQRT(massfac)
       this%accel%data1d(:) = this%accel%data1d(:) * massfac
       this%pot%data1d(:) = this%pot%data1d(:) * massfac
       this%omega%data1d(:) = this%omega%data1d(:) * sqrmassfac
       this%mdot = bflux(Physics%DENSITY)
    END IF

    ! modify acceleration during switchon phase
    IF (time.GE.0.0.AND.time.LE.this%switchon) THEN
      ! compute new scaling
      scaling = SIN(0.5*PI*time/this%switchon)**2
      ! compute acceleration and scale it;
      ! during the switchon phase accretion is disabled
!NEC$ COLLAPSE
      DO k=Mesh%KGMIN,Mesh%KGMAX
        DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
          DO i=Mesh%IGMIN,Mesh%IGMAX
!NEC$ SHORTLOOP
            DO l=1,Physics%VDIM
              this%accel%data4d(i,j,k,l) = -scaling * this%omega2%data3d(i,j,k) * this%posvec_prim(i,j,k,Physics%VIDX(l))
            END DO
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE UpdateGravity_single


  !> \public computes pressure scale height of a geometrically thin Keplerian disk
  !!
  !! (see \ref getdiskheight )
  PURE SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    CLASS(marray_compound),   INTENT(INOUT) :: pvar
    TYPE(marray_base),        INTENT(INOUT) :: bccsound,h_ext,height
    !------------------------------------------------------------------------!
    ! compute disk height
    h_ext%data1d(:) = GetDiskHeight(bccsound%data1d(:),this%omega%data1d(:))
  END SUBROUTINE CalcDiskHeight_single

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(gravity_pointmass), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=128) :: buffer
    !------------------------------------------------------------------------!
    IF (this%GetRank().EQ.0) THEN
       CALL this%Info("-------------------------------------------------------------------")
       WRITE(buffer, "(A,(ES11.3))") " central mass: ", this%mass
       CALL this%Info(buffer)
       IF (this%acclimit.GT.0.0) THEN
          WRITE(buffer, "(A,(ES11.3))") "    mass loss: ", this%massloss
          CALL this%Info(buffer)
       END IF
    END IF

    DEALLOCATE(this%omega,this%omega2,this%r_prim,this%fr_prim, &
               this%posvec_prim,this%posvec_prim_tmp,this%fposvec_prim,this%mass, &
               this%accrate,this%massloss,this%pos)

    IF (ALLOCATED(this%potential)) DEALLOCATE(this%potential)
    IF (ALLOCATED(this%pot)) DEALLOCATE(this%pot)
    CALL this%Finalize_base()
  END SUBROUTINE Finalize

  !> \public return pressure scale height of a geometrically thin Keplerian disk
  !!
  !! The pressure scale height for an geometrically thin Keplerian disk
  !! is given by
  !! \f[
  !!     h = c_s / \Omega
  !! \f]
  !! where \f$ c_s \f$ is the local speed of sound and \f$ \Omega \f$ is
  !! the local Keplerian angular velocity (see \cite pringle1981 ).
  ELEMENTAL FUNCTION GetDiskHeight(cs,omega) RESULT(height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN) :: cs,omega
    REAL             :: height
    !------------------------------------------------------------------------!
    height = cs / omega
  END FUNCTION GetDiskHeight

END MODULE gravity_pointmass_mod
