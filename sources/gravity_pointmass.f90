!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_pointmass.f03                                             #
!#                                                                           #
!# Copyright (C) 2007-2018                                                   #
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
!! - parameters of \link gravity_pointmass \endlink as key-values
!! \key{potential,INTEGER,type of the potential}
!! \key{mass,REAL,mass of the point mass, 1.0}
!! \key{x,REAL,cartesian x-position of the point mass,0.0}
!! \key{y,REAL,cartesian y-position of the point mass,0.0}
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
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_pointmass_mod
  USE gravity_base_mod
  USE boundary_base_mod
  USE fluxes_base_mod
  USE physics_base_mod
  USE mesh_base_mod
  USE logging_base_mod
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
    CHARACTER(LEN=32) :: gravity_name = "central point mass"
    CLASS(logging_base), ALLOCATABLE    :: potential    !< newton or wiita
    REAL, POINTER                       :: mass         !< mass of pointmass
    REAL, POINTER                       :: accrate      !< true (limited) accretion rate
    REAL, POINTER                       :: massloss     !< mass loss due to acc. limit
    REAL                                :: mdot         !< mass flux at the boundary
    REAL                                :: acclimit     !< accretion limit
    REAL, DIMENSION(:,:),       POINTER :: pos          !< 3D cart. positions
    REAL, DIMENSION(:),         POINTER :: r0           !< center of mass
    REAL, DIMENSION(:,:,:),     POINTER :: r_prim       !< distance to primary point mass
    REAL, DIMENSION(:,:,:,:),   POINTER :: fr_prim
    REAL, DIMENSION(:,:,:,:),   POINTER :: posvec_prim  !< pos. vectors from primary
    REAL, DIMENSION(:,:,:,:,:), POINTER :: fposvec_prim !< face pos.
  CONTAINS
    PROCEDURE :: InitGravity_pointmass
    PROCEDURE :: CalcPotential
    PROCEDURE :: UpdateGravity_single
    PROCEDURE :: InfoGravity
    PROCEDURE :: CalcDiskHeight_single
    PROCEDURE :: Finalize
  END TYPE

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       gravity_pointmass, &
       ! constants
       NEWTON, WIITA
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
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    TYPE(Dict_TYP),           POINTER       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER :: potential_number, valwrite, gtype
    INTEGER :: err,i,j,k
    REAL    :: r,r_schwarzschild=0.0,eps
    REAL    :: invdt_x, invdt_y
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", gtype)
    ! allocate memory for new gravity term
    CALL this%InitLogging(gtype,this%gravity_name)
    ! type of the potential
    CALL GetAttr(config, "potential", potential_number, NEWTON)
    ALLOCATE(logging_base::this%potential)
    CALL this%potential%InitLogging(potential_number,potential_name(potential_number))

    SELECT CASE(potential_number)
    CASE(NEWTON)
       r_schwarzschild = 0.0
    CASE(WIITA)
       ! Schwarzschild radius
       r_schwarzschild = 2.*Physics%constants%GN * this%mass / Physics%constants%C**2
    CASE DEFAULT
       CALL this%Error("InitGravity_pointmass", "potential must be either NEWTON or WIITA")
    END SELECT


    !\todo last dimensions are absolutely not clear if right.. What are the standing for?
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%DIM), &
             this%pot(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,4), &
             this%omega(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
             this%omega2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,1), &
             this%r_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),&
             this%fr_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
             this%posvec_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
             this%fposvec_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3,Mesh%NDIMS),& ! last two entries (EAST,NORTH)x(dim1,dim2)
             this%mass, this%accrate, this%massloss, this%pos(1,3), &
         STAT = err)
    IF (err.NE.0) CALL this%Error("InitGravity_pointmass", "Unable allocate memory!")

    ! set mass
    CALL GetAttr(config, "mass", this%mass, 1.0)
    CALL SetAttr(IO, "mass", Ref(this%mass))

    ! cartesian position of the point mass
    CALL GetAttr(config, "x", this%pos(1,1), 0.0)
    CALL GetAttr(config, "y", this%pos(1,2), 0.0)
    CALL GetAttr(config, "z", this%pos(1,3), 0.0)

    ! output cartesian positions of the point mass
    CALL GetAttr(config, "output/position", valwrite, 0)
    IF (valwrite .EQ. 1) &
       CALL SetAttr(IO, "position", this%pos)

    ! Softening (e.g. for planets inside the computational domain)
    CALL GetAttr(config, "softening", eps, 0.0)

    ! soft switch on
    CALL GetAttr(config, "switchon", this%switchon, -1.0)
    IF (this%switchon.GT.0.0) this%scaling = 1.0

    ! accretion limit (e.g. eddington limit) divided by central mass
    ! the units are therefore (e.g. SI): kg/s / central_mass = 1 / s
    ! -1.0 == off
    CALL GetAttr(config, "acclimit", this%acclimit, -1.0)

    ! reset mass flux and massloss and register for output
    this%mdot = 0.0
    this%massloss = 0.0
    this%accrate = 0.0
    CALL SetAttr(IO, "accrate", Ref(this%accrate))
    IF (this%acclimit.GT.0.0) &
       CALL SetAttr(IO, "massloss", Ref(this%massloss))

    ! define position vector from the central object to all cell bary centers
    ! and the corresponding distances
    ! TODO Here should be a different evaluation for either 2D or 3D simulation
    IF (ABS(this%pos(1,1)).LE.TINY(this%pos(1,1)).AND.ABS(this%pos(1,2)).LE.TINY(this%pos(1,2))) THEN
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
       this%posvec_prim(:,:,:,1) = this%pos(1,1)
       this%posvec_prim(:,:,:,2) = this%pos(1,2)
       this%posvec_prim(:,:,:,3) = this%pos(1,3)

       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,EAST,:),this%posvec_prim(:,:,:,:),this%fposvec_prim(:,:,:,1,:))
       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,NORTH,:),this%posvec_prim(:,:,:,:),this%fposvec_prim(:,:,:,2,:))
       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,TOP,:),this%posvec_prim(:,:,:,:),this%fposvec_prim(:,:,:,3,:))
       CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,this%posvec_prim,this%posvec_prim)

       ! subtract the result from the position vector:
       ! this gives you the curvilinear components of all vectors pointing
       ! from the point mass to the bary center of any cell on the mesh
       this%posvec_prim(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - this%posvec_prim(:,:,:,:)
       this%fposvec_prim(:,:,:,1,:) = Mesh%posvec%faces(:,:,:,EAST,:) - this%fposvec_prim(:,:,:,1,:)  ! EAST
       this%fposvec_prim(:,:,:,2,:) = Mesh%posvec%faces(:,:,:,NORTH,:) - this%fposvec_prim(:,:,:,2,:) ! NORTH
       this%fposvec_prim(:,:,:,2,:) = Mesh%posvec%faces(:,:,:,TOP,:) - this%fposvec_prim(:,:,:,3,:)   ! TOP

       this%r_prim(:,:,:) = SQRT(this%posvec_prim(:,:,:,1)**2+this%posvec_prim(:,:,:,2)**2+this%posvec_prim(:,:,:,3)**2)
       this%fr_prim(:,:,:,1) = SQRT(this%fposvec_prim(:,:,:,1,1)**2+this%fposvec_prim(:,:,:,1,2)**2+this%fposvec_prim(:,:,:,1,3)**2) ! shifted EAST-faces
       this%fr_prim(:,:,:,2) = SQRT(this%fposvec_prim(:,:,:,2,1)**2+this%fposvec_prim(:,:,:,2,2)**2+this%fposvec_prim(:,:,:,2,3)**2) ! shifted NORTH-faces
       this%fr_prim(:,:,:,3) = SQRT(this%fposvec_prim(:,:,:,3,1)**2+this%fposvec_prim(:,:,:,3,2)**2+this%fposvec_prim(:,:,:,3,3)**2) ! shifted TOP-faces
    END IF

    ! initialize gravitational acceleration and Keplerian angular velocity
    ! loop over ghost cells to ensure that all values are well defined
    ! and remain finite; this is of particular importance for the calculation
    ! of the disk scale height (see below)
!NEC$ COLLAPSE
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
        DO i=Mesh%IGMIN,Mesh%IGMAX
          ! compute local Keplerian angular velocity
          ! Since omega usually becomes infinite in the limit r -> r_prim (or a, if a > 0 )
          ! (unless softening is enabled, i. e. eps > 0) we avoid devision by zero by first
          ! limiting the inverse of omega to some very small value.
          this%omega2(i,j,k,1) = 1./ MAX(10*TINY(eps), &
             (( this%r_prim(i,j,k) * (this%r_prim(i,j,k)-r_schwarzschild)**2) + eps**3) / (Physics%Constants%GN*this%mass) )
          this%omega(i,j,k) = SQRT(this%omega2(i,j,k,1))
          ! curvilinear components of the gravitational acceleration
          this%accel(i,j,k,1:2) = -this%omega2(i,j,k,1) * this%posvec_prim(i,j,k,1:2)
        END DO
      END DO
    END DO

    !! \todo implement computation for pseudo-Newton Paczynski-Wiita potential
    CALL this%CalcPotential(Mesh,Physics,this%mass,this%r_prim,this%fr_prim,this%pot)

    ! enable mass accretion by setting "outbound" to one of the four boundaries
    ! of the computational domain (depends on mesh geometry)
    SELECT CASE(Mesh%geometry%GetType())
    CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR,SPHERICAL,OBLATE_SPHEROIDAL,SINHSPHERICAL)
       CALL GetAttr(config, "outbound", this%outbound, WEST)
    CASE(CYLINDRICAL,TANCYLINDRICAL,BIPOLAR)
       CALL GetAttr(config, "outbound", this%outbound, SOUTH)
    CASE DEFAULT
       CALL GetAttr(config, "outbound", this%outbound, 0)
       CALL this%WARNING("GravitySources_pointmass","geometry does not support accretion")
    END SELECT

    CALL this%InitGravity(Mesh,Physics,config,IO)
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

  SUBROUTINE InfoGravity(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(IN) :: this
    CLASS(mesh_base),         INTENT(IN) :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: param_str
    !------------------------------------------------------------------------!
    WRITE (param_str,'(ES8.2)') this%mass
    CALL this%Info("            potential:         " // &
         TRIM(this%potential%GetName()))
    CALL this%Info("            mass:              " // TRIM(param_str))
    IF (this%acclimit.GT.0.0) THEN
       WRITE (param_str,'(ES8.2)') this%acclimit
       CALL this%Info("            accretion limit:   " // TRIM(param_str))
    END IF
    IF (this%switchon.GT.0.0) THEN
       WRITE (param_str,'(ES8.2)') this%switchon
       CALL this%Info("            switchon time:     " // TRIM(param_str))
    END IF
  END SUBROUTINE InfoGravity


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
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                              INTENT(IN)    :: pvar
    REAL,                     INTENT(IN)    :: time,dt
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Physics%VNUM) :: bflux
    REAL                          :: oldmass,massfac,sqrmassfac,dmass,dmasslim
    !------------------------------------------------------------------------!
    ! update accel and omega only in case of accretion
    IF (this%outbound.NE.0) THEN
       ! get boundary flux
#ifdef PARALLEL
!NEC$ IEXPAND
       bflux(:)  = Fluxes%GetBoundaryFlux(Mesh,Physics,this%outbound,MPI_COMM_WORLD)
#else
!NEC$ IEXPAND
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
       this%accel(:,:,:,:) = this%accel(:,:,:,:) * massfac
       this%pot(:,:,:,:) = this%pot(:,:,:,:) * massfac
       this%omega(:,:,:) = this%omega(:,:,:) * sqrmassfac
       this%mdot = bflux(Physics%DENSITY)
    END IF

    ! modify acceleration during switchon phase
    IF (time.GE.0.0.AND.time.LE.this%switchon) THEN
      ! compute new scaling
      this%scaling = SIN(0.5*PI*time/this%switchon)**2
      ! compute acceleration and scale it;
      ! during the switchon phase accretion is disabled
      this%accel(:,:,:,1) = -this%scaling * this%omega2(:,:,:,1) * this%posvec_prim(:,:,:,1)
      this%accel(:,:,:,2) = -this%scaling * this%omega2(:,:,:,1) * this%posvec_prim(:,:,:,2)
    END IF

  END SUBROUTINE UpdateGravity_single


  !> \public computes pressure scale height of a geometrically thin Keplerian disk
  !!
  !! (see \link getdiskheight_pointmass \endlink )
  PURE SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
    CLASS(mesh_base),         INTENT(IN)    :: Mesh
    CLASS(physics_base),      INTENT(IN)    :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM), &
                              INTENT(IN)    :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(IN)    :: bccsound
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(OUT)   :: h_ext
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
                              INTENT(INOUT) :: height
    !------------------------------------------------------------------------!
    ! compute disk height
    h_ext(:,:,:) = GetDiskHeight(bccsound(:,:,:),this%omega(:,:,:))
  END SUBROUTINE CalcDiskHeight_single

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_pointmass), INTENT(INOUT) :: this
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

    DEALLOCATE(this%accel,this%pot,this%omega,this%omega2,this%r_prim,this%fr_prim, &
               this%posvec_prim,this%fposvec_prim,this%mass,this%accrate,this%massloss,this%pos)

    DEALLOCATE(this%potential)

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
