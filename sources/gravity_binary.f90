!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: gravity_binary.f90                                                #
!#                                                                           #
!# Copyright (C) 2010-2024                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Anna Feiler      <afeiler@astrophysik.uni-kiel.de>                        #
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
!! - parameters of \link gravity_binary_mod gravity_binary \endlink as key-values
!! \key{mass1,REAL,mass of primary component,1.0}
!! \key{mass2,REAL,mass of secondary component,1.0}
!! \key{excentricity,REAL,excentricity,0.0}
!! \key{semimayoraxis,REAL,semi major axis,1.0}
!! \key{softening1,REAL,softening parameter of primary component,0.0}
!! \key{softening2,REAL,softening parameter of secondary component,0.0}
!! \key{switchon1,REAL,soft switch on,-1.0}
!! \key{omega_rot,REAL,convert binary positions to a rotating reference frame
!! \key{mesh_rot,INTEGER,convert binary positions to a rotating mesh (0|1),0}
!! with angular velocity omega ,0.0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Anna Feiler
!! \author Manuel Jung
!! \author Jannes Klee
!!
!! \brief source terms module for gravitational acceleration due to two
!! pointmasses
!!
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_binary_mod
  USE gravity_pointmass_mod
  USE gravity_base_mod
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
  TYPE, EXTENDS(gravity_pointmass) :: gravity_binary
    CLASS(marray_base), ALLOCATABLE     :: pot_prim,pot_sec !< potential of primary & secondary
    REAL, DIMENSION(:,:,:), POINTER     :: r_sec        !<    and to secondary point mass
    REAL, DIMENSION(:,:,:,:), POINTER   :: posvec_sec   !<   secondary to all cell bary centers
    REAL, DIMENSION(:,:,:,:), POINTER   :: posvec_sec_tmp !<  tmp. secondary to all cell bary centers
    REAL, DIMENSION(:,:,:,:), POINTER   :: fr_sec
    REAL, DIMENSION(:,:,:,:,:), POINTER :: fposvec_sec
    REAL, POINTER                       :: mass2        !< 2nd mass for binaries
    REAL                                :: excent       !< excentricity
    REAL                                :: semaaxis     !< semi major axis
    REAL                                :: period       !< period of binaries
    REAL                                :: eps1,eps2    !< softening parameter
    REAL                                :: switchon2    !< same for secondary in binary systems
    REAL                                :: omega_rot    !< angular velocity of rotating reference frame
  CONTAINS
    PROCEDURE :: InitGravity_binary
    PROCEDURE :: UpdateGravity_single
    PROCEDURE :: UpdatePositions
    PROCEDURE :: InfoGravity
    PROCEDURE :: SetOutput
    PROCEDURE :: GetMass_primary
    PROCEDURE :: GetMass_secondary
    PROCEDURE :: CalcDiskHeight_single
    FINAL :: Finalize
  END TYPE

  PUBLIC :: &
       ! types
       gravity_binary
  !--------------------------------------------------------------------------!

CONTAINS
  !> \public Constructor of gravity binary module
  !!
  !! This subroutine reads the necessary config data for the binary.
  !! It initializes the gravity type and various mesh data arrays. Some of those
  !! are marked for output.
  SUBROUTINE InitGravity_binary(this,Mesh,Physics,config,IO)
    USE physics_euler_mod, ONLY: physics_euler
    USE physics_eulerisotherm_mod, ONLY: physics_eulerisotherm
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(INOUT) :: this  !< \param [in,out] this all gravity data
    CLASS(mesh_base),      INTENT(IN) :: Mesh     !< \param [in] Mesh mesh type
    CLASS(physics_base),   INTENT(IN) :: Physics  !< \param [in] Physics physics type
    TYPE(Dict_TYP),POINTER            :: config   !< \param [in,out] config sub-dictionary
                                                  !! with binary configuration data
    TYPE(Dict_TYP),POINTER            :: IO       !< \param [in,out] IO output sub-dictionary
    !------------------------------------------------------------------------!
    INTEGER           :: mesh_rot
    INTEGER           :: err
    !------------------------------------------------------------------------!
    ! init base class
    CALL this%InitGravity(Mesh,Physics,"binary point masses",config,IO)

    SELECT TYPE(Physics)
    TYPE IS(physics_euler)
       ! do nothing
    TYPE IS(physics_eulerisotherm)
      ! do nothing
    CLASS DEFAULT
       CALL this%Error("InitGravity_binary", &
            "Physics not supported for binary gravity term")
    END SELECT

    ALLOCATE(&
         this%mass,this%massloss,this%accrate,this%mass2,this%pos(3,2),this%r0(3), &
         this%omega2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,2), &
         this%omega(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%r_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),&
         this%r_sec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),&
         this%posvec_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
         this%posvec_prim_tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
         this%posvec_sec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
         this%posvec_sec_tmp(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3),&
         this%pot,this%pot_prim,this%pot_sec, &
         this%fr_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
         this%fr_sec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3), &
         this%fposvec_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3,3),&
         this%fposvec_sec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,3,3),&
           ! last two entries (EAST,NORTH)x(dim1,dim2)
         STAT = err)
    IF (err.NE.0) CALL this%Error("gravity_pointmass::InitGravity_pointmass", "Unable allocate memory!")

    ! initialize potentials
    this%pot = marray_base(4)
    this%pot_prim = marray_base(4)
    this%pot_sec  = marray_base(4)
    this%pot%data1d(:) = 0.0
    this%pot_prim%data1d(:) = 0.0
    this%pot_sec%data1d(:) = 0.0

    ! mass of primary component
    CALL GetAttr(config, "mass1", this%mass, 1.0)
    CALL SetAttr(IO, "mass1", Ref(this%mass))

    ! mass of secondary component
    CALL GetAttr(config, "mass2", this%mass2, 1.0)
    CALL SetAttr(IO, "mass2", Ref(this%mass2))

    ! excentricity
    CALL GetAttr(config, "excentricity", this%excent, 0.0)

    ! semi major axis
    CALL GetAttr(config, "semimayoraxis", this%semaaxis, 1.0)

    ! Softening parameter
    CALL GetAttr(config, "softening1", this%eps1, 0.0)

    CALL GetAttr(config, "softening2", this%eps2, 0.0)

    ! soft switch on for primary component
    CALL GetAttr(config, "switchon1", this%switchon, -1.0)

    ! soft switch on for secondary component
    CALL GetAttr(config, "switchon2", this%switchon2, -1.0)

    ! rotating reference frame
    CALL GetAttr(config, "omega_rot", this%omega_rot, 0.0)

    ! check whether rotating mesh is activated
    CALL GetAttr(config, "mesh_rot", mesh_rot, 0)
    IF (mesh_rot.EQ.1) THEN
       this%omega_rot = Mesh%OMEGA
       this%r0(:) = Mesh%rotcent(:)
       IF (Mesh%omega.EQ.0.0) &
          CALL this%Warning("InitGravity_binary","Rotating mesh enabled but OMEGA set to zero in mesh")
    ELSE
       this%r0(:) = 0.0
    END IF

    ! reset mass flux
    this%mdot = 0.0
    this%massloss = 0.0
    this%accrate = 0.0

    ! set period
    this%period = 2.*PI*SQRT(this%semaaxis/(this%mass+this%mass2)/Physics%constants%GN)*this%semaaxis

    ! reset omega and omega**2
    this%omega  = 0.0
    this%omega2 = 0.0
    
    ! set ghost cell data to one to avoid zero division
    this%omega(Mesh%IGMIN:Mesh%IMIN+Mesh%IM1,:,:) = 1.0
    this%omega(Mesh%IMAX+Mesh%IP1:Mesh%IGMAX,:,:) = 1.0
    this%omega(:,Mesh%JGMIN:Mesh%JMIN+Mesh%JM1,:) = 1.0
    this%omega(:,Mesh%JMAX+Mesh%JP1:Mesh%JGMAX,:) = 1.0
    this%omega(:,:,Mesh%KGMIN:Mesh%KMIN+Mesh%KM1) = 1.0
    this%omega(:,:,Mesh%KMAX+Mesh%KP1:Mesh%KGMAX) = 1.0

    !initialise output
    CALL this%SetOutput(Mesh,Physics,config,IO)
    CALL this%InfoGravity()

  END SUBROUTINE InitGravity_binary


  !> \public write binary parameters to screen
  SUBROUTINE InfoGravity(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(IN) :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: mass1_str,mass2_str,excent_str,sema_str,omega_str
    !------------------------------------------------------------------------!
    WRITE (mass1_str,'(ES8.2)')   this%mass
    WRITE (mass2_str,'(ES8.2)')   this%mass2
    WRITE (excent_str,'(ES8.2)')  this%excent
    WRITE (sema_str,'(ES8.2)')    this%semaaxis
    CALL this%Info("            primary mass:      " // TRIM(mass1_str))
    CALL this%Info("            secondary mass:    " // TRIM(mass2_str))
    CALL this%Info("            excentricity:      " // TRIM(excent_str))
    CALL this%Info("            semi major axis:   " // TRIM(sema_str))
    IF(this%period.NE.0.0) THEN
        WRITE(omega_str,'(ES8.2)') 2.0*PI/this%period
        CALL this%Info("            pattern speed:     " // TRIM(omega_str))
    END IF
  END SUBROUTINE InfoGravity


  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(INOUT) :: this
    CLASS(mesh_base),      INTENT(IN)    :: Mesh
    CLASS(physics_base),   INTENT(IN)    :: Physics
    TYPE(Dict_TYP),        POINTER       :: config,IO
    !------------------------------------------------------------------------!
    INTEGER :: valwrite
    !------------------------------------------------------------------------!

    CALL this%gravity_pointmass%SetOutput(Mesh,Physics,config,IO)

    ! output cartesian positions of both components of the binary system
    CALL GetAttr(config, "output/binpos", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "binpos", this%pos)

  END SUBROUTINE SetOutput


  !> \public calculates the resultant gravitational accerleration of a binary system
  !!
  !! The acceleration is computed with respect to the curvilinear local frame of
  !! reference. Each component of the binary system contributes according to
  !! \f[
  !!     \vec{g} = -\frac{d\Phi}{d r} = -r \Omega^2 \hat{e}_r = -\Omega^2 \vec{r}
  !! \f]
  !! to the total acceleration. Thereby \f$ \Omega^2 = GM/r^3 \f$ is the square of
  !! the local Keplerian velocity.
  SUBROUTINE UpdateGravity_single(this,Mesh,Physics,Fluxes,pvar,time,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(INOUT) :: this
    CLASS(mesh_base),      INTENT(IN)    :: Mesh
    CLASS(physics_base),   INTENT(IN)    :: Physics
    CLASS(fluxes_base),    INTENT(IN)    :: Fluxes
    CLASS(marray_compound),INTENT(INOUT) :: pvar
    REAL,                  INTENT(IN)    :: time,dt
    !------------------------------------------------------------------------!
    REAL              :: GM1, GM2
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    ! update the positions of the binary system
    CALL this%UpdatePositions(time)

    ! compute curvilinear components of the position vectors
    ! 1. primary component
    this%posvec_prim_tmp(:,:,:,1) = this%pos(1,1)
    this%posvec_prim_tmp(:,:,:,2) = this%pos(2,1)
    this%posvec_prim_tmp(:,:,:,3) = this%pos(3,1)
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,EAST,:),this%posvec_prim_tmp(:,:,:,:),this%fposvec_prim(:,:,:,1,:))
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,NORTH,:),this%posvec_prim_tmp(:,:,:,:),this%fposvec_prim(:,:,:,2,:))
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,TOP,:),this%posvec_prim_tmp(:,:,:,:),this%fposvec_prim(:,:,:,3,:))
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,this%posvec_prim_tmp,this%posvec_prim)
    ! 2. secondary component
    this%posvec_sec_tmp(:,:,:,1) = this%pos(1,2)
    this%posvec_sec_tmp(:,:,:,2) = this%pos(2,2)
    this%posvec_sec_tmp(:,:,:,3) = this%pos(3,2)
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,EAST,:),this%posvec_sec_tmp(:,:,:,:),this%fposvec_sec(:,:,:,1,:))
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,NORTH,:),this%posvec_sec_tmp(:,:,:,:),this%fposvec_sec(:,:,:,2,:))
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%curv%faces(:,:,:,TOP,:),this%posvec_sec_tmp(:,:,:,:),this%fposvec_sec(:,:,:,3,:))
    CALL Mesh%Geometry%Convert2Curvilinear(Mesh%bcenter,this%posvec_sec_tmp,this%posvec_sec)

    ! subtract the result from the position vectors of all cell bary centers
    ! this gives you the curvilinear components of all vectors pointing
    ! from the point mass of each component of the binary system to the bary
    ! centers of all cells on the mesh
    this%posvec_prim(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - this%posvec_prim(:,:,:,:)
    this%fposvec_prim(:,:,:,1,:) = Mesh%posvec%faces(:,:,:,EAST,:) - this%fposvec_prim(:,:,:,1,:)   ! EAST
    this%fposvec_prim(:,:,:,2,:) = Mesh%posvec%faces(:,:,:,NORTH,:) - this%fposvec_prim(:,:,:,2,:)  ! NORTH
    this%fposvec_prim(:,:,:,3,:) = Mesh%posvec%faces(:,:,:,TOP,:) - this%fposvec_prim(:,:,:,3,:)    ! TOP

    this%posvec_sec(:,:,:,:) = Mesh%posvec%bcenter(:,:,:,:) - this%posvec_sec(:,:,:,:)
    this%fposvec_sec(:,:,:,1,:) = Mesh%posvec%faces(:,:,:,EAST,:) - this%fposvec_sec(:,:,:,1,:)     ! EAST
    this%fposvec_sec(:,:,:,2,:) = Mesh%posvec%faces(:,:,:,NORTH,:) - this%fposvec_sec(:,:,:,2,:)    ! NORTH
    this%fposvec_sec(:,:,:,3,:) = Mesh%posvec%faces(:,:,:,TOP,:) - this%fposvec_sec(:,:,:,3,:)      ! TOP

    ! compute the distances between each component of the binary system
    ! and all cell bary centers
    this%r_prim(:,:,:) = SQRT(this%posvec_prim(:,:,:,1)**2 &
      + this%posvec_prim(:,:,:,2)**2 + this%posvec_prim(:,:,:,3)**2)
    this%fr_prim(:,:,:,1) = SQRT(this%fposvec_prim(:,:,:,1,1)**2 &
      + this%fposvec_prim(:,:,:,1,2)**2 + this%fposvec_prim(:,:,:,1,3)**2) ! shifted EAST-faces
    this%fr_prim(:,:,:,2) = SQRT(this%fposvec_prim(:,:,:,2,1)**2 &
      + this%fposvec_prim(:,:,:,2,2)**2 + this%fposvec_prim(:,:,:,2,3)**2) ! shifted NORTH-faces
    this%fr_prim(:,:,:,3) = SQRT(this%fposvec_prim(:,:,:,3,1)**2 &
      + this%fposvec_prim(:,:,:,3,2)**2 + this%fposvec_prim(:,:,:,3,3)**2) ! shifted TOP-faces
    this%r_sec(:,:,:) = SQRT(this%posvec_sec(:,:,:,1)**2 &
      + this%posvec_sec(:,:,:,2)**2 + this%posvec_sec(:,:,:,3)**2)
    this%fr_sec(:,:,:,1) = SQRT(this%fposvec_sec(:,:,:,1,1)**2 &
      + this%fposvec_sec(:,:,:,1,2)**2 + this%fposvec_sec(:,:,:,1,3)**2)   ! shifted EAST-faces
    this%fr_sec(:,:,:,2) = SQRT(this%fposvec_sec(:,:,:,2,1)**2 &
      + this%fposvec_sec(:,:,:,2,2)**2 + this%fposvec_sec(:,:,:,2,3)**2)   ! shifted NORTH-faces
    this%fr_sec(:,:,:,3) = SQRT(this%fposvec_sec(:,:,:,3,1)**2 &
      + this%fposvec_sec(:,:,:,3,2)**2 + this%fposvec_sec(:,:,:,3,3)**2)   ! shifted NORTH-faces

    ! compute square of Keplerian velocities for updated time value
    ! with respect to PRIMARY star
    GM1 = Physics%Constants%GN*this%GetMass_primary(time)
    this%omega2(:,:,:,1) = GM1 / (this%r_prim(:,:,:)**3 + this%eps1**3)

    ! and SECONDARY star
    GM2 = Physics%Constants%GN*this%GetMass_secondary(time)
    this%omega2(:,:,:,2) = GM2 / (this%r_sec(:,:,:)**3 + this%eps2**3)

    ! potential calculation for binary and use pointmassroutine for it
    CALL this%gravity_pointmass%CalcPotential(Mesh,Physics,this%mass,this%r_prim,this%fr_prim,this%pot_prim%data4d)
    CALL this%gravity_pointmass%CalcPotential(Mesh,Physics,this%mass2,this%r_sec,this%fr_sec,this%pot_sec%data4d)

    this%pot%data1d(:) = this%pot_prim%data1d(:) + this%pot_sec%data1d(:)

    ! set curvilinear components of the gravitational acceleration
!NEC$ collapse
    DO k=Mesh%KGMIN,Mesh%KGMAX
      DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ ivdep
        DO i=Mesh%IGMIN,Mesh%IGMAX
          this%accel%data4d(i,j,k,1:Physics%VDIM) = -this%omega2(i,j,k,1) * this%posvec_prim(i,j,k,1:Physics%VDIM)&
                           -this%omega2(i,j,k,2) * this%posvec_sec(i,j,k,1:Physics%VDIM)
        END DO
      END DO
    END DO

  END SUBROUTINE UpdateGravity_single


  !> \private compute new positions for both components of the binary system
  !!
  !! This subroutine solves the Kepler equation to get the positions at a
  !! given time. A detailed explanation and a derivation of the relevant
  !! equations can be found in \cite taff1985 .
  !!
  !! First we map the dimensionless orbital time into the interval \f$[0,2\pi]\f$
  !! to calculate the mean anomaly \f$ \tau \f$.
  !! Then we have to solve the Kepler equation:
  !! \f[
  !!    \tau= E-\epsilon\sin{E}
  !! \f]
  !! to calculate the eccentric anomaly \f$ E \f$.
  !! \f$ \epsilon \f$ is the eccentricity of the orbit which is given as an
  !! input parameter of the binary.
  !! We can solve the equation using the function \ref roots.getroot_regulafalsi .
  !! To do this, we have to know an Intervall that includes the correct solution
  !! for \f$ E \f$. This should be as small as possible to speed up the calculation.
  !! Looking at the Kepler equation we can see, that if \f$ \tau\f$ is in
  !! the interval \f$[0,\pi]\f$ \f$E\f$ also has a value in  \f$[0,\pi]\f$.
  !! The same is true for the other possible intervall of values for \f$ \tau \f$
  !! and \f$ E \f$ \f$]\pi,2\pi[\f$.
  !! We can further constrain the Intervall for the correct solution for \f$ E \f$
  !! by looking at these separate cases, keeping in mind that the eccentricity
  !! \f$ \epsilon \f$ can only have values of:
  !! \f[
  !!    0 \leq \epsilon \leq 1
  !! \f]
  !! Case 1 : \f$ E,\tau \in [0,\pi] \f$
  !! \f{eqnarray*}{
  !!    E &=& \tau + \epsilon\sin{E} \leq \tau + \epsilon \\
  !!    E &=& \tau + \epsilon\sin{E} \geq \tau  \\
  !! \f}
  !! Case 2 : \f$ E,\tau \in ]\pi,2\pi[ \f$
  !! \f{eqnarray*}{
  !!    E &=& \tau + \epsilon\sin{E} < \tau  \\
  !!    E &=& \tau + \epsilon\sin{E} > \tau - \epsilon \\
  !! \f}
  !! If we know the mean anomaly \f$ E \f$ we can calculate the distance
  !! between the stars
  !! \f[
  !!    r = a (1-\epsilon \cos{E}).
  !! \f]
  !! \f$ a \f$ is the semimajor axis of the binary.
  !! We can also calculate the angle \f$ \phi \f$ between the semimajoraxis of
  !! the ellipse and the line connecting the primary component to the center of mass
  !! \f[
  !!    \phi = 2\arctan\left(\sqrt{\frac{1+\epsilon}{1-\epsilon}} \tan{\frac{E}{2}} \right)
  !! \f]
  !! This equation can be written as
  !! \f[
  !!    \phi = 2\arctan\left(\pm\sqrt{\frac{1+\epsilon}{1-\epsilon} \frac{1-\cos{E}}{1+\cos{E}}} \right)
  !! \f]
  !! to avoid one evaluation of \f$ \tan \f$.
  !! To choose the coorect sign in this expression we have to look at the sign of
  !! \f$ \tan{(E/2)} \f$ which is positive for  \f$ E\in[0,\pi] \f$ and negative for \f$ E\in[\pi,2\pi]\f$.
  !! When we know the angle \f$ \phi \f$ and the distance between the stars \f$ r \f$
  !! we can calculate the positions of the stars \f$ \vec{r}_1,\vec{r}_2 \f$ using
  !! the equations:
  !! \f{eqnarray*}{
  !!    \vec{r}_1 &=& \frac{m_2}{m_1+m_2} \vec{r}\\
  !!    \vec{r}_2 &=& \vec{r}-\vec{r}_1 \\
  !! \f}
  !! with the masses of the binary components \f$ m_1, m_2 \f$ and \f$ r_1=|\vec{r}_1|\f$.
  SUBROUTINE UpdatePositions(this,time)
    USE roots
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(INOUT) :: this
    REAL,                  INTENT(IN)    :: time
    !------------------------------------------------------------------------!
    REAL              :: E,cosE,r,phi,r1
    REAL              :: tau,excent
    REAL, DIMENSION(2):: plist
    INTEGER           :: err
    !------------------------------------------------------------------------!
    ! mean anomaly
    tau = 2.*PI*MODULO(time,this%period)/this%period
    excent = this%excent

    plist(1)=excent
    plist(2)=tau
    ! solve the Kepler equation: E-tau-excent*sin(E) = 0
    !   i.e. compute the intersection of f(E) = E-tau and g(E) = excent*sin(E)
    IF ((tau.GE.0.0) .AND. (tau.LE.PI)) THEN
       ! intersection lies ABOVE abscissa
       CALL GetRoot(func,MAX(0.0,tau),MIN(PI,excent+tau),E,err,plist)
    ELSE
       ! intersection lies BELOW abscissa
       CALL GetRoot(func,MAX(PI,tau-excent),MIN(2.*PI,tau),E,err,plist)
    END IF
    IF(err.NE.0) &
      CALL this%Error("UpdatePositions_binary","Error in root finding!")

    ! compute positions
    cosE = COS(E)
    ! distance to the origin of primary star position
    r = this%semaaxis * (1.0 - excent*cosE)
    ! position angle of primary star
    !   phi1 = 2.0*ATAN(SQRT((1.0+excent)/(1.0-excent))*TAN(0.5*E))
    !   A little bit more complicated but this avoids one evaluation of TAN;
    !   uses TAN(x/2) = +/- SQRT((1-cos(x))/(1+cos(x)))
    IF (cosE.NE.-1.0) THEN
       phi = 2.0*ATAN(SIGN(SQRT(((1.+excent)*(1.-cosE)) &
            /((1.-excent)*(1.+cosE))),PI-E))
    ELSE
       ! since lim_{x->0} atan(q/x) = pi/2 for q>0
       phi = PI
    END IF

    !transform to rotating reference frame if necessary
    phi = phi - this%omega_rot * time

    ! cartesian coordinates of SECONDARY component (right at t=0)
    ! \todo This implementation only works in 2D
    r1               = this%mass/(this%mass+this%mass2)*r
    this%pos(1,2) = r1*COS(phi)
    this%pos(2,2) = r1*SIN(phi)
    this%pos(3,2) = 0.0
    ! and PRIMARY COMPONENT
    this%pos(:,1) = -this%pos(:,2) * this%mass2/this%mass
    ! shift with respect to center of rotation
    this%pos(:,1) = this%pos(:,1) + this%r0(:)
    this%pos(:,2) = this%pos(:,2) + this%r0(:)

  END SUBROUTINE UpdatePositions

  !> \public compute the pressure scale height of a disk for the binary potential
  !!
  !! Since the potential of a binary system is given by the superposition
  !! of the potential of two individual stars we get
  !! \f[
  !!     h = c_s \left(\sum_{i=1}^2\frac{\partial^2\Phi_i}{\partial z^2}\right)^{-1/2}
  !!       = c_s \left(\Omega_1^2 + \Omega_2^2\right)^{-1/2}
  !! \f]
  !! where \f$ \Phi_i = -GM_i / r_i \f$ is the pointmass potential in the
  !! equatorial plane with respect to the \f$ i \f$ th component of the binary system
  !! and \f$ \Omega_i \f$ is the associated local Keplerian angular velocity.
  PURE SUBROUTINE CalcDiskHeight_single(this,Mesh,Physics,pvar,bccsound,h_ext,height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary),  INTENT(INOUT) :: this
    CLASS(mesh_base),       INTENT(IN)    :: Mesh
    CLASS(physics_base),    INTENT(IN)    :: Physics
    CLASS(marray_compound), INTENT(INOUT) :: pvar
    TYPE(marray_base),      INTENT(INOUT) :: bccsound,h_ext,height
    !------------------------------------------------------------------------!
    ! compute the effective omega
    this%omega(:,:,:) = SQRT(this%omega2(:,:,:,1) + this%omega2(:,:,:,2))
    ! compute the disk height
    height%data3d(:,:,:) = GetDiskHeight(bccsound%data3d(:,:,:),this%omega(:,:,:))
  END SUBROUTINE CalcDiskHeight_single


  !> get mass of the primary component, eventually scaled by switch on factor
  FUNCTION GetMass_primary(this,time) RESULT(mass)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(IN) :: this
    REAL,                  INTENT(IN) :: time
    !------------------------------------------------------------------------!
    REAL              :: mass
    !------------------------------------------------------------------------!
    IF(time.LE.this%switchon) THEN
        mass = this%mass * SIN(0.5*PI*time/this%switchon)**2
    ELSE
        mass = this%mass
    END IF
  END FUNCTION GetMass_primary

  !> get mass of the primary component, eventually scaled by switch on factor
  FUNCTION GetMass_secondary(this,time) RESULT(mass)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(gravity_binary), INTENT(IN) :: this
    REAL,                  INTENT(IN) :: time
    !------------------------------------------------------------------------!
    REAL              :: mass
    !------------------------------------------------------------------------!
    IF(time.LE.this%switchon2) THEN
        mass = this%mass2 * SIN(0.5*PI*time/this%switchon2)**2
    ELSE
        mass = this%mass2
    END IF
  END FUNCTION GetMass_secondary

  !> \public Closes the binary source term
  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(gravity_binary), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%r0,this%r_sec,this%posvec_sec,this%posvec_sec_tmp,&
               this%pot_prim,this%pot_sec,this%fr_sec,this%fposvec_sec,this%mass2)

    CALL this%Finalize_base()
  END SUBROUTINE Finalize


  !> find root of the Kepler equation
  PURE SUBROUTINE func(x,fx,plist)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist
    !plist[1]=excent, plist[2]=tau
    REAL, INTENT(OUT) :: fx
    !------------------------------------------------------------------------!
    fx  = x - plist(1)*SIN(x) - plist(2)
  END SUBROUTINE func

END MODULE gravity_binary_mod
