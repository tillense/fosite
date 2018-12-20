!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: sources_viscosity.f90                                             #
!#                                                                           #
!# Copyright (C) 2008-2012                                                   #
!# Bjoern Sperling <sperling@astrophysik.uni-kiel.de>                        #
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
!# along with Mesh program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!> \addtogroup sources
!! - parameters of \link sources_viscosity \endlink as key-values
!! \key{vismodel,INTEGER,viscosity model}
!! \key{dynconst,REAL,dynamic viscosity constant,0.1}
!! \key{bulkconst,REAL,bulk viscosity constant,-2/3*dynconst}
!! \key{cvis,REAL,viscous courant number,0.5}
!! \key{output/stress,INTEGER,enable(=1) output of the stress tensor,0}
!! \key{output/dynvis,INTEGER,enable(=1) output of dynamic viscosity,0}
!! \key{output/kinvis,INTEGER,enable(=1) output of kinematic viscosity,0}
!! \key{output/bulkvis,INTEGER,enable(=1) output of bulk viscosity,0}
!----------------------------------------------------------------------------!
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief viscosity of Newtonian fluid
!----------------------------------------------------------------------------!
MODULE sources_viscosity_mod
  USE sources_base_mod
  USE physics_base_mod
  USE fluxes_base_mod
  USE mesh_base_mod
  USE marray_compound_mod
  USE logging_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
    ! flags for viscosity model
    INTEGER,PARAMETER :: MOLECULAR = 1     ! constant viscosity
    INTEGER,PARAMETER :: ALPHA     = 2     ! Shakura-Sunyaev prescription
    INTEGER,PARAMETER :: BETA      = 3     ! Duschl prescription
    INTEGER,PARAMETER :: PRINGLE   = 4     ! constant kinematic viscosity
    INTEGER,PARAMETER :: ALPHA_ALT = 5     ! alternative Shakura-Sunyaev
    CHARACTER(LEN=32),PARAMETER, DIMENSION(5) :: viscosity_name = (/ &
                                     "constant viscosity              ", &
                                     "turbulent Shakura-Sunyaev       ", &
                                     "turbulent Duschl                ", &
                                     "const. kinematic viscosity      ", &
                                     "alternative Shakura-Sunyaev     "/)


  TYPE, EXTENDS(sources_base) :: sources_viscosity
    CHARACTER(LEN=32) :: source_name = "viscosity of Newtonian fluid"
    CLASS(logging_base), ALLOCATABLE :: viscosity
    REAL :: dynconst,bulkconst                      !< viscosity const.
    REAL, DIMENSION(:,:,:), POINTER  :: dynvis, &   !< dynamic viscosity
                                        kinvis, &   !< kinematic viscosity
                                        bulkvis, &  !< bulk viscosity
                                        envelope
     REAL, DIMENSION(:,:,:), POINTER :: &           !< comp. of stress tensor
           btxx,btyy,btzz,btxy,btxz,btyz,tmp,tmp2,tmp3
     REAL, DIMENSION(:,:), POINTER   :: &
           Sxx,Syy,Szz,Sxy,Sxz,Syz
  CONTAINS
    PROCEDURE :: InitSources_viscosity
    PROCEDURE :: InfoSources
    PROCEDURE :: SetOutput
    PROCEDURE :: UpdateViscosity
    PROCEDURE :: ExternalSources_single
    PROCEDURE :: CalcTimestep_single

    PROCEDURE :: Finalize
  END TYPE
  !--------------------------------------------------------------------------!
 PUBLIC :: &
       ! types
       Sources_viscosity, &
       ! constants
       MOLECULAR,ALPHA,BETA,PRINGLE,ALPHA_ALT,viscosity_name
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources_viscosity(this,Mesh,Physics,Fluxes,config,IO)
    USE geometry_cylindrical_mod, ONLY: geometry_cylindrical
    USE geometry_logcylindrical_mod, ONLY: geometry_logcylindrical
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_viscosity) :: this
    CLASS(Mesh_base)    :: Mesh
    CLASS(Physics_base) :: Physics
    CLASS(Fluxes_base)  :: Fluxes
    TYPE(Dict_TYP),POINTER :: config,IO
    INTEGER           :: stype
    !------------------------------------------------------------------------!
    INTEGER           :: err, viscosity_number
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,Fluxes
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "stype", stype)
    CALL this%InitLogging(stype,this%source_name)
    ! viscosity model
    CALL GetAttr(config, "vismodel", viscosity_number)
    ALLOCATE(logging_base::this%viscosity)
    CALL this%viscosity%InitLogging(viscosity_number,viscosity_name(viscosity_number))


    IF (.NOT.Fluxes%Initialized()) &
         CALL this%Error("InitSources_viscosity","fluxes module uninitialized")

    ! check mesh
    IF (Mesh%GetType().NE.MIDPOINT) &
         CALL this%Error("InitSources_viscosity", &
         "only midpoint rule is currently supported")


    ! dynamic viscosity constant
    CALL GetAttr(config, "dynconst", this%dynconst, 0.1)
    ! bulk viscosity constant
    ! At the moment this only affects the molecular viscosity. The default
    ! value refers to a trace free viscosity tensor.
    ! All other viscosity models should always be trace free. Therefore this
    ! is hard coded and cannot be changed.
    CALL GetAttr(config, "bulkconst", this%bulkconst, -2./3.*this%dynconst)

    ! set viscous courant number
    CALL GetAttr(config, "cvis", this%cvis, 0.5)

    ! enable the envelope mask field, which scales the viscosity with [0,1]
    CALL GetAttr(config, "envelope", this%use_envelope, 0)
    ALLOCATE(this%dynvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%kinvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%bulkvis(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%btxx(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%btyy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%btzz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%btxy(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%btxz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         this%btyz(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
         STAT=err)
    IF (err.NE.0) &
         CALL this%Error("InitSources_viscosity","Memory allocation failed.")
    ! do the model specific initialization
    SELECT CASE(this%viscosity%GetType())
    CASE(MOLECULAR,PRINGLE)
       ! do nothing
    CASE(ALPHA)
       ALLOCATE(this%invr(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX), &
            STAT=err)
       IF (err.NE.0) CALL this%Error("InitSources_viscosity",&
               "Unable to allocate memory.")
       ! check geometry
       SELECT TYPE(geo=>Mesh%geometry)
       TYPE IS(geometry_cylindrical)
          ! compute inverse of distance to origin
          this%invr(:,:,:) = 1./(TINY(1.0)+Mesh%radius%bcenter(:,:,:))
       TYPE IS(geometry_logcylindrical)
          this%invr(:,:,:) = 1./(TINY(1.0)+Mesh%radius%bcenter(:,:,:))
!      TYPE IS(TANCYLINDRICAL,SPHERICAL,OBLATE_SPHEROIDAL)
!          ! compute inverse of distance to axis
!          this%invr(:,:,:) = 1./(TINY(1.0)+Mesh%hz%bcenter(:,:,:))
       CLASS DEFAULT
          CALL this%Error("InitSources_viscosity",&
               "Geometry not supported for alpha-viscosity")
       END SELECT
    CASE(BETA)
!       SELECT CASE(Physics%GetType())
!       CASE (EULER2D,EULER2D_ISOTHERM,&
!             EULER2D_IAMT,EULER2D_ISOIAMT,&
!             EULER3D_ROTSYM,EULER3D_ROTAMT,EULER3D_ROTSYMSGS,EULER3D,EULER3D_ISOTHERM)
!          ! do nothing
!       CASE DEFAULT
          CALL this%Error("InitSources_viscosity",&
               "Physics not supported for beta-viscosity")
!       END SELECT
    CASE(ALPHA_ALT)
       !> \todo check if this is really sufficient
       IF (Physics%VDIM.NE.2) &
          CALL this%Error("InitSources_viscosity",&
               "alternative alpha-viscosity works only for flat disks")
   !    IF (.NOT.Timedisc%always_update_bccsound) &
   !       CALL Error(this,"InitSources_viscosity","always_update_bccsound must be enabled in timedisc")
    CASE DEFAULT
       CALL this%Error("InitSources_viscosity",&
            "Viscosity prescription not supported.")
    END SELECT
    ! set initial time to negative value;
    ! this guarantees update of viscosity for time t=0.0
    this%time = -1.0
    ! initialize viscosity arrays
    this%dynvis(:,:,:)  = this%dynconst
    this%bulkvis(:,:,:) = this%bulkconst
    ! set stress tensor arrays to zero
    this%btxx(:,:,:) = 0.0
    this%btyy(:,:,:) = 0.0
    this%btzz(:,:,:) = 0.0
    this%btxy(:,:,:) = 0.0
    this%btxz(:,:,:) = 0.0
    this%btyz(:,:,:) = 0.0
    IF(this%use_envelope.EQ.1) THEN
      ALLOCATE(this%envelope(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX),&
               STAT=err)
      IF (err.NE.0) &
        CALL this%Error("InitSources_viscosity","Memory allocation failed.")
      this%envelope = 1.
    END IF
    CALL this%SetOutput(Mesh,Physics,config,IO)
    IF (this%GetType().EQ.ALPHA_ALT) this%update_disk_height = .True.

    CALL this%InitSources(Mesh,Fluxes,Physics,config,IO)

  END SUBROUTINE InitSources_viscosity

  SUBROUTINE SetOutput(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_viscosity),INTENT(IN) :: this
    CLASS(Mesh_base) ,INTENT(IN)      :: Mesh
    CLASS(Physics_base),INTENT(IN)    :: Physics
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    INTEGER              :: valwrite
    !------------------------------------------------------------------------!
   ! Hier wird der Stresstensor per default in 3D ausgegeben! Vorher nur in 2D plus abhängig
   ! von der Physik auch die Komponenten der dritten Dimension
    CALL GetAttr(config, "output/stress", valwrite, 0)
    IF (valwrite .EQ. 1) THEN
       CALL SetAttr(IO, "stress:xx", this%btxx(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
       CALL SetAttr(IO, "stress:xy", this%btxy(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
       CALL SetAttr(IO, "stress:xz", this%btxz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
       CALL SetAttr(IO, "stress:yy", this%btyy(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
       CALL SetAttr(IO, "stress:yz", this%btyz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))
       CALL SetAttr(IO, "stress:zz", this%btzz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))


!     SELECT CASE(GetType(Physics))
!       CASE (EULER3D_ROTSYM,EULER3D_ROTAMT,EULER3D_ROTSYMSGS)
!       CALL IO%SetAttr("stress:xz", this%btxz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAXX,Mesh%KMIN:Mesh%KMAX))
!       CALL IO%SetAttr("stress:yz", this%btyz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAXX,Mesh%KMIN:Mesh%KMAX))
!       CALL IO%SetAttr("stress:zz", this%btzz(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX))
!      END SELECT
    END IF

    CALL GetAttr(config, "output/dynvis", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO,"dynvis", this%dynvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    CALL GetAttr(config, "output/kinvis", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "kinvis", this%kinvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

    CALL GetAttr(config, "output/bulkvis", valwrite, 0)
    IF (valwrite .EQ. 1) &
         CALL SetAttr(IO, "bulkvis", this%bulkvis(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX))

  END SUBROUTINE SetOutput


  SUBROUTINE InfoSources(this,Mesh)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_viscosity),INTENT(IN) :: this
    CLASS(Mesh_base),INTENT(IN)         :: Mesh
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: dynconst_str,bulkconst_str
    !------------------------------------------------------------------------!
    WRITE (dynconst_str,'(ES9.2)') this%dynconst
    WRITE (bulkconst_str,'(ES9.2)') this%bulkconst
    CALL this%Info("            viscosity model:   " // this%viscosity%GetName())
    SELECT CASE(this%viscosity%GetType())
    CASE(MOLECULAR)
       CALL this%Info("            dynamic viscosity: " // TRIM(dynconst_str))
       CALL this%Info("            bulk viscosity:    " // TRIM(bulkconst_str))
    CASE(ALPHA,ALPHA_ALT)
       CALL this%Info("            alpha:             " // TRIM(dynconst_str))
    CASE(BETA)
       CALL this%Info("            beta:              " // TRIM(dynconst_str))
    CASE(PRINGLE)
       CALL this%Info("            kinemat. viscosity:" // TRIM(dynconst_str))
    END SELECT
  END SUBROUTINE InfoSources


  SUBROUTINE UpdateViscosity(this,Mesh,Physics,Fluxes,time,pvar,cvar)
    USE physics_eulerisotherm_mod, ONLY : physics_eulerisotherm
    USE physics_euler_mod, ONLY : physics_euler
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_viscosity),INTENT(INOUT) :: this
    CLASS(Mesh_base),INTENT(IN)         :: Mesh
    CLASS(Physics_base),INTENT(INOUT)   :: Physics
    CLASS(Fluxes_base),INTENT(IN)       :: Fluxes
    REAL,INTENT(IN)                     :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM) &
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    INTEGER           :: kv=0
    !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,cvar
    !------------------------------------------------------------------------!
    ! only update viscosity if time has changed
    IF ((time.NE.this%time) .OR. (time .EQ. 0.)) THEN
       ! Set the dynamic viscosity
       SELECT CASE(this%viscosity%GetType())
       CASE(MOLECULAR)
          ! do nothing, since it has already been initialized
          ! in InitSources_viscosity
       CASE(ALPHA)
          ! Shakura-Sunyaev type alpha viscosity
          ! standard alpha prescription: nu = alpha*cs*h
          ! or nu = alpha / A * cs**2 / omega with A = -d ln(omega)/d ln(r)
          ! (see Kato, Fukue & Minishige: Black Hole Accretion Disks, 2008;
          ! equation (3.46))
          !
          ! this is a rough estimation assuming that the logarithmic derivative
          ! of the angular velocity (d ln(omega)/d ln(r)) is of the order of -1
          !
          ! check for 2D / 3D
          !> \todo check if this is really sufficient
          SELECT CASE(Physics%VDIM)
          CASE(2)
             kv  = Physics%YVELOCITY
          CASE(3)
             kv  = Physics%ZVELOCITY
          END SELECT

          ! check for non-isothermal physics
          SELECT TYPE(phys => Physics)
          CLASS IS(physics_euler)
            ! compute alpha viscosity
            this%dynvis(:,:,:) = etafkt_alpha(this%dynconst,pvar(:,:,:,phys%PRESSURE), &
                                    pvar(:,:,:,kv)*this%invr(:,:,:) + Mesh%Omega)
          CLASS IS(physics_eulerisotherm)
            ! compute alpha viscosity
            this%dynvis(:,:,:) = etafkt_alpha(this%dynconst, &
                                    phys%bccsound%data3d(:,:,:)**2*pvar(:,:,:,phys%DENSITY), &
                                    pvar(:,:,:,kv)*this%invr(:,:,:) + Mesh%Omega)            
          END SELECT

       CASE(BETA)
          ! Duschl type beta viscosity
          !TODO Beta-viscosity is NOT converted to 3D!!!!!
          !an usage of beta-viscosity with a cartesian grid will be wrong!!!!!
          CALL this%ERROR("UpdateViscosity","Beta viscosity not supported at the moment")
!          SELECT CASE(Physics%GetType())
!          CASE (EULER2D,EULER2D_ISOTHERM)
!             this%dynvis(:,:,:) = etafkt_beta(this%dynconst, &
!                  Mesh%hy%bcenter(:,:,:) &
!                    * (cvar(:,:,:,Physics%YMOMENTUM) &
!                       +cvar(:,:,:,Physics%DENSITY)*Mesh%hy%bcenter(:,:,:)*Mesh%omega))
!          CASE (EULER2D_IAMT,EULER2D_ISOIAMT)
!             this%dynvis(:,:,:) = etafkt_beta(this%dynconst, &
!                  cvar(:,:,:,Physics%YMOMENTUM))
!          CASE (EULER3D_ROTSYM,EULER3D_ROTSYMSGS)
!             this%dynvis(:,:,:) = etafkt_beta(this%dynconst, &
!                  Mesh%hz%bcenter(:,:,:)*cvar(:,:,:,Physics%ZMOMENTUM))
!          CASE (EULER3D_ROTAMT)
!             this%dynvis(:,:,:) = etafkt_beta(this%dynconst, &
!                  cvar(:,:,:,Physics%ZMOMENTUM))
!          END SELECT
       CASE(PRINGLE) ! constant kinematic viscosity
          this%dynvis(:,:,:) = etafkt_pringle(this%dynconst,pvar(:,:,:,Physics%DENSITY))
! TODO: height of the disk is calculated in another source module, which is not translated yet. Therefore, ALPHA_ALT is not
! working
!        CASE(ALPHA_ALT)
!          ! eta = nu*rho = alpha*cs*h * rho
!          ! compute dynamic viscosity
!          this%dynvis(:,:,:) = this%dynconst * Physics%bccsound%data3d(:,:,:) &
!                * Physics%sources%height(:,:,:) * pvar(:,:,:,Physics%DENSITY)
       END SELECT

       IF(this%use_envelope.EQ.1) &
         this%dynvis(:,:,:) = this%envelope(:,:,:) * this%dynvis(:,:,:)

       ! Set the bulk viscosity
       SELECT CASE(this%viscosity%GetType())
       CASE(MOLECULAR)
         ! do nothing, since it has already been initialized
         ! in InitSources_viscosity
       CASE(ALPHA,BETA,PRINGLE,ALPHA_ALT)
         ! All available parametrized viscosities should be trace free.
         ! This can be achieved by setting the following:
         this%bulkvis(:,:,:) = -2./3.*this%dynvis(:,:,:)
       END SELECT

       ! update time
       this%time = time
    END IF
  CONTAINS
    ! some elemental functions for computation of the viscosity

    ! Alpha viscosity
    ELEMENTAL FUNCTION etafkt_alpha(alpha,P,omega) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: alpha,P,omega
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = alpha * P / ABS(omega)
    END FUNCTION etafkt_alpha

    ! Beta viscosity
    ELEMENTAL FUNCTION etafkt_beta(beta,L) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: beta,L
      REAL :: eta
      !----------------------------------------------------------------------!
      ! L = r * rho * v_phi is the angular momentum
      eta = beta * abs(L)
    END FUNCTION etafkt_beta

    ! Pringle disk (kinematic viscosity is constant)
    ELEMENTAL FUNCTION etafkt_pringle(nu,rho) RESULT(eta)
      IMPLICIT NONE
      !----------------------------------------------------------------------!
      REAL, INTENT(IN) :: nu,rho
      REAL :: eta
      !----------------------------------------------------------------------!
      eta = nu * rho   !nu is the kinematic viscosity
    END FUNCTION etafkt_pringle

  END SUBROUTINE UpdateViscosity


  SUBROUTINE ExternalSources_single(this,Mesh,Physics,Fluxes,Sources,time,dt,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(sources_viscosity),INTENT(INOUT) :: this
    CLASS(mesh_base),INTENT(IN)            :: Mesh
    CLASS(physics_base),INTENT(INOUT)      :: Physics
    CLASS(fluxes_base),INTENT(IN)          :: Fluxes
    CLASS(sources_base), INTENT(INOUT)  :: Sources
    REAL,INTENT(IN)                        :: time, dt
    CLASS(marray_compound),INTENT(INOUT)   :: pvar,cvar,sterm
    !------------------------------------------------------------------------!
    CALL this%UpdateViscosity(Mesh,Physics,Fluxes,time,pvar%data4d,cvar%data4d)
    CALL Physics%CalcStresses_euler(Mesh,pvar%data4d,this%dynvis,this%bulkvis, &
             this%btxx,this%btxy,this%btxz,this%btyy,this%btyz,this%btzz)
    CALL Physics%ViscositySources(Mesh,pvar%data4d,this%btxx,this%btxy,this%btxz, &
             this%btyy,this%btyz,this%btzz,sterm%data4d)
  END SUBROUTINE ExternalSources_single


  SUBROUTINE CalcTimestep_single(this,Mesh,Physics,Fluxes,time,pvar,cvar,dt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_viscosity),INTENT(INOUT) :: this
    CLASS(Mesh_base),INTENT(IN)         :: Mesh
    CLASS(Physics_base),INTENT(INOUT)   :: Physics
    CLASS(Fluxes_base),INTENT(IN)       :: Fluxes
    REAL,INTENT(IN)                     :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM) &
                      :: pvar,cvar
    REAL              :: dt
    !------------------------------------------------------------------------!
    REAL              :: invdt_x, invdt_y,invdt_z
    !------------------------------------------------------------------------!
    INTENT(IN)        :: pvar,cvar
    INTENT(OUT)       :: dt
    !------------------------------------------------------------------------!
    ! update dynamic and bulk viscosity
    CALL this%UpdateViscosity(Mesh,Physics,Fluxes,time,pvar,cvar)

    ! compute kinematic viscosity
    this%kinvis(:,:,:) = this%dynvis(:,:,:) / pvar(:,:,:,Physics%DENSITY)

    ! x-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_x = MAXVAL(this%kinvis(:,:,:) / Mesh%dlx%data3d(:,:,:)**2, &
                        MASK=Mesh%without_ghost_zones%mask3d(:,:,:))
    ELSE
       ! set to zero, i.e. no limit in x-direction due to diffusion
       invdt_x = 0.0
    END IF
    ! y-direction
    IF (Mesh%JNUM.GT.1) THEN
       invdt_y = MAXVAL(this%kinvis(:,:,:)  / Mesh%dly%data3d(:,:,:)**2, &
                        MASK=Mesh%without_ghost_zones%mask3d(:,:,:))
    ELSE
       ! set to zero, i.e. no limit in y-direction due to diffusion
       invdt_y = 0.0
    END IF
    ! z-direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_z = MAXVAL(this%kinvis(:,:,:)  / Mesh%dlz%data3d(:,:,:)**2, &
                        MASK=Mesh%without_ghost_zones%mask3d(:,:,:))
    ELSE
       ! set to zero, i.e. no limit in z-direction due to diffusion
       invdt_z = 0.0
    END IF
    ! largest time step due to diffusion
    invdt_x = MAX(invdt_x, invdt_y,invdt_z)
    IF (invdt_x.GT.TINY(invdt_x)) THEN
       dt = this%cvis / invdt_x
    ELSE
       dt = HUGE(dt)
    END IF
  END SUBROUTINE CalcTimestep_single


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Sources_viscosity),INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%dynvis,this%kinvis,this%bulkvis, &
               this%btxx,this%btyy,this%btzz,this%btxy,this%btxz,this%btyz)
    IF(this%viscosity%GetType().EQ.ALPHA) DEALLOCATE(this%invr)
    IF(this%use_envelope.EQ.1) DEALLOCATE(this%envelope)
    CALL this%Finalize_base()
    IF(ASSOCIATED(this%next)) CALL this%next%Finalize()
  END SUBROUTINE Finalize


END MODULE sources_viscosity_mod
