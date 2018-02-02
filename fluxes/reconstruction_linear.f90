!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: reconstruction_linear.f90                                         #
!#                                                                           #
!# Copyright (C) 2007-2017                                                   #
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

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Jannes Klee
!!
!! \brief module for linear (first order) TVD reconstruction using slope limiters
!!
!! Currentely supported limiter functions:
!! - minmod
!! - mc (monotonized central)
!! - sweby
!! - superbee
!! - ospre
!! - pp (positivity preserving)
!! - vanleer
!! - nolimit: for testing purposes only, disables limiting
!!
!! \extends reconstruction_common
!! \ingroup reconstruction
!----------------------------------------------------------------------------!
MODULE reconstruction_linear_mod
  USE reconstruction_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS (reconstruction_base)  :: reconstruction_linear
  CONTAINS
    PROCEDURE :: InitReconstruction_linear
    PROCEDURE :: CalculateStates
    PROCEDURE :: CalculateSlopes
    FINAL     :: Finalize
  END TYPE reconstruction_linear
  !--------------------------------------------------------------------------!
  INTEGER, PARAMETER :: MINMOD   = 1
  INTEGER, PARAMETER :: MONOCENT = 2
  INTEGER, PARAMETER :: SWEBY    = 3
  INTEGER, PARAMETER :: SUPERBEE = 4
  INTEGER, PARAMETER :: OSPRE    = 5
  INTEGER, PARAMETER :: PP       = 6
  INTEGER, PARAMETER :: VANLEER  = 7
  INTEGER, PARAMETER :: NOLIMIT  = 8
  CHARACTER(LEN=32), PARAMETER  :: recontype_name = "linear"
  CHARACTER(LEN=32), DIMENSION(8), PARAMETER :: limitertype_name = (/ &
         "minmod    ", "mc        ", "sweby     ", "superbee  ", "ospre     ", &
         "pp        ", "van leer  ", "no limit  "/)
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       reconstruction_linear, &
       ! constants
       MINMOD, MONOCENT, SWEBY, SUPERBEE, OSPRE, PP, VANLEER, NOLIMIT
  !--------------------------------------------------------------------------!

CONTAINS


  !> \public Constructor of linear reconstruction module
  SUBROUTINE InitReconstruction_linear(this,Mesh,Physics,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_linear), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_base),          INTENT(IN)    :: Physics
    TYPE(DICT_TYP),               POINTER       :: config
    TYPE(DICT_TYP),               POINTER       :: IO
    !------------------------------------------------------------------------!
    INTEGER                                     :: err, limiter
    REAL                                        :: theta
    !------------------------------------------------------------------------!
    ! allocate memory for all arrays used in reconstruction_linear
    ALLOCATE( &
         this%xslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                      Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),         &
         this%yslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                      Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),         &
         this%zslopes(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                      Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM),         &
         STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("InitReconstruction_linear",  "Unable to allocate memory.")
    END IF


    ! initialize parent
    CALL this%InitReconstruction(Mesh,Physics,config,IO,LINEAR,recontype_name)!Mesh,Physics,config,rtype,rname)

    ! set defaults for the limiter
    CALL GetAttr(config, "limiter", limiter, MINMOD)
    CALL GetAttr(config, "theta", theta, 1.0)

    ! limiter settings
    SELECT CASE(limiter)
    CASE(MINMOD,MONOCENT,SWEBY,SUPERBEE,OSPRE,PP,VANLEER,NOLIMIT)
       CALL this%InitLogging(limiter,limitertype_name(limiter))
       this%limiter_param= theta
    CASE DEFAULT
       CALL this%Error("InitReconstruction_linear", "Unknown limiter")
    END SELECT

    ! check parameter settings for PP limiter
    IF (limiter.EQ.PP) THEN
       IF (this%limiter_param.LT.0.0) &
          CALL this%Error("InitReconstruction_linear", &
               "Parameter of PP limiter must be greater than 0")
       IF (this%limiter_param.GT.EPSILON(this%limiter_param)) &
          CALL this%Warning("InitReconstruction_linear", &
               "Parameter of PP limiter should be less than machine precision.")
    END IF

    ! zero the slopes
    this%xslopes(:,:,:,:) = 0.0
    this%yslopes(:,:,:,:) = 0.0
    this%zslopes(:,:,:,:) = 0.0


    CALL this%Info("            limiter:           " // TRIM(this%GetName()))
  END SUBROUTINE InitReconstruction_linear


  !> \public Reconstructes states at cell boundaries
  PURE SUBROUTINE CalculateStates(this,Mesh,Physics,npos,dx,dy,dz,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_linear), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_base),          INTENT(IN)    :: Physics
    INTEGER                                     :: npos
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,npos) :: dx,dy,dz
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM)          &
                                                :: rvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,npos,Physics%VNUM)     &
                                                :: rstates
    !------------------------------------------------------------------------!
    INTEGER                                     :: i,j,k,l,n
    !------------------------------------------------------------------------!
    INTENT(IN)                                  :: npos,dx,dy,dz,rvar
    INTENT(OUT)                                 :: rstates
    !------------------------------------------------------------------------!
    ! calculate slopes first
    CALL this%CalculateSlopes(Mesh,Physics,rvar)

    ! reconstruct states at positions pos
    DO l=1,Physics%VNUM
      DO n=1,npos
!CDIR COLLAPSE
        DO k=Mesh%KGMIN,Mesh%KGMAX
          DO j=Mesh%JGMIN,Mesh%JGMAX
            DO i=Mesh%IGMIN,Mesh%IGMAX
!CDIR IEXPAND
              rstates(i,j,k,n,l) = reconstruct(rvar(i,j,k,l), &
                        this%xslopes(i,j,k,l),this%yslopes(i,j,k,l),this%zslopes(i,j,k,l), &
                        dx(i,j,k,n),dy(i,j,k,n),dz(i,j,k,n))
            END DO
          END DO
        END DO
      END DO
    END DO

    CONTAINS

      ELEMENTAL FUNCTION reconstruct(cvar0,xslope0,yslope0,zslope0,dx,dy,dz) RESULT(rstate)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: rstate
        REAL, INTENT(IN) :: cvar0,xslope0,yslope0,zslope0,dx,dy,dz
        !--------------------------------------------------------------------!
        rstate = cvar0 + xslope0*dx + yslope0*dy + zslope0*dz
      END FUNCTION reconstruct

  END SUBROUTINE CalculateStates


  !> \public computes the limited slopes
  PURE SUBROUTINE CalculateSlopes(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_linear), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_base),          INTENT(IN)    :: Physics
    REAL    :: rvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                    Mesh%KGMIN:Mesh%KGMAX,Physics%vnum)
    !------------------------------------------------------------------------!
    INTEGER                                     :: i,j,k,l
    !------------------------------------------------------------------------!
    INTENT(IN)                                  :: rvar
    !------------------------------------------------------------------------!

    ! choose limiter & check for cases where dimension needs to be neclected
!CDIR IEXPAND
    SELECT CASE(this%GetType())
    CASE(MINMOD)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR UNROLL=8
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * minmod2_limiter(&
                   rvar(i,j,k,l) - rvar(i-1,j,k,l), rvar(i+1,j,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%yslopes(i,j,k,l) = Mesh%invdy * minmod2_limiter(&
                  rvar(i,j,k,l) - rvar(i,j-1,k,l), rvar(i,j+1,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%zslopes(i,j,k,l) = Mesh%invdz * minmod2_limiter(&
                  rvar(i,j,k,l) - rvar(i,j,k-1,l), rvar(i,j,k+1,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
    CASE(MONOCENT)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR UNROLL=8
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * minmod3_limiter(&
                   this%limiter_param*(rvar(i,j,k,l) - rvar(i-1,j,k,l)),&
                   this%limiter_param*(rvar(i+1,j,k,l) - rvar(i,j,k,l)),&
                   0.5*(rvar(i+1,j,k,l) - rvar(i-1,j,k,l)))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
!CDIR COLLAPSE
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%yslopes(i,j,k,l) = Mesh%invdy * minmod3_limiter(&
                   this%limiter_param*(rvar(i,j,k,l)- rvar(i,j-1,k,l)),&
                   this%limiter_param*(rvar(i,j+1,k,l) - rvar(i,j,k,l)),&
                   0.5*(rvar(i,j+1,k,l) - rvar(i,j-1,k,l)))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%zslopes(i,j,k,l) = Mesh%invdz * minmod3_limiter(&
                   this%limiter_param*(rvar(i,j,k,l)- rvar(i,j,k-1,l)),&
                   this%limiter_param*(rvar(i,j,k+1,l) - rvar(i,j,k,l)),&
                   0.5*(rvar(i,j,k+1,l) - rvar(i,j,k-1,l)))
              END DO
            END DO
          END DO
        END DO
      END IF
    CASE(SWEBY)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR UNROLL=8
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * sweby_limiter(&
                   rvar(i,j,k,l) - rvar(i-1,j,k,l), rvar(i+1,j,k,l) - rvar(i,j,k,l),&
                   this%limiter_param)
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
!CDIR COLLAPSE
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                 this%yslopes(i,j,k,l) = Mesh%invdy * sweby_limiter(&
                    rvar(i,j,k,l) - rvar(i,j-1,k,l), rvar(i,j+1,k,l) - rvar(i,j,k,l),&
                    this%limiter_param)
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%zslopes(i,j,k,l) = Mesh%invdz * sweby_limiter(&
                   rvar(i,j,k,l) - rvar(i,j,k-1,l), rvar(i,j,k+1,l) - rvar(i,j,k,l),&
                   this%limiter_param)
              END DO
            END DO
          END DO
        END DO
      END IF
    CASE(SUPERBEE)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
!CDIR UNROLL=8
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * sweby_limiter(&
                   rvar(i,j,k,l) - rvar(i-1,j,k,l), rvar(i+1,j,k,l) - rvar(i,j,k,l), 2.0)
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                 this%yslopes(i,j,k,l) = Mesh%invdy * sweby_limiter(&
                    rvar(i,j,k,l) - rvar(i,j-1,k,l), rvar(i,j+1,k,l) - rvar(i,j,k,l), 2.0)
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%zslopes(i,j,k,l) = Mesh%invdz * sweby_limiter(&
                   rvar(i,j,k,l) - rvar(i,j,k-1,l), rvar(i,j,k+1,l) - rvar(i,j,k,l), 2.0)
              END DO
            END DO
          END DO
        END DO
      END IF
    CASE(OSPRE)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR UNROLL=8
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR COLLAPSE
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * ospre_limiter(&
                   rvar(i,j,k,l) - rvar(i-1,j,k,l), rvar(i+1,j,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                 this%yslopes(i,j,k,l) = Mesh%invdy * ospre_limiter(&
                    rvar(i,j,k,l) - rvar(i,j-1,k,l), rvar(i,j+1,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                 this%zslopes(i,j,k,l) = Mesh%invdz * ospre_limiter(&
                    rvar(i,j,k,l) - rvar(i,j,k-1,l), rvar(i,j,k+1,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
!      CASE(PP)
!        ! calculate slopes in both-directions
!        DO l=1,Physics%VNUM
!!CDIR UNROLL=8
!          DO k=Mesh%KGMIN,Mesh%KGMAX
!            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
!!CDIR NODEP
!              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
!                 CALL pp_limiter(this%xslopes(i,j,k,l),&
!                                 this%yslopes(i,j,k,l),&
!                    rvar(i-1,j+1,k)-rvar(i,j,k),&
!                    rvar(i  ,j+1,k)-rvar(i,j,k),&
!                    rvar(i+1,j+1,k)-rvar(i,j,k),&
!                    rvar(i-1,j  ,k)-rvar(i,j,k),&
!                    this%limiter_param,&
!                    rvar(i+1,j  ,k)-rvar(i,j,k),&
!                    rvar(i-1,j-1,k)-rvar(i,j,k),&
!                    rvar(i  ,j-1,k)-rvar(i,j,k),&
!                    rvar(i+1,j-1,k)-rvar(i,j,k))
!              END DO
!            END DO
!          END DO
!        END DO
    CASE(VANLEER)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR UNROLL=8
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR COLLAPSE
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * vanleer_limiter(&
                   rvar(i,j,k,l) - rvar(i-1,j,k,l), rvar(i+1,j,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
       ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%yslopes(i,j,k,l) = Mesh%invdy * vanleer_limiter(&
                   rvar(i,j,k,l) - rvar(i,j-1,k,l), rvar(i,j+1,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%zslopes(i,j,k,l) = Mesh%invdz * vanleer_limiter(&
                   rvar(i,j,k,l) - rvar(i,j,k-1,l), rvar(i,j,k+1,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
    CASE(NOLIMIT)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR UNROLL=8
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
              DO i=Mesh%IGMIN+1,Mesh%IGMAX-1
                this%xslopes(i,j,k,l) = Mesh%invdx * nolimit_limiter(&
                   rvar(i,j,k,l) - rvar(i-1,j,k,l), rvar(i+1,j,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN+1,Mesh%JGMAX-1
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%yslopes(i,j,k,l) = Mesh%invdy * nolimit_limiter(&
                   rvar(i,j,k,l) - rvar(i,j-1,k,l), rvar(i,j+1,k,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        DO l=1,Physics%VNUM
!CDIR COLLAPSE
          DO k=Mesh%KGMIN+1,Mesh%KGMAX-1
            DO j=Mesh%JGMIN,Mesh%JGMAX
              DO i=Mesh%IGMIN,Mesh%IGMAX
                this%zslopes(i,j,k,l) = Mesh%invdz * nolimit_limiter(&
                   rvar(i,j,k,l) - rvar(i,j,k-1,l), rvar(i,j,k+1,l) - rvar(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
      END IF
    END SELECT

    CONTAINS

      ELEMENTAL FUNCTION minmod2_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        IF (SIGN(1.0,arg1)*SIGN(1.0,arg2).GT.0) THEN
           limarg = SIGN(MIN(ABS(arg1),ABS(arg2)),arg1)
        ELSE
           limarg = 0.
        END IF
      END FUNCTION minmod2_limiter

      ELEMENTAL FUNCTION minmod3_limiter(arg1,arg2,arg3) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2, arg3
        !--------------------------------------------------------------------!
        IF (((SIGN(1.0,arg1)*SIGN(1.0,arg2)).GT.0).AND.&
            ((SIGN(1.0,arg2)*SIGN(1.0,arg3)).GT.0)) THEN
           limarg = SIGN(MIN(ABS(arg1),ABS(arg2),ABS(arg3)),arg1)
        ELSE
           limarg = 0.
        END IF
      END FUNCTION minmod3_limiter

      ELEMENTAL FUNCTION sweby_limiter(arg1,arg2,param) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2, param
        !--------------------------------------------------------------------!
        IF (SIGN(1.0,arg1)*SIGN(1.0,arg2).GT.0) THEN
           limarg = SIGN(MAX(MIN(param*ABS(arg1),ABS(arg2)), &
                MIN(ABS(arg1),param*ABS(arg2))),arg1)
        ELSE
           limarg = 0.
        END IF
      END FUNCTION sweby_limiter

      ELEMENTAL FUNCTION ospre_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        limarg = 1.5*arg1*arg2*(arg1 + arg2) / (arg1*(arg1 + 0.5*arg2) &
             + arg2*(arg2 + 0.5*arg1) + TINY(1.0))
      END FUNCTION ospre_limiter

!      ELEMENTAL SUBROUTINE pp_limiter(xslope,yslope,arg1,arg2,arg3,arg4,param,arg6,arg7,arg8,arg9)  ! TODO: 3D
!        IMPLICIT NONE
!        !--------------------------------------------------------------------!
!        REAL, INTENT(OUT):: xslope,yslope,zslope
!        REAL, INTENT(IN) :: arg1,arg2,arg3,arg4,param,arg6,arg7,arg8,arg9
!        !--------------------------------------------------------------------!
!        REAL             :: Vmin,Vmax,V
!        !--------------------------------------------------------------------!
!        xslope = (arg6-arg4)*0.5
!        yslope = (arg2-arg8)*0.5 ! TODO: 3D add zslope
!        Vmin = MIN(arg1,arg2,arg3,arg4,-param,arg6,arg7,arg8,arg9)
!        Vmax = MAX(arg1,arg2,arg3,arg4,+param,arg6,arg7,arg8,arg9)
!        V = MIN(1.0,2.0*MIN(ABS(Vmin),ABS(Vmax))/(ABS(xslope)+ABS(yslope)))
!        xslope = V*xslope*Mesh%invdx
!        yslope = V*yslope*Mesh%invdy ! TODO: 3D add zslope
!      END SUBROUTINE pp_limiter

      ELEMENTAL FUNCTION vanleer_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        REAL             :: a1,a2
        !--------------------------------------------------------------------!
        a1 = ABS(arg1)
        a2 = ABS(arg2)
        limarg = (arg1 * a2 + arg2 * a1)/(a1 + a2 + TINY(a1))
      END FUNCTION vanleer_limiter

      ELEMENTAL FUNCTION nolimit_limiter(arg1,arg2) RESULT(limarg)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        REAL             :: limarg
        REAL, INTENT(IN) :: arg1, arg2
        !--------------------------------------------------------------------!
        limarg = 0.5*(arg1+arg2)
      END FUNCTION nolimit_limiter

  END SUBROUTINE CalculateSlopes


  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(reconstruction_linear), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%xslopes,this%yslopes,this%zslopes)

    CALL this%FinalizeReconstruction()
  END SUBROUTINE Finalize

END MODULE reconstruction_linear_mod
