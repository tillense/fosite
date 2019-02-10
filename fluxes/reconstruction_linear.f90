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
!! - pp (positivity preserving), NOT SUPPORTED in Fosite3D
!! - vanleer
!! - nolimit: for testing purposes only, disables limiting
!!
!! \todo implement 1D/2D/3D PP limiter
!! \todo fix broken output of slopes, uncomment the lines below yields segfault
!!
!! \extends reconstruction_common
!! \ingroup reconstruction
!----------------------------------------------------------------------------!
MODULE reconstruction_linear_mod
  USE reconstruction_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE marray_base_mod
  USE marray_compound_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, EXTENDS (reconstruction_base)  :: reconstruction_linear
  CLASS(marray_base), ALLOCATABLE :: slopes,dx
  CONTAINS
    PROCEDURE :: InitReconstruction_linear
    PROCEDURE :: CalculateStates
    PROCEDURE :: CalculateSlopes
    PROCEDURE :: Finalize
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
    CHARACTER(LEN=60)                           :: key
    INTEGER                                     :: err,limiter,valwrite,l,m
    REAL                                        :: theta
    !------------------------------------------------------------------------!
    ! allocate memory for all arrays used in reconstruction_linear
    ALLOCATE(this%slopes,this%dx, &
         STAT = err)
    IF (err.NE.0) THEN
       CALL this%Error("InitReconstruction_linear",  "Unable to allocate memory.")
    END IF


    ! initialize parent
    CALL this%InitReconstruction(Mesh,Physics,config,IO,LINEAR,recontype_name)

    ! set defaults for the limiter
    CALL GetAttr(config, "limiter", limiter, MINMOD)
    CALL GetAttr(config, "theta", theta, 1.0)

    ! limiter settings
    SELECT CASE(limiter)
    CASE(MINMOD,MONOCENT,SWEBY,SUPERBEE,OSPRE,PP,VANLEER,NOLIMIT)
       CALL this%InitLogging(limiter,limitertype_name(limiter))
       SELECT CASE(limiter)
       CASE(SWEBY,MONOCENT)
         this%limiter_param = theta
       END SELECT
    CASE DEFAULT
       CALL this%Error("InitReconstruction_linear", "Unknown limiter")
    END SELECT

    ! check parameter settings for PP limiter
    IF (limiter.EQ.PP) THEN
      CALL this%Error("InitReconstruction_linear","PP limiter currently not supported.")
!       IF (this%limiter_param.LT.0.0) &
!         CALL this%Error("InitReconstruction_linear", &
!                "Parameter of PP limiter must be greater than 0")
!       IF (this%limiter_param.GT.EPSILON(this%limiter_param)) &
!         CALL this%Warning("InitReconstruction_linear", &
!                "Parameter of PP limiter should be less than machine precision.")
    END IF

    ! create new rank 2 mesh array for slopes
    this%slopes = marray_base(Mesh%NDIMS,Physics%VNUM)

    ! create new rank 1 mesh array for reconstruction points
    this%dx = marray_base(Mesh%NFACES)

    ! zero the slopes
    this%slopes%data1d(:) = 0.0

    ! initialize coordinate differences for reconstruction
    m = 1
    IF (Mesh%INUM.GT.1) THEN
      ! reconstruct along x-direction
      this%dx%data4d(:,:,:,m) = Mesh%curv%faces(:,:,:,1,1) - Mesh%curv%bcenter(:,:,:,1)
      this%dx%data4d(:,:,:,m+1) = Mesh%curv%faces(:,:,:,2,1) - Mesh%curv%bcenter(:,:,:,1)
      m = m + 2
    END IF
    IF (Mesh%JNUM.GT.1) THEN
      this%dx%data4d(:,:,:,m) = Mesh%curv%faces(:,:,:,3,2) - Mesh%curv%bcenter(:,:,:,2)
      this%dx%data4d(:,:,:,m+1) = Mesh%curv%faces(:,:,:,4,2) - Mesh%curv%bcenter(:,:,:,2)
      m = m + 2
    END IF
    IF (Mesh%KNUM.GT.1) THEN
      this%dx%data4d(:,:,:,m) = Mesh%curv%faces(:,:,:,5,3) - Mesh%curv%bcenter(:,:,:,3)
      this%dx%data4d(:,:,:,m+1) = Mesh%curv%faces(:,:,:,6,3) - Mesh%curv%bcenter(:,:,:,3)
    END IF

     CALL GetAttr(config, "output/slopes", valwrite, 0)
     IF(valwrite.EQ.1) THEN
       m = 1
       IF (Mesh%INUM.GT.1) THEN
!NEC$ NOVECTOR
         DO l=1,Physics%VNUM
           key = TRIM(Physics%pvarname(l)) // "_xslope"
           CALL SetAttr(IO, TRIM(key), &
             this%slopes%data5d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m,l))
         END DO
         m = m + 1
       END IF
       IF (Mesh%JNUM.GT.1) THEN
!NEC$ NOVECTOR
         DO l=1,Physics%VNUM
           key = TRIM(Physics%pvarname(l)) // "_yslope"
           CALL SetAttr(IO, TRIM(key), &
             this%slopes%data5d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m,l))
         END DO
         m = m + 1
       END IF
       IF (Mesh%KNUM.GT.1) THEN
!NEC$ NOVECTOR
         DO l=1,Physics%VNUM
           key = TRIM(Physics%pvarname(l)) // "_zslope"
           CALL SetAttr(IO, TRIM(key), &
             this%slopes%data5d(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,m,l))
         END DO
       END IF
     END IF

    CALL this%Info("            limiter:           " // TRIM(this%GetName()))
  END SUBROUTINE InitReconstruction_linear


  !> \public Reconstructes states at cell boundaries
  PURE SUBROUTINE CalculateStates(this,Mesh,Physics,rvar,rstates)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_linear), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_base),          INTENT(IN)    :: Physics
    CLASS(marray_compound),       INTENT(INOUT) :: rvar
    CLASS(marray_compound),       INTENT(INOUT) :: rstates
    !------------------------------------------------------------------------!
    INTEGER                                     :: l,n
    !------------------------------------------------------------------------!
    ! calculate slopes first
    CALL this%CalculateSlopes(Mesh,Physics,rvar)

    ! reconstruct cell face values
!NEC$ SHORTLOOP
    DO l=1,Physics%VNUM
!NEC$ SHORTLOOP
      DO n=1,Mesh%NFACES
        rstates%data3d(:,n,l) = rvar%data2d(:,l) + &
           this%slopes%data3d(:,((n+1)/2),l)*this%dx%data2d(:,n)
      END DO
    END DO
  END SUBROUTINE CalculateStates

  !> \public computes the limited slopes
  PURE SUBROUTINE CalculateSlopes(this,Mesh,Physics,rvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_linear), INTENT(INOUT) :: this
    CLASS(mesh_base),             INTENT(IN)    :: Mesh
    CLASS(physics_base),          INTENT(IN)    :: Physics
    CLASS(marray_compound),       INTENT(INOUT) :: rvar
    !------------------------------------------------------------------------!
    INTEGER                                     :: i,j,k,l,m
    !------------------------------------------------------------------------!
    m = 1 ! this counts the dimensions used for computation (last index in slopes)
    ! choose limiter & check for cases where dimension needs to be neclected
    SELECT CASE(this%GetType())
    CASE(MINMOD)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * minmod2_limiter(&
                   rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l), &
                   rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
!NEC$ NOVECTOR
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ NOVECTOR
            DO j=i+1,SIZE(rvar%data3d,1)-i
                this%slopes%data4d(j,k,m,l) = Mesh%invdy * minmod2_limiter(&
                  rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l), &
                  rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l))
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
              this%slopes%data3d(j,m,l) = Mesh%invdz * minmod2_limiter(&
                rvar%data2d(j,l) - rvar%data2d(j-i,l), &
                rvar%data2d(j+i,l) - rvar%data2d(j,l))
          END DO
        END DO
      END IF
    CASE(MONOCENT)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * minmod3_limiter(&
                   this%limiter_param*(rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l)),&
                   this%limiter_param*(rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l)),&
                   0.5*(rvar%data4d(i+1,j,k,l) - rvar%data4d(i-1,j,k,l)))
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ IVDEP
            DO j=i+1,SIZE(rvar%data3d,1)-i
              this%slopes%data4d(j,k,m,l) = Mesh%invdy * minmod3_limiter(&
                  this%limiter_param*(rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l)),&
                  this%limiter_param*(rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l)),&
                  0.5*(rvar%data3d(j+i,k,l) - rvar%data3d(j-i,k,l)))
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
            this%slopes%data3d(j,m,l) = Mesh%invdz * minmod3_limiter(&
                this%limiter_param*(rvar%data2d(j,l) - rvar%data2d(j-i,l)),&
                this%limiter_param*(rvar%data2d(j+i,l) - rvar%data2d(j,l)),&
                0.5*(rvar%data2d(j+i,l) - rvar%data2d(j-i,l)))
          END DO
        END DO
      END IF
    CASE(SWEBY)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * sweby_limiter(&
                   rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l), &
                   rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l),&
                   this%limiter_param)
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ IVDEP
            DO j=i+1,SIZE(rvar%data3d,1)-i
              this%slopes%data4d(j,k,m,l) = Mesh%invdy * sweby_limiter(&
                  rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l),&
                  rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l),&
                  this%limiter_param)
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
            this%slopes%data3d(j,m,l) = Mesh%invdz * sweby_limiter(&
                rvar%data2d(j,l) - rvar%data2d(j-i,l),&
                rvar%data2d(j+i,l) - rvar%data2d(j,l),&
                this%limiter_param)
          END DO
        END DO
      END IF
    CASE(SUPERBEE)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * sweby_limiter(&
                   rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l), &
                   rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l), 2.0)
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ IVDEP
            DO j=i+1,SIZE(rvar%data3d,1)-i
              this%slopes%data4d(j,k,m,l) = Mesh%invdy * sweby_limiter(&
                  rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l),&
                  rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l), 2.0)
            END DO
          END DO
        END DO
         m = m + 1
     END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
            this%slopes%data3d(j,m,l) = Mesh%invdz * sweby_limiter(&
                rvar%data2d(j,l) - rvar%data2d(j-i,l),&
                rvar%data2d(j+i,l) - rvar%data2d(j,l), 2.0)
          END DO
        END DO
      END IF
    CASE(OSPRE)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * ospre_limiter(&
                   rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l), &
                   rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ IVDEP
            DO j=i+1,SIZE(rvar%data3d,1)-i
              this%slopes%data4d(j,k,m,l) = Mesh%invdy * ospre_limiter(&
                  rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l),&
                  rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l))
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
            this%slopes%data3d(j,m,l) = Mesh%invdz * ospre_limiter(&
                rvar%data2d(j,l) - rvar%data2d(j-i,l),&
                rvar%data2d(j+i,l) - rvar%data2d(j,l))
          END DO
        END DO
      END IF
!      CASE(PP)
!        ! calculate slopes in both-directions
! !NEC$ SHORTLOOP
!         DO l=1,Physics%VNUM
!          DO k=Mesh%KGMIN,Mesh%KGMAX
!            DO j=Mesh%JGMIN+Mesh%JP1,Mesh%JGMAX+Mesh%JM1
!!NEC$ IVDEP
!              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
!                 CALL pp_limiter(this%slopes%data5d(i,j,k,m,l),&
!                                 this%slopes%data5d(i,j,k,m,l),&
!                    rvar%data4d(i-1,j+1,k)-rvar%data4d(i,j,k),&
!                    rvar%data4d(i  ,j+1,k)-rvar%data4d(i,j,k),&
!                    rvar%data4d(i+1,j+1,k)-rvar%data4d(i,j,k),&
!                    rvar%data4d(i-1,j  ,k)-rvar%data4d(i,j,k),&
!                    this%limiter_param,&
!                    rvar%data4d(i+1,j  ,k)-rvar%data4d(i,j,k),&
!                    rvar%data4d(i-1,j-1,k)-rvar%data4d(i,j,k),&
!                    rvar%data4d(i  ,j-1,k)-rvar%data4d(i,j,k),&
!                    rvar%data4d(i+1,j-1,k)-rvar%data4d(i,j,k))
!              END DO
!            END DO
!          END DO
!        END DO
    CASE(VANLEER)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * vanleer_limiter(&
                   rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l), &
                   rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
       ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ IVDEP
            DO j=i+1,SIZE(rvar%data3d,1)-i
              this%slopes%data4d(j,k,m,l) = Mesh%invdy * vanleer_limiter(&
                  rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l),&
                  rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l))
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
            this%slopes%data3d(j,m,l) = Mesh%invdz * vanleer_limiter(&
                rvar%data2d(j,l) - rvar%data2d(j-i,l),&
                rvar%data2d(j+i,l) - rvar%data2d(j,l))
          END DO
        END DO
      END IF
    CASE(NOLIMIT)
      ! calculate slopes in x-direction
      IF (Mesh%INUM.GT.1) THEN
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
            DO j=Mesh%JGMIN,Mesh%JGMAX
!NEC$ IVDEP
              DO i=Mesh%IGMIN+Mesh%IP1,Mesh%IGMAX+Mesh%IM1
                this%slopes%data5d(i,j,k,m,l) = Mesh%invdx * nolimit_limiter(&
                   rvar%data4d(i,j,k,l) - rvar%data4d(i-1,j,k,l), &
                   rvar%data4d(i+1,j,k,l) - rvar%data4d(i,j,k,l))
              END DO
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in y-direction
      IF (Mesh%JNUM.GT.1) THEN
        i = Mesh%IGMAX-Mesh%IGMIN+1
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
          DO k=Mesh%KGMIN,Mesh%KGMAX
! use collapsed arrays for better vectorization
!NEC$ IVDEP
            DO j=i+1,SIZE(rvar%data3d,1)-i
              this%slopes%data4d(j,k,m,l) = Mesh%invdy * nolimit_limiter(&
                  rvar%data3d(j,k,l) - rvar%data3d(j-i,k,l),&
                  rvar%data3d(j+i,k,l) - rvar%data3d(j,k,l))
            END DO
          END DO
        END DO
        m = m + 1
      END IF
      ! calculate slopes in z-direction
      IF (Mesh%KNUM.GT.1) THEN
        i = (Mesh%IGMAX-Mesh%IGMIN+1)*(Mesh%JGMAX-Mesh%JGMIN+1)
!NEC$ SHORTLOOP
        DO l=1,Physics%VNUM
! use collapsed arrays for better vectorization
!NEC$ IVDEP
          DO j=i+1,SIZE(rvar%data2d,1)-i
            this%slopes%data3d(j,m,l) = Mesh%invdz * nolimit_limiter(&
                rvar%data2d(j,l) - rvar%data2d(j-i,l),&
                rvar%data2d(j+i,l) - rvar%data2d(j,l))
          END DO
        END DO
      END IF
    END SELECT
  END SUBROUTINE CalculateSlopes

  SUBROUTINE Finalize(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_linear), INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    CALL this%slopes%Destroy()
    CALL this%dx%Destroy()
    DEALLOCATE(this%slopes,this%dx)

    CALL this%Finalize_base()
  END SUBROUTINE Finalize

!----------------------------------------------------------------------------!
!> \par elemental non-class subroutines / functions

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

END MODULE reconstruction_linear_mod
