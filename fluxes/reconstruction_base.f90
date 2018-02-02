!############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: reconstruction_base.f90                                           #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
!> \author Tobias Illenseer
!> \author Jannes Klee
!!
!! \brief base module for reconstruction process
!!
!! \ingroup reconstruction
!----------------------------------------------------------------------------!
MODULE reconstruction_base_mod
  USE logging_base_mod
  USE mesh_base_mod
  USE physics_base_mod
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  TYPE, ABSTRACT, EXTENDS (logging_base)   ::  reconstruction_base
     !> \name Variables
     CLASS(logging_base), ALLOCATABLE      :: limiter                 !< limiter (linear, const)
     LOGICAL                               :: primcons                !< true if primitive
     REAL                                  :: limiter_param           !< limiter parameter
     REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: xslopes,yslopes,zslopes !< limited slopes
  CONTAINS
    PROCEDURE                              :: InitReconstruction
    PROCEDURE (CalculateStates), DEFERRED  :: CalculateStates
    PROCEDURE                              :: PrimRecon
    PROCEDURE                              :: FinalizeReconstruction
  END TYPE reconstruction_base

  ABSTRACT INTERFACE
    PURE SUBROUTINE CalculateStates(this,Mesh,Physics,npos,dx,dy,dz,rvar,rstates)
      IMPORT reconstruction_base, mesh_base, physics_base
      IMPLICIT NONE
      CLASS(reconstruction_base), INTENT(INOUT)   :: this
      CLASS(mesh_base),           INTENT(IN)      :: Mesh
      CLASS(physics_base),        INTENT(IN)      :: Physics
      INTEGER                                     :: npos
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                      Mesh%KGMIN:Mesh%KGMAX,npos) :: dx,dy,dz
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                      Mesh%KGMIN:Mesh%KGMAX,Physics%VNUM) &
                                                  :: rvar
      REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                      Mesh%KGMIN:Mesh%KGMAX,npos,Physics%VNUM) &
                                                  :: rstates
      INTENT(IN)                                  :: npos,dx,dy,dz,rvar
      INTENT(OUT)                                 :: rstates
    END SUBROUTINE
  END INTERFACE

  !--------------------------------------------------------------------------!

  !\todo{here something is still not right! (constant, linear, primitive, conservative all at this place?)}
  !> \name Public Attributes
  INTEGER, PARAMETER :: PRIMITIVE    = 1           ! True  ! not type LOGICAL!
  INTEGER, PARAMETER :: CONSERVATIVE = 0           ! False ! because it cant !

  INTEGER, PARAMETER :: CONSTANT     = 1
  INTEGER, PARAMETER :: LINEAR       = 2
  !--------------------------------------------------------------------------!
  PUBLIC ::                 &
       ! types
       reconstruction_base, &
       ! constants
       CONSTANT, LINEAR,    &
       PRIMITIVE, CONSERVATIVE
  !--------------------------------------------------------------------------!

CONTAINS


  !> \public Constructor of base reconstruction module
  SUBROUTINE InitReconstruction(this,Mesh,Physics,config,IO,rtype,rname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_base), INTENT(INOUT) :: this
    CLASS(mesh_base),           INTENT(IN)    :: Mesh
    CLASS(physics_base),        INTENT(IN)    :: Physics
    TYPE(Dict_TYP),             POINTER       :: config,IO
    INTEGER                                   :: rtype
    CHARACTER(LEN=32)                         :: rname
    !------------------------------------------------------------------------!
    INTEGER                                   :: variables
    CHARACTER(LEN=32)                         :: infostr
    CHARACTER(LEN=60)                         :: key
    INTEGER                                   :: valwrite,i
    INTEGER                                   :: order
    !------------------------------------------------------------------------!
    INTENT(IN)                                :: config,rtype,rname
    !------------------------------------------------------------------------!
    CALL this%InitLogging(rtype,rname)

    ! check initialization of Mesh and Physics
    IF (.NOT.Mesh%Initialized().OR..NOT.Physics%Initialized()) & ! TODO: Why should this not be seperated?
         CALL this%Error("InitFluxes","mesh and/or physics module uninitialized")

    ! set general reconstruction defaults
    CALL GetAttr(config, "order", order, LINEAR)

    CALL GetAttr(config, "variables", variables, CONSERVATIVE)

    IF(variables.EQ.0) THEN
        this%primcons = .FALSE.
    ELSE
        this%primcons = .TRUE.
    END IF

    ! print some information
    CALL this%Info(" RECONSTR-> order:             " // TRIM(this%GetName()))
    IF (PrimRecon(this)) THEN
       WRITE (infostr,'(A)') "primitive"
    ELSE
       WRITE (infostr,'(A)') "conservative"
    END IF
    CALL this%Info("            variables:         " // TRIM(infostr))

    CALL GetAttr(config, "output/slopes", valwrite, 0)
    IF(valwrite.EQ.1) THEN
      DO i=1, Physics%VNUM
        key = TRIM(Physics%pvarname(i)) // "_xslope"
        CALL SetAttr(IO, TRIM(key), &
          this%xslopes(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
        key = TRIM(Physics%pvarname(i)) // "_yslope"
        CALL SetAttr(IO, TRIM(key), &
          this%yslopes(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
        key = TRIM(Physics%pvarname(i)) // "_zslope"
        CALL SetAttr(IO, TRIM(key), &
          this%zslopes(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,Mesh%KMIN:Mesh%KMAX,i))
      END DO
    END IF
  END SUBROUTINE InitReconstruction


  PURE FUNCTION PrimRecon(this) RESULT(pc)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_base), INTENT(IN) :: this
    LOGICAL                                :: pc
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    pc = this%primcons
  END FUNCTION PrimRecon


  SUBROUTINE FinalizeReconstruction(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(reconstruction_base), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    IF (.NOT.this%Initialized()) &
         CALL this%Error("CloseReconstruction","not initialized")
  END SUBROUTINE FinalizeReconstruction

END MODULE reconstruction_base_mod
