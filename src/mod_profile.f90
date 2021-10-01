!
!     PostProcessing module
!
      MODULE mProfile

      USE mCommon
      IMPLICIT NONE
!
!     Types
!
!
! ----------------------------------------------------------------------
      TYPE lineProfile
        
        REAL(rk)                        :: firstLine(2,3)
        REAL(rk)                        :: secndLine(2,3)
        REAL(rk), POINTER               :: ptsOnLine(:,:)
        REAL(rk), POINTER               :: ptsOnProj(:,:)
        REAL(rk), POINTER               :: ptsOnBody(:,:)
        REAL(rk), POINTER               :: results(:,:)
        INTEGER                     :: numOfDiv
        CHARACTER(255)              :: name
        CHARACTER(255),POINTER      :: WallName(:)
        INTEGER                     :: numWalls
      END TYPE
! ----------------------------------------------------------------------
!

!
!     Variables
!
      INTEGER noLines ! number of profiles
      TYPE(lineProfile), POINTER :: prof(:)

!
! ----------------------------------------------------------------------
!

! ----------------------------------------------------------------------
!

!
!     Subroutines
!
      CONTAINS
!
! ----------------------------------------------------------------------
!
      SUBROUTINE prof_init(nOfProf)

      INTEGER     ::  nOfProf

      noLines = nOfProf
      ALLOCATE (prof(noLines))

      END SUBROUTINE

!
! ----------------------------------------------------------------------
!
      SUBROUTINE lines_init(numNeq)
  
        INTEGER i, numNeq
        DO i = 1, noLines
            ALLOCATE (prof(i)%ptsOnLine(prof(i)%numOfDiv + 1, 3))
            ALLOCATE (prof(i)%ptsOnProj(prof(i)%numOfDiv + 1, 3))
            ALLOCATE (prof(i)%ptsOnBody(prof(i)%numOfDiv + 1, 3))
            ALLOCATE (prof(i)%results( numNeq + 1, prof(i)%numOfDiv + 1 ))
        END DO
  
      END SUBROUTINE
! ----------------------------------------------------------------------
!
      SUBROUTINE walls_init(wN,which)

      INTEGER     :: which, wN

      prof(which)%numWalls = wN
      ALLOCATE (prof(which)%WallName(wN))

      END SUBROUTINE

      END MODULE