!
!     domain mesh module
!     (c) Jure Ravnik, 2025
!
MODULE mDomainMesh
      
    USE mCommon
    IMPLICIT NONE
  
!
!     Types
!
!
! ----------------------------------------------------------------------
!
    TYPE DMElementType
        INTEGER nno  ! number of nodes in element
        INTEGER type ! element type
        INTEGER, POINTER :: con(:)  ! id's of nodes
    END TYPE

    TYPE DMNodeType
        REAL(rk) :: x(3)  ! node coordiantes, x,y,z
    END TYPE    
!
!     Variables
!
    INTEGER DMnnodes
    TYPE(DMNodeType), POINTER :: DMnode(:)

    INTEGER DMnelem
    TYPE (DMElementType), POINTER :: DMelement(:)
 
    CHARACTER(255) DMmeshName
  
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
    SUBROUTINE DMinitElement(n)

        INTEGER i,n

        DMnelem = n

        ALLOCATE ( DMelement(DMnelem) )

        DO i=1,DMnelem
          DMelement(i)%type = -1
          DMelement(i)%nno  = -1
        END DO

    END SUBROUTINE
  
!
! ----------------------------------------------------------------------
!

!
! ----------------------------------------------------------------------
!
    SUBROUTINE DMinitNodes(n)

        INTEGER i,j,n
  
        DMnnodes = n
  
        ALLOCATE ( DMnode(DMnnodes) )
  
        DO i=1,DMnnodes
          DO j=1,3
            DMnode(i)%x(j)=i*j
          END DO
        END DO
  
    END SUBROUTINE


END MODULE