!
!
!
!     domain data module
!     (c) Jure Ravnik, 2025
!
      MODULE mDomainData
      
        USE mCommon
        IMPLICIT NONE
  
  !
  !     Types
  !
  
  !
  ! ----------------------------------------------------------------------
  !
        TYPE DataType
  
          CHARACTER(255) name
          INTEGER n
          REAL(rk), POINTER :: val(:)
  
        END TYPE
  
  
  
!
! ----------------------------------------------------------------------
!
  
!
!     Variables
!
    INTEGER ddnN ! number of nodal data sets
    INTEGER ddnE ! number of element data sets
    TYPE(DataType), POINTER :: ddN(:),ddE(:)

  
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
    SUBROUTINE dd_init(nnodes,nelem)
  
        INTEGER i,nnodes,nelem
  
        ddnN = 1
        ddnE = 4
      
        allocate( ddN(ddnN) )
        allocate( ddE(ddnE) )

        ddN(1)%name = "node"
        do i = 1,ddnN
            ddN(i)%n = nnodes
            allocate(ddN(i)%val(ddN(i)%n))
        end do

        do i =1,ddN(1)%n
            ddN(1)%val(i) = i
        end do

        ddE(1)%name = "u"
        ddE(2)%name = "v"
        ddE(3)%name = "w"
        ddE(4)%name = "p"
        do i = 1,ddnE
            ddE(i)%n = nelem
            allocate(ddE(i)%val(ddE(i)%n))
            ddE(i)%val = 0.0_rk
        end do

    END SUBROUTINE
  
  
  
  !
  ! ----------------------------------------------------------------------
  !
  
        END MODULE
  
  