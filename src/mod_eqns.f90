!
!
!
!     equations module
!     (c) Jure Ravnik, 2017
!
      MODULE mEqns
      
      USE mCommon
      USE plis
      IMPLICIT NONE

!
!     Types
!
!
! ----------------------------------------------------------------------
!
      TYPE ConditionType

!       Boundary condition
        INTEGER known ! is function of flux known at the boundary
        INTEGER type ! const,lin,quad distribution at the boundary
        REAL(rk) params(10) ! boundary condition parameters

      END TYPE
!
! ----------------------------------------------------------------------
!
      TYPE SolverType

        INTEGER type ! solver type
        INTEGER pret ! preconditioner type
        INTEGER prep ! preconditioner placement
        INTEGER maxit! maximal number of iterations
        INTEGER stopt! stop type
        REAL(rk) eps  ! stop epsilon
!
! ----------------------------------------------------------------------
!
      END TYPE

      TYPE CRSmatrixType
        INTEGER    :: neq
        INTEGER(8) :: nnz
        INTEGER(8), POINTER :: i(:)
        INTEGER, POINTER :: j(:)
        REAL(rk), POINTER :: v(:)
      END TYPE CRSmatrixType


      TYPE DivideNodesType
            INTEGER, POINTER :: rowTG(:)
            INTEGER, POINTER :: rowRHS(:)
            INTEGER, POINTER :: rowSYS(:)
            INTEGER, POINTER :: UQ(:)
            INTEGER, POINTER :: een(:)
            INTEGER, POINTER :: irow(:)
      END TYPE

!
! ----------------------------------------------------------------------
!
      TYPE EquationType

        CHARACTER(255) name
        REAL(rk), POINTER :: u(:),q(:),p(:)
        TYPE(ConditionType), POINTER :: boundary(:)
        TYPE(ConditionType), POINTER :: initial(:)

        TYPE(SolverType) :: slv  ! definition of SLE solver

        INTEGER nx  ! number of unknowns (=cols in sysM)
        INTEGER nb  ! number of rows in rhs vector (=cols in rhs matrix)
        INTEGER neq ! number of equations (=rows in sys matrix)
        REAL(rk), POINTER :: x(:),b(:),pivot(:)

        TYPE(CRSmatrixType) :: sysMcrs ! system matrix
        TYPE(CRSmatrixType) :: rhsMcrs ! right hand side matrix

        INTEGER, POINTER :: col(:) ! ! position in x or b vector
        INTEGER, POINTER :: qcol(:) ! ! position in x or b vector

        TYPE(systemLinEq) sle ! system to be solved with LIS library
        TYPE(DivideNodesType) dn ! nodes divided between processes

      END TYPE



!
! ----------------------------------------------------------------------
!

!
!     Variables
!
      INTEGER neq ! number of equations
      TYPE(EquationType), POINTER :: eqn(:)
      TYPE(EquationType) stk

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
      SUBROUTINE eqn_init(nofw,nnodes,nqnodes)

      INTEGER i,j,nofw,nnodes,nqnodes


      ALLOCATE (eqn(neq))
      DO i=1,neq
        ALLOCATE (eqn(i)%boundary(nofw))
        ALLOCATE (eqn(i)%initial(nofw))
        ALLOCATE (eqn(i)%u(nnodes))
        ALLOCATE (eqn(i)%q(nqnodes))
        ALLOCATE (eqn(i)%p(nqnodes))

        DO j=1,nqnodes
            eqn(i)%p(j)=0.0_rk
        end do

      END DO
      END SUBROUTINE



!
! ----------------------------------------------------------------------
!

      END MODULE

