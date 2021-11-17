!
!     mesh module
!     (c) Jure Ravnik, 2017
!
      MODULE mMesh
      
      USE mCommon
      IMPLICIT NONE

!
!     Types
!
!
! ----------------------------------------------------------------------
!
      TYPE ElementType

        INTEGER type ! element type
        INTEGER nno  ! number of nodes in element
        INTEGER bcid ! which wall the element belongs to
        INTEGER, POINTER :: con(:)  ! id's of nodes
        REAL(rk) area ! area of element
        REAL(rk) normal(3) ! normal vector

      END TYPE
!
! ----------------------------------------------------------------------
!

      TYPE NodeType

        REAL(rk) :: x(3)  ! node coordiantes, x,y,z
        INTEGER bcid  ! wall to which the node belongs

      END TYPE
!
! ----------------------------------------------------------------------
!
      TYPE WallType

        CHARACTER(100) name ! wall name
        INTEGER id
        REAL(rk) area

      END TYPE

!
! ----------------------------------------------------------------------
!
      TYPE SubdomainType

        CHARACTER(100) name ! wall name
        REAL(rk) diff ! diffusivity of the material in the subdomain
        INTEGER nofw
        INTEGER, POINTER :: loWalls(:) ! list of walls in subdomain (nofw)
        REAL(rk), POINTER :: normMul(:) ! normal multiplyer for wall (nofw)
        REAL(rk) fspCerr,QspCerr


        INTEGER nnodes,nqnodes,nelem,nsp
        INTEGER, POINTER :: elementList(:) ! (nelem)
        INTEGER, POINTER :: nodeList(:) ! (nnodes)
        INTEGER, POINTER :: BCidList(:) !  (nnodes) ! wall to which node belongs
        INTEGER, POINTER :: qnodeList(:) ! (nqnodes)
        INTEGER, POINTER :: qBCidList(:) ! (nqnodes)

        ! Laplace matrices
        REAL(rk), POINTER :: Hmat(:,:),Gmat(:,:)

        ! Stokes matrices
        REAL(rk), POINTER :: TxxMat(:,:), GxxMat(:,:)
        REAL(rk), POINTER :: TyyMat(:,:), GyyMat(:,:)
        REAL(rk), POINTER :: TzzMat(:,:), GzzMat(:,:)
        REAL(rk), POINTER :: TxyMat(:,:), GxyMat(:,:)
        REAL(rk), POINTER :: TxzMat(:,:), GxzMat(:,:)
        REAL(rk), POINTER :: TyzMat(:,:), GyzMat(:,:)
        REAL(rk), POINTER :: prTxMat(:,:),prTyMat(:,:),prTzMat(:,:)
        REAL(rk), POINTER :: prGxMat(:,:),prGyMat(:,:),prGzMat(:,:)

      END TYPE



!
! ----------------------------------------------------------------------
!

!
!     Variables
!
      INTEGER nnodes
      TYPE (NodeType), POINTER :: node(:)

      INTEGER nqnodes
      TYPE (NodeType), POINTER :: qnode(:)

      INTEGER nelem
      TYPE (ElementType), POINTER :: element(:)

      INTEGER nofw
      TYPE (WallType), POINTER :: wall(:)

      INTEGER nosd ! number of subdomains
      TYPE(SubdomainType), POINTER :: subdomain(:)

      CHARACTER(255) meshName

!     Constants to make coding easier
      INTEGER iFunction,iFlux,iContact,iConst,iLin,iQuad
      PARAMETER (iFunction=1,iFlux=2,iContact=3)
      PARAMETER (iConst=1, iLin=2, iQuad=3)

!     Mesh extents
      REAL(rk) mextMin(3),mextMax(3)

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
      SUBROUTINE initWall(n)

      INTEGER i,n

      nofw = n

      ALLOCATE ( wall(nofw) )

      DO i=1,nofw
        wall(i)%name = ""
        wall(i)%id  = -1
      END DO

      END SUBROUTINE


!
! ----------------------------------------------------------------------
!
      SUBROUTINE initElement(n)

      INTEGER i,n

      nelem = n

      ALLOCATE ( element(nelem) )

      DO i=1,nelem
        element(i)%type = -1
        element(i)%nno  = -1
        element(i)%area  = 0.0_rk
      END DO

      END SUBROUTINE


!
! ----------------------------------------------------------------------
!
      SUBROUTINE PrintElements()

      INTEGER i,j

      DO i=1,nelem
        print *,i,element(i)%type,(element(i)%con(j),j=1,element(i)%nno)
      END DO

      END SUBROUTINE

!
! ----------------------------------------------------------------------
!
      SUBROUTINE initNodes(n)

      INTEGER i,j,n

      nnodes = n

      ALLOCATE ( node(nnodes) )

      DO i=1,nnodes
        DO j=1,3
          node(i)%x(j)=i*j
        END DO
      END DO

      END SUBROUTINE

!
! ----------------------------------------------------------------------
!
      SUBROUTINE PrintNodes()

      INTEGER i

      DO i=1,nnodes
        print *,i,node(i)%x(1),node(i)%x(2),node(i)%x(3)
      END DO

      END SUBROUTINE

!
! ----------------------------------------------------------------------
!

      END MODULE





