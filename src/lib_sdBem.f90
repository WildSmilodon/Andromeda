!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdSetUpXandB()

      USE mMesh
      USE mEqns
      USE mPar

      IMPLICIT NONE

      INTEGER en,isp,mywall,isd,i
      INTEGER notSet

      PARAMETER (notSet = -999)

      DO en=1,neq

        eqn(en)%nx=0  ! number of unknowns (=cols in sysM)
        eqn(en)%nb=0  ! number of rows in rhs vector  (=cols in rhs matrix)
        eqn(en)%neq=0 ! number of equations (=rows in sys matrix)

        ALLOCATE (eqn(en)%col(nnodes))
        ALLOCATE (eqn(en)%qcol(nqnodes))

        eqn(en)%col = notSet
        eqn(en)%qcol= notSet

!
!       Loop over subdomains
!
        DO isd=1,nosd
!
!         function source points
!
          DO i = 1,subdomain(isd)%nnodes
            isp = subdomain(isd)%nodeList(i)  ! source point
            myWall = subdomain(isd)%BCidList(i)

            IF (eqn(en)%boundary(myWall)%known.EQ.iFunction) THEN
              eqn(en)%nb=eqn(en)%nb+1
              eqn(en)%col(isp)=eqn(en)%nb
            ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iFlux) THEN
              eqn(en)%nx=eqn(en)%nx+1      ! add an unknown
              eqn(en)%col(isp)=eqn(en)%nx
              eqn(en)%neq=eqn(en)%neq+1    ! add an equation

            ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iContact) THEN
              eqn(en)%neq=eqn(en)%neq+1    ! add an equation

              IF (eqn(en)%col(isp).EQ.notSet) THEN ! do not add unknown twice
                eqn(en)%nx=eqn(en)%nx+1      ! add an unknown
                eqn(en)%col(isp)=eqn(en)%nx
              END IF

            END IF

          END DO ! u source point
!
!         Flux source points
!
          DO i = 1,subdomain(isd)%nqnodes
            isp = subdomain(isd)%qnodeList(i)  ! source point
            myWall = subdomain(isd)%qBCidList(i)

            IF (eqn(en)%boundary(myWall)%known.EQ.iFunction) THEN
              eqn(en)%nx=eqn(en)%nx+1      ! add an unknown
              eqn(en)%qcol(isp)=eqn(en)%nx
              eqn(en)%neq=eqn(en)%neq+1    ! add an equation

            ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iFlux) THEN
              eqn(en)%nb=eqn(en)%nb+1
              eqn(en)%qcol(isp)=eqn(en)%nb
            ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iContact) THEN
              eqn(en)%neq=eqn(en)%neq+1    ! add an equation

              IF (eqn(en)%qcol(isp).EQ.notSet) THEN ! do not add unknown twice
                eqn(en)%nx=eqn(en)%nx+1      ! add an unknown
                eqn(en)%qcol(isp)=eqn(en)%nx
              END IF
            END IF

          END DO ! q source point

        END DO ! subdomains
!
!       Alocate vectors of unknowns and knowns
!
        ALLOCATE (eqn(en)%x(eqn(en)%nx))
        ALLOCATE (eqn(en)%Pivot(eqn(en)%nx))
        ALLOCATE (eqn(en)%b(eqn(en)%nb))
!
!       Report to log file
!
        CALL WriteToLog("")
        CALL WriteToLog("--------------------")
        WRITE (parLogTekst,'(A,A)') "Equation = ", TRIM(eqn(en)%name)
        CALL WriteToLog(parLogTekst)
        WRITE (parLogTekst,'(A,I0,A,I0,A)') "System matrix = ",eqn(en)%neq, " x ",eqn(en)%nx
        CALL WriteToLog(parLogTekst)
        WRITE (parLogTekst,'(A,I0,A,I0,A)') "RHS    matrix = ",eqn(en)%neq, " x ",eqn(en)%nb
        CALL WriteToLog(parLogTekst)


      END DO ! equations

      END SUBROUTINE


! ----------------------------------------------------------------------
!
      SUBROUTINE HandleCornersEdges(fname)
!
! ----------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon
      IMPLICIT NONE

      CHARACTER*(*) fname
      CHARACTER(255) OneLine,KeyWord
      CHARACTER(255) WallName1,WallName2,WallName3,TargetWallName
      CHARACTER(255) sdName      
      INTEGER lun
      INTEGER i,j,k,l,jj,ll,wid1,wid2,twid,wid3,isd,ie,je,m,jje
      INTEGER Cstring

      INTEGER, ALLOCATABLE :: tmp(:)

      ALLOCATE (tmp(100))


!
!     Read BiC file
!
      lun=11

      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')

        READ(Oneline,*) KeyWord
!
!       Node placement
!
        IF (Cstring(KeyWord,"EDGE").EQ.parYes) THEN


          READ(Oneline,*) KeyWord,sdName,WallName1,WallName2,TargetWallName

!         which subdomain
          isd=-1
          DO i=1,nosd
            IF (TRIM(sdName).EQ.TRIM(subdomain(i)%name)) THEN
              isd=i
            END IF
          END DO
          IF (isd.EQ.-1) CALL WriteToLog("BiC error! - EDGE subdomain name")

!         which wall 1
          wid1=-1
          DO i=1,nofw
            IF (TRIM(WallName1).EQ.TRIM(wall(i)%name)) THEN
              wid1=i
            END IF
          END DO
          IF (wid1.EQ.-1) CALL WriteToLog("BiC error! - EDGE wall name")

!         which wall 2
          wid2=-1
          DO i=1,nofw
            IF (TRIM(WallName2).EQ.TRIM(wall(i)%name)) THEN
              wid2=i
            END IF
          END DO
          IF (wid2.EQ.-1) CALL WriteToLog("BiC error! - EDGE wall name")

!         which wall TARGET
          twid=-1
          DO i=1,nofw
            IF (TRIM(TargetWallName).EQ.TRIM(wall(i)%name)) THEN
              twid=i
            END IF
          END DO
          IF (twid.EQ.-1) CALL WriteToLog("BiC error! - EDGE wall name")

!         apply edge rule to node%bcid

          DO i=1,subdomain(isd)%nelem
            ie = subdomain(isd)%elementList(i)
            IF (element(ie)%bcid.EQ.wid1) THEN
              DO j=1,subdomain(isd)%nelem
                je=subdomain(isd)%elementList(j)
                IF (element(je)%bcid.EQ.wid2) THEN
!                 Do elemnts share nodes?
                  DO k=1,element(ie)%nno
                    DO l=1,element(je)%nno
                      IF (element(ie)%con(k).EQ.element(je)%con(l)) THEN

                        DO m=1,subdomain(isd)%nnodes ! find node in my list
                          IF (subdomain(isd)%nodeList(m).EQ.element(ie)%con(k)) THEN
                            subdomain(isd)%BCidList(m)= twid
                          END IF
                        END DO
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END IF
          END DO

        END IF

        IF (Cstring(KeyWord,"CORNER").EQ.parYes) THEN


          READ(Oneline,*) KeyWord,sdName,WallName1,WallName2,WallName3,TargetWallName


!         which subdomain
          isd=-1
          DO i=1,nosd
            IF (TRIM(sdName).EQ.TRIM(subdomain(i)%name)) THEN
              isd=i
            END IF
          END DO
          IF (isd.EQ.-1) CALL WriteToLog("BiC error! - CORNER subdomain name")

!         which wall 1
          wid1=-1
          DO i=1,nofw
            IF (TRIM(WallName1).EQ.TRIM(wall(i)%name)) THEN
              wid1=i
            END IF
          END DO
          IF (wid1.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         which wall 2
          wid2=-1
          DO i=1,nofw
            IF (TRIM(WallName2).EQ.TRIM(wall(i)%name)) THEN
              wid2=i
            END IF
          END DO
          IF (wid2.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         which wall 3
          wid3=-1
          DO i=1,nofw
            IF (TRIM(WallName3).EQ.TRIM(wall(i)%name)) THEN
              wid3=i
            END IF
          END DO
          IF (wid3.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         which wall TARGET
          twid=-1
          DO i=1,nofw
            IF (TRIM(TargetWallName).EQ.TRIM(wall(i)%name)) THEN
              twid=i
            END IF
          END DO
          IF (twid.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         apply edge rule to BCidList

          DO i=1,subdomain(isd)%nelem
            ie = subdomain(isd)%elementList(i)
            IF (element(ie)%bcid.EQ.wid1) THEN
              DO j=1,nelem
                je = subdomain(isd)%elementList(j)
                IF (element(je)%bcid.EQ.wid2) THEN
                  DO jj=1,nelem
                    jje = subdomain(isd)%elementList(jj)
                    IF (element(jje)%bcid.EQ.wid3) THEN
!                     Do elemnts share nodes?
                      DO k=1,element(ie)%nno
                        DO l=1,element(je)%nno
                          DO ll=1,element(jje)%nno
                            IF (element(ie)%con(k).EQ.element(je)%con(l).AND.element(ie)%con(k).EQ.element(jje)%con(ll)) THEN

                              DO m=1,subdomain(isd)%nnodes ! find node in my list
                                IF (subdomain(isd)%nodeList(m).EQ.element(ie)%con(k)) THEN
                                  subdomain(isd)%BCidList(m)= twid
                                END IF
                              END DO
                            END IF
                          END DO
                        END DO
                      END DO
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END DO

        END IF

        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
      DEALLOCATE (tmp)
      RETURN

10    CALL WriteToLog("BiC error! - Can not open BiC file!")


      END SUBROUTINE


!
!     ------------------------------------------------------------------
!
      SUBROUTINE SetUpSubdomains()

      USE mMesh
      USE mEqns
      USE mPar

      IMPLICIT NONE

      INTEGER isd,iWall,jj,ie,i,j
      INTEGER, ALLOCATABLE :: tmp(:)


!
!     For each subdomain, get a list of elements
!
      DO isd=1,nosd ! loop over subdomains
        subdomain(isd)%nelem = 0
        DO jj=1,subdomain(isd)%nofw ! loop over wall in a subdomain
          iWall = subdomain(isd)%loWalls(jj) ! current wall number
          DO ie=1,nelem
            IF (element(ie)%bcid.EQ.iWall) THEN
              subdomain(isd)%nelem = subdomain(isd)%nelem + 1
            END IF
          END DO
        END DO

        ALLOCATE (subdomain(isd)%elementList(subdomain(isd)%nelem))

        subdomain(isd)%nelem = 0
        DO jj=1,subdomain(isd)%nofw ! loop over wall in a subdomain
          iWall = subdomain(isd)%loWalls(jj) ! current wall number
          DO ie=1,nelem
            IF (element(ie)%bcid.EQ.iWall) THEN
              subdomain(isd)%nelem = subdomain(isd)%nelem + 1
              subdomain(isd)%elementList(subdomain(isd)%nelem)=ie
            END IF
          END DO
        END DO

      END DO
!
!     For each subdomain, get a list of q-nodes
!
      DO isd=1,nosd ! loop over subdomains
        subdomain(isd)%nqnodes = subdomain(isd)%nelem
        ALLOCATE (subdomain(isd)%qnodeList(subdomain(isd)%nqnodes))
        ALLOCATE (subdomain(isd)%qbcidList(subdomain(isd)%nqnodes))
        DO ie=1,subdomain(isd)%nelem
          subdomain(isd)%qnodeList(ie)=subdomain(isd)%elementList(ie)
          subdomain(isd)%qbcidList(ie)=element(subdomain(isd)%elementList(ie))%bcid
        END DO
      END DO
!
!     For each subdomain, get a list of u-nodes
!
      ALLOCATE (tmp(nnodes))
      DO isd=1,nosd ! loop over subdomains
        tmp=0
        DO i=1,subdomain(isd)%nelem ! loop over elements in subdomain
          ie = subdomain(isd)%elementList(i)
          DO j=1,element(ie)%nno
            tmp(element(ie)%con(j))=element(ie)%bcid
          END DO
        END DO
        subdomain(isd)%nnodes = 0
        DO i=1,nnodes
          IF (tmp(i).NE.0) subdomain(isd)%nnodes = subdomain(isd)%nnodes + 1
        END DO
        ALLOCATE (subdomain(isd)%nodeList(subdomain(isd)%nnodes))
        ALLOCATE (subdomain(isd)%BCidList(subdomain(isd)%nnodes))
        subdomain(isd)%nnodes = 0
        DO i=1,nnodes
          IF (tmp(i).NE.0) THEN
            subdomain(isd)%nnodes = subdomain(isd)%nnodes + 1
            subdomain(isd)%nodeList(subdomain(isd)%nnodes) = i
            subdomain(isd)%BCidList(subdomain(isd)%nnodes) = tmp(i) ! wall ID for node
          END IF
        END DO
      END DO
      DEALLOCATE(tmp)


!
!     For each subdomain, get number of source points for integration
!
      DO isd=1,nosd ! loop over subdomains
        subdomain(isd)%nsp = subdomain(isd)%nnodes+subdomain(isd)%nqnodes ! number of source points (u+q)
      END DO
!
!       Report to log file
!
      CALL WriteToLog("")
      CALL WriteToLog("------- Subdomain list (name,nofw,nelem,nnodes,nqnodes) ----------")
      DO isd=1,nosd ! loop over subdomains
        WRITE (parLogTekst,'(A15,4(1X,I0))') TRIM(subdomain(isd)%name),subdomain(isd)%nofw,subdomain(isd)%nelem, &
                                       subdomain(isd)%nnodes,subdomain(isd)%nqnodes
        CALL WriteToLog(parLogTekst)
      END DO

      END




!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdIntRowU(isd,irow,rowH,rowG,err)

      USE mMesh
      USE Triangle ! provides tipw
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER isd ! subdomain number
      INTEGER irow ! source point
      INTEGER ie,i,j,IntRecDepth,isrc,jj,iWall
      REAL(rk) sx,sy,sz
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,c
      REAL(rk) integ,err,multi

      REAL(rk) rowH(nnodes),rowG(nqnodes)
      REAL(rk), ALLOCATABLE :: inteH(:)
!
!     Set integrals to zero
!
      rowH=0.0_rk
      rowG=0.0_rk

!
!     source point location
!
      sx=node(irow)%x(1)
      sy=node(irow)%x(2)
      sz=node(irow)%x(3)

!
!     Integrate over boudnary elements in subdomain
!
      DO jj=1,subdomain(isd)%nofw ! loop over walls in a subdomain
        iWall = subdomain(isd)%loWalls(jj) ! current wall number

        DO ie=1,nelem
          IF (element(ie)%bcid.EQ.iWall) THEN ! loop over elements in a wall
!
!           Find if element is singular
!
            isrc=0
            DO j=1,element(ie)%nno
              IF (irow.EQ.element(ie)%con(j)) isrc=j
            END DO
!
!           Perform intergration over a triangle
!           (due to recursive nature of integration, set to zero here)
            ALLOCATE (inteH(element(ie)%nno))
            inteG=0.0_rk
            DO j=1,element(ie)%nno
              inteH(j)=0.0_rk
            END DO
!
!           Element corners
!
            x1=node(element(ie)%con(1))%x(1)  !p%x(p%ibc(ie,1),1)
            y1=node(element(ie)%con(1))%x(2)  !p%x(p%ibc(ie,1),2)
            z1=node(element(ie)%con(1))%x(3)  !p%x(p%ibc(ie,1),3)

            x2=node(element(ie)%con(2))%x(1)  !p%x(p%ibc(ie,2),1)
            y2=node(element(ie)%con(2))%x(2)  !p%x(p%ibc(ie,2),2)
            z2=node(element(ie)%con(2))%x(3)  !p%x(p%ibc(ie,2),3)

            x3=node(element(ie)%con(3))%x(1)  !p%x(p%ibc(ie,3),1)
            y3=node(element(ie)%con(3))%x(2)  !p%x(p%ibc(ie,3),2)
            z3=node(element(ie)%con(3))%x(3)  !p%x(p%ibc(ie,3),3)

            IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
!
!             Set recursion depth for singular triangles
!
              IntRecDepth=parTriRecur

              CALL Triangle_BEMInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                              subdomain(isd)%normMul(jj)*element(ie)%normal(1), &
                              subdomain(isd)%normMul(jj)*element(ie)%normal(2), & 
                              subdomain(isd)%normMul(jj)*element(ie)%normal(3), &
                              element(ie)%area,integ,inteH(1),inteH(2),inteH(3),isrc,IntRecDepth)


            ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad

              x4=node(element(ie)%con(4))%x(1)
              y4=node(element(ie)%con(4))%x(2)
              z4=node(element(ie)%con(4))%x(3)

              CALL Quad_BEMInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sx,sy,sz,isrc, &
                            subdomain(isd)%normMul(jj),integ,inteh)

            ELSE
              CALL WriteToLog("Error :: Element type not supported!")
            END IF
!
!           Distribute results to integral matrices
!
            DO j=1,element(ie)%nno
              rowH(element(ie)%con(j))=rowH(element(ie)%con(j))+inteH(j)
            END DO
            rowG(ie)=integ

            DEALLOCATE (inteH)
          END IF
        END DO ! nbelem in wall
      END DO ! walls in subdomain
!
!     Find c parameter : u=1, q=(0,0,0)*(nx,ny,nz)=0
!
      c=0.0_rk
      DO i=1,nnodes
        c=c-rowH(i)
      END DO
!
!     Add c to diagonal of H matrix since c+Hu=Gq
!
      rowH(irow)=rowH(irow)+c
!
!       Check solution : u=x, q=(1,0,0)*(nx,ny,nz)=nx
!
      DO j=1,3
          c=0.0_rk
          DO i=1,nnodes
            c=c+rowH(i)*node(i)%x(j)
          END DO
          DO i=1,nelem
            multi=1.0_rk
            DO jj=1,subdomain(isd)%nofw ! find my wall in list of subdomain walls
              IF (element(i)%bcid.EQ.subdomain(isd)%loWalls(jj)) THEN
                multi = subdomain(isd)%normMul(jj)
              END IF
            END DO
            c=c-rowG(i)*element(i)%normal(j)*multi
          END DO
          err=MAX(err,ABS(c))
      END DO


      END SUBROUTINE

!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdIntRowQ(isd,irow,rowH,rowG,err)


      USE mMesh
      USE Triangle ! provides tipw
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER isd ! subdomain number
      INTEGER irow ! source point
      INTEGER ie,i,j,IntRecDepth,isrc,jj,iWall

      REAL(rk) sx,sy,sz
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,c
      REAL(rk) integ,err,shpf

      REAL(rk) rowH(nnodes),rowG(nqnodes)
      REAL(rk), ALLOCATABLE :: inteH(:)
!
!     Set integrals to zero
!
      rowH=0.0_rk
      rowG=0.0_rk

!
!     source point location
!
      sx=qnode(irow)%x(1)
      sy=qnode(irow)%x(2)
      sz=qnode(irow)%x(3)

!
!     Integrate over boudnary elements in subdomain
!
      DO jj=1,subdomain(isd)%nofw ! loop over walls in a subdomain
        iWall = subdomain(isd)%loWalls(jj) ! current wall number

        DO ie=1,nelem
          IF (element(ie)%bcid.EQ.iWall) THEN ! loop over elements in a wall
!
!           Find if element is singular
!
            isrc=0
            IF (ie.EQ.irow) THEN
              IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
                isrc=4
              ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad
                isrc=5
              END IF
            END IF
!
!           Perform intergration over a triangle
!           (due to recursive nature of integration, set to zero here)
            ALLOCATE (inteH(element(ie)%nno))
            inteG=0.0_rk
            DO j=1,element(ie)%nno
              inteH(j)=0.0_rk
            END DO
!
!           Element corners
!
            x1=node(element(ie)%con(1))%x(1)  !p%x(p%ibc(ie,1),1)
            y1=node(element(ie)%con(1))%x(2)  !p%x(p%ibc(ie,1),2)
            z1=node(element(ie)%con(1))%x(3)  !p%x(p%ibc(ie,1),3)

            x2=node(element(ie)%con(2))%x(1)  !p%x(p%ibc(ie,2),1)
            y2=node(element(ie)%con(2))%x(2)  !p%x(p%ibc(ie,2),2)
            z2=node(element(ie)%con(2))%x(3)  !p%x(p%ibc(ie,2),3)

            x3=node(element(ie)%con(3))%x(1)  !p%x(p%ibc(ie,3),1)
            y3=node(element(ie)%con(3))%x(2)  !p%x(p%ibc(ie,3),2)
            z3=node(element(ie)%con(3))%x(3)  !p%x(p%ibc(ie,3),3)

            IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
!
!             Set recursion depth for singular triangles
!
              IntRecDepth=parTriRecur

              CALL Triangle_BEMInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                              subdomain(isd)%normMul(jj)*element(ie)%normal(1), &
                              subdomain(isd)%normMul(jj)*element(ie)%normal(2), &
                              subdomain(isd)%normMul(jj)*element(ie)%normal(3), &
                              element(ie)%area,integ,inteH(1),inteH(2),inteH(3),isrc,IntRecDepth)

            ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad

              x4=node(element(ie)%con(4))%x(1)
              y4=node(element(ie)%con(4))%x(2)
              z4=node(element(ie)%con(4))%x(3)

              CALL Quad_BEMInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sx,sy,sz,isrc, &
                            subdomain(isd)%normMul(jj),integ,inteh)
            ELSE
              CALL WriteToLog("Error :: Element type not supported!")
            END IF
!
!           Distribute results to integral matrices
!
            DO j=1,element(ie)%nno
              rowH(element(ie)%con(j))=rowH(element(ie)%con(j))+inteH(j)
            END DO
            rowG(ie)=integ

            DEALLOCATE (inteH)

          END IF
        END DO ! nbelem in wall
      END DO ! walls in subdomain

!
!       Find c parameter : must be 0.5 since source point is in the middle of element
!
      c=0.0_rk
      DO i=1,nnodes
        c=c-rowH(i)
      END DO
      err = MAX(err,ABS(0.5_rk-c))
!       dodamo del ! ja k vsaki izmed 3 oz. 4 tock
!       element v katerem je izvorna tocka
      ie=irow
!       shape functions for center point (constant flux approximation)
      IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
          shpf=1.0_rk/3.0_rk
      ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad
          shpf=0.25_rk
      END IF
!       razdelim c med vozlisca
      DO j=1,element(ie)%nno
          rowH(element(ie)%con(j))=rowH(element(ie)%con(j))+shpf*c
      END DO

      END SUBROUTINE






!
!     ******************************************************************
!
      RECURSIVE SUBROUTINE Triangle_BEMInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                                nx,ny,nz,area,integ,integ1,integ2,integ3,isrc,n)
!
!     Recursively divides triangle towards the source point
!
!
!                3
!              /   \
!             6 --- 5
!           /   \ /   \
!          1 --- 4 --- 2
!
!     ******************************************************************
      USE Triangle ! provides tipw
      USE mCommon
      IMPLICIT NONE

      REAL(rk) sx,sy,sz ! source point
      REAL(rk) nx,ny,nz ! unit normal on triangle surface
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
      REAL(rk) x4,y4,z4,x5,y5,z5,x6,y6,z6 ! small triangle vertexes
      REAL(rk) integ,area,areaS,integ1,integ2,integ3
      INTEGER isrc,isrc1,isrc2,isrc3,isrc4,n

      IF (isrc.EQ.0.OR.n.LT.1) THEN
!
!       non-singular integrals ( or singular at the end of recursion)
!
        CALL Triangle_3DLapInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                              nx,ny,nz,area,integ,integ1,integ2,integ3)
      ELSE
!
!       singular integrals
!
        x4=0.5_rk*(x1+x2)
        y4=0.5_rk*(y1+y2)
        z4=0.5_rk*(z1+z2)

        x5=0.5_rk*(x3+x2)
        y5=0.5_rk*(y3+y2)
        z5=0.5_rk*(z3+z2)

        x6=0.5_rk*(x1+x3)
        y6=0.5_rk*(y1+y3)
        z6=0.5_rk*(z1+z3)
!
!       Detremine new singular integrals
!
        isrc1=0
        isrc2=0
        isrc3=0
        isrc4=0
        IF (isrc.EQ.1) isrc1=1
        IF (isrc.EQ.2) isrc2=2
        IF (isrc.EQ.3) isrc3=3
        IF (isrc.EQ.4) isrc4=4
!
!       Integrate four smaller pieces
!
        CALL   Triangle_Area(x1,y1,z1,x4,y4,z4,x6,y6,z6,areaS)
        CALL Triangle_BEMInt(x1,y1,z1,x4,y4,z4,x6,y6,z6,sx,sy,sz, &
                            nx,ny,nz,areaS,integ,integ1,integ2,integ3,isrc1,n-1)

        CALL   Triangle_Area(x4,y4,z4,x2,y2,z2,x5,y5,z5,areaS)
        CALL Triangle_BEMInt(x4,y4,z4,x2,y2,z2,x5,y5,z5,sx,sy,sz, &
                            nx,ny,nz,areaS,integ,integ1,integ2,integ3,isrc2,n-1)

        CALL   Triangle_Area(x6,y6,z6,x5,y5,z5,x3,y3,z3,areaS)
        CALL Triangle_BEMInt(x6,y6,z6,x5,y5,z5,x3,y3,z3,sx,sy,sz, &
                            nx,ny,nz,areaS,integ,integ1,integ2,integ3,isrc3,n-1)

        CALL   Triangle_Area(x4,y4,z4,x5,y5,z5,x6,y6,z6,areaS)
        CALL Triangle_BEMInt(x4,y4,z4,x5,y5,z5,x6,y6,z6,sx,sy,sz, &
                            nx,ny,nz,areaS,integ,integ1,integ2,integ3,isrc4,n-1)

      END IF

      END
!


!
!     ******************************************************************
!
      SUBROUTINE Triangle_3DLapInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                                nx,ny,nz,area,integ,integ1,integ2,integ3)
!
!     Calculates integral over a triangle
!
!     ******************************************************************
      USE Triangle ! provides tipw
      USE mCommon
      IMPLICIT NONE

      REAL(rk) sx,sy,sz ! source point
      REAL(rk) nx,ny,nz ! unit normal on triangle surface
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
      REAL(rk) x,y,z    ! integration point
      REAL(rk) FUNG,FUNH
      REAL(rk) integ,area,integ1,integ2,integ3
      INTEGER i

      !PI=4.0_rk*ATAN(1.0_rk)

!
!     Loop over integration points
!
      DO i=1,tipw%n
!
!       Calculate point in R^3 space, where function must be evaluated
!
        x=tipw%l1(i)*x1+tipw%l2(i)*x2+tipw%l3(i)*x3
        y=tipw%l1(i)*y1+tipw%l2(i)*y2+tipw%l3(i)*y3
        z=tipw%l1(i)*z1+tipw%l2(i)*z2+tipw%l3(i)*z3

!
!       3D Laplace fundamental solution 1/(4pi r)
!
        CALL BemKernel(sx,sy,sz,x,y,z,nx,ny,nz,funG,funH)
!
!       Sum up integral
!
        integ=integ+tipw%w(i)*FUNG*area ! constant interpolation of flux

        integ1=integ1+tipw%w(i)*FUNH*area*tipw%l1(i) ! linear
        integ2=integ2+tipw%w(i)*FUNH*area*tipw%l2(i) ! interpolation of
        integ3=integ3+tipw%w(i)*FUNH*area*tipw%l3(i) ! function

      END DO ! integration points

      END


!
!     ******************************************************************
!
      SUBROUTINE BemKernel(sx,sy,sz,fx,fy,fz,nx,ny,nz,funG,funH)
!
!     Calculates BEM kernel
!
!     ******************************************************************
      USE mCommon
      IMPLICIT NONE

      REAL(rk) sx,sy,sz ! source point
      REAL(rk) nx,ny,nz ! unit normal on boundary element surface
      REAL(rk) fx,fy,fz ! field (integration) point
      REAL(rk) RX,RY,RZ,RA,IRA2,FUNG,FUNH

!
!     ELLIPTIC LAPLACE FUNDAMENTAL SOLUTION
!
      RX=-fX+sx
      RY=-fY+sy
      RZ=-fZ+sz
      RA=SQRT(RX*RX+RY*RY+RZ*RZ)
      IRA2=1.0_rk/(RA*RA)
      FUNG=0.25_rk/(PI*RA)
      FUNH=( RX*nx+RY*ny+RZ*nz )*IRA2*FUNG


      END
