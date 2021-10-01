!
!     ------------------------------------------------------------------
!
      SUBROUTINE OutputDomainProfile()

      USE mPar
      USE mEqns
      USE mMesh
      USE mCommon
      IMPLICIT NONE

      INTEGER i,j,k,en,ierr,lun
      REAL(rk) Point(3),dist,p
      REAL(rk), ALLOCATABLE :: res(:)
      CHARACTER(255) FileName,FMT

      ALLOCATE (res(neq))
      lun=65

      p = 0.0_rk

      CALL WriteToLog("Calculating domain line profile!")

      DO i=1,parNoDPLE
        WRITE(FileName,'(3A)') "and.dple.",TRIM(parDomLine(i)%name),".txt"
        OPEN (lun,FILE=TRIM(FileName),ERR=10,STATUS='UNKNOWN')

        WRITE(FMT,'(A,I0,A,A)') "(",1+neq,"(A,1X)",",1X,A)"
        WRITE (lun,FMT) "x y z l isd",("u"//TRIM(eqn(k)%name),k=1,neq),"p"

        DO j=1,parDomLine(i)%nno

          DO k=1,3
            Point(k) = parDomLine(i)%startPoint(k) + (parDomLine(i)%endPoint(k) - parDomLine(i)%startPoint(k) ) &
                      / DBLE(parDomLine(i)%nno - 1) * DBLE (j-1)
          END DO

          IF (parPrType .EQ. parLaplace) THEN
            DO en=1,neq
              CALL GetDomainValue(en,Point,res(en),ierr)
            END DO
          ELSE
            CALL GetStokesDV(Point,res,p,ierr)
          END IF

          CALL Distance(Point,parDomLine(i)%startPoint,dist)

          IF (ierr.GT.0) THEN
            WRITE(FMT,'(A,I0,A)') "(4G25.15,1X,I0,1X,",neq+1,"G25.15)"
            WRITE (lun,FMT) Point,dist,ierr,(res(k),k=1,neq),p
          END IF

        END DO

        CLOSE(lun)

      END DO

      DEALLOCATE (res)

      RETURN

10    CONTINUE ! error when opening input file
      CALL WriteToLog("Could not open domain profile export file!")

      END SUBROUTINE




!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetDomainTEST(en)

      USE mEqns
      USE mMesh
      USE mCommon
      IMPLICIT NONE

      INTEGER en,ierr
      REAL(rk) res
      REAL(rk) SourcePoint(3)

      SourcePoint(1)=0.0
      SourcePoint(2)=0.0
      SourcePoint(3)=0.5

      CALL GetDomainValue(en,SourcePoint,res,ierr)


      END SUBROUTINE

!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetDomainValue(en,SourcePoint,res,ierr)

      USE mEqns
      USE mMesh
      USE mCommon
      IMPLICIT NONE

      INTEGER en,isd
      REAL(rk) res,err
      INTEGER ierr
      REAL(rk) SourcePoint(3)

      ierr = 0
      res = -99.99_rk

      CALL FindSubdomain(SourcePoint,isd)
      IF (isd.GT.0) THEN
        CALL GetSubDomainValue(SourcePoint,isd,eqn(en)%u,eqn(en)%q,res,err)
        ierr=isd
        IF (err.GT.1.0E-1) ierr = -1
      ELSE
        ierr = -1
      END IF


      END SUBROUTINE


!
!     ------------------------------------------------------------------
!
      SUBROUTINE FindSubdomain(Point,isdRes)

      USE mMesh
      USE mCommon
      IMPLICIT NONE

      REAL(rk) Point(3)
      REAL(rk) ElementCenter(3),minEC(3),r(3)
      INTEGER isd,ie,in,j,cel,isdRes
      REAL(rk) dist,minDist,multi,dotp

!
!     Find closest element
!
      minDist = 1.0D10
      cel=-1
      DO ie=1,nelem
        ElementCenter = 0.0_rk
        DO in=1,element(ie)%nno
          DO j=1,3
            ElementCenter(j)=ElementCenter(j)+node(element(ie)%con(in))%x(j)
          END DO
        END DO
        DO j=1,3
          ElementCenter(j)=ElementCenter(j)/DBLE(element(ie)%nno)
        END DO
        dist = (Point(1) - ElementCenter(1) ) ** 2 + &
              (Point(2) - ElementCenter(2) ) ** 2 + &
              (Point(3) - ElementCenter(3) ) ** 2
        IF (dist.LT.minDist) THEN
          minDist = dist
          cel = ie
          minEC = ElementCenter
        END IF
      END DO
!
!     Loop over subdomains, do dot product
!
      isdRes = -1
      DO isd = 1,nosd
        CALL GetNormalMulti(isd,cel,multi)
        IF (multi.NE.0.0_rk) THEN  ! element is in this subdomain

          r = minEC - Point

          CALL DotProduct(r,element(cel)%normal,dotp)
          dotp = dotp * multi

          IF (dotp.GE.-1.0D-14) isdRes = isd

        END IF

      END DO

      END SUBROUTINE


!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetNormalMulti(isd,el,multi)

      USE mMesh
      USE mCommon
      IMPLICIT NONE

      INTEGER isd ! subdomain number
      INTEGER el ! element in subdomian
      REAL(rk) multi
      INTEGER jj

      multi=0.0_rk
      DO jj=1,subdomain(isd)%nofw ! find my wall in list of subdomain walls
        IF (element(el)%bcid.EQ.subdomain(isd)%loWalls(jj)) THEN
          multi = subdomain(isd)%normMul(jj)
        END IF
      END DO

      END SUBROUTINE
!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetSubDomainValue(SourcePoint,isd,u,q,res,err)

      USE mMesh
      USE Triangle ! provides tipw
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER isd ! subdomain number
      INTEGER ie,i,j,IntRecDepth,isrc,jj,iWall
      REAL(rk) SourcePoint(3)
      REAL(rk) sx,sy,sz ! source point location
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,c
      REAL(rk) integ,err,multi
      REAL(rk) u(nnodes),q(nqnodes)

      REAL(rk) res


      REAL(rk), ALLOCATABLE :: inteH(:)
      REAL(rk), ALLOCATABLE :: rowH(:),rowG(:)
!
!     Set result to zero
!
      res = 0.0_rk
!
!     Set integrals to zero
!
      ALLOCATE (rowH(nnodes),rowG(nqnodes))
      rowH = 0.0_rk
      rowG = 0.0_rk
!
!     source point location
!
      sx = SourcePoint(1)
      sy = SourcePoint(2)
      sz = SourcePoint(3)
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
!     Estimate error of integration
!
      CALL VerifyIntegrals(SourcePoint,c,isd,rowH,rowG,err)
!
!     Get results
!
      DO i=1,nnodes
        res=res-rowH(i)*u(i)
      END DO
      DO i=1,nelem
        CALL GetNormalMulti(isd,i,multi)
        res=res-rowG(i)*q(i)*multi / subdomain(isd)%diff
      END DO
      res=res/c


      DEALLOCATE (rowH,rowG)

      END SUBROUTINE





!
!     ------------------------------------------------------------------
!
      SUBROUTINE VerifyIntegrals(SourcePoint,cc,isd,rowH,rowG,err)
!
!     ------------------------------------------------------------------
      USE mMesh
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER isd ! subdomain number
      REAL(rk) rowH(nnodes),rowG(nqnodes) ! integrals
      REAL(rk) err,c,multi,cc
      REAL(rk) SourcePoint(3)
      INTEGER i,j,jj

!
!       Check solution : u=x, q=(1,0,0)*(nx,ny,nz)=nx
!
      err=0.0_rk
      DO j=1,3
          c=0.0_rk
          DO i=1,nnodes
            c=c+rowH(i)*node(i)%x(j)
          END DO
          c=c+SourcePoint(j)*cc
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
