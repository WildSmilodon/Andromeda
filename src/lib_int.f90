!
!     ------------------------------------------------------------------
!
      SUBROUTINE CoR_Integrals()

      USE mMesh
      USE mPar
      USE mCommon

      IMPLICIT NONE

      INTEGER ierr,isd,i
      REAL(rk) uSize,qSize,tSize

!
!     Report sizes to log file
!
      tSize = 0.0_rk
      DO isd=1,nosd

        uSize = DBLE(subdomain(isd)%nsp)*DBLE(subdomain(isd)%nnodes)*8.0_rk/1024.0_rk/1024.0_rk
        WRITE (parLogTekst,'(A,A,I0,1X,I0,A,F10.2,A)') TRIM(subdomain(isd)%name)," H matrix size = ",subdomain(isd)%nsp &
          ,subdomain(isd)%nnodes," = ",uSize," Mb"
        CALL WriteToLog(parLogTekst)

        qSize = DBLE(subdomain(isd)%nsp)*DBLE(subdomain(isd)%nqnodes)*8.0/1024.0/1024.0
        WRITE (parLogTekst,'(A,A,I0,1X,I0,A,F10.2,A)') TRIM(subdomain(isd)%name)," G matrix size = ",subdomain(isd)%nsp &
          ,subdomain(isd)%nqnodes," = ",qSize," Mb"
        CALL WriteToLog(parLogTekst)

        tSize = tSize + uSize + qSize
      END DO

      WRITE (parLogTekst,'(A,F10.2,A)') "total size = ",tSize," Mb"
      CALL WriteToLog(parLogTekst)


!
!     Read
!
      CALL WriteToLog("Verifying integrals file!")
      CALL VerifyIntegralsFile(ierr)
!
!     Are integrals on disk
!
      IF (ierr.NE.0) THEN
!
!       No
!
        CALL WriteToLog("Verify failed - computing integrals!")
        CALL sdFormIntegralMatrices()
        IF (parWriteIntegrals.EQ.parYes) THEN
          CALL WriteToLog("Writing integrals to disk!")
          CALL WriteIntegralsToDisk()
          IF (parConserveMemory.EQ.parYes) THEN
            DO i=1,nosd
              DEALLOCATE (subdomain(i)%Hmat)
              DEALLOCATE (subdomain(i)%Gmat)
            END DO
            parfromDisk = parYes
          ELSE
            parfromDisk = parNo
          END IF
        ELSE
          parfromDisk = parNo
        END IF
      ELSE
!
!       Yes
!
        CALL WriteToLog("Verify successful - no need to compute integrals!")
        IF (parConserveMemory.EQ.parYes) THEN
          parfromDisk = parYes
        ELSE
          CALL WriteToLog("Reading integrals from disk!")
          CALL ReadIntegralsFile(ierr)
          parfromDisk = parNo
        END IF
      END IF

!
!     Report errors to log file
!
      DO isd=1,nosd
        WRITE (parLogTekst,'(A,A,E18.12)') TRIM(subdomain(isd)%name)," fSP c err = ", subdomain(isd)%fspCerr
        CALL WriteToLog(parLogTekst)
        WRITE (parLogTekst,'(A,A,E18.12)') TRIM(subdomain(isd)%name)," qSP c err = ", subdomain(isd)%QspCerr
        CALL WriteToLog(parLogTekst)
      END DO

      END SUBROUTINE


!
!     ------------------------------------------------------------------
!
      SUBROUTINE VerifyIntegralsFile(ierr)


      USE mMesh
      USE mPar

      IMPLICIT NONE

      INTEGER isd,lun,ierr,a,b,c

      ierr=1
      lun=12

      OPEN (lun,FILE=TRIM(parIntegralsFileName),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
!
!     Loop over subdomains
!
      DO isd=1,nosd
        READ (lun,ERR=10) a,b,c
        IF ( (a.NE.subdomain(isd)%nsp) .OR. (subdomain(isd)%nnodes.NE.b) .OR. (subdomain(isd)%nqnodes.NE.c) ) THEN
          GOTO 10
        END IF
      END DO
!
!     Errors
!
      DO isd=1,nosd
        READ (lun,ERR=10) subdomain(isd)%FspCerr,subdomain(isd)%QspCerr
      END DO

      ierr=0
10    CONTINUE
      CLOSE(lun)

      END SUBROUTINE

!
!     ------------------------------------------------------------------
!
      SUBROUTINE ReadIntegralsFile(ierr)


      USE mMesh
      USE mPar

      IMPLICIT NONE

      INTEGER isd,lun,ierr,a,b,c

      ierr=1
      lun=12

      OPEN (lun,FILE=TRIM(parIntegralsFileName),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
!
!     Loop over subdomains
!
      DO isd=1,nosd
        READ (lun,ERR=10) a,b,c
        IF ( (a.NE.subdomain(isd)%nsp) .OR. (subdomain(isd)%nnodes.NE.b) .OR. (subdomain(isd)%nqnodes.NE.c) ) THEN
          GOTO 10
        END IF
      END DO
!
!     Errors
!
      DO isd=1,nosd
        READ (lun,ERR=10) subdomain(isd)%FspCerr,subdomain(isd)%QspCerr
      END DO
!
!     Matrices
!
      DO isd=1,nosd
        ALLOCATE (subdomain(isd)%Hmat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
        ALLOCATE (subdomain(isd)%Gmat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
        CALL RdSMat(subdomain(isd)%Hmat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%Gmat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
      END DO


      ierr=0
10    CONTINUE
      CLOSE(lun)

      END SUBROUTINE

!
!     ------------------------------------------------------------------
!
      SUBROUTINE WriteIntegralsToDisk()

      USE mMesh
      USE mPar

      IMPLICIT NONE

      INTEGER isd,lun

      lun=12

      OPEN (lun,FILE=TRIM(parIntegralsFileName),FORM='UNFORMATTED',STATUS='UNKNOWN')
!
!     Sizes
!
      DO isd=1,nosd
        WRITE (lun) subdomain(isd)%nsp,subdomain(isd)%nnodes,subdomain(isd)%nqnodes
      END DO
!
!     Errors
!
      DO isd=1,nosd
        WRITE (lun) subdomain(isd)%FspCerr,subdomain(isd)%QspCerr
      END DO
!
!     Matrices
!
      DO isd=1,nosd
        CALL WrSMat(subdomain(isd)%Hmat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%Gmat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
      END DO


      CLOSE(lun)


      END SUBROUTINE
!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdFormIntegralMatrices()

      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon
      IMPLICIT NONE

      REAL(rk), ALLOCATABLE :: rowH(:),rowG(:)
      REAL(rk) ura,ura0,err,cptime
      INTEGER isd,irow,icol,i,j,row
!
!     Allocate space for integrals for a single source point
!
      ALLOCATE (rowH(nnodes))
      ALLOCATE (rowG(nqnodes))

      ura0=cptime(0.0)
!
!     Loop over subdomains
!
      DO isd=1,nosd
        row=0
        err=0.0

        WRITE (parLogTekst,'(A,A,A)') "working on: ",TRIM(subdomain(isd)%name)," function source points."
        CALL WriteToLog(parLogTekst)


        ALLOCATE (subdomain(isd)%Hmat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
        ALLOCATE (subdomain(isd)%Gmat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
!
!       function source points
!
        DO i = 1,subdomain(isd)%nnodes
          irow = subdomain(isd)%nodeList(i)  ! source point
!
!         Integrate row
!
          row=row+1
          CALL sdIntRowU(isd,irow,rowH,rowG,err) ! irow = source point ID, samo po subdomainu
!
!         Distribute to columns
!
          DO j = 1,subdomain(isd)%nnodes
            icol = subdomain(isd)%nodeList(j)  ! where to take column value
            subdomain(isd)%Hmat(row,j)=rowH(icol)
          END DO

          DO j = 1,subdomain(isd)%nqnodes
            icol = subdomain(isd)%qnodeList(j)  ! where to take column value
            subdomain(isd)%Gmat(row,j)=rowG(icol)
          END DO

        END DO ! u source point

        subdomain(isd)%fspCerr=err

        WRITE (parLogTekst,'(A,A,A)') "working on: ",TRIM(subdomain(isd)%name)," flux source points."
        CALL WriteToLog(parLogTekst)

!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          irow = subdomain(isd)%qnodeList(i)  ! source point
!
!         Integrate row
!
          row=row+1
          CALL sdIntRowQ(isd,irow,rowH,rowG,err)
!
!         Distribute to columns
!
          DO j = 1,subdomain(isd)%nnodes
            icol = subdomain(isd)%nodeList(j)  ! where to take column value
            subdomain(isd)%Hmat(row,j)=rowH(icol)
          END DO

          DO j = 1,subdomain(isd)%nqnodes
            icol = subdomain(isd)%qnodeList(j)  ! where to take column value
            subdomain(isd)%Gmat(row,j)=rowG(icol)
          END DO


        END DO ! q source point

        subdomain(isd)%QspCerr=err

      END DO ! subdomains

      ura=cptime(ura0)
      WRITE (parLogTekst,'(A,F14.4,A4)') "Integration time = ", ura/60.0_rk," min"
      CALL WriteToLog(parLogTekst)


      DEALLOCATE (rowH,rowG)

      END SUBROUTINE


!
! -----------------------------------------------------------------------------------------
!
      SUBROUTINE Quad_BEMInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,isrc, &
                            multi,integ,inteh)
!
!     $: Integracija Robni 4 tockovni element s 4 tockovno geometrijo
!
! -----------------------------------------------------------------------------
      USE GaussIntegration
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER i,j,k,isrc,isip,ng1s,ng2s,ng1,ng2
!
      REAL(rk) xp,yp,zp
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(rk) xc0,yc0,zc0,xet,yet,zet,xks,yks,zks
      REAL(rk) xx1,yy1
      REAL(rk) ro,roth,th,ajac
      REAL(rk) pi2
      REAL(rk) eti,etj,eta1m,eta2m,eta1p,eta2p
      REAL(rk) anx,any,anz,anx1,any1,anz1
      REAL(rk) fung,funh

      REAL(rk) integ,inteh(4)
      REAL(rk) multi ! flips normals

      REAL(rk) fig4(4)
      REAL(rk) al(4),fii(4),th0(4),th1(4)
      REAL(rk) ksi(5),eta(5) ! prve 4 za funkcijo, peta v sredini za fluks SP
!
!      integral divison
      REAL(rk) a,b,c,d,dex,gii,gij,d1,d2,d3,d4,d13max,d24max,minedge
      INTEGER idivXi,ndivXi,idivEt,ndivEt,ising

!
      DATA ksi / 0.0_rk,  1.0_rk, 1.0_rk, 0.0_rk, 0.5_rk/
      DATA eta / 0.0_rk,  0.0_rk, 1.0_rk, 1.0_rk, 0.5_rk/
!
!*** SET NUMBER OF INTEGRATION POINTS
!
!     singular
      ng1s=gaus%ng1(parQuadIntegSing)
      ng2s=gaus%ng2(parQuadIntegSing)
!     regular
      ng1=gaus%ng1(parQuadIntegRegu)
      ng2=gaus%ng2(parQuadIntegRegu)
!
!***  4 NODE CONTINUOUS BOUNDARY ELEMENT
!
      !PI=2.0_rk*ASIN(1.0_rk)
      PI2=2.0_rk*PI
!
      integ=0.0_rk
      inteh(1)=0.0_rk
      inteh(2)=0.0_rk
      inteh(3)=0.0_rk
      inteh(4)=0.0_rk

      minedge=1.0_rk ! to popravi

!
!     Integral razdelimo tako, da integriramo po kvdadratkih
!
      d1=SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      d2=SQRT((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
      d3=SQRT((x3-x4)**2+(y3-y4)**2+(z3-z4)**2)
      d4=SQRT((x4-x1)**2+(y4-y1)**2+(z4-z1)**2)
      d13max=max(d1,d3)
      d24max=max(d2,d4)

      ndivXi=max(2,INT(d13max/minedge+0.5_rk))
      ndivEt=max(2,INT(d24max/minedge+0.5_rk))
!      print *,ndivXi,ndivEt,INT(d13max/minedge+0.5_rk),INT(d24max/minedge+0.5_rk)
!
!     glavna integracijska zanka
!
      DO idivXi=1,ndivXi
        DO idivEt=1,ndivEt
          a=-1.0_rk+(idivXi-1)*2.0_rk/ndivXi
          b=-1.0_rk+(idivXi)*2.0_rk/ndivXi
          c=-1.0_rk+(idivEt-1)*2.0_rk/ndivEt
          d=-1.0_rk+(idivEt)*2.0_rk/ndivEt
          dex=0.25_rk*(b-a)*(d-c)

          IF (isrc.NE.0) THEN
            XX1=ksi(isrc)  ! med 0 in 1
            XX1=-1.0_rk+2.0_rk*XX1  ! med -1 in 1
            XX1=(XX1-a)/(b-a)  ! med 0 in 1 v intervalu a,b

            YY1=eta(isrc)  ! med 0 in 1
            YY1=-1.0_rk+2.0_rk*YY1  ! med -1 in 1
            YY1=(YY1-c)/(d-c)  ! med 0 in 1 v intervalu c,d

            IF (XX1.GE.0.0_rk.AND.XX1.LE.1.0_rk.AND.YY1.GE.0.0_rk.AND.YY1.LE.1.0_rk) THEN
              ising=1
            ELSE
              ising=0
            END IF
          ELSE
            ising=0
          ENDIF
!
!*** SINGULAR INTEGRALS
!
      IF (ising.NE.0) THEN
!        XX1=ksi(isrc)
!        YY1=eta(isrc)
!
        DO K=1,4
!
          IF(K.EQ.1) THEN
            IF (1.0_rk-XX1.EQ.0.0_rk) GOTO 1000
            AL(K)=SQRT(YY1**2+(1.0_rk-XX1)**2)
            FII(K)=ACOS(YY1/AL(K))
            TH0(K)=1.5_rk*PI+ACOS(YY1/AL(K))
            TH1(K)=PI2+ATAN((1.0_rk-YY1)/(1.0_rk-XX1))
!
          ELSE IF(K.EQ.2) THEN
            IF (1.0_rk-YY1.EQ.0.0_rk) GOTO 1000
            AL(K)=SQRT((1.0_rk-YY1)**2+(1.0_rk-XX1)**2)
            FII(K)=ACOS((1.0_rk-XX1)/AL(K))
            TH0(K)=ASIN((1.0_rk-YY1)/AL(K))
            TH1(K)=0.5_rk*PI+ATAN(XX1/(1.0_rk-YY1))
!
          ELSE IF(K.EQ.3) THEN
            IF (XX1.EQ.0.0_rk) GOTO 1000
            AL(K)=SQRT(XX1**2+(1.0_rk-YY1)**2)
            FII(K)=ACOS((1.0_rk-YY1)/AL(K))
            TH0(K)=0.5_rk*PI+ACOS((1.0_rk-YY1)/AL(K))
            TH1(K)=PI+ATAN(YY1/XX1)

          ELSE
            IF (YY1.EQ.0.0_rk) GOTO 1000
            AL(K)=SQRT(YY1**2+XX1**2)
            FII(K)=ACOS(XX1/AL(K))
            TH0(K)=PI+ACOS(XX1/AL(K))
            TH1(K)=1.5_rk*PI+ATAN((1.0_rk-XX1)/YY1)
          END IF
!
!***      GAUSS INTEGRATION (48 points)
!
          DO I=ng1s,ng2s
            DO J=ng1s,ng2s

              TH=(TH1(K)-TH0(K))*gaus%GI(I)/2.0_rk+(TH1(K)+TH0(K))/2.0_rk
              ROTH=AL(K)*SIN(FII(K))/SIN(TH-TH0(K)+FII(K))
              RO=ROTH*gaus%GI(J)/2.0_rk+ROTH/2.0_rk

              ETI=2.0_rk*(XX1+RO*COS(TH))-1.0_rk
              ETJ=2.0_rk*(YY1+RO*SIN(TH))-1.0_rk

              ETI=a+0.5_rk*(ETI+1.0_rk)*(b-a)
              ETJ=c+0.5_rk*(ETJ+1.0_rk)*(d-c)

              ETA1M=1._rk-ETI
              ETA1P=1._rk+ETI
              ETA2M=1._rk-ETJ
              ETA2P=1._rk+ETJ

              FIG4(1)=0.25_rk*ETA1M*ETA2M
              FIG4(2)=0.25_rk*ETA1P*ETA2M
              FIG4(3)=0.25_rk*ETA1P*ETA2P
              FIG4(4)=0.25_rk*ETA1M*ETA2P

              XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
              YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
              ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4

              XKS=0.25_rk*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
              YKS=0.25_rk*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
              ZKS=0.25_rk*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)

              XET=0.25_rk*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
              YET=0.25_rk*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
              ZET=0.25_rk*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)

              ANX=YKS*ZET-YET*ZKS
              ANY=XET*ZKS-XKS*ZET
              ANZ=XKS*YET-XET*YKS

              AJAC=1.0_rk/SQRT(ANX**2+ANY**2+ANZ**2)

              ANX1=multi*ANX*AJAC
              ANY1=multi*ANY*AJAC
              ANZ1=multi*ANZ*AJAC

              AJAC=dex*(TH1(K)-TH0(K))*ROTH*RO*gaus%OME(I)*gaus%OME(J)/AJAC
!
!             Calculate BEM kernel
!
              CALL BemKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,FUNG,FUNH)

!             calculate H continous shape functions
              DO isip=1,4
                inteh(isip)=inteh(isip)+FUNH*AJAC*FIG4(isip)
              END DO

!             calculate G constant shape functions
              integ = integ + FUNG * AJAC

            END DO
          END DO
 1000   END DO

      ELSE
!
!*** REGULAR INTEGRALS
!
        DO i=ng1,ng2
          gii=a+0.5_rk*(gaus%GI(I)+1.0_rk)*(b-a)
          ETA1M=1.0_rk-gii
          ETA1P=1.0_rk+gii
          DO j=ng1,ng2
            gij=c+0.5_rk*(gaus%GI(J)+1.0_rk)*(d-c)
            ETA2M=1.0_rk-gij
            ETA2P=1.0_rk+gij
!
            FIG4(1)=0.25_rk*ETA1M*ETA2M
            FIG4(2)=0.25_rk*ETA1P*ETA2M
            FIG4(3)=0.25_rk*ETA1P*ETA2P
            FIG4(4)=0.25_rk*ETA1M*ETA2P

            XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
            YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
            ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4

            XKS=0.25_rk*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
            YKS=0.25_rk*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
            ZKS=0.25_rk*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)

            XET=0.25_rk*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
            YET=0.25_rk*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
            ZET=0.25_rk*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)

            ANX=YKS*ZET-YET*ZKS
            ANY=XET*ZKS-XKS*ZET
            ANZ=XKS*YET-XET*YKS

            AJAC=1.0_rk/SQRT(ANX**2+ANY**2+ANZ**2)

            ANX1=multi*ANX*AJAC
            ANY1=multi*ANY*AJAC
            ANZ1=multi*ANZ*AJAC

            AJAC=dex*gaus%OME(I)*gaus%OME(J)/AJAC
!
!           Calculate BEM kernel
!
            CALL BemKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,FUNG,FUNH)

!             calculate H continous shape functions
              DO isip=1,4
                inteh(isip)=inteh(isip)+FUNH*AJAC*FIG4(isip)
              END DO

!             calculate G constant shape functions
              integ = integ + FUNG * AJAC


          END DO
        END DO
        CONTINUE
      END IF

        END DO
      END DO

      END

!
!
! -----------------------------------------------------------------------------------------
!
      SUBROUTINE CalElementArea(e)

      USE mMesh
      USE mPar
      USE Triangle
      USE GaussIntegration
      USE mCommon
      IMPLICIT NONE

      INTEGER i,j,isip,ng1,ng2
!
      REAL(rk)   x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
             xc0,yc0,zc0,xet,yet,zet,xks,yks,zks, &
             ajac, &
             eta1m,eta2m,eta1p,eta2p, &
             anx,any,anz
      REAL(rk)  integral
      REAL(rk)  fig4(4)
!
!      integral divison
      REAL(rk) a,b,c,d,dex,gii,gij

      INTEGER n1,n2,n3
      TYPE(ElementType) e


      IF (e%type.EQ.2) THEN ! 3 node trangle

!         list of nodes in triangle
          n1 = e%con(1)
          n2 = e%con(2)
          n3 = e%con(3)

          CALL Triangle_Area(node(n1)%x(1),node(n1)%x(2),node(n1)%x(3), &
                            node(n2)%x(1),node(n2)%x(2),node(n2)%x(3), & 
                            node(n3)%x(1),node(n3)%x(2),node(n3)%x(3), &
                            e%area) 


      ELSE IF (e%type.EQ.3) THEN ! 4 node quad

        ng1=gaus%ng1(4)
        ng2=gaus%ng2(4)
!
        integral=0.00_rk

        X1 = node(e%con(1))%x(1)
        Y1 = node(e%con(1))%x(2)
        Z1 = node(e%con(1))%x(3)

        X2 = node(e%con(2))%x(1)
        Y2 = node(e%con(2))%x(2)
        Z2 = node(e%con(2))%x(3)

        X3 = node(e%con(3))%x(1)
        Y3 = node(e%con(3))%x(2)
        Z3 = node(e%con(3))%x(3)

        X4 = node(e%con(4))%x(1)
        Y4 = node(e%con(4))%x(2)
        Z4 = node(e%con(4))%x(3)

        a=-1.0_rk !+(idivXi-1)*2.0_rk/ndivXi
        b=1.0_rk !-1.0_rk+(idivXi)*2.0_rk/ndivXi
        c=-1.0_rk !+(idivEt-1)*2.0_rk/ndivEt
        d=1.0_rk !-1.0_rk+(idivEt)*2.0_rk/ndivEt
        dex=0.25_rk*(b-a)*(d-c)
!
        DO i=ng1,ng2
          gii=a+0.5_rk*(gaus%GI(I)+1.0_rk)*(b-a)
          ETA1M=1.0_rk-gii
          ETA1P=1.0_rk+gii
          DO j=ng1,ng2
            gij=c+0.5_rk*(gaus%GI(J)+1.0_rk)*(d-c)
            ETA2M=1.0_rk-gij
            ETA2P=1.0_rk+gij
!
            FIG4(1)=0.25_rk*ETA1M*ETA2M
            FIG4(2)=0.25_rk*ETA1P*ETA2M
            FIG4(3)=0.25_rk*ETA1P*ETA2P
            FIG4(4)=0.25_rk*ETA1M*ETA2P

            XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
            YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
            ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4

            XKS=0.25_rk*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
            YKS=0.25_rk*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
            ZKS=0.25_rk*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)

            XET=0.25_rk*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
            YET=0.25_rk*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
            ZET=0.25_rk*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)

            ANX=YKS*ZET-YET*ZKS
            ANY=XET*ZKS-XKS*ZET
            ANZ=XKS*YET-XET*YKS

            AJAC=dex*gaus%OME(I)*gaus%OME(J)*SQRT(ANX**2+ANY**2+ANZ**2)

            DO isip=1,4
              integral=integral+AJAC*FIG4(isip)
            END DO

          END DO
        END DO

        e%area = integral


      ELSE
        CALL WriteToLog("Error :: Element type not supported!")
      END IF

      END


!
!     ------------------------------------------------------------------
!
      SUBROUTINE IntegrateFluxes()

        USE mEqns
        USE mPar
        IMPLICIT NONE
        integer en
        character(255) vrsta
      
        call WriteToLog("")
        call WriteToLog("Integration of flux through the walls")
        DO en=1,neq
          write(vrsta,'(A,A)') "Equation: ", trim(eqn(en)%name)
          call WriteToLog(vrsta)
          call CalFluxIntegralWalls(eqn(en)%q,eqn(en)%name)
        END DO ! en
      
        call WriteToLog("")
      
      end subroutine
      


!
!     ------------------------------------------------------------------
!
SUBROUTINE IntegrateTorques()
  USE mEqns
  USE mPar
  USE mMesh
  IMPLICIT NONE
  integer j
  character(255) vrsta
  real(rk) T(3)

  call WriteToLog("")
  call WriteToLog("Integration of cross product of flux through the walls and r")

  DO j=1,nofw
    call getTorqueIntegrateFluxes(j,T)
    write(vrsta,'(A,3G20.10)') trim(wall(j)%name),T
    call WriteToLog(vrsta)
  END DO
  call WriteToLog("")

end subroutine
      

      !
      !     ------------------------------------------------------------------
      !
      SUBROUTINE IntegrateFunction()
      
        USE mEqns
        USE mPar
        IMPLICIT NONE
        integer en
        character(255) vrsta
      
        call WriteToLog("")
        call WriteToLog("Integration of function over the walls")
        DO en=1,neq
          write(vrsta,'(A,A)') "Equation: ", trim(eqn(en)%name)
          call WriteToLog(vrsta)
          call CalFunctionIntegralWalls(eqn(en)%u,eqn(en)%name)
        END DO ! en
      
        call WriteToLog("")
      
      end subroutine
      

!
!
! -----------------------------------------------------------------------------------------
!
      SUBROUTINE Quad_FunctionInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,f1,f2,f3,f4,integral)

        USE mMesh
        USE mPar
        USE Triangle
        USE GaussIntegration
        USE mCommon
        IMPLICIT NONE
  
        INTEGER i,j,ng1,ng2
  !
        REAL(rk)   x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,f1,f2,f3,f4,integral, &
               xc0,yc0,zc0,xks,yks,zks,xet,yet,zet,ajac, &
               eta1m,eta2m,eta1p,eta2p, &
               anx,any,anz
        REAL(rk)  fig4(4)
  !
  !      integral divison
        REAL(rk) a,b,c,d,dex,gii,gij
  
 
          ng1=gaus%ng1(4)
          ng2=gaus%ng2(4)
  !
          integral=0.00_rk

          a=-1.0_rk !+(idivXi-1)*2.0_rk/ndivXi
          b=1.0_rk !-1.0_rk+(idivXi)*2.0_rk/ndivXi
          c=-1.0_rk !+(idivEt-1)*2.0_rk/ndivEt
          d=1.0_rk !-1.0_rk+(idivEt)*2.0_rk/ndivEt
          dex=0.25_rk*(b-a)*(d-c)
  !
          DO i=ng1,ng2
            gii=a+0.5_rk*(gaus%GI(I)+1.0_rk)*(b-a)
            ETA1M=1.0_rk-gii
            ETA1P=1.0_rk+gii
            DO j=ng1,ng2
              gij=c+0.5_rk*(gaus%GI(J)+1.0_rk)*(d-c)
              ETA2M=1.0_rk-gij
              ETA2P=1.0_rk+gij
  !
              FIG4(1)=0.25_rk*ETA1M*ETA2M
              FIG4(2)=0.25_rk*ETA1P*ETA2M
              FIG4(3)=0.25_rk*ETA1P*ETA2P
              FIG4(4)=0.25_rk*ETA1M*ETA2P
  
              XC0=FIG4(1)*X1+FIG4(2)*X2+FIG4(3)*X3+FIG4(4)*X4
              YC0=FIG4(1)*Y1+FIG4(2)*Y2+FIG4(3)*Y3+FIG4(4)*Y4
              ZC0=FIG4(1)*Z1+FIG4(2)*Z2+FIG4(3)*Z3+FIG4(4)*Z4
  
              XKS=0.25_rk*(-ETA2M*X1+ETA2M*X2+ETA2P*X3-ETA2P*X4)
              YKS=0.25_rk*(-ETA2M*Y1+ETA2M*Y2+ETA2P*Y3-ETA2P*Y4)
              ZKS=0.25_rk*(-ETA2M*Z1+ETA2M*Z2+ETA2P*Z3-ETA2P*Z4)
  
              XET=0.25_rk*(-ETA1M*X1-ETA1P*X2+ETA1P*X3+ETA1M*X4)
              YET=0.25_rk*(-ETA1M*Y1-ETA1P*Y2+ETA1P*Y3+ETA1M*Y4)
              ZET=0.25_rk*(-ETA1M*Z1-ETA1P*Z2+ETA1P*Z3+ETA1M*Z4)
  
              ANX=YKS*ZET-YET*ZKS
              ANY=XET*ZKS-XKS*ZET
              ANZ=XKS*YET-XET*YKS
  
              AJAC=dex*gaus%OME(I)*gaus%OME(J)*SQRT(ANX**2+ANY**2+ANZ**2)
  
              integral = integral + AJAC * FIG4(1) * f1
              integral = integral + AJAC * FIG4(2) * f2
              integral = integral + AJAC * FIG4(3) * f3
              integral = integral + AJAC * FIG4(4) * f4
  
            END DO
          END DO
  
  
        END
  
  


! *********************************************************************
! **                                                                 **
! ** Squeeze multiple blanks into single blank                       **
! **                                                                 **
! *********************************************************************
!---------------------------------------------------------------------!
      SUBROUTINE sqblnk(lun,line)
!---------------------------------------------------------------------!
      CHARACTER*(*) line
      LOGICAL flag
      INTEGER i,ii,j,lun !,lnblnk

      j=1
      flag=.FALSE.
      DO i=1,LNBLNK(line)+1
        IF (line(i:i).NE.' ' .AND. .NOT.flag) THEN
          flag=.TRUE.
          ii=i
        ELSE IF (line(i:i).EQ.' ' .AND. flag) THEN
          flag=.FALSE.
          line(j:j+i-ii)=line(ii:i)
          j=j+i-ii+1
        END IF
      END DO
      WRITE(lun,'(A)') line(1:j-2)
      RETURN
      END
!      INTEGER FUNCTION LNBLNK (string) - je del fortrana
!!
!!     LNBLNK returns the index of the last non-blank character in string
!!
!      CHARACTER*(*) string
!      INTEGER i
!
!      DO i=LEN(string),1,-1
!        IF (string(i:i).NE.' ') THEN
!          lnblnk=i
!          RETURN
!        END IF
!      END DO
!      lnblnk=0
!      RETURN
!      END

