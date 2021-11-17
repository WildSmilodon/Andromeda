!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdSolveCRS(en)

      USE mMesh
      USE mCommon
      USE mEqns
      USE mPar
      IMPLICIT NONE

      REAL(rk), ALLOCATABLE :: b(:)
      INTEGER en,i,isp,isd,myWall

      INTEGER nits,ierr
      REAL(rk)    cput

      WRITE (parLogTekst,'(A,A)') "Solving for ",TRIM(eqn(en)%name)
      CALL WriteToLog(parLogTekst)

!
!     Set up x, b and rhs-vector vectors
!
      ALLOCATE (b(eqn(en)%neq))
!
!     Loop over subdomains
!
      DO isd=1,nosd
!
!       function source points
!
        DO i = 1,subdomain(isd)%nnodes
          isp = subdomain(isd)%nodeList(i)  ! source point
          myWall = subdomain(isd)%BCidList(i)

          IF (eqn(en)%boundary(myWall)%known.EQ.iFunction) THEN
              eqn(en)%b(eqn(en)%col(isp))=eqn(en)%u(isp)
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iFlux) THEN
              eqn(en)%x(eqn(en)%col(isp))=eqn(en)%u(isp)
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iContact) THEN
              eqn(en)%x(eqn(en)%col(isp))=eqn(en)%u(isp)
          END IF

        END DO ! u source point
!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          isp = subdomain(isd)%qnodeList(i)  ! source point
          myWall = subdomain(isd)%qBCidList(i)

          IF (eqn(en)%boundary(myWall)%known.EQ.iFunction) THEN
            eqn(en)%x(eqn(en)%qcol(isp))=eqn(en)%q(isp)
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iFlux) THEN
            eqn(en)%b(eqn(en)%qcol(isp))=eqn(en)%q(isp)
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iContact) THEN
            eqn(en)%x(eqn(en)%qcol(isp))=eqn(en)%q(isp)
          END IF

        END DO ! q source point
      END DO ! subdomains

!
!     Calculate rhs vector
!
      CALL CRSxV(eqn(en)%rhsMcrs,-eqn(en)%b,eqn(en)%nb,b)
!
!     Solve
!
      ierr=0
      cput=0.0_rk
      nits=0

!
!     preconditioner (calculate : eqn(en)%Pivot )
!
      CALL lsqr_precon(eqn(en)%sysMcrs%neq,eqn(en)%nx,eqn(en)%sysMcrs%nnz,eqn(en)%Pivot, &
                       eqn(en)%sysMcrs%i,eqn(en)%sysMcrs%j,eqn(en)%sysMcrs%v)

!
!     solve overdetermined system of linear equations
!
      CALL lsqr_solve(eqn(en)%sysMcrs%neq,eqn(en)%nx,eqn(en)%sysMcrs%nnz, &
                      eqn(en)%slv%maxit,eqn(en)%slv%eps,nits,ierr,eqn(en)%Pivot, &
                      eqn(en)%sysMcrs%i,eqn(en)%sysMcrs%j,eqn(en)%sysMcrs%v,b,eqn(en)%x)

!
!       Report to log file
!
      WRITE (parLogTekst,'(A,A,I0,1X,I0)') TRIM(eqn(en)%name)," LSQR Solver : ", ierr,nits
      CALL WriteToLog(parLogTekst)
      IF (eqn(en)%slv%maxit.LE.nits) THEN
        WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
        CALL WriteToLog(parLogTekst)
      END IF
!
!     Distribute SLE results to u and q vectors
!
!
!
!     Loop over subdomains
!
      DO isd=1,nosd
!
!       function source points
!
        DO i = 1,subdomain(isd)%nnodes
          isp = subdomain(isd)%nodeList(i)  ! source point
          myWall = subdomain(isd)%BCidList(i)

          IF (eqn(en)%boundary(myWall)%known.EQ.iFunction) THEN
              eqn(en)%u(isp)=eqn(en)%b(eqn(en)%col(isp))
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iFlux) THEN
              eqn(en)%u(isp)=eqn(en)%x(eqn(en)%col(isp))
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iContact) THEN
              eqn(en)%u(isp)=eqn(en)%x(eqn(en)%col(isp))
          END IF

        END DO ! u source point
!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          isp = subdomain(isd)%qnodeList(i)  ! source point
          myWall = subdomain(isd)%qBCidList(i)

          IF (eqn(en)%boundary(myWall)%known.EQ.iFunction) THEN
            eqn(en)%q(isp)=eqn(en)%x(eqn(en)%qcol(isp))
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iFlux) THEN
            eqn(en)%q(isp)=eqn(en)%b(eqn(en)%qcol(isp))
          ELSE IF (eqn(en)%boundary(myWall)%known.EQ.iContact) THEN
            eqn(en)%q(isp)=eqn(en)%x(eqn(en)%qcol(isp))
          END IF

        END DO ! q source point
      END DO ! subdomains
!
!     Free memory
!
      DEALLOCATE (b)

      END


!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdSetUpSysMrhsM_CRS(en)

      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon

      IMPLICIT NONE

      INTEGER en,row,mywall,isd,i,j,r,col,mywallC,k,bc
      REAL(rk) val,multi
      INTEGER(8) nnzSys,nnzRHS

      nnzSys=0
      nnzRHS=0
!
!     Count NNZ in CRS matrices
!
      CALL sdSetUpSysMrhsM_CRScount(en,eqn(en)%sysMcrs%nnz,eqn(en)%rhsMcrs%nnz)
      eqn(en)%sysMcrs%neq = eqn(en)%neq
      eqn(en)%rhsMcrs%neq = eqn(en)%neq

      WRITE (parLogTekst,'(A,A,F10.2,A)') TRIM(eqn(en)%name)," system matrix size = ", &
      (DBLE(eqn(en)%sysMcrs%nnz)+0.5*DBLE(eqn(en)%sysMcrs%neq))*8.0/1024.0/1024.0," Mb"
      CALL WriteToLog(parLogTekst)

      WRITE (parLogTekst,'(A,A,F10.2,A)') TRIM(eqn(en)%name)," RHS    matrix size = ", &
      (DBLE(eqn(en)%rhsMcrs%nnz)+0.5*DBLE(eqn(en)%rhsMcrs%neq))*8.0/1024.0/1024.0," Mb"
      CALL WriteToLog(parLogTekst)



!
!     allocate memory for CRS matrix
!
      ALLOCATE (eqn(en)%sysMcrs%v(eqn(en)%sysMcrs%nnz),eqn(en)%sysMcrs%i(eqn(en)%sysMcrs%neq+1))
      ALLOCATE (eqn(en)%sysMcrs%j(eqn(en)%sysMcrs%nnz))

      ALLOCATE (eqn(en)%rhsMcrs%v(eqn(en)%rhsMcrs%nnz),eqn(en)%rhsMcrs%i(eqn(en)%rhsMcrs%neq+1))
      ALLOCATE (eqn(en)%rhsMcrs%j(eqn(en)%rhsMcrs%nnz))
!
!     Allocate memory for row
!

      row = 0  ! row in sysM and rhsM
!
!     Loop over subdomains
!
      DO isd=1,nosd
!
!       function source points
!
        r=0 ! row in H and G matrix

        DO i = 1,subdomain(isd)%nnodes
          myWall = subdomain(isd)%BCidList(i)
          bc = eqn(en)%boundary(myWall)%known

          r=r+1
!
!         Determine if I need equation for this source point
!
          IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN

            row=row+1
            eqn(en)%sysMcrs%i(row) = nnzSys + 1
            eqn(en)%rhsMcrs%i(row) = nnzRHS + 1
!
!           H values
!
            DO j = 1,subdomain(isd)%nnodes ! loop over H matrix columns
              col =eqn(en)%col(subdomain(isd)%nodeList(j))
              myWallC = subdomain(isd)%BCidList(j)
              bc = eqn(en)%boundary(myWallC)%known
              val = subdomain(isd)%Hmat(r,j)

              IF (bc.EQ.iFunction) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF


              IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF

            END DO
!
!           G values
!
            DO j = 1,subdomain(isd)%nqnodes ! loop over G matrix columns
              col =eqn(en)%qcol(subdomain(isd)%qnodeList(j))
              myWallC = subdomain(isd)%qBCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              DO k=1,subdomain(isd)%nofw
                IF (myWallC.EQ.subdomain(isd)%loWalls(k)) multi = subdomain(isd)%normMul(k)
              END DO
              val = multi / subdomain(isd)%diff * subdomain(isd)%Gmat(r,j)

              IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF
              IF (bc.EQ.iFlux) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF

            END DO

          END IF
        END DO
!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          myWall = subdomain(isd)%qBCidList(i)
          bc = eqn(en)%boundary(myWall)%known
!
!         Integrate row
!
          r=r+1
!
!         Determine if I need equation for this source point
!
          IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN

            row=row+1
            eqn(en)%sysMcrs%i(row) = nnzSys + 1
            eqn(en)%rhsMcrs%i(row) = nnzRHS + 1
!
!           H values
!
            DO j = 1,subdomain(isd)%nnodes
              col =eqn(en)%col(subdomain(isd)%nodeList(j))
              myWallC = subdomain(isd)%BCidList(j)
              bc = eqn(en)%boundary(myWallC)%known
              val = subdomain(isd)%Hmat(r,j)

              IF (bc.EQ.iFunction) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF


              IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF

            END DO
!
!           G values
!
            DO j = 1,subdomain(isd)%nqnodes
              col =eqn(en)%qcol(subdomain(isd)%qnodeList(j))
              myWallC = subdomain(isd)%qBCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              DO k=1,subdomain(isd)%nofw ! define q based on sub. normal
                IF (myWallC.EQ.subdomain(isd)%loWalls(k)) multi = subdomain(isd)%normMul(k)
              END DO
              val = multi / subdomain(isd)%diff * subdomain(isd)%Gmat(r,j)

              IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF
              IF (bc.EQ.iFlux) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF

            END DO

          END IF
        END DO

      END DO ! subdomains
!
!     Additional nonexistent start of line
!
      eqn(en)%sysMcrs%i(row+1) = nnzSys + 1
      eqn(en)%rhsMcrs%i(row+1) = nnzRHS + 1

      END SUBROUTINE



!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdSetUpSysMrhsM_CRS_fromDisk(en)

      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon

      IMPLICIT NONE

      INTEGER en,row,mywall,isd,i,j,r,col,mywallC,k,bc
      REAL(rk) val,multi
      INTEGER(8) nnzSys,nnzRHS
      INTEGER lun,a,b,c

      lun=12
!
!     Open integrals file
!
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

      nnzSys=0
      nnzRHS=0
!
!     Count NNZ in CRS matrices
!
      CALL sdSetUpSysMrhsM_CRScount(en,eqn(en)%sysMcrs%nnz,eqn(en)%rhsMcrs%nnz)
      eqn(en)%sysMcrs%neq = eqn(en)%neq
      eqn(en)%rhsMcrs%neq = eqn(en)%neq

      WRITE (parLogTekst,'(A,A,F10.2,A)') TRIM(eqn(en)%name)," system matrix size = ", &
      (DBLE(eqn(en)%sysMcrs%nnz)+0.5*DBLE(eqn(en)%sysMcrs%neq))*8.0/1024.0/1024.0," Mb"
      CALL WriteToLog(parLogTekst)

      WRITE (parLogTekst,'(A,A,F10.2,A)') TRIM(eqn(en)%name)," RHS    matrix size = ", &
      (DBLE(eqn(en)%rhsMcrs%nnz)+0.5*DBLE(eqn(en)%rhsMcrs%neq))*8.0/1024.0/1024.0," Mb"
      CALL WriteToLog(parLogTekst)



!
!     allocate memory for CRS matrix
!
      ALLOCATE (eqn(en)%sysMcrs%v(eqn(en)%sysMcrs%nnz),eqn(en)%sysMcrs%i(eqn(en)%sysMcrs%neq+1))
      ALLOCATE (eqn(en)%sysMcrs%j(eqn(en)%sysMcrs%nnz))

      ALLOCATE (eqn(en)%rhsMcrs%v(eqn(en)%rhsMcrs%nnz),eqn(en)%rhsMcrs%i(eqn(en)%rhsMcrs%neq+1))
      ALLOCATE (eqn(en)%rhsMcrs%j(eqn(en)%rhsMcrs%nnz))
!
!     Allocate memory for row
!

      row = 0  ! row in sysM and rhsM
!
!     Loop over subdomains
!
      DO isd=1,nosd

!
!       Read integrals from disk
!
        ALLOCATE (subdomain(isd)%Hmat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
        ALLOCATE (subdomain(isd)%Gmat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
        CALL RdSMat(subdomain(isd)%Hmat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%Gmat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)


!
!       function source points
!
        r=0 ! row in H and G matrix

        DO i = 1,subdomain(isd)%nnodes
          myWall = subdomain(isd)%BCidList(i)
          bc = eqn(en)%boundary(myWall)%known

          r=r+1
!
!         Determine if I need equation for this source point
!
          IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN

            row=row+1
            eqn(en)%sysMcrs%i(row) = nnzSys + 1
            eqn(en)%rhsMcrs%i(row) = nnzRHS + 1
!
!           H values
!
            DO j = 1,subdomain(isd)%nnodes ! loop over H matrix columns
              col =eqn(en)%col(subdomain(isd)%nodeList(j))
              myWallC = subdomain(isd)%BCidList(j)
              bc = eqn(en)%boundary(myWallC)%known
              val = subdomain(isd)%Hmat(r,j)

              IF (bc.EQ.iFunction) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF


              IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF

            END DO
!
!           G values
!
            DO j = 1,subdomain(isd)%nqnodes ! loop over G matrix columns
              col =eqn(en)%qcol(subdomain(isd)%qnodeList(j))
              myWallC = subdomain(isd)%qBCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              DO k=1,subdomain(isd)%nofw
                IF (myWallC.EQ.subdomain(isd)%loWalls(k)) multi = subdomain(isd)%normMul(k)
              END DO
              val = multi / subdomain(isd)%diff * subdomain(isd)%Gmat(r,j)

              IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF
              IF (bc.EQ.iFlux) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF

            END DO

          END IF
        END DO
!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          myWall = subdomain(isd)%qBCidList(i)
          bc = eqn(en)%boundary(myWall)%known
!
!         Integrate row
!
          r=r+1
!
!         Determine if I need equation for this source point
!
          IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN

            row=row+1
            eqn(en)%sysMcrs%i(row) = nnzSys + 1
            eqn(en)%rhsMcrs%i(row) = nnzRHS + 1
!
!           H values
!
            DO j = 1,subdomain(isd)%nnodes
              col =eqn(en)%col(subdomain(isd)%nodeList(j))
              myWallC = subdomain(isd)%BCidList(j)
              bc = eqn(en)%boundary(myWallC)%known
              val = subdomain(isd)%Hmat(r,j)

              IF (bc.EQ.iFunction) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF


              IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF

            END DO
!
!           G values
!
            DO j = 1,subdomain(isd)%nqnodes
              col =eqn(en)%qcol(subdomain(isd)%qnodeList(j))
              myWallC = subdomain(isd)%qBCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              DO k=1,subdomain(isd)%nofw ! define q based on sub. normal
                IF (myWallC.EQ.subdomain(isd)%loWalls(k)) multi = subdomain(isd)%normMul(k)
              END DO
              val = multi / subdomain(isd)%diff * subdomain(isd)%Gmat(r,j)

              IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                eqn(en)%sysMcrs%v(nnzSys) = val
                eqn(en)%sysMcrs%j(nnzSys) = col
              END IF
              IF (bc.EQ.iFlux) THEN
                nnzRHS = nnzRHS + 1
                eqn(en)%rhsMcrs%v(nnzRHS) = val
                eqn(en)%rhsMcrs%j(nnzRHS) = col
              END IF

            END DO

          END IF
        END DO
!
!       Free memory (deallocate integrals for this subdomain)
!
        DEALLOCATE (subdomain(isd)%Hmat)
        DEALLOCATE (subdomain(isd)%Gmat)

      END DO ! subdomains

      CLOSE(lun)

!
!     Additional nonexistent start of line
!
      eqn(en)%sysMcrs%i(row+1) = nnzSys + 1
      eqn(en)%rhsMcrs%i(row+1) = nnzRHS + 1

      RETURN
10    CONTINUE

      CALL WriteToLog("ERROR in sdSetUpSysMrhsM_CRS_fromDisk!")
      CALL StopProgram

      END SUBROUTINE


!
!     ------------------------------------------------------------------
!
      SUBROUTINE sdSetUpSysMrhsM_CRScount(en,nnzSys,nnzRHS)

      USE mMesh
      USE mEqns
      USE mPar

      IMPLICIT NONE

      INTEGER en,row,mywall,isd,i,j,r,col,mywallC,bc
      INTEGER(8) nnzSys,nnzRHS

      nnzSys=0
      nnzRHS=0

      row = 0  ! row in sysM and rhsM
!
!     Loop over subdomains
!
      DO isd=1,nosd
!
!       function source points
!
        r=0 ! row in H and G matrix

        DO i = 1,subdomain(isd)%nnodes
          myWall = subdomain(isd)%BCidList(i)
          bc = eqn(en)%boundary(myWall)%known

          r=r+1
!
!         Determine if I need equation for this source point
!
          IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN

            row=row+1
!
!           H values
!
            DO j = 1,subdomain(isd)%nnodes ! loop over H matrix columns
              col =eqn(en)%col(subdomain(isd)%nodeList(j))
              myWallC = subdomain(isd)%BCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzRHS = nnzRHS + 1
              IF (bc.EQ.iFlux)     nnzSys = nnzSys + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1

            END DO
!
!           G values
!
            DO j = 1,subdomain(isd)%nqnodes ! loop over G matrix columns
              col =eqn(en)%qcol(subdomain(isd)%qnodeList(j))
              myWallC = subdomain(isd)%qBCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzSys = nnzSys + 1
              IF (bc.EQ.iFlux)     nnzRHS = nnzRHS + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1

            END DO

          END IF
        END DO
!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          myWall = subdomain(isd)%qBCidList(i)
          bc = eqn(en)%boundary(myWall)%known
!
!         Integrate row
!
          r=r+1
!
!         Determine if I need equation for this source point
!
          IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN

            row=row+1
!
!           H values
!
            DO j = 1,subdomain(isd)%nnodes
              col =eqn(en)%col(subdomain(isd)%nodeList(j))
              myWallC = subdomain(isd)%BCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzRHS = nnzRHS + 1
              IF (bc.EQ.iFlux)     nnzSys = nnzSys + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1

            END DO
!
!           G values
!
            DO j = 1,subdomain(isd)%nqnodes
              col =eqn(en)%qcol(subdomain(isd)%qnodeList(j))
              myWallC = subdomain(isd)%qBCidList(j)
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzSys = nnzSys + 1
              IF (bc.EQ.iFlux)     nnzRHS = nnzRHS + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1

            END DO

          END IF
        END DO

      END DO ! subdomains

      END SUBROUTINE


!
!     ------------------------------------------------------------------
!
      SUBROUTINE TransformToCRS(mat,nrow,ncol,crsM)

      USE mEqns
      USE mCommon

      IMPLICIT NONE

      INTEGER nrow,ncol,i,j
      REAL(rk) mat(nrow,ncol)
      TYPE(CRSmatrixType) crsM

!
!     count
!
      crsM%nnz = 0
      crsM%neq = nrow

      DO i=1,nrow
        DO j=1,ncol
          IF (mat(i,j).NE.0.0_rk) THEN
            crsM%nnz = crsM%nnz + 1
          END IF
        END DO
      END DO
!
!     allocate memory for CRS matrix
!
      ALLOCATE(crsM%v(crsM%nnz),crsM%i(crsM%neq+1),crsM%j(crsM%nnz))
!
!     transform
!
      crsM%nnz = 0
      DO i=1,nrow
        crsM%i(i)=crsM%nnz + 1 ! line start
        DO j=1,ncol
          IF (mat(i,j).NE.0.0_rk) THEN
            crsM%nnz = crsM%nnz + 1
            crsM%v(crsM%nnz) = mat(i,j)  ! values
            crsM%j(crsM%nnz) = j         ! column
          END IF
        END DO
      END DO
      crsM%i(crsM%neq + 1)=crsM%nnz + 1 ! start of non-existent line


      END SUBROUTINE


! ----------------------------------------------------------------------
      SUBROUTINE CRSxV(mat,vek,DolVek,rez)
!
!     $: mnozi CRS matriko z vektorjem
!
! ----------------------------------------------------------------------
      USE mEqns
      USE mCommon

      TYPE(CRSmatrixType) :: mat
      INTEGER DolVek
      INTEGER i
      INTEGER(8) j

      REAL(rk) vek(DolVek)
      REAL(rk) rez(mat%neq)

      DO i=1,mat%neq
        rez(i)=0.0_rk
        DO j=mat%i(i),mat%i(i+1)-1
          rez(i)=rez(i)+mat%v(j)*vek(mat%j(j))
        END DO
      END DO

      END SUBROUTINE

