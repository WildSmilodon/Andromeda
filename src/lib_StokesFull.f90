!
!     ------------------------------------------------------------------
!
SUBROUTINE formQrow(isd,een,r,sysMrow,rhsMrow)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE
  
  REAL(rk) sysMrow(stk%neq),rhsMrow(stk%nb)
  integer j,k,isd,isp,mywallC,en,ispL,col,bc,een,r
  real(rk) val,multi

  sysMrow = 0.0_rk
  rhsMrow = 0.0_rk
  !
  ! T matrix values
  !
  DO j = 1,subdomain(isd)%nnodes ! loop over T matrix columns
    myWallC = subdomain(isd)%BCidList(j)
    isp = subdomain(isd)%nodeList(j)
    DO en=1,neq
      ispL = isp + (en-1)*nnodes ! source point number in the large system     
      col  = stk%col(ispL)
      bc   = eqn(en)%boundary(myWallC)%known

      IF (een.EQ.1.AND.en.EQ.1) val = subdomain(isd)%TxxMat(r,j)
      IF (een.EQ.1.AND.en.EQ.2) val = subdomain(isd)%TxyMat(r,j)
      IF (een.EQ.1.AND.en.EQ.3) val = subdomain(isd)%TxzMat(r,j)
      IF (een.EQ.2.AND.en.EQ.1) val = subdomain(isd)%TxyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.2) val = subdomain(isd)%TyyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.3) val = subdomain(isd)%TyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.1) val = subdomain(isd)%TxzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.2) val = subdomain(isd)%TyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.3) val = subdomain(isd)%TzzMat(r,j)
      
      IF (bc.EQ.iFunction) THEN
        rhsMrow(col) = val
      END IF
      IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
        sysMrow(col) =  - val                ! MINUS
      END IF

    END DO
  END DO
  !
  ! G matrix values
  !
  DO j = 1,subdomain(isd)%nqnodes ! loop over G matrix columns
    myWallC = subdomain(isd)%qBCidList(j)
    isp = subdomain(isd)%qnodeList(j)
    DO en=1,neq
      ispL = isp + (en-1)*nqnodes ! source point number in the large system     
      col =stk%qcol(ispL)
      bc = eqn(en)%boundary(myWallC)%known

      DO k=1,subdomain(isd)%nofw
        IF (myWallC.EQ.subdomain(isd)%loWalls(k)) multi = subdomain(isd)%normMul(k)
      END DO


      IF (een.EQ.1.AND.en.EQ.1) val = subdomain(isd)%GxxMat(r,j)
      IF (een.EQ.1.AND.en.EQ.2) val = subdomain(isd)%GxyMat(r,j)
      IF (een.EQ.1.AND.en.EQ.3) val = subdomain(isd)%GxzMat(r,j)
      IF (een.EQ.2.AND.en.EQ.1) val = subdomain(isd)%GxyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.2) val = subdomain(isd)%GyyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.3) val = subdomain(isd)%GyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.1) val = subdomain(isd)%GxzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.2) val = subdomain(isd)%GyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.3) val = subdomain(isd)%GzzMat(r,j)
 
      val = multi / subdomain(isd)%diff * val

      IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
        sysMrow(col) = - val  ! XXXXXX
      END IF

      IF (bc.EQ.iFlux) THEN
        rhsMrow(col) =  val
      END IF

    END DO
  END DO  


end subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE formUrow(isd,een,r,sysMrow,rhsMrow)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE
  
  REAL(rk) sysMrow(stk%neq),rhsMrow(stk%nb)
  integer j,k,isd,isp,mywallC,en,ispL,col,bc,een,r
  real(rk) val,multi

  sysMrow = 0.0_rk
  rhsMrow = 0.0_rk
  !
  ! T matrix values
  !
  DO j = 1,subdomain(isd)%nnodes ! loop over T matrix columns
    myWallC = subdomain(isd)%BCidList(j)
    isp = subdomain(isd)%nodeList(j)
    DO en=1,neq
      ispL = isp + (en-1)*nnodes ! source point number in the large system     
      col  = stk%col(ispL)
      bc   = eqn(en)%boundary(myWallC)%known
      IF (een.EQ.1.AND.en.EQ.1) val = subdomain(isd)%TxxMat(r,j)
      IF (een.EQ.1.AND.en.EQ.2) val = subdomain(isd)%TxyMat(r,j)
      IF (een.EQ.1.AND.en.EQ.3) val = subdomain(isd)%TxzMat(r,j)
      IF (een.EQ.2.AND.en.EQ.1) val = subdomain(isd)%TxyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.2) val = subdomain(isd)%TyyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.3) val = subdomain(isd)%TyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.1) val = subdomain(isd)%TxzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.2) val = subdomain(isd)%TyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.3) val = subdomain(isd)%TzzMat(r,j)            
      IF (bc.EQ.iFunction) THEN
        rhsMrow(col) =   val
      END IF
      IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
        sysMrow(col) = - val                ! MINUS    
      END IF
    END DO
  END DO
  !
  ! G matrix values
  !
  DO j = 1,subdomain(isd)%nqnodes ! loop over G matrix columns
    myWallC = subdomain(isd)%qBCidList(j)
    isp = subdomain(isd)%qnodeList(j)
    DO en=1,neq
      ispL = isp + (en-1)*nqnodes ! source point number in the large system     
      col =stk%qcol(ispL)          
      bc = eqn(en)%boundary(myWallC)%known
      DO k=1,subdomain(isd)%nofw
        IF (myWallC.EQ.subdomain(isd)%loWalls(k)) multi = subdomain(isd)%normMul(k)
      END DO
      IF (een.EQ.1.AND.en.EQ.1) val = subdomain(isd)%GxxMat(r,j)
      IF (een.EQ.1.AND.en.EQ.2) val = subdomain(isd)%GxyMat(r,j)
      IF (een.EQ.1.AND.en.EQ.3) val = subdomain(isd)%GxzMat(r,j)
      IF (een.EQ.2.AND.en.EQ.1) val = subdomain(isd)%GxyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.2) val = subdomain(isd)%GyyMat(r,j)
      IF (een.EQ.2.AND.en.EQ.3) val = subdomain(isd)%GyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.1) val = subdomain(isd)%GxzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.2) val = subdomain(isd)%GyzMat(r,j)
      IF (een.EQ.3.AND.en.EQ.3) val = subdomain(isd)%GzzMat(r,j)
 
      val = multi / subdomain(isd)%diff * val
      IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
        sysMrow(col) = - val  ! XXXXXX
      END IF
      IF (bc.EQ.iFlux) THEN
        rhsMrow(col) =   val
      END IF
    END DO
  END DO
end subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE formStokesSysRhsMatrices(sysM,rhsM)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER row,mywall,isd,i,r,col,bc,een
  REAL(rk), ALLOCATABLE :: sysMrow(:),rhsMrow(:)
  real(rk) sysM(stk%neq,stk%neq),rhsM(stk%neq,stk%nb)

  !
  !  Form single rows for system and rhs matrices
  !
  ALLOCATE(sysMrow(stk%neq))
  ALLOCATE(rhsMrow(stk%nb))
 
  row = 0  ! row in sysM and rhsM
  !
  ! Loop over subdomains
  !

  isd = 1 ! single subdomain only

  !
  ! Function source points
  !  
  r=0 ! row in T and G matrices
  DO i = 1,subdomain(isd)%nnodes
    myWall = subdomain(isd)%BCidList(i)
    r = r + 1
    !
    ! Loop over all three small systems
    !
    DO een=1,neq
      bc = eqn(een)%boundary(myWall)%known 
      !
      ! Use this equation only, if function is unknown
      !
      IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
        row = row + 1        
        !
        !  For a single row in matrices
        !
        CALL formUrow(isd,een,r,sysMrow,rhsMrow)
        !
        ! Copy row to full matrix
        !
        DO col = 1,stk%neq
          sysM(row,col) = sysMrow(col)
        END DO
        DO col = 1,stk%nb
          rhsM(row,col) = rhsMrow(col)
        END DO            
      END IF ! I need this equation
    END DO ! en
  END DO ! u source point

  !      
  ! Flux source points
  !        
  DO i = 1,subdomain(isd)%nqnodes
    myWall = subdomain(isd)%qBCidList(i)
    r=r+1
    !
    ! Loop over all three small systems
    !
    DO een=1,neq             
      bc = eqn(een)%boundary(myWall)%known
      !
      ! Use this equation only, if flux is unknown
      !
      IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
        row=row+1      
        !
        !  For a single row in matrices
        !
        CALL formQrow(isd,een,r,sysMrow,rhsMrow)
        !
        ! Copy row to full matrix
        !
        DO col = 1,stk%neq
          sysM(row,col) = sysMrow(col)
        END DO
        DO col = 1,stk%nb
          rhsM(row,col) = rhsMrow(col)
        END DO              
      END IF
    END DO ! en        
  END DO ! q source point

  DEALLOCATE(sysMrow,rhsMrow)

end subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE countStokesSysRhsMatrices(sysMnnz,rhsMnnz)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER row,mywall,isd,i,r,col,bc,een
  REAL(rk), ALLOCATABLE :: sysMrow(:),rhsMrow(:)
  integer sysMnnz,rhsMnnz

  sysMnnz = 0
  rhsMnnz = 0
  !
  !  Form single rows for system and rhs matrices
  !
  ALLOCATE(sysMrow(stk%neq))
  ALLOCATE(rhsMrow(stk%nb))
 
  row = 0  ! row in sysM and rhsM
  !
  ! Loop over subdomains
  !

  isd = 1 ! single subdomain only

  !
  ! Function source points
  !  
  r=0 ! row in T and G matrices
  DO i = 1,subdomain(isd)%nnodes
    myWall = subdomain(isd)%BCidList(i)
    r = r + 1
    !
    ! Loop over all three small systems
    !
    DO een=1,neq
      bc = eqn(een)%boundary(myWall)%known 
      !
      ! Use this equation only, if function is unknown
      !
      IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
        row = row + 1        
        !
        !  For a single row in matrices
        !
        CALL formUrow(isd,een,r,sysMrow,rhsMrow)
        !
        ! Copy row to full matrix
        !
        DO col = 1,stk%neq
          if (sysMrow(col).NE.0.0_rk) sysMnnz = sysMnnz + 1
        END DO
        DO col = 1,stk%nb
          if (rhsMrow(col).NE.0.0_rk) rhsMnnz = rhsMnnz + 1 
        END DO            
      END IF ! I need this equation
    END DO ! en
  END DO ! u source point

  !      
  ! Flux source points
  !        
  DO i = 1,subdomain(isd)%nqnodes
    myWall = subdomain(isd)%qBCidList(i)
    r=r+1
    !
    ! Loop over all three small systems
    !
    DO een=1,neq             
      bc = eqn(een)%boundary(myWall)%known
      !
      ! Use this equation only, if flux is unknown
      !
      IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
        row=row+1      
        !
        !  For a single row in matrices
        !
        CALL formQrow(isd,een,r,sysMrow,rhsMrow)
        !
        ! Copy row to full matrix
        !
        DO col = 1,stk%neq
          if (sysMrow(col).NE.0.0_rk) sysMnnz = sysMnnz + 1
        END DO
        DO col = 1,stk%nb
          if (rhsMrow(col).NE.0.0_rk) rhsMnnz = rhsMnnz + 1 
        END DO              
      END IF
    END DO ! en        
  END DO ! q source point

  DEALLOCATE(sysMrow,rhsMrow)

end subroutine



!
!     ------------------------------------------------------------------
!
subroutine formStokesLeftRightHandSideVectors()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER en,mywall,isd,i,bc,isp,ispL

  !
  !     Set up lhs and rhs-vectors
  !
  !
  ! Loop over subdomains
  !
  isd=1
  !
  ! Function source points
  !
  DO i = 1,subdomain(isd)%nnodes
    isp = subdomain(isd)%nodeList(i)  ! source point
    myWall = subdomain(isd)%BCidList(i)

    DO en=1,neq
      ispL = isp + (en-1)*nnodes ! source point number in the large system       
      bc = eqn(en)%boundary(myWall)%known 
      IF (bc.EQ.iFunction) THEN
        stk%b(stk%col(ispL))=eqn(en)%u(isp)
      ELSE IF (bc.EQ.iFlux) THEN
        stk%x(stk%col(ispL))=eqn(en)%u(isp)
      ELSE IF (bc.EQ.iContact) THEN
        stk%x(stk%col(ispL))=eqn(en)%u(isp)
      END IF
    END DO ! en
  END DO ! u source point
  !
  ! Flux source points
  !
  DO i = 1,subdomain(isd)%nqnodes
    isp = subdomain(isd)%qnodeList(i)  ! source point
    myWall = subdomain(isd)%qBCidList(i)
    !
    ! Loop over all three small systems
    !
    DO en=1,neq
      ispL = isp + (en-1)*nqnodes! source point number in the large system 
      bc = eqn(en)%boundary(myWall)%known
      IF (bc.EQ.iFunction) THEN
        stk%x(stk%qcol(ispL))=eqn(en)%q(isp)
      ELSE IF (bc.EQ.iFlux) THEN
        stk%b(stk%qcol(ispL))=eqn(en)%q(isp)
      ELSE IF (bc.EQ.iContact) THEN
        stk%x(stk%qcol(ispL))=eqn(en)%q(isp)
      END IF
    END DO ! en
  END DO ! q source point
!      END DO ! subdomains

end subroutine
      



!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesDistributeUnknowns()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER en,mywall,isd,i,bc,isp,ispL
  !
  ! Loop over subdomains
  !
  isd=1
  !
  ! Function source points
  !
  DO i = 1,subdomain(isd)%nnodes
    isp = subdomain(isd)%nodeList(i)  ! source point
    myWall = subdomain(isd)%BCidList(i)

    DO en=1,neq
      ispL = isp + (en-1)*nnodes ! source point number in the large system       
      bc = eqn(en)%boundary(myWall)%known 
      IF (bc.EQ.iFunction) THEN
          eqn(en)%u(isp) = stk%b(stk%col(ispL))
      ELSE IF (bc.EQ.iFlux) THEN
          eqn(en)%u(isp) = stk%x(stk%col(ispL))
      ELSE IF (bc.EQ.iContact) THEN
          eqn(en)%u(isp) = stk%x(stk%col(ispL))
      END IF
    END DO ! en
  END DO ! u source point
  !
  ! Flux source points
  !
  DO i = 1,subdomain(isd)%nqnodes
    isp = subdomain(isd)%qnodeList(i)  ! source point
    myWall = subdomain(isd)%qBCidList(i)
    !
    ! Loop over all three small systems
    !
    DO en=1,neq
      ispL = isp + (en-1)*nqnodes ! source point number in the large system 
      bc = eqn(en)%boundary(myWall)%known

      IF (bc.EQ.iFunction) THEN
        eqn(en)%q(isp) = stk%x(stk%qcol(ispL))
      ELSE IF (bc.EQ.iFlux) THEN
        eqn(en)%q(isp) = stk%b(stk%qcol(ispL))
      ELSE IF (bc.EQ.iContact) THEN
        eqn(en)%q(isp) = stk%x(stk%qcol(ispL))
      END IF
    END DO ! en
  END DO ! q source point

end subroutine



!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesBigSystemSOLVE()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: sysM(:,:),rhsM(:,:),b(:)

  ALLOCATE(sysM(stk%neq,stk%neq))
  ALLOCATE(rhsM(stk%neq,stk%nb))
  !
  ! Form matrices
  !
  CALL WriteToLog("Stokes: forming system and rhs matrices.")
  call formStokesSysRhsMatrices(sysM,rhsM)
  !
  ! Free memory
  !
  do isd=1,nosd
    call dellocateStokesMatrices(isd)
  end do

  !
  ! Set up left and right hand side vectors
  !
  call formStokesLeftRightHandSideVectors()

  ALLOCATE (b(stk%neq))              
  b = MATMUL(rhsM,stk%b)
  DEALLOCATE(rhsM)

  !
  ! solve
  !

  ierr=0
  cput=0.0
  nits=0

  !
  ! solve full system of linear equations
  !
  stk%slv%type  = eqn(1)%slv%type ! 0 = direct solver, 1 = lsqr
  stk%slv%maxit = eqn(1)%slv%maxit
  stk%slv%stopt = eqn(1)%slv%stopt
  stk%slv%eps   = eqn(1)%slv%eps

  CALL WriteToLog("Stokes: solving.")
  IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver

    stk%slv%pret = 2 ! 2 = LU faktorizacija
    stk%slv%prep = 2 ! 2 = preconditioner calculated outside 
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL FormPREfm(stk%slv%pret,stk%neq,stk%Pivot,sysM,cput,ierr) 
    !
    ! Solve with direct solver
    !
    CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                   stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)

    DEALLOCATE(sysM)

  ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver

    CALL  TransformToCRS(sysM,stk%neq,stk%neq,stk%sysMcrs)
    DEALLOCATE(sysM)
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)
    !
    ! Solve with LSQR solver
    !
    CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
    stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 
    DEALLOCATE(stk%sysMcrs%v)


  END IF
  !
  !       Report to log file
  !
  WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
  CALL WriteToLog(parLogTekst)
  IF (stk%slv%maxit.LE.nits) THEN
    WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
    CALL WriteToLog(parLogTekst)
  END IF
  !
  !     Distribute SLE results to u and q vectors
  !
  call StokesDistributeUnknowns()
  !
  !     Free memory
  !      
  DEALLOCATE(b)

end

