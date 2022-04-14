!
! ------------------------------------------------------------------
!
SUBROUTINE StokesBigSystemSOLVElis()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  USE parallel
  USE plis
  IMPLICIT NONE

  INTEGER nits,ierr
  REAL ts,te

  CALL CPU_TIME(ts)
  !
  ! Set up left and right hand side vectors
  !
  !call formStokesLeftRightHandSideVectors()
  !
  ! Calculate b = rhsM * rhsVector
  !
  call plis_rhs_matvec(stk%sle,stk%b)
  !
  ! Solve with LIS solver
  !
  call plis_solve(stk%sle)
  ierr = plis_get_solver_status(stk%sle)
  nits = plis_get_solver_iterations(stk%sle)
  !
  ! Get result to all processes
  !
  call plis_getX(stk%sle,stk%x) 
  !
  ! Report to log file
  !
  WRITE (parLogTekst,'(A,A,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : LIS ",ierr,nits
  CALL WriteToLog(parLogTekst)
  IF (ierr.NE.0) THEN
    WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: solver error!"
    CALL WriteToLog(parLogTekst)
  END IF
  !
  !     Distribute SLE results to u and q vectors
  !
  call StokesDistributeUnknowns()
  !
  !     Free memory
  !      
  
  !
  !  Write to log
  ! 
  CALL CPU_TIME(te)
  WRITE (parLogTekst,'(A,F10.4)') "TIMER :: StokesBigSystemSOLVElis [s] = ",te-ts
  CALL WriteToLog(parLogTekst)

end


!
!     ------------------------------------------------------------------
!
SUBROUTINE DivideRows()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  USE parallel
  USE plis
  IMPLICIT NONE
    
  
  INTEGER row,mywall,isd,i,j,r,bc,een,rowSYS,irow

  ALLOCATE (stk%dn%rowTG(stk%sle%n))
  ALLOCATE (stk%dn%rowRHS(stk%sle%n))
  ALLOCATE (stk%dn%rowSYS(stk%sle%n))
  ALLOCATE (stk%dn%UQ(stk%sle%n))
  ALLOCATE (stk%dn%een(stk%sle%n))
  ALLOCATE (stk%dn%irow(stk%sle%n))

  !
  ! Loop over subdomains
  ! 
  isd = 1 ! single subdomain only  
  !
  ! Function source points
  !  
  r=0 ! row in T and G matrices
  j=0
  DO i = 1,subdomain(isd)%nnodes
    irow = subdomain(isd)%nodeList(i)
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
        row = stk%col(irow + (een-1)*nnodes)
        if (plis_is_row_mine(stk%sle,row)) then
          j=j+1
          rowSYS = ( row - stk%sle%is ) * stk%sle%neq - 1
          !   r = row in T & G matrix
          !   row = row in rhs matrix
          !   rowSYS = row in system matrix
          stk%dn%rowTG(j) = r
          stk%dn%rowRHS(j) = row
          stk%dn%rowSYS(j) = rowSYS
          stk%dn%een(j) = een
          stk%dn%irow(j) = irow
          stk%dn%UQ(j) = parU
        END IF
      END IF ! I need this equation
    END DO ! en
  END DO ! u source point
 
  !      
  !        
  ! Flux source points
  DO i = 1,subdomain(isd)%nqnodes
    irow = subdomain(isd)%qnodeList(i)
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
        row = stk%qcol(irow + (een-1)*nqnodes) 
        if (plis_is_row_mine(stk%sle,row)) then
          j=j+1
          rowSYS = ( row - stk%sle%is ) * stk%sle%neq - 1
          !   r = row in T & G matrix
          !   row = row in rhs matrix
          !   rowSYS = row in system matrix
          stk%dn%rowTG(j) = r
          stk%dn%rowRHS(j) = row
          stk%dn%rowSYS(j) = rowSYS
          stk%dn%een(j) = een
          stk%dn%irow(j) = irow
          stk%dn%UQ(j) = parQ

        END IF ! my row
      END IF ! need this equation
    END DO ! en        
  END DO ! q source point

END subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE verifyLISintegrals(integralsAvailable)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  USE parallel
  USE plis
  IMPLICIT NONE
      

  INTEGER, PARAMETER :: lun = 12
  CHARACTER(255) fileName
  LOGICAL integralsAvailable
  INTEGER myIntegralsAvailable,totIntegralsAvailable
  INTEGER ver(8)
  INTEGER ierr
    
  CALL WriteToLog("Looking for integrals files on hard disk!")
  WRITE(fileName,'(A,A,I0,A,I0)')  TRIM(parIntegralsFileName),".",env%myRank,"-",env%nProcs  
  OPEN (lun,FILE=TRIM(fileName),FORM='UNFORMATTED',STATUS='OLD', ERR = 10)
  READ(lun) ver(1),ver(2),ver(3),ver(4),ver(5),ver(6),ver(7),ver(8)
  IF (ver(1).EQ.parStokes.AND. &
      ver(2).EQ.env%myRank.AND. &
      ver(3).EQ.env%nProcs.AND. &
      ver(4).EQ.stk%neq.AND. &
      ver(5).EQ.stk%nb.AND. &
      ver(6).EQ.stk%sle%is.AND. &
      ver(7).EQ.stk%sle%ie.AND. &
      ver(8).EQ.stk%sle%n) THEN
        myIntegralsAvailable = 1
        GOTO 20      
  END IF  
10 CONTINUE ! error occured
  myIntegralsAvailable = 0
20 CONTINUE  
  CLOSE(lun)
  !
  !  Send info on integrals on all prcessors
  !
  CALL MPI_ALLREDUCE(myIntegralsAvailable, totIntegralsAvailable, 1, MPI_INTEGER, MPI_SUM, env%comm, ierr)
  IF (totIntegralsAvailable.EQ.env%nProcs) THEN 
    integralsAvailable = .TRUE.
    CALL WriteToLog("Integrals found!")
  ELSE
    integralsAvailable = .FALSE.
    CALL WriteToLog("Integrals not found!")
  END IF
END subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE stokesFormLISsysMrhsMcsrDiv()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  USE parallel
  USE plis
  IMPLICIT NONE
      
  INTEGER isd,i,r,c,col,myRow
  REAL(rk), ALLOCATABLE :: sysMrow(:),rhsMrow(:)
  REAL ts,te
  INTEGER nnz
  INTEGER dummy
  INTEGER, PARAMETER :: lun = 12
  CHARACTER(255) fileName
  LOGICAL integralsAvailable
    
  CALL CPU_TIME(ts)
  !
  ! Form matrices
  !
  CALL WriteToLog("Stokes: forming LIS system and rhs matrices in CSR format.")     
  !
  !  Do integrals exist on hard disk?
  ! 
  CALL verifyLISintegrals(integralsAvailable)  
  !
  !  Am I writing to disk?
  !
  IF ( (.NOT.integralsAvailable) .AND. (parWriteIntegrals.EQ.parYes) ) THEN
    CALL WriteToLog("Will write integrals to disk!")
    WRITE(fileName,'(A,A,I0,A,I0)')  TRIM(parIntegralsFileName),".",env%myRank,"-",env%nProcs

    OPEN (lun,FILE=TRIM(fileName),FORM='UNFORMATTED',STATUS='UNKNOWN')
    WRITE (lun) parStokes,env%myRank,env%nProcs,stk%neq,stk%nb,stk%sle%is,stk%sle%ie,stk%sle%n
  END IF
  !
  !  Reading from disk
  !
  IF (integralsAvailable) THEN
    CALL WriteToLog("Reading integrals from disk!")
    WRITE(fileName,'(A,A,I0,A,I0)')  TRIM(parIntegralsFileName),".",env%myRank,"-",env%nProcs
    OPEN (lun,FILE=TRIM(fileName),FORM='UNFORMATTED',STATUS='OLD')
    READ (lun) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
  ELSE
    CALL WriteToLog("Calculating integrals!")
  END IF

  !
  !  Form single rows for system and rhs matrices
  !
  ALLOCATE(sysMrow(stk%neq))
  ALLOCATE(rhsMrow(stk%nb))
   
  call plis_setup_bx(stk%sle)
  call plis_init_rhs_matrix(stk%sle,stk%nb) 

  nnz = stk%sle%n * stk%sle%neq
  
  ALLOCATE (  stk%sle%sm_val(0:nnz-1) )
  ALLOCATE ( stk%sle%sm_colI(0:nnz-1) )
  ALLOCATE ( stk%sle%sm_rowS(0:stk%sle%n  ) )

  ! matrix is dense, so I know rowS + colI in advance
  ! pazi max integer value !!!! INTEGER(4) max = 2147483647

  DO r = stk%sle%is,stk%sle%ie-1  !0,stk%sle%n-1  
    myRow = ( r - stk%sle%is ) * stk%sle%neq 
    DO c = 0,stk%sle%neq
      stk%sle%sm_colI(myRow+c) = c 
    END DO
    stk%sle%sm_rowS(r- stk%sle%is) = (r  - stk%sle%is )* stk%sle%neq
  END DO
  stk%sle%sm_rowS(stk%sle%ie-stk%sle%is) = nnz

  !
  ! Loop over subdomains
  ! 
  sysMrow = 0.0_rk
  rhsMrow = 0.0_rk

  isd = 1 ! single subdomain only  
  
  DO i = 1,stk%sle%n
    IF (integralsAvailable) THEN
      ! Read
      CALL rdvec(lun,stk%neq,sysMrow)
      CALL rdvec(lun,stk%nb,rhsMrow)      
    ELSE    
      ! Calculate  
      IF (stk%dn%UQ(i) .EQ. parU) THEN
        !
        ! Function source points
        !
        CALL formUrowFly(isd,stk%dn%een(i),stk%dn%irow(i),sysMrow,rhsMrow)
      ELSE 
        !        
        ! Flux source points
        !
        CALL formQrowFly(isd,stk%dn%een(i),stk%dn%irow(i),sysMrow,rhsMrow)
      END IF
      IF (parWriteIntegrals.EQ.parYes) THEN
        CALL wrvec(lun,stk%neq,sysMrow)
        CALL wrvec(lun,stk%nb,rhsMrow)
      END IF
    END IF
    !
    ! Copy row to system matrix
    !
    DO col = 1,stk%neq
      stk%sle%sm_val(stk%dn%rowSYS(i)  + col ) = sysMrow(col)
    END DO
    !
    ! Copy row to rhs matrix
    !
    DO col = 1,stk%nb
      call plis_set_RHSmatrix_element(stk%sle,stk%dn%rowRHS(i),col,rhsMrow(col))
    END DO 
  END DO
  !
  ! Free memory
  !
  DEALLOCATE(sysMrow,rhsMrow)
  !
  !  Assemble CRS structure of the system matrix
  !
  CALL plis_assemble_system_matrix_csr(nnz,stk%sle)  
  !
  !  Close integrals file
  !
  IF ( (integralsAvailable).OR.( (.NOT.integralsAvailable) .AND. (parWriteIntegrals.EQ.parYes) ) ) THEN
    CALL WriteToLog("Closing integrals files!")
    CLOSE(lun)
  END IF  
  !
  !  Write to log
  ! 
  CALL CPU_TIME(te)
  WRITE (parLogTekst,'(A,F10.4)') "TIMER :: stokesFormLISsysMrhsMcsrDiv [s] = ",te-ts
  CALL WriteToLog(parLogTekst)

END subroutine




!
!     ------------------------------------------------------------------
!
SUBROUTINE stokesFormLISsysMrhsMcsr()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  USE parallel
  USE plis
  IMPLICIT NONE
      
  INTEGER row,mywall,isd,i,r,c,col,bc,een,myRow
  REAL(rk), ALLOCATABLE :: sysMrow(:),rhsMrow(:)
  REAL ts,te
  INTEGER nnz
    
  CALL CPU_TIME(ts)
  !
  ! Form matrices
  !
  CALL WriteToLog("Stokes: forming LIS system and rhs matrices in CSR format.")      
  !
  !  Form single rows for system and rhs matrices
  !
  ALLOCATE(sysMrow(stk%neq))
  ALLOCATE(rhsMrow(stk%nb))
   
  call plis_setup_bx(stk%sle)
  call plis_init_rhs_matrix(stk%sle,stk%nb) 

  nnz = stk%sle%n * stk%sle%neq
  
  ALLOCATE (  stk%sle%sm_val(0:nnz-1) )
  ALLOCATE ( stk%sle%sm_colI(0:nnz-1) )
  ALLOCATE ( stk%sle%sm_rowS(0:stk%sle%n  ) )

  ! matrix is dense, so I know rowS + colI in advance
  ! pazi max integer value !!!! INTEGER(4) max = 2147483647
  DO r = stk%sle%is,stk%sle%ie-1  !0,stk%sle%n-1  
    myRow = ( r - stk%sle%is ) * stk%sle%neq 
    DO c = 0,stk%sle%neq
      stk%sle%sm_colI(myRow+c) = c 
    END DO
    stk%sle%sm_rowS(r- stk%sle%is) = (r  - stk%sle%is )* stk%sle%neq
  END DO
  stk%sle%sm_rowS(stk%sle%ie-stk%sle%is) = nnz

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
        row = stk%col(subdomain(isd)%nodeList(i) + (een-1)*nnodes)
        if (plis_is_row_mine(stk%sle,row)) then
          !
          !  For a single row in matrices
          !
          CALL formUrow(isd,een,r,sysMrow,rhsMrow)
          !
          ! Copy row to system matrix
          !
          myRow = ( row - stk%sle%is ) * stk%sle%neq - 1
          DO col = 1,stk%neq
            stk%sle%sm_val(myRow  + col ) = sysMrow(col)
          END DO
          !
          ! Copy row to rhs matrix
          !
          DO col = 1,stk%nb
            call plis_set_RHSmatrix_element(stk%sle,row,col,rhsMrow(col))
          END DO            
        END IF
      END IF ! I need this equation
    END DO ! en
  END DO ! u source point
 
  !      
  !        
  ! Flux source points
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
        row = stk%qcol(subdomain(isd)%qnodeList(i) + (een-1)*nqnodes) 
        if (plis_is_row_mine(stk%sle,row)) then
          !
          !  For a single row in matrices
          !
          CALL formQrow(isd,een,r,sysMrow,rhsMrow)
          !
          ! Copy row to system matrix
          !
          myRow = ( row - stk%sle%is ) * stk%sle%neq - 1
          DO col = 1,stk%neq
            stk%sle%sm_val(myRow  + col ) = sysMrow(col)
          END DO
          !
          ! Copy row to rhs matrix
          !
          DO col = 1,stk%nb
            call plis_set_RHSmatrix_element(stk%sle,row,col,rhsMrow(col))
          END DO                        
        END IF ! my row
      END IF ! need this equation
    END DO ! en        
  END DO ! q source point
 
  DEALLOCATE(sysMrow,rhsMrow)
  !
  ! Free memory
  !
  do isd=1,nosd
    call dellocateStokesMatrices(isd)
  end do
  !
  !  Assemble CRS structure of the system matrix
  !
  CALL plis_assemble_system_matrix_csr(nnz,stk%sle)  

  !
  !  Write to log
  ! 
  CALL CPU_TIME(te)
  WRITE (parLogTekst,'(A,F10.4)') "TIMER :: stokesFormLISsysMrhsMcsr [s] = ",te-ts
  CALL WriteToLog(parLogTekst)

END subroutine



!
!     ------------------------------------------------------------------
!
SUBROUTINE stokesFormLISsysMrhsMveryslow()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  USE parallel
  USE plis
  IMPLICIT NONE
      
  INTEGER row,mywall,isd,i,r,col,bc,een
  REAL(rk), ALLOCATABLE :: sysMrow(:),rhsMrow(:)   
  !
  ! Form matrices
  !
  CALL WriteToLog("Stokes: forming LIS system and rhs matrices.")      
  !
  !  Form single rows for system and rhs matrices
  !
  ALLOCATE(sysMrow(stk%neq))
  ALLOCATE(rhsMrow(stk%nb))

   
  call plis_setup_Abx(stk%sle,stk%neq)
  call plis_init_rhs_matrix(stk%sle,stk%nb)

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
        row = stk%col(subdomain(isd)%nodeList(i) + (een-1)*nnodes)
        if (plis_is_row_mine(stk%sle,row)) then
          !
          !  For a single row in matrices
          !
          CALL formUrow(isd,een,r,sysMrow,rhsMrow)
          !
          ! Copy row to system matrix
          !
          DO col = 1,stk%neq
            !if (sysMrow(col).NE.0.0_rk) 
            call plis_set_matrix_element(stk%sle,row,col,sysMrow(col))
          END DO
          !
          ! Copy row to rhs matrix
          !
          DO col = 1,stk%nb
            call plis_set_RHSmatrix_element(stk%sle,row,col,rhsMrow(col))
          END DO            
        END IF
      END IF ! I need this equation
    END DO ! en
  END DO ! u source point
 
  !      
  !        
  ! Flux source points
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
        row = stk%qcol(subdomain(isd)%qnodeList(i) + (een-1)*nqnodes) 
        if (plis_is_row_mine(stk%sle,row)) then
          !
          !  For a single row in matrices
          !
          CALL formQrow(isd,een,r,sysMrow,rhsMrow)
          !
          ! Copy row to system matrix
          !
          DO col = 1,stk%neq
            !if (sysMrow(col).NE.0.0_rk) 
            call plis_set_matrix_element(stk%sle,row,col,sysMrow(col))
          END DO
          !
          ! Copy row to rhs matrix
          !
          DO col = 1,stk%nb
            call plis_set_RHSmatrix_element(stk%sle,row,col,rhsMrow(col))
          END DO                        
        END IF ! my row
      END IF ! need this equation
    END DO ! en        
  END DO ! q source point
 
  DEALLOCATE(sysMrow,rhsMrow)

  !
  ! Free memory
  !
  do isd=1,nosd
    call dellocateStokesMatrices(isd)
  end do
  !
  !  Assemble CRS structure of the system matrix
  !
  call plis_assemble_system_matrix(stk%sle) 

END subroutine





!
!     ------------------------------------------------------------------
!
SUBROUTINE formUrowFly(isd,een,irow,sysMrow,rhsMrow)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE


  REAL(rk), ALLOCATABLE :: rowTxx(:), rowTxy(:), rowTxz(:), rowTyy(:), rowTyz(:), rowTzz(:)
  REAL(rk), ALLOCATABLE :: rowGxx(:),rowGxy(:),rowGxz(:),rowGyy(:),rowGyz(:),rowGzz(:)

  REAL(rk), ALLOCATABLE :: rowprTx(:), rowprTy(:), rowprTz(:)
  REAL(rk), ALLOCATABLE :: rowprGx(:),rowprGy(:),rowprGz(:)    
  
  REAL(rk) sysMrow(stk%neq),rhsMrow(stk%nb)
  integer j,k,isd,isp,mywallC,en,ispL,col,bc,een,irow
  real(rk) val,multi,err

!
!   Allocate space for integrals for a single source point
!
  
  ALLOCATE (rowTxx(nnodes), rowTxy(nnodes), rowTxz(nnodes), rowTyy(nnodes), rowTyz(nnodes), rowTzz(nnodes))
  ALLOCATE (rowGxx(nqnodes),rowGxy(nqnodes),rowGxz(nqnodes),rowGyy(nqnodes),rowGyz(nqnodes),rowGzz(nqnodes))

  ALLOCATE (rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes))
  ALLOCATE (rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes))  

  isd = 1

  CALL sdIntRowUStokes(isd,irow, & ! irow = source point ID, samo po subdomainu
        rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz, &
        rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz, &
        rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz, &
        err)

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
      IF (een.EQ.1.AND.en.EQ.1) val = rowTxx(isp)
      IF (een.EQ.1.AND.en.EQ.2) val = rowTxy(isp)
      IF (een.EQ.1.AND.en.EQ.3) val = rowTxz(isp)
      IF (een.EQ.2.AND.en.EQ.1) val = rowTxy(isp)
      IF (een.EQ.2.AND.en.EQ.2) val = rowTyy(isp)
      IF (een.EQ.2.AND.en.EQ.3) val = rowTyz(isp)
      IF (een.EQ.3.AND.en.EQ.1) val = rowTxz(isp)
      IF (een.EQ.3.AND.en.EQ.2) val = rowTyz(isp)
      IF (een.EQ.3.AND.en.EQ.3) val = rowTzz(isp)
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

      IF (een.EQ.1.AND.en.EQ.1) val = rowGxx(isp)
      IF (een.EQ.1.AND.en.EQ.2) val = rowGxy(isp)
      IF (een.EQ.1.AND.en.EQ.3) val = rowGxz(isp)
      IF (een.EQ.2.AND.en.EQ.1) val = rowGxy(isp)
      IF (een.EQ.2.AND.en.EQ.2) val = rowGyy(isp)
      IF (een.EQ.2.AND.en.EQ.3) val = rowGyz(isp)
      IF (een.EQ.3.AND.en.EQ.1) val = rowGxz(isp)
      IF (een.EQ.3.AND.en.EQ.2) val = rowGyz(isp)
      IF (een.EQ.3.AND.en.EQ.3) val = rowGzz(isp)
 
      val = multi / subdomain(isd)%diff * val
      IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
        sysMrow(col) = - val  ! XXXXXX
      END IF
      IF (bc.EQ.iFlux) THEN
        rhsMrow(col) =   val
      END IF
    END DO
  END DO

  DEALLOCATE (rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz)
  DEALLOCATE (rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz)
  DEALLOCATE (rowprTx,rowprTy,rowprTz)
  DEALLOCATE (rowprGx,rowprGy,rowprGz)

end subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE formQrowFly(isd,een,irow,sysMrow,rhsMrow)

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE
  
  REAL(rk), ALLOCATABLE :: rowTxx(:), rowTxy(:), rowTxz(:), rowTyy(:), rowTyz(:), rowTzz(:)
  REAL(rk), ALLOCATABLE :: rowGxx(:),rowGxy(:),rowGxz(:),rowGyy(:),rowGyz(:),rowGzz(:)

  REAL(rk), ALLOCATABLE :: rowprTx(:), rowprTy(:), rowprTz(:)
  REAL(rk), ALLOCATABLE :: rowprGx(:),rowprGy(:),rowprGz(:)    
  
  REAL(rk) sysMrow(stk%neq),rhsMrow(stk%nb)
  integer j,k,isd,isp,mywallC,en,ispL,col,bc,een,irow
  real(rk) val,multi,err

!
!   Allocate space for integrals for a single source point
!
  
  ALLOCATE (rowTxx(nnodes), rowTxy(nnodes), rowTxz(nnodes), rowTyy(nnodes), rowTyz(nnodes), rowTzz(nnodes))
  ALLOCATE (rowGxx(nqnodes),rowGxy(nqnodes),rowGxz(nqnodes),rowGyy(nqnodes),rowGyz(nqnodes),rowGzz(nqnodes))

  ALLOCATE (rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes))
  ALLOCATE (rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes))  

  isd = 1

  CALL sdIntRowQStokes(isd,irow, & ! irow = source point ID, samo po subdomainu
  rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz, &
  rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz, &
  rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz, &
  err)     

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
      
      IF (een.EQ.1.AND.en.EQ.1) val = rowTxx(isp)
      IF (een.EQ.1.AND.en.EQ.2) val = rowTxy(isp)
      IF (een.EQ.1.AND.en.EQ.3) val = rowTxz(isp)
      IF (een.EQ.2.AND.en.EQ.1) val = rowTxy(isp)
      IF (een.EQ.2.AND.en.EQ.2) val = rowTyy(isp)
      IF (een.EQ.2.AND.en.EQ.3) val = rowTyz(isp)
      IF (een.EQ.3.AND.en.EQ.1) val = rowTxz(isp)
      IF (een.EQ.3.AND.en.EQ.2) val = rowTyz(isp)
      IF (een.EQ.3.AND.en.EQ.3) val = rowTzz(isp)

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

      IF (een.EQ.1.AND.en.EQ.1) val = rowGxx(isp)
      IF (een.EQ.1.AND.en.EQ.2) val = rowGxy(isp)
      IF (een.EQ.1.AND.en.EQ.3) val = rowGxz(isp)
      IF (een.EQ.2.AND.en.EQ.1) val = rowGxy(isp)
      IF (een.EQ.2.AND.en.EQ.2) val = rowGyy(isp)
      IF (een.EQ.2.AND.en.EQ.3) val = rowGyz(isp)
      IF (een.EQ.3.AND.en.EQ.1) val = rowGxz(isp)
      IF (een.EQ.3.AND.en.EQ.2) val = rowGyz(isp)
      IF (een.EQ.3.AND.en.EQ.3) val = rowGzz(isp)
 
      val = multi / subdomain(isd)%diff * val

      IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
        sysMrow(col) = - val  ! XXXXXX
      END IF

      IF (bc.EQ.iFlux) THEN
        rhsMrow(col) =  val
      END IF

    END DO
  END DO  


  DEALLOCATE (rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz)
  DEALLOCATE (rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz)
  DEALLOCATE (rowprTx,rowprTy,rowprTz)
  DEALLOCATE (rowprGx,rowprGy,rowprGz)  

end subroutine

