

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesInlet3sidesTorquePiCRS()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: b(:)

  integer i,j,k,nid,outside,particle,lun
  real(rk) flowVelocity(3)

  real(rk) c,rho,nu,resT(3),cdre(3),rotT(3),cmre(3)
  REAL(rk), ALLOCATABLE :: Fnum(:,:),Tnum(:,:),Fdrag(:,:),Flift(:,:)

  rho=1.0
  nu=1.0
  call getSupElc(parSuelLam1,parSuelLam2,parSuelE1,parSuelE2,c)
  !WRITE (parLogTekst,'(A,5G20.10)') " SEL (a,b,c,e1,e2) : ",parSuelLam1*c,parSuelLam2*c,c,parSuelE1,parSuelE2
  !CALL WriteToLog(parLogTekst)

  ALLOCATE (Fnum(3,3),Tnum(3,3))
  ALLOCATE (Fdrag(3,3),Flift(3,3))
  Fnum=0.0_rk
  Tnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk


  CALL WriteToLog("Stokes: solving.")


  ierr=0
  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1


  !
  ! TORQUE PI tensor part
  !
  

  do i=1,3
    if (i.eq.1) CALL WriteToLog("Torque PI tensor: solving for X")
    if (i.eq.2) CALL WriteToLog("Torque PI tensor: solving for Y")
    if (i.eq.3) CALL WriteToLog("Torque PI tensor: solving for Z")

    !
    ! Get flow velocity
    !
    flowVelocity(i)=0.0_rk

    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   

            if (i.eq.1) then ! 
              flowVelocity(2) =   node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(2)
            end if

            if (i.eq.2) then ! 
              flowVelocity(1) =   node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(1)
            end if

            if (i.eq.3) then ! 
              flowVelocity(2) =   node(nid)%x(1)
              flowVelocity(1) =   node(nid)%x(2)
            end if


            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)                
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

  !
  ! Calculate b = rhsM * rhsVector
  !
    ALLOCATE (b(stk%neq))              
    call CRSxV(stk%rhsMcrs,stk%b,stk%nb,b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 

        
    !
    !  Report to log file
    !
    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
    CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    ! force coeficient
    cdre(i)=4.0_rk*resT(i)
    !
    !  Calculate torque numerically
    !
    call getTorqueIntegrateFluxes(particle,Tnum(i,:))
    !Tnum(i,i) = - Tnum(i,i) ! MINUS !!!
    ! Calculate rotation tensor diagonal components
    rotT(i)=Tnum(i,i)/(pi*c*c*c*rho*nu)
    ! moment coeficient
    cmre(i)=4.0_rk*rotT(i)

    WRITE (parLogTekst,'(3G20.10)') Tnum(i,i),rotT(i),cmre(i)
    CALL WriteToLog(parLogTekst)


  end do


  lun=58
  OPEN (unit=lun,file=TRIM("and.piT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 PIxx PIyy PIzz cmreX cmreY cmreZ Tx Ty Tz"
  WRITE (lun,'(16G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                       rotT(1),rotT(2),rotT(3),cmre(1),cmre(2),cmre(3),Tnum(1,1),Tnum(2,2),Tnum(3,3)
  CLOSE (lun)

end subroutine




!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesInlet3sidesTorqueOmegaCRS()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: b(:)

  integer i,j,k,nid,outside,particle,lun
  real(rk) flowVelocity(3)

  real(rk) c,rho,nu,resT(3),cdre(3),rotT(3),cmre(3)
  REAL(rk), ALLOCATABLE :: Fnum(:,:),Tnum(:,:),Fdrag(:,:),Flift(:,:)

  rho=1.0
  nu=1.0
  call getSupElc(parSuelLam1,parSuelLam2,parSuelE1,parSuelE2,c)
  !WRITE (parLogTekst,'(A,5G20.10)') " SEL (a,b,c,e1,e2) : ",parSuelLam1*c,parSuelLam2*c,c,parSuelE1,parSuelE2
  !CALL WriteToLog(parLogTekst)

  ALLOCATE (Fnum(3,3),Tnum(3,3))
  ALLOCATE (Fdrag(3,3),Flift(3,3))
  Fnum=0.0_rk
  Tnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk


  CALL WriteToLog("Stokes: solving.")


  ierr=0
  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1


  !
  ! TORQUE
  !

  do i=1,3
    if (i.eq.1) CALL WriteToLog("Torque: solving for X")
    if (i.eq.2) CALL WriteToLog("Torque: solving for Y")
    if (i.eq.3) CALL WriteToLog("Torque: solving for Z")

    !
    ! Get flow velocity
    !
    flowVelocity(i)=0.0_rk

    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   

            if (i.eq.1) then ! rotation around x axis
              flowVelocity(2) = - node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(2)
            end if

            if (i.eq.2) then ! rotation around y axis
              flowVelocity(1) =   node(nid)%x(3)
              flowVelocity(3) = - node(nid)%x(1)
            end if

            if (i.eq.3) then ! rotation around z axis ! PAZI NAROBE PREDZNAK !!! (-y,x,0) bi moralo biti
              flowVelocity(1) =   node(nid)%x(2)
              flowVelocity(2) = - node(nid)%x(1)              
            end if

            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)                
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

  !
  ! Calculate b = rhsM * rhsVector
  !
    ALLOCATE (b(stk%neq))              
    call CRSxV(stk%rhsMcrs,stk%b,stk%nb,b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 

        
    !
    !  Report to log file
    !
    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
    CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    ! force coeficient
    cdre(i)=4.0_rk*resT(i)
    !
    !  Calculate torque numerically
    !
    call getTorqueIntegrateFluxes(particle,Tnum(i,:))
    !Tnum(i,i) = - Tnum(i,i) ! MINUS !!!
    ! Calculate rotation tensor diagonal components
    rotT(i)=Tnum(i,i)/(pi*c*c*c*rho*nu) 
    ! moment coeficient
    cmre(i)=4.0_rk*rotT(i)

    WRITE (parLogTekst,'(3G20.10)') Tnum(i,i),rotT(i),cmre(i)
    CALL WriteToLog(parLogTekst)


  end do


  lun=58
  OPEN (unit=lun,file=TRIM("and.rotT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 Omegaxx Omegayy Omegazz cmreX cmreY cmreZ Tx Ty Tz"
  WRITE (lun,'(16G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                       rotT(1),rotT(2),rotT(3),cmre(1),cmre(2),cmre(3),Tnum(1,1),Tnum(2,2),Tnum(3,3)
  CLOSE (lun)

end subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesInlet3sidesForceCRS()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: b(:)

  integer i,j,k,nid,outside,particle,lun
  real(rk) flowVelocity(3)

  real(rk) c,rho,nu,resT(3),cdre(3)
  REAL(rk), ALLOCATABLE :: Fnum(:,:),Tnum(:,:),Fdrag(:,:),Flift(:,:)

  rho=1.0
  nu=1.0
  call getSupElc(parSuelLam1,parSuelLam2,parSuelE1,parSuelE2,c)
  WRITE (parLogTekst,'(A,5G20.10)') " SEL (a,b,c,e1,e2) : ",parSuelLam1*c,parSuelLam2*c,c,parSuelE1,parSuelE2
  CALL WriteToLog(parLogTekst)

  ALLOCATE (Fnum(3,3),Tnum(3,3))
  ALLOCATE (Fdrag(3,3),Flift(3,3))
  Fnum=0.0_rk
  Tnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk


  CALL WriteToLog("Stokes: solving.")


  ierr=0
  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1

  !
  ! FORCE
  !

  do i=1,3
    if (i.eq.1) CALL WriteToLog("Force: solving for X")
    if (i.eq.2) CALL WriteToLog("Force: solving for Y")
    if (i.eq.3) CALL WriteToLog("Force: solving for Z")
    !
    ! Get flow velocity
    !
    flowVelocity = 0.0_rk
    flowVelocity(i)=1.0_rk
    call normalize(flowVelocity)
    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   
            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

  !
  ! Calculate b = rhsM * rhsVector
  !
    ALLOCATE (b(stk%neq))              
    call CRSxV(stk%rhsMcrs,stk%b,stk%nb,b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 
        
    !
    !  Report to log file
    !
    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
    CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    !
    ! Calculate drag and lift from total force
    !
    call getDragLift(flowVelocity,Fnum(i,:),Fdrag(i,:),Flift(i,:))
    ! 
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    cdre(i)=8.0_rk*c*resT(i)

    WRITE (parLogTekst,'(3G20.10)') Fnum(i,i),resT(i),cdre(i)
    CALL WriteToLog(parLogTekst)

  end do

  lun=58
  OPEN (unit=lun,file=TRIM("and.resT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 Kxx Kyy Kzz cdreX cdreY cdreZ Fx Fy Fz"
  WRITE (lun,'(19G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                           resT(1),resT(2),resT(3),cdre(1),cdre(2),cdre(3),Fnum(1,1),Fnum(2,2),Fnum(3,3)
  CLOSE (lun)

end subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesInlet3sidesForceTorque()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: sysM(:,:),rhsM(:,:),b(:)

  integer i,j,k,nid,outside,particle,lun
  real(rk) flowVelocity(3)

  real(rk) c,rho,nu,resT(3),cdre(3),rotT(3),cmre(3)
  REAL(rk), ALLOCATABLE :: Fnum(:,:),Tnum(:,:),Fdrag(:,:),Flift(:,:)

  rho=1.0
  nu=1.0
  call getSupElc(parSuelLam1,parSuelLam2,parSuelE1,parSuelE2,c)
  WRITE (parLogTekst,'(A,5G20.10)') " SEL (a,b,c,e1,e2) : ",parSuelLam1*c,parSuelLam2*c,c,parSuelE1,parSuelE2
  CALL WriteToLog(parLogTekst)

  ALLOCATE(sysM(stk%neq,stk%neq))
  ALLOCATE(rhsM(stk%neq,stk%nb))

  ALLOCATE (Fnum(3,3),Tnum(3,3))
  ALLOCATE (Fdrag(3,3),Flift(3,3))
  Fnum=0.0_rk
  Tnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk

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


  CALL WriteToLog("Stokes: solving.")
  !
  ! Calculate preconditioner
  !          
  ierr=0
  !
  ! solve full system of linear equations
  !
  stk%slv%type  = eqn(1)%slv%type ! 0 = direct solver, 1 = lsqr
  stk%slv%maxit = eqn(1)%slv%maxit
  stk%slv%stopt = eqn(1)%slv%stopt
  stk%slv%eps   = eqn(1)%slv%eps 

  IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver

    stk%slv%pret = 2 ! 2 = LU faktorizacija
    stk%slv%prep = 2 ! 2 = preconditioner calculated outside 
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL FormPREfm(stk%slv%pret,stk%neq,stk%Pivot,sysM,cput,ierr)   
    DEALLOCATE(sysM)

  ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver

    CALL  TransformToCRS(sysM,stk%neq,stk%neq,stk%sysMcrs)
    DEALLOCATE(sysM)
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)
  END IF

  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1

  !
  ! FORCE
  !

  do i=1,3
    if (i.eq.1) CALL WriteToLog("Force: solving for X")
    if (i.eq.2) CALL WriteToLog("Force: solving for Y")
    if (i.eq.3) CALL WriteToLog("Force: solving for Z")
    !
    ! Get flow velocity
    !
    flowVelocity = 0.0_rk
    flowVelocity(i)=1.0_rk
    call normalize(flowVelocity)
    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   
            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

    ALLOCATE (b(stk%neq))              
    b = MATMUL(rhsM,stk%b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
      !
      ! Solve with direct solver
      !
      CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                     stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
            
    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 
        
    END IF
    !
    !  Report to log file
    !
    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
    CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    !
    ! Calculate drag and lift from total force
    !
    call getDragLift(flowVelocity,Fnum(i,:),Fdrag(i,:),Flift(i,:))
    ! 
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    cdre(i)=8.0_rk*c*resT(i)

    WRITE (parLogTekst,'(3G20.10)') Fnum(i,i),resT(i),cdre(i)
    CALL WriteToLog(parLogTekst)

  end do

  lun=58
  OPEN (unit=lun,file=TRIM("and.resT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 Kxx Kyy Kzz cdreX cdreY cdreZ Fx Fy Fz"
  WRITE (lun,'(19G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                           resT(1),resT(2),resT(3),cdre(1),cdre(2),cdre(3),Fnum(1,1),Fnum(2,2),Fnum(3,3)
  CLOSE (lun)

  !
  ! TORQUE
  !

  do i=1,3
    if (i.eq.1) CALL WriteToLog("Torque: solving for X")
    if (i.eq.2) CALL WriteToLog("Torque: solving for Y")
    if (i.eq.3) CALL WriteToLog("Torque: solving for Z")

    !
    ! Get flow velocity
    !
    flowVelocity(i)=0.0_rk

    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   

            if (i.eq.1) then ! rotation around x axis
              flowVelocity(2) = - node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(2)
            end if

            if (i.eq.2) then ! rotation around y axis
              flowVelocity(1) =   node(nid)%x(3)
              flowVelocity(3) = - node(nid)%x(1)
            end if

            if (i.eq.3) then ! rotation around z axis
              flowVelocity(1) =   node(nid)%x(2)
              flowVelocity(2) = - node(nid)%x(1)              
            end if

            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)                
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

    ALLOCATE (b(stk%neq))              
    b = MATMUL(rhsM,stk%b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
      !
      ! Solve with direct solver
      !
      CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                     stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
            
    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 

        
    END IF
    !
    !  Report to log file
    !
    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
    CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    ! force coeficient
    cdre(i)=4.0_rk*resT(i)
    !
    !  Calculate torque numerically
    !
    call getTorqueIntegrateFluxes(particle,Tnum(i,:))
    !Tnum(i,i) = - Tnum(i,i) ! MINUS !!!
    ! Calculate rotation tensor diagonal components
    rotT(i)=Tnum(i,i)/(pi*c*c*c*rho*nu) 
    ! moment coeficient
    cmre(i)=4.0_rk*rotT(i)

    WRITE (parLogTekst,'(3G20.10)') Tnum(i,i),rotT(i),cmre(i)
    CALL WriteToLog(parLogTekst)


  end do


  lun=58
  OPEN (unit=lun,file=TRIM("and.rotT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 Omegaxx Omegayy Omegazz cmreX cmreY cmreZ Tx Ty Tz"
  WRITE (lun,'(16G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                       rotT(1),rotT(2),rotT(3),cmre(1),cmre(2),cmre(3),Tnum(1,1),Tnum(2,2),Tnum(3,3)
  CLOSE (lun)




  !
  ! TORQUE PI tensor part
  !
  

  do i=1,3
    if (i.eq.1) CALL WriteToLog("Torque PI tensor: solving for X")
    if (i.eq.2) CALL WriteToLog("Torque PI tensor: solving for Y")
    if (i.eq.3) CALL WriteToLog("Torque PI tensor: solving for Z")

    !
    ! Get flow velocity
    !
    flowVelocity(i)=0.0_rk

    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   

            if (i.eq.1) then ! 
              flowVelocity(2) =   node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(2)
            end if

            if (i.eq.2) then ! 
              flowVelocity(1) =   node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(1)
            end if

            if (i.eq.3) then ! 
              flowVelocity(2) =   node(nid)%x(1)
              flowVelocity(1) =   node(nid)%x(2)
            end if


            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)                
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

    ALLOCATE (b(stk%neq))              
    b = MATMUL(rhsM,stk%b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
      !
      ! Solve with direct solver
      !
      CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                     stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
            
    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 

        
    END IF
    !
    !  Report to log file
    !
    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
    CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    ! force coeficient
    cdre(i)=4.0_rk*resT(i)
    !
    !  Calculate torque numerically
    !
    call getTorqueIntegrateFluxes(particle,Tnum(i,:))
    !Tnum(i,i) = - Tnum(i,i) ! MINUS !!!
    ! Calculate rotation tensor diagonal components
    rotT(i)=Tnum(i,i)/(pi*c*c*c*rho*nu)
    ! moment coeficient
    cmre(i)=4.0_rk*rotT(i)

    WRITE (parLogTekst,'(3G20.10)') Tnum(i,i),rotT(i),cmre(i)
    CALL WriteToLog(parLogTekst)


  end do


  lun=58
  OPEN (unit=lun,file=TRIM("and.piT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 PIxx PIyy PIzz cmreX cmreY cmreZ Tx Ty Tz"
  WRITE (lun,'(16G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                       rotT(1),rotT(2),rotT(3),cmre(1),cmre(2),cmre(3),Tnum(1,1),Tnum(2,2),Tnum(3,3)
  CLOSE (lun)




!  !
!  ! TEST DELETE FOR PRODUCTION RUN
!  !
!  
!  do i=1,3
!    if (i.eq.1) CALL WriteToLog("TEST RUN 1")
!    if (i.eq.2) CALL WriteToLog("TEST RUN 2")
!    if (i.eq.3) CALL WriteToLog("TEST RUN 3")
!
!    !
!    ! Get flow velocity
!    !
!    flowVelocity(i)=0.0_rk
!
!    !
!    ! Apply flow velocity as boundary condition
!    !
!    do k=1,subdomain(isd)%nnodes
!        if (subdomain(isd)%BCidList(k).eq.outside) then
!            nid = subdomain(isd)%nodeList(k)   
!
!            if (i.eq.1) then ! 
!              flowVelocity(1) =   node(nid)%x(2)
!              flowVelocity(2) =   node(nid)%x(3)
!            end if
!
!            if (i.eq.2) then ! 
!              flowVelocity(1) =   node(nid)%x(2)
!              flowVelocity(2) =   3.0_rk * node(nid)%x(3)
!              flowVelocity(3) =   5.0_rk * node(nid)%x(1)
!            end if
!
!
!            !
!            if (i.eq.3) then ! 
!              flowVelocity(1) =   1.0_rk * node(nid)%x(3)
!              flowVelocity(2) =   1.3_rk * node(nid)%x(1)
!              flowVelocity(3) =   2.7_rk * node(nid)%x(2)
!            end if
!
!
!            DO j=1,3
!                eqn(j)%u(nid) = flowVelocity(j)                
!            END DO
!        end if
!    end do
!    !
!    ! Set up left and right hand side vectors
!    !
!    call formStokesLeftRightHandSideVectors()
!
!    ALLOCATE (b(stk%neq))              
!    b = MATMUL(rhsM,stk%b)
! 
!    !
!    ! solve
!    !
!  
!    ierr=0
!    cput=0.0_rk
!    nits=0
!              
!    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
!      !
!      ! Solve with direct solver
!      !
!      CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
!                     stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
!            
!    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
!      !
!      ! Solve with LSQR solver
!      !
!      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
!      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 
!
!        
!    END IF
!    !
!    !  Report to log file
!    !
!    WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
!    CALL WriteToLog(parLogTekst)
!    IF (stk%slv%maxit.LE.nits) THEN
!      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
!      CALL WriteToLog(parLogTekst)
!    END IF
!    !
!    !  Distribute SLE results to u and q vectors
!    !
!    call StokesDistributeUnknowns()
!    !
!    !  Free memory
!    !      
!    DEALLOCATE(b)
!    !
!    !  Calculate force numerically
!    !
!    !call getForceIntegrateFluxes(particle,Fnum(i,:))
!    ! Calculate resistance tensor diagonal components
!    !resT(i)=Fnum(i,i)/(pi*c*rho*nu)
!    ! force coeficient
!    !cdre(i)=4.0_rk*resT(i)
!    !
!    !  Calculate torque numerically
!    !
!    call getTorqueIntegrateFluxes(particle,Tnum(i,:))
!    !Tnum(i,1) = - Tnum(i,1) ! MINUS !!!
!    !Tnum(i,2) = - Tnum(i,2) ! MINUS !!!
!    !Tnum(i,3) = - Tnum(i,3) ! MINUS !!!
!    ! Calculate total torque PI + OMEGA parts
!    rotT(1)=Tnum(i,1)/(pi*c*c*c*rho*nu)
!    rotT(2)=Tnum(i,2)/(pi*c*c*c*rho*nu)
!    rotT(3)=Tnum(i,3)/(pi*c*c*c*rho*nu)
!
!    WRITE (parLogTekst,'(3G20.10)') rotT(1),rotT(2),rotT(3)
!    CALL WriteToLog(parLogTekst)
!
!
!  end do







  DEALLOCATE(stk%sysMcrs%v)
  DEALLOCATE(rhsM)

  DEALLOCATE(Fnum,Fdrag,Flift)


end subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesInlet3sidesTorque()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: sysM(:,:),rhsM(:,:),b(:)

  integer i,j,k,nid,outside,particle,lun
  real(rk) flowVelocity(3)
  real(rk) c,rho,nu,resT(3),cdre(3),rotT(3),cmre(3)

  REAL(rk), ALLOCATABLE :: Fnum(:,:),Tnum(:,:),Fdrag(:,:),Flift(:,:)

  rho=1.0
  nu=1.0
  call getSupElc(parSuelLam1,parSuelLam2,parSuelE1,parSuelE2,c)
  

  ALLOCATE(sysM(stk%neq,stk%neq))
  ALLOCATE(rhsM(stk%neq,stk%nb))

  ALLOCATE (Fnum(3,3),Tnum(3,3))
  ALLOCATE (Fdrag(3,3),Flift(3,3))
  Fnum=0.0_rk
  Tnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk

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


  CALL WriteToLog("Stokes: solving.")
  !
  ! Calculate preconditioner
  !          
  ierr=0
  !
  ! solve full system of linear equations
  !
  stk%slv%type  = eqn(1)%slv%type ! 0 = direct solver, 1 = lsqr
  stk%slv%maxit = eqn(1)%slv%maxit
  stk%slv%stopt = eqn(1)%slv%stopt
  stk%slv%eps   = eqn(1)%slv%eps 

  IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver

    stk%slv%pret = 2 ! 2 = LU faktorizacija
    stk%slv%prep = 2 ! 2 = preconditioner calculated outside 
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL FormPREfm(stk%slv%pret,stk%neq,stk%Pivot,sysM,cput,ierr)   
    DEALLOCATE(sysM)

  ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver

    CALL  TransformToCRS(sysM,stk%neq,stk%neq,stk%sysMcrs)
    DEALLOCATE(sysM)
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)

  END IF

  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1


  do i=1,3
    !
    ! Get flow velocity
    !
    flowVelocity(i)=0.0_rk

    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   

            if (i.eq.1) then ! rotation around x axis
              flowVelocity(2) =   node(nid)%x(3)
              flowVelocity(3) = - node(nid)%x(2)
            end if

            if (i.eq.2) then ! rotation around y axis
              flowVelocity(1) = - node(nid)%x(3)
              flowVelocity(3) =   node(nid)%x(1)
            end if

            if (i.eq.3) then ! rotation around z axis
              flowVelocity(2) = - node(nid)%x(1)
              flowVelocity(1) =   node(nid)%x(2)
            end if


            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)                
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

    ALLOCATE (b(stk%neq))              
    b = MATMUL(rhsM,stk%b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0_rk
    nits=0
              
    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
      !
      ! Solve with direct solver
      !
      CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                     stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
            
    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 

        
    END IF
    !
    !  Report to log file
    !
!            WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
!            CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    ! force coeficient
    cdre(i)=4.0_rk*resT(i)
    !
    !  Calculate torque numerically
    !
    call getTorqueIntegrateFluxes(particle,Tnum(i,:))
    ! Calculate rotation tensor diagonal components
    rotT(i)=Tnum(i,i)/(pi*c*c*c*rho*nu)
    ! moment coeficient
    cmre(i)=4.0_rk*rotT(i)

    WRITE (parLogTekst,'(I0,1X,6G20.10)') i,flowVelocity,Fnum(i,i),resT(i),cdre(i)
    CALL WriteToLog(parLogTekst)

  end do


  DEALLOCATE(stk%sysMcrs%v)
  DEALLOCATE(rhsM)

  lun=58
  OPEN (unit=lun,file=TRIM("and.rotT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 Omegaxx Omegayy Omegazz cmreX cmreY cmreZ Tx Ty Tz"
  WRITE (lun,'(16G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                       rotT(1),rotT(2),rotT(3),cmre(1),cmre(2),cmre(3),Tnum(1,1),Tnum(2,2),Tnum(3,3)
  CLOSE (lun)


  DEALLOCATE(Fnum,Tnum,Fdrag,Flift)


end subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesInlet3sides()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr,isd
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: sysM(:,:),rhsM(:,:),b(:)

  integer i,j,k,nid,outside,particle,lun
  real(rk) flowVelocity(3)
  real(rk) c,rho,nu,resT(3),cdre(3)

  REAL(rk), ALLOCATABLE :: Fnum(:,:),Fdrag(:,:),Flift(:,:)

  rho=1.0
  nu=1.0
  call getSupElc(parSuelLam1,parSuelLam2,parSuelE1,parSuelE2,c)
  

  ALLOCATE(sysM(stk%neq,stk%neq))
  ALLOCATE(rhsM(stk%neq,stk%nb))

  ALLOCATE (Fnum(3,3))
  ALLOCATE (Fdrag(3,3),Flift(3,3))
  Fnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk

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


  CALL WriteToLog("Stokes: solving.")
  !
  ! Calculate preconditioner
  !          
  ierr=0
  !
  ! solve full system of linear equations
  !
  stk%slv%type  = eqn(1)%slv%type ! 0 = direct solver, 1 = lsqr
  stk%slv%maxit = eqn(1)%slv%maxit
  stk%slv%stopt = eqn(1)%slv%stopt
  stk%slv%eps   = eqn(1)%slv%eps 

  IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver

    stk%slv%pret = 2 ! 2 = LU faktorizacija
    stk%slv%prep = 2 ! 2 = preconditioner calculated outside 
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL FormPREfm(stk%slv%pret,stk%neq,stk%Pivot,sysM,cput,ierr)   
    DEALLOCATE(sysM)

  ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver

    CALL  TransformToCRS(sysM,stk%neq,stk%neq,stk%sysMcrs)
    DEALLOCATE(sysM)
    !
    ! preconditioner (calculate : stk%Pivot )
    !
    CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)

  END IF

  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1


  do i=1,3
    !
    ! Get flow velocity
    !
    flowVelocity = 0.0_rk
    flowVelocity(i)=1.0_rk
    call normalize(flowVelocity)
    !
    ! Apply flow velocity as boundary condition
    !
    do k=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(k).eq.outside) then
            nid = subdomain(isd)%nodeList(k)   
            DO j=1,3
                eqn(j)%u(nid) = flowVelocity(j)
            END DO
        end if
    end do
    !
    ! Set up left and right hand side vectors
    !
    call formStokesLeftRightHandSideVectors()

    ALLOCATE (b(stk%neq))              
    b = MATMUL(rhsM,stk%b)
 
    !
    ! solve
    !
  
    ierr=0
    cput=0.0
    nits=0
              
    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
      !
      ! Solve with direct solver
      !
      CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                     stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
            
    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
      !
      ! Solve with LSQR solver
      !
      CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
      stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 

        
    END IF
    !
    !  Report to log file
    !
!            WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
!            CALL WriteToLog(parLogTekst)
    IF (stk%slv%maxit.LE.nits) THEN
      WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
      CALL WriteToLog(parLogTekst)
    END IF
    !
    !  Distribute SLE results to u and q vectors
    !
    call StokesDistributeUnknowns()
    !
    !  Free memory
    !      
    DEALLOCATE(b)
    !
    !  Calculate force numerically
    !
    call getForceIntegrateFluxes(particle,Fnum(i,:))
    !
    ! Calculate drag and lift from total force
    !
    call getDragLift(flowVelocity,Fnum(i,:),Fdrag(i,:),Flift(i,:))
    ! 
    ! Calculate resistance tensor diagonal components
    resT(i)=Fnum(i,i)/(pi*c*rho*nu)
    cdre(i)=8.0_rk*c*resT(i)

    WRITE (parLogTekst,'(I0,1X,6G20.10)') i,flowVelocity,Fnum(i,i),resT(i),cdre(i)
    CALL WriteToLog(parLogTekst)

  end do


  DEALLOCATE(stk%sysMcrs%v)
  DEALLOCATE(rhsM)

  lun=58
  OPEN (unit=lun,file=TRIM("and.resT.dat"),status="UNKNOWN")
  WRITE (lun,*) "# a b c l1 l2 e1 e2 Kxx Kyy Kzz cdreX cdreY cdreZ Fx Fy Fz"
  WRITE (lun,'(19G20.10)') parSuelLam1*c,parSuelLam2*c,c,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2, &
                           resT(1),resT(2),resT(3),cdre(1),cdre(2),cdre(3),Fnum(1,1),Fnum(2,2),Fnum(3,3)
  CLOSE (lun)


  DEALLOCATE(Fnum,Fdrag,Flift)


end subroutine


!
!
! -----------------------------------------------------------------------------------------
subroutine getSupElc(l1,l2,e1,e2,c)
  USE mPar
  implicit NONE
  real(rk) l1,l2,e1,e2,c,d,V,b1,b2

  ! assuming sphere with d=1
  d=1
  V=4_rk*pi/3.0_rk*(d/2.0_rk)**3

  call beta(e1/2.0_rk+1.0_rk,e1,b1)
  call beta(e2/2.0_rk,e2/2.0_rk,b2)

  c=(V/(2.0_rk*l1*l2*e1*e2*b1*b2))**(1.0_rk/3.0_rk)

end subroutine

! -----------------------------------------------------------------------------------------
subroutine beta(x,y,b)
  USE mPar
  implicit none
  real(rk) x,y,b
  b = Gamma(x)*Gamma(y)/Gamma(x+y)
end subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesEllipsoidValidationSOLVE()

    USE mMesh
    USE mEqns
    USE mPar
    USE mCommon
    IMPLICIT NONE

    INTEGER nits,ierr,isd
    REAL(rk)    cput
    REAL(rk), ALLOCATABLE :: sysM(:,:),rhsM(:,:),b(:)

    integer i,j,k,nid,outside,particle,noid
    real(rk) flowVelocity(3)
    real(rk) st,im

    REAL(rk), ALLOCATABLE :: Fnum(:,:),Fana(:,:),err(:),Fdrag(:,:),Flift(:,:)

    ALLOCATE(sysM(stk%neq,stk%neq))
    ALLOCATE(rhsM(stk%neq,stk%nb))

    ALLOCATE (Fnum(nnodes,3),Fana(nnodes,3),err(nnodes))
    ALLOCATE (Fdrag(nnodes,3),Flift(nnodes,3))
    Fnum=0.0_rk
    Fdrag=0.0_rk
    Flift=0.0_rk
    Fana=0.0_rk
    err=0.0_rk

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


    CALL WriteToLog("Stokes: solving.")
    !
    ! Calculate preconditioner
    !          
    ierr=0
    !
    ! solve full system of linear equations
    !
    stk%slv%type  = eqn(1)%slv%type ! 0 = direct solver, 1 = lsqr
    stk%slv%maxit = eqn(1)%slv%maxit
    stk%slv%stopt = eqn(1)%slv%stopt
    stk%slv%eps   = eqn(1)%slv%eps 

    IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver
  
      stk%slv%pret = 2 ! 2 = LU faktorizacija
      stk%slv%prep = 2 ! 2 = preconditioner calculated outside 
      !
      ! preconditioner (calculate : stk%Pivot )
      !
      CALL FormPREfm(stk%slv%pret,stk%neq,stk%Pivot,sysM,cput,ierr)   
      DEALLOCATE(sysM)
  
    ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver
  
      CALL  TransformToCRS(sysM,stk%neq,stk%neq,stk%sysMcrs)
      DEALLOCATE(sysM)
      !
      ! preconditioner (calculate : stk%Pivot )
      !
      CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)
  
    END IF

    !
    !    Loop over outside nodes
    !
    outside = 2
    particle = 1
    isd=1

    CALL WriteToLog("noid nits flowVelocity Ftot Fdrag Flift")


    do i=1,subdomain(isd)%nnodes
        if (subdomain(isd)%BCidList(i).eq.outside) then
            !
            ! Get flow velocity
            !
            noid = subdomain(isd)%nodeList(i)            
            flowVelocity = node(noid)%x
            call normalize(flowVelocity)
!            d=SQRT(flowVelocity(1)**2+flowVelocity(2)**2+flowVelocity(3)**2)
!            DO j=1,3
!                flowVelocity(j)=flowVelocity(j)/d
!            END DO

!            print *,"outside node =",noid
!            print *,flowVelocity

            !
            ! Apply flow velocity as boundary condition
            !
            do k=1,subdomain(isd)%nnodes
                if (subdomain(isd)%BCidList(k).eq.outside) then
                    nid = subdomain(isd)%nodeList(k)   
                    DO j=1,3
                        eqn(j)%u(nid) = flowVelocity(j)
                    END DO
                end if
            end do
            !
            ! Set up left and right hand side vectors
            !
            call formStokesLeftRightHandSideVectors()

            ALLOCATE (b(stk%neq))              
            b = MATMUL(rhsM,stk%b)
            !DEALLOCATE(rhsM)
          
            !
            ! solve
            !
          
            ierr=0
            cput=0.0_rk
            nits=0
                      
            IF (stk%slv%type.EQ.0) THEN ! 0 = direct solver         
              !
              ! Solve with direct solver
              !
              CALL SolvEQNfm(stk%slv%type,stk%slv%pret,stk%slv%prep,stk%slv%maxit,stk%slv%stopt,stk%slv%eps, &
                             stk%neq,stk%Pivot,sysM,b,stk%x,nits,cput,ierr)
                    
            ELSE IF (stk%slv%type.EQ.1) THEN ! 1 = lsqr solver          
              !
              ! Solve with LSQR solver
              !
              CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
              stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 
!              DEALLOCATE(stk%sysMcrs%v)
          
            END IF
            !
            !  Report to log file
            !
!            WRITE (parLogTekst,'(A,A,I0,1X,I0,1X,I0)') TRIM(stk%name)," solver (type,ierr,nits) : ",stk%slv%type,ierr,nits
!            CALL WriteToLog(parLogTekst)
            IF (stk%slv%maxit.LE.nits) THEN
              WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
              CALL WriteToLog(parLogTekst)
            END IF
            !
            !  Distribute SLE results to u and q vectors
            !
            call StokesDistributeUnknowns()
            !
            !  Free memory
            !      
            DEALLOCATE(b)
            !
            !  Calculate force numerically
            !
            call getForceIntegrateFluxes(particle,Fnum(noid,:))
            !
            !  Calculate force analytically
            !
            call calForceOnParticle(flowVelocity,Fana(noid,:))
            !
            ! Calculate drag and lift from total force
            !
            call getDragLift(flowVelocity,Fnum(noid,:),Fdrag(noid,:),Flift(noid,:))
            !call getDragLift(flowVelocity,Fana(noid,:),Fdrag(noid,:),Flift(noid,:))

            !
            !  Estimate error
            !
            st = 0.0_rk
            im = 0.0_rk
            do j=1,3
                st = st + (Fana(noid,j)-Fnum(noid,j))**2
                im = im + (Fana(noid,j))**2
            end do
            err(noid) = sqrt(st/im)

            !WRITE (parLogTekst,'(I0,1X,4G20.10)') noid,flowVelocity,err(noid)
            WRITE (parLogTekst,'(I0,1X,I0,1X,12G20.10)') noid,nits,flowVelocity,Fnum(noid,:),Fdrag(noid,:),Flift(noid,:)
            CALL WriteToLog(parLogTekst)
            !print *, err(noid)
        end if
    end do


    call OutputForceErrorParaview(outside,Fnum,Fana,err,Fdrag,Flift)


    DEALLOCATE(stk%sysMcrs%v)
    DEALLOCATE(rhsM)
    DEALLOCATE(Fnum,Fana,err,Fdrag,Flift)


end subroutine


subroutine getDragLift(v,F,drag,lift)
  use mPar
  implicit none
  real(rk) v(3),F(3),drag(3),lift(3),dir(3),dp

  ! velocity direction
  dir = v
  call normalize(dir)

  call DotProduct(F,dir,dp)
  drag = dir * dp
  lift = F - drag

end subroutine


subroutine getLatLon(pointOnSphere,lat,lon)
  USE mPar
  implicit NONE
  real(rk) r(3),pointOnSphere(3),lat,lon

  r = pointOnSphere
  call normalize(r)

  lat = acos(r(3)); !theta (sever - jug), 0 -- pi
  lon = atan2(r(2),r(1));   !phi (vzhod zahod), 0 -- 2pi
  if (lon.lt.0.0_rk) lon = 2.0_rk*pi + lon

  ! (1,0,0) -> lat = pi/2, lon = 0
  ! (0,1,0) -> lat = pi/2, lon = pi/2

end subroutine


subroutine normalize(v)
  USE mPar
  implicit none 
  real(rk) v(3),d
  integer j

  d=SQRT(v(1)**2+v(2)**2+v(3)**2)

  if (d.ne.0.0_rk) then
    DO j=1,3
        v(j)=v(j)/d
    END DO
  end if

end subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE getForceIntegrateFluxes(wallID,F)

    USE mEqns
    USE mPar
    IMPLICIT NONE
    integer en,wallID
    real(rk) F(3)
  
    DO en=1,neq
      call CalFluxIntegralSingleWall(eqn(en)%q,wallID,F(en))
    END DO ! en
 
  
  end subroutine

!
!     ------------------------------------------------------------------
!
SUBROUTINE getElementCenter(ie,elementCenter)
  
  USE mCommon
  use mMesh
  IMPLICIT NONE
  
  integer ie,in,j
  real(rk) elementCenter(3)
   
  elementCenter = 0.0_rk
  DO in=1,element(ie)%nno
    DO j=1,3
      elementCenter(j)=elementCenter(j)+node(element(ie)%con(in))%x(j)
    END DO
  END DO
  DO j=1,3
    elementCenter(j)=elementCenter(j)/DBLE(element(ie)%nno)
  END DO

end subroutine

!
!     ------------------------------------------------------------------
!
  SUBROUTINE getTorqueIntegrateFluxes(wallID,T)

    USE mEqns
    USE mPar
    use mMesh
    IMPLICIT NONE
    integer wallID,i,j
    real(rk) T(3),cp(3),r(3),q(3)
  
    T = 0.0_rk

    DO i=1,nelem
      IF (element(i)%bcid.EQ.wallID) THEN
        !
        !  Set up traction vector
        !
        DO j=1,3
          q(j) = eqn(j)%q(i)
        END DO
        !
        ! Radij vector to mesh element
        !
        CALL getElementCenter(i,r)
        ! 
        ! Cross product
        !
        CALL CrossProduct(cp,r,q)
        !
        ! Add to integral (constant interpolation)
        !
        DO j=1,3
          T(j) = T(j) + element(i)%area * cp(j)
        END DO

      END IF
    END DO
     
  
  end subroutine



!
!
! ----------------------------------------------------------------------
!
  SUBROUTINE CalFluxIntegralSingleWall(q,wallID,integral)
!
! ----------------------------------------------------------------------
    USE mMesh
    USE mPar
  
  
    implicit none
  
    real(rk) integral
    real(rk) q(nqnodes)
    integer i,wallID
            
    !
    !     Calculate integral of flux
    !
    integral = 0.0_rk
    DO i=1,nelem
      IF (element(i)%bcid.EQ.wallID) THEN
        integral = integral + element(i)%area * q(i)
      END IF
    END DO
    
    
    END SUBROUTINE
            
    


! -----------------------------------------------------------------------------------------
SUBROUTINE CalResTensorPrime(lambda,ResTprime)
!
!     Calculate resistance tensor in particle frame of reference
!
! -----------------------------------------------------------------------------------------
USE mPar
IMPLICIT NONE
REAL(rk) lambda,ResTprime(3),l2

ResTprime=0.0_rk
IF (lambda.LT.1.0_rk) THEN
  CALL WriteToLog("ERROR :: CalResTensorPrime :: Particle aspect ratio < 1 !")   
ELSE IF (lambda.EQ.1.0_rk) THEN
  ResTprime(1)=6.0_rk
  ResTprime(2)=6.0_rk
  ResTprime(3)=6.0_rk
ELSE
  l2=Lambda**2 - 1.0_rk
  ResTprime(2)=(16.0_rk*l2**(1.5))/( (2*l2 - 1.0_rk)*Log(Lambda + Sqrt(l2)) + Lambda*Sqrt(l2) )
  ResTprime(3)=ResTprime(2)
  ResTprime(1)=( 8.0_rk*l2**(1.5))/( (2*l2 + 1.0_rk)*Log(Lambda + Sqrt(l2)) - Lambda*Sqrt(l2) )
END IF

END

subroutine calForceOnParticle(u,F)

    use mPar
    implicit none
    real(rk) b,ResTprime(3),u(3),F(3),lambda,nu,rho

    lambda = 2.0_rk
    b = 0.5_rk
    nu = 1.0_rk
    rho = 1.0_rk
!    u(1) = 0.3_rk
!    u(2) = -0.5_rk
!    u(3) = 0.812404_rk

    call CalResTensorPrime(lambda,ResTprime)

    F(1)=resTprime(1)*u(1)*b*pi*nu*rho
    F(2)=resTprime(2)*u(2)*b*pi*nu*rho
    F(3)=resTprime(3)*u(3)*b*pi*nu*rho

!    print *,resTprime
!    print *,F

end subroutine


subroutine calForceOnSphere(u,F,V)
  use mPar
  implicit none

  real(rk) R,u(3),F(3),nu,rho,V

  !V = pi/6.0_rk
  R = (3.0_rk*V/4.0_rk/pi)**(1.0_rk/3.0_rk)
  nu = 1.0_rk
  rho = 1.0_rk

  F(1) = 6.0_rk *u(1)*R*pi*nu*rho
  F(2) = 6.0_rk *u(2)*R*pi*nu*rho
  F(3) = 6.0_rk *u(3)*R*pi*nu*rho

end subroutine


subroutine calForceOnEllipsoid(u,F,lambda,V)

  use mPar
  implicit none
  real(rk) R,ResTprime(3),u(3),F(3),lambda,nu,rho,V

  ! semi-minor axis
  R = (3.0_rk*V/4.0_rk/pi/lambda)**(1.0_rk/3.0_rk)
  nu = 1.0_rk
  rho = 1.0_rk

  call CalResTensorPrime(lambda,ResTprime)

  F(1)=resTprime(1)*u(1)*R*pi*nu*rho
  F(2)=resTprime(2)*u(2)*R*pi*nu*rho
  F(3)=resTprime(3)*u(3)*R*pi*nu*rho

  !print *,R
  !print *,resTprime
  !print *,F

end subroutine


! -----------------------------------------------------------------------------

SUBROUTINE OutputForceErrorParaview(iWall,Fnum,Fana,err,Fdrag,Flift)

    !
    !     $: Each wall in separate file
    !
    ! -----------------------------------------------------------------------------
          USE mMesh
          USE mEqns
          USE mPar
          USE mCommon
          IMPLICIT NONE
    
          INTEGER itp,i,of,k,j
          CHARACTER*255 filename,tmp,vrstica
    
          INTEGER iWall,noe,non
          INTEGER, ALLOCATABLE :: nlist(:)

          real(rk) Fnum(nnodes,3),Fana(nnodes,3),err(nnodes),lat,lon
          real(rk) Fdrag(nnodes,3),Flift(nnodes,3)
    
    
          itp=96
        
          WRITE (filename,'(3A)') TRIM(parResultsForceErr),TRIM(wall(iWall)%name),".vtu"
          OPEN (itp,FILE=TRIM(filename),STATUS='UNKNOWN')
    
    !     count number of nodes in wall
    
          noe = 0
          ALLOCATE (nlist(nnodes))
          nlist=0
          DO  i=1,nelem
            IF (element(i)%bcid.EQ.iWall) THEN
              noe = noe + 1
              DO j=1,element(i)%nno
                nlist(element(i)%con(j))=1
              END DO
            END IF
          END DO
          non=SUM(nlist)
    
          j=0
          DO  i=1,nnodes
            IF (nlist(i).EQ.1) THEN
              j=j+1
              nlist(i)=j
            END IF
          END DO
    
    
    
    
    
          WRITE (itp,'(A)') '<?xml version="1.0"?>'
          WRITE (itp,'(A)')&
         '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
          WRITE (itp,'(A)') '<UnstructuredGrid>'
          WRITE (itp,'(A,I10,A,I10,A)') '<Piece NumberOfPoints="',non,'" NumberOfCells="',noe,'">'
    
          WRITE (tmp,'("(",I1,"A)")') neq+2
          WRITE (itp,tmp) '<PointData Scalars="node,err,lat,lon" Vectors="Fana,Fnum,Fdrag,Flift">'
    !
    !     PODATKI V VOZLISCIH
    !
          WRITE (itp,'(A)') '<DataArray type="Float32" Name="node" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) WRITE (itp,*) nlist(i)
          END DO
          WRITE (itp,'(A)') '</DataArray>'
    
          WRITE (itp,'(A)') '<DataArray type="Float32" Name="err" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) WRITE (itp,*) err(i)
          END DO
          WRITE (itp,'(A)') '</DataArray>'          

          WRITE (itp,'(A)') '<DataArray type="Float32" Name="lat" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) THEN
              call getLatLon(node(i)%x,lat,lon)
              WRITE (itp,*) lat
            END IF
          END DO
          WRITE (itp,'(A)') '</DataArray>'    
          
          WRITE (itp,'(A)') '<DataArray type="Float32" Name="lon" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) THEN
              call getLatLon(node(i)%x,lat,lon)
              WRITE (itp,*) lon
            END IF
          END DO
          WRITE (itp,'(A)') '</DataArray>'              


          WRITE (itp,'(A)') '<DataArray type="Float32" Name="Fana" NumberOfComponents="3" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) THEN 
              WRITE (vrstica,*) Fana(i,1),Fana(i,2),Fana(i,3)
              CALL sqblnk(itp,vrstica)
            END IF
          END DO
          WRITE (itp,'(A)') '</DataArray>'

          WRITE (itp,'(A)') '<DataArray type="Float32" Name="Fnum" NumberOfComponents="3" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) THEN 
              WRITE (vrstica,*) Fnum(i,1),Fnum(i,2),Fnum(i,3)
              CALL sqblnk(itp,vrstica)
            END IF
          END DO
          WRITE (itp,'(A)') '</DataArray>'

          WRITE (itp,'(A)') '<DataArray type="Float32" Name="Fdrag" NumberOfComponents="3" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) THEN 
              WRITE (vrstica,*) Fdrag(i,1),Fdrag(i,2),Fdrag(i,3)
              CALL sqblnk(itp,vrstica)
            END IF
          END DO
          WRITE (itp,'(A)') '</DataArray>'

          WRITE (itp,'(A)') '<DataArray type="Float32" Name="Flift" NumberOfComponents="3" format="ascii">'
          DO  i=1,nnodes
            IF (nlist(i).NE.0) THEN 
              WRITE (vrstica,*) Flift(i,1),1E-15_rk+Flift(i,2),Flift(i,3)                  
              CALL sqblnk(itp,vrstica)
            END IF
          END DO
          WRITE (itp,'(A)') '</DataArray>'



          WRITE (itp,'(A)') '</PointData>'
    
    
    !      KOORDINATE VOZLISC
          WRITE (itp,'(A)') '<Points>'
          WRITE (itp,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
    
          DO  i=1,nnodes
             IF (nlist(i).NE.0) THEN
               WRITE (vrstica,'(3G18.9)') node(i)%x(1),node(i)%x(2),node(i)%x(3)
               CALL sqblnk(itp,vrstica)
             END IF
          END DO
    
          WRITE (itp,'(A)') '</DataArray>'
          WRITE (itp,'(A)') '</Points>'
    
          WRITE (itp,'(A)') '<Cells>'
            WRITE (itp,'(A)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
    
            DO i=1,nelem
              IF (element(i)%bcid.EQ.iWall) THEN
                WRITE(vrstica,*) (nlist(element(i)%con(k))-1,k=1,element(i)%nno)
                CALL sqblnk(itp,vrstica)
              END IF
            END DO
    
    
            WRITE (itp,'(A)') '</DataArray>'
            WRITE (itp,'(A)') '<DataArray type="Int32" Name="offsets" format="ascii">'
              of = 0
              DO i=1,nelem
                IF (element(i)%bcid.EQ.iWall) THEN
                  of = of + element(i)%nno
                  WRITE (itp,*) of
                END IF
              END DO
            WRITE (itp,'(A)') '</DataArray>'
            WRITE (itp,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
              DO i=1,nelem
              IF (element(i)%bcid.EQ.iWall) THEN
                IF (element(i)%type.EQ.2) THEN ! three node triangle
                  WRITE (itp,*) "5" ! 5 tikotnik, 9 quuad, 12 heksaeder
                ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle
                 WRITE (itp,*) "9" ! 5 tikotnik, 9 quuad, 12 heksaeder
                ELSE
                  CALL WriteToLog("Error :: Element type not supported!")
                END IF
              END IF
              END DO
            WRITE (itp,'(A)') '</DataArray>'
          WRITE (itp,'(A)') '</Cells>'
    
    
          WRITE (itp,'(A)') '</Piece>'
          WRITE (itp,'(A)') '</UnstructuredGrid>'
          WRITE (itp,'(A)') '</VTKFile>'
    
          CLOSE (itp)
    
          DEALLOCATE (nlist)
    
    
          END
    

subroutine sphere2superE()
    use mPar
    use mMesh
    implicit none
    real(rk) selX(3),r
    integer i

    do i=1,nnodes
        r = node(i)%x(1)**2 + node(i)%x(2)**2 + node(i)%x(3)**2

        if (r.LT.100.0_rk) then
            call mapSphereToSuperE(parMTsel(1),parMTsel(2),parMTsel(3),parMTsel(4),parMTsel(5),selX,node(i)%x)
            node(i)%x = selX
            !print *,i,parMTsel(1),parMTsel(2),parMTsel(3),parMTsel(4),parMTsel(5)
        end if
    end do
        

end subroutine


! -----------------------------------------------------------------------------------------
SUBROUTINE mapSphereToSuperE(a,b,c,e1,e2,t,r)
!
!     t = point on super ellipsoid (local ks)
!     r = point on sphere (|r|=1,c=(0,0,0)) (local ks)
!     (eq 2.10 in Jaklic : Segmentation and Recovery of Superquadrics)
!
! -----------------------------------------------------------------------------------------
    USE mPar
    IMPLICIT NONE

    REAL(rk) a,b,c,e1,e2
    REAL(rk) t(3),r(3)
    REAL(rk) eta,om,xy

    om=ATAN2(r(2),r(1))
    xy=SQRT(r(1)**2+r(2)**2)
    eta=ATAN2(r(3),xy)

    ! SIGN(A,B) returns the value of A with the sign of B
    t(1)=a * SIGN( (ABS(COS(eta)))**e1 , COS(eta) ) * SIGN( (ABS(COS(om)))**e2 , COS(om) )
    t(2)=b * SIGN( (ABS(COS(eta)))**e1 , COS(eta) ) * SIGN( (ABS(SIN(om)))**e2 , SIN(om) )
    t(3)=c * SIGN( (ABS(SIN(eta)))**e1 , SIN(eta) )

END SUBROUTINE
    

! -----------------------------------------------------------------------------------------
subroutine projectPoints2superE()
!
!     Project mesh points to supere surface
!
! -----------------------------------------------------------------------------------------  
  use mPar
  use mMesh
  implicit none
  real(rk) r
  integer i,k,maxk
  
  real(rk) u0,u,res,np(3)
  real(rk) a,b,c,e1,e2,f,df
  
  !PRINT *, "Entering Subroutine :: projectPoints2superE"
  
  DO i=1,nnodes
    r = node(i)%x(1)**2 + node(i)%x(2)**2 + node(i)%x(3)**2
    IF (r.LT.10.0_rk) THEN
  
        np(1) = node(i)%x(1)
        np(2) = node(i)%x(2)
        np(3) = node(i)%x(3)
  
        CALL NormVector(np)
  
        a = parMTsel(1)
        b = parMTsel(2)
        c = parMTsel(3)
        e1 = parMTsel(4)
        e2 = parMTsel(5)
        
        u0 = 1.0D0
        res = 1.0D0
      
        k = 0
        maxk = 1000
        DO WHILE ( res > 1.0e-3 )
      
          call calF(a,b,c,e1,e2,u0,np,f)
          call calDF(a,b,c,e1,e2,u0,np,df)
          u = u0 - f / df
          !PRINT *, "u=",u,f(u0,np),df(u0,np)
      
          call calF(a,b,c,e1,e2,u,np,f)
          res = abs( f )
          u0 = u
      
          IF (k > maxk) EXIT
          k = k+1
      
        END DO
      
        node(i)%x(1) = u * np(1)
        node(i)%x(2) = u * np(2)
        node(i)%x(3) = u * np(3)
      
        !PRINT *, "u=",u
        !PRINT *, "res=",res
        !PRINT *, "k=",k
        !PRINT *, "p=",node(i)%x 
  
    END IF
  END DO
  
end subroutine
        
subroutine calF(a,b,c,e1,e2,t,p,f)
  USE mPar
  IMPLICIT NONE
   
   real(rk) t,p(3),f
   real(rk) x,y,z,a,b,c,e1,e2
  !real(rk) sign_x, sign_y, sign_z
  !sign_x = SIGN( real(1,rk), x )
  !sign_y = SIGN( real(1,rk), y )
  !sign_z = SIGN( real(1,rk), z )
   x = abs( p(1) )
   y = abs( p(2) )
   z = abs( p(3) )
   f = ( (t*x/a)**(2/e2) + (t*y/b)**(2/e2) )**(e2/e1) + (t*z/c)**(2/e1) - 1
   !print *, "f+1=",f+1
   
   RETURN
 end subroutine
 
subroutine  calDF(a,b,c,e1,e2,t,p,df)
 use mPar
 implicit none 
 real(rk) df,t,p(3)
 real(rk) x,y,z,a,b,c,e1,e2

  x = abs( p(1) )
  y = abs( p(2) )
  z = abs( p(3) )
  
  df = 2 * ( ((t*x/a)**(2/e2) + (t*y/b)**(2/e2))**(e2/e1) + (t*z/c)**(2/e1) )
  df = df / (e1*t)

  RETURN
end subroutine   