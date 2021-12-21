!
!     ------------------------------------------------------------------
!
SUBROUTINE pressureStokes()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER row,isd,icol,j,i
  REAL(rk) osemPi
  REAL(rk), ALLOCATABLE :: vx(:),vy(:),vz(:)


  osemPi = ATAN(1.0)*4.0_rk*8.0_rk


  ALLOCATE (vx(nqnodes),vy(nqnodes),vz(nqnodes))

  CALL InterpolateToQMesh(eqn(1)%u,vx)
  CALL InterpolateToQMesh(eqn(2)%u,vy)
  CALL InterpolateToQMesh(eqn(3)%u,vz)
  !
  ! Loop over subdomains
  !
  isd=1
  row = 0

  DO i = nnodes+1,nnodes+nqnodes
    row = row + 1

    eqn(1)%p(row) = 0.0_rk
    DO j = 1,subdomain(isd)%nnodes
      icol = subdomain(isd)%nodeList(j)  ! where to take column value

      ! Ingber 1991, Surface pressure solution ...
      eqn(1)%p(row) = eqn(1)%p(row) + &
                        subdomain(isd)%diff * ( subdomain(isd)%prTxMat(i,j) * ( eqn(1)%u(icol)-vx(row) ) &
                                              + subdomain(isd)%prTyMat(i,j) * ( eqn(2)%u(icol)-vy(row) ) &
                                              + subdomain(isd)%prTzMat(i,j) * ( eqn(3)%u(icol)-vz(row) ) )  


    END DO 
    DO j = 1,subdomain(isd)%nqnodes
      icol = subdomain(isd)%qnodeList(j)  ! where to take column value
      eqn(1)%p(row) = eqn(1)%p(row)         - ( subdomain(isd)%prGxMat(i,j) * eqn(1)%q(icol) &
                                              + subdomain(isd)%prGyMat(i,j) * eqn(2)%q(icol) &
                                              + subdomain(isd)%prGzMat(i,j) * eqn(3)%q(icol) )
    END DO 
    eqn(1)%p(row)  = - eqn(1)%p(row) / osemPi * 2.0_rk ! *2 ker na ravni steni

  END DO

  DEALLOCATE(vx,vy,vz)

!  print *,"u tocke"
!  ALLOCATE (p(nnodes))
!  !
!  ! Loop over subdomains
!  !
!  isd=1
!  row = 0
!
!  DO i = 1,nnodes
!    row = row + 1
!    p(row) = 0.0_rk
!    DO j = 1,subdomain(isd)%nnodes
!      icol = subdomain(isd)%nodeList(j)  ! where to take column value
!!      p(row) = p(row) + subdomain(isd)%diff * ( subdomain(isd)%prTxMat(i,j) * eqn(1)%u(icol) &
!!                                              + subdomain(isd)%prTyMat(i,j) * eqn(2)%u(icol) &
!!                                              + subdomain(isd)%prTzMat(i,j) * eqn(3)%u(icol) ) 
!      ! Ingber 1991, Surface pressure solution ...
!      p(row) = p(row) + subdomain(isd)%diff * ( subdomain(isd)%prTxMat(i,j) * (eqn(1)%u(icol)-eqn(1)%u(row)) &
!                                              + subdomain(isd)%prTyMat(i,j) * (eqn(2)%u(icol)-eqn(2)%u(row)) &
!                                              + subdomain(isd)%prTzMat(i,j) * (eqn(3)%u(icol)-eqn(3)%u(row)) )  
!
!    END DO 
!    DO j = 1,subdomain(isd)%nqnodes
!      icol = subdomain(isd)%qnodeList(j)  ! where to take column value
!      p(row) = p(row)                       - ( subdomain(isd)%prGxMat(i,j) * eqn(1)%q(icol) &
!                                              + subdomain(isd)%prGyMat(i,j) * eqn(2)%q(icol) &
!                                              + subdomain(isd)%prGzMat(i,j) * eqn(3)%q(icol) )
!    END DO 
!    p(row)  = - p(row) / osemPi * 2.0_rk
!    print *,row,p(row),node(row)%x(1) !,1.0-8.0*qnode(row)%x(1)
!  END DO
!  eqn(3)%u = p
!  DEALLOCATE(p)



end subroutine


!
!  get force on particle from all sides
!
subroutine StokesDragFromAllSidesCRS()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  implicit none


  integer i,j,k,nid,outside,particle,noid,isd,itp
  real(rk) flowVelocity(3)
  real(rk) lambda,V,Fsphere
  !REAL(rk)    R(3,3),RT(3,3),vr(3),sukanje(3),zasukano(3)

  REAL(rk), ALLOCATABLE :: Fnum(:,:),Fana(:,:),err(:),Fdrag(:,:),Flift(:,:)
  REAL(rk), ALLOCATABLE :: Fs(:), Fbd(:,:),Tnum(:,:)

  ALLOCATE (Fnum(nnodes,3),Fana(nnodes,3),err(nnodes))
  ALLOCATE (Fdrag(nnodes,3),Flift(nnodes,3))
  ALLOCATE (Fs(3),Fbd(nnodes,3),Tnum(nnodes,3))

  Fnum=0.0_rk
  Fdrag=0.0_rk
  Flift=0.0_rk
  Fana=0.0_rk
  err=0.0_rk
  Tnum = 0.0_rk
  !
  !    Loop over outside nodes
  !
  outside = 2
  particle = 1
  isd=1
  lambda = 1.0_rk
  V = pi/6.0_rk
  !
  ! Volume of particle
  !
  call calculateVolumeInsideWall(1,particle,V)
  WRITE (parLogTekst,'(A,G20.10)') "Particle volume = ",V
  CALL WriteToLog(parLogTekst)
  !
  ! Force on a sphere of the same volume
  !
  flowVelocity(1) = 1.0_rk
  flowVelocity(2) = 0.0_rk
  flowVelocity(3) = 0.0_rk
  call calForceOnSphere(flowVelocity,Fs,V)
  Fsphere = Fs(1)
  WRITE (parLogTekst,'(A,G20.10)') "Force on vol. eq. sphere = ",Fsphere
  CALL WriteToLog(parLogTekst)

  itp=96
  OPEN (itp,FILE=TRIM(parFromAllSidesFileName),STATUS='UNKNOWN')

  write(itp,'(A)') "#"
  write(itp,'(A)') "# Andromeda results file, drag force from all sides"
  write(itp,'(A)') "#"
  write(itp,'(A)') "#noid flowVelocity drag drag/dragOnSphere torque"


  !
  ! Loop over outside nodes
  !
  do i=1,subdomain(isd)%nnodes
    if (subdomain(isd)%BCidList(i).eq.outside) then
      !
      ! Get flow velocity
      !
      noid = subdomain(isd)%nodeList(i)            
      flowVelocity = node(noid)%x
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

      ! TEST - vse možne smeri sukanje krogle
      ! v flow velocity je smer okoli katere želim, da se tekočina suče
!      print *,"jure"
!      CALL GetRotationMatrix(R,RT,flowVelocity)
!      do k=1,subdomain(isd)%nnodes
!        if (subdomain(isd)%BCidList(k).eq.outside) then
!            nid = subdomain(isd)%nodeList(k)   
!            sukanje(1) =   0.0_rk ! v x,y,z sukam okoli osi x
!            sukanje(2) = - node(nid)%x(3)
!            sukanje(3) =   node(nid)%x(2)
!            ! zasukam
!            zasukano = MATMUL(RT,sukanje)
!            call normalize(zasukano)
!            print *,zasukano
!            DO j=1,3
!                eqn(j)%u(nid) = zasukano(j)
!            END DO
!        end if
!      end do
!      call OutputInitialParaview()
!      stop


      ! 
      ! Solve
      !
      CALL StokesBigSystemSOLVEcrs()
      !
      ! Postprocess
      !
      !
      !  Calculate force numerically
      !
      call getForceIntegrateFluxes(particle,Fnum(noid,:))
      !
      !  Nondimensional force
      !      
      do j=1,3
        Fbd(noid,j) = Fnum(noid,j) / Fsphere
      end do

      !
      !  Calculate torque numerically
      !
      call getTorqueIntegrateFluxes(particle,Tnum(noid,:))

      !
      ! Calculate drag and lift from total force
      !
      !call getDragLift(flowVelocity,Fnum(noid,:),Fdrag(noid,:),Flift(noid,:))
      !
      !  Calculate force analytically for spheres or ellipsoids
      !
      !call calForceOnEllipsoid(flowVelocity,Fana(noid,:),lambda,V)
      !
      !  Estimate error
      !
      !st = 0.0_rk
      !im = 0.0_rk
      !do j=1,3
      !    st = st + (Fana(noid,j)-Fnum(noid,j))**2
      !    im = im + (Fana(noid,j))**2
      !end do
      !err(noid) = sqrt(st/im)

      write (itp,'(I0,1X,12G18.10)') noid,flowVelocity,Fnum(noid,:),Fbd(noid,:),Tnum(noid,:)
      call flush(itp)
      !WRITE (parLogTekst,'(I0,1X,4G20.10)') noid,flowVelocity,err(noid)
      !WRITE (parLogTekst,'(I0,1X,13G18.10)') noid,err(noid),flowVelocity,Fnum(noid,:),Fdrag(noid,:),Flift(noid,:)
      !CALL WriteToLog(parLogTekst)

    end if
  end do

  close(itp)

  deallocate(Fnum,Fana,err,Fdrag,Flift,Fs,Fbd,Tnum)

end subroutine




!     ------------------------------------------------------------------
!
SUBROUTINE StokesFLOPcrs()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  use linFlowFields
  implicit none


  integer j,k,nid,isd,particle,outside
  real(rk) flowVelocity(3)
  real(rk) V,R,mu,v0,L

  REAL(rk), ALLOCATABLE :: Fnum(:),Tbd(:)
  REAL(rk), ALLOCATABLE :: Fbd(:),Tnum(:)

  ALLOCATE (Fnum(3),Tnum(3))
  ALLOCATE (Fbd(3),Tbd(3))

  Fnum=0.0_rk
  Tnum = 0.0_rk
  Fbd=0.0_rk
  Tbd = 0.0_rk
  mu = 1.0_rk
  v0 = 1.0_rk
  L = 1.0_rk

  !
  ! Volume of particle
  !
  particle = 1
  outside = 2
  isd = 1
  call calculateVolumeInsideWall(1,particle,V)
  R=(3.0_rk*V/(4*pi))**(1.0_rk/3.0_rk)
  WRITE (parLogTekst,'(2(A,G20.10))') "Particle volume = ",V, " Raduis vol. eq. sphere = ",R
  CALL WriteToLog(parLogTekst)
  !
  ! Init particle and its rotation
  !
  call lff_init()
  !
  ! Apply flow velocity as boundary condition
  !
  do k=1,subdomain(isd)%nnodes
      if (subdomain(isd)%BCidList(k).eq.outside) then
          nid = subdomain(isd)%nodeList(k)   
          call lff_getFlowVelocityPFR( node(nid)%x(1),node(nid)%x(2),node(nid)%x(3),flowVelocity)
          DO j=1,3
              eqn(j)%u(nid) = flowVelocity(j)
          END DO
      end if
  end do
  ! 
  ! Solve
  !
  CALL StokesBigSystemSOLVEcrs()
  !
  ! Postprocess
  !
  !
  !  Calculate force numerically
  !
  call getForceIntegrateFluxes(particle,Fnum)
  !
  !  Nondimensional force
  !      
  do j=1,3
    Fbd(j) = Fnum(j) / (pi*R*mu*v0)
  end do
  !
  !  Calculate torque numerically
  !
  call getTorqueIntegrateFluxes(particle,Tnum)
  !
  !  Nondimensional force
  !      
  do j=1,3
    Tbd(j) = Tnum(j) / (pi*R*R*R*mu*v0/L)
  end do  
  !
  !  Write force and torque results to log file
  !
  WRITE (parLogTekst,'(A)') "Force / (pi*R*mu*v0)"
  CALL WriteToLog(parLogTekst)  
  WRITE (parLogTekst,'(3G18.10)') Fbd
  CALL WriteToLog(parLogTekst)
  WRITE (parLogTekst,'(A)') "Torque / (pi*R^3*mu*v0/L)"
  CALL WriteToLog(parLogTekst)    
  WRITE (parLogTekst,'(3G18.10)') Tbd
  CALL WriteToLog(parLogTekst)
  WRITE (parLogTekst,'(A)') "Force"
  CALL WriteToLog(parLogTekst)  
  WRITE (parLogTekst,'(3G18.10)') Fnum
  CALL WriteToLog(parLogTekst)
  WRITE (parLogTekst,'(A)') "Torque"
  CALL WriteToLog(parLogTekst)    
  WRITE (parLogTekst,'(3G18.10)') Tnum
  CALL WriteToLog(parLogTekst)

  deallocate(Fnum,Tbd,Fbd,Tnum)  

end



!     ------------------------------------------------------------------
!
SUBROUTINE StokesBigSystemSOLVEcrs()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER nits,ierr
  REAL(rk)    cput
  REAL(rk), ALLOCATABLE :: b(:)

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
  ! Solve with LSQR solver
  !

  ierr=0
  cput=0.0
  nits=0
  CALL WriteToLog("Stokes: solving.")

  CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
    stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x) 
  !
  ! Report to log file
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
  !DEALLOCATE(b,stk%sysMcrs%v,stk%rhsMcrs%v)
  DEALLOCATE(b)

end


!
!     ------------------------------------------------------------------
!
SUBROUTINE stokesFormCRSsysMrhsM()

  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  integer sysMnnz,rhsMnnz

  INTEGER row,mywall,isd,i,r,col,bc,een
  REAL(rk), ALLOCATABLE :: sysMrow(:),rhsMrow(:)

  !
  ! Form matrices
  !
  CALL WriteToLog("Stokes: counting CRS system and rhs matrices.")
  call countStokesSysRhsMatrices(sysMnnz,rhsMnnz)
  CALL WriteToLog("Stokes: forming CRS system and rhs matrices.")


  stk%sysMcrs%nnz = sysMnnz
  stk%sysMcrs%neq = stk%neq

  stk%rhsMcrs%nnz = rhsMnnz
  stk%rhsMcrs%neq = stk%neq

  ALLOCATE(stk%sysMcrs%v(stk%sysMcrs%nnz),stk%sysMcrs%i(stk%sysMcrs%neq+1),stk%sysMcrs%j(stk%sysMcrs%nnz))
  ALLOCATE(stk%rhsMcrs%v(stk%rhsMcrs%nnz),stk%rhsMcrs%i(stk%rhsMcrs%neq+1),stk%rhsMcrs%j(stk%rhsMcrs%nnz))

  stk%sysMcrs%nnz = 0
  stk%rhsMcrs%nnz = 0
  
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
          stk%sysMcrs%i(row)=stk%sysMcrs%nnz + 1 ! line start        
          stk%rhsMcrs%i(row)=stk%rhsMcrs%nnz + 1 ! line start        
          !
          !  For a single row in matrices
          !
          CALL formUrow(isd,een,r,sysMrow,rhsMrow)
          !
          ! Copy row to full matrix
          !
          DO col = 1,stk%neq
            if (sysMrow(col).NE.0.0_rk) then
              stk%sysMcrs%nnz = stk%sysMcrs%nnz + 1
              stk%sysMcrs%v(stk%sysMcrs%nnz) = sysMrow(col)  ! values
              stk%sysMcrs%j(stk%sysMcrs%nnz) = col         ! column
            end if
          END DO
          DO col = 1,stk%nb
            if (rhsMrow(col).NE.0.0_rk) then
              stk%rhsMcrs%nnz = stk%rhsMcrs%nnz + 1
              stk%rhsMcrs%v(stk%rhsMcrs%nnz) = rhsMrow(col)  ! values
              stk%rhsMcrs%j(stk%rhsMcrs%nnz) = col         ! column
            end if
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
          stk%sysMcrs%i(row)=stk%sysMcrs%nnz + 1 ! line start        
          stk%rhsMcrs%i(row)=stk%rhsMcrs%nnz + 1 ! line start          
          !
          !  For a single row in matrices
          !
          CALL formQrow(isd,een,r,sysMrow,rhsMrow)
          !
          ! Copy row to full matrix
          !
          DO col = 1,stk%neq
            if (sysMrow(col).NE.0.0_rk) then
              stk%sysMcrs%nnz = stk%sysMcrs%nnz + 1
              stk%sysMcrs%v(stk%sysMcrs%nnz) = sysMrow(col)  ! values
              stk%sysMcrs%j(stk%sysMcrs%nnz) = col         ! column
            end if
          END DO
          DO col = 1,stk%nb
            if (rhsMrow(col).NE.0.0_rk) then
              stk%rhsMcrs%nnz = stk%rhsMcrs%nnz + 1
              stk%rhsMcrs%v(stk%rhsMcrs%nnz) = rhsMrow(col)  ! values
              stk%rhsMcrs%j(stk%rhsMcrs%nnz) = col         ! column
            end if
          END DO                        
        END IF
      END DO ! en        
    END DO ! q source point
  
    DEALLOCATE(sysMrow,rhsMrow)

    ! Zacetek prve neobstojece vrstice
    row=row+1      
    stk%sysMcrs%i(row)=stk%sysMcrs%nnz + 1 ! line start        
    stk%rhsMcrs%i(row)=stk%rhsMcrs%nnz + 1 ! line start 

!
! Free memory
!
  do isd=1,nosd
    call dellocateStokesMatrices(isd)
  end do


  !
  ! Solver settings
  !
  stk%slv%type  = eqn(1)%slv%type ! 0 = direct solver, 1 = lsqr
  stk%slv%maxit = eqn(1)%slv%maxit
  stk%slv%stopt = eqn(1)%slv%stopt
  stk%slv%eps   = eqn(1)%slv%eps
  !
  ! preconditioner (calculate : stk%Pivot )
  !
  CALL WriteToLog("Stokes: calculating preconditioner.")
  CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)

  
end subroutine


!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesSolveCRS()

  USE mMesh
  USE mCommon
  USE mEqns
  USE mPar
  IMPLICIT NONE

  REAL(rk), ALLOCATABLE :: b(:)
  INTEGER en,i,isp,ispL,bc,isd,myWall
  INTEGER nits,ierr
  REAL(rk)    cput

  WRITE (parLogTekst,'(A,A)') "Solving for ",TRIM(stk%name)
  CALL WriteToLog(parLogTekst)
  
  !
  !     Set up x, b and rhs-vector vectors
  !
  ALLOCATE (b(stk%neq))
  !
  ! Loop over subdomains
  !
  DO isd=1,nosd
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
  END DO ! subdomains
  

  !
  !     Calculate rhs vector
  !
        CALL CRSxV(stk%rhsMcrs,-stk%b,stk%nb,b)
  !
  !     Solve
  !
        ierr=0
        cput=0.0
        nits=0
  !
  !     preconditioner (calculate : eqn(en)%Pivot )
  !
        print *,"precon"
        CALL lsqr_precon(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v)
  !
  !     solve overdetermined system of linear equations
  !
        stk%slv%eps = 1.0E-7
        print *,"solve"
        CALL lsqr_solve(stk%sysMcrs%neq,stk%nx,stk%sysMcrs%nnz, &
        stk%slv%maxit,stk%slv%eps,nits,ierr,stk%Pivot,stk%sysMcrs%i,stk%sysMcrs%j,stk%sysMcrs%v,b,stk%x)
  
  !
  !       Report to log file
  !
        print *,ierr,nits
        WRITE (parLogTekst,'(A,A,I0,1X,I0)') TRIM(stk%name)," LSQR Solver : ", ierr,nits
        CALL WriteToLog(parLogTekst)
        IF (stk%slv%maxit.LE.nits) THEN
          WRITE (parLogTekst,'(A,I0,1X,I0)') "WARNING :: MAX solver iterations reached!"
          CALL WriteToLog(parLogTekst)
        END IF
  !
  !     Distribute SLE results to u and q vectors
  !
  !

  !
  ! Loop over subdomains
  !
  DO isd=1,nosd
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
        ispL = isp + (en-1)*nqnodes! source point number in the large system 
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
  END DO ! subdomains

  !
  !     Free memory
  !
  DEALLOCATE (b)
  
END
  


!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesBigSystem_RHS()
  USE mMesh
  USE mEqns
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER en,row,mywall,isd,i,j,r,col,mywallC,k,bc,een,isp,ispL
  REAL(rk) val,multi
  INTEGER(8) nnzSys,nnzRHS
  
  nnzSys = 0
  nnzRhs = 0
  !
  ! Count NNZ in CRS matrices
  !
  CALL StokesBigSystemCount(stk%sysMcrs%nnz,stk%rhsMcrs%nnz)

  stk%sysMcrs%neq = stk%neq
  stk%rhsMcrs%neq = stk%neq
  
  WRITE (parLogTekst,'(A,A,F10.2,A)') TRIM(stk%name)," system matrix size = ", &
  (DBLE(stk%sysMcrs%nnz)+0.5*DBLE(stk%sysMcrs%neq))*8.0/1024.0/1024.0," Mb"
  CALL WriteToLog(parLogTekst)
  
  WRITE (parLogTekst,'(A,A,F10.2,A)') TRIM(stk%name)," RHS    matrix size = ", &
  (DBLE(stk%rhsMcrs%nnz)+0.5*DBLE(stk%rhsMcrs%neq))*8.0/1024.0/1024.0," Mb"
  CALL WriteToLog(parLogTekst)

  
  !
  ! Allocate memory for CRS matrix
  !
  ALLOCATE (stk%sysMcrs%v(stk%sysMcrs%nnz),stk%sysMcrs%i(stk%sysMcrs%neq+1))
  ALLOCATE (stk%sysMcrs%j(stk%sysMcrs%nnz))

  ALLOCATE (stk%rhsMcrs%v(stk%rhsMcrs%nnz),stk%rhsMcrs%i(stk%rhsMcrs%neq+1))
  ALLOCATE (stk%rhsMcrs%j(stk%rhsMcrs%nnz))
  
  row = 0  ! row in sysM and rhsM
  !
  ! Loop over subdomains
  !
  DO isd=1,nosd
    !
    ! Function source points
    !
    
    r=0 ! row in T and G matrices

    DO i = 1,subdomain(isd)%nnodes
      myWall = subdomain(isd)%BCidList(i)
      r=r+1
      !
      ! Loop over all three small systems
      !
      DO een=1,neq
        bc = eqn(een)%boundary(myWall)%known 
        !
        ! Use this equation only, if function is unknown
        !
        IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN

          row=row+1
          stk%sysMcrs%i(row) = nnzSys + 1
          stk%rhsMcrs%i(row) = nnzRHS + 1          
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
                nnzRHS = nnzRHS + 1
                stk%rhsMcrs%v(nnzRHS) = val
                stk%rhsMcrs%j(nnzRHS) = col
              END IF


              IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                stk%sysMcrs%v(nnzSys) = val
                stk%sysMcrs%j(nnzSys) = col
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
                nnzSys = nnzSys + 1
                stk%sysMcrs%v(nnzSys) = val
                stk%sysMcrs%j(nnzSys) = col
              END IF
              IF (bc.EQ.iFlux) THEN
                nnzRHS = nnzRHS + 1
                stk%rhsMcrs%v(nnzRHS) = val
                stk%rhsMcrs%j(nnzRHS) = col
              END IF

            END DO
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
          stk%sysMcrs%i(row) = nnzSys + 1
          stk%rhsMcrs%i(row) = nnzRHS + 1           

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
                nnzRHS = nnzRHS + 1
                stk%rhsMcrs%v(nnzRHS) = val
                stk%rhsMcrs%j(nnzRHS) = col
              END IF


              IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
                nnzSys = nnzSys + 1
                stk%sysMcrs%v(nnzSys) = val
                stk%sysMcrs%j(nnzSys) = col
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
                nnzSys = nnzSys + 1
                stk%sysMcrs%v(nnzSys) = val
                stk%sysMcrs%j(nnzSys) = col
              END IF
              IF (bc.EQ.iFlux) THEN
                nnzRHS = nnzRHS + 1
                stk%rhsMcrs%v(nnzRHS) = val
                stk%rhsMcrs%j(nnzRHS) = col
              END IF

            END DO
          END DO


        END IF
      END DO ! en        
    END DO ! q source point
  END DO ! subdomains
 
END 

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesBigSystemCount(nnzSys,nnzRHS)
  USE mMesh
  USE mEqns
  USE mPar

  IMPLICIT NONE
  INTEGER i,en,isp,isd,myWall,bc
  INTEGER j,myWallC,een

  INTEGER(8) nnzSys,nnzRHS

  nnzSys = 0
  nnzRhs = 0
  !
  ! Loop over subdomains
  !
  DO isd=1,nosd
    !
    ! Function source points
    !
    DO i = 1,subdomain(isd)%nnodes
      isp = subdomain(isd)%nodeList(i)  ! source point in the small system
      myWall = subdomain(isd)%BCidList(i)
      !
      ! Loop over all three small systems
      !
      DO een=1,neq
        bc = eqn(een)%boundary(myWall)%known 
        !
        ! Use this equation only, if function is unknown
        !
        IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
          !
          ! T matrix values
          !
          DO j = 1,subdomain(isd)%nnodes ! loop over T matrix columns
            myWallC = subdomain(isd)%BCidList(j)
            DO en=1,neq
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzRHS = nnzRHS + 1
              IF (bc.EQ.iFlux)     nnzSys = nnzSys + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1
            END DO
          END DO
          !
          ! G matrix values
          !
          DO j = 1,subdomain(isd)%nqnodes ! loop over G matrix columns
            myWallC = subdomain(isd)%qBCidList(j)
            DO en=1,neq          
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzSys = nnzSys + 1
              IF (bc.EQ.iFlux)     nnzRHS = nnzRHS + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1
            END DO
          END DO
        END IF ! I need this equation
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
      DO een=1,neq          
        bc = eqn(een)%boundary(myWall)%known
        !
        ! Use this equation only, if flux is unknown
        !
        IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
          !
          ! T matrix values
          !
          DO j = 1,subdomain(isd)%nnodes            
            myWallC = subdomain(isd)%BCidList(j)
            DO en=1,neq
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzRHS = nnzRHS + 1
              IF (bc.EQ.iFlux)     nnzSys = nnzSys + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1
            END DO
          END DO
          !
          ! G matrix values
          !
          DO j = 1,subdomain(isd)%nqnodes
            myWallC = subdomain(isd)%qBCidList(j)
            DO en=1,neq
              bc = eqn(en)%boundary(myWallC)%known

              IF (bc.EQ.iFunction) nnzSys = nnzSys + 1
              IF (bc.EQ.iFlux)     nnzRHS = nnzRHS + 1
              IF (bc.EQ.iContact)  nnzSys = nnzSys + 1
            END DO
          END DO
        END IF
      END DO ! en        
    END DO ! q source point
  END DO ! subdomains

END
  
  

!
!     ------------------------------------------------------------------
!
SUBROUTINE StokesBigSystemXB()
  USE mMesh
  USE mEqns
  USE mPar

  IMPLICIT NONE
  INTEGER i,en,isp,ispL,isd,myWall,bc

  ALLOCATE (stk%boundary(nofw))
  ALLOCATE (stk%initial(nofw))
  ALLOCATE (stk%u(3*nnodes))
  ALLOCATE (stk%q(3*nqnodes))
  ALLOCATE (stk%col(3*nnodes))
  ALLOCATE (stk%qcol(3*nqnodes))


  WRITE(stk%name,'(3A)') TRIM(eqn(1)%name),TRIM(eqn(2)%name),TRIM(eqn(3)%name)



  !
  !   Shift column values from 3 small SLE (eqn%) to one large (stk%)
  !              
  DO isd=1,nosd
    !
    ! Function source points
    !
    DO i = 1,subdomain(isd)%nnodes
      isp = subdomain(isd)%nodeList(i)  ! source point in the small system
      myWall = subdomain(isd)%BCidList(i)
      !
      ! Loop over all three small systems
      !
      stk%nx = 0
      stk%nb = 0
      stk%neq = 0
      DO en=1,neq
        ispL = isp + (en-1)*nnodes ! source point number in the large system       
        bc = eqn(en)%boundary(myWall)%known 
        IF (bc.EQ.iFunction) THEN           
          ! function is known
          stk%col(ispL) = eqn(en)%col(isp) + stk%nb
        ELSE IF (bc.EQ.iFlux.OR.bc.EQ.iContact) THEN
          ! function is unknown 
          stk%col(ispL) = eqn(en)%col(isp) + stk%nx
        END IF
        stk%nx = stk%nx + eqn(en)%nx ! number of unknowns (=cols in sysM)
        stk%nb = stk%nb + eqn(en)%nb ! number of equations (=rows in sys matrix)
        stk%neq = stk%neq + eqn(en)%neq ! number of rows in rhs vector  (=cols in rhs matrix)           
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
      stk%nx = 0
      stk%nb = 0
      stk%neq = 0
      DO en=1,neq        
        ispL = isp + (en-1)*nqnodes! source point number in the large system        
        bc = eqn(en)%boundary(myWall)%known
        IF (bc.EQ.iFunction.OR.bc.EQ.iContact) THEN
          ! flux is unknown
          stk%qcol(ispL) = eqn(en)%qcol(isp) + stk%nx
        ELSE IF (bc.EQ.iFlux) THEN
          ! flux is known
          stk%qcol(ispL) = eqn(en)%qcol(isp) + stk%nb
        END IF
        stk%nx = stk%nx + eqn(en)%nx ! number of unknowns (=cols in sysM)
        stk%nb = stk%nb + eqn(en)%nb ! number of equations (=rows in sys matrix)
        stk%neq = stk%neq + eqn(en)%neq ! number of rows in rhs vector  (=cols in rhs matrix)    
      END DO ! en        
    END DO ! q source point
  END DO ! subdomains
   
  
!
! Alocate vectors of unknowns and knowns
!
  ALLOCATE (stk%x(stk%nx))
  ALLOCATE (stk%Pivot(stk%nx))
  ALLOCATE (stk%b(stk%nb))
!
! Report to log file
!
  CALL WriteToLog("")
  CALL WriteToLog("--------------------")
  WRITE (parLogTekst,'(A,A)') "Equation = ", TRIM(stk%name)
  CALL WriteToLog(parLogTekst)
  WRITE (parLogTekst,'(A,I0,A,I0,A)') "System matrix = ",stk%neq, " x ",stk%nx
  CALL WriteToLog(parLogTekst)
  WRITE (parLogTekst,'(A,I0,A,I0,A)') "RHS    matrix = ",stk%neq, " x ",stk%nb
  CALL WriteToLog(parLogTekst)


END
  
  


!
!     ------------------------------------------------------------------
!
subroutine sdFormIntegralMatricesStokes()
!
    USE mMesh
    USE mEqns
    USE mPar
    USE mCommon
    IMPLICIT NONE

    REAL(rk), ALLOCATABLE :: rowTxx(:), rowTxy(:), rowTxz(:), rowTyy(:), rowTyz(:), rowTzz(:)
    REAL(rk), ALLOCATABLE :: rowGxx(:),rowGxy(:),rowGxz(:),rowGyy(:),rowGyz(:),rowGzz(:)

    REAL(rk), ALLOCATABLE :: rowprTx(:), rowprTy(:), rowprTz(:)
    REAL(rk), ALLOCATABLE :: rowprGx(:),rowprGy(:),rowprGz(:)    
 
    REAL(rk) ura,ura0,err,cptime
    INTEGER isd,irow,icol,i,j,row
!
!   Allocate space for integrals for a single source point
!
    
    ALLOCATE (rowTxx(nnodes), rowTxy(nnodes), rowTxz(nnodes), rowTyy(nnodes), rowTyz(nnodes), rowTzz(nnodes))
    ALLOCATE (rowGxx(nqnodes),rowGxy(nqnodes),rowGxz(nqnodes),rowGyy(nqnodes),rowGyz(nqnodes),rowGzz(nqnodes))

    ALLOCATE (rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes))
    ALLOCATE (rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes))
  
    ura0=cptime(0.0_rk)
    
!
!   Loop over subdomains
!
    DO isd=1,nosd
        row=0

        WRITE (parLogTekst,'(A,A,A)') "Stokes working on: ",TRIM(subdomain(isd)%name)," function source points."
        CALL WriteToLog(parLogTekst)

!
!       Allocate memory for matrices
!                
        CALL AllocateStokesMatrices(isd) 
!
!       function source points
!       set error indicator        
        subdomain(isd)%fspCerr = 0.0_rk
!
        DO i = 1,subdomain(isd)%nnodes 
            irow = subdomain(isd)%nodeList(i)  ! source point
!
!           Integrate row
!
            row=row+1

            CALL sdIntRowUStokes(isd,irow, & ! irow = source point ID, samo po subdomainu
                                 rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz, &
                                 rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz, &
                                 rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz, &
                                 err)
!
!           Distribute to columns
!
            DO j = 1,subdomain(isd)%nnodes
              icol = subdomain(isd)%nodeList(j)  ! where to take column value
              subdomain(isd)%TxxMat(row,j)=rowTxx(icol)
              subdomain(isd)%TyyMat(row,j)=rowTyy(icol)
              subdomain(isd)%TzzMat(row,j)=rowTzz(icol)
              subdomain(isd)%TxyMat(row,j)=rowTxy(icol)
              subdomain(isd)%TxzMat(row,j)=rowTxz(icol)
              subdomain(isd)%TyzMat(row,j)=rowTyz(icol)
              ! pressure
              subdomain(isd)%prTxMat(row,j)=rowprTx(icol)
              subdomain(isd)%prTyMat(row,j)=rowprTy(icol)
              subdomain(isd)%prTzMat(row,j)=rowprTz(icol)
            END DO
!  
            DO j = 1,subdomain(isd)%nqnodes
              icol = subdomain(isd)%qnodeList(j)  ! where to take column value
              subdomain(isd)%GxxMat(row,j)=rowGxx(icol)
              subdomain(isd)%GyyMat(row,j)=rowGyy(icol)
              subdomain(isd)%GzzMat(row,j)=rowGzz(icol)
              subdomain(isd)%GxyMat(row,j)=rowGxy(icol)
              subdomain(isd)%GxzMat(row,j)=rowGxz(icol)
              subdomain(isd)%GyzMat(row,j)=rowGyz(icol)
              ! pressure
              subdomain(isd)%prGxMat(row,j)=rowprGx(icol)
              subdomain(isd)%prGyMat(row,j)=rowprGy(icol)
              subdomain(isd)%prGzMat(row,j)=rowprGz(icol)              
            END DO
!
!           Remember error
!
            subdomain(isd)%fspCerr = MAX(subdomain(isd)%fspCerr,err)
!  
        END DO ! u source point
!  
        
!  
        WRITE (parLogTekst,'(A,A,A)') "Stokes working on: ",TRIM(subdomain(isd)%name)," flux source points."
        CALL WriteToLog(parLogTekst)
        subdomain(isd)%QspCerr = 0.0_rk 
!
!       Flux source points
!
        DO i = 1,subdomain(isd)%nqnodes
          irow = subdomain(isd)%qnodeList(i)  ! source point
!
!         Integrate row
!
            row=row+1
            CALL sdIntRowQStokes(isd,irow, & ! irow = source point ID, samo po subdomainu
                                 rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz, &
                                 rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz, &
                                 rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz, &
                                 err)     
!
!           Distribute to columns
!
            DO j = 1,subdomain(isd)%nnodes
               icol = subdomain(isd)%nodeList(j)  ! where to take column value
               subdomain(isd)%TxxMat(row,j)=rowTxx(icol)
               subdomain(isd)%TyyMat(row,j)=rowTyy(icol)
               subdomain(isd)%TzzMat(row,j)=rowTzz(icol)
               subdomain(isd)%TxyMat(row,j)=rowTxy(icol)
               subdomain(isd)%TxzMat(row,j)=rowTxz(icol)
               subdomain(isd)%TyzMat(row,j)=rowTyz(icol)
              ! pressure
               subdomain(isd)%prTxMat(row,j)=rowprTx(icol)
               subdomain(isd)%prTyMat(row,j)=rowprTy(icol)
               subdomain(isd)%prTzMat(row,j)=rowprTz(icol)               
             END DO
!  
             DO j = 1,subdomain(isd)%nqnodes
               icol = subdomain(isd)%qnodeList(j)  ! where to take column value
               subdomain(isd)%GxxMat(row,j)=rowGxx(icol)
               subdomain(isd)%GyyMat(row,j)=rowGyy(icol)
               subdomain(isd)%GzzMat(row,j)=rowGzz(icol)
               subdomain(isd)%GxyMat(row,j)=rowGxy(icol)
               subdomain(isd)%GxzMat(row,j)=rowGxz(icol)
               subdomain(isd)%GyzMat(row,j)=rowGyz(icol)
              ! pressure
               subdomain(isd)%prGxMat(row,j)=rowprGx(icol)
               subdomain(isd)%prGyMat(row,j)=rowprGy(icol)
               subdomain(isd)%prGzMat(row,j)=rowprGz(icol)                
             END DO
                      
!
!           Remember error
!
            subdomain(isd)%QspCerr = MAX(subdomain(isd)%QspCerr,err)
!  
!  
          END DO ! q source point
    
    END DO ! subdomains
!  
    ura=cptime(ura0)
    WRITE (parLogTekst,'(A,F14.4,A4)') "Integration time = ", ura/60.0_rk," min"    
    CALL WriteToLog(parLogTekst)

    DEALLOCATE (rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz)
    DEALLOCATE (rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz)
    DEALLOCATE (rowprTx,rowprTy,rowprTz)
    DEALLOCATE (rowprGx,rowprGy,rowprGz)
    
    end subroutine



!
!     ------------------------------------------------------------------
!
SUBROUTINE sdIntRowQStokes(isd,irow, & 
    rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz, &
    rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz, &
    rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz, &
    err)    


    USE mMesh
    USE Triangle ! provides tipw
    USE mPar
    USE mCommon
    IMPLICIT NONE
  
    INTEGER isd ! subdomain number
    INTEGER irow ! source point
    INTEGER ie,i,j,IntRecDepth,isrc,jj,iWall
    REAL(rk) sx,sy,sz
    REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    REAL(rk) err,shpf

    REAL(rk) cx,cy,cz,osemPi

    REAL(rk) rowTxx(nnodes), rowTxy(nnodes), rowTxz(nnodes), rowTyy(nnodes), rowTyz(nnodes), rowTzz(nnodes)
    REAL(rk) rowGxx(nqnodes),rowGxy(nqnodes),rowGxz(nqnodes),rowGyy(nqnodes),rowGyz(nqnodes),rowGzz(nqnodes)

    REAL(rk) rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes)
    REAL(rk) rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes)
    
    REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
    REAL(rk), ALLOCATABLE :: Txx(:),Txy(:),Txz(:),Tyy(:),Tyz(:),Tzz(:)
    
    REAL(rk) prGx,prGy,prGz
    REAL(rk), ALLOCATABLE ::  prTx(:),prTy(:),prTz(:)


    osemPi = ATAN(1.0)*4.0_rk*8.0_rk
    
!   Set integrals to zero
    rowTxx = 0.0_rk
    rowTxy = 0.0_rk
    rowTxz = 0.0_rk
    rowTyy = 0.0_rk
    rowTyz = 0.0_rk
    rowTzz = 0.0_rk

    rowprTx = 0.0_rk
    rowprTy = 0.0_rk
    rowprTz = 0.0_rk

    rowGxx = 0.0_rk
    rowGxy = 0.0_rk
    rowGxz = 0.0_rk
    rowGyy = 0.0_rk
    rowGyz = 0.0_rk
    rowGzz = 0.0_rk

  
!
!   source point location
!
    sx=qnode(irow)%x(1)
    sy=qnode(irow)%x(2)
    sz=qnode(irow)%x(3)
  
!
!   Integrate over boudnary elements in subdomain
!
    DO jj=1,subdomain(isd)%nofw ! loop over walls in a subdomain
       iWall = subdomain(isd)%loWalls(jj) ! current wall number
  
          DO ie=1,nelem
            IF (element(ie)%bcid.EQ.iWall) THEN ! loop over elements in a wall
  
  
!
!               Find if element is singular
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
!               Perform intergration over a triangle
!               (due to recursive nature of integration, set to zero here)
                ALLOCATE (Txx(element(ie)%nno))
                ALLOCATE (Txy(element(ie)%nno))
                ALLOCATE (Txz(element(ie)%nno))
                ALLOCATE (Tyy(element(ie)%nno))
                ALLOCATE (Tyz(element(ie)%nno))
                ALLOCATE (Tzz(element(ie)%nno))
                ALLOCATE (prTx(element(ie)%nno))
                ALLOCATE (prTy(element(ie)%nno))
                ALLOCATE (prTz(element(ie)%nno))
                Gxx=0.0_rk
                Gxy=0.0_rk
                Gxz=0.0_rk
                Gyy=0.0_rk
                Gyz=0.0_rk
                Gzz=0.0_rk
                prGx=0.0_rk
                prGy=0.0_rk
                prGz=0.0_rk
                DO j=1,element(ie)%nno
                    Txx(j)=0.0_rk
                    Txy(j)=0.0_rk
                    Txz(j)=0.0_rk
                    Tyy(j)=0.0_rk
                    Tyz(j)=0.0_rk
                    Tzz(j)=0.0_rk
                    prTx(j)=0.0_rk
                    prTy(j)=0.0_rk
                    prTz(j)=0.0_rk
                END DO
!
!               Element corners
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
  
                CALL Triangle_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                    subdomain(isd)%normMul(jj)*element(ie)%normal(1), &
                    subdomain(isd)%normMul(jj)*element(ie)%normal(2), & 
                    subdomain(isd)%normMul(jj)*element(ie)%normal(3), &
                    element(ie)%area,isrc,IntRecDepth, &                           
                    Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
  
              ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad
  

                    x4=node(element(ie)%con(4))%x(1)
                    y4=node(element(ie)%con(4))%x(2)
                    z4=node(element(ie)%con(4))%x(3)
       
                    CALL Quad_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sx,sy,sz,isrc,subdomain(isd)%normMul(jj), &
                               Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
 
              ELSE
                CALL WriteToLog("Error :: Element type not supported!")
              END IF
  !
  !           Distribute results to integral matrices
  !
              DO j=1,element(ie)%nno
                rowTxx(element(ie)%con(j)) = rowTxx(element(ie)%con(j)) + Txx(j)
                rowTxy(element(ie)%con(j)) = rowTxy(element(ie)%con(j)) + Txy(j)
                rowTxz(element(ie)%con(j)) = rowTxz(element(ie)%con(j)) + Txz(j)
                rowTyy(element(ie)%con(j)) = rowTyy(element(ie)%con(j)) + Tyy(j)
                rowTyz(element(ie)%con(j)) = rowTyz(element(ie)%con(j)) + Tyz(j)
                rowTzz(element(ie)%con(j)) = rowTzz(element(ie)%con(j)) + Tzz(j)

                rowprTx(element(ie)%con(j)) = rowprTx(element(ie)%con(j)) + prTx(j)
                rowprTy(element(ie)%con(j)) = rowprTy(element(ie)%con(j)) + prTy(j)
                rowprTz(element(ie)%con(j)) = rowprTz(element(ie)%con(j)) + prTz(j)
            END DO

            rowGxx(ie) = Gxx
            rowGxy(ie) = Gxy
            rowGxz(ie) = Gxz
            rowGyy(ie) = Gyy
            rowGyz(ie) = Gyz
            rowGzz(ie) = Gzz

            rowprGx(ie) = prGx
            rowprGy(ie) = prGy
            rowprGz(ie) = prGz

            DEALLOCATE (Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz)
  
            END IF
        END DO ! nbelem in wall
    END DO ! walls in subdomain
  
!
!       Find c parameter : must be 0.5 since source point is in the middle of element
!


!    Check solution : integral n_k T_ijk naj bi bil ..., Pozdrikis en. 7.2.20
!                                        = 0 če i != j  in če je na ravni steni, če vogal ali rob ni nič.
 
    cx = 0.0_rk
    cy = 0.0_rk
    cz = 0.0_rk
    DO i=1,nnodes
      cx = cx - ( rowTxy(i) )
      cy = cy - ( rowTyz(i) )
      cz = cz - ( rowTxz(i) )
    END DO
 
    err=MAX(err,ABS(cx))
    err=MAX(err,ABS(cy))
    err=MAX(err,ABS(cz))
 
!    if (cx>5.0E-5_rk.or.cy>5.0E-5_rk.or.cz>5.0E-5_rk) then
!      write(*,*) "Qsp: TxyTxzTyz integrals verify",cx,cy,cz
!      write(*,*) "sp=",sx,sy,sz 
!    end if

!
!   Get singular integral : integral n_k T_ijk naj bi bil ..., Pozdrikis en. 7.2.20
!                                       = 4*pi če i=j  
    cx = 0.0_rk
    cy = 0.0_rk
    cz = 0.0_rk
    DO i=1,nnodes
        cx = cx - ( rowTxx(i) )
        cy = cy - ( rowTyy(i) )
        cz = cz - ( rowTzz(i) )
    END DO

    err=MAX(err,ABS(cx - 0.5_rk*osemPi))
    err=MAX(err,ABS(cy - 0.5_rk*osemPi))
    err=MAX(err,ABS(cz - 0.5_rk*osemPi))
    
!    write(*,*) "Qsp: cji diagonala",cx/osemPi,cy/osemPi,cz/osemPi,err

!   dodamo del ! ja k vsaki izmed 3 oz. 4 tock
!   element v katerem je izvorna tocka
    ie=irow
!   shape functions for center point (constant flux approximation)
    IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
        shpf=1.0_rk/3.0_rk
    ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad
        shpf=0.25_rk
    END IF
!   razdelim c med vozlisca
    DO j=1,element(ie)%nno
        rowTxx(element(ie)%con(j)) = rowTxx(element(ie)%con(j)) + shpf * cx
        rowTyy(element(ie)%con(j)) = rowTyy(element(ie)%con(j)) + shpf * cy
        rowTzz(element(ie)%con(j)) = rowTzz(element(ie)%con(j)) + shpf * cz
    END DO
  

!    cx=0.0_rk
!    cy=0.0_rk
!    cz=0.0_rk
!    DO i=1,nnodes
!      cx = cx +  rowTxx(i)
!      cy = cy +  rowTyy(i)
!      cz = cz +  rowTzz(i)
!    END DO
!    print *,"qSP: Vsota Tii_jev ",cx,cy,cz
!    if (abs(cx)>1.0E-10) stop
!    if (abs(cy)>1.0E-10) stop
!    if (abs(cz)>1.0E-10) stop

! Correct pressure Q


    
!    cx = 0.0_rk
!    cy = 0.0_rk
!    cz = 0.0_rk
!    DO i=1,nnodes
!      cx = cx - ( rowprTx(i) )!* node(i)%x(3) )
!      cy = cy - ( rowprTy(i) )!* node(i)%x(2) )
!      cz = cz - ( rowprTz(i) )!* node(i)%x(3) )
!    END DO
!
!    write(*,*) "PR Qsp: cji diagonala",cx/osemPi,cy/osemPi,cz/osemPi
!
!   Add c to diagonal of H matrix since c+Hu=Gq, za ta del nisem prepričan, da mora biti res.
!
   ! print *,rowprTx(irow),rowprTy(irow),rowprTz(irow)
   ! rowprTx(irow) = rowprTx(irow) + cx 
   ! rowprTy(irow) = rowprTy(irow) + cy 
   ! rowprTz(irow) = rowprTz(irow) + cz 
   ! print *,rowprTx(irow),rowprTy(irow),rowprTz(irow)

END SUBROUTINE
  
  

!
!     ------------------------------------------------------------------
!
subroutine sdIntRowUStokes(isd,irow, &
                           rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz, &
                           rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz, &
                           rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz, &
                           err)
    USE mMesh
    USE Triangle ! provides tipw
    USE mPar
    USE mCommon
    IMPLICIT NONE

    INTEGER isd ! subdomain number
    INTEGER irow ! source point
    INTEGER ie,i,j,IntRecDepth,isrc,jj,iWall

    REAL(rk) sx,sy,sz
    REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    REAL(rk) err,multi

    REAL(rk) cx,cy,cz,osemPi

    REAL(rk) rowTxx(nnodes), rowTxy(nnodes), rowTxz(nnodes), rowTyy(nnodes), rowTyz(nnodes), rowTzz(nnodes)
    REAL(rk) rowGxx(nqnodes),rowGxy(nqnodes),rowGxz(nqnodes),rowGyy(nqnodes),rowGyz(nqnodes),rowGzz(nqnodes)

    REAL(rk) rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes)
    REAL(rk) rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes)
    
    REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
    REAL(rk), ALLOCATABLE :: Txx(:),Txy(:),Txz(:),Tyy(:),Tyz(:),Tzz(:)
    
    REAL(rk) prGx,prGy,prGz
    REAL(rk), ALLOCATABLE ::  prTx(:),prTy(:),prTz(:)


    osemPi = ATAN(1.0)*4.0_rk*8.0_rk

!   Set integrals to zero
    rowTxx = 0.0_rk
    rowTxy = 0.0_rk
    rowTxz = 0.0_rk
    rowTyy = 0.0_rk
    rowTyz = 0.0_rk
    rowTzz = 0.0_rk

    rowprTx = 0.0_rk
    rowprTy = 0.0_rk
    rowprTz = 0.0_rk

    rowGxx = 0.0_rk
    rowGxy = 0.0_rk
    rowGxz = 0.0_rk
    rowGyy = 0.0_rk
    rowGyz = 0.0_rk
    rowGzz = 0.0_rk
!
!   source point location
!
    sx=node(irow)%x(1)
    sy=node(irow)%x(2)
    sz=node(irow)%x(3)
    
!
!   Integrate over boudnary elements in subdomain
!
    DO jj=1,subdomain(isd)%nofw ! loop over walls in a subdomain
        iWall = subdomain(isd)%loWalls(jj) ! current wall number

        DO ie=1,nelem
            IF (element(ie)%bcid.EQ.iWall) THEN ! loop over elements in a wall
!
!               Find if element is singular
!
                isrc=0
                DO j=1,element(ie)%nno
                    IF (irow.EQ.element(ie)%con(j)) isrc=j
                END DO
!
!               Perform intergration over a triangle
!               (due to recursive nature of integration, set to zero here)
                ALLOCATE (Txx(element(ie)%nno))
                ALLOCATE (Txy(element(ie)%nno))
                ALLOCATE (Txz(element(ie)%nno))
                ALLOCATE (Tyy(element(ie)%nno))
                ALLOCATE (Tyz(element(ie)%nno))
                ALLOCATE (Tzz(element(ie)%nno))
                ALLOCATE (prTx(element(ie)%nno))
                ALLOCATE (prTy(element(ie)%nno))
                ALLOCATE (prTz(element(ie)%nno))
                Gxx=0.0_rk
                Gxy=0.0_rk
                Gxz=0.0_rk
                Gyy=0.0_rk
                Gyz=0.0_rk
                Gzz=0.0_rk
                prGx=0.0_rk
                prGy=0.0_rk
                prGz=0.0_rk
                DO j=1,element(ie)%nno
                    Txx(j)=0.0_rk
                    Txy(j)=0.0_rk
                    Txz(j)=0.0_rk
                    Tyy(j)=0.0_rk
                    Tyz(j)=0.0_rk
                    Tzz(j)=0.0_rk
                    prTx(j)=0.0_rk
                    prTy(j)=0.0_rk
                    prTz(j)=0.0_rk
                END DO
!
!               Element corners
!
                x1=node(element(ie)%con(1))%x(1)  
                y1=node(element(ie)%con(1))%x(2)  
                z1=node(element(ie)%con(1))%x(3)  
                
                x2=node(element(ie)%con(2))%x(1)  
                y2=node(element(ie)%con(2))%x(2)  
                z2=node(element(ie)%con(2))%x(3)  
                
                x3=node(element(ie)%con(3))%x(1)  
                y3=node(element(ie)%con(3))%x(2)  
                z3=node(element(ie)%con(3))%x(3)  

                IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
!
!                   Set recursion depth for singular triangles
!
                    IntRecDepth=parTriRecur

                    CALL Triangle_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                            subdomain(isd)%normMul(jj)*element(ie)%normal(1), &
                            subdomain(isd)%normMul(jj)*element(ie)%normal(2), & 
                            subdomain(isd)%normMul(jj)*element(ie)%normal(3), &
                            element(ie)%area,isrc,IntRecDepth, &                           
                            Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )

               ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad
     
                  x4=node(element(ie)%con(4))%x(1)
                  y4=node(element(ie)%con(4))%x(2)
                  z4=node(element(ie)%con(4))%x(3)
     
                  CALL Quad_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sx,sy,sz,isrc,subdomain(isd)%normMul(jj), &
                             Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
     
                ELSE
                  CALL WriteToLog("Error :: Element type not supported!")
                END IF
!
!               Distribute results to integral matrices
!

                DO j=1,element(ie)%nno
                    rowTxx(element(ie)%con(j)) = rowTxx(element(ie)%con(j)) + Txx(j)
                    rowTxy(element(ie)%con(j)) = rowTxy(element(ie)%con(j)) + Txy(j)
                    rowTxz(element(ie)%con(j)) = rowTxz(element(ie)%con(j)) + Txz(j)
                    rowTyy(element(ie)%con(j)) = rowTyy(element(ie)%con(j)) + Tyy(j)
                    rowTyz(element(ie)%con(j)) = rowTyz(element(ie)%con(j)) + Tyz(j)
                    rowTzz(element(ie)%con(j)) = rowTzz(element(ie)%con(j)) + Tzz(j)

                    rowprTx(element(ie)%con(j)) = rowprTx(element(ie)%con(j)) + prTx(j)
                    rowprTy(element(ie)%con(j)) = rowprTy(element(ie)%con(j)) + prTy(j)
                    rowprTz(element(ie)%con(j)) = rowprTz(element(ie)%con(j)) + prTz(j)
                END DO

                rowGxx(ie) = Gxx
                rowGxy(ie) = Gxy
                rowGxz(ie) = Gxz
                rowGyy(ie) = Gyy
                rowGyz(ie) = Gyz
                rowGzz(ie) = Gzz

                rowprGx(ie) = prGx
                rowprGy(ie) = prGy
                rowprGz(ie) = prGz

                DEALLOCATE (Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz)
            END IF
        END DO ! nbelem in wall
    END DO ! walls in subdomain


!
!   Set integration error
!        
    err = -9.99_rk
!
!    write(*,*) ""
!    write(*,*) "source point",sx,sy,sz    


!
!   Check solution : integral n_k T_ijk naj bi bil ..., Pozdrikis en. 7.2.20
!                                       = 0 če i != j  in če je na ravni steni, če vogal ali rob ni nič.
!!
!    cx = 0.0_rk
!    cy = 0.0_rk
!    cz = 0.0_rk
!    DO i=1,nnodes
!      cx = cx - ( rowTxy(i) )
!      cy = cy - ( rowTxz(i) )
!      cz = cz - ( rowTyz(i) )
!    END DO
!
!    write(*,*) "TxyTxzTyz integrals verify",cx,cy,cz
!
!
!    cx = 0.0_rk
!    cy = 0.0_rk
!    cz = 0.0_rk
!    DO i=1,nqnodes
!      cx = cx - ( rowGxy(i) )
!      cy = cy - ( rowGxz(i) )
!      cz = cz - ( rowGyz(i) )
!    END DO
!
!    write(*,*) "GxyGxzGyz integrals verify",cx,cy,cz
!

!
!   Get singular integral : integral n_k T_ijk naj bi bil ..., Pozdrikis en. 7.2.20
!                                       = 4*pi če i=j  in če je na ravni steni ??      
    cx = 0.0_rk
    cy = 0.0_rk
    cz = 0.0_rk
    DO i=1,nnodes
      cx = cx - ( rowTxx(i) )
      cy = cy - ( rowTyy(i) )
      cz = cz - ( rowTzz(i) )
    END DO

    !write(*,*) "Usp: cji diagonala",cx/osemPi,cy/osemPi,cz/osemPi
!
!   Add c to diagonal of H matrix since c+Hu=Gq, za ta del nisem prepričan, da mora biti res.
!
    rowTxx(irow) = rowTxx(irow) + cx 
    rowTyy(irow) = rowTyy(irow) + cy 
    rowTzz(irow) = rowTzz(irow) + cz 

!
!   Check solution : integral n_i G_ij naj bi bil nič za vse j, Pozdrikis en. 7.2.19
!
    cx=0.0_rk
    cy=0.0_rk
    cz=0.0_rk
    DO i=1,nelem
        multi=1.0_rk
        DO jj=1,subdomain(isd)%nofw ! find my wall in list of subdomain walls
            IF (element(i)%bcid.EQ.subdomain(isd)%loWalls(jj)) THEN
                multi = subdomain(isd)%normMul(jj)
            END IF
        END DO
        cx = cx - rowGxx(i)*element(i)%normal(1)*multi &
                - rowGxy(i)*element(i)%normal(2)*multi &
                - rowGxz(i)*element(i)%normal(3)*multi
        cy = cy - rowGxy(i)*element(i)%normal(1)*multi &
                - rowGyy(i)*element(i)%normal(2)*multi &
                - rowGyz(i)*element(i)%normal(3)*multi
        cz = cz - rowGxz(i)*element(i)%normal(1)*multi &
                - rowGyz(i)*element(i)%normal(2)*multi &
                - rowGzz(i)*element(i)%normal(3)*multi
    END DO

    err=MAX(err,ABS(cx))
    err=MAX(err,ABS(cy))
    err=MAX(err,ABS(cz))

    !write(*,*) "Usp: G integrals verify",cx,cy,cz

!    cx=0.0_rk
!    cy=0.0_rk
!    cz=0.0_rk
!    DO i=1,nnodes
!      cx = cx +  rowTxx(i)
!      cy = cy +  rowTyy(i)
!      cz = cz +  rowTzz(i)
!    END DO
!    print *,"uSP: Vsota Tii_jev ",cx,cy,cz



! Correct pressure
    
!    cx = 0.0_rk
!    cy = 0.0_rk
!    cz = 0.0_rk
!    DO i=1,nnodes
!      cx = cx - ( rowprTx(i) )
!      cy = cy - ( rowprTy(i) )
!      cz = cz - ( rowprTz(i) )
!    END DO

    !write(*,*) "PR Usp: cji diagonala",cx/osemPi,cy/osemPi,cz/osemPi
!
!   Add c to diagonal of H matrix since c+Hu=Gq, za ta del nisem prepričan, da mora biti res.
!
    !print *,rowprTx(irow),rowprTy(irow),rowprTz(irow)
    !rowprTx(irow) = rowprTx(irow) + cx 
    !rowprTy(irow) = rowprTy(irow) + cy 
    !rowprTz(irow) = rowprTz(irow) + cz 
    !print *,rowprTx(irow),rowprTy(irow),rowprTz(irow)

!    cx = 0.0_rk
!    cy = 0.0_rk
!    cz = 0.0_rk
!    DO i=1,nnodes
!      cx = cx - ( rowprTx(i) )
!      cy = cy - ( rowprTy(i) )
!      cz = cz - ( rowprTz(i) )
!    END DO
!
!!    write(*,*) "PR Usp: cji verify",cx/osemPi,cy/osemPi,cz/osemPi
!
!    cx = 0.0_rk
!    cy = 0.0_rk
!    cz = 0.0_rk
!    DO i=1,nqnodes
!      cx = cx - ( rowprGx(i) )
!      cy = cy - ( rowprGy(i) )
!      cz = cz - ( rowprGz(i) )
!    END DO

    !write(*,*) "PR Qsp: cji diagonala",cx/osemPi,cy/osemPi,cz/osemPi
!
!   Add c to diagonal of H matrix since c+Hu=Gq, za ta del nisem prepričan, da mora biti res.
!
    !rowprGx(irow) = rowprGx(irow) + cx 
    !rowprGy(irow) = rowprGy(irow) + cy 
    !rowprGz(irow) = rowprGz(irow) + cz 



end subroutine


!
!     ******************************************************************
!
RECURSIVE SUBROUTINE Triangle_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area,isrc,n, &
                                        Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
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
    REAL(rk) area,areaS
    INTEGER isrc,isrc1,isrc2,isrc3,isrc4,n

    ! Results
    REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
    REAL(rk) prGx,prGy,prGz
    REAL(rk) Txx(3),Txy(3),Txz(3),Tyy(3),Tyz(3),Tzz(3)
    REAL(rk) prTx(3),prTy(3),prTz(3)

    IF (isrc.EQ.0.OR.n.LT.1) THEN
!
!       non-singular integrals ( or singular at the end of recursion)
!
        CALL Triangle_3DStokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area, &
                                  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
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
        CALL        Triangle_Area(x1,y1,z1,x4,y4,z4,x6,y6,z6,areaS)
        CALL Triangle_StokesInt(x1,y1,z1,x4,y4,z4,x6,y6,z6,sx,sy,sz,nx,ny,nz,areaS,isrc1,n-1, &
                                  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )

        CALL        Triangle_Area(x4,y4,z4,x2,y2,z2,x5,y5,z5,areaS)
        CALL Triangle_StokesInt(x4,y4,z4,x2,y2,z2,x5,y5,z5,sx,sy,sz,nx,ny,nz,areaS,isrc2,n-1, &
                                  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )    

        CALL        Triangle_Area(x6,y6,z6,x5,y5,z5,x3,y3,z3,areaS)
        CALL Triangle_StokesInt(x6,y6,z6,x5,y5,z5,x3,y3,z3,sx,sy,sz,nx,ny,nz,areaS,isrc3,n-1, &
                                  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )    
        

        CALL        Triangle_Area(x4,y4,z4,x5,y5,z5,x6,y6,z6,areaS)
        CALL Triangle_StokesInt(x4,y4,z4,x5,y5,z5,x6,y6,z6,sx,sy,sz,nx,ny,nz,areaS,isrc4,n-1, &
                                  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )    
        

    END IF

END

!
!     ******************************************************************
!
SUBROUTINE Triangle_3DStokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area, &
                                Gxx,Gxy,Gxz,Gyy,Gyz,Gzz, &
                                prGx,prGy,prGz, &
                                Txx,Txy,Txz,Tyy,Tyz,Tzz, &
                                prTx,prTy,prTz )
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
    REAL(rk) area,w,w1,w2,w3

    ! Kernels
    REAL(rk) G(3,3) ! Stokeslet
    REAL(rk) T(3,3) ! Stresslet times normal
    REAL(rk) prG(3) ! Pressure 
    REAL(rk) prT(3) ! Pressure 

    ! Results
    REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
    REAL(rk) prGx,prGy,prGz
    REAL(rk) Txx(3),Txy(3),Txz(3),Tyy(3),Tyz(3),Tzz(3)
    REAL(rk) prTx(3),prTy(3),prTz(3)

    INTEGER i
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
!       3D Stokes flow fundamental solutions
!
        CALL StokesKernel(sx,sy,sz,x,y,z,nx,ny,nz,G,T,prG,prT)
!
!       Sum up integral
!
        w = tipw%w(i) * area
        w1 = w * tipw%l1(i)
        w2 = w * tipw%l2(i)
        w3 = w * tipw%l3(i)

        Gxx = Gxx + w * G(1,1)  ! constant interpolation of flux
        Gxy = Gxy + w * G(1,2)  ! constant interpolation of flux
        Gxz = Gxz + w * G(1,3)  ! constant interpolation of flux

        Gyy = Gyy + w * G(2,2)  ! constant interpolation of flux
        Gyz = Gyz + w * G(2,3)  ! constant interpolation of flux
        
        Gzz = Gzz + w * G(3,3)  ! constant interpolation of flux

        prGx = prGx + w * prG(1)  ! constant interpolation of flux               
        prGy = prGy + w * prG(2)  ! constant interpolation of flux               
        prGz = prGz + w * prG(3)  ! constant interpolation of flux               

        Txx(1) = Txx(1) + w1 * T(1,1)   ! linear
        Txx(2) = Txx(2) + w2 * T(1,1)   ! interpolation of
        Txx(3) = Txx(3) + w3 * T(1,1)   ! function

        Txy(1) = Txy(1) + w1 * T(1,2)   ! linear
        Txy(2) = Txy(2) + w2 * T(1,2)   ! interpolation of
        Txy(3) = Txy(3) + w3 * T(1,2)   ! function

        Txz(1) = Txz(1) + w1 * T(1,3)   ! linear
        Txz(2) = Txz(2) + w2 * T(1,3)   ! interpolation of
        Txz(3) = Txz(3) + w3 * T(1,3)   ! function

        Tyy(1) = Tyy(1) + w1 * T(2,2)   ! linear
        Tyy(2) = Tyy(2) + w2 * T(2,2)   ! interpolation of
        Tyy(3) = Tyy(3) + w3 * T(2,2)   ! function

        Tyz(1) = Tyz(1) + w1 * T(2,3)   ! linear
        Tyz(2) = Tyz(2) + w2 * T(2,3)   ! interpolation of
        Tyz(3) = Tyz(3) + w3 * T(2,3)   ! function

        Tzz(1) = Tzz(1) + w1 * T(3,3)   ! linear
        Tzz(2) = Tzz(2) + w2 * T(3,3)   ! interpolation of
        Tzz(3) = Tzz(3) + w3 * T(3,3)   ! function

        prTx(1) = prTx(1) + w1 * prT(1)   ! linear
        prTx(2) = prTx(2) + w2 * prT(1)   ! interpolation of
        prTx(3) = prTx(3) + w3 * prT(1)   ! function

        prTy(1) = prTy(1) + w1 * prT(2)   ! linear
        prTy(2) = prTy(2) + w2 * prT(2)   ! interpolation of
        prTy(3) = prTy(3) + w3 * prT(2)   ! function

        prTz(1) = prTz(1) + w1 * prT(3)   ! linear
        prTz(2) = prTz(2) + w2 * prT(3)   ! interpolation of
        prTz(3) = prTz(3) + w3 * prT(3)   ! function

    END DO ! integration points

END
!
!     ******************************************************************
!
SUBROUTINE StokesKernel(sx,sy,sz,fx,fy,fz,nx,ny,nz,G,T,prG,prT)
!
!     Calculates BEM kernel
!
!     ******************************************************************
    USE mCommon
    IMPLICIT NONE

    REAL(rk) sx,sy,sz ! source point
    REAL(rk) nx,ny,nz ! unit normal on boundary element surface
    REAL(rk) fx,fy,fz ! field (integration) point    
    ! Kernels
    REAL(rk) G(3,3) ! Stokeslet
    REAL(rk) T(3,3) ! Stresslet times normal
    REAL(rk) prG(3) ! Pressure 
    REAL(rk) prT(3) ! Pressure 

    REAL(rk) rx,ry,rz,ir,ir2,ir3,c
    REAL(rk) rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz

    rx = fx - sx
    ry = fy - sy
    rz = fz - sz

    rxx = rx * rx
    rxy = rx * ry
    rxz = rx * rz

    ryx = rxy
    ryy = ry * ry
    ryz = ry * rz

    rzx = rxz
    rzy = ryz
    rzz = rz * rz

    ir2 = 1.0_rk / (rxx + ryy + rzz)
    ir = dsqrt(ir2)
    ir3 = ir2 * ir

!
!     STOKESLET
!
    G(1,1) = ir * (1.0_rk + rxx * ir2 )
    G(1,2) =                rxy * ir3
    G(1,3) =                rxz * ir3

    G(2,1) = G(1,2)
    G(2,2) = ir * (1.0_rk + ryy * ir2 )
    G(2,3) =                ryz * ir3

    G(3,1) = G(1,3)
    G(3,2) = G(2,3)
    G(3,3) = ir * (1.0_rk + rzz * ir2 )

!
!   Stresslet times normal
!            
    c = - 6.0_rk * ir2 * ir3 * ( rx * nx + ry * ny + rz * nz ) ! tu ze mnozim z normalo

    T(1,1) = c * rxx
    T(1,2) = c * rxy
    T(1,3) = c * rxz

    T(2,1) = T(1,2)
    T(2,2) = c * ryy
    T(2,3) = c * ryz

    T(3,1) = T(1,3)
    T(3,2) = T(2,3)
    T(3,3) = c * rzz

!
!   Pressure G
!            
    c = 2.0_rk * ir3
    prG(1) = c * rx
    prG(2) = c * ry
    prG(3) = c * rz

!
!   Pressure T times normal
!            
    c = 4.0_rk * ir3
    prT(1) = c * ( ( - 1.0_rk + 3.0_rk * rxx * ir2 ) * nx + &
                                3.0_rk * rxy * ir2   * ny + & 
                                3.0_rk * rxz * ir2   * nz )
    prT(2) = c * (              3.0_rk * ryx * ir2   * nx + &
                   ( - 1.0_rk + 3.0_rk * ryy * ir2 ) * ny + & 
                                3.0_rk * ryz * ir2   * nz )
    prT(3) = c * (              3.0_rk * rzx * ir2   * nx + &
                                3.0_rk * rzy * ir2   * ny + &
                   ( - 1.0_rk + 3.0_rk * rzz * ir2 ) * nz )

    END


!
!      subroutine sgf_3d_fs 
!        +
!        +   (Iopt
!        +   ,x,y,z
!        +   ,x0,y0,z0
!        +   ,Gxx,Gxy,Gxz
!        +   ,Gyx,Gyy,Gyz
!        +   ,Gzx,Gzy,Gzz
!        +   ,px,py,pz
!        +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
!        +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
!        +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
!        +   )
!   
!   c-----------------------------------------
!   c FDLIB, BEMLIB
!   c
!   c Copyright by C. Pozrikidis, 1999
!   c All rights reserved.
!   c
!   c This program is to be used only under the
!   c stipulations of the licensing agreement.
!   c----------------------------------------
!   
!   c---------------------------------------
!   c Free-space Green's function: Stokeslet
!   c
!   c Pozrikidis (1992, p. 23)
!   c
!   c Iopt =  1 generates only the Green's function
!   c      ne 1 generates the Green's function,
!   c           pressure, and stress
!   c---------------------------------------
!   
!         Implicit Double Precision (a-h,o-z)
!   
!         dx = x-x0
!         dy = y-y0
!         dz = z-z0
!   
!         dxx = dx*dx
!         dxy = dx*dy
!         dxz = dx*dz
!         dyy = dy*dy
!         dyz = dy*dz
!         dzz = dz*dz
!   
!         r  = Dsqrt(dxx+dyy+dzz)
!         r3 = r**3
!   
!         ri  = 1.0D0/r
!         ri3 = 1.0D0/r3
!   
!         Gxx = ri + dxx*ri3
!         Gxy =      dxy*ri3
!         Gxz =      dxz*ri3
!         Gyy = ri + dyy*ri3
!         Gyz =      dyz*ri3
!         Gzz = ri + dzz*ri3
!   
!         Gyx = Gxy
!         Gzx = Gxz
!         Gzy = Gyz
!   
!         If(Iopt.eq.1) Go to 99
!   
!   c--------------
!   c stress tensor
!   c--------------
!   
!         cf = -6.0D0/r**5
!   
!         Txxx = dxx*dx * cf
!         Txxy = dxy*dx * cf
!         Txxz = dxz*dx * cf
!         Tyxy = dyy*dx * cf
!         Tyxz = dyz*dx * cf
!         Tzxz = dzz*dx * cf
!   
!         Txyx = Txxy
!         Txyy = Tyxy
!         Txyz = Tyxz
!         Tyyy = dyy*dy * cf
!         Tyyz = dyz*dy * cf
!         Tzyz = dzz*dy * cf
!   
!         Txzx = Txxz
!         Txzy = Tyxz
!         Txzz = Tzxz
!         Tyzy = dyy*dz * cf
!         Tyzz = dyz*dz * cf
!         Tzzz = dzz*dz * cf
!   
!   c---------
!   c pressure
!   c---------
!   
!         cf = 2.0D0*ri3
!         px = dx * cf
!         py = dy * cf
!         pz = dz * cf
!   
!     99  Continue
!   
!   c-----
!   c Done
!   c-----
!   
!         Return
!         End
!   


!
!     ------------------------------------------------------------------
!
SUBROUTINE CoR_Integrals_Stokes()

      USE mMesh
      USE mPar
      USE mCommon

      IMPLICIT NONE

      INTEGER ierr,isd
      REAL(rk) uSize,qSize,tSize

!
!     Report sizes to log file
!
      tSize = 0.0_rk
      DO isd=1,nosd

        uSize = DBLE(subdomain(isd)%nsp)*DBLE(subdomain(isd)%nnodes)*8.0_rk/1024.0_rk/1024.0_rk
        WRITE (parLogTekst,'(A,A,I0,1X,I0,A,F10.2,A)') TRIM(subdomain(isd)%name)," Tij & pT matrix size = ",subdomain(isd)%nsp &
          ,subdomain(isd)%nnodes," = ",uSize," Mb"
        CALL WriteToLog(parLogTekst)

        qSize = DBLE(subdomain(isd)%nsp)*DBLE(subdomain(isd)%nqnodes)*8.0/1024.0/1024.0
        WRITE (parLogTekst,'(A,A,I0,1X,I0,A,F10.2,A)') TRIM(subdomain(isd)%name)," Gij & pG matrix size = ",subdomain(isd)%nsp &
          ,subdomain(isd)%nqnodes," = ",qSize," Mb"
        CALL WriteToLog(parLogTekst)

        tSize = tSize + 9.0_rk*(uSize + qSize)
      END DO

      WRITE (parLogTekst,'(A,F10.2,A)') "total size = ",tSize," Mb"
      CALL WriteToLog(parLogTekst)


!
!     Read
!
      CALL WriteToLog("Verifying integrals file!")
      CALL VerifyIntegralsFileStokes(ierr)
!
!     Are integrals on disk
!
      IF (ierr.NE.0) THEN
!
!       No
!
        CALL WriteToLog("Verify failed - computing integrals!")
        CALL sdFormIntegralMatricesStokes()
        IF (parWriteIntegrals.EQ.parYes) THEN
          CALL WriteToLog("Writing integrals to disk!")
          CALL WriteIntegralsToDiskStokes()
        END IF
        parfromDisk = parNo
      ELSE
!
!       Yes
!
        CALL WriteToLog("Verify successful - no need to compute integrals!")
        CALL WriteToLog("Reading integrals from disk!")
        CALL ReadIntegralsFileStokes(ierr)
        parfromDisk = parNo
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
SUBROUTINE WriteIntegralsToDiskStokes()
      USE mMesh
      USE mPar

      IMPLICIT NONE

      INTEGER isd,lun

      lun=12

      OPEN (lun,FILE=TRIM(parIntegralsFileName),FORM='UNFORMATTED',STATUS='UNKNOWN')

      WRITE (lun) parStokes
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
        CALL WrSMat(subdomain(isd)%TxxMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%TyyMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%TzzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%TxyMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%TxzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%TyzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)

        CALL WrSMat(subdomain(isd)%GxxMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%GyyMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%GzzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%GxyMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%GxzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%GyzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)

        ! pressure
        CALL WrSMat(subdomain(isd)%prTxMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%prTyMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL WrSMat(subdomain(isd)%prTzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)

        CALL WrSMat(subdomain(isd)%prGxMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%prGyMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL WrSMat(subdomain(isd)%prGzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)

      END DO


      CLOSE(lun)
end subroutine        


!
!     ------------------------------------------------------------------
!
    SUBROUTINE ReadIntegralsFileStokes(ierr)
      USE mMesh
      USE mPar

      IMPLICIT NONE

      INTEGER isd,lun,ierr,a,b,c,id

      ierr=1
      lun=12

      OPEN (lun,FILE=TRIM(parIntegralsFileName),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ (lun,ERR=10) id
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
!
!       Allocate memory for matrices
!                
        CALL AllocateStokesMatrices(isd) 
!
!       Read data from disk
!
        CALL RdSMat(subdomain(isd)%TxxMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%TyyMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%TzzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%TxyMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%TxzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%TyzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)

        CALL RdSMat(subdomain(isd)%GxxMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%GyyMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%GzzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%GxyMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%GxzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%GyzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)  

        ! pressure
        CALL RdSMat(subdomain(isd)%prTxMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%prTyMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)
        CALL RdSMat(subdomain(isd)%prTzMat,subdomain(isd)%nsp,subdomain(isd)%nnodes,lun)

        CALL RdSMat(subdomain(isd)%prGxMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%prGyMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)
        CALL RdSMat(subdomain(isd)%prGzMat,subdomain(isd)%nsp,subdomain(isd)%nqnodes,lun)        
        
      END DO

      ierr=0
10    CONTINUE
      CLOSE(lun)

      END SUBROUTINE
  

!
!     ------------------------------------------------------------------
!
    SUBROUTINE VerifyIntegralsFileStokes(ierr)
      USE mMesh
      USE mPar

      IMPLICIT NONE

      INTEGER isd,lun,ierr,a,b,c,id

      ierr=1
      lun=12

      OPEN (lun,FILE=TRIM(parIntegralsFileName),FORM='UNFORMATTED',STATUS='OLD',ERR=10)

      READ (lun,ERR=10) id
      IF (id.NE.parStokes) GOTO 10
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
subroutine AllocateStokesMatrices(isd)
    USE mMesh
    USE mEqns
    USE mCommon
    IMPLICIT NONE

    INTEGER isd
 
    ALLOCATE (subdomain(isd)%TxxMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%TyyMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%TzzMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%TxyMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%TxzMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%TyzMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))

    ALLOCATE (subdomain(isd)%GxxMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%GyyMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%GzzMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%GxyMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%GxzMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%GyzMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))

!   Pressure
    ALLOCATE (subdomain(isd)%prTxMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%prTyMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))
    ALLOCATE (subdomain(isd)%prTzMat(subdomain(isd)%nsp,subdomain(isd)%nnodes))

    ALLOCATE (subdomain(isd)%prGxMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%prGyMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))
    ALLOCATE (subdomain(isd)%prGzMat(subdomain(isd)%nsp,subdomain(isd)%nqnodes))

end subroutine

!
!     ------------------------------------------------------------------
!
subroutine dellocateStokesMatrices(isd)
  USE mMesh
  USE mEqns
  USE mCommon
  IMPLICIT NONE

  INTEGER isd

  DEALLOCATE (subdomain(isd)%TxxMat)
  DEALLOCATE (subdomain(isd)%TyyMat)
  DEALLOCATE (subdomain(isd)%TzzMat)
  DEALLOCATE (subdomain(isd)%TxyMat)
  DEALLOCATE (subdomain(isd)%TxzMat)
  DEALLOCATE (subdomain(isd)%TyzMat)

  DEALLOCATE (subdomain(isd)%GxxMat)
  DEALLOCATE (subdomain(isd)%GyyMat)
  DEALLOCATE (subdomain(isd)%GzzMat)
  DEALLOCATE (subdomain(isd)%GxyMat)
  DEALLOCATE (subdomain(isd)%GxzMat)
  DEALLOCATE (subdomain(isd)%GyzMat)

end subroutine


!
!     ------------------------------------------------------------------
!
subroutine dellocateStokesPressureMatrices(isd)
  USE mMesh
  USE mEqns
  USE mCommon
  IMPLICIT NONE

  INTEGER isd

  ! pressure
  DEALLOCATE (subdomain(isd)%prTxMat)
  DEALLOCATE (subdomain(isd)%prTyMat)
  DEALLOCATE (subdomain(isd)%prTzMat)

  DEALLOCATE (subdomain(isd)%prGxMat)
  DEALLOCATE (subdomain(isd)%prGyMat)
  DEALLOCATE (subdomain(isd)%prGzMat)

end subroutine


subroutine getStokesDVtest()
  USE mMesh
  USE mCommon
  USE mEqns
  IMPLICIT NONE
  real(rk) SourcePoint(3),rx,ry,rz,err,p,diff

  SourcePoint(1)=0.99_rk
  SourcePoint(2)=0.4_rk
  SourcePoint(3)=0.6_rk
  diff = subdomain(1)%diff

  call GetStokesDomainValue(SourcePoint,eqn(1)%u,eqn(2)%u,eqn(3)%u,eqn(1)%q,eqn(2)%q,eqn(3)%q,rx,ry,rz,p,diff,err)
  print *,"jre"
  print *,rx,ry,rz,p,err
  stop
end subroutine

subroutine getStokesDV(SourcePoint,res,p,ierr)
  USE mMesh
  USE mCommon
  USE mEqns
  IMPLICIT NONE
  integer ierr
  real(rk) SourcePoint(3),res(3),err,p
  call GetStokesDomainValue(SourcePoint,eqn(1)%u,eqn(2)%u,eqn(3)%u,eqn(1)%q,eqn(2)%q,eqn(3)%q,res(1),res(2),res(3), &
                            p,subdomain(1)%diff,err)
  
  if (err.gt.1E-2) then
    ierr = -1
  else
    ierr = 1 ! ierr = stevilka subdomaina (pri Stokes samo 1 subdomain)
  end if


end subroutine


!
! --------------------------------------
!
SUBROUTINE GetStokesDomainValue(SourcePoint,ux,uy,uz,qx,qy,qz,rx,ry,rz,p,diff,err)

  USE mMesh
  USE Triangle ! provides tipw
  USE mPar
  USE mCommon
  IMPLICIT NONE

  REAL(rk) ux(nnodes),qx(nqnodes)
  REAL(rk) uy(nnodes),qy(nqnodes)
  REAL(rk) uz(nnodes),qz(nqnodes)
  REAL(rk) rx,ry,rz,p,diff
  REAL(rk) SourcePoint(3)

  INTEGER isd ! subdomain number
  INTEGER ie,i,j,IntRecDepth,isrc,jj,iWall

  REAL(rk) sx,sy,sz
  REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  REAL(rk) err,multi

  REAL(rk) cx,cy,cz,osemPi

  REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
  REAL(rk), ALLOCATABLE :: Txx(:),Txy(:),Txz(:),Tyy(:),Tyz(:),Tzz(:)

  REAL(rk) prGx,prGy,prGz
  REAL(rk), ALLOCATABLE ::  prTx(:),prTy(:),prTz(:)

  REAL(rk), ALLOCATABLE :: rowTxx(:), rowTxy(:), rowTxz(:), rowTyy(:), rowTyz(:), rowTzz(:)
  REAL(rk), ALLOCATABLE :: rowGxx(:),rowGxy(:),rowGxz(:),rowGyy(:),rowGyz(:),rowGzz(:)

  REAL(rk), ALLOCATABLE :: rowprTx(:),rowprTy(:),rowprTz(:)
  REAL(rk), ALLOCATABLE :: rowprGx(:),rowprGy(:),rowprGz(:)    

  !
  !   Allocate space for integrals for a single source point
  !

  ALLOCATE (rowTxx(nnodes), rowTxy(nnodes), rowTxz(nnodes), rowTyy(nnodes), rowTyz(nnodes), rowTzz(nnodes))
  ALLOCATE (rowGxx(nqnodes),rowGxy(nqnodes),rowGxz(nqnodes),rowGyy(nqnodes),rowGyz(nqnodes),rowGzz(nqnodes))

  ALLOCATE (rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes))
  ALLOCATE (rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes))


  osemPi = ATAN(1.0)*4.0_rk*8.0_rk

  !   Set integrals to zero
  rowTxx = 0.0_rk
  rowTxy = 0.0_rk
  rowTxz = 0.0_rk
  rowTyy = 0.0_rk
  rowTyz = 0.0_rk
  rowTzz = 0.0_rk

  rowprTx = 0.0_rk
  rowprTy = 0.0_rk
  rowprTz = 0.0_rk

  rowGxx = 0.0_rk
  rowGxy = 0.0_rk
  rowGxz = 0.0_rk
  rowGyy = 0.0_rk
  rowGyz = 0.0_rk
  rowGzz = 0.0_rk
  !
  !   source point location
  !
  sx=SourcePoint(1)
  sy=SourcePoint(2)
  sz=SourcePoint(3)
  !
  !   Stokes flow, only 1 subdomain supported
  !
  isd = 1


  !
  !   Integrate over boudnary elements in subdomain
  !
  DO jj=1,subdomain(isd)%nofw ! loop over walls in a subdomain
    iWall = subdomain(isd)%loWalls(jj) ! current wall number

    DO ie=1,nelem
      IF (element(ie)%bcid.EQ.iWall) THEN ! loop over elements in a wall
        !
        ! Element is not singular, since the source point is in the domain
        !
        isrc=0
        !
        ! Perform intergration over a triangle
        ! (due to recursive nature of integration, set to zero here)
        ALLOCATE (Txx(element(ie)%nno))
        ALLOCATE (Txy(element(ie)%nno))
        ALLOCATE (Txz(element(ie)%nno))
        ALLOCATE (Tyy(element(ie)%nno))
        ALLOCATE (Tyz(element(ie)%nno))
        ALLOCATE (Tzz(element(ie)%nno))
        ALLOCATE (prTx(element(ie)%nno))
        ALLOCATE (prTy(element(ie)%nno))
        ALLOCATE (prTz(element(ie)%nno))
        Gxx=0.0_rk
        Gxy=0.0_rk
        Gxz=0.0_rk
        Gyy=0.0_rk
        Gyz=0.0_rk
        Gzz=0.0_rk
        prGx=0.0_rk
        prGy=0.0_rk
        prGz=0.0_rk
        DO j=1,element(ie)%nno
          Txx(j)=0.0_rk
          Txy(j)=0.0_rk
          Txz(j)=0.0_rk
          Tyy(j)=0.0_rk
          Tyz(j)=0.0_rk
          Tzz(j)=0.0_rk
          prTx(j)=0.0_rk
          prTy(j)=0.0_rk
          prTz(j)=0.0_rk
        END DO
        !
        !               Element corners
        !
        x1=node(element(ie)%con(1))%x(1)  
        y1=node(element(ie)%con(1))%x(2)  
        z1=node(element(ie)%con(1))%x(3)  

        x2=node(element(ie)%con(2))%x(1)  
        y2=node(element(ie)%con(2))%x(2)  
        z2=node(element(ie)%con(2))%x(3)  

        x3=node(element(ie)%con(3))%x(1)  
        y3=node(element(ie)%con(3))%x(2)  
        z3=node(element(ie)%con(3))%x(3)  

        IF (element(ie)%type.EQ.2) THEN ! 3 node trangle
        !
        ! Set recursion depth for singular triangles
        !
          IntRecDepth=parTriRecur

          CALL Triangle_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
           subdomain(isd)%normMul(jj)*element(ie)%normal(1), &
           subdomain(isd)%normMul(jj)*element(ie)%normal(2), & 
           subdomain(isd)%normMul(jj)*element(ie)%normal(3), &
           element(ie)%area,isrc,IntRecDepth, &                           
           Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )

        ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad
          
          x4=node(element(ie)%con(4))%x(1)
          y4=node(element(ie)%con(4))%x(2)
          z4=node(element(ie)%con(4))%x(3)
          
          CALL Quad_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sx,sy,sz,isrc,subdomain(isd)%normMul(jj), &
          Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )

        ELSE
          CALL WriteToLog("Error :: Element type not supported!")
        END IF
        !
        !               Distribute results to integral matrices
        !

        DO j=1,element(ie)%nno
          rowTxx(element(ie)%con(j)) = rowTxx(element(ie)%con(j)) + Txx(j)
          rowTxy(element(ie)%con(j)) = rowTxy(element(ie)%con(j)) + Txy(j)
          rowTxz(element(ie)%con(j)) = rowTxz(element(ie)%con(j)) + Txz(j)
          rowTyy(element(ie)%con(j)) = rowTyy(element(ie)%con(j)) + Tyy(j)
          rowTyz(element(ie)%con(j)) = rowTyz(element(ie)%con(j)) + Tyz(j)
          rowTzz(element(ie)%con(j)) = rowTzz(element(ie)%con(j)) + Tzz(j)

          rowprTx(element(ie)%con(j)) = rowprTx(element(ie)%con(j)) + prTx(j)
          rowprTy(element(ie)%con(j)) = rowprTy(element(ie)%con(j)) + prTy(j)
          rowprTz(element(ie)%con(j)) = rowprTz(element(ie)%con(j)) + prTz(j)
        END DO

        rowGxx(ie) = Gxx
        rowGxy(ie) = Gxy
        rowGxz(ie) = Gxz
        rowGyy(ie) = Gyy
        rowGyz(ie) = Gyz
        rowGzz(ie) = Gzz

        rowprGx(ie) = prGx
        rowprGy(ie) = prGy
        rowprGz(ie) = prGz

        DEALLOCATE (Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz)
      END IF
    END DO ! nbelem in wall
  END DO ! walls in subdomain


  !
  !   Set integration error
  !        
  err = -9.99_rk

  !
  !   Get singular integral : integral n_k T_ijk naj bi bil ..., Pozdrikis en. 7.2.20
  !                                       = 4*pi če i=j  in če je na ravni steni ??      
  cx = 0.0_rk
  cy = 0.0_rk
  cz = 0.0_rk
  DO i=1,nnodes
    cx = cx - ( rowTxx(i) )
    cy = cy - ( rowTyy(i) )
    cz = cz - ( rowTzz(i) )
  END DO


  err=MAX(err,ABS(cx/osemPi-1.0_rk))
  err=MAX(err,ABS(cy/osemPi-1.0_rk))
  err=MAX(err,ABS(cz/osemPi-1.0_rk))

  !write(*,*) "cji diagonala",cx/osemPi,cy/osemPi,cz/osemPi,err
  !
  !   Add c to diagonal of H matrix since c+Hu=Gq, za ta del nisem prepričan, da mora biti res.
  !
  !rowTxx(irow) = rowTxx(irow) + cx 
  !rowTyy(irow) = rowTyy(irow) + cy 
  !rowTzz(irow) = rowTzz(irow) + cz 

  !
  !   Check solution : integral n_i G_ij naj bi bil nič za vse j, Pozdrikis en. 7.2.19
  !
  cx=0.0_rk
  cy=0.0_rk
  cz=0.0_rk
  DO i=1,nelem
    multi=1.0_rk
    DO jj=1,subdomain(isd)%nofw ! find my wall in list of subdomain walls
      IF (element(i)%bcid.EQ.subdomain(isd)%loWalls(jj)) THEN
        multi = subdomain(isd)%normMul(jj)
      END IF
    END DO
    cx = cx - rowGxx(i)*element(i)%normal(1)*multi &
    - rowGxy(i)*element(i)%normal(2)*multi &
    - rowGxz(i)*element(i)%normal(3)*multi
    cy = cy - rowGxy(i)*element(i)%normal(1)*multi &
    - rowGyy(i)*element(i)%normal(2)*multi &
    - rowGyz(i)*element(i)%normal(3)*multi
    cz = cz - rowGxz(i)*element(i)%normal(1)*multi &
    - rowGyz(i)*element(i)%normal(2)*multi &
    - rowGzz(i)*element(i)%normal(3)*multi
  END DO

  err=MAX(err,ABS(cx))
  err=MAX(err,ABS(cy))
  err=MAX(err,ABS(cz))

!  write(*,*) "G integrals verify",cx,cy,cz
!  print *,err

!    cx=0.0_rk
!    cy=0.0_rk
!    cz=0.0_rk
!    DO i=1,nnodes
!      cx = cx +  rowTxx(i)
!      cy = cy +  rowTyy(i)
!      cz = cz +  rowTzz(i)
!    END DO
!    print *,"uSP: Vsota Tii_jev ",cx,cy,cz

!  s1 = 0.0_rk
!  s2 = 0.0_rk
!  s3 = 0.0_rk
!  DO i=1,nnodes
!    s1 = s1 + rowprTx(i)
!    s2 = s2 + rowprTy(i)
!    s3 = s3 + rowprTz(i)
!  END DO
!  print *,s1,s2,s3
!
!  s1 = 0.0_rk
!  s2 = 0.0_rk
!  s3 = 0.0_rk 
!  DO i=1,nqnodes
!
!    s1 = s1 + rowprGx(i)
!    s2 = s2 + rowprGy(i)
!    s3 = s3 + rowprGz(i)
!
!  END DO
!  print *,s1,s2,s3
!

  !
  !  Calculate velocity & pressure values
  !
  rx = 0.0_rk
  ry = 0.0_rk
  rz = 0.0_rk
  p = 0.0_rk

  DO i=1,nnodes
    rx = rx + rowTxx(i) * ux(i)
    rx = rx + rowTxy(i) * uy(i)
    rx = rx + rowTxz(i) * uz(i)

    ry = ry + rowTxy(i) * ux(i)
    ry = ry + rowTyy(i) * uy(i)
    ry = ry + rowTyz(i) * uz(i)

    rz = rz + rowTxz(i) * ux(i)
    rz = rz + rowTyz(i) * uy(i)
    rz = rz + rowTzz(i) * uz(i)

    p = p + diff * ( rowprTx(i) * ux(i) + rowprTy(i) * uy(i) + rowprTz(i) * uz(i) )

  END DO
  
  DO i=1,nqnodes

    rx = rx + rowGxx(i) * qx(i)
    rx = rx + rowGxy(i) * qy(i)
    rx = rx + rowGxz(i) * qz(i)

    ry = ry + rowGxy(i) * qx(i)
    ry = ry + rowGyy(i) * qy(i)
    ry = ry + rowGyz(i) * qz(i)

    rz = rz + rowGxz(i) * qx(i)
    rz = rz + rowGyz(i) * qy(i)
    rz = rz + rowGzz(i) * qz(i)    

    p = p - ( rowprGx(i) * qx(i) + rowprGy(i) * qy(i) + rowprGz(i) * qz(i) )

  END DO


  rx = - rx/osemPi
  ry = - ry/osemPi
  rz = - rz/osemPi  
  p  = -  p/osemPi

  DEALLOCATE (rowTxx,rowTxy,rowTxz,rowTyy,rowTyz,rowTzz)
  DEALLOCATE (rowGxx,rowGxy,rowGxz,rowGyy,rowGyz,rowGzz)
  DEALLOCATE (rowprTx,rowprTy,rowprTz)
  DEALLOCATE (rowprGx,rowprGy,rowprGz)

  !print *,"p=",p
  !stop

end subroutine




!
! -----------------------------------------------------------------------------
!
SUBROUTINE Quad_StokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,isrc,multi, &
  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
  !
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
    REAL(rk) xp,yp,zp ! source point
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 ! quad vertexes
  
  
    REAL(rk) xc0,yc0,zc0,xet,yet,zet,xks,yks,zks
    REAL(rk) xx1,yy1
    REAL(rk) ro,roth,th,ajac
    REAL(rk) pi2
    REAL(rk) eti,etj,eta1m,eta2m,eta1p,eta2p
    REAL(rk) anx,any,anz,anx1,any1,anz1
  
    REAL(rk) multi ! flips normals
  
    REAL(rk) fig4(4),w
    REAL(rk) al(4),fii(4),th0(4),th1(4)
    REAL(rk) ksi(5),eta(5) ! prve 4 za funkcijo, peta v sredini za fluks SP
    !
    !      integral divison
    REAL(rk) a,b,c,d,dex,gii,gij,d1,d2,d3,d4,d13max,d24max,minedge
    INTEGER idivXi,ndivXi,idivEt,ndivEt,ising
  
    ! Kernels
    REAL(rk) G(3,3) ! Stokeslet
    REAL(rk) T(3,3) ! Stresslet times normal
    REAL(rk) prG(3) ! Pressure 
    REAL(rk) prT(3) ! Pressure 
  
    ! Results
    REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
    REAL(rk) prGx,prGy,prGz
    REAL(rk) Txx(4),Txy(4),Txz(4),Tyy(4),Tyz(4),Tzz(4)
    REAL(rk) prTx(4),prTy(4),prTz(4)
  
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
    ! Set to zero
    !
    Gxx = 0.0_rk
    Gxy = 0.0_rk
    Gxz = 0.0_rk
    Gyy = 0.0_rk
    Gyz = 0.0_rk
    Gzz = 0.0_rk
    prGx = 0.0_rk
    prGy = 0.0_rk
    prGz = 0.0_rk
    DO i=1,4
      Txx(i) = 0.0_rk
      Txy(i) = 0.0_rk
      Txz(i) = 0.0_rk
      Tyy(i) = 0.0_rk
      Tyz(i) = 0.0_rk
      Tzz(i) = 0.0_rk
      prTx(i) = 0.0_rk
      prTy(i) = 0.0_rk
      prTz(i) = 0.0_rk
    END DO    
  
  
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
    ndivXi=3
    ndivEt=3

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
            ! Calculate Stokes kernel
            !
            CALL StokesKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,G,T,prG,prT)
            !
            !       Sum up integral
            !
            w = AJAC
            
            Gxx = Gxx + w * G(1,1)  ! constant interpolation of flux
            Gxy = Gxy + w * G(1,2)  ! constant interpolation of flux
            Gxz = Gxz + w * G(1,3)  ! constant interpolation of flux
    
            Gyy = Gyy + w * G(2,2)  ! constant interpolation of flux
            Gyz = Gyz + w * G(2,3)  ! constant interpolation of flux
            
            Gzz = Gzz + w * G(3,3)  ! constant interpolation of flux
    
            prGx = prGx + w * prG(1)  ! constant interpolation of flux               
            prGy = prGy + w * prG(2)  ! constant interpolation of flux               
            prGz = prGz + w * prG(3)  ! constant interpolation of flux               
            
            DO isip=1,4
              w = AJAC * FIG4(isip)
  
              Txx(isip) = Txx(isip) + w * T(1,1)   ! linear interpolation of function
              Txy(isip) = Txy(isip) + w * T(1,2)   ! linear interpolation of function
              Txz(isip) = Txz(isip) + w * T(1,3)   ! linear interpolation of function
  
              Tyy(isip) = Tyy(isip) + w * T(2,2)   ! linear interpolation of function
              Tyz(isip) = Tyz(isip) + w * T(2,3)   ! linear interpolation of function
              Tzz(isip) = Tzz(isip) + w * T(3,3)   ! linear interpolation of function
                        
              prTx(isip) = prTx(isip) + w * prT(1)   ! linear interpolation of function
              prTy(isip) = prTy(isip) + w * prT(2)   ! linear interpolation of function
              prTz(isip) = prTz(isip) + w * prT(3)   ! linear interpolation of function
            END DO
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
              ! Calculate Stokes kernel
              !
          CALL StokesKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,G,T,prG,prT)
          !
          !       Sum up integral
          !
          w = AJAC
  
          Gxx = Gxx + w * G(1,1)  ! constant interpolation of flux
          Gxy = Gxy + w * G(1,2)  ! constant interpolation of flux
          Gxz = Gxz + w * G(1,3)  ! constant interpolation of flux
  
          Gyy = Gyy + w * G(2,2)  ! constant interpolation of flux
          Gyz = Gyz + w * G(2,3)  ! constant interpolation of flux
  
          Gzz = Gzz + w * G(3,3)  ! constant interpolation of flux
  
          prGx = prGx + w * prG(1)  ! constant interpolation of flux               
          prGy = prGy + w * prG(2)  ! constant interpolation of flux               
          prGz = prGz + w * prG(3)  ! constant interpolation of flux               
  
          DO isip=1,4
            w = AJAC * FIG4(isip)
  
            Txx(isip) = Txx(isip) + w * T(1,1)   ! linear interpolation of function
            Txy(isip) = Txy(isip) + w * T(1,2)   ! linear interpolation of function
            Txz(isip) = Txz(isip) + w * T(1,3)   ! linear interpolation of function
  
            Tyy(isip) = Tyy(isip) + w * T(2,2)   ! linear interpolation of function
            Tyz(isip) = Tyz(isip) + w * T(2,3)   ! linear interpolation of function
            Tzz(isip) = Tzz(isip) + w * T(3,3)   ! linear interpolation of function
  
            prTx(isip) = prTx(isip) + w * prT(1)   ! linear interpolation of function
            prTy(isip) = prTy(isip) + w * prT(2)   ! linear interpolation of function
            prTz(isip) = prTz(isip) + w * prT(3)   ! linear interpolation of function
          END DO
  
        END DO
      END DO
      CONTINUE
    END IF
  
    END DO
  END DO
  
  END
  
  