!
!     ------------------------------------------------------------------
!
SUBROUTINE pressureStokesLIS()

   USE mMesh
   USE mEqns
   USE mPar
   USE mCommon
   IMPLICIT NONE
#include "lisf.h"    
 
   INTEGER row,isd,icol
   REAL(rk) osemPi,val
   REAL(rk), ALLOCATABLE :: vx(:),vy(:),vz(:)
   REAL(rk), ALLOCATABLE :: rowprTx(:),rowprTy(:),rowprTz(:)
   REAL(rk), ALLOCATABLE :: rowprGx(:),rowprGy(:),rowprGz(:) 
 
   LIS_VECTOR pressure 
   LIS_INTEGER is,ie
   REAL ts,te

   CALL CPU_TIME(ts)

   CALL WriteToLog("Solving for pressure!")

   osemPi = ATAN(1.0)*4.0_rk*8.0_rk
 

   ALLOCATE ( rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes) )
   ALLOCATE ( rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes) )
 
   ALLOCATE (vx(nqnodes),vy(nqnodes),vz(nqnodes))
 
   CALL InterpolateToQMesh(eqn(1)%u,vx)
   CALL InterpolateToQMesh(eqn(2)%u,vy)
   CALL InterpolateToQMesh(eqn(3)%u,vz)

   !
   ! Loop over subdomains (only 1 since Stokes!)
   !
   isd = 1
   !
   ! Divide qnodes between processors
   !
   call lis_vector_create(env%comm,pressure,env%ierr)      
   call lis_vector_set_size(pressure,0,subdomain(isd)%nqnodes,env%ierr)
   call lis_vector_set_all(0.0_rk,pressure,env%ierr)
   call lis_vector_get_range(pressure, is, ie, env%ierr)
   !
   ! Loop over my boundary elements
   !
   DO row = is,ie-1
      !
      !  Calculate integrals
      !               
      CALL sdIntRowQStokesPressure(isd,row,rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz)  
      !
      !  Calculate pressure, Ingber 1991, Surface pressure solution ...
      ! 
      val = 0.0_rk
      DO icol = 1,subdomain(isd)%nnodes
                
         val =  val +    subdomain(isd)%diff * ( rowprTx(icol) * ( eqn(1)%u(icol)-vx(row) ) & 
                                               + rowprTy(icol) * ( eqn(2)%u(icol)-vy(row) ) &
                                               + rowprTz(icol) * ( eqn(3)%u(icol)-vz(row) ) )  
      END DO 

      DO icol = 1,subdomain(isd)%nqnodes

         val =  val  - ( rowprGx(icol) * eqn(1)%q(icol) &
                       + rowprGy(icol) * eqn(2)%q(icol) &
                       + rowprGz(icol) * eqn(3)%q(icol) )         
     END DO      
     !
     ! Upload value to parallel vector
     !
     call lis_vector_set_value(LIS_INS_VALUE, row, val, pressure, env%ierr)
   END DO   
   call lis_vector_scale( - 2.0_rk /osemPi, pressure, env%ierr) ! *2 ker na ravni steni
   !
   ! Gather result from all processes
   !
   call lis_vector_gather(pressure,eqn(1)%p,env%ierr)
 
   DEALLOCATE ( vx,vy,vz )
   DEALLOCATE ( rowprTx,rowprTy,rowprTz )
   DEALLOCATE ( rowprGx,rowprGy,rowprGz )
   !
   !  Deallocate parallel vector
   !   
   call lis_vector_destroy(pressure,env%ierr)
   !
   !  Write to log
   ! 
   CALL CPU_TIME(te)
   WRITE (parLogTekst,'(A,F10.4)') "TIMER :: pressureStokesLIS [s] = ",te-ts
   CALL WriteToLog(parLogTekst)

 end subroutine
 



!
!     ------------------------------------------------------------------
!
SUBROUTINE sdIntRowQStokesPressure(isd,irow,rowprGx,rowprGy,rowprGz,rowprTx,rowprTy,rowprTz)    


   USE mMesh
   USE Triangle ! provides tipw
   USE mPar
   USE mCommon
   IMPLICIT NONE
 
   INTEGER isd ! subdomain number
   INTEGER irow ! source point
   INTEGER ie,j,IntRecDepth,isrc,jj,iWall
   REAL(rk) sx,sy,sz
   REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4

   REAL(rk) rowprTx(nnodes), rowprTy(nnodes), rowprTz(nnodes)
   REAL(rk) rowprGx(nqnodes),rowprGy(nqnodes),rowprGz(nqnodes)  
   
   REAL(rk) prGx,prGy,prGz
   REAL(rk), ALLOCATABLE ::  prTx(:),prTy(:),prTz(:)
   REAL(rk) cx,cy,cz

!   
!  Set integrals to zero
!   
   rowprTx = 0.0_rk
   rowprTy = 0.0_rk
   rowprTz = 0.0_rk

   rowprGx = 0.0_rk
   rowprGy = 0.0_rk
   rowprGz = 0.0_rk
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
            ALLOCATE (prTx(element(ie)%nno))
            ALLOCATE (prTy(element(ie)%nno))
            ALLOCATE (prTz(element(ie)%nno))

            prGx=0.0_rk
            prGy=0.0_rk
            prGz=0.0_rk
            DO j=1,element(ie)%nno
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
 !                Set recursion depth for singular triangles
 !
               IntRecDepth=parTriRecur
 
               CALL Triangle_StokesPressureInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz, &
                  subdomain(isd)%normMul(jj)*element(ie)%normal(1), &
                  subdomain(isd)%normMul(jj)*element(ie)%normal(2), & 
                  subdomain(isd)%normMul(jj)*element(ie)%normal(3), &
                  element(ie)%area,isrc,IntRecDepth, &                           
                  prGx,prGy,prGz,prTx,prTy,prTz )
 
            ELSE IF (element(ie)%type.EQ.3) THEN ! 4 node quad

               x4=node(element(ie)%con(4))%x(1)
               y4=node(element(ie)%con(4))%x(2)
               z4=node(element(ie)%con(4))%x(3)
      
               CALL Quad_StokesPressureInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sx,sy,sz,isrc,subdomain(isd)%normMul(jj), &
                                           prGx,prGy,prGz,prTx,prTy,prTz )
            ELSE
               CALL WriteToLog("Error :: Element type not supported!")
            END IF
 !
 !          Distribute results to integral matrices
 !
         if (isrc.eq.0) then ! ignore singular integral
            DO j=1,element(ie)%nno
               rowprTx(element(ie)%con(j)) = rowprTx(element(ie)%con(j)) + prTx(j)
               rowprTy(element(ie)%con(j)) = rowprTy(element(ie)%con(j)) + prTy(j)
               rowprTz(element(ie)%con(j)) = rowprTz(element(ie)%con(j)) + prTz(j)
            END DO
         
            rowprGx(ie) = prGx
            rowprGy(ie) = prGy
            rowprGz(ie) = prGz
         end if

            DEALLOCATE (prTx,prTy,prTz)
 
         END IF
       END DO ! nbelem in wall
   END DO ! walls in subdomain


   !
   !  Estimate singular integral
   !
   cx = 0.0_rk
   cy = 0.0_rk
   cz = 0.0_rk
   do j = 1,nnodes
      cx = cx + rowprTx(j)
      cy = cy + rowprTy(j)
      cz = cz + rowprTz(j)
   end do
   DO j=1,element(irow)%nno
      rowprTx(element(irow)%con(j)) = rowprTx(element(irow)%con(j)) - cx / element(irow)%nno      
      rowprTy(element(irow)%con(j)) = rowprTy(element(irow)%con(j)) - cy / element(irow)%nno   
      rowprTz(element(irow)%con(j)) = rowprTz(element(irow)%con(j)) - cz / element(irow)%nno         
   END DO   


   !cx = 0.0_rk
   !cy = 0.0_rk
   !cz = 0.0_rk
   !do j = 1,nqnodes
   !   cx = cx + rowprGx(j)
   !   cy = cy + rowprGy(j)
   !   cz = cz + rowprGz(j)
   !end do
   !rowprGx(irow) = - cx
   !rowprGy(irow) = - cy
   !owprGz(irow) = - cz
   !print *,irow,cx,cy,cz,qnode(irow)%x
   !print *,irow,cy,qnode(irow)%x,"G"
   !print *,irow,cz,qnode(irow)%x,"G"



END SUBROUTINE
 
 






!
!     ******************************************************************
!
RECURSIVE SUBROUTINE Triangle_StokesPressureInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area,isrc,n,&
                                                prGx,prGy,prGz,prTx,prTy,prTz )
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

!  Results
   REAL(rk) prGx,prGy,prGz
   REAL(rk) prTx(3),prTy(3),prTz(3)

   IF (isrc.EQ.0.OR.n.LT.1) THEN
!
!       non-singular integrals ( or singular at the end of recursion)
!
      CALL Triangle_3DStokesPressureInt(  x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area, &
                                          prGx,prGy,prGz,prTx,prTy,prTz )
   ELSE
!
!     singular integrals 
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
!     Detremine new singular integrals
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
!     Integrate four smaller pieces
!
      CALL Triangle_Area(x1,y1,z1,x4,y4,z4,x6,y6,z6,areaS)
      CALL Triangle_StokesPressureInt( x1,y1,z1,x4,y4,z4,x6,y6,z6,sx,sy,sz,nx,ny,nz,areaS,isrc1,n-1, &
                                       prGx,prGy,prGz,prTx,prTy,prTz )

      CALL Triangle_Area(x4,y4,z4,x2,y2,z2,x5,y5,z5,areaS)
      CALL Triangle_StokesPressureInt( x4,y4,z4,x2,y2,z2,x5,y5,z5,sx,sy,sz,nx,ny,nz,areaS,isrc2,n-1, &
                                       prGx,prGy,prGz,prTx,prTy,prTz )    

      CALL Triangle_Area(x6,y6,z6,x5,y5,z5,x3,y3,z3,areaS)
      CALL Triangle_StokesPressureInt( x6,y6,z6,x5,y5,z5,x3,y3,z3,sx,sy,sz,nx,ny,nz,areaS,isrc3,n-1, &
                                       prGx,prGy,prGz,prTx,prTy,prTz )    


      CALL Triangle_Area(x4,y4,z4,x5,y5,z5,x6,y6,z6,areaS)
      CALL Triangle_StokesPressureInt( x4,y4,z4,x5,y5,z5,x6,y6,z6,sx,sy,sz,nx,ny,nz,areaS,isrc4,n-1, &
                                       prGx,prGy,prGz,prTx,prTy,prTz )    
   END IF

END

!
!  ******************************************************************
!
SUBROUTINE Triangle_3DStokesPressureInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area,prGx,prGy,prGz,prTx,prTy,prTz )
!
!  Calculates integral over a triangle
!
!  ******************************************************************
   USE Triangle ! provides tipw
   USE mCommon
   IMPLICIT NONE

   REAL(rk) sx,sy,sz ! source point
   REAL(rk) nx,ny,nz ! unit normal on triangle surface
   REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
   REAL(rk) x,y,z    ! integration point
   REAL(rk) area,w,w1,w2,w3

!  Kernels
   REAL(rk) prG(3) ! Pressure 
   REAL(rk) prT(3) ! Pressure times normal

!  Results
   REAL(rk) prGx,prGy,prGz
   REAL(rk) prTx(3),prTy(3),prTz(3)

   INTEGER i
!
!     Loop over integration points
!
   DO i=1,tipw%n
!
!     Calculate point in R^3 space, where function must be evaluated
!
      x=tipw%l1(i)*x1+tipw%l2(i)*x2+tipw%l3(i)*x3
      y=tipw%l1(i)*y1+tipw%l2(i)*y2+tipw%l3(i)*y3
      z=tipw%l1(i)*z1+tipw%l2(i)*z2+tipw%l3(i)*z3       
!
!     3D Stokes flow fundamental solutions
!
      CALL StokesPressureKernel(sx,sy,sz,x,y,z,nx,ny,nz,prG,prT)
!
!     Sum up integral
!
      w = tipw%w(i) * area
      w1 = w * tipw%l1(i)
      w2 = w * tipw%l2(i)
      w3 = w * tipw%l3(i)

      prGx = prGx + w * prG(1)  ! constant interpolation of flux               
      prGy = prGy + w * prG(2)  ! constant interpolation of flux               
      prGz = prGz + w * prG(3)  ! constant interpolation of flux               

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
! -----------------------------------------------------------------------------
!
SUBROUTINE Quad_StokesPressureInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,isrc,multi,prGx,prGy,prGz,prTx,prTy,prTz )
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
   REAL(rk) prG(3) ! Pressure 
   REAL(rk) prT(3) ! Pressure times normal
  
    ! Results
   REAL(rk) prGx,prGy,prGz
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
   PI2=2.0_rk*PI
   !
   ! Set to zero
   !
   prGx = 0.0_rk
   prGy = 0.0_rk
   prGz = 0.0_rk
   DO i=1,4
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
            CALL StokesPressureKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,prG,prT)
            !
            !       Sum up integral
            !
            w = AJAC           
    
            prGx = prGx + w * prG(1)  ! constant interpolation of flux               
            prGy = prGy + w * prG(2)  ! constant interpolation of flux               
            prGz = prGz + w * prG(3)  ! constant interpolation of flux               
            
            DO isip=1,4
              w = AJAC * FIG4(isip)
                          
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
          CALL StokesPressureKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,prG,prT)
          !
          !       Sum up integral
          !
          w = AJAC
  

          prGx = prGx + w * prG(1)  ! constant interpolation of flux               
          prGy = prGy + w * prG(2)  ! constant interpolation of flux               
          prGz = prGz + w * prG(3)  ! constant interpolation of flux               
  
          DO isip=1,4
            w = AJAC * FIG4(isip)
  
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
  
  
!
!  ******************************************************************
!
SUBROUTINE StokesPressureKernel(sx,sy,sz,fx,fy,fz,nx,ny,nz,prG,prT)
!
!  Calculates Stokes BEM kernel for pressure
!
!  ******************************************************************
   USE mCommon
   IMPLICIT NONE

   REAL(rk) sx,sy,sz ! source point
   REAL(rk) nx,ny,nz ! unit normal on boundary element surface
   REAL(rk) fx,fy,fz ! field (integration) point    
   ! Kernels
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
