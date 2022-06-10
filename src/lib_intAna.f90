
!
! -----------------------------------------------------------------------------
!
SUBROUTINE Quad_StokesIntNumAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,isrc,multi, &
  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
  !
  !
  !     $: Integracija Robni 4 tockovni element s 4 tockovno geometrijo
  !
  ! -----------------------------------------------------------------------------
    USE mCommon
    IMPLICIT NONE
  
    INTEGER isrc
    !
    REAL(rk) xp,yp,zp ! source point
    REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 ! quad vertexes
  
    REAL(rk) multi ! flips normals
  
    ! Results
    REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
    REAL(rk) prGx,prGy,prGz
    REAL(rk) Txx(4),Txy(4),Txz(4),Tyy(4),Tyz(4),Tzz(4)
    REAL(rk) prTx(4),prTy(4),prTz(4)

    IF (isrc.EQ.0) THEN ! non singular integral
      CALL Quad_StokesIntNum(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,multi, &
      Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
    ELSE
      CALL Quad_StokesIntAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,isrc, &
      Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
    END IF

END SUBROUTINE


!
! -----------------------------------------------------------------------------------------
!
SUBROUTINE Quad_StokesIntAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,isrc, &
      Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
!
!     $: Integracija Robni 4 tockovni element s 4 tockovno geometrijo
!     Analiticno singularni
!
!
!     4 ----- 3
!     |       |
!     |   5   |
!     |       |
!     1 ----- 2    
!
! -----------------------------------------------------------------------------
  USE mCommon
  use singInt
  IMPLICIT NONE
  
  REAL(rk), INTENT(IN)  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  INTEGER,  INTENT(IN)  :: isrc

  REAL(rk), INTENT(OUT) ::  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
  REAL(rk), INTENT(OUT) ::  prGx,prGy,prGz
  REAL(rk), INTENT(OUT) ::  Txx(4),Txy(4),Txz(4),Tyy(4),Tyz(4),Tzz(4)
  REAL(rk), INTENT(OUT) ::  prTx(4),prTy(4),prTz(4)
  

  INTEGER i


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

  prGx = 0.0_rk
  prGy = 0.0_rk
  prGz = 0.0_rk

  call si_qua_stk(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,isrc)

END SUBROUTINE


!
!     ******************************************************************
!
SUBROUTINE Triangle_StokesIntAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area,isrc, &
          Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
!
!  isrc = 0 - regular integral, isrc=1..4 singular
!
!
!                3
!              /   \
!             /     \
!            /   4   \
!           /         \
!          1 --------- 2
!
!     ******************************************************************
  USE mCommon
  use singInt
  IMPLICIT NONE

  REAL(rk) sx,sy,sz ! source point
  REAL(rk) nx,ny,nz ! unit normal on triangle surface
  REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
  REAL(rk) area
  INTEGER isrc

! Results
  REAL(rk) Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
  REAL(rk) prGx,prGy,prGz
  REAL(rk) Txx(3),Txy(3),Txz(3),Tyy(3),Tyz(3),Tzz(3)
  REAL(rk) prTx(3),prTy(3),prTz(3)


  select case (isrc) 
    case (0)          
      !
      !       non-singular integrals
      !
      CALL Triangle_3DStokesInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,sx,sy,sz,nx,ny,nz,area, &
           Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,prGx,prGy,prGz,Txx,Txy,Txz,Tyy,Tyz,Tzz,prTx,prTy,prTz )
      !
      !       singular integrals
      !
    case default  
      Txx = 0.0_rk
      Txy = 0.0_rk
      Tyy = 0.0_rk
      Tyz = 0.0_rk
      Tzz = 0.0_rk
      prTx = 0.0_rk
      prTy = 0.0_rk
      prTz = 0.0_rk        
      call si_tri_stk(x1,y1,z1,x2,y2,z2,x3,y3,z3,area,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,isrc)
  end select

end subroutine



!
! -----------------------------------------------------------------------------------------
!
SUBROUTINE Quad_BEMIntNumAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,isrc, &
                          multi,integ,inteh)
!
!     $: Integracija Robni 4 tockovni element s 4 tockovno geometrijo
!     Analiticno singularni, numericno regularni        
!
! -----------------------------------------------------------------------------
  USE mCommon
  IMPLICIT NONE

  REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  REAL(rk) xp,yp,zp
  INTEGER isrc
  REAL(rk) integ,inteh(4)
  REAL(rk) multi ! flips normals
  
  IF (isrc.EQ.0) THEN
    ! numeric integration
    CALL Quad_BEMIntNum(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp, &
    multi,integ,inteh)
  ELSE 
    ! singular integral - analytic integration
    CALL Quad_BEMIntAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,isrc,integ,inteh) 
  END IF
END subroutine       

!
! -----------------------------------------------------------------------------------------
!
SUBROUTINE Quad_BEMIntAna(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,isrc,integ,inteh)
!
!     $: Integracija Robni 4 tockovni element s 4 tockovno geometrijo
!     Analiticno singularni
!
!
!     4 ----- 3
!     |       |
!     |   5   |
!     |       |
!     1 ----- 2    
!
! -----------------------------------------------------------------------------
  USE mCommon
  use singInt
  IMPLICIT NONE

  REAL(rk), INTENT(IN)  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  INTEGER,  INTENT(IN)  :: isrc
  REAL(rk), INTENT(OUT) :: integ,inteh(4)

  inteh(1) = 0.0_rk
  inteh(2) = 0.0_rk
  inteh(3) = 0.0_rk
  inteh(4) = 0.0_rk

  call si_qua_lap(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,integ,isrc)

END subroutine       


!
! -----------------------------------------------------------------------------------------
!
SUBROUTINE Quad_BEMIntNum(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp, &
                          multi,integ,inteh)
!
!     $: Integracija Robni 4 tockovni element s 4 tockovno geometrijo
!
! -----------------------------------------------------------------------------
  USE GaussIntegration
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER i,j,isip,ng1,ng2

  REAL(rk) xp,yp,zp
  REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  REAL(rk) xc0,yc0,zc0,xet,yet,zet,xks,yks,zks
  REAL(rk) ajac
  REAL(rk) eta1m,eta2m,eta1p,eta2p
  REAL(rk) anx,any,anz,anx1,any1,anz1
  REAL(rk) fung,funh

  REAL(rk) integ,inteh(4)
  REAL(rk) multi ! flips normals

  REAL(rk) fig4(4)
!
!      integral divison
  REAL(rk) a,b,c,d,dex,gii,gij,d1,d2,d3,d4,d13max,d24max,minedge
  INTEGER idivXi,ndivXi,idivEt,ndivEt
!
!*** SET NUMBER OF INTEGRATION POINTS
!
!     regular
  ng1=gaus%ng1(parQuadIntegRegu)
  ng2=gaus%ng2(parQuadIntegRegu)
!
!***  4 NODE CONTINUOUS BOUNDARY ELEMENT
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
!
!***  REGULAR INTEGRALS
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
!         Calculate BEM kernel
!
          CALL BemKernel(XP,YP,ZP,XC0,YC0,ZC0,ANX1,ANY1,ANZ1,FUNG,FUNH)

!         calculate H continous shape functions
          DO isip=1,4
            inteh(isip)=inteh(isip)+FUNH*AJAC*FIG4(isip)
          END DO

!         calculate G constant shape functions
          integ = integ + FUNG * AJAC
        END DO
      END DO
    END DO
  END DO

END subroutine

!
! ******************************************************************
!
SUBROUTINE Triangle_3DSingLapInt(x1,y1,z1,x2,y2,z2,x3,y3,z3, &
                                 area,integ,integ1,integ2,integ3,isrc)
!
! Analytical integration over a triangle, when source points is within triangle
!
!                3
!              /   \
!             /     \
!            /   4   \
!           1 ------- 2
!                                 
! ******************************************************************
  USE mCommon
  USE singInt
  IMPLICIT NONE
    
  INTEGER, INTENT(IN) :: isrc ! source point location
  REAL(rk),INTENT(IN) :: x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
  REAL(rk),INTENT(IN) :: area ! area of the triangle
  REAL(rk),INTENT(OUT) ::  integ,integ1,integ2,integ3 ! result
    
!  REAL(rk) a,b,c,f,jac
    
  ! since source point is in the triangle dot(field_p,source_p) = 0 -> H integrals = 0
  integ1 = 0.0_rk
  integ2 = 0.0_rk
  integ3 = 0.0_rk
  call si_tri_lap(x1,y1,z1,x2,y2,z2,x3,y3,z3,area,integ,isrc)  

!  ! translate, rotate and map triangle to reference triangle
!  call getRenChanABC(x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c)
!
!  select case (isrc)
!    case (1) ! xi = (1,0,0)
!      call getXi100(a,b,c,f)
!    case (2) ! xi = (0,1,0)
!      call getXi010(a,b,c,f)
!    case (3) ! xi = (0,0,0)
!      call getXi000(a,b,c,f)            
!    case (4) ! xi = (1/3,1/3,0)
!      call getRenChanFF2(a,b,c,f) ! center of element 
!    case default
!      print *,"Error in Triangle_3DSingLapInt"
!      stop
!  end select  
!  
!  jac = 2.0_rk * area
!  integ = 0.25_rk/pi*jac*f
   
end subroutine
