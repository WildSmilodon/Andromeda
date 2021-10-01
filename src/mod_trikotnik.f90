!
!     Triangle integration and interpolation library
!     (c) Jure Ravnik, 2017
!
      MODULE Triangle
!
!     Provides the follwing routines
!
!      -> Triangle_IntegrationInit(n)             <- provides weights and abcissas
!         Valid n
!            =1 =>  7 points (5  order polynomials accuarate)
!            =2 => 25 points (10 order polynomials accuarate)
!            =3 => 54 points (15 order polynomials accuarate)
!            =4 => 85 points (20 order polynomials accuarate)
!            =5 =>126 points (25 order polynomials accuarate)
!
!      -> Triangle_GetBarycentricCoordinates
!      -> Triangle_Area
!      -> Triangle_FunctionInt
!      -> Triangle_TestIntegration
!
!
      USE mCommon
      IMPLICIT NONE
!
! ----------------------------------------------------------------------
!
      TYPE IntegrationPointsType

        INTEGER :: n
        REAL(rk), POINTER :: l1(:),l2(:),l3(:),w(:)

      END TYPE
!
! ----------------------------------------------------------------------
!
      TYPE WandzuratXiaoType

        INTEGER :: n
        INTEGER, POINTER :: m(:)
        REAL(rk), POINTER :: x(:),y(:),w(:)

      END TYPE
!
! ----------------------------------------------------------------------
!
!
!     Variables
!
      TYPE(IntegrationPointsType) :: tipw
!
!     Subroutines
!
      CONTAINS


!
!     ******************************************************************
!
      SUBROUTINE Triangle_FunctionInt(x1,y1,z1,x2,y2,z2,x3,y3,z3,f1,f2,f3,integ)
!
!     Calculates integral over a triangle
!
!     ******************************************************************
      INTEGER i
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3
      REAL(rk) f1,f2,f3,f
      REAL(rk) area,integ

!
!     Calculate area of triangle
!
      CALL Triangle_Area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
!
!     Set result to zero
!
      integ=0.0_rk
!
!     Loop over integration points
!
      DO i=1,tipw%n
!
!       Calculate point in R^3 space, where function must be evaluated
!
!        x=tipw%l1(i)*x1+tipw%l2(i)*x2+tipw%l3(i)*x3
!        y=tipw%l1(i)*y1+tipw%l2(i)*y2+tipw%l3(i)*y3
!        z=tipw%l1(i)*z1+tipw%l2(i)*z2+tipw%l3(i)*z3
!
!       interpolate function to integration point
!
        f=tipw%l1(i)*f1+tipw%l2(i)*f2+tipw%l3(i)*f3
!
!       Sum up integral
!
        integ=integ+tipw%w(i)*f

      END DO
!
!     Multiply final result by triangle area
!
      integ=integ*area

      END SUBROUTINE

!
!     ******************************************************************
!
      SUBROUTINE Triangle_TestIntegration()
!
!     Comapres numerical integration with analytical results
!
!     ******************************************************************
      INTEGER i
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z,sx,sy,sz
      REAL(rk) area
      REAL(rk) integ1,anal1
      REAL(rk) integ2,anal2
      REAL(rk) integ3,anal3
      REAL(rk) integ4,anal4

!
!     Define test triangle
!
      x1=0.0_rk
      y1=0.0_rk
      z1=1.0_rk
      x2=2.0_rk
      y2=0.0_rk
      z2=1.0_rk
      x3=2.0_rk
      y3=2.0_rk
      z3=1.0_rk

      sx=0.0_rk
      sy=0.0_rk
      sz=1.0_rk
!
!     Calculate area of triangle
!
      CALL Triangle_Area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)

!
!     Set result to zero
!
      anal1=2097152.0_rk/121.0_rk  ! x^10*y^10
      integ1=0.0_rk

      anal2=1.0522458619541769455_rk ! sin(x*y)
      integ2=0.0_rk

      anal3=0.98364468921563619296_rk ! exp(-x*y)
      integ3=0.0_rk

      anal4=0.14027496308479503178_rk ! 1/4 pi r -> Log[1+sqrt(2)]/2pi
      integ4=0.0_rk
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
!       Sum up integral
!
        integ1=integ1+tipw%w(i)*x**10*y**10
        integ2=integ2+tipw%w(i)*Sin(x*y)
        integ3=integ3+tipw%w(i)*exp(-x*y)
        integ4=integ4+tipw%w(i)/(SQRT((x-sx)**2+(y-sy)**2+(z-sz)**2))


      END DO
!
!     Multiply final result by triangle area
!
      integ1=integ1*area
      integ2=integ2*area
      integ3=integ3*area
      integ4=integ4*area/(16.0_rk*ATAN(1.0_rk))


      print *,tipw%n,abs(integ1-anal1)/anal1,integ1,anal1
      print *,tipw%n,abs(integ2-anal2)/anal2,integ2,anal2
      print *,tipw%n,abs(integ3-anal3)/anal3,integ3,anal3
      print *,tipw%n,abs(integ4-anal4)/anal4,integ4,anal4


      END SUBROUTINE




!
!     ******************************************************************
!
      SUBROUTINE Triangle_Int(x1,y1,z1,x2,y2,z2,x3,y3,z3,integ)
!
!     Calculates integral over a triangle
!
!     ******************************************************************
      INTEGER i
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z
      REAL(rk) area,integ

!
!     Calculate area of triangle
!
      CALL Triangle_Area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
!
!     Set result to zero
!
      integ=0.0_rk
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
!       Sum up integral
!
        integ=integ+tipw%w(i)*x**10*y**10

      END DO
!
!     Multiply final result by triangle area
!
      integ=integ*area

      END SUBROUTINE


!
!     ******************************************************************
!
      SUBROUTINE Triangle_GetBarycentricCoordinates(x1,y1,z1,x2,y2,z2,x3,y3,z3,px,py,pz,l1,l2,l3)
!
!     Calculates Barycentric coordinates for a triangle in 3D
!
!     ******************************************************************
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertex
      REAL(rk) px,py,pz ! point in triangle in 3D space
      REAL(rk) l1,l2,l3 ! Barycentric coordinates of point p

      REAL(rk) v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z
      REAL(rk) d00,d01,d11,d20,d21,denom

      v0x=x2-x1
      v0y=y2-y1
      v0z=z2-z1

      v1x=x3-x1
      v1y=y3-y1
      v1z=z3-z1

      v2x=px-x1
      v2y=py-y1
      v2z=pz-z1

      d00 = v0x*v0x + v0y*v0y + v0z*v0z
      d01 = v0x*v1x + v0y*v1y + v0z*v1z
      d11 = v1x*v1x + v1y*v1y + v1z*v1z
      d20 = v2x*v0x + v2y*v0y + v2z*v0z
      d21 = v2x*v1x + v2y*v1y + v2z*v1z
      denom = d00 * d11 - d01 * d01
      l1 = (d11 * d20 - d01 * d21) / denom
      l2 = (d00 * d21 - d01 * d20) / denom
      l3 = 1.0_rk - l1 - l2

      END SUBROUTINE

!
!     ******************************************************************
!
      SUBROUTINE Triangle_IntegrationInit(n)
!
!     Calculates positions and weights for integration based on
!     Wandzurat, S., & Xiao, H. (2003). Symmetric quadrature rules on a triangle.
!     Computers & Mathematics with Applications, 45(12), 1829–1840.
!     https://doi.org/10.1016/S0898-1221(03)90004-6
!
!     Valid n
!            =1 =>  7 points (5  order polynomials accuarate)
!            =2 => 25 points (10 order polynomials accuarate)
!            =3 => 54 points (15 order polynomials accuarate)
!            =4 => 85 points (20 order polynomials accuarate)
!            =5 =>126 points (25 order polynomials accuarate)
!
!     ******************************************************************
      INTEGER n
      TYPE(WandzuratXiaoType) :: waxi

      IF (n.EQ.1) THEN
        waxi%n=3
      ELSE IF (n.EQ.2) THEN
        waxi%n=7
      ELSE IF (n.EQ.3) THEN
        waxi%n=12
      ELSE IF (n.EQ.4) THEN
        waxi%n=19
      ELSE IF (n.EQ.5) THEN
        waxi%n=26
      ELSE
        WRITE (*,*) "ERROR !"
      END IF

      ALLOCATE (waxi%x(waxi%n))
      ALLOCATE (waxi%y(waxi%n))
      ALLOCATE (waxi%w(waxi%n))
      ALLOCATE (waxi%m(waxi%n))

!
      IF (waxi%n.EQ.3) THEN ! accurate for poliynomials up to order 5

        waxi%m(1)= 1
        waxi%x(1)= 0.00_rk
        waxi%y(1)= 0.00_rk
        waxi%w(1)= 0.2250_rk


        waxi%m(2)= 3
        waxi%x(2)=-0.4104261923153453_rk ! to so v bistvu 3, zasukane za 120 stopinj
        waxi%y(2)= 0.00_rk
        waxi%w(2)= 0.1323941527885062_rk

        waxi%m(3)= 3
        waxi%x(3)= 0.6961404780296310_rk
        waxi%y(3)= 0.00_rk
        waxi%w(3)= 0.1259391805448271_rk

      ELSE IF (waxi%n.EQ.7) THEN ! accurate for poliynomials up to order 10

        waxi%m(1)= 1
        waxi%x(1)= 0.00_rk
        waxi%y(1)= 0.00_rk
        waxi%w(1)= 0.8352339980519638D-01

        waxi%m(2)= 3
        waxi%x(2)=-0.4935962988634245_rk
        waxi%y(2)= 0.00_rk
        waxi%w(2)= 0.7229850592056743D-02

        waxi%m(3)= 3
        waxi%x(3)=-0.2840373491871686_rk
        waxi%y(3)= 0.00_rk
        waxi%w(3)= 0.7449217792098051D-01

        waxi%m(4)= 3
        waxi%x(4)= 0.4457307617703263_rk
        waxi%y(4)= 0.00_rk
        waxi%w(4)= 0.7864647340310853D-01

        waxi%m(5)= 3
        waxi%x(5)= 0.9385563442849673_rk
        waxi%y(5)= 0.00_rk
        waxi%w(5)= 0.6928323087107504D-02

        waxi%m(6)= 6
        waxi%x(6)=-0.4474955151540920_rk
        waxi%y(6)=-0.5991595522781586_rk
        waxi%w(6)= 0.2951832033477940D-01

        waxi%m(7)= 6
        waxi%x(7)=-0.4436763946123360_rk
        waxi%y(7)=-0.2571781329392130_rk
        waxi%w(7)= 0.3957936719606124D-01


      ELSE IF (waxi%n.EQ.12) THEN ! accurate for poliynomials up to order 15

        waxi%m(1)= 3
        waxi%x(1)=-0.3748423891073751_rk
        waxi%y(1)= 0.00_rk
        waxi%w(1)= 0.3266181884880529D-01

        waxi%m(2)= 3
        waxi%x(2)=-0.2108313937373917_rk
        waxi%y(2)= 0.00_rk
        waxi%w(2)= 0.2741281803136436D-01

        waxi%m(3)= 3
        waxi%x(3)= 0.1204084962609239_rk
        waxi%y(3)= 0.00_rk
        waxi%w(3)= 0.2651003659870330D-01

        waxi%m(4)= 3
        waxi%x(4)= 0.5605966391716812_rk
        waxi%y(4)= 0.00_rk
        waxi%w(4)= 0.2921596213648611D-01

        waxi%m(5)= 3
        waxi%x(5)= 0.8309113970031897_rk
        waxi%y(5)= 0.00_rk
        waxi%w(5)= 0.1058460806624399D-01

        waxi%m(6)= 3
        waxi%x(6)= 0.9502746194248890_rk
        waxi%y(6)= 0.00_rk
        waxi%w(6)= 0.3614643064092035D-02

        waxi%m(7)= 6
        waxi%x(7)=-0.4851316950361628_rk
        waxi%y(7)=-0.4425551659467111_rk
        waxi%w(7)= 0.8527748101709436D-02

        waxi%m(8)= 6
        waxi%x(8)=-0.4762943440546580_rk
        waxi%y(8)=-0.1510682717598242_rk
        waxi%w(8)= 0.1391617651669193D-01

        waxi%m(9)= 6
        waxi%x(9)=-0.4922845867745440_rk
        waxi%y(9)=-0.6970224211436132_rk
        waxi%w(9)= 0.4291932940734835D-02

        waxi%m(10)= 6
        waxi%x(10)=-0.4266165113705168_rk
        waxi%y(10)=-0.5642774363966393_rk
        waxi%w(10)= 0.1623532928177489D-01

        waxi%m(11)= 6
        waxi%x(11)=-0.3968468770512212_rk
        waxi%y(11)=-0.3095105740458471_rk
        waxi%w(11)= 0.2560734092126239D-01

        waxi%m(12)= 6
        waxi%x(12)=-0.2473933728129512_rk
        waxi%y(12)=-0.2320292030461791_rk
        waxi%w(12)= 0.3308819553164567D-01


      ELSE IF (waxi%n.EQ.19) THEN ! accurate for poliynomials up to order 20

        waxi%m(1)= 1
        waxi%x(1)=0.0000000000000000_rk
        waxi%y(1)=0.0_rk
        waxi%w(1)=0.2761042699769952D-01

        waxi%m(2)= 3
        waxi%x(2)=-0.4977490260133565_rk
        waxi%y(2)=0.0_rk
        waxi%w(2)=0.1779029547326740D-02

        waxi%m(3)= 3
        waxi%x(3)=-0.3587903720915737_rk
        waxi%y(3)=0.0_rk
        waxi%w(3)=0.2011239811396117D-01

        waxi%m(4)= 3
        waxi%x(4)=-0.1932918138657104_rk
        waxi%y(4)=0.0_rk
        waxi%w(4)=0.2681784725933157D-01

        waxi%m(5)= 3
        waxi%x(5)=0.2064993924016380_rk
        waxi%y(5)=0.0_rk
        waxi%w(5)=0.2452313380150201D-01

        waxi%m(6)= 3
        waxi%x(6)=0.3669431077237697_rk
        waxi%y(6)=0.0_rk
        waxi%w(6)=0.1639457841069539D-01

        waxi%m(7)= 3
        waxi%x(7)=0.6767931784861860_rk
        waxi%y(7)=0.0_rk
        waxi%w(7)=0.1479590739864960D-01

        waxi%m(8)= 3
        waxi%x(8)=0.8827927364865920_rk
        waxi%y(8)=0.0_rk
        waxi%w(8)=0.4579282277704251D-02

        waxi%m(9)= 3
        waxi%x(9)=0.9664768608120111_rk
        waxi%y(9)=0.0_rk
        waxi%w(9)=0.1651826515576217D-02

        waxi%m(10)= 6
        waxi%x(10)=-0.4919755727189941_rk
        waxi%y(10)=-0.7513212483763635_rk
        waxi%w(10)=0.2349170908575584D-02

        waxi%m(11)= 6
        waxi%x(11)=-0.4880677744007016_rk
        waxi%y(11)=-0.5870191642967427_rk
        waxi%w(11)=0.4465925754181793D-02

        waxi%m(12)= 6
        waxi%x(12)=-0.4843664025781043_rk
        waxi%y(12)=-0.1717270984114328_rk
        waxi%w(12)=0.6099566807907972D-02

        waxi%m(13)= 6
        waxi%x(13)=-0.4835533778058150_rk
        waxi%y(13)=-0.3833898305784408_rk
        waxi%w(13)=0.6891081327188203D-02

        waxi%m(14)= 6
        waxi%x(14)=-0.4421499318718065_rk
        waxi%y(14)=-0.6563281974461070_rk
        waxi%w(14)=0.7997475072478163D-02

        waxi%m(15)= 6
        waxi%x(15)=-0.4466292382741727_rk
        waxi%y(15)=-0.6157647932662624D-01
        waxi%w(15)=0.7386134285336024D-02

        waxi%m(16)= 6
        waxi%x(16)=-0.4254937754558538_rk
        waxi%y(16)=-0.4783124082660027_rk
        waxi%w(16)=0.1279933187864826D-01

        waxi%m(17)= 6
        waxi%x(17)=-0.4122204123735024_rk
        waxi%y(17)=-0.2537089901614676_rk
        waxi%w(17)=0.1725807117569655D-01

        waxi%m(18)= 6
        waxi%x(18)=-0.3177533194934086_rk
        waxi%y(18)=-0.3996183176834929_rk
        waxi%w(18)=0.1867294590293547D-01

        waxi%m(19)= 6
        waxi%x(19)=-0.2889337325840919_rk
        waxi%y(19)=-0.1844183967233982_rk
        waxi%w(19)=0.2281822405839526D-01

      ELSE IF (waxi%n.EQ.26) THEN ! accurate for poliynomials up to order 25

        waxi%m(1)= 3
        waxi%x(1)=-0.4580802753902387_rk
        waxi%y(1)=0.0_rk
        waxi%w(1)=0.8005581880020417D-02

        waxi%m(2)= 3
        waxi%x(2)=-0.3032320980085228_rk
        waxi%y(2)=0.0_rk
        waxi%w(2)=0.1594707683239050D-01

        waxi%m(3)= 3
        waxi%x(3)=-0.1696674057318916_rk
        waxi%y(3)=0.0_rk
        waxi%w(3)=0.1310914123079553D-01

        waxi%m(4)= 3
        waxi%x(4)=0.1046702979405866_rk
        waxi%y(4)=0.0_rk
        waxi%w(4)=0.1958300096563562D-01

        waxi%m(5)= 3
        waxi%x(5)=0.2978674829878846_rk
        waxi%y(5)=0.0_rk
        waxi%w(5)=0.1647088544153727D-01

        waxi%m(6)= 3
        waxi%x(6)=0.5455949961729473_rk
        waxi%y(6)=0.0_rk
        waxi%w(6)=0.8547279074092100D-02

        waxi%m(7)= 3
        waxi%x(7)=0.6617983193620190_rk
        waxi%y(7)=0.0_rk
        waxi%w(7)=0.8161885857226492D-02

        waxi%m(8)= 3
        waxi%x(8)=0.7668529237254211_rk
        waxi%y(8)=0.0_rk
        waxi%w(8)=0.6121146539983779D-02

        waxi%m(9)= 3
        waxi%x(9)=0.8953207191571090_rk
        waxi%y(9)=0.0_rk
        waxi%w(9)=0.2908498264936665D-02

        waxi%m(10)= 3
        waxi%x(10)=0.9782254461372029_rk
        waxi%y(10)=0.0_rk
        waxi%w(10)=0.6922752456619963D-03

        waxi%m(11)= 6
        waxi%x(11)=-0.4980614709433367_rk
        waxi%y(11)=-0.4713592181681879_rk
        waxi%w(11)= 0.1248289199277397D-02

        waxi%m(12)= 6
        waxi%x(12)=-0.4919004480918257_rk
        waxi%y(12)=-0.1078887424748246_rk
        waxi%w(12)= 0.3404752908803022D-02

        waxi%m(13)= 6
        waxi%x(13)=-0.4904239954490375_rk
        waxi%y(13)=-0.3057041948876942_rk
        waxi%w(13)= 0.3359654326064051D-02

        waxi%m(14)= 6
        waxi%x(14)=-0.4924576827470104_rk
        waxi%y(14)=-0.7027546250883238_rk
        waxi%w(14)= 0.1716156539496754D-02

        waxi%m(15)= 6
        waxi%x(15)=-0.4897598620673272_rk
        waxi%y(15)=-0.7942765584469995_rk
        waxi%w(15)= 0.1480856316715606D-02

        waxi%m(16)= 6
        waxi%x(16)=-0.4849757005401057_rk
        waxi%y(16)=-0.5846826436376921_rk
        waxi%w(16)= 0.3511312610728685D-02

        waxi%m(17)= 6
        waxi%x(17)=-0.4613632802399150_rk
        waxi%y(17)=-0.4282174042835178_rk
        waxi%w(17)= 0.7393550149706484D-02

        waxi%m(18)= 6
        waxi%x(18)=-0.4546581528201263_rk
        waxi%y(18)=-0.2129434060653430_rk
        waxi%w(18)= 0.7983087477376558D-02

        waxi%m(19)= 6
        waxi%x(19)=-0.4542425148392569_rk
        waxi%y(19)=-0.6948910659636692_rk
        waxi%w(19)= 0.4355962613158041D-02

        waxi%m(20)= 6
        waxi%x(20)=-0.4310651789561460_rk
        waxi%y(20)=-0.5691146659505208_rk
        waxi%w(20)= 0.7365056701417832D-02


        waxi%m(21)= 6
        waxi%x(21)=-0.3988357991895837_rk
        waxi%y(21)=-0.3161666335733065_rk
        waxi%w(21)= 0.1096357284641955D-01

        waxi%m(22)= 6
        waxi%x(22)=-0.3949323628761341_rk
        waxi%y(22)=-0.1005941839340892_rk
        waxi%w(22)= 0.1174996174354112D-01

        waxi%m(23)= 6
        waxi%x(23)=-0.3741327130398251_rk
        waxi%y(23)=-0.4571406037889341_rk
        waxi%w(23)= 0.1001560071379857D-01

        waxi%m(24)= 6
        waxi%x(24)=-0.3194366964842710_rk
        waxi%y(24)=-0.2003599744104858_rk
        waxi%w(24)= 0.1330964078762868D-01

        waxi%m(25)= 6
        waxi%x(25)=-0.2778996512639500_rk
        waxi%y(25)=-0.3406754571040736_rk
        waxi%w(25)= 0.1415444650522614D-01

        waxi%m(26)= 6
        waxi%x(26)=-0.2123422011990124_rk
        waxi%y(26)=-0.1359589640107579_rk
        waxi%w(26)= 0.1488137956116801D-01



      END IF
!
!     Calculate Barycentric coordinates
!
      CALL WandzuratXiaoBary(waxi)

      END SUBROUTINE

!
!     ******************************************************************
!
      SUBROUTINE WandzuratXiaoBary(waxi)
!
!     Calculates positions and weights for integration based on
!     Wandzurat, S., & Xiao, H. (2003). Symmetric quadrature rules on a triangle.
!     Computers & Mathematics with Applications, 45(12), 1829–1840.
!     https://doi.org/10.1016/S0898-1221(03)90004-6
!
!     ******************************************************************
      TYPE(WandzuratXiaoType) :: waxi
      INTEGER i,j

      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertex
      REAL(rk) c,s,pi,c4,s4

      REAL(rk), ALLOCATABLE :: x(:),y(:)

!     unit triangle vertices
      x1=-0.5_rk
      y1=+0.5_rk*sqrt(3.0_rk)
      z1= 0.0_rk
      x2=-0.5_rk
      y2=-0.5_rk*sqrt(3.0_rk)
      z2= 0.0_rk
      x3= 1.0_rk
      y3= 0.0_rk
      z3= 0.0_rk

      pi=4.0_rk*DATAN(1.00_rk)
      s=sin(2.0_rk/3.0_rk*pi) ! sin(120)
      c=cos(2.0_rk/3.0_rk*pi) ! cos(120)
      s4=sin(4.0_rk/3.0_rk*pi) ! sin(240)
      c4=cos(4.0_rk/3.0_rk*pi) ! cos(240)


      tipw%n=0
      DO i=1,waxi%n
        tipw%n=tipw%n+waxi%m(i)
      END DO

      ALLOCATE (x(tipw%n))
      ALLOCATE (y(tipw%n))
      ALLOCATE (tipw%l1(tipw%n))
      ALLOCATE (tipw%l2(tipw%n))
      ALLOCATE (tipw%l3(tipw%n))
      ALLOCATE (tipw%w(tipw%n))
!
!     Transform the waxi points with 3 or 6 rotations (see Table 1 in paper)
!
      i=0
      DO j=1,waxi%n
        IF (waxi%m(j).EQ.1) THEN
          i=i+1
          x(i)=waxi%x(j)
          y(i)=waxi%y(j)
          tipw%w(i)=waxi%w(j)
        ELSE IF (waxi%m(j).EQ.3) THEN
          i=i+1
          x(i)=waxi%x(j)
          y(i)=waxi%y(j)
          tipw%w(i)=waxi%w(j)
          x(i+1)= x(i)*c-y(i)*s
          y(i+1)=-x(i)*s-y(i)*c
          tipw%w(i+1)= tipw%w(i)
          x(i+2)= x(i+1)*c+y(i+1)*s
          y(i+2)=-x(i+1)*s+y(i+1)*c
          tipw%w(i+2)= tipw%w(i+1)
          i=i+2
        ELSE IF (waxi%m(j).EQ.6) THEN
          i=i+1
          x(i)=waxi%x(j)
          y(i)=waxi%y(j)
          tipw%w(i)=waxi%w(j)
          x(i+1)= x(i)*c-y(i)*s
          y(i+1)=-x(i)*s-y(i)*c
          tipw%w(i+1)= tipw%w(i)
          x(i+2)= x(i+1)*c+y(i+1)*s
          y(i+2)=-x(i+1)*s+y(i+1)*c
          tipw%w(i+2)= tipw%w(i+1)
          i=i+2
          i=i+1
          x(i)= waxi%x(j)
          y(i)=-waxi%y(j)
          tipw%w(i)= waxi%w(j)
          x(i+1)= x(i)*c4-y(i)*s4
          y(i+1)=-x(i)*s4-y(i)*c4
          tipw%w(i+1)= tipw%w(i)
          x(i+2)= x(i+1)*c4+y(i+1)*s4
          y(i+2)=-x(i+1)*s4+y(i+1)*c4
          tipw%w(i+2)= tipw%w(i+1)
          i=i+2
        END IF
      END DO

      DO i=1,tipw%n
        CALL Triangle_GetBarycentricCoordinates(x1,y1,z1,x2,y2,z2,x3,y3,z3,x(i),y(i),0.0_rk,tipw%l1(i),tipw%l2(i),tipw%l3(i))
      END DO

      DEALLOCATE (x,y)

      END SUBROUTINE

!
!     ******************************************************************
!
      SUBROUTINE Triangle_Area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
!
!     Calculates area of triangle
!
!     ******************************************************************
      REAL(rk) x1,y1,z1,x2,y2,z2,x3,y3,z3,area,a,b,c,p

        a=SQRT( (x1-x2)**2.0_rk + (y1-y2)**2.0_rk +(z1-z2)**2.0_rk )
        b=SQRT( (x3-x2)**2.0_rk + (y3-y2)**2.0_rk +(z3-z2)**2.0_rk )
        c=SQRT( (x1-x3)**2.0_rk + (y1-y3)**2.0_rk +(z1-z3)**2.0_rk )

!       Heron's formula

        p=0.5_rk*(a+b+c)
        area=SQRT( p*(p-a)*(p-b)*(p-c) )

      END SUBROUTINE





      END MODULE

