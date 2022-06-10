!
! Singular integration for Laplace and Stokes kernels
! for BEM using linear interpolation of function and constant interpolation
! of flux using triangular or quadrilateral boundary elements
!
! Module provides 4 subroutines calculating singular integrals for boundary element
! Needed  : - coordinates of element vertices
!           - source point location (vertex or center of element)
!
!
! Module is free to use, provided that citation is given to
!
! J. Ravnik:    , 2022, EABE
!
!
module singInt
   implicit none

   private

   public :: si_tri_lap  ! Laplace kernel, triangular boundary element
   public :: si_tri_stk  ! Stokes kernel, triangular boundary element
   public :: si_qua_lap  ! Laplace kernel, quadrilateral boundary element
   public :: si_qua_stk  ! Stokes kernel, quadrilateral boundary element

   integer,  parameter :: rk = selected_real_kind(8)
   real(rk), parameter :: pi = 3.141592653589793238462643383_rk

   contains

!
! ------------------------------------------------------------------------------
!
   subroutine si_tri_stk(x1,y1,z1,x2,y2,z2,x3,y3,z3,area,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,isrc)
!
!
!  Analytical integration over a triangle, when source points is within triangle
!  Stokes kernel
!
!                3
!              /   \
!             /     \
!            /   4   \
!           1 ------- 2
!
! ------------------------------------------------------------------------------
!
      real(rk), intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
      real(rk), intent(in) :: area   ! area of the triangle
      integer,  intent(in) ::  isrc  ! source point location(=1,2,3 or 4)

      real(rk), intent(out) :: Gxx,Gxy,Gxz,Gyy,Gyz,Gzz

      ! Temporary
      real(rk) a,b,c,dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy
      real(rk) ff,ggxx,ggxy,ggyy

      real(rk), allocatable :: Mrot(:,:) ! rotation matrix (3D)
      real(rk), allocatable :: MrotT(:,:) ! transpose of rotation matrix (3D)
      real(rk), allocatable :: mat(:,:) ! my matrix


      allocate ( Mrot(3,3),MrotT(3,3),mat(3,3) )  

!     translate, rotate and map triangle to reference triangle    
      call getABCDEF(x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy,Mrot)
!     get transpose of rotation matrix
      MrotT = TRANSPOSE(Mrot)

      select case (isrc) 

         case (1) ! xi = (1,0,0)
            call getXi100(a,b,c,ff)
            call getXi100G(a,b,c,dxx,exx,fxx,ggxx)
            call getXi100G(a,b,c,dxy,exy,fxy,ggxy)
            call getXi100G(a,b,c,dyy,eyy,fyy,ggyy)      
         
         case (2) ! xi = (0,1,0)
            call getXi010(a,b,c,ff)
            call getXi010G(a,b,c,dxx,exx,fxx,ggxx)
            call getXi010G(a,b,c,dxy,exy,fxy,ggxy)
            call getXi010G(a,b,c,dyy,eyy,fyy,ggyy)
         
         case (3) ! xi = (0,0,0)
            call getXi000(a,b,c,ff)            
            call getXi000G(a,b,c,dxx,exx,fxx,ggxx)
            call getXi000G(a,b,c,dxy,exy,fxy,ggxy)
            call getXi000G(a,b,c,dyy,eyy,fyy,ggyy)
         
         case (4) ! xi = (1/3,1/3,0)
            call getXi13130(a,b,c,ff) ! center of element   
            call getXi13130G(a,b,c,dxx,exx,fxx,ggxx)
            call getXi13130G(a,b,c,dxy,exy,fxy,ggxy)
            call getXi13130G(a,b,c,dyy,eyy,fyy,ggyy)
         
         case default          
            print *,"Error in si_tri_stk"
            stop
      end select


      ! set up my matrix
      mat = 0.0_rk
      mat(1,1) = 2.0_rk * area * ( ff + ggxx )
      mat(2,2) = 2.0_rk * area * ( ff + ggyy )
      mat(3,3) = 2.0_rk * area * ( ff )

      mat(1,2) = 2.0_rk * area * ( ggxy )
      mat(2,1) = mat(1,2)

      ! rotate to original frame of refence
      mat = MATMUL(MATMUL(MrotT,mat),Mrot)

      Gxx = mat(1,1)
      Gxy = mat(1,2)
      Gxz = mat(1,3)
      Gyy = mat(2,2)
      Gyz = mat(2,3)
      Gzz = mat(3,3)

      deallocate (Mrot,MrotT,mat)

   end subroutine
!
! ------------------------------------------------------------------------------
!
   subroutine si_tri_lap(x1,y1,z1,x2,y2,z2,x3,y3,z3,area,G,isrc)
!
!  Analytical integration over a triangle, when source points is within triangle
!  Laplace kernel
!
!                3
!              /   \
!             /     \
!            /   4   \
!           1 ------- 2
!
! ------------------------------------------------------------------------------
!
      integer, intent(in) :: isrc ! source point location(=1,2,3 or 4)
      real(rk),intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
      real(rk),intent(in) :: area ! area of the triangle
      real(rk),intent(out) ::  G ! result

      real(rk) a,b,c,f,jac

      ! translate, rotate and map triangle to reference triangle
      call getrcABC(x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c)

      select case (isrc)
         case (1) ! xi = (1,0,0)
            call getXi100(a,b,c,f)
         case (2) ! xi = (0,1,0)
            call getXi010(a,b,c,f)
         case (3) ! xi = (0,0,0)
            call getXi000(a,b,c,f)            
         case (4) ! xi = (1/3,1/3,0)
            call getXi13130(a,b,c,f) ! center of element 
         case default
            print *,"Error in si_tri_lap"
         stop
      end select  

      jac = 2.0_rk * area
      G = 0.25_rk/pi*jac*f

   end subroutine

!
! ------------------------------------------------------------------------------
!
   subroutine si_qua_lap(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,G,isrc)
!
!     Singular quadrilateral integration, isrc = 1,2,3,4 or 5
!     Laplace kernel
!
!
!     4 ----- 3
!     |       |
!     |   5   |
!     |       |
!     1 ----- 2    
!
! -----------------------------------------------------------------------------
      
      real(rk), intent(in)  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      integer,  intent(in)  :: isrc
      real(rk), intent(out) :: G
      
      real(rk) int1,int2,int3,int4,x5,y5,z5
      real(rk) area
      
      select case (isrc) 
         case (1)          
            
            call tri_area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
            call si_tri_lap(x1,y1,z1,x2,y2,z2,x3,y3,z3,area,int1,1)
      
            call tri_area(x1,y1,z1,x3,y3,z3,x4,y4,z4,area)
            call si_tri_lap(x1,y1,z1,x3,y3,z3,x4,y4,z4,area,int2,1)
      
            G = int1 + int2
      
         case (2)          
      
            call tri_area(x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
            call si_tri_lap(x2,y2,z2,x3,y3,z3,x4,y4,z4,area,int1,1)
      
            call tri_area(x2,y2,z2,x4,y4,z4,x1,y1,z1,area)
            call si_tri_lap(x2,y2,z2,x4,y4,z4,x1,y1,z1,area,int2,1)
      
            G = int1 + int2
            
         case (3)          
      
            call tri_area(x3,y3,z3,x4,y4,z4,x1,y1,z1,area)
            call si_tri_lap(x3,y3,z3,x4,y4,z4,x1,y1,z1,area,int1,1)
      
            call tri_area(x3,y3,z3,x1,y1,z1,x2,y2,z2,area)
            call si_tri_lap(x3,y3,z3,x1,y1,z1,x2,y2,z2,area,int2,1)
      
            G = int1 + int2
      
      
         case (4)          
      
            call tri_area(x4,y4,z4,x1,y1,z1,x2,y2,z2,area)
            call si_tri_lap(x4,y4,z4,x1,y1,z1,x2,y2,z2,area,int1,1)
      
            call tri_area(x4,y4,z4,x2,y2,z2,x3,y3,z3,area)
            call si_tri_lap(x4,y4,z4,x2,y2,z2,x3,y3,z3,area,int2,1)
      
            G = int1 + int2
      
         case (5)      
            
            x5 = 0.25_rk * (x1+x2+x3+x4)
            y5 = 0.25_rk * (y1+y2+y3+y4)
            z5 = 0.25_rk * (z1+z2+z3+z4)
      
            call tri_area(x5,y5,z5,x1,y1,z1,x2,y2,z2,area)
            call si_tri_lap(x5,y5,z5,x1,y1,z1,x2,y2,z2,area,int1,1)
      
            call tri_area(x5,y5,z5,x2,y2,z2,x3,y3,z3,area)
            call si_tri_lap(x5,y5,z5,x2,y2,z2,x3,y3,z3,area,int2,1)
      
            call tri_area(x5,y5,z5,x3,y3,z3,x4,y4,z4,area)
            call si_tri_lap(x5,y5,z5,x3,y3,z3,x4,y4,z4,area,int3,1)
      
            call tri_area(x5,y5,z5,x4,y4,z4,x1,y1,z1,area)
            call si_tri_lap(x5,y5,z5,x4,y4,z4,x1,y1,z1,area,int4,1)
      
            G = int1 + int2 + int3 + int4
      
         case default          
            print *,"ERROR in si_qua_lap"
            stop
      end select
      
   end subroutine       
      

!
! ------------------------------------------------------------------------------
!
   subroutine si_qua_stk(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,isrc)
!
!     Singular quadrilateral integration, isrc = 1,2,3,4 or 5
!     Stokes kernel
!
!
!     4 ----- 3
!     |       |
!     |   5   |
!     |       |
!     1 ----- 2    
!
!
! ------------------------------------------------------------------------------
      real(rk), intent(in)  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      integer,  intent(in)  :: isrc

      real(rk), intent(out) ::  Gxx,Gxy,Gxz,Gyy,Gyz,Gzz
  
      real(rk) area,x5,y5,z5
      real(rk) Gxx1,Gxy1,Gxz1,Gyy1,Gyz1,Gzz1
      real(rk) Gxx2,Gxy2,Gxz2,Gyy2,Gyz2,Gzz2

      real(rk) Gxx3,Gxy3,Gxz3,Gyy3,Gyz3,Gzz3
      real(rk) Gxx4,Gxy4,Gxz4,Gyy4,Gyz4,Gzz4

      select case (isrc) 
         case (1)          
    
            call tri_area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
            call si_tri_stk(x1,y1,z1,x2,y2,z2,x3,y3,z3,area,Gxx1,Gxy1,Gxz1,Gyy1,Gyz1,Gzz1,1 )

            call tri_area(x1,y1,z1,x3,y3,z3,x4,y4,z4,area)
            call si_tri_stk(x1,y1,z1,x3,y3,z3,x4,y4,z4,area,Gxx2,Gxy2,Gxz2,Gyy2,Gyz2,Gzz2,1 )

            Gxx = Gxx1 + Gxx2
            Gxy = Gxy1 + Gxy2
            Gxz = Gxz1 + Gxz2

            Gyy = Gyy1 + Gyy2
            Gyz = Gyz1 + Gyz2
            Gzz = Gzz1 + Gzz2

         case (2)          

            call tri_area(x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
            call si_tri_stk(x2,y2,z2,x3,y3,z3,x4,y4,z4,area,Gxx1,Gxy1,Gxz1,Gyy1,Gyz1,Gzz1,1)
                  
            call tri_area(x2,y2,z2,x4,y4,z4,x1,y1,z1,area)
            call si_tri_stk(x2,y2,z2,x4,y4,z4,x1,y1,z1,area,Gxx2,Gxy2,Gxz2,Gyy2,Gyz2,Gzz2,1)
                  
            Gxx = Gxx1 + Gxx2
            Gxy = Gxy1 + Gxy2
            Gxz = Gxz1 + Gxz2
                  
            Gyy = Gyy1 + Gyy2
            Gyz = Gyz1 + Gyz2
            Gzz = Gzz1 + Gzz2

         case (3)          

            call tri_area(x3,y3,z3,x4,y4,z4,x1,y1,z1,area)
            call si_tri_stk(x3,y3,z3,x4,y4,z4,x1,y1,z1,area,Gxx1,Gxy1,Gxz1,Gyy1,Gyz1,Gzz1,1)

            call tri_area(x3,y3,z3,x1,y1,z1,x2,y2,z2,area)
            call si_tri_stk(x3,y3,z3,x1,y1,z1,x2,y2,z2,area,Gxx2,Gxy2,Gxz2,Gyy2,Gyz2,Gzz2,1)

            Gxx = Gxx1 + Gxx2
            Gxy = Gxy1 + Gxy2
            Gxz = Gxz1 + Gxz2

            Gyy = Gyy1 + Gyy2
            Gyz = Gyz1 + Gyz2
            Gzz = Gzz1 + Gzz2

         case (4)          

            call tri_area(x4,y4,z4,x1,y1,z1,x2,y2,z2,area)
            call si_tri_stk(x4,y4,z4,x1,y1,z1,x2,y2,z2,area,Gxx1,Gxy1,Gxz1,Gyy1,Gyz1,Gzz1,1)

            call tri_area(x4,y4,z4,x2,y2,z2,x3,y3,z3,area)
            call si_tri_stk(x4,y4,z4,x2,y2,z2,x3,y3,z3,area,Gxx2,Gxy2,Gxz2,Gyy2,Gyz2,Gzz2,1)

            Gxx = Gxx1 + Gxx2
            Gxy = Gxy1 + Gxy2
            Gxz = Gxz1 + Gxz2

            Gyy = Gyy1 + Gyy2
            Gyz = Gyz1 + Gyz2
            Gzz = Gzz1 + Gzz2
   
         case (5)      
    
            x5 = 0.25_rk * (x1+x2+x3+x4)
            y5 = 0.25_rk * (y1+y2+y3+y4)
            z5 = 0.25_rk * (z1+z2+z3+z4)

            call tri_area(x5,y5,z5,x1,y1,z1,x2,y2,z2,area)
            call si_tri_stk(x5,y5,z5,x1,y1,z1,x2,y2,z2,area,Gxx1,Gxy1,Gxz1,Gyy1,Gyz1,Gzz1,1)

            call tri_area(x5,y5,z5,x2,y2,z2,x3,y3,z3,area)
            call si_tri_stk(x5,y5,z5,x2,y2,z2,x3,y3,z3,area,Gxx2,Gxy2,Gxz2,Gyy2,Gyz2,Gzz2,1)

            call tri_area(x5,y5,z5,x3,y3,z3,x4,y4,z4,area)
            call si_tri_stk(x5,y5,z5,x3,y3,z3,x4,y4,z4,area,Gxx3,Gxy3,Gxz3,Gyy3,Gyz3,Gzz3,1)

            call tri_area(x5,y5,z5,x4,y4,z4,x1,y1,z1,area)
            call si_tri_stk(x5,y5,z5,x4,y4,z4,x1,y1,z1,area,Gxx4,Gxy4,Gxz4,Gyy4,Gyz4,Gzz4,1)

            Gxx = Gxx1 + Gxx2 + Gxx3 + Gxx4
            Gxy = Gxy1 + Gxy2 + Gxy3 + Gxy4
            Gxz = Gxz1 + Gxz2 + Gxz3 + Gxz4

            Gyy = Gyy1 + Gyy2 + Gyy3 + Gyy4
            Gyz = Gyz1 + Gyz2 + Gyz3 + Gyz4
            Gzz = Gzz1 + Gzz2 + Gzz3 + Gzz4

         case default          
            print *,"ERROR in si_qua_stk"
            stop
      end select
   end subroutine



!
! ------------------------------------------------------------------------------
!
   subroutine tri_area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
!
!     Calculates area of triangle
!
! ------------------------------------------------------------------------------
!
      real(rk), intent(in)  :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(rk), intent(out) :: area
      real(rk) a,b,c,p

      a=SQRT( (x1-x2)**2.0_rk + (y1-y2)**2.0_rk +(z1-z2)**2.0_rk )
      b=SQRT( (x3-x2)**2.0_rk + (y3-y2)**2.0_rk +(z3-z2)**2.0_rk )
      c=SQRT( (x1-x3)**2.0_rk + (y1-y3)**2.0_rk +(z1-z3)**2.0_rk )

!     Heron's formula

      p=0.5_rk*(a+b+c)
      area=SQRT( p*(p-a)*(p-b)*(p-c) )
   
   end subroutine
         
!
! ------------------------------------------------------------------------------
!
   subroutine getXi100(a,b,c,f)
      
      real(rk),intent(in) :: a,b,c ! a,b,c - Eq. (23)
      real(rk),intent(out) ::  f 
      real(rk) d,e
    
      d = 2.0_rk * sqrt(b)*sqrt(a+b-c)+2.0_rk*b-c
      e = 2.0_rk * sqrt(a)*sqrt(b)-c
    
      f = log(d/e)/sqrt(b)
    
    end subroutine
    
!
! ------------------------------------------------------------------------------
!
   subroutine getXi010(a,b,c,f)

      
      real(rk),intent(in) :: a,b,c ! a,b,c - Eq. (23)
      real(rk),intent(out) ::  f 
      real(rk) d,e
    
      d = 2.0_rk * sqrt(a) * sqrt(a+b-c) + 2.0_rk * a - c
      e = 2.0_rk * sqrt(a) * sqrt(b) - c
    
      f = log(d/e)/sqrt(a)
    
    end subroutine
    
    
!
! ------------------------------------------------------------------------------
!
    subroutine getXi000(a,b,c,f)
      
      real(rk),intent(in) :: a,b,c ! a,b,c - Eq. (23)
      real(rk),intent(out) ::  f 
      real(rk) d,e
    
      d = 2.0_rk * sqrt(a) * sqrt(a+b-c) + 2.0_rk * a - c
      e = 2.0_rk * sqrt(b) * sqrt(a+b-c) - 2.0_rk * b + c
    
      f = log(d/e)/sqrt(a+b-c)
    
    end subroutine
    
!
! ------------------------------------------------------------------------------
!
   subroutine getXi13130G(a,b,c,d,e,f,G)
!
!     Implements Eq. (34) in the paper:
!     Calculates integral over a triangle, using analitical method of
!     Ren, Q., & Chan, C. L. (2015). Analytical evaluation of the BEM singular integrals 
!     for 3D Laplace and Stokes flow equations using coordinate transformation. 
!     Engineering Analysis with Boundary Elements, 53, 1–8. 
!     https://doi.org/10.1016/j.enganabound.2014.11.018
!
!
! ------------------------------------------------------------------------------
!
      
      real(rk),intent(in) :: a,b,c,d,e,f 
      real(rk),intent(out) ::  G
     
      real(rk) st,im
    
      G = 0.0_rk
    
      st = 2.0_rk*sqrt(b)*(b*c*d+a*c*e-2.0_rk*a*b*f)*(sqrt(a+b+c)-sqrt(4.0_rk*a+b-2.0_rk*c))
      im = a*(4.0_rk*a*b-c*c)
    
      G = G + st/im
    
      st = 2.0_rk*sqrt(b)*(  &
        b*(2.0_rk*b-c)*d-a*(2.0_rk*b*e+c*e-2.0_rk*b*f)+c*(c*e-b*f) &
        )*( sqrt(a+4.0_rk*b-2.0_rk*c) - sqrt(4.0_rk*a+b-2.0_rk*c) ) 
      im = (a+b-c)*(4.0_rk*a*b-c*c)
    
      G = G + st/im
    
      st = 4.0_rk*a+2.0_rk*b+2.0_rk*sqrt(4.0_rk*a+b-2.0_rk*c)*sqrt(a+b-c)-3.0_rk*c
      im =-2.0_rk*a-4.0_rk*b+2.0_rk*sqrt(a+4.0_rk*b-2.0_rk*c)*sqrt(a+b-c)+3.0_rk*c
    
      G = G + e*sqrt(b)*log(st/im)/sqrt(a+b-c)
      G = G + sqrt(b)*(b*d-a*e+c*e-b*f)*log(st/im)/(a+b-c)**1.5_rk
    
      st = 4.0_rk*a+2*sqrt(a)*sqrt(4.0_rk*a+b-2.0_rk*c)-c
      im =-2.0_rk*a-c+2*sqrt(a)*sqrt(a+b+c) 
    
      G = G + e*sqrt(b)*log(st/im)/sqrt(a)
      G = G + sqrt(b)*(b*d-a*e)*log(st/im)/a**1.5_rk
    
      st = (4.0_rk*b+2.0_rk*sqrt(b)*sqrt(a+4.0_rk*b-2.0_rk*c)-c)* &
           (2.0_rk*b+2.0_rk*sqrt(b)*sqrt(a+b+c)              +c)
      im = 4.0_rk* ( b+sqrt(b)*sqrt(4.0_rk*a+b-2.0_rk*c)-c )*     &
                   (-b+sqrt(b)*sqrt(4.0_rk*a+b-2.0_rk*c)+c ) 
    
      G = G + e*log(st/im)
      G = G + e*log(4.0_rk)
      G = G / (3.0_rk*b**1.5_rk)
    
   end subroutine
!
! ------------------------------------------------------------------------------
!
   subroutine getXi13130Gmathematica(a,b,c,d,e,f,G)
!
! Direct FortanForm output from Mathematica, resuls the same as getXi13130G
!            
      
      real(rk),intent(in) :: a,b,c,d,e,f 
      real(rk),intent(out) ::  G
      
    
      G =  &
      (8.0_rk*b*d*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (2.0_rk*c*d*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (8.0_rk*a*e*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (2.0_rk*a*c*e*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*b*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) +   &
        (4.0_rk*c**2.0_rk*e*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*b*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) +   &
        (4.0_rk*a*f*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (4.0_rk*c*f*(Log(2.0_rk*a + 4.0_rk*b +   &
               2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*Sqrt(a + b - c) -   &
               3.0_rk*c) - Log(-4.0_rk*a - 2.0_rk*b +   &
               2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*Sqrt(a + b - c) +   &
               3.0_rk*c)))/  &
         (3.0_rk*Sqrt(a + b - c)*(4.0_rk*a*b - c**2.0_rk)) +   &
        (e*(-3 + Log(4.0_rk) +   &
             (Sqrt(b)*Log(4.0_rk*a + 2.0_rk*b +   &
                  2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                   Sqrt(a + b - c) - 3.0_rk*c))/  &
              Sqrt(a + b - c) +   &
             Log(4.0_rk*b + 2.0_rk*Sqrt(b)*Sqrt(a + 4.0_rk*b - 2.0_rk*c) -   &
               c) + 2.0_rk*Log(-b +   &
                Sqrt(b)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) -   &
             (Sqrt(b)*Log(-2.0_rk*a - 4.0_rk*b +   &
                  2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                   Sqrt(a + b - c) + 3.0_rk*c))/  &
              Sqrt(a + b - c)))/(3.0_rk*b**1.5_rk) -   &
        (2.0_rk*b*d*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) +   &
        (c*d*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) +   &
        (2.0_rk*a*e*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) +   &
        (a*c*e*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*b*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) -   &
        (c**2.0_rk*e*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*b*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) -   &
        (2.0_rk*a*f*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) +   &
        (c*f*(2.0_rk*(Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
                Sqrt(a + 4.0_rk*b - 2.0_rk*c))*Sqrt(a + b - c) +   &
             (2.0_rk*a + 4.0_rk*b - 3.0_rk*c)*  &
              Log(4.0_rk*a + 2.0_rk*b +   &
                2.0_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) - 3.0_rk*c) +   &
             (-2.0_rk*a - 4.0_rk*b + 3.0_rk*c)*  &
              Log(-2.0_rk*a - 4.0_rk*b +   &
                2.0_rk*Sqrt(a + 4.0_rk*b - 2.0_rk*c)*  &
                 Sqrt(a + b - c) + 3.0_rk*c)))/  &
         (3.0_rk*(a + b - c)**1.5_rk*(4.0_rk*a*b - c**2.0_rk)) -   &
        (c*d*(-2.0_rk*Sqrt(a)*  &
              (a*Sqrt(4.0_rk*a + b - 2.0_rk*c) +   &
                b*Sqrt(4.0_rk*a + b - 2.0_rk*c) +   &
                Sqrt(4.0_rk*a + b - 2.0_rk*c)*c -   &
                4.0_rk*a*Sqrt(a + b + c) -   &
                b*Sqrt(a + b + c) + 2.0_rk*c*Sqrt(a + b + c)  &
                ) + Sqrt(4.0_rk*a + b - 2.0_rk*c)*(2.0_rk*a + c)*  &
              Sqrt(a + b + c)*  &
              Log(4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) - c) -   &
             Sqrt(4.0_rk*a + b - 2.0_rk*c)*(2.0_rk*a + c)*  &
              Sqrt(a + b + c)*  &
              Log(-2.0_rk*a - c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))  &
             ))/  &
         (3.0_rk*a**1.5_rk*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
           Sqrt(a + b + c)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (c*e*(-2.0_rk*Sqrt(a)*  &
              (a*Sqrt(4.0_rk*a + b - 2.0_rk*c) +   &
                b*Sqrt(4.0_rk*a + b - 2.0_rk*c) +   &
                Sqrt(4.0_rk*a + b - 2.0_rk*c)*c -   &
                4.0_rk*a*Sqrt(a + b + c) -   &
                b*Sqrt(a + b + c) + 2.0_rk*c*Sqrt(a + b + c)  &
                ) + Sqrt(4.0_rk*a + b - 2.0_rk*c)*(2.0_rk*a + c)*  &
              Sqrt(a + b + c)*  &
              Log(4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) - c) -   &
             Sqrt(4.0_rk*a + b - 2.0_rk*c)*(2.0_rk*a + c)*  &
              Sqrt(a + b + c)*  &
              Log(-2.0_rk*a - c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))  &
             ))/  &
         (3.0_rk*Sqrt(a)*b*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
           Sqrt(a + b + c)*(4.0_rk*a*b - c**2.0_rk)) +   &
        (2.0_rk*f*(-2.0_rk*Sqrt(a)*  &
              (a*Sqrt(4.0_rk*a + b - 2.0_rk*c) +   &
                b*Sqrt(4.0_rk*a + b - 2.0_rk*c) +   &
                Sqrt(4.0_rk*a + b - 2.0_rk*c)*c -   &
                4.0_rk*a*Sqrt(a + b + c) -   &
                b*Sqrt(a + b + c) + 2.0_rk*c*Sqrt(a + b + c)  &
                ) + Sqrt(4.0_rk*a + b - 2.0_rk*c)*(2.0_rk*a + c)*  &
              Sqrt(a + b + c)*  &
              Log(4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) - c) -   &
             Sqrt(4.0_rk*a + b - 2.0_rk*c)*(2.0_rk*a + c)*  &
              Sqrt(a + b + c)*  &
              Log(-2.0_rk*a - c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))  &
             ))/  &
         (3.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c)*  &
           Sqrt(a + b + c)*(4.0_rk*a*b - c**2.0_rk)) +   &
        (4.0_rk*b*d*(-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*Sqrt(a)*(4.0_rk*a*b - c**2.0_rk)) +   &
        (2.0_rk*c*d*(-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*Sqrt(a)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (4.0_rk*Sqrt(a)*e*(-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*(4.0_rk*a*b - c**2.0_rk)) +   &
        (2.0_rk*Sqrt(a)*c*e*  &
           (-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*b*(4.0_rk*a*b - c**2.0_rk)) +   &
        (2.0_rk*c**2.0_rk*e*(-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*Sqrt(a)*b*(4.0_rk*a*b - c**2.0_rk)) -   &
        (4.0_rk*Sqrt(a)*f*(-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*(4.0_rk*a*b - c**2.0_rk)) -   &
        (2.0_rk*c*f*(-Log(-4.0_rk*a +   &
                2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             Log(2.0_rk*a + c + 2.0_rk*Sqrt(a)*Sqrt(a + b + c))))  &
          /(3.0_rk*Sqrt(a)*(4.0_rk*a*b - c**2.0_rk)) -   &
        (e*(-3.0_rk + Log(4.0_rk) -   &
             (Sqrt(b)*Log(4.0_rk*a +   &
                  2.0_rk*Sqrt(a)*Sqrt(4.0_rk*a + b - 2.0_rk*c) - c))/  &
              Sqrt(a) +   &
             Log(b + Sqrt(b)*Sqrt(4.0_rk*a + b - 2.0_rk*c) -   &
               c) + 3.0_rk*Log(-b +   &
                Sqrt(b)*Sqrt(4.0_rk*a + b - 2.0_rk*c) + c) +   &
             (Sqrt(b)*Log(-2.0_rk*a - c +   &
                  2.0_rk*Sqrt(a)*Sqrt(a + b + c)))/Sqrt(a)  &
              - Log(2.0_rk*b + c +   &
               2.0_rk*Sqrt(b)*Sqrt(a + b + c))))/(3.0_rk*b**1.5_rk)  
    
   end subroutine
    
!
! ------------------------------------------------------------------------------
!
   subroutine getXi000G(a,b,c,d,e,f,G)
      
      real(rk),intent(in) :: a,b,c,d,e,f 
      real(rk),intent(out) ::  G
      
      G =  &
          (-4*a*b*e + c**2*e + (2*b*(2*b**2*d + c*(-a + c)*e - b*(2*a*(e - f) + c*(d + f))))/ &
          (a + b - c) + (2*Sqrt(a)*Sqrt(b)* &
          (-2*b**2*d + (a - c)*c*e + b*(2*a*(e - f) + c*(d + f))))/(a + b - c) +  &
          (b**1.5*(4*a*b - c**2)*(d + e - f)*Log(2*a + 2*Sqrt(a)*Sqrt(a + b - c) - c))/ &
          (a + b - c)**1.5 + (-2*b**1.5*c*d - 2*a*Sqrt(b)*c*e + 4*a*b**1.5*f +  &
          Sqrt(a)*(-4*a*b + c**2)*e*(-1 + Log(2*Sqrt(a)*Sqrt(b) + c)))/Sqrt(a) +  &
          4*a*b*e*Log(2*Sqrt(a)*Sqrt(b) + c) - c**2*e*Log(2*Sqrt(a)*Sqrt(b) + c) - &
          (b**1.5*(4*a*b - c**2)*(d + e - f)*Log(-2*b + 2*Sqrt(b)*Sqrt(a + b - c) + c))/ &
          (a + b - c)**1.5)/(b**1.5*(4*a*b - c**2))  
    
   end subroutine
    
!
! ------------------------------------------------------------------------------
!
   subroutine getXi100G(a,b,c,d,e,f,G)

      real(rk),intent(in) :: a,b,c,d,e,f 
      real(rk),intent(out) ::  G
      
      G =  &
        (2*Sqrt(b)*(b*Sqrt(a + b - c)*c*d - a**1.5*(c*e + 2*b*(e - f)) + &
        a*Sqrt(a + b - c)*(c*e - 2*b*f) + Sqrt(a)*(2*b**2*d + c**2*e - b*c*(d + f))) - &
        Sqrt(a)*Sqrt(a + b - c)*(4*a*b - c**2)*e* &
        (Log(2*Sqrt(a)*Sqrt(b) - c) - Log(2*b + 2*Sqrt(b)*Sqrt(a + b - c) - c)))/ &
        (Sqrt(a)*b**1.5*Sqrt(a + b - c)*(4*a*b - c**2))
    
   end subroutine
    
    
!
! ------------------------------------------------------------------------------
!
   subroutine getXi010G(a,b,c,d,e,f,G)

      
      real(rk),intent(in) :: a,b,c,d,e,f 
      real(rk),intent(out) ::  G
      
      G =  &
          ((2*Sqrt(b)*(-2*b**2*d + (a - c)*c*e + b*(2*a*(e - f) + c*(d + f))) + &
             Sqrt(a + b - c)*(4*a*b - c**2)*e* &
              (-1 + Log(-2*b + 2*Sqrt(b)*Sqrt(a + b - c) + c)))/Sqrt(a + b - c) + &
          (b**1.5*Sqrt(a + b - c)*(-4*a*b + c**2)*d*Log(2*Sqrt(a)*Sqrt(b) - c) + &
             b**1.5*Sqrt(a + b - c)*(4*a*b - c**2)*d* &
              Log(2*a + 2*Sqrt(a)*Sqrt(a + b - c) - c) +  &
             Sqrt(a)*(2*b**1.5*c*(-b + Sqrt(b)*Sqrt(a + b - c) + c)*d +  &
                a**2*(4*b*Sqrt(a + b - c)*e - 2*Sqrt(b)*c*e + 4*b**1.5*f) +  &
                a*(2*b*Sqrt(a + b - c)*c*e + 2*Sqrt(b)*c**2*e -  &
                   Sqrt(a + b - c)*c**2*e + 4*b**2.5*f - 4*b**2*Sqrt(a + b - c)*f - &
                   2*b**1.5*c*(d + e + 2*f)) +  &
                a*Sqrt(a + b - c)*(-4*a*b + c**2)*e* &
                 Log(-2*b + 2*Sqrt(b)*Sqrt(a + b - c) + c)))/(a**1.5*Sqrt(a + b - c)) &
          )/(b**1.5*(4*a*b - c**2)) 
    
    
   end subroutine
    
!
! ------------------------------------------------------------------------------
!
   subroutine getXi13130(a,b,c,FF2)
!
!     Implements FF2 in Eq. (24) in the paper:
!     Calculates integral over a triangle, using analitical method of
!     Ren, Q., & Chan, C. L. (2015). Analytical evaluation of the BEM singular integrals 
!     for 3D Laplace and Stokes flow equations using coordinate transformation. 
!     Engineering Analysis with Boundary Elements, 53, 1–8. 
!     https://doi.org/10.1016/j.enganabound.2014.11.018
!
!
! ------------------------------------------------------------------------------
!
          
         real(rk),intent(in) :: a,b,c ! a,b,c - Eq. (23)
         real(rk),intent(out) ::  FF2 ! result - Eq. (24)
        
         real(rk) d,e,t1,t2,t3
        
         d = 4.0_rk*a + 2.0_rk*b + 2.0_rk* &
              SQRT( 4.0_rk*a + b - 2.0_rk*c )* &
              SQRT( a + b - c )  &
              -3.0_rk*c
        
         e = -2.0_rk*a - 4.0_rk*b + 2.0_rk* &
              SQRT( a + 4.0_rk*b - 2.0_rk*c )* &
              SQRT( a + b - c )  &
              +3.0_rk*c
        
         t1 = log(d/e)/SQRT(a+b-c)
        
        
         d = 4.0_rk*a + 2.0_rk* &
              SQRT( 4.0_rk*a + b - 2.0_rk*c )* &
              SQRT( a )  &
              -c
        
         e = -2.0_rk*a - c + 2.0_rk* &
              SQRT( a + b + c )* &
              SQRT( a )  
        
         t2 = log(d/e)/SQRT(a)
        
        
         d = (                                &
             4.0_rk*b + 2.0_rk*               &
             SQRT( a + 4.0_rk*b - 2.0_rk*c )* &
             SQRT( b )                        &
             -c                               &
         ) * (                                &
             2.0_rk*b + 2.0_rk*               &
             SQRT( a + b + c )*               &
             SQRT( b )                        &
             +c                               &
         )
        
         e = (                                &
             b +                              &
             SQRT( 4.0_rk*a + b - 2.0_rk*c )* &
             SQRT( b )                        &
             -c                               &
         ) * (                                &
             -b +                             &
             SQRT( 4.0_rk*a + b - 2.0_rk*c )* &
             SQRT( b )                        &
             +c                               &
         )
        
         t3 = log(d/e)/SQRT(b)
        
        
         FF2 = ( t1 + t2 + t3 ) / 3.0_rk
        
   end subroutine
    
    
!
! ------------------------------------------------------------------------------
!  
   subroutine getMrot(rr1,rr2,rr3,Mrot)
            
      real(rk),intent(in) :: rr1(3),rr2(3),rr3(3) ! vectors in local FR
      real(rk),intent(out) ::  Mrot(3,3) ! result, rotation matrix (3D)   
   
      real(rk), allocatable :: ov1(:),ov2(:),ov3(:),ov4(:) ! ortogonal vectors
      real(rk), allocatable :: e1(:),e2(:),e3(:) ! unit vectors
   
      allocate ( ov1(3),ov2(3),ov3(3),ov4(3) )
      allocate ( e1(3), e2(3), e3(3) )

!     create 2 ortogonal vectors in the plane of triangle
      ov1 = rr2 - rr1
      ov4 = rr3 - rr1 
      call siCrossProduct(ov3,ov4,ov1) ! ov3 = ov4 x ov1
      call siCrossProduct(ov2,ov3,ov1) ! ov2 = ov3 x ov1
   
      call siNormVector(ov1)
      call siNormVector(ov2)
      call siNormVector(ov3)
   
   
!     unit vectors in new coordiante system
      e1(1) = 1.0_rk; e1(2) = 0.0_rk; e1(3) = 0.0_rk
      e2(1) = 0.0_rk; e2(2) = 1.0_rk; e2(3) = 0.0_rk
      e3(1) = 0.0_rk; e3(2) = 0.0_rk; e3(3) = 1.0_rk
   
!     rotation matrix 
      call siDotProduct( ov1, e1, Mrot(1,1) )
      call siDotProduct( ov1, e2, Mrot(1,2) )
      call siDotProduct( ov1, e3, Mrot(1,3) )
   
      call siDotProduct( ov2, e1, Mrot(2,1) )
      call siDotProduct( ov2, e2, Mrot(2,2) )
      call siDotProduct( ov2, e3, Mrot(2,3) )
   
      call siDotProduct( ov3, e1, Mrot(3,1) )
      call siDotProduct( ov3, e2, Mrot(3,2) )
      call siDotProduct( ov3, e3, Mrot(3,3) )

      deallocate(ov1,ov2,ov3,ov4,e1,e2,e3)
    
   end subroutine
    
!
! ------------------------------------------------------------------------------
!  
   subroutine getABC(rrr1,rrr2,rrr3,a,b,c)
    
      real(rk),intent(in)  ::  rrr1(3),rrr2(3),rrr3(3) ! vectors
      real(rk),intent(out) ::  a,b,c ! result

      real(rk), allocatable :: Amat(:,:) ! rotation matrix (2D)
      allocate ( Amat(2,2) )  

      Amat(1,1) = rrr1(1) - rrr3(1)      
      Amat(1,2) = rrr2(1) - rrr3(1)
      Amat(2,1) = rrr1(2) - rrr3(2)
      Amat(2,2) = rrr2(2) - rrr3(2)

      a = Amat(1,1)**2 + Amat(2,1)**2
      b = Amat(1,2)**2 + Amat(2,2)**2 
      c = 2.0_rk * (Amat(1,1)*Amat(1,2) + Amat(2,1)*Amat(2,2))   

      deallocate(Amat)
        
   end subroutine
!
! ------------------------------------------------------------------------------
!  
   subroutine getDEF(rrr1,rrr2,rrr3,dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy)
    
    
      real(rk),intent(in)  ::  rrr1(3),rrr2(3),rrr3(3) ! vectors
      real(rk),intent(out) ::  dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy ! result

      real(rk), allocatable :: Amat(:,:) ! rotation matrix (2D)

      ! map triangle on (x,y) plane to reference triangle
      ! 2D rotation matrix to rotate to reference triangle
        
      allocate ( Amat(2,2) )  

      Amat(1,1) = rrr1(1) - rrr3(1)      
      Amat(1,2) = rrr2(1) - rrr3(1)
      Amat(2,1) = rrr1(2) - rrr3(2)
      Amat(2,2) = rrr2(2) - rrr3(2)

      dxx = Amat(1,1)**2
      exx = Amat(1,2)**2 
      fxx = 2.0_rk * Amat(1,1)*Amat(1,2)
      dyy = Amat(2,1)**2
      eyy = Amat(2,2)**2
      fyy = 2.0_rk * Amat(2,1)*Amat(2,2)
      dxy = Amat(1,1)*Amat(2,1) 
      exy = Amat(1,2)*Amat(2,2) 
      fxy = Amat(1,1)*Amat(2,2)+Amat(2,1)*Amat(1,2)
      
      deallocate(Amat)
        
   end subroutine
    
!
! ------------------------------------------------------------------------------
!  
   subroutine getrcABC(x1,y1,z1,x2,y2,z2,x3,y3,z3,fa,fb,fc)
!
! ------------------------------------------------------------------------------
!  
        
      real(rk),intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
      real(rk),intent(out) ::  fa,fb,fc ! result
   
      real(rk), allocatable :: r1(:),r2(:),r3(:),c(:) ! vertexes and center
      real(rk), allocatable :: rr1(:),rr2(:),rr3(:) ! vertexes after translation
      real(rk), allocatable :: rrr1(:),rrr2(:),rrr3(:) ! vertexes after rotation
      real(rk), allocatable :: Mrot(:,:) ! rotation matrix (3D)
   
      allocate ( r1(3),r2(3),r3(3),c(3) )
      allocate ( rr1(3),rr2(3),rr3(3) )
      allocate ( rrr1(3),rrr2(3),rrr3(3) )  

      allocate ( Mrot(3,3) )  
   
   
      ! store vertexes in vector form
      r1(1) = x1 ; r2(1) = x2 ; r3(1) = x3
      r1(2) = y1 ; r2(2) = y2 ; r3(2) = y3
      r1(3) = z1 ; r2(3) = z2 ; r3(3) = z3
   
      ! element barycenter 
      c = ( r1 + r2 + r3 ) / 3.0_rk
   
      ! translation - set origin to Barycentric coordinates = 1/3
      rr1 = r1 - c
      rr2 = r2 - c
      rr3 = r3 - c

      ! get 3D rotation matrix
      call getMrot(rr1,rr2,rr3,Mrot)
   
      ! rotate - all vertexes in the rotated FR shoule be in (x,y) plane
      rrr1 = MATMUL(Mrot,rr1)
      rrr2 = MATMUL(Mrot,rr2)
      rrr3 = MATMUL(Mrot,rr3)
   
      ! map triangle on (x,y) plane to reference triangle
      ! 2D rotation matrix to rotate to reference triangle
      call getABC(rrr1,rrr2,rrr3,fa,fb,fc)
        
      deallocate (r1,r2,r3,c,rr1,rr2,rr3,rrr1,rrr2,rrr3)
      deallocate (Mrot)
        
        
   end subroutine
        
!
! ------------------------------------------------------------------------------
!  
   subroutine getABCDEF(x1,y1,z1,x2,y2,z2,x3,y3,z3,fa,fb,fc,dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy,Mrot)
!
! ------------------------------------------------------------------------------
!  
      real(rk),intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3 ! triangle vertexes
      real(rk),intent(out) :: Mrot(3,3) ! rotation matrix (3D)
      real(rk),intent(out) :: fa,fb,fc,dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy ! result
    
      real(rk), allocatable :: r1(:),r2(:),r3(:),c(:) ! vertexes and center
      real(rk), allocatable :: rr1(:),rr2(:),rr3(:) ! vertexes after translation
      real(rk), allocatable :: rrr1(:),rrr2(:),rrr3(:) ! vertexes after rotation
    
      allocate ( r1(3),r2(3),r3(3),c(3) )
      allocate ( rr1(3),rr2(3),rr3(3) )
      allocate ( rrr1(3),rrr2(3),rrr3(3) )   
    
    
      ! store vertexes in vector form
      r1(1) = x1 ; r2(1) = x2 ; r3(1) = x3
      r1(2) = y1 ; r2(2) = y2 ; r3(2) = y3
      r1(3) = z1 ; r2(3) = z2 ; r3(3) = z3
    
      ! element barycenter 
      c = ( r1 + r2 + r3 ) / 3.0_rk
    
      ! translation - set origin to Barycentric coordinates = 1/3
      rr1 = r1 - c
      rr2 = r2 - c
      rr3 = r3 - c
      ! get 3D rotation matrix
      call getMrot(rr1,rr2,rr3,Mrot)
    
      ! rotate - all vertexes in the rotated FR shoule be in (x,y) plane
      rrr1 = MATMUL(Mrot,rr1)
      rrr2 = MATMUL(Mrot,rr2)
      rrr3 = MATMUL(Mrot,rr3)
    
      ! map triangle on (x,y) plane to reference triangle
      ! 2D rotation matrix to rotate to reference triangle
      call getABC(rrr1,rrr2,rrr3,fa,fb,fc)
      call getDEF(rrr1,rrr2,rrr3,dxx,exx,fxx,dxy,exy,fxy,dyy,eyy,fyy)
        
      deallocate (r1,r2,r3,c,rr1,rr2,rr3,rrr1,rrr2,rrr3)
      
   end subroutine
          


!
! ------------------------------------------------------------------------------
!  
   subroutine siCrossProduct(a,b,c)
!  a = b x c

      real(rk) a(3),b(3),c(3) 

      a(1)= b(2)*c(3)-b(3)*c(2)
      a(2)=-b(1)*c(3)+b(3)*c(1)
      a(3)= b(1)*c(2)-b(2)*c(1)
      
   end subroutine
!
! ------------------------------------------------------------------------------
!  
   subroutine siDotProduct(a,b,dp)
!   dp = a \cdot b
           
      real(rk) a(3),b(3),dp
      
      dp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      
   end subroutine
      

!
! ------------------------------------------------------------------------------
! 
   subroutine siNormVector(a)
!  a=a / norm(a)
               
      real(rk) a(3),norm
      integer i

      do i=1,3
        norm=norm+a(i)**2
      end do
      norm=SQRT(norm)
      do i=1,3
        a(i)=a(i)/norm
      end do

   end subroutine
               

end module