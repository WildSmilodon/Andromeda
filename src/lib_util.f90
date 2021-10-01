! -----------------------------------------------------------------------------
      SUBROUTINE DatumInUra(cas)
!
!     $: writes date and time to cas
!
! -----------------------------------------------------------------------------
      CHARACTER*50 D,T,Z
      INTEGER Value(8)
      CHARACTER*39 cas

      CALL DATE_AND_TIME(D,T,Z,Value)
      WRITE (cas,33) &
         'Date and time = ',Value(3),'.',Value(2),'.',Value(1),' ', &
         Value(5),':',Value(6),':',Value(7),'.',Value(8)
33    FORMAT (A16,I2,A1,I2,A1,I4,A1,I2,A1,I2,A1,I2,A1,I3)

      END



! -----------------------------------------------------------------------------
      function Cstring(strIn,cStr) result(iOut)
!
!     $: comapres two strings
!
! -----------------------------------------------------------------------------

      USE mPar

      implicit none

      character(len=*), intent(in) :: strIn
      character(len=*), intent(in) :: cstr
      character(len=len(strIn)) :: strOut
      integer :: i,j,iOut

      do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
      end do

      IF (TRIM(strOut).EQ.TRIM(cStr)) THEN
        iOut=parYes
      ELSE
        iOut=parNo
      END IF


      end function Cstring



! ********************************************************************
      SUBROUTINE GetRotationMatrix(R,RT,n1)
!     Creates a rotation matrix R, so that
!     vector n1 is rotated into e1=(1,0,0)
!     e1=R*n1
!     n1=RT*e1
            USE mCommon
            implicit none
      REAL(rk) e1(3),e2(3),n1(3),n2(3),tmp(3),dp
      REAL(rk) e3(3),e4(3),n3(3),n4(3)
      REAL(rk) Re(3,3),RnT(3,3),R(3,3),RT(3,3)
!
      e1(1)=1.0_rk
      e1(2)=0.0_rk
      e1(3)=0.0_rk
!
      e2(1)=0.0_rk
      e2(2)=1.0_rk
      e2(3)=0.0_rk

!     Create a vector n2 perpendicular to vector n1
!     tmp is a vector, which is not colilear with n1
      tmp(1)=n1(2)+3.345345
      tmp(2)=n1(3)-7.345344
      tmp(3)=n1(1)+6.243234

      CALL CrossProduct(n2,n1,tmp)
      CALL NormVector(n2)
!     Dot product of perpendicular vectors should be zero
      CALL DotProduct(n1,n2,dp)
      IF (ABS(dp).GT.1.0D-14) THEN
        PRINT *,"ERROR IN GetRotationMatrix"
        STOP
      END IF
!     Use TRIAD algorithm http://en.wikipedia.org/wiki/Triad_Method
      CALL CrossProduct(n3,n1,n2)
      CALL NormVector(n3)
      CALL CrossProduct(n4,n1,n3)
      CALL NormVector(n4)
      CALL CrossProduct(e3,e1,e2)
      CALL NormVector(e3)
      CALL CrossProduct(e4,e1,e3)
      CALL NormVector(e4)

      RnT(1,1)=n1(1)
      RnT(1,2)=n1(2)
      RnT(1,3)=n1(3)

      RnT(2,1)=n3(1)
      RnT(2,2)=n3(2)
      RnT(2,3)=n3(3)

      RnT(3,1)=n4(1)
      RnT(3,2)=n4(2)
      RnT(3,3)=n4(3)

      Re(1,1)=e1(1)
      Re(2,1)=e1(2)
      Re(3,1)=e1(3)

      Re(1,2)=e3(1)
      Re(2,2)=e3(2)
      Re(3,2)=e3(3)

      Re(1,3)=e4(1)
      Re(2,3)=e4(2)
      Re(3,3)=e4(3)

!     Rotational matrix
      R=MATMUL(Re,RnT)

!     Transpose rotation matrix
      RT(1,1)=R(1,1)
      RT(1,2)=R(2,1)
      RT(1,3)=R(3,1)
      RT(2,1)=R(1,2)
      RT(2,2)=R(2,2)
      RT(2,3)=R(3,2)
      RT(3,1)=R(1,3)
      RT(3,2)=R(2,3)
      RT(3,3)=R(3,3)

      END

! ********************************************************************
      SUBROUTINE CrossProduct(a,b,c)
!     a=b x c
      USE mCommon
      implicit none
      REAL(rk) a(3),b(3),c(3)

      a(1)= b(2)*c(3)-b(3)*c(2)
      a(2)=-b(1)*c(3)+b(3)*c(1)
      a(3)= b(1)*c(2)-b(2)*c(1)

      END

! ********************************************************************
      SUBROUTINE DotProduct(a,b,dp)
!     dp = a \cdot b
      USE mCommon
      implicit none            
      REAL(rk) a(3),b(3),dp

      dp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

      END

! ********************************************************************
      SUBROUTINE Distance(a,b,d)
!     d = euclidian distance from a to b
      USE mCommon
      implicit none            
      REAL(rk) a(3),b(3),d

      d= SQRT ( (a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2 )

      END




! ********************************************************************
      SUBROUTINE NormVector(a)
!     a=a / norm(a)

      USE mCommon
      implicit none            
      REAL(rk) a(3),norm
      INTEGER i

      DO i=1,3
        norm=norm+a(i)**2
      END DO
      norm=SQRT(norm)
      DO i=1,3
        a(i)=a(i)/norm
      END DO

      END


!C______________________________________________________________________C
!C______________________________________________________________________C
      SUBROUTINE WrSMat(mat,nrow,ncol,io)
!        __    _      ___
!        Write Single Matrix
!C______________________________________________________________________C
!C______________________________________________________________________C
            USE mCommon
            implicit none            
      INTEGER nrow,ncol,io,j
      REAL(rk) mat(nrow,ncol)

      DO j=1,ncol
        CALL wrvec(io,nrow,mat(1,j))
      END DO

      END
!C______________________________________________________________________C
!C______________________________________________________________________C
      SUBROUTINE RdSMat(mat,nrow,ncol,io)
!        _  _ _      ___
!        Read Single Matrix
!C______________________________________________________________________C
!C______________________________________________________________________C
            USE mCommon
            implicit none            
      INTEGER nrow,ncol,io,j
      REAL(rk) mat(nrow,ncol)

      DO j=1,ncol
        CALL rdvec(io,nrow,mat(1,j))
      END DO

      END

!----------------------------------------------------------------------!
!c----------------------------------------------------------------------c
!c                                                                      c
      SUBROUTINE rdvec(ifr,nnx,vec)
!c                                                                      c
!c----------------------------------------------------------------------c
!----------------------------------------------------------------------!
!c.......................................................................
!c..                                                                   ..
!c..   REad VECtor using MAXSIZE chunks                                ..
!c..   --   ---                                                        ..
!c.......................................................................
            USE mCommon
            implicit none            
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      REAL(rk)  vec(nnx)
      PARAMETER (maxsize=8192)

      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        READ(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      READ(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END
!----------------------------------------------------------------------!
!c----------------------------------------------------------------------c
!c                                                                      c
      SUBROUTINE wrvec(ifr,nnx,vec)
!c                                                                      c
!c----------------------------------------------------------------------c
!----------------------------------------------------------------------!
!c.......................................................................
!c..                                                                   ..
!c..   WRite VECtor using MAXSIZE chunks                               ..
!c..   --    ---                                                       ..
!c.......................................................................
            USE mCommon
            implicit none            
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      REAL(rk)  vec(nnx)
      PARAMETER (maxsize=8192)
!c....
      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        WRITE(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      WRITE(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END



!----------------------------------------------------------------------!

      REAL(rk) FUNCTION cptime (rvar)

!----------------------------------------------------------------------!
! **********************************************************************
! **                                                                  **
! ** Returns the CPU TIME in seconds                                  **
! **             --  ----                                             **
! **********************************************************************
!
!     with Intel compiler Makefile must have -fpp
!     with gfortran compiler Makefile must have -cpp
!
      USE mCommon
      implicit none
      REAL(rk) rvar,tmp
      CALL CPU_TIME(tmp)
      cptime = tmp - rvar

      RETURN
      END
