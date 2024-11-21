C----------------------------------------------------------------------C
c
      SUBROUTINE SolvEQNfm(slvt,pret,prep,maxit,stopt,eps,neq,q,a,b,x,nit,cpu,ierr)
C
C----------------------------------------------------------------------C
C **********************************************************************
C **                                                                  **
C ** Solve system of equations using various solvers (full matrix)    **
C ** ----            --     -                         -    -          **
C **                                                                  **
C ** slvt ... solver type                                             **
C **          0 direct                                                **
C ** pret ... preconditioner Q type                                   **
C **          0 none                                                  **
C **          1 diag                                                  **
C **          2 lu                                                    **
C ** prep ... preconditioner Q place of calculation                   **
C **          1 inside                                                **
C **          2 outside                                               **
C ** neq  ... number of equations                                     **
C ** maxit .. maximum number of iterations                            **
C ** stopt .. stopping criterum type                                  **
C ** eps  ... stopping criterium value                                **
C ** q  ..... preconditioner vector Q                                 **
C ** a  ..... matrix A                                                **
C ** b  ..... right hand side vector                                  **
C ** x  ..... vector of unknowns                                      **
C ** nit .... number of iterations                                    **
C ** cpu .... elapsed CPU time                                        **
C ** ierr ... error status (0=no error)                               **
C **                                                                  **
C **********************************************************************
      USE mPar  
      USE mCommon
      IMPLICIT NONE
C Arguments 
      INTEGER   slvt,pret,prep,maxit,stopt,neq,nit,ierr
      REAL(rk)  eps
      REAL(rk)  q(neq),a(neq,neq),b(neq),x(neq)
C Internal
      REAL(rk)  cpu,cpu0,cptime
*     .. External Functions ..
      REAL(rk)  DNRM2
      EXTERNAL DNRM2
C
C Measure CPU time
C
      cpu0=cptime(0.0_rk)
C
C Check for b=0
C
       IF (DNRM2(neq,b,1).EQ.0.0_rk) THEN
         CALL dinit(neq,0.0D0,x,1)
         nit=0
         ierr=0
         RETURN
       END IF
c
c Calculate preconditioner
c
      IF (pret.NE.0 .AND. prep.EQ.1) THEN
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdiafm(neq,1,a,q)
        ELSE IF (pret.EQ.2) THEN
c lu
          CALL dgefa(a,neq,neq,q,ierr)
        ELSE
c error
          WRITE (parLogTekst,'(I0,A,A,I0)') 2,'SolvEQNfm','Wrong preconditioner type!',0
          CALL WriteToLog(parLogTekst)
        END IF
      END IF
c
c Run solver
c
      IF (slvt.EQ.0) THEN
c Direct solver
        x=b
        CALL dgesl(a,neq,neq,q,x,0)
      ELSE
c error
        WRITE (parLogTekst,'(I0,A,A,I0)') 2,'SolvEQNfm','Wrong solver type!',0
        CALL WriteToLog(parLogTekst)
      END IF
c
      IF (ierr.EQ.-1) THEN
        WRITE (parLogTekst,'(I0,A,A,I0)') 2,'SolvEQNfm','Solver reached maximum number of iterations!',0
        CALL WriteToLog(parLogTekst)        
c     ELSE IF (ierr.LE.-2) THEN
c       WRITE(*,*) 'Warning: SolvSLE: solver finished with error',ierr,'.'
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END


!----------------------------------------------------------------------C
!
      SUBROUTINE FormPREfm(pret,neq,q,a,cpu,ierr)
!
!----------------------------------------------------------------------C
! **********************************************************************
! **                                                                  **
! ** Form Preconditioning Matrix Q (full matrix A)                    **
! ** ---- ---                       -    -                            **
! **                                                                  **
! ** pret ... preconditioner Q type                                   **
! **          0 none                                                  **
! **          1 diag                                                  **
! **          2 lu                                                    **
! ** neq  ... number of equations                                     **
! ** q    ... preconditioner matrix Q                                 **
! ** a    ... matrix A                                                **
! ** cpu .... elapsed CPU time                                        **
! ** ierr ... error status (0=no error)                               **
! **                                                                  **
! **********************************************************************
      USE mPar
      USE mCommon
      IMPLICIT NONE
! Arguments
      INTEGER pret,neq,ierr
      REAL    cpu
      REAL(rk)  q(neq),a(neq,neq)
! Internal
      REAL(rk)    cpu0,cptime
c
c Measure CPU time
c
      cpu0=cptime(0.0_rk)
c
c Calculate preconditioner
c
      IF (pret.NE.0) THEN
        IF (pret.EQ.1) THEN
c diagonal
          CALL frmdiafm(neq,1,a,q)
        ELSE IF (pret.EQ.2) THEN
c lu
          CALL dgefa(a,neq,neq,q,ierr)
        ELSE
c error
          WRITE (parLogTekst,'(I0,A,A,I0)') 2,'FormPREfm','Wrong preconditioner type!',0
          CALL WriteToLog(parLogTekst)
       
        END IF
      END IF
c
      cpu=cpu+cptime(cpu0)
c
      RETURN
      END
C^L

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE frmdiafm(neq,pres,a,q)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Form Preconditioner with Diagonal matrix (full matrix)          **
c **  - --                     ---              -    -                **
c **                                                                  **
c **  neq  ... number of equations                                    **
c **  pres ... preconditioner side (1=left,2=right,3=both)            **
c **  a    ... matrix A                                               **
c **  q  ..... preconditioner vector Q                                **
c **                                                                  **
c **********************************************************************
      USE mCommon
      IMPLICIT NONE
c Arguments
      INTEGER neq,pres
      REAL(rk)  a(neq,neq),q(neq)
c Internal
      INTEGER i
c
      IF (pres.EQ.1.OR.pres.EQ.2) THEN
        DO i = 1,neq
          q(i) = 1.D0/a(i,i)
        END DO
      ELSE IF (pres.EQ.3) THEN
        DO i = 1,neq
           q(i) = 1.D0/SQRT(a(i,i))
        END DO
      END IF
c
      RETURN
      END


C----------------------------------------------------------------------C
c
      subroutine dgefa(a,lda,n,ipvt,info)
c
C----------------------------------------------------------------------C
c***begin prologue  dgefa
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  double precision,factor,linear algebra,linpack,matrix
c***author  moler, c. b., (u. of new mexico)
c***purpose  factors a double precision matrix by gaussian elimination.
c***description
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    double precision(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,dscal,idamax
c***end prologue  dgefa
      USE mCommon
      IMPLICIT NONE
      integer lda,n,info
      real(rk) ipvt(*),a(lda,*)
c
      real(rk) t
      integer idamax,j,k,kp1,l,nm1
c
c     gaussian elimination with partial pivoting
c
c***first executable statement  dgefa
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30       continue
         go to 50
40    continue
            info = k
50    continue
60    continue
70    continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end      

C----------------------------------------------------------------------C
c
      subroutine dgesl(a,lda,n,ipvt,b,job)
c
C----------------------------------------------------------------------C
c***begin prologue  dgesl
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  double precision,linear algebra,linpack,matrix,solve
c***author  moler, c. b., (u. of new mexico)
c***purpose  solves the double precision system  a*x=b or  trans(a)*x=b
c            using the factors computed by dgeco or dgefa.
c***description
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    double precision(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,ddot
c***end prologue  dgesl
      USE mCommon
      IMPLICIT NONE
      integer lda,n,job
      real(rk) ipvt(*),a(lda,*),b(*)
c
      real(rk) ddot,t
      integer k,kb,l,nm1
c***first executable statement  dgesl
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
20    continue
30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
40         continue
      go to 100
50    continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
70       continue
80    continue
90    continue
100   continue
      return
      end
   

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE dinit(n,alpha,dx,incx)
c                                                                      c
C----------------------------------------------------------------------C
c
*
*     Initialises a vector x with a scalar alpha.
*     Modified from dcopy, BLAS Level 1.
*     Rudnei Dias da Cunha, 14/6/93.
*

*     copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*
      USE mCommon
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      REAL(rk)  ALPHA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL(rk)  DX(n)
*     ..
*     .. Local Scalars ..
      INTEGER I,IX,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
          DX(IX) = ALPHA
          IX = IX + INCX
10    CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
20    M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = ALPHA
30    CONTINUE
      IF (N.LT.7) RETURN
40    MP1 = M + 1
      DO 50 I = MP1,N,7
          DX(I) = ALPHA
          DX(I+1) = ALPHA
          DX(I+2) = ALPHA
          DX(I+3) = ALPHA
          DX(I+4) = ALPHA
          DX(I+5) = ALPHA
          DX(I+6) = ALPHA
50    CONTINUE
      RETURN
      END      

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE prediafm(neq,nnz,iq,jq,dq,q,u,v,ipar)

      USE mCommon
      IMPLICIT NONE

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq),ipar(13)
      REAL(rk)  q(neq),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with diagonal preconditioner                **
c **  ---                      ---                                    **
c **********************************************************************
      v=q*u
c
      RETURN
      END      


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE progress(loclen,itno,normres,x,res,trueres)
c                                                                      c
C----------------------------------------------------------------------C
c
      USE mCommon
      IMPLICIT NONE

*     .. Scalar Arguments ..
      INTEGER ITNO,LOCLEN
      REAL(rk)  NORMRES
*     ..
*     .. Array Arguments ..
      REAL(rk)  RES(loclen),TRUERES(loclen),X(loclen)
*     ..
*     .. External Subroutines ..
      EXTERNAL PRINTV
*     ..
*     WRITE (6,FMT=9000) ITNO,NORMRES
*     WRITE (6,FMT=9010) 'X:'
*     CALL PRINTV(LOCLEN,X)
*     WRITE (6,FMT=9010) 'RES:'
*     CALL PRINTV(LOCLEN,RES)
*     WRITE (6,FMT=9010) 'TRUE RES:'
*     CALL PRINTV(LOCLEN,TRUERES)
      RETURN
 9000 FORMAT (I5,1X,D16.10)
 9010 FORMAT (A)
      END      

C----------------------------------------------------------------------C
c                                                                      c
      REAL(rk)  FUNCTION pdnrm2(loclen,u,ipar)
c                                                                      c
C----------------------------------------------------------------------C
c
      USE mCommon
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER LOCLEN
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL(rk)  U(loclen)
*     ..
*     .. External Functions ..
      REAL(rk)  DNRM2
      EXTERNAL DNRM2
*     ..
      PDNRM2 = DNRM2(LOCLEN,U,1)
*     ..
      RETURN
      END      

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pdsum(isize,x,ipar)
c                                                                      c
C----------------------------------------------------------------------C
c
      USE mCommon
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER ISIZE
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL(rk)  X(isize)
*     ..
      RETURN
      END      

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE matvecfm(neq,nnz,ia,ja,a,u,v,ipar)
c                                                                      c
C----------------------------------------------------------------------C
      USE mCommon
      IMPLICIT NONE

      INTEGER neq,nnz,ia(neq+1),ja(nnz),ipar(13)
      REAL(rk)  a(neq,neq),u(neq),v(neq)

      INTEGER i,j,n
c
      n = ipar(2)
      v = 0.0D0
      DO j = 1,n
        DO i = 1,n
          v(i) = v(i) + a(i,j)*u(j)
        END DO
      END DO
c
      RETURN
      END      

C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pimRBiCGSTAB
     &   (neq,nnz,pret,maxit,stopt,eps,nit,ierr,iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using RBi-CGSTAB iterative method      **
c ** - --                            --- -----                        **
c **********************************************************************
      USE mCommon
      IMPLICIT NONE
c Parameters
      INTEGER  neq,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(neq+1),jq(nnz),dq(neq),
     &         ia(neq+1),ja(nnz),da(neq)
      REAL(rk)   q(nnz),a(nnz),b(neq),x(neq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,pdsum,pdnrm,progress
c Internal variables
      INTEGER  nwrk,ipar(13),basis,pres
      REAL(rk)   wrk,dpar(6)
      EXTERNAL dinit,dvprod,rbicgstab
c
      ALLOCATABLE wrk(:)
c
c Set initial values
c
c if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c basis dimension for RBi-CGSTAB vectors
      basis=2
c
c Allocate memory for matrix and vectors
c
      nwrk = (6+2*basis)*neq
      ALLOCATE(wrk(nwrk))
      CALL pimdsetpar(ipar,dpar,neq,neq,neq,neq,basis,-1,-1,
     &                pres,stopt,maxit,eps)
c
c RBi-CGSTAB
c
      CALL rbicgstab(neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &               matvec,preconl,preconr,pdsum,pdnrm,progress)
c
c Set output values
c
      nit=ipar(11)
      ierr=ipar(12)
c
c Free memory
c
      DEALLOCATE(wrk)
c
      RETURN
      END
C^L      

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pimdsetpar
     &   (IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,PROCID,
     &    PRECONT,STOPT,MAXIT,EPSILON)
c                                                                      c
C----------------------------------------------------------------------C
C     ..
      USE mCommon
      IMPLICIT NONE
C     .. Parameters ..
      REAL(rk)  ONE
      PARAMETER (ONE=1.0_rk)
C     ..
C     .. Arguments ..
      INTEGER IPAR(13),BASIS,BLKSZ,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,
     &        PROCID,STOPT
      REAL*8  DPAR(6),EPSILON
C     ..
      ipar(1) = LDA
      ipar(2) = N
      ipar(3) = BLKSZ
      ipar(4) = LOCLEN
      ipar(5) = BASIS
      ipar(6) = NPROCS
      ipar(7) = PROCID
      ipar(8) = PRECONT
      ipar(9) = STOPT
      ipar(10) = MAXIT
      ipar(11) = -1
      ipar(12) = -1
      ipar(13) = -1

      dpar(1) = EPSILON
      dpar(2) = -ONE

      RETURN
      END

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE pimdgetpar
     &   (IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,PROCID,
     &    PRECONT,STOPT,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)
c                                                                      c
C----------------------------------------------------------------------C
C     ..
      USE mCommon
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER IPAR(13),BASIS,BLKSZ,ITNO,LDA,LOCLEN,MAXIT,N,NPROCS,
     &        PRECONT,PROCID,STATUS,STEPERR,STOPT
      REAL(rk) DPAR(6),EPSILON,EXITNORM
C     ..
      LDA = ipar(1)
      N = ipar(2)
      BLKSZ = ipar(3)
      LOCLEN = ipar(4)
      BASIS = ipar(5)
      NPROCS = ipar(6)
      PROCID = ipar(7)
      PRECONT = ipar(8)
      STOPT = ipar(9)
      MAXIT = ipar(10)
      ITNO = ipar(11)
      STATUS = ipar(12)
      STEPERR = ipar(13)

      EPSILON = dpar(1)
      EXITNORM = dpar(2)

      RETURN
      END      

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE rbicgstab
     &   (neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c
C----------------------------------------------------------------------C
c
      USE mCommon
      IMPLICIT NONE
*     .. Parameters ..
      REAL(rk)  ZERO
      PARAMETER (ZERO=0.0_rk)
      REAL(rk)  ONE
      PARAMETER (ONE=1.0_rk)
      INTEGER IBDIM
      PARAMETER (IBDIM=8)
*     ..
*     .. Matrices Q,A ..
      INTEGER neq,nnz,nwrk,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz)
      REAL(rk)  q(nnz),a(nnz)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL(rk)  b(neq),x(neq),wrk(nwrk),dpar(6)
*     ..
*     .. Function Arguments ..
      REAL(rk)   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PDSUM,PRECONL,PRECONR,PROGRESS
*     ..
*     .. Local Scalars ..
      INTEGER BASIS,BLKSZ,CNVRTX,I,I0,I1,I2,I3,I4,IR,IRTILDE,ITNO,IU,
     +        IW,IXOLD,IZ,J,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,PROCID,
     +        STATUS,STEPERR,STOPT
      REAL*8  ALPHA,BETA,EPSILON,EXITNORM,KSI,OMEGA,RHO0,RHO1,RHSSTOP,S
*     ..
*     .. Local Arrays ..
      REAL*8  DOTS(IBDIM),GAMMA(IBDIM),GAMMA1(IBDIM),
     +        GAMMA2(IBDIM),SIGMA(IBDIM),TAU(IBDIM,IBDIM)
*     ..
*     .. External Functions ..
      REAL*8   DDOT,DSETRHSSTOP
      EXTERNAL DDOT,DSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DINIT,PIMDGETPAR,STOPCRIT
*     ..
      CALL PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,
     +                PROCID,PRECONT,STOPT,MAXIT,ITNO,STATUS,
     +                STEPERR,EPSILON,EXITNORM)

*  Check consistency of preconditioning and stop types
      IF (((PRECONT.EQ.0).OR. (PRECONT.EQ.2)) .AND. (STOPT.EQ.6)) THEN
          ITNO = 0
          STATUS = -4
          STEPERR = 0
          GO TO 9999
      END IF

*  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

*  Set indices for mapping local vectors into wrk
      IRTILDE = 1
      IW = IRTILDE + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN
      IR = IXOLD + LOCLEN
      IU = IR + (BASIS+1)*LOCLEN

*  Set rhs of stopping criteria
      RHSSTOP = DSETRHSSTOP(NEQ,NNZ,IQ,JQ,DQ,Q,B,WRK(IR),EPSILON,IPAR,
     +                      PRECONL,PDNRM)

*  1. r=Q1(b-AQ2x)
      IF (PRECONT.EQ.0) THEN
*     r=b-Ax
          CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IR),1)

      ELSE IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
          CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IR),ipar)

      ELSE IF (PRECONT.EQ.2) THEN
*     r=b-AQ2x
          CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),ipar)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

      ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
          CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IR),ipar)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IR),WRK(IZ),ipar)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR),ipar)
      END IF

*  2. rtilde=r
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)

*  3. u0=0
      CALL DINIT(LOCLEN,ZERO,WRK(IU),1)

*  4. rho0=1, alpha=0, omega=1
      RHO0 = ONE
      ALPHA = ZERO
      OMEGA = ONE

*  Loop
      STATUS = 0
      STEPERR = -1
      EXITNORM = -ONE
      DO 120 ITNO = 1,MAXIT

*  5. rho0=-omega*rho0
          RHO0 = -OMEGA*RHO0

*  BiCG loop
          I1 = 0
          I2 = LOCLEN
          DO 30 J = 0,BASIS - 1

*  6. rho1=r(j)^{T}rtilde
              DOTS(1) = DDOT(LOCLEN,WRK(IR+I1),1,WRK(IRTILDE),1)
              CALL PDSUM(1,DOTS,ipar)
              RHO1 = DOTS(1)

*  7. beta=alpha*rho1/rho0
              IF (RHO0.EQ.ZERO) THEN
                  STATUS = -3
                  STEPERR = 7
                  GO TO 9999

              END IF

              BETA = ALPHA*RHO1/RHO0

*  8. rho0=rho1
              RHO0 = RHO1

*  9. u(i)=r(i)-beta*u(i), i=0:j
              I3 = 0
              DO 10 I = 0,J
                  CALL DCOPY(LOCLEN,WRK(IU+I3),1,WRK(IZ),1)
                  CALL DCOPY(LOCLEN,WRK(IR+I3),1,WRK(IU+I3),1)
                  CALL DAXPY(LOCLEN,-BETA,WRK(IZ),1,WRK(IU+I3),1)
                  I3 = I3 + LOCLEN
   10         CONTINUE

* 10. u(j+1)=Q1AQ2u(j)
              IF (PRECONT.EQ.0) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IU+I1),WRK(IU+I2),ipar)

              ELSE IF (PRECONT.EQ.1) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IU+I1),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU+I2),ipar)

              ELSE IF (PRECONT.EQ.2) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IU+I1),WRK(IW),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IU+I2),ipar)

              ELSE IF (PRECONT.EQ.3) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IU+I1),WRK(IZ),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU+I2),ipar)
              END IF

* 11. ksi=u(j+1)^{T}rtilde
              DOTS(1) = DDOT(LOCLEN,WRK(IU+I2),1,WRK(IRTILDE),1)
              CALL PDSUM(1,DOTS,ipar)
              KSI = DOTS(1)

* 12. alpha=rho0/ksi
              IF (KSI.EQ.ZERO) THEN
                  STATUS = -3
                  STEPERR = 12
                  GO TO 9999
              END IF

              ALPHA = RHO0/KSI

* 13. r(i)=r(i)-alpha*u(i+1), i=0:j
              I3 = 0
              I4 = LOCLEN
              DO 20 I = 0,J
                  CALL DAXPY(LOCLEN,-ALPHA,WRK(IU+I4),1,WRK(IR+I3),1)
                  I3 = I3 + LOCLEN
                  I4 = I4 + LOCLEN
   20         CONTINUE

* 14. r(j+1)=Q1AQ2r(j)
              IF (PRECONT.EQ.0) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IR+I1),WRK(IR+I2),ipar)

              ELSE IF (PRECONT.EQ.1) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IR+I1),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR+I2),ipar)

              ELSE IF (PRECONT.EQ.2) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IR+I1),WRK(IW),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IR+I2),ipar)

              ELSE IF (PRECONT.EQ.3) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IR+I1),WRK(IZ),ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),ipar)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IR+I2),ipar)
              END IF

* 15. x0=x0+alpha*u0
              CALL DCOPY(LOCLEN,X,1,WRK(IXOLD),1)
              CALL DAXPY(LOCLEN,ALPHA,WRK(IU),1,X,1)
              I1 = I1 + LOCLEN
              I2 = I2 + LOCLEN
   30     CONTINUE

* 16. check stopping criterion
          CALL STOPCRIT(NEQ,NNZ,IQ,JQ,DQ,Q,IA,JA,A,B,
     +                  WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),
     +                  RHSSTOP,CNVRTX,EXITNORM,STATUS,IPAR,
     +                  MATVEC,MATVEC,PRECONR,PDSUM,PDNRM)

*  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

          IF (STATUS.EQ.-5) THEN
              STEPERR = 16
              GO TO 9999
          ELSE IF (STATUS.EQ.0) THEN
              GO TO 9999
          END IF
*  MR loop

* 17. sigma(1)=r(1)^{T}r(1), gamma'(1)=r(0)^{T}r(1)/sigma(1)
          DOTS(1) = DDOT(LOCLEN,WRK(IR+LOCLEN),1,WRK(IR+LOCLEN),1)
          DOTS(2) = DDOT(LOCLEN,WRK(IR),1,WRK(IR+LOCLEN),1)
          CALL PDSUM(2,DOTS,ipar)
          SIGMA(1) = DOTS(1)

          IF (SIGMA(1).EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 17
              GO TO 9999
          END IF

          GAMMA1(1) = DOTS(2)/SIGMA(1)

          I0 = LOCLEN + LOCLEN
          DO 60 J = 2,BASIS

* 18. tau(i,j)=r(j)^{T}r(i)/sigma(i), r(j)=r(j)-tau(i,j)r(i)
              I1 = LOCLEN
              DO 40 I = 1,J - 1
                  DOTS(I) = DDOT(LOCLEN,WRK(IR+I0),1,WRK(IR+I1),1)
                  I1 = I1 + LOCLEN
   40         CONTINUE
              CALL PDSUM(J-1,DOTS,ipar)
              I1 = LOCLEN
              DO 50 I = 1,J - 1
                  TAU(I,J) = DOTS(I)/SIGMA(I)
                  CALL DAXPY(LOCLEN,-TAU(I,J),WRK(IR+I1),1,WRK(IR+I0),1)
   50         CONTINUE

* 19. sigma(j)=r(j)^{T}r(j), gamma'(j)=r(0)^{T}r(j)/sigma(j)
              DOTS(1) = DDOT(LOCLEN,WRK(IR+I0),1,WRK(IR+I0),1)
              DOTS(2) = DDOT(LOCLEN,WRK(IR),1,WRK(IR+I0),1)
              CALL PDSUM(2,DOTS,ipar)
              SIGMA(J) = DOTS(1)

              IF (SIGMA(J).EQ.ZERO) THEN
                  STATUS = -3
                  STEPERR = 19
                  GO TO 9999
              END IF

              GAMMA1(J) = DOTS(2)/SIGMA(J)
              I0 = I0 + LOCLEN
   60     CONTINUE

* 20. gamma_{l}=omega=gamma'_{l}
*     gamma_{j}=gamma'_{j}-\sum_{i=j+1}^{l}{tau_{j,i}gamma_{i}}
          GAMMA(BASIS) = GAMMA1(BASIS)
          OMEGA = GAMMA(BASIS)
          DO 80 J = BASIS - 1,1,-1
              S = ZERO
              DO 70 I = J + 1,BASIS
                  S = S + TAU(J,I)*GAMMA(I)
   70         CONTINUE
              GAMMA(J) = GAMMA1(J) - S
   80     CONTINUE

* 21. gamma''=gamma_{j+1}+\sum_{i=j+1}^{l-1}{tau_{j,i}gamma_{i+1}}
          DO 100 J = 1,BASIS - 1
              S = ZERO
              DO 90 I = J + 1,BASIS - 1
                  S = S + TAU(J,I)*GAMMA(I+1)
   90         CONTINUE
              GAMMA2(J) = GAMMA(J+1) + S
  100     CONTINUE

*  Update

* 22. x(0)=x(0)+gamma(1)r(0)
          CALL DAXPY(LOCLEN,GAMMA(1),WRK(IR),1,X,1)

* 23. r(0)=r(0)-gamma'(l)r(l)
          CALL DAXPY(LOCLEN,-GAMMA1(BASIS),WRK(IR+BASIS*LOCLEN),1,
     +               WRK(IR),1)

* 24. u(0)=u(0)-gamma(l)u(l)
          CALL DAXPY(LOCLEN,-GAMMA(BASIS),WRK(IU+BASIS*LOCLEN),1,
     +               WRK(IU),1)

          I0 = LOCLEN
          DO 110 J = 1,BASIS - 1

* 25. u(0)=u(0)-gamma(j)u(j), j=1:l-1
              CALL DAXPY(LOCLEN,-GAMMA(J),WRK(IU+I0),1,WRK(IU),1)

* 26. x(0)=x(0)+gamma''(j)r(j), j=1:l-1
              CALL DAXPY(LOCLEN,GAMMA2(J),WRK(IR+I0),1,X,1)

* 27. r(0)=r(0)-gamma'(j)r(j), j=1:l-1
              CALL DAXPY(LOCLEN,-GAMMA1(J),WRK(IR+I0),1,WRK(IR),1)
              I0 = I0 + LOCLEN
  110     CONTINUE

  120 CONTINUE

      IF (ITNO.GT.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),X,ipar)
      END IF

*  Set output parameters
      ipar(11) = ITNO
      ipar(12) = STATUS
      ipar(13) = STEPERR
      dpar(2) = EXITNORM

      RETURN
      END
C^L      


C----------------------------------------------------------------------C
c                                                                      c
      REAL(rk) FUNCTION dsetrhsstop(neq,nnz,iq,jq,dq,q,b,r,epsilon,ipar,preconl,pdnrm)
c                                                                      c
C----------------------------------------------------------------------C
c
      USE mCommon
      IMPLICIT NONE
*     .. Scalar Arguments ..
      REAL(rk)  EPSILON
*     ..
*     .. Array Arguments ..
      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq)
      REAL(rk)  q(neq)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(13)
      REAL(rk)  b(neq),r(neq)
*     ..
*     .. Function Arguments ..
      REAL(rk)   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL PRECONL
*     ..
*     .. Local Scalars ..
      INTEGER LOCLEN,STOPT
*     ..
      LOCLEN = ipar(4)
      STOPT = ipar(9)
*     ..
      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.4) .OR. (STOPT.EQ.7)) THEN

*  ||r||<epsilon or ||Q1r||<epsilon ||x(k)-x(k-1)||<epsilon
          DSETRHSSTOP = EPSILON

      ELSE IF ((STOPT.EQ.2) .OR. (STOPT.EQ.3) .OR. (STOPT.EQ.5)) THEN

*  ||r||<epsilon||b|| or sqrt(r(Q1r))<epsilon||b|| or ||Q1r||<epsilon||b||
          DSETRHSSTOP = EPSILON*PDNRM(LOCLEN,B,ipar)

      ELSE IF (STOPT.EQ.6) THEN
*  ||Q1r||<epsilon||Q1b||
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,B,R,ipar)
          DSETRHSSTOP = EPSILON*PDNRM(LOCLEN,R,ipar)
      END IF

      RETURN
      END

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE stopcrit
     &   (neq,nnz,iq,jq,dq,q,ia,ja,a,b,r,rtrue,x,xold,wrk,
     &    rhsstop,cnvrtx,exitnorm,status,ipar,
     &    matvec,tmatvec,preconr,pdsum,pdnrm)
c                                                                      c
C----------------------------------------------------------------------C
*     ..
      USE mCommon
      IMPLICIT NONE
*     .. Arguments ..
      INTEGER  neq,nnz,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz),
     &         cnvrtx,status,ipar(13)
      REAL(rk)   q(nnz),a(nnz),
     &         b(neq),r(neq),rtrue(neq),x(neq),xold(neq),wrk(neq),
     &         rhsstop,exitnorm,pdnrm
      EXTERNAL MATVEC,TMATVEC,PRECONR,PDSUM,PDNRM
*     ..
*     .. Local ..
      INTEGER  LOCLEN,PRECONT,STOPT
      REAL*8   DOTS(1),DDOT
      EXTERNAL DAXPY,DCOPY,DDOT
      INTRINSIC SQRT
*     ..
*     .. Parameters ..
      REAL*8   ZERO, ONE
      PARAMETER (ZERO=0.0_rk, ONE=1.0_rk)
*     ..
      LOCLEN = ipar(4)
      PRECONT = ipar(8)
      STOPT = ipar(9)

      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.2) .OR. (STOPT.EQ.3)) THEN

*  Compute true residual if needed
          CALL DCOPY(LOCLEN,B,1,RTRUE,1)

          IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK,ipar)
              IF (CNVRTX.EQ.1) THEN
*    r=b-AATQ2x
                  CALL TMATVEC(NEQ,NNZ,IA,JA,A,WRK,XOLD,ipar)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,XOLD,WRK,ipar)
                  CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
              ELSE
*    r=b-AQ2x
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK,XOLD,ipar)
                  CALL DAXPY(LOCLEN,-ONE,XOLD,1,RTRUE,1)
              END IF
          ELSE IF (CNVRTX.EQ.1) THEN
*    r=b-AATx
              CALL TMATVEC(NEQ,NNZ,IA,JA,A,X,XOLD,ipar)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,XOLD,WRK,ipar)
              CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
          ELSE
*    r=b-Ax
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK,ipar)
              CALL DAXPY(LOCLEN,-ONE,WRK,1,RTRUE,1)
          END IF
      END IF

      IF ((STOPT.EQ.1) .OR. (STOPT.EQ.2)) THEN

*  ||r||<epsilon or ||r||<epsilon||b||
          EXITNORM = PDNRM(LOCLEN,RTRUE,ipar)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF (STOPT.EQ.3) THEN

*  sqrt(rT(Q1r))<epsilon||b||
          DOTS(1) = DDOT(LOCLEN,RTRUE,1,R,1)
          CALL PDSUM(1,DOTS(1),ipar)
          IF (DOTS(1).LT.ZERO) THEN
              STATUS = -5
              RETURN
          END IF
          EXITNORM = SQRT(DOTS(1))
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF ((STOPT.EQ.4) .OR. (STOPT.EQ.5) .OR. (STOPT.EQ.6)) THEN

*  ||Q1r||<epsilon or ||Q1r||<epsilon||b|| or ||Q1r||<epsilon||Q1b||
          EXITNORM = PDNRM(LOCLEN,R,ipar)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF

      ELSE IF (STOPT.EQ.7) THEN

*  ||x-x0||<epsilon
          CALL DCOPY(LOCLEN,X,1,WRK,1)
          CALL DAXPY(LOCLEN,-ONE,XOLD,1,WRK,1)
          EXITNORM = PDNRM(LOCLEN,WRK,ipar)
          IF (EXITNORM.LT.RHSSTOP) THEN
              STATUS = 0
          ELSE
              STATUS = -99
          END IF
      END IF

      RETURN
      END
C^L      

C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pimRGMRES
     &   (neq,nnz,pret,maxit,stopt,eps,nit,ierr,iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using RGMRES iterative method          **
c ** - --                            ------                           **
c **********************************************************************
      USE mCommon
      IMPLICIT NONE
c Parameters
      INTEGER  neq,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(neq+1),jq(nnz),dq(neq),
     &         ia(neq+1),ja(nnz),da(neq)
      REAL(rk)   q(nnz),a(nnz),b(neq),x(neq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,pdsum,pdnrm,progress
c Internal variables
      INTEGER  nwrk,ipar(13),pres,rst
      REAL(rk)   wrk,dpar(6)
      EXTERNAL dinit,dvprod,drgmres
c
      ALLOCATABLE wrk(:)
c
c Set initial values
c
      rst=10
c if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c
c Allocate memory for matrix and vectors
c
      nwrk = (4+rst)*neq
      ALLOCATE(wrk(nwrk))
      CALL pimdsetpar(ipar,dpar,neq,neq,neq,neq,rst,-1,-1,
     &                pres,stopt,maxit,eps)
c
c RGMRES
c
      CALL drgmres(neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &             matvec,preconl,preconr,pdsum,pdnrm,progress)
c
c Set output values
c
      nit=ipar(11)
      ierr=ipar(12)
c
c Free memory
c
      DEALLOCATE(wrk)
c
      RETURN
      END

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE prenon(neq,nnz,iq,jq,dq,q,u,v,ipar)

      USE mCommon
      IMPLICIT NONE

      INTEGER neq,nnz,iq(neq+1),jq(nnz),dq(neq),ipar(13)
      REAL(rk)  q(neq),u(neq),v(neq)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Precondition system with none preconditioner                    **
c **  ---                      ---                                    **
c **********************************************************************
c
      RETURN
      END


C----------------------------------------------------------------------C
c                                                                      c      
      SUBROUTINE pimCGS
     &   (neq,nnz,pret,maxit,stopt,eps,nit,ierr,iq,jq,dq,q,ia,ja,a,b,x,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c      
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using CGS iterative method             **
c ** - --                            ---                              **
c **********************************************************************
      USE mCommon
      IMPLICIT NONE
c Parameters
      INTEGER  neq,nnz,pret,maxit,stopt,nit,ierr,
     &         iq(neq+1),jq(nnz),dq(neq),
     &         ia(neq+1),ja(nnz),da(neq)
      REAL(rk)   q(nnz),a(nnz),b(neq),x(neq),eps,pdnrm
      EXTERNAL matvec,preconl,preconr,pdsum,pdnrm,progress
c Internal variables
      INTEGER  nwrk,ipar(13),pres
      REAL(rk)   wrk,dpar(6)
      EXTERNAL dinit,dvprod,dcgs
c
      ALLOCATABLE wrk(:)
c
c Set initial values
c
c if pret>0 use left preconditioning
      pres=0
      IF (pret.GT.0) pres=1
c
c Allocate memory for matrix and vectors
c
      nwrk = (9)*neq
      ALLOCATE(wrk(nwrk))
      CALL pimdsetpar(ipar,dpar,neq,neq,neq,neq,-1,-1,-1,
     &                pres,stopt,maxit,eps)
c
c CGS
c
      CALL dcgs(neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &         matvec,preconl,preconr,pdsum,pdnrm,progress)
c
c Set output values
c
      nit=ipar(11)
      ierr=ipar(12)
c
c Free memory
c
      DEALLOCATE(wrk)
c
      RETURN
      END      

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE dcgs
     &   (neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &    matvec,preconl,preconr,pdsum,pdnrm,progress)
c                                                                      c
C----------------------------------------------------------------------C
*     ..
      USE mCommon
      IMPLICIT NONE
*     .. Parameters ..
      REAL(rk)  ZERO
      PARAMETER (ZERO=0.0_rk)
      REAL(rk)  ONE
      PARAMETER (ONE=1.0_rk)
      INTEGER IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER DPARSIZ
      PARAMETER (DPARSIZ=6)
*     ..
*     .. Matrices A,Q ..
      INTEGER neq,nnz,nwrk,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz)
      REAL*8  q(nnz),a(nnz)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(IPARSIZ)
      REAL*8  b(neq),x(neq),dpar(DPARSIZ),wrk(nwrk)
*     ..
*     .. Function Arguments ..
      REAL*8   PDNRM
      EXTERNAL PDNRM
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PDSUM,PRECONL,PRECONR,PROGRESS
*     ..
*     .. Local Scalars ..
      INTEGER BASIS,BLKSZ,CNVRTX,IP,IR,IRTILDE,IS,IT,ITNO,IU,IW,
     +        IXOLD,IZ,LDA,LOCLEN,MAXIT,N,NPROCS,PRECONT,PROCID,
     +        STATUS,STEPERR,STOPT
      REAL*8  ALPHA,BETA,EPSILON,EXITNORM,RHO,RHO0,RHSSTOP,XI
*     ..
*     .. Local Arrays ..
      REAL*8  DOTS(1)
*     ..
*     .. External Functions ..
      REAL*8   DDOT,DSETRHSSTOP
      EXTERNAL DDOT,DSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,PIMDGETPAR,STOPCRIT
*     ..
      CALL PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASIS,NPROCS,
     +                PROCID,PRECONT,STOPT,MAXIT,ITNO,STATUS,
     +                STEPERR,EPSILON,EXITNORM)

*  Check consistency of preconditioning and stop types
      IF (((PRECONT.EQ.0).OR. (PRECONT.EQ.2)) .AND. (STOPT.EQ.6)) THEN
          ITNO = 0
          STATUS = -4
          STEPERR = 0
          GO TO 9999
      END IF

*  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

*  Set indices for mapping local vectors into wrk
      IR = 1
      IRTILDE = IR + LOCLEN
      IP = IRTILDE + LOCLEN
      IS = IP + LOCLEN
      IT = IS + LOCLEN
      IU = IT + LOCLEN
      IW = IU + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN

*  Set rhs of stopping criteria
      RHSSTOP = DSETRHSSTOP(NEQ,NNZ,IQ,JQ,DQ,Q,B,WRK(IR),EPSILON,IPAR,
     +                      PRECONL,PDNRM)

*  1. r=Q1(b-AQ2x)
      IF (STOPT.NE.6) THEN
          IF (PRECONT.EQ.0) THEN
*     r=b-Ax
              CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IR),IPAR)

          ELSE IF (PRECONT.EQ.2) THEN
*     r=b-AQ2x
              CALL DCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL DCOPY(LOCLEN,B,1,WRK(IP),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IP),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IR),IPAR)
          END IF

      ELSE
*     r has been set to Qb in the call to dsetrhsstop
          IF (PRECONT.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)

          ELSE IF (PRECONT.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IZ),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IR),1)
          END IF

      END IF

*  2. p=s=rtilde=r
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
      CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)

*  3. rho=dot(rtilde,r)
      DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
      CALL PDSUM(1,DOTS,IPAR)
      RHO = DOTS(1)

*  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO 20 ITNO = 1,MAXIT

*  4. w=Q1AQ2p
          IF (PRECONT.EQ.0) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IP),WRK(IW),IPAR)

          ELSE IF (PRECONT.EQ.1) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IP),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IW),IPAR)

          ELSE IF (PRECONT.EQ.2) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IZ),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),IPAR)

          ELSE IF (PRECONT.EQ.3) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IP),WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IW),IPAR)
          END IF

*  5. xi=dot(rtilde,w)
          DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IW),1)
          CALL PDSUM(1,DOTS,IPAR)
          XI = DOTS(1)

*  6. alpha=rho/xi
          IF (XI.EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 6
              GO TO 9999

          END IF

          ALPHA = RHO/XI

*  7. t=s-alpha*w
          CALL DCOPY(LOCLEN,WRK(IS),1,WRK(IT),1)
          CALL DAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IT),1)

*  8. w=s+t
          CALL DCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
          CALL DAXPY(LOCLEN,ONE,WRK(IT),1,WRK(IW),1)

*  9. x=x+alpha*w
          CALL DCOPY(LOCLEN,X,1,WRK(IXOLD),1)
          CALL DAXPY(LOCLEN,ALPHA,WRK(IW),1,X,1)

* 10. r=r-alpha*Q1AQ2w
          IF (PRECONT.EQ.0) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IU),IPAR)

          ELSE IF (PRECONT.EQ.1) THEN
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IU),IPAR)

          ELSE IF (PRECONT.EQ.2) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IU),IPAR)

          ELSE IF (PRECONT.EQ.3) THEN
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IU),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IU),WRK(IZ),IPAR)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IU),IPAR)
          END IF

          CALL DAXPY(LOCLEN,-ALPHA,WRK(IU),1,WRK(IR),1)

* 11. check stopping criterion
          CALL STOPCRIT(NEQ,NNZ,IQ,JQ,DQ,Q,IA,JA,A,B,
     +                  WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IU),RHSSTOP,
     +                  CNVRTX,EXITNORM,STATUS,IPAR,MATVEC,MATVEC,
     +                  PRECONR,PDSUM,PDNRM)

*  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

          IF (STATUS.EQ.-5) THEN
              STEPERR = 11
              GO TO 9999
          ELSE IF (STATUS.EQ.0) THEN
              GO TO 9999
          END IF

* 12. rho=dot(rtilde0,r)
          RHO0 = RHO
          DOTS(1) = DDOT(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
          CALL PDSUM(1,DOTS,IPAR)
          RHO = DOTS(1)

* 13. beta=rho/rho0
          IF (RHO0.EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 13
              GO TO 9999

          END IF

          BETA = RHO/RHO0

* 14. s=r+beta*t
          CALL DCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)
          CALL DAXPY(LOCLEN,BETA,WRK(IT),1,WRK(IS),1)

* 15. w=t+beta*p
          CALL DCOPY(LOCLEN,WRK(IT),1,WRK(IW),1)
          CALL DAXPY(LOCLEN,BETA,WRK(IP),1,WRK(IW),1)

* 16. p=s+beta*w
          CALL DCOPY(LOCLEN,WRK(IS),1,WRK(IP),1)
          CALL DAXPY(LOCLEN,BETA,WRK(IW),1,WRK(IP),1)

   20 CONTINUE

      IF (ITNO.GT.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONT.EQ.2) .OR. (PRECONT.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),X,IPAR)
      END IF

*  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      DPAR(2) = EXITNORM

      RETURN

      END
C^L      


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE drgmres
     &   (neq,nnz,nwrk,iq,jq,dq,q,ia,ja,a,b,x,wrk,ipar,dpar,
     &    matvec,preconl,preconr,pdsum,pdnrm2,progress)
c                                                                      c
C----------------------------------------------------------------------C
      USE mCommon
      IMPLICIT NONE
*     ..
*     .. Parameters ..
      REAL(rk)  ZERO
      PARAMETER (ZERO=0.0_rk)
      REAL(rk)  ONE
      PARAMETER (ONE=1.0_rk)
      INTEGER IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER DPARSIZ
      !PARAMETER (DPARSIZ=2) ! original
      PARAMETER (DPARSIZ=6)
      INTEGER IBDIM
      PARAMETER (IBDIM=50)
      INTEGER LDR
      PARAMETER (LDR=IBDIM+1)
      INTEGER LDG
      PARAMETER (LDG=IBDIM+2)
*     ..
*     .. Matrices A,Q ..
      INTEGER neq,nnz,nwrk,iq(neq+1),jq(nnz),dq(neq),ia(neq+1),ja(nnz)
      REAL(rk)  q(nnz),a(nnz)
*     ..
*     .. Array Arguments ..
      INTEGER ipar(IPARSIZ)
      REAL(rk)  b(neq),x(neq),dpar(DPARSIZ),wrk(nwrk)
*     ..
*     .. Function Arguments ..
      REAL(rk)   PDNRM2
      EXTERNAL PDNRM2
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL MATVEC,PDSUM,PRECONL,PRECONR,PROGRESS
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION BETA,EPSILON,ETA,EXITNORM,KSI,RHSSTOP,TAU1,TAU2
      INTEGER BASISDIM,BLKSZ,I,IRES,ITNO,IV,IW,IZ,J,K0,K1,LDA,LOCLEN,
     +        MAXIT,N,NPROCS,PRECONTYPE,PROCID,STATUS,STEPERR,STOPTYPE
      LOGICAL ENDED
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION G(LDG),R(LDR,LDR),RHO(IBDIM)
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DDOT,DSETRHSSTOP
      EXTERNAL DDOT,DSETRHSSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DECODE,DINIT,DSCAL,DTRSV,ENCODE,GIVENS,
     +         PIMDGETPAR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
      CALL PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,
     +                PROCID,PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,
     +                STEPERR,EPSILON,EXITNORM)

*  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE.EQ.0).OR. (PRECONTYPE.EQ.2)) .AND.
     +    (STOPTYPE.EQ.6)) THEN
          ITNO = 0
          STATUS = -4
          STEPERR = 0
          GO TO 9999

      END IF

*  Set indices for mapping local vectors into wrk
      IRES = 1
      IZ = IRES + LOCLEN
      IW = IZ + LOCLEN
      IV = IW + LOCLEN

*  Set rhs of stopping criteria
      RHSSTOP = DSETRHSSTOP(NEQ,NNZ,IQ,JQ,DQ,Q,B,WRK(IRES),EPSILON,IPAR,
     +                      PRECONL,PDNRM2)

*  1. r=Q1(b-AQ2x)
      IF (PRECONTYPE.EQ.0) THEN
*     r=b-Ax
          CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IRES),1)

      ELSE IF (PRECONTYPE.EQ.1) THEN
*     r=Q1(b-Ax)
          CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IRES),IPAR)

      ELSE IF (PRECONTYPE.EQ.2) THEN
*     r=b-AQ2x
          CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IRES),1)

      ELSE IF (PRECONTYPE.EQ.3) THEN
*     r=Q1(b-AQ2x)
          CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IRES),IPAR)
          CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IRES),WRK(IZ),IPAR)
          CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
          CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IRES),IPAR)
      END IF

*  2. beta=||r||_2
      BETA = PDNRM2(LOCLEN,WRK(IRES),IPAR)

*  Loop
      STATUS = 0
      STEPERR = -1
      EXITNORM = -ONE
      ENDED = .FALSE.
      DO 20 ITNO = 1,MAXIT

*  3. g=(beta,beta,...)
          G(1) = BETA
          G(2) = BETA

*  4. V(1)=r/beta
          IF (BETA.EQ.ZERO) THEN
              STATUS = -3
              STEPERR = 4
              GO TO 9999

          END IF

          CALL DCOPY(LOCLEN,WRK(IRES),1,WRK(IV),1)
          CALL DSCAL(LOCLEN,ONE/BETA,WRK(IV),1)

          K0 = 0
          DO 40 J = 1,BASISDIM

*     z=Q1AQ2V(j)
              IF (PRECONTYPE.EQ.0) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IV+K0),WRK(IZ),IPAR)

              ELSE IF (PRECONTYPE.EQ.1) THEN
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IV+K0),WRK(IW),IPAR)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)

              ELSE IF (PRECONTYPE.EQ.2) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IV+K0),WRK(IW),IPAR)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)

              ELSE IF (PRECONTYPE.EQ.3) THEN
                  CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IV+K0),WRK(IZ),IPAR)
                  CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IZ),WRK(IW),IPAR)
                  CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IZ),IPAR)
              END IF


*  5. R(i,j)=dot(V(i),Q1AQ2V(j))
              K1 = 0
              DO 50 I = 1,J
                  R(I,J) = DDOT(LOCLEN,WRK(IV+K1),1,WRK(IZ),1)
                  K1 = K1 + LOCLEN
   50         CONTINUE
              CALL PDSUM(J,R(1,J),IPAR)

*  6. Vhat(j)=Q1AQ2V(j)-sum_{i=1}^{j}{R(i,j)V(i)}
              K1 = 0
              CALL DINIT(LOCLEN,ZERO,WRK(IW),1)
              DO 60 I = 1,J
                  CALL DAXPY(LOCLEN,R(I,J),WRK(IV+K1),1,WRK(IW),1)
                  K1 = K1 + LOCLEN
   60         CONTINUE
              CALL DSCAL(LOCLEN,-ONE,WRK(IW),1)
              CALL DAXPY(LOCLEN,ONE,WRK(IZ),1,WRK(IW),1)

*  From this point, w holds the (j+1)-st column of vhat

*  7. R(j+1,j)=||Vhat(j)||_2
              R(J+1,J) = PDNRM2(LOCLEN,WRK(IW),IPAR)

*  8. V(j+1)=Vhat(j)/R(j+1,j)
              IF (R(J+1,J).EQ.ZERO) THEN
                  STATUS = -2
                  STEPERR = 8
                  GO TO 9999

              END IF

              K0 = K0 + LOCLEN
              CALL DSCAL(LOCLEN,ONE/R(J+1,J),WRK(IW),1)
              CALL DCOPY(LOCLEN,WRK(IW),1,WRK(IV+K0),1)

*  9. Apply previous Givens' rotations to column j of R
              DO 70 I = 1,J - 1
                  CALL DECODE(RHO(I),KSI,ETA)
                  TAU1 = R(I,J)
                  TAU2 = R(I+1,J)
                  R(I,J) = KSI*TAU1 - ETA*TAU2
                  R(I+1,J) = ETA*TAU1 + KSI*TAU2
   70         CONTINUE

* 10. Compute Givens' rotation to zero element R(j+1,j)
              CALL GIVENS(R(J,J),R(J+1,J),KSI,ETA)
              TAU1 = R(J,J)
              TAU2 = R(J+1,J)
              R(J,J) = KSI*TAU1 - ETA*TAU2
              R(J+1,J) = ETA*TAU1 + KSI*TAU2
              CALL ENCODE(RHO(J),KSI,ETA)

*  11. Update g
              G(J) = G(J)*KSI
              G(J+1) = G(J+1)*ETA
              G(J+2) = G(J+1)

*  12. If |g(j+1)|<rhsstop stop
              EXITNORM = ABS(G(J+1))
              IF (EXITNORM.LT.RHSSTOP) THEN
                  BASISDIM = J
                  ENDED = .TRUE.
                  GO TO 80

              END IF

   40     CONTINUE
   80     CONTINUE

*  13. Solve Ry=g (solution to least-squares problem)
          CALL DTRSV('U','N','N',BASISDIM,R,LDR,G,1)

*  14. x=x+Vy (Form approximate solution after a c-cycle)
          K1 = 0
          DO 100 I = 1,BASISDIM
              CALL DAXPY(LOCLEN,G(I),WRK(IV+K1),1,X,1)
              K1 = K1 + LOCLEN
  100     CONTINUE

*  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IRES),WRK(IRES))

          IF (ENDED) GO TO 9999

*  15. r=Q1(b-AQ2x)
          IF (PRECONTYPE.EQ.0) THEN
*     r=b-Ax
              CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IRES),1)

          ELSE IF (PRECONTYPE.EQ.1) THEN
*     r=Q1(b-Ax)
              CALL DCOPY(LOCLEN,B,1,WRK(IZ),1)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,X,WRK(IW),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IW),1,WRK(IZ),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),WRK(IRES),IPAR)

          ELSE IF (PRECONTYPE.EQ.2) THEN
*     r=b-AQ2x
              CALL DCOPY(LOCLEN,B,1,WRK(IRES),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IW),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IW),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IRES),1)

          ELSE IF (PRECONTYPE.EQ.3) THEN
*     r=Q1(b-AQ2x)
              CALL DCOPY(LOCLEN,B,1,WRK(IW),1)
              CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,X,WRK(IRES),IPAR)
              CALL MATVEC(NEQ,NNZ,IA,JA,A,WRK(IRES),WRK(IZ),IPAR)
              CALL DAXPY(LOCLEN,-ONE,WRK(IZ),1,WRK(IW),1)
              CALL PRECONL(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IW),WRK(IRES),IPAR)
          END IF

*  16. beta=||r||_2
          BETA = PDNRM2(LOCLEN,WRK(IRES),IPAR)

   20 CONTINUE

      IF (ITNO.GT.MAXIT) THEN
          STATUS = -1
          ITNO = MAXIT
      END IF

 9999 CONTINUE

      IF ((PRECONTYPE.EQ.2) .OR. (PRECONTYPE.EQ.3)) THEN
          CALL DCOPY(LOCLEN,X,1,WRK(IZ),1)
          CALL PRECONR(NEQ,NNZ,IQ,JQ,DQ,Q,WRK(IZ),X,IPAR)
      END IF

*  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      DPAR(2) = EXITNORM

      RETURN

      END

      SUBROUTINE ENCODE(RHO,C,S)
      USE mCommon
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,RHO,S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0_rk)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0_rk)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0_rk)
*     ..
      IF (C.EQ.ZERO) THEN
          RHO = ONE

      ELSE IF (ABS(S).LT.ABS(C)) THEN
          RHO = SIGN(ONE,C)*S/TWO

      ELSE
          RHO = TWO*SIGN(ONE,S)/C
      END IF

      RETURN

      END

      SUBROUTINE DECODE(RHO,C,S)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,RHO,S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
*     ..
      IF (RHO.EQ.ONE) THEN
          C = ZERO
          S = ONE

      ELSE IF (ABS(RHO).LT.ONE) THEN
          S = TWO*RHO
          C = SQRT(ONE-S**2)

      ELSE
          C = TWO/RHO
          S = SQRT(ONE-C**2)
      END IF

      RETURN

      END

      SUBROUTINE GIVENS(A,B,C,S)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,S
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TAU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
*     ..
      IF (B.EQ.ZERO) THEN
          C = ONE
          S = ZERO

      ELSE IF (ABS(B).GT.ABS(A)) THEN
          TAU = -A/B
          S = ONE/SQRT(ONE+TAU**2)
          C = S*TAU

      ELSE
          TAU = -B/A
          C = ONE/SQRT(ONE+TAU**2)
          S = C*TAU
      END IF

      RETURN

      END
C^L