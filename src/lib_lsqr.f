C
C     Least sqares overdetermined system solver
C
C     Usage :
C
C     preconditioner (calculate : Pivot )
C
c      CALL lsqr_precon(neq,nx,nnz,Pivot,crs%i,crs%j,crs%v)

C
C     solve overdetermined system of linear equations (calculate : x )
C
c      CALL lsqr_solve(neq,nx,nnz,maxit,eps,nits,ierr,Pivot,crs%i,crs%j,crs%v,b,x)

C
C     Sizes
C
C      - x(nx)
C      - b(neq)
C      - Pivot(neq)
C      - crs%i(neq+1)
C      - crs%j(nnz)
C      - crs%v(nnz)
C
C     Types
C
C      - INTEGER(8) nnz,crs%i
C      - INTEGER neq,nx,crs%j,maxit,nits,ierr
C      - REAL(rk) x,b,Pivot,crs%v,eps
C
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_solve(m,n,nnz,maxit,eps,nit,ierr,q,ia,ja,a,b,x)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using Least Square iterative method    **
c ** - --                            -     --  -                      **
c ** !!! With built-in diagonal preconditioning !!!                   **
c ** a call to prelsqr2 must preceede a call to slvlsqr2              **
c **********************************************************************
      USE mCommon
      IMPLICIT NONE
c Parameters
      INTEGER(8) nnz
      INTEGER  m,n,maxit,nit,ierr,ja(nnz)
      INTEGER(8) ia(m+1)
      REAL(rk)   eps,q(n),a(nnz),b(m),x(n)
c Internal variables
      LOGICAL  wantse,extra
      INTEGER  istop,iout,i
      REAL(rk)   x0,wm,u,se,v,w
      REAL(rk)   atol,btol,nrm_b,nrm_r
      REAL(rk)   zero,damp,conlim,anorm,acond,rnorm,arnorm,xnorm
      REAL(rk)   dnrm2
c
      ALLOCATABLE x0(:),wm(:),u(:),se(:),v(:),w(:)
c
c Allocate memory for matrix and vectors
c
      ALLOCATE(x0(n),wm(m),u(m),se(1),v(n),w(n))
c
c Set initial values
c
      wantse=.FALSE.
      extra=.FALSE.
c
      zero=0.0
      damp=zero
      conlim=1.0D8
c
      IF (maxit.EQ.0) THEN
        maxit = 10*(m + n + 50)
        damp=zero
      END IF
c
      IF (extra) THEN
        iout=66
      ELSE
        iout=-1
      END IF
c
      CALL dcopy(n,x,1,x0,1)                       !x0=x
      wm=zero
      DO i=1,n
        x(i)=x(i)/q(i)                             !x=Q^-1*x
      END DO
      CALL lsqr_aprod(1,m,n,nnz,x,wm,iA,jA,A)           !wm=A*x
      CALL dcopy(m,b,1,u,1)                        !
      CALL daxpy(m,-1.0_rk,wm,1,u,1)                !u=u-wm=b-Ax
      nrm_r=dnrm2(m,u,1)
      nrm_b=dnrm2(m,b,1)
      atol=eps
      btol=eps*nrm_b/nrm_r
c      print *,btol,nrm_b,nrm_r
c
c LSQR
c
      CALL lsqr_LSQR(
     &  m,n,nnz,damp,wantse,
     &  iA,jA,A,
     &  u,x,se,v,w,
     &  atol,btol,conlim,maxit,
     &  istop,nit,anorm,acond,rnorm,arnorm,xnorm)
c
c Preconditioninig
c
      CALL lsqr_dvprod(n,Q,1,x,1)                         !x=Q*x
c
c Hotstart
c
      CALL daxpy(n,1.0_rk,x0,1,x,1)                   !x=x0+dx
c
      IF (istop.GT.3) THEN
c       WRITE (*,*) 'Warning: LSQR: stoping condition ISTOP= ',istop
        ierr=istop
      ELSE
        ierr=0
      END IF
c
c Free memory
c
      DEALLOCATE(x0,wm,u,se,v,w)
c
      RETURN
      END

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_precon(m,n,nnz,q,ia,ja,a)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c ** Solve system of equations using Least Square iterative method    **
c ** - --                            -     --  -                      **
c ** !!! With built-in diagonal preconditioning !!!                   **
c **********************************************************************
      USE mCommon
      IMPLICIT NONE
c Parameters
      INTEGER(8) nnz
      INTEGER m,n
      INTEGER(8) ia(m+1)
      INTEGER ja(nnz)
      REAL(rk)   q(n),a(nnz)
c Internal variables
      REAL(rk)   wm
c
      ALLOCATABLE wm(:)
c
c Allocate memory for matrix and vectors
c
      ALLOCATE(wm(m))
c
c Compute preconditioner
c
      CALL lsqr_aprod(3,m,n,nnz,Q,wm,iA,jA,A)            !Q=M^-1=1/diag(sqrt(AT*A))
      CALL lsqr_aprod(4,m,n,nnz,Q,wm,iA,jA,A)            !A=A*M^-1
c
c Free memory
c
      DEALLOCATE(wm)
c
      RETURN
      END


C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_aprod(MODE,M,N,NNZ,X,Y,IA,JA,A )
      
      USE mCommon            
      IMPLICIT NONE

      INTEGER(8) nnz
      INTEGER mode,m,n,ja(nnz)
      INTEGER(8) ia(m+1)
      REAL(rk)  x(n),y(m),a(nnz)
c                                                                      c
C----------------------------------------------------------------------C
c                If mode = 1, compute  y = y + A*x.
c                If mode = 2, compute  x = x + AT*y.
c                If mode = 3, compute  x = 1/diag(sqrt((AT*A)).
c                If mode = 4, compute  A_ij = A_ij*x_j.
C----------------------------------------------------------------------C
      IF (mode.EQ.1) THEN
        CALL lsqr_matvec1(m,n,nnz,ia,ja,a,x,y)
      ELSE IF (mode.EQ.2) THEN
        CALL lsqr_tmatvec2(m,n,nnz,ia,ja,a,y,x)
      ELSE IF (mode.EQ.3) THEN
        CALL lsqr_tmatvec3(m,n,nnz,ia,ja,a,x)
      ELSE IF (mode.EQ.4) THEN
        CALL lsqr_matvec4(m,n,nnz,ia,ja,a,x)
      END IF
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_matvec1(m,n,nnz,ia,ja,a,u,v)
      
      USE mCommon
      IMPLICIT NONE

      INTEGER(8) nnz
      INTEGER m,n,ja(nnz)
      INTEGER(8) ia(m+1)
      REAL(rk)  a(nnz),u(n),v(m)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Matrix vector product: v = v + A*u                              **
c **                                                                  **
c **********************************************************************
      INTEGER i,ii
      INTEGER(8) j
      REAL(rk)  sum
c
      DO i = 1,m
        sum=0.0_rk
        DO j = ia(i),ia(i+1) - 1
          ii=ja(j)
          sum=sum+a(j)*u(ii)
        END DO
        v(i)=v(i)+sum
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_tmatvec2(m,n,nnz,ia,ja,a,u,v)
      USE mCommon
      IMPLICIT NONE

      INTEGER(8) nnz
      INTEGER m,n,ja(nnz)
      INTEGER(8) ia(m+1)
      REAL(rk)  a(nnz),u(m),v(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Transpose matrix vector product: v = v + AT*u                   **
c **                                                                  **
c **********************************************************************
      INTEGER ii,j
      INTEGER(8) i
c
      DO j = 1,m
        DO i = ia(j),ia(j+1) - 1
          ii=ja(i)
          v(ii) = v(ii)+a(i)*u(j)
        END DO
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_tmatvec3(m,n,nnz,ia,ja,a,v)
      USE mCommon
      IMPLICIT NONE

      INTEGER(8) nnz
      INTEGER m,n,ja(nnz)
      INTEGER(8) ia(m+1)
      REAL(rk)  a(nnz),v(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Create diagonal preconditioner: v = 1 / diag(sqrt(AT * A))      **
c **                                                                  **
c **********************************************************************
      INTEGER ii,j
      INTEGER(8) i
c
      CALL lsqr_dinit(n,0.0_rk,v,1)
c
      DO j = 1,m
        DO i = ia(j),ia(j+1) - 1
          ii=ja(i)
          v(ii) = v(ii)+a(i)**2
        END DO
      END DO
c
      DO i=1,n
        v(i)=1.0_rk/sqrt(v(i))
      END DO
c
      RETURN
      END
c
C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_matvec4(m,n,nnz,ia,ja,a,v)
      USE mCommon
      IMPLICIT NONE

      INTEGER(8) nnz
      INTEGER m,n,ja(nnz)
      INTEGER(8) ia(m+1)
      REAL(rk)  a(nnz),v(n)
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c **  Apply diagonal preconditioner: A_ij = A_ij*v_j                  **
c **                                                                  **
c **********************************************************************
      INTEGER i,j
      INTEGER(8) ij
c
      DO i = 1,m
        DO ij = ia(i),ia(i+1) - 1
          j=ja(ij)
          a(ij)=a(ij)*v(j)
        END DO
      END DO
c
      RETURN
      END

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_dinit(n,alpha,dx,incx)
      USE mCommon            
      IMPLICIT NONE
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
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = ALPHA
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          DX(I) = ALPHA
          DX(I+1) = ALPHA
          DX(I+2) = ALPHA
          DX(I+3) = ALPHA
          DX(I+4) = ALPHA
          DX(I+5) = ALPHA
          DX(I+6) = ALPHA
   50 CONTINUE
      RETURN
      END
C^L

C----------------------------------------------------------------------C
c                                                                      c
      SUBROUTINE lsqr_dvprod(n,dx,incx,dy,incy)
      USE mCommon            
      IMPLICIT NONE
c                                                                      c
C----------------------------------------------------------------------C
*
*     Modified from daxpy level 1 BLAS
*     element-wise vector multiplication, y<-x*y
*     Rudnei Dias da Cunha, 16/6/93
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL(rk)  DX(n),DY(n)
*     ..
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY)*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I)*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I)*DX(I)
          DY(I+1) = DY(I+1)*DX(I+1)
          DY(I+2) = DY(I+2)*DX(I+2)
          DY(I+3) = DY(I+3)*DX(I+3)
   50 CONTINUE
      RETURN
      END


      SUBROUTINE lsqr_LSQR  (
     &  m, n, nnz, damp, wantse,
     &  ia, ja, a,
     &  u, x, se, v, w,
     &  atol, btol, conlim, itnlim,
     &  istop, itn, anorm, acond, rnorm, arnorm, xnorm)

      USE mCommon            
      IMPLICIT NONE

      LOGICAL            wantse
      INTEGER(8)         nnz
      INTEGER            m, n, itnlim, nout, istop, itn
      INTEGER(8)         ia(m+1)
      INTEGER            ja(nnz)
      REAL(rk)   a(nnz), u(m), x(n), se(*), v(n), w(n),
     &                   atol, btol, conlim, damp,
     &                   anorm, acond, rnorm, arnorm, xnorm
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c
c     LSQR  finds a solution x to the following problems:
c
c     1. Unsymmetric equations --    solve  A*x = b
c
c     2. Linear least squares  --    solve  A*x = b
c                                    in the least-squares sense
c
c     3. Damped least squares  --    solve  (   A    )*x = ( b )
c                                           ( damp*I )     ( 0 )
c                                    in the least-squares sense
c
c     where A is a matrix with m rows and n columns, b is an
c     m-vector, and damp is a scalar.  (All quantities are REAL(rk).)
c     The matrix A is intended to be large and sparse.  It is accessed
c     by means of SUBROUTINE CALLs of the form
c
c                CALL aprod ( mode, m, n, nnz, x, y, ia, ja, a )
c
c     which must perform the following functions:
c
c                If mode = 1, compute  y = y + A*x.
c                If mode = 2, compute  x = x + A(transpose)*y.
c
c     The vectors x and y are input parameters in both cases.
c     If  mode = 1,  y should be altered without changing x.
c     If  mode = 2,  x should be altered without changing y.
c     The parameters, ia, ja, a  may be used for workspace
c     as described below.
c
c     The rhs vector b is input via u, and subsequently overwritten.
c
c
c     Note:  LSQR uses an iterative method to approximate the solution.
c     The number of iterations required to reach a certain accuracy
c     depends strongly on the scaling of the problem.  Poor scaling of
c     the rows or columns of A should therefore be avoided where
c     possible.
c
c     For example, in problem 1 the solution is unaltered by
c     row-scaling.  If a row of A is very small or large compared to
c     the other rows of A, the corresponding row of ( A  b ) should be
c     scaled up or down.
c
c     In problems 1 and 2, the solution x is easily recovered
c     following column-scaling.  Unless better information is known,
c     the nonzero columns of A should be scaled so that they all have
c     the same Euclidean norm (e.g., 1.0).
c
c     In problem 3, there is no freedom to re-scale if damp is
c     nonzero.  However, the value of damp should be assigned only
c     after attention has been paid to the scaling of A.
c
c     The parameter damp is intended to help regularize
c     ill-conditioned systems, by preventing the true solution from
c     being very large.  Another aid to regularization is provided by
c     the parameter acond, which may be used to terminate iterations
c     before the computed solution becomes very large.
c
c     Note that x is not an input parameter.
c     If some initial estimate x0 is known and if damp = 0,
c     one could proceed as follows:
c
c       1. Compute a residual vector     r0 = b - A*x0.
c       2. Use LSQR to solve the system  A*dx = r0.
c       3. Add the correction dx to obtain a final solution x = x0 + dx.
c
c     This requires that x0 be available before and after the CALL
c     to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
c     to solve A*x = b and k2 iterations to solve A*dx = r0.
c     If x0 is "good", norm(r0) will be smaller than norm(b).
c     If the same stopping tolerances atol and btol are used for each
c     system, k1 and k2 will be similar, but the final solution x0 + dx
c     should be more accurate.  The only way to reduce the total work
c     is to use a larger stopping tolerance for the second system.
c     If some value btol is suitable for A*x = b, the larger value
c     btol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
c
c     Preconditioning is another way to reduce the number of iterations.
c     If it is possible to solve a related system M*x = b efficiently,
c     where M approximates A in some helpful way
c     (e.g. M - A has low rank or its elements are small relative to
c     those of A), LSQR may converge more rapidly on the system
c           A*M(inverse)*z = b,
c     after which x can be recovered by solving M*x = z.
c
c     NOTE: If A is symmetric, LSQR should not be used!
c     Alternatives are the symmetric conjugate-gradient method (cg)
c     and/or SYMMLQ.
c     SYMMLQ is an implementation of symmetric cg that applies to
c     any symmetric A and will converge more rapidly than LSQR.
c     If A is positive definite, there are other implementations of
c     symmetric cg that require slightly less work per iteration
c     than SYMMLQ (but will take the same number of iterations).
c
c
c     Notation
c     --------
c
c     The following quantities are used in discussing the SUBROUTINE
c     parameters:
c
c     Abar   =  (   A    ),          bbar  =  ( b )
c               ( damp*I )                    ( 0 )
c
c     r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
c
c     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
c            =  norm( rbar )
c
c     relpr  =  the relative precision of floating-point arithmetic
c               on the machine being used.  On most machines,
c               relpr is about 1.0e-7 and 1.0D-16 in single and double
c               precision respectively.
c
c     LSQR  minimizes the function rnorm with respect to x.
c
c
c     Parameters
c     ----------
c
c     m       input      m, the number of rows in A.
c
c     n       input      n, the number of columns in A.
c
c     damp    input      The damping parameter for problem 3 above.
c                        (damp should be 0.0 for problems 1 and 2.)
c                        If the system A*x = b is incompatible, values
c                        of damp in the range 0 to sqrt(relpr)*norm(A)
c                        will probably have a negligible effect.
c                        Larger values of damp will tend to decrease
c                        the norm of x and reduce the number of
c                        iterations required by LSQR.
c
c                        The work per iteration and the storage needed
c                        by LSQR are the same for all values of damp.
c
c     wantse  input      A logical variable to say if the array se(*)
c                        of standard error estimates should be computed.
c                        If m .GT. n  or  damp .GT. 0,  the system is
c                        overdetermined and the standard errors may be
c                        useful.  (See the first LSQR reference.)
c                        Otherwise (m .LE. n  and  damp = 0) they DO not
c                        mean much.  Some time and storage can be saved
c                        by setting wantse = .FALSE. and using any
c                        convenient array for se(*), which won't be
c                        touched.
c
c     nnz     input      The length of the workspace array rw.
c     ia      workspace  An INTEGER array of length m+1.
c     ja      workspace  An INTEGER array of length nnz.
c     a       workspace  A REAL(rk) array of length nnz.
c
c             Note:  LSQR  does not explicitly use the previous four
c             PARAMETERs, but passes them to SUBROUTINE aprod for
c             possible use as workspace.
c
c     u(m)    input      The rhs vector b.  Beware that u is
c                        over-written by LSQR.
c
c     v(n)    workspace
c
c     w(n)    workspace
c
c     x(n)    output     Returns the computed solution x.
c
c     se(*)   output     If wantse is true, the DIMENSION of se must be
c             (maybe)    n or more.  se(*) THEN RETURNs standard error
c                        estimates for the components of x.
c                        For each i, se(i) is set to the value
c                           rnorm * sqrt( sigma(i,i) / t ),
c                        where sigma(i,i) is an estimate of the i-th
c                        diagonal of the inverse of Abar(transpose)*Abar
c                        and  t = 1      if  m .LE. n,
c                             t = m - n  if  m .GT. n  and  damp = 0,
c                             t = m      if  damp .NE. 0.
c
c                        If wantse is false, se(*) will not be touched.
c                        The actual parameter can be any suitable array
c                        of any length.
c
c     atol    input      An estimate of the relative error in the data
c                        defining the matrix A.  For example,
c                        if A is accurate to about 6 digits, set
c                        atol = 1.0e-6 .
c
c     btol    input      An estimate of the relative error in the data
c                        defining the rhs vector b.  For example,
c                        if b is accurate to about 6 digits, set
c                        btol = 1.0e-6 .
c
c     conlim  input      An upper limit on cond(Abar), the apparent
c                        condition number of the matrix Abar.
c                        Iterations will be terminated if a computed
c                        estimate of cond(Abar) exceeds conlim.
c                        This is intended to prevent certain small or
c                        zero singular values of A or Abar from
c                        coming into effect and causing unwanted growth
c                        in the computed solution.
c
c                        conlim and damp may be used separately or
c                        together to regularize ill-conditioned systems.
c
c                        Normally, conlim should be in the range
c                        1000 to 1/relpr.
c                        Suggested value:
c                        conlim = 1/(100*relpr)  for compatible systems,
c                        conlim = 1/(10*sqrt(relpr)) for least squares.
c
c             Note:  If the user is not concerned about the parameters
c             atol, btol and conlim, any or all of them may be set
c             to zero.  The effect will be the same as the values
c             relpr, relpr and 1/relpr respectively.
c
c     itnlim  input      An upper limit on the number of iterations.
c                        Suggested value:
c                        itnlim = n/2   for well-conditioned systems
c                                       with clustered singular values,
c                        itnlim = 4*n   otherwise.
c
c     nout    input      File number for printed output.  If positive,
c                        a summary will be printed on file nout.
c
c     istop   output     An integer giving the reason for termination:
c
c                0       x = 0  is the exact solution.
c                        No iterations were performed.
c
c                1       The equations A*x = b are probably
c                        compatible.  Norm(A*x - b) is sufficiently
c                        small, given the values of atol and btol.
c
c                2       damp is zero.  The system A*x = b is probably
c                        not compatible.  A least-squares solution has
c                        been obtained that is sufficiently accurate,
c                        given the value of atol.
c
c                3       damp is nonzero.  A damped least-squares
c                        solution has been obtained that is sufficiently
c                        accurate, given the value of atol.
c
c                4       An estimate of cond(Abar) has exceeded
c                        conlim.  The system A*x = b appears to be
c                        ill-conditioned.  Otherwise, there could be an
c                        error in SUBROUTINE aprod.
c
c                5       The iteration limit itnlim was reached.
c
c     itn     output     The number of iterations performed.
c
c     anorm   output     An estimate of the Frobenius norm of  Abar.
c                        This is the square-root of the sum of squares
c                        of the elements of Abar.
c                        If damp is small and if the columns of A
c                        have all been scaled to have length 1.0,
c                        anorm should increase to roughly sqrt(n).
c                        A radiCALLy different value for anorm may
c                        indicate an error in SUBROUTINE aprod (there
c                        may be an inconsistency between modes 1 and 2).
c
c     acond   output     An estimate of cond(Abar), the condition
c                        number of Abar.  A very high value of acond
c                        may again indicate an error in aprod.
c
c     rnorm   output     An estimate of the final value of norm(rbar),
c                        the function being minimized (see notation
c                        above).  This will be small if A*x = b has
c                        a solution.
c
c     arnorm  output     An estimate of the final value of
c                        norm( Abar(transpose)*rbar ), the norm of
c                        the residual for the usual normal equations.
c                        This should be small in all cases.  (arnorm
c                        will often be smaller than the true value
c                        computed from the output vector x.)
c
c     xnorm   output     An estimate of the norm of the final
c                        solution vector x.
c
c
c     Subroutines and functions used
c     ------------------------------
c
c     USER               aprod
c     LSQR               d2norm
c     BLAS               dcopy, dnrm2, dscal (see Lawson et al. below)
c
c
c     Precision
c     ---------
c
c     The number of iterations required by LSQR will usually decrease
c     if the computation is performed in higher precision.
c     At least 15-digit arithmetic should normally be used.
c     To convert LSQR and D2NORM between single and REAL(rk),
c     change
c                        REAL(rk)
c                        dcopy, dnrm2, dscal
c     to the appropriate FORTRAN and BLAS equivalents.
c     Also change 'd+' or 'e+' in the parameter statement.
c
c
c     References
c     ----------
c
c     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
c          linear equations and sparse least squares,
c          ACM Transactions on Mathematical Software 8, 1 (March 1982),
c          pp. 43-71.
c
c     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
c          linear equations and least-squares problems,
c          ACM Transactions on Mathematical Software 8, 2 (June 1982),
c          pp. 195-209.
c
c     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
c          Basic linear algebra subprograms for Fortran usage,
c          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
c          pp. 308-323 and 324-325.
c     ------------------------------------------------------------------
c
c
c     LSQR development:
c     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
c     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
c     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
c                     IF ( (one + ABS(t)) .LE. one ) GO TO 200
c                  from loop 200.  The test was an attempt to reduce
c                  underflows, but caused w(i) not to be updated.
c     17 Mar 1989: First F77 version.
c     04 May 1989: Bug (David Gay, AT&T).  When the second beta is zero,
c                  rnorm = 0 and
c                  test2 = arnorm / (anorm * rnorm) overflows.
c                  Fixed by testing for rnorm = 0.
c     05 May 1989: Sent to "misc" in netlib.
c     14 Mar 1990: Bug (John Tomlin via IBM OSL testing).
c                  Setting rhbar2 = rhobar**2 + dampsq can give zero
c                  if rhobar underflows and damp = 0.
c                  Fixed by testing for damp = 0 specially.
c     15 Mar 1990: Converted to lower case.
c     21 Mar 1990: d2norm introduced to avoid overflow in numerous
c                  items like  c = sqrt( a**2 + b**2 ).
c     04 Sep 1991: wantse added as an argument to LSQR, to make
c                  standard errors optional.  This saves storage and
c                  time when se(*) is not wanted.
c     13 Feb 1992: istop now RETURNs a value in [1,5], not [1,7].
c                  1, 2 or 3 means that x solves one of the problems
c                  Ax = b,  min norm(Ax - b)  or  damped least squares.
c                  4 means the limit on cond(A) was reached.
c                  5 means the limit on iterations was reached.
c     07 Dec 1994: Keep track of dxmax = max_k norm( phi_k * d_k ).
c                  So far, this is just printed at the end.
c                  A large value (relative to norm(x)) indicates
c                  significant cancellation in forming
c                  x  =  D*f  =  sum( phi_k * d_k ).
c                  A large column of D need NOT be serious if the
c                  corresponding phi_k is small.
c     27 Dec 1994: Include estimate of alfa_opt in iteration log.
c                  alfa_opt is the optimal scale factor for the
c                  residual in the "augmented system", as described by
c                  A. Bjorck (1992),
c                  Pivoting and stability in the augmented system method,
c                  in D. F. Griffiths and G. A. Watson (eds.),
c                  "Numerical Analysis 1991",
c                  Proceedings of the 14th Dundee Conference,
c                  Pitman Research Notes in Mathematics 260,
c                  Longman Scientific and Technical, Harlow, Essex, 1992.
c
c
c     Michael A. Saunders                  mike@sol-michael.stanford.edu
c     Dept of Operations Research          na.Msaunders@na-net.ornl.gov
c     Stanford University
c     Stanford, CA 94305-4022              (415) 723-1875
c **                                                                  **
c **********************************************************************

c     EXTERNAL           d2norm, dnrm2, dcopy, dscal
      REAL(rk)   lsqr_d2norm, dnrm2

c     Local variables

      LOGICAL            damped, extra
      INTEGER            i, maxdx, nconv, nstop
      REAL(rk)   alpha, beta, bnorm,
     &                   cs, cs1, cs2, ctol,
     &                   delta, dknorm, dnorm, dxk, dxmax,
     &                   gamma, gambar, phi, phibar, psi,
     &                   res2, rho, rhobar, rhbar1,
     &                   rhs, rtol, sn, sn1, sn2,
     &                   t, tau, temp, test1, test2, test3,
     &                   theta, t1, t2, t3, xnorm1, z, zbar, alfopt

      REAL(rk)   zero,           one
      PARAMETER        ( zero = 0.0,  one = 1.0 )

      CHARACTER*14       enter, exit
      CHARACTER*53       msg(0:5)

      DATA               enter /' Enter LSQR.  '/
      DATA               exit  /' Exit  LSQR.  '/
      DATA               msg
     &/ 'The exact solution is  x = 0',
     &'A solution to Ax = b was found, given atol, btol',
     &'A least-squares solution was found, given atol',
     &'A damped least-squares solution was found, given atol',
     &'Cond(Abar) seems to be too large, given conlim',
     &'The iteration limit was reached' /
c-----------------------------------------------------------------------

c     Initialize.
      nout=-1
      IF (nout .GT. 0) THEN
        WRITE(nout, 1000) enter, m, n, damp, wantse,
     &  atol, conlim, btol, itnlim
      END IF

      damped =   damp .GT. zero
      extra  =  .true.     ! true for extra printing below.
      itn    =   0
      istop  =   0
      nstop  =   0
      maxdx  =   0
      ctol   =   zero
      IF (conlim .GT. zero) ctol = one / conlim
      anorm  =   zero
      acond  =   zero
      dnorm  =   zero
      dxmax  =   zero
      res2   =   zero
      psi    =   zero
      xnorm  =   zero
      xnorm1 =   zero
      cs2    = - one
      sn2    =   zero
      z      =   zero

c     ------------------------------------------------------------------
c     Set up the first vectors u and v for the bidiagonalization.
c     These satisfy  beta*u = b,  alpha*v = A(transpose)*u.
c     ------------------------------------------------------------------
      DO 10  i = 1, n
        v(i)  =  zero
        x(i)  =  zero
   10 CONTINUE

      IF ( wantse ) THEN
        DO 20  i = 1, n
          se(i) =  zero
   20   CONTINUE
      END IF

      alpha  =   zero
      beta   =   dnrm2 ( m, u, 1 )  ! to paraleliziraj u(zac:kon), racuna normo

      IF (beta .GT. zero) THEN
        CALL dscal ( m, (one / beta), u, 1 ) ! u=u*(one/beta)
        CALL lsqr_aprod ( 2, m, n, nnz, v, u, ia, ja, a ) ! paraleliziraj
        alpha  =   dnrm2 ( n, v, 1 )
      END IF

      IF (alpha .GT. zero) THEN
        CALL dscal ( n, (one / alpha), v, 1 )
        CALL dcopy ( n, v, 1, w, 1 )
      END IF

      arnorm =   alpha * beta
      IF (arnorm .EQ. zero) GOTO 800

      rhobar =   alpha
      phibar =   beta
      bnorm  =   beta
      rnorm  =   beta

      IF (nout   .GT.  0  ) THEN
        IF ( damped ) THEN
          WRITE(nout, 1300)
        ELSE
          WRITE(nout, 1200)
        END IF
        test1  = one
        test2  = alpha / beta

        IF ( extra ) THEN
          WRITE(nout, 1400)
        END IF
        WRITE(nout, 1500) itn, x(1), rnorm, test1, test2
        WRITE(nout, 1600)
      END IF

c     ==================================================================
c     Main iteration loop.
c     ==================================================================
  100 itn    = itn + 1

c     ------------------------------------------------------------------
c     Perform the next step of the bidiagonalization to obtain the
c     next  beta, u, alpha, v.  These satisfy the relations
c                beta*u  =  A*v  -  alpha*u,
c               alpha*v  =  A(transpose)*u  -  beta*v.
c     ------------------------------------------------------------------
      CALL dscal ( m, (- alpha), u, 1 )
      CALL lsqr_aprod ( 1, m, n, nnz, v, u, ia, ja, a )
      beta   =   dnrm2 ( m, u, 1 )

c     Accumulate  anorm = || Bk ||
c                       =  sqrt( sum of  alpha**2 + beta**2 + damp**2 ).

      temp   =   lsqr_d2norm( alpha, beta )
      temp   =   lsqr_d2norm( temp , damp )
      anorm  =   lsqr_d2norm( anorm, temp )

      IF (beta .GT. zero) THEN
        CALL dscal ( m, (one / beta), u, 1 )
        CALL dscal ( n, (- beta), v, 1 )
        CALL lsqr_aprod ( 2, m, n, nnz, v, u, ia, ja, a )
        alpha  =   dnrm2 ( n, v, 1 )
        IF (alpha .GT. zero) THEN
          CALL dscal ( n, (one / alpha), v, 1 )
        END IF
      END IF

c     ------------------------------------------------------------------
c     Use a plane rotation to eliminate the damping parameter.
c     This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
c     ------------------------------------------------------------------
      rhbar1 = rhobar
      IF ( damped ) THEN
        rhbar1 = lsqr_d2norm( rhobar, damp )
        cs1    = rhobar / rhbar1
        sn1    = damp   / rhbar1
        psi    = sn1 * phibar
        phibar = cs1 * phibar
      END IF

c     ------------------------------------------------------------------
c     Use a plane rotation to eliminate the subdiagonal element (beta)
c     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
c     ------------------------------------------------------------------
      rho    =   lsqr_d2norm( rhbar1, beta )
      cs     =   rhbar1 / rho
      sn     =   beta   / rho
      theta  =   sn * alpha
      rhobar = - cs * alpha
      phi    =   cs * phibar
      phibar =   sn * phibar
      tau    =   sn * phi

c     ------------------------------------------------------------------
c     Update  x, w  and (perhaps) the standard error estimates.
c     ------------------------------------------------------------------
      t1     =   phi   / rho
      t2     = - theta / rho
      t3     =   one   / rho
      dknorm =   zero

      IF ( wantse ) THEN
        DO 200  i =  1, n
          t      =  w(i)
          x(i)   =  t1*t  +  x(i)
          w(i)   =  t2*t  +  v(i)
          t      = (t3*t)**2
          se(i)  =  t     +  se(i)
          dknorm =  t     +  dknorm
  200   CONTINUE
      ELSE
        DO 220  i =  1, n
          t      =  w(i)
          x(i)   =  t1*t  +  x(i)
          w(i)   =  t2*t  +  v(i)
          dknorm = (t3*t)**2  +  dknorm
  220   CONTINUE
      END IF

c     ------------------------------------------------------------------
c     Monitor the norm of d_k, the update to x.
c     dknorm = norm( d_k )
c     dnorm  = norm( D_k ),        where   D_k = (d_1, d_2, ..., d_k )
c     dxk    = norm( phi_k d_k ),  where new x = x_k + phi_k d_k.
c     ------------------------------------------------------------------
      dknorm = SQRT( dknorm )
      dnorm  = lsqr_d2norm( dnorm, dknorm )
      dxk    = ABS( phi * dknorm )
      IF (dxmax .LT. dxk ) THEN
        dxmax   =  dxk
        maxdx   =  itn
      END IF

c     ------------------------------------------------------------------
c     Use a plane rotation on the right to eliminate the
c     super-diagonal element (theta) of the upper-bidiagonal matrix.
c     Then use the result to estimate  norm(x).
c     ------------------------------------------------------------------
      delta  =   sn2 * rho
      gambar = - cs2 * rho
      rhs    =   phi    - delta * z
      zbar   =   rhs    / gambar
      xnorm  =   lsqr_d2norm( xnorm1, zbar  )
      gamma  =   lsqr_d2norm( gambar, theta )
      cs2    =   gambar / gamma
      sn2    =   theta  / gamma
      z      =   rhs    / gamma
      xnorm1 =   lsqr_d2norm( xnorm1, z     )

c     ------------------------------------------------------------------
c     Test for convergence.
c     First, estimate the norm and condition of the matrix  Abar,
c     and the norms of  rbar  and  Abar(transpose)*rbar.
c     ------------------------------------------------------------------
      acond  =   anorm * dnorm
      res2   =   lsqr_d2norm( res2 , psi    )
      rnorm  =   lsqr_d2norm( res2 , phibar )
      arnorm =   alpha * ABS( tau )

c     Now use these norms to estimate certain other quantities,
c     some of which will be small near a solution.

      alfopt =   SQRT( rnorm / (dnorm * xnorm) )
      test1  =   rnorm /  bnorm
      test2  =   zero
      IF (rnorm .GT. zero) test2 = arnorm / (anorm * rnorm)
      test3  =   one   /  acond
      t1     =   test1 / (one  +  anorm * xnorm / bnorm)
      rtol   =   btol  +  atol *  anorm * xnorm / bnorm

c     The following tests guard against extremely small values of
c     atol, btol  or  ctol.  (The user may have set any or all of
c     the parameters  atol, btol, conlim  to zero.)
c     The effect is equivalent to the normal tests using
c     atol = relpr,  btol = relpr,  conlim = 1/relpr.

      t3     =   one + test3
      t2     =   one + test2
      t1     =   one + t1
      IF (itn .GE. itnlim) istop = 5
      IF (t3  .LE. one   ) istop = 4
      IF (t2  .LE. one   ) istop = 2
      IF (t1  .LE. one   ) istop = 1

c     Allow for tolerances set by the user.

      IF (test3 .LE. ctol) istop = 4
      IF (test2 .LE. atol) istop = 2
      IF (test1 .LE. rtol) istop = 1

c     ------------------------------------------------------------------
c     See if it is time to print something.
c     ------------------------------------------------------------------
      IF (nout  .LE.  0       ) GOTO 600
      IF (n     .LE. 40       ) GOTO 400
      IF (itn   .LE. 10       ) GOTO 400
      IF (itn   .GE. itnlim-10) GOTO 400
      IF (MOD(itn,10) .EQ. 0  ) GOTO 400
      IF (test3 .LE.  2.0*ctol) GOTO 400
      IF (test2 .LE. 10.0*atol) GOTO 400
      IF (test1 .LE. 10.0*rtol) GOTO 400
      IF (istop .NE.  0       ) GOTO 400
      GOTO 600

c     Print a line for this iteration.
c     "extra" is for experimental purposes.

  400 IF ( extra ) THEN
        WRITE(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
     &  , phi, dknorm, dxk, alfopt
      ELSE
        WRITE(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
      END IF
      IF (MOD(itn,10) .EQ. 0) WRITE(nout, 1600)

c     ------------------------------------------------------------------
c     Stop if appropriate.
c     The convergence criteria are required to be met on  nconv
c     consecutive iterations, where  nconv  is set below.
c     Suggested value:  nconv = 1, 2  or  3.
c     ------------------------------------------------------------------
  600 IF (istop .EQ. 0) THEN
        nstop  = 0
      ELSE
        nconv  = 1
        nstop  = nstop + 1
        IF (nstop .LT. nconv  .AND.  itn .LT. itnlim) istop = 0
      END IF
      IF (istop .EQ. 0) GOTO 100

c     ==================================================================
c     End of iteration loop.
c     ==================================================================

c     Finish off the standard error estimates.

      IF ( wantse ) THEN
        t    =   one
        IF (m .GT. n)  t = m - n
        IF ( damped )  t = m
        t    =   rnorm / SQRT( t )

        DO 700  i = 1, n
          se(i)  = t * SQRT( se(i) )
  700   CONTINUE
      END IF

c     Decide if istop = 2 or 3.
c     Print the stopping condition.

  800 IF (damped  .AND.  istop .EQ. 2) istop = 3
      IF (nout .GT. 0) THEN
        WRITE(nout, 2000) EXIT, istop, itn,
     &  EXIT, anorm, acond,
     &  EXIT, bnorm, xnorm,
     &  EXIT, rnorm, arnorm
        WRITE(nout, 2100) EXIT, dxmax, maxdx,
     &  EXIT, dxmax/(xnorm + 1.0D-20)
        WRITE(nout, 3000) EXIT, msg(istop)
      END IF

  900 RETURN

c     ------------------------------------------------------------------
 1000 format(// 1p, a, '     Least-squares solution of  Ax = b'
     &/ ' The matrix  A  has', i7, ' rows   and', i7, ' columns'
     &/ ' damp   =', e22.14, 3x,        'wantse =', l10
     &/ ' atol   =', e10.2, 15x,        'conlim =', e10.2
     &/ ' btol   =', e10.2, 15x,        'itnlim =', i10)
 1200 format(// '   Itn       x(1)           Function',
     &'     Compatible   LS        Norm A    Cond A')
 1300 format(// '   Itn       x(1)           Function',
     &'     Compatible   LS     Norm Abar Cond Abar')
 1400 format(80x, '    phi    dknorm   dxk  alfa_opt')
 1500 format(1p, i6, 2e17.9, 4e10.2, e9.1, 3e8.1)
 1600 format(1x)
 2000 format(/ 1p, a, 5x, 'istop  =', i2,   15x, 'itn    =', i8
     &/     a, 5x, 'anorm  =', e12.5, 5x, 'acond  =', e12.5
     &/     a, 5x, 'bnorm  =', e12.5, 5x, 'xnorm  =', e12.5
     &/     a, 5x, 'rnorm  =', e12.5, 5x, 'arnorm =', e12.5)
 2100 format(  1p, a, 5x, 'max dx =', e8.1 , ' occurred at itn ', i8,
     &/     a, 5x, '       =', e8.1 , '*xnorm' )
 3000 format( a, 5x, a )

c     End of LSQR
      END


C----------------------------------------------------------------------C
c                                                                      c
      FUNCTION lsqr_d2norm( a, b )
      USE mCommon            
      IMPLICIT NONE
      REAL(rk)  lsqr_d2norm, a, b
c                                                                      c
C----------------------------------------------------------------------C
c **********************************************************************
c **                                                                  **
c     d2norm  returns  sqrt( a**2 + b**2 )  with precautions
c     to avoid overflow.
c
c     21 Mar 1990: First version.
c **                                                                  **
c **********************************************************************
      INTRINSIC         abs, sqrt
      REAL(rk)  scale
      REAL(rk)  zero
      PARAMETER       ( zero = 0.0d+0 )

      scale  = abs( a ) + abs( b )
      IF (scale .EQ. zero) THEN
        lsqr_d2norm = zero
      ELSE
        lsqr_d2norm = scale * sqrt( (a/scale)**2   +  (b/scale)**2 )
      END IF
c
      RETURN
      END
