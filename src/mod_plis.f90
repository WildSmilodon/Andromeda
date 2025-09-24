module plis

   use parallel
   use mCommon
   use mPar
   implicit none
 
   private
#include "lisf.h" 
    
   TYPE systemLinEq
      LIS_INTEGER neq  ! number of equations in SLE
      LIS_INTEGER n    ! number of rows on processor
      LIS_INTEGER is,ie! start row and end row + 1      

      LIS_INTEGER nnz    
      LIS_MATRIX A      ! system matrix
      LIS_SOLVER solver ! pointer to solver
      LIS_VECTOR b,x    ! r.h.s. and solution vectors

      REAL(rk), POINTER :: rhsM(:,:) ! parallel version of rectangular RHS matrix
      INTEGER nb

      REAL(rk), ALLOCATABLE :: sm_val(:) ! element values array
      INTEGER, ALLOCATABLE :: sm_rowS(:) ! start of row in CSR format
      INTEGER, ALLOCATABLE :: sm_colI(:) ! column indices in CSR format

   END TYPE

   public :: systemLinEq
   public :: plis_setup_A, plis_setup_bx, plis_setup_Abx
   public :: plis_set_matrix_element, plis_set_matrix_element_check
   public :: plis_assemble_system_matrix, plis_assemble_system_matrix_csr
   public :: plis_set_RHSmatrix_element, plis_init_rhs_matrix
   public :: plis_solve
   public :: plis_get_solver_status, plis_get_solver_iterations
   public :: plis_get_solver_name, plis_get_solver_time
   public :: plis_createSolver
   public :: plis_get_nnz
   public :: plis_cleanup
   public :: plis_getX
   public :: plis_is_row_mine
   public :: plis_fill_121
   public :: plis_rhs_matvec
   public :: plis_print_system_matrix

   

   contains 




   ! ----------------------------------------------------------------------------      
   subroutine plis_print_system_matrix(sle)

      TYPE(systemLinEq) sle 
      integer i,j

      DO i=0,sle%n
         DO j=sle%sm_rowS(i),sle%sm_rowS(i+1)-1
           print *,i+1,sle%sm_colI(j)+1,sle%sm_val(j)
         END DO
      END DO

   end subroutine




   ! ----------------------------------------------------------------------------      
   subroutine plis_assemble_system_matrix(sle)

      TYPE(systemLinEq) sle 
      call lis_matrix_set_type(sle%A,LIS_MATRIX_DNS,env%ierr)
      call lis_matrix_assemble(sle%A,env%ierr)

   end subroutine


   ! ----------------------------------------------------------------------------      
   subroutine plis_assemble_system_matrix_csr(nnz,sle)

      TYPE(systemLinEq) sle 
      INTEGER nnz
      call lis_matrix_set_csr(nnz,sle%sm_rowS,sle%sm_colI,sle%sm_val,sle%A,env%ierr)
      call lis_matrix_assemble(sle%A,env%ierr)
   end subroutine


   

   ! ----------------------------------------------------------------------------      
   subroutine plis_rhs_matvec(sle,v)
      INTEGER r,c,row
      TYPE(systemLinEq) sle 
      REAL(rk) v(sle%nb)

      ! set result vector to 0
      call lis_vector_set_all(0.0_rk,sle%b,env%ierr)

      DO c = 1, sle%nb
         DO r = 1, sle%n
            row = r + sle%is - 1
            call lis_vector_set_value(LIS_ADD_VALUE,row, sle%rhsM(r,c)*v(c), sle%b,env%ierr)
         END DO
      END DO

   end subroutine

   ! ----------------------------------------------------------------------------      
   subroutine plis_init_rhs_matrix(sle,nb)
      INTEGER nb,r,c
      TYPE(systemLinEq) sle 

      sle%nb = nb
      ALLOCATE (sle%rhsM(sle%n,sle%nb))

      DO c=1,sle%nb
         DO r=1,sle%n
            sle%rhsM(r,c) = 0.0_rk
         end do
      end do

   end subroutine

   ! ----------------------------------------------------------------------------   
   subroutine plis_set_RHSmatrix_element(sle,row,col,value)
   
      INTEGER row,col ! row is global
      REAL(8) value
      TYPE(systemLinEq) sle ! system to be solved

      sle%rhsM(row - sle%is + 1, col) = value

   end subroutine     

   ! ----------------------------------------------------------------------------      
   logical function plis_is_row_mine(sle,row)
      integer row
      TYPE(systemLinEq) sle 

      logical result
      result = .false.
      if (row.ge.sle%is.and.row.lt.sle%ie) result = .true.
      plis_is_row_mine = result
   end function

   ! ----------------------------------------------------------------------------   
   subroutine plis_set_matrix_element(sle,row,col,val)
   
      INTEGER row,col
      REAL(8) val
      TYPE(systemLinEq) sle ! system to be solved

      call lis_matrix_set_value(LIS_INS_VALUE,row,col,val,sle%A,env%ierr)

   end subroutine     


   ! ----------------------------------------------------------------------------   
   subroutine plis_set_matrix_element_check(sle,row,col,value)
   
      INTEGER row,col
      REAL(8) value
      TYPE(systemLinEq) sle ! system to be solved

      if (row.ge.sle%is.and.row.lt.sle%ie) call lis_matrix_set_value(LIS_INS_VALUE,row,col,value,sle%A,env%ierr)

   end subroutine     


   ! ----------------------------------------------------------------------------
   subroutine plis_solve(sle)

      TYPE(systemLinEq) sle ! system to be solved

      call lis_solve(sle%A,sle%b,sle%x,sle%solver,env%ierr)      

   end subroutine


   ! ----------------------------------------------------------------------------
   INTEGER function plis_get_nnz(sle)
      LIS_INTEGER nnz
      TYPE(systemLinEq) sle ! system to be solved
      
      call lis_matrix_get_nnz(sle%A, nnz, env%ierr)
      plis_get_nnz = nnz
   end function

   ! ----------------------------------------------------------------------------
   INTEGER function plis_get_solver_status(sle)
      LIS_INTEGER ierr
      TYPE(systemLinEq) sle ! system to be solved
      
      call lis_solver_get_status(sle%solver, ierr, env%ierr)
      plis_get_solver_status = ierr
   end function

   ! ----------------------------------------------------------------------------
   INTEGER function plis_get_solver_iterations(sle)
      LIS_INTEGER iter
      TYPE(systemLinEq) sle ! system to be solved
      
      call lis_solver_get_iter(sle%solver, iter, env%ierr)
      plis_get_solver_iterations = iter
   end function

   ! ----------------------------------------------------------------------------
   CHARACTER(255) function plis_get_solver_name(sle)
      CHARACTER(255) name
      TYPE(systemLinEq) sle ! system to be solved
      call lis_solver_get_solvername(sle%solver, name, env%ierr)
      plis_get_solver_name = name
   end function   

   ! ----------------------------------------------------------------------------
   REAL(8) function plis_get_solver_time(sle)
      REAL(8) time
      TYPE(systemLinEq) sle ! system to be solved
      
      call lis_solver_get_time(sle%solver, time, env%ierr)
      plis_get_solver_time = time
   end function


   ! ----------------------------------------------------------------------------
   subroutine plis_createSolver(sle)

      TYPE(systemLinEq) sle ! system to be solved

      call lis_solver_create(sle%solver,env%ierr)

      call lis_solver_set_option(parLISslvSet,sle%solver,env%ierr)
      !call lis_solver_set_option("-i bicg -p none -tol 1.0e-8 -maxiter 10000",sle%solver,env%ierr)
      !call lis_solver_set_option("-i bicg -p jacobi -tol 1.0e-12 -print all",sle%solver,env%ierr)

   end subroutine


   ! ----------------------------------------------------------------------------
   subroutine plis_setUpVector(vector,size)

      LIS_VECTOR vector 
      LIS_INTEGER size
 
      call lis_vector_create(env%comm,vector,env%ierr)      
      call lis_vector_set_size(vector,0,size,env%ierr)
      call lis_vector_set_all(0.0D0,vector,env%ierr)

   end subroutine

   ! ----------------------------------------------------------------------------
   subroutine plis_vector_print(vector)
 
      INTEGER i
      LIS_VECTOR vector
      LIS_INTEGER nl,ng
      REAL(8), ALLOCATABLE :: v(:)

      call lis_vector_get_size(vector,nl,ng,env%ierr)
      ALLOCATE (v(ng))
      call lis_vector_gather(vector,v,env%ierr)

      if (env%myRank.EQ.0) then
         do i = 1,ng
            print *,i,v(i)
         end do
      end if

      DEALLOCATE(v)      

   end subroutine


   ! ----------------------------------------------------------------------------
   subroutine plis_setup_A(sle,numberOfEquations)

      integer numberOfEquations
      TYPE(systemLinEq) sle ! system to be solved

      sle%neq = numberOfEquations

      call lis_matrix_create(env%comm,sle%A,env%ierr)
      call lis_matrix_set_size(sle%A,0,sle%neq,env%ierr)
      call lis_matrix_get_size(sle%A,sle%n,sle%neq,env%ierr)
      call lis_matrix_get_range(sle%A,sle%is,sle%ie,env%ierr)
   
   end subroutine

   ! ----------------------------------------------------------------------------
   subroutine plis_setup_bx(sle)

      TYPE(systemLinEq) sle ! system to be solved

      call plis_setUpVector(sle%x,sle%neq)
      call plis_setUpVector(sle%b,sle%neq)
   
   end subroutine


   ! ----------------------------------------------------------------------------
   subroutine plis_setup_Abx(sle,numberOfEquations)

      integer numberOfEquations
      TYPE(systemLinEq) sle ! system to be solved

      sle%neq = numberOfEquations

      call lis_matrix_create(env%comm,sle%A,env%ierr)
      call lis_matrix_set_size(sle%A,0,sle%neq,env%ierr)
      call lis_matrix_get_size(sle%A,sle%n,sle%neq,env%ierr)
      call lis_matrix_get_range(sle%A,sle%is,sle%ie,env%ierr)
      call plis_setUpVector(sle%x,sle%neq)
      call plis_setUpVector(sle%b,sle%neq)
   
   end subroutine

   ! ----------------------------------------------------------------------------
   subroutine plis_fill_121(sle)

      INTEGER idx,i
      TYPE(systemLinEq) sle ! system to be solved

      ! system matrix
      do i=sle%is,sle%ie-1
         if( i>1   ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,1.0d0,sle%A,env%ierr)
         if( i<sle%neq ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,1.0d0,sle%A,env%ierr)
         call lis_matrix_set_value(LIS_INS_VALUE,i,i,-2.0d0,sle%A,env%ierr)
      end do

      call lis_matrix_set_type(sle%A,LIS_MATRIX_CSR,env%ierr)
      call lis_matrix_assemble(sle%A,env%ierr)
      call lis_matrix_get_nnz(sle%A,sle%nnz,env%ierr)      

      ! right hand side
      idx = 1
      if (sle%is.le.idx.and.sle%ie.gt.idx) call lis_vector_set_value(LIS_INS_VALUE,idx,-10.0D0,sle%b,env%ierr)
      idx = sle%neq
      if (sle%is.le.idx.and.sle%ie.gt.idx) call lis_vector_set_value(LIS_INS_VALUE,idx,-20.0D0,sle%b,env%ierr)


   end subroutine


   ! ----------------------------------------------------------------------------
   subroutine plis_cleanup(sle)  

      TYPE(systemLinEq) sle ! system to be solved
      call lis_matrix_destroy(sle%A,env%ierr)
      call lis_vector_destroy(sle%b,env%ierr)
      call lis_vector_destroy(sle%x,env%ierr)

   end subroutine


   ! ----------------------------------------------------------------------------
   subroutine plis_getX(sle,x)  

      TYPE(systemLinEq) sle ! system to be solved
      real(8) x(sle%neq)

      call lis_vector_gather(sle%x,x,env%ierr)

   end subroutine



end module plis