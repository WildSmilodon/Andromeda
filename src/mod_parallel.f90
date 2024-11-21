module parallel

   implicit none
   
   include "mpif.h"
#include "lisf.h"       
      
   TYPE paralelEnv
      INTEGER ierr
      LIS_INTEGER comm
      INTEGER nProcs
      INTEGER myRank
   END TYPE
    
   TYPE(paralelEnv) env  ! paralel environment
   LOGICAL amIroot
   

   contains 
    
   ! ----------------------------------------------------------------------------
   subroutine par_init()
  
      call lis_initialize(env%ierr)
      env%comm = LIS_COMM_WORLD 
      call MPI_Comm_size(env%comm,env%nprocs,env%ierr)
      call MPI_Comm_rank(env%comm,env%myRank,env%ierr)

      if (env%myRank.eq.0) then 
         amIroot = .true.
      else
         amIroot = .false.
      end if

   end subroutine
 
   ! ----------------------------------------------------------------------------
   subroutine par_finalize()  
      !call lis_finalize(env%ierr)
   end subroutine   

   
end module parallel