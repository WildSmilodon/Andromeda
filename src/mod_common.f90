module mCommon 
    implicit none
!
!   REAL(rk) number representation (4=single, 8=double)
!    
    integer,  parameter :: rk = selected_REAL_kind(8)
    REAL(rk), parameter :: pi = 3.141592653589793238462643383_rk

end module  