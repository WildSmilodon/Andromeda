!
!     Linearised flow fields module
!     (c) Jure Ravnik, 2021
!
module linFlowFields
      
  use mCommon
  use mPar
  implicit none
      
  !
  ! Types
  !

  !
  ! Parameters
  !
  integer, parameter :: pipeFlow = 1
  integer, parameter :: poiseuilleFlow = 2
  integer, parameter :: ABCflow = 3
  integer, parameter :: Xvortex = 4
  integer, parameter :: Yvortex = 5
  integer, parameter :: Zvortex = 6
  integer, parameter :: Xstrain = 7
  integer, parameter :: Ystrain = 8
  integer, parameter :: Zstrain = 9
  integer, parameter :: XplugFlow = 10
  integer, parameter :: YplugFlow = 11
  integer, parameter :: ZplugFlow = 12
  !
  ! Variables
  !
  real(rk), private :: v0 ! velocity length scale
  real(rk), private :: A,B,C,H,L,D ! length scales
  real(rk), private :: nu ! frequency in ABC flow
  real(rk), private :: r(3) ! particle position in IFR 
  real(rk), private :: ep(4) ! particle orientation in Euler paramters in IFR
  real(rk), private :: velocityIFR(3) ! flow velocity at particle in IFR : v_0
  real(rk), private :: velGradTenIFR(3,3) ! flow velocity gradient tensor at particle in IFR 
  real(rk), private :: RM(3,3) ! rotation matrix
  real(rk), private :: RMT(3,3) ! transpose of the rotation matrix
  integer, private :: flowType 

  !
  ! Subroutines
  !
  contains

  !
  ! ----------------------------------------------------------------------
  !
  subroutine lff_getFlowVelocityPFR(x,y,z,u,debug)
    real(rk), intent(in) :: x,y,z
    real(rk), intent(out) :: u(3)
    logical, intent(in), optional :: debug
    real(rk) rr(3)

    if (present(debug)) then
      print *,"velocityIFR"
      print *,velocityIFR  ! hitrost, ki jo ima tekocina na mestu delca
    end if

    u(1) = velocityIFR(1)
    u(2) = velocityIFR(2)
    u(3) = velocityIFR(3)

    if (present(debug)) then
      print *,"velGradTenIFR"
      print *,velGradTenIFR
      print *,"xyz"
      print *,x,y,z
    end if

    ! x,y,z is a point in PFR, rr is the same point in IFR
    rr(1)=x
    rr(2)=y
    rr(3)=z
    rr = MATMUL(RMT,rr)

    ! here w is calculated
    u(1) = u(1) + velGradTenIFR(1,1) * rr(1)  + velGradTenIFR(2,1) * rr(2)  + velGradTenIFR(3,1) * rr(3)
    u(2) = u(2) + velGradTenIFR(1,2) * rr(1)  + velGradTenIFR(2,2) * rr(2)  + velGradTenIFR(3,2) * rr(3)
    u(3) = u(3) + velGradTenIFR(1,3) * rr(1)  + velGradTenIFR(2,3) * rr(2)  + velGradTenIFR(3,3) * rr(3)

    if (present(debug)) then
      print *,"w"
      print *,u
    end if

    u = MATMUL(RM,u)    

    if (present(debug)) then
      print *,"u"
      print *,u
    end if

  end subroutine

  !
  ! ----------------------------------------------------------------------
  !
  subroutine lff_setFlowVelocityAndGradIFR()

    real(rk) x,y,z
    x = r(1)
    y = r(2)
    z = r(3)
    !
    ! Laminar flow in a pipe
    !
    if (flowType.eq.pipeFlow) then

      !v0 = k*D*D/(32.0_rk*mu)

      velocityIFR(1) = 2.0_rk * v0 * (1.0_rk - 4.0_rk*(y*y+z*z)/(D*D))
      velocityIFR(2) = 0.0_rk
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = - 16.0_rk * v0 / (D*D) * y 
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = - 16.0_rk * v0 / (D*D) * z
      velGradTenIFR(3,2) = 0.0_rk
      velGradTenIFR(3,3) = 0.0_rk     

    end if
    !
    ! Poissuille 2D flow
    !    
    if (flowType.eq.poiseuilleFlow) then


      !v0 = k*H*H/(12.0_rk*mu)

      velocityIFR(1) = 6.0_rk * v0 * z/H * (1.0_rk - z/H)
      velocityIFR(2) = 0.0_rk
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 6.0_rk * v0 /  H * ( 1.0_rk - 2.0_rk * z /H)
      velGradTenIFR(3,2) = 0.0_rk
      velGradTenIFR(3,3) = 0.0_rk     


    end if
    !
    ! Arnold–Beltrami–Childress flow
    !
    if (flowType.eq.ABCflow) then

      velocityIFR(1) = v0 * ( A * sin(nu*z/L) + C * cos(nu*y/L) )
      velocityIFR(2) = v0 * ( B * sin(nu*x/L) + A * cos(nu*z/L) )
      velocityIFR(3) = v0 * ( C * sin(nu*y/L) + B * cos(nu*x/L) )

      velGradTenIFR(1,1) =   0.0_rk
      velGradTenIFR(1,2) =   B*nu/L * cos(nu*x/L)
      velGradTenIFR(1,3) = - B*nu/L * sin(nu*x/L)

      velGradTenIFR(2,1) = - C*nu/L * sin(nu*y/L)
      velGradTenIFR(2,2) =   0.0_rk
      velGradTenIFR(2,3) =   C*nu/L * cos(nu*y/L)

      velGradTenIFR(3,1) =   A*nu/L * cos(nu*z/L)
      velGradTenIFR(3,2) = - A*nu/L * sin(nu*z/L)
      velGradTenIFR(3,3) =   0.0_rk  

    end if
    !
    ! Plug flow X
    !
    if (flowType.eq.XplugFlow) then
      velocityIFR(1) = 1.0_rk
      velocityIFR(2) = 0.0_rk
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! Plug flow Y
    !
    if (flowType.eq.YplugFlow) then
      velocityIFR(1) = 0.0_rk
      velocityIFR(2) = 1.0_rk
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! Plug flow X
    !
    if (flowType.eq.ZplugFlow) then
      velocityIFR(1) = 0.0_rk
      velocityIFR(2) = 1.0_rk
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if

    !
    ! Vortex around x axis
    !
    if (flowType.eq.Xvortex) then
      velocityIFR(1) = 0.0_rk
      velocityIFR(2) = -y
      velocityIFR(3) = z

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 1.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) =-1.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! Vortex around y axis
    !
    if (flowType.eq.Yvortex) then
      velocityIFR(1) = z
      velocityIFR(2) = 0.0_rk
      velocityIFR(3) = -x

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) =-1.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 1.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! Vortex around z axis
    !
    if (flowType.eq.Zvortex) then
      velocityIFR(1) = -y
      velocityIFR(2) = x
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 1.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) =-1.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! Strain around x axis
    !
    if (flowType.eq.Xstrain) then
      velocityIFR(1) = 0.0_rk
      velocityIFR(2) = y
      velocityIFR(3) = z

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 1.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) = 1.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! strain around y axis
    !
    if (flowType.eq.Ystrain) then
      velocityIFR(1) = z
      velocityIFR(2) = 0.0_rk
      velocityIFR(3) = x

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 0.0_rk
      velGradTenIFR(1,3) = 1.0_rk

      velGradTenIFR(2,1) = 0.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 1.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if
    !
    ! Strain around z axis
    !
    if (flowType.eq.Zstrain) then
      velocityIFR(1) = -y
      velocityIFR(2) = x
      velocityIFR(3) = 0.0_rk

      velGradTenIFR(1,1) = 0.0_rk
      velGradTenIFR(1,2) = 1.0_rk
      velGradTenIFR(1,3) = 0.0_rk

      velGradTenIFR(2,1) = 1.0_rk
      velGradTenIFR(2,2) = 0.0_rk
      velGradTenIFR(2,3) = 0.0_rk

      velGradTenIFR(3,1) = 0.0_rk
      velGradTenIFR(3,2) = 0.0_rk 
      velGradTenIFR(3,3) = 0.0_rk  
    end if

    

  end subroutine

  !
  ! ----------------------------------------------------------------------
  !
  subroutine lff_setParticlePosition(partX, partY, partZ)
    real(rk) partX, partY, partZ
    r(1) = partX
    r(2) = partY
    r(3) = partZ
  end subroutine

  !
  ! ----------------------------------------------------------------------
  !
  subroutine lff_setParticleOrientation(e0,e1,e2,e3)
    real(rk) e0,e1,e2,e3
    ep(1) = e0
    ep(2) = e1
    ep(3) = e2
    ep(4) = e3
  end subroutine

  !
  ! ----------------------------------------------------------------------
  !
  subroutine lff_setRotationMatrix()

    RM(1,1)=  ep(1)**2 + ep(2)**2 - ep(3)**2 - ep(4)**2
    RM(1,2)=  2.0D0*(ep(2)*ep(3) + ep(1)*ep(4))
    RM(1,3)=  2.0D0*(ep(2)*ep(4) - ep(1)*ep(3))

    RM(2,1)=  2.0D0*(ep(2)*ep(3) - ep(1)*ep(4))
    RM(2,2)=  ep(1)**2 - ep(2)**2 + ep(3)**2 - ep(4)**2
    RM(2,3)=  2.0D0*(ep(3)*ep(4) + ep(1)*ep(2))

    RM(3,1)=  2.0D0*(ep(2)*ep(4) + ep(1)*ep(3))
    RM(3,2)=  2.0D0*(ep(3)*ep(4) - ep(1)*ep(2))
    RM(3,3)=  ep(1)**2 - ep(2)**2 - ep(3)**2 + ep(4)**2

    RMT(1,1) = RM(1,1)
    RMT(1,2) = RM(2,1)
    RMT(1,3) = RM(3,1)

    RMT(2,1) = RM(1,2)
    RMT(2,2) = RM(2,2)
    RMT(2,3) = RM(3,2)

    RMT(3,1) = RM(1,3)
    RMT(3,2) = RM(2,3)
    RMT(3,3) = RM(3,3)



  end subroutine

  !
  ! ----------------------------------------------------------------------
  !
  subroutine lff_init()

    flowType = parFlopFlowType
    !
    ! Laminar flow in a pipe
    !
    if (flowType.eq.pipeFlow) then
      v0 = 1.0_rk
      D = 1.0_rk
    end if
    !
    ! Poissuille 2D flow
    !    
    if (flowType.eq.poiseuilleFlow) then
      v0 = 1.0_rk
      H = 1.0_rk
    end if
    !
    ! Arnold–Beltrami–Childress flow
    !
    if (flowType.eq.ABCflow) then
      A = 1.0_rk
      B = 1.0_rk
      C = 1.0_rk
      L = 1.0_rk
      nu = 2.0_rk * pi
    end if
    
    if (flowType.gt.0) then ! use special predefinded flow types
      call lff_setParticlePosition( parFlopUgU(5), parFlopUgU(6), parFlopUgU(7) ) 
      call lff_setFlowVelocityAndGradIFR()
    else ! use velocity + grad(vel) from input file
      velocityIFR(1) = parFlopUgU(5)
      velocityIFR(2) = parFlopUgU(6)
      velocityIFR(3) = parFlopUgU(7)

      velGradTenIFR(1,1) = parFlopUgU(8)
      velGradTenIFR(1,2) = parFlopUgU(9)
      velGradTenIFR(1,3) = parFlopUgU(10)

      velGradTenIFR(2,1) = parFlopUgU(11)
      velGradTenIFR(2,2) = parFlopUgU(12)
      velGradTenIFR(2,3) = parFlopUgU(13)

      velGradTenIFR(3,1) = parFlopUgU(14)
      velGradTenIFR(3,2) = parFlopUgU(15)
      velGradTenIFR(3,3) = parFlopUgU(16)     
    end if

    WRITE (parLogTekst,'(A)') "Creating BC's for flow over particle:"
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(A)') "Flow velocity at particle position in IFR"
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(3G18.10)') velocityIFR
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(A)') "Flow velocity gradient at particle position in IFR"
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(3G18.10)') velGradTenIFR(1,1),velGradTenIFR(1,2),velGradTenIFR(1,3)
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(3G18.10)') velGradTenIFR(2,1),velGradTenIFR(2,2),velGradTenIFR(2,3)
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(3G18.10)') velGradTenIFR(3,1),velGradTenIFR(3,2),velGradTenIFR(3,3)
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(A)') "Particle orientation - Euler parameters"
    CALL WriteToLog(parLogTekst)
    WRITE (parLogTekst,'(4G18.10)') parFlopUgU(1),parFlopUgU(2), parFlopUgU(3), parFlopUgU(4)
    CALL WriteToLog(parLogTekst)


    call lff_setParticleOrientation(parFlopUgU(1),parFlopUgU(2), parFlopUgU(3), parFlopUgU(4))
    call lff_setRotationMatrix()
  
  end subroutine lff_init
 

end module