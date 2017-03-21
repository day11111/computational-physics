!      This program solves driven morse oscillator model. 
!      The hamiltonian is H(x,p) = p**2/2 + (1-exp(-x))**2 - a*x*sin(wt)
!      where a = 0.04, w = 0.25
!
!      We use a two-dimensional array to represent x and p
!      as functions of t
!      y[0] == x
!      y[1] == p
!      dy[0]/dt = p
!      dy[1]/dt = (-2)*(exp(-x)-exp(-2*x))+0.04_dp*sin(0.25_dp*t)
!
!      Start	99	trajectories	at	equally	spaced	energies	between	
!      E = 0.01 and E = 0.99 with x(t=0) = 0 and momentum chosen
!      appropriately. 	Compute	the	time	evolution	for	each	of		
!      these	trajectories	until	a	maximum	time	of t = 400.
!      	Make	a	plot	of	ionization	time	vs.	starting	energy	
!
!  this is the number of differential equations  as a global parameter

MODULE parameters
  INTEGER, PARAMETER, PUBLIC :: number_differential_eqs =2 
END MODULE parameters

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
END MODULE constants
!
!     Main function begins here 
!
PROGRAM driven_Morse
  USE constants
  USE parameters
  IMPLICIT NONE
  REAL(DP), DIMENSION(number_differential_eqs) :: y, dydt, yout
  REAL(DP) ::  t, E, initial_x = 0.00_dp, initial_E = 0.01_dp, initial_p, h = 0.10_dp

  do while (initial_E <= 0.99_dp)
    t=0.0_dp                                      ! initial time      
    initial_p = sqrt(2*initial_E)              ! the initial p
    y(1) = initial_x                          ! initial position  
    y(2) = initial_p                          ! initial momentum  
    E = initial_E                             ! initial energy

    ! now we start solving the differential equations using the RK4 method 
    OPEN(6,FILE='driven.dat', position='append')
    yout = 0.0_dp; dydt = 0.0_dp
    DO WHILE (E <= 1 .and. t <= 400)
      ! initial derivatives 
      CALL derivatives(t, y, dydt)                
      ! here we call the runge-kutta method and get the new y-value in yout                     
      CALL runge_kutta_4(y, t, h, yout, dydt)
      y = yout  
      t = t + h
      E = 0.5_dp*(y(2)**2) + (1-exp(-y(1)))**2 - 0.04_dp*y(1)*sin(0.25_dp*t)
      !print *, E, t
    ENDDO
    ! writing energy and time
    WRITE(6,'(5(E12.6,1X))') initial_E, t
    CLOSE (6)
    initial_E = initial_E + 0.01_dp
  end do

END PROGRAM driven_Morse


!
!   this function sets up the derivatives for this special case  
!
SUBROUTINE derivatives(t, y, dydt)
  USE constants
  USE parameters
  IMPLICIT NONE
  REAL(DP), DIMENSION(number_differential_eqs) :: y, dydt
  REAL(DP) :: t

  dydt(1) = y(2)    ! derivative of x
  dydt(2) = (-2)*(exp(-y(1))-exp(-2*y(1)))+0.04_dp*sin(0.25_dp*t) ! derivative of p

END SUBROUTINE derivatives

!
!     Runge-kutta procedure 
!
SUBROUTINE runge_kutta_4(y,x,diff_eq_step,yout,dydx)
  USE constants
  USE parameters
  IMPLICIT NONE
  REAL(DP), DIMENSION(number_differential_eqs) :: yt, dyt, dym
  REAL(DP), DIMENSION(number_differential_eqs), INTENT(IN) :: y, dydx
  REAL(DP), DIMENSION(number_differential_eqs), INTENT(OUT) :: yout
  REAL(DP) :: hh, h6, xh
  REAL(DP), INTENT(IN) :: x, diff_eq_step 

  hh=diff_eq_step*0.5_dp; h6=diff_eq_step/6.00_dp; xh=x+hh
  !     first rk-step
  yt=y+hh*dydx
  CALL derivatives(xh,yt,dyt)
  !     second rk-step
  yt=y+hh*dyt
  CALL derivatives(xh,yt,dym)
  !     third rk-step
  yt=y+diff_eq_step*dym;  dym=dyt+dym
  CALL derivatives(x+diff_eq_step,yt,dyt)
  !     fourth rk-step
  yout=y+h6*(dydx+dyt+2.00_dp*dym)

END SUBROUTINE runge_kutta_4

