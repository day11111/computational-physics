program trapezoidal
implicit none

    ! ----------------------------------------------------------------------------------
    ! Calculating the normal probability function integral using the trapezoidal rule
    ! ----------------------------------------------------------------------------------

    !the number of partitions is n
    !a and b are the beginning point and the ending point of the integral
    !the function is even, so let a = 0.0

    integer :: i, j
    real :: h, sum, x, f
    integer :: n = 10 , m
    real, parameter :: a = 0.0, b = 1.0

    do j = 1, 5
      m = n ** j
      h = (b-a)/real(m)  
      sum = 0.5*(f(a) + f(b))     
      do  i=1,m-1
        x = a + i*h
        sum = sum + f(x)   
      end do
    sum = 2*h*sum
    print *,"The trapezoidal rule : ", "n = ", m ,"the integral is : ", sum
    end do
end program trapezoidal
            
    function f (t)
    implicit none
    real, intent(in) :: t
    real :: f
    ! pi = 4.*atan(1.)
    f = sqrt(1./(8.*atan(1.))) /exp(t*t)
    end function f
    


