program simpson
implicit none

    ! ----------------------------------------------------------------------------------
    ! Calculating the normal probability function integral using the simpson's rule
    ! ----------------------------------------------------------------------------------

    !the number of partitions is n
    !a and b are the beginning point and the ending point of the integral
    !the function is even, so let a = 0.0
    !we only calculate integral from 0.0 to 1.0

    integer :: i, j
    real :: h, sum, x, f
    integer :: n = 10, m
    real, parameter :: a = 0.0, b = 1.0

    do j = 1, 5
        m = n ** j
        h = (b-a)/real(m)  
        sum = f(a) + f(b)     
        do  i=1,m-1,2
        x = a + i*h
        sum = sum + 4*f(x)   
        end do
        do  i=2,m-1,2
        x = a + i*h
        sum = sum + 2*f(x)   
        end do
        sum = 2*h/3*sum
        print *,"The simpson's rule : ", "n = ", m ,"the integral is : ", sum
    end do
end program simpson
            
    function f (t)
    implicit none
    real, intent(in) :: t
    real :: f
    ! pi = 4.*atan(1.)
    f = sqrt(1./(8.*atan(1.))) /exp(t*t)
    end function f
    
