program monteCarlo
implicit none

    ! ----------------------------------------------------------------------------------
    ! Calculating the normal probability function integral mc
    ! ----------------------------------------------------------------------------------

    !the number of partitions is n
    !a and b are the beginning point and the ending point of the integral
    !the function is even, so let a = 0.0
    !we only calculate integral from 0.0 to 1.0

    integer :: n = 10,  i, j, m
    real :: harvest, sum =0, f
    do j = 1, 5
        m = n ** j
        do i = 1, m
            call random_number(harvest)
            sum = sum + f(harvest)
        end do 
        sum = 2*sum / m
        print *, "Monte Carlo : ", "n = ", m ,"the integral is : ", sum
    end do

end program monteCarlo

    function f (t)
    implicit none
    real, intent(in) :: t
    real :: f
    ! pi = 4.*atan(1.)
    f = sqrt(1./(8.*atan(1.))) /exp(t*t)
    end function f
