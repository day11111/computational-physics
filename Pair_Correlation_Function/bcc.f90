program bcc

    implicit none
    real :: i = 0.0, j = 0.0, k = 0.0, distance, dr = 0.1, r = 0.1, gr
    integer :: count = 0

    do while(r <= 10.1)
        ! count the points in simple cubic
        do while(i <= 20.0)
            do while ( j <= 20.0 ) 
                do while(k <= 20.0)
                    if(distance(i, j, k) >= r .and. distance(i, j, k) < (r+dr)) then 
                        ! print *, "i, j, k = ", i, j, k
                        ! print *, "distance = ", distance(i, j, k)
                        count = count + 1
                    end if
                k = k + 1.
                end do
            k = 0.0
            j = j + 1.
            end do
        j = 0.0
        i = i + 1.
        end do

        ! count the points in body center of the cubic
        i = 0.5
        j = 0.5
        k = 0.5
        do while(i <= 19.5)
            do while ( j <= 19.5 ) 
                do while(k <= 19.5)
                    if(distance(i, j, k) >= r .and. distance(i, j, k) < (r+dr)) then 
                        ! print *, "i, j, k = ", i, j, k
                        ! print *, "distance = ", distance(i, j, k)
                        count = count + 1
                    end if
                k = k + 1.
                end do
            k = 0.5
            j = j + 1.
            end do
        j = 0.5
        i = i + 1.
        end do

        gr = count/(4*(4*atan(1.)))/r**2
        ! print *, " count is ", count
        ! print *, "r = ", r
        ! print *, " gr = ", gr

        open(6, file = 'bcc.dat', position = 'append')
        write(6,'(5(E12.6,1X))') r, gr
        close(6)

        i = 0.0
        count = 0
        r = r + 0.1
    end do

end program bcc

function distance(i1, j1, k1)
    implicit none
    real, intent(in) :: i1, j1, k1
    real :: distance
    distance = sqrt((i1 - 10.0)**2 + (j1 - 10.0)**2 + (k1 - 10.0)**2)
    ! print *, "distance = ", distance
end function distance