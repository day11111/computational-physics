! Name： Wang Jiping  ID: 11510049
! This program calculates the surface energy for given k in Pb(111) and Ag(111)

program fermi

    implicit none
    real :: k = 1.0, N = 1.0, S, EsPb, EsAg
    real, parameter :: kFpb = 1.57, EFpb = 9.37, kfAg = 1.20, EFAg = 5.48

    do while(k <= 20.0)
        N = 1.0   ! initial N

        ! calculate N for infinite barrier
        do while(N**2/(k**2) > (2*k/(3*N) + (N+1)*(2*N+1)/(6*k**2)) &
            .or. (2*k/(3*N) + (N+1)*(2*N+1)/(6*k**2)) >= (N+1)**2/(k**2))
            N = N + 1.0
        end do

        S = 4*k/(9*N) + 2*(N+1)*(2*N+1)/(9*k**2) &
            - N*(N+1)*(2*N+1)*(8*N**2+3*N-11)/(180*k**5)
       
        EsPb = kFpb**2*EFpb/(16*atan(1.))*k*(S - 4./5.)
        EsAg = kfAg**2*EFAg/(16*atan(1.))*k*(S - 4./5.)

        open(6,FILE='fermi.dat', position='append')
        write(6,'(5(E12.6,1X))') k, EsPb, EsAg
        close(6)

        k = k + 0.1
    end do

end program fermi
