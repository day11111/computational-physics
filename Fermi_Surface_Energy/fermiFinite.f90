! This program is designed to calculate the surface energy for given k in Pb(111) and Ag(111)
! For now this program can calculate Es-k for Pb(111). For Ag(111), three variables should change values.
! V0Pb = 1.45 * 9.37 -> V0Pb = 1.78 * 5.48, kFPb = 1.57 -> kFPb = 1.20, EFPb = 9.37-> EFPb = 5.48 
! in line 11, 13 and 83.
!
! This program is coded by Jinlong Huang. Copyright (C) 2016. All right reserved.

program fermiFinite

    implicit none
    real :: k = 1.0, N = 1.0, S, EsPb = 0., EsAg = 0., yita, V0Pb = 1.45 * 9.37, fPb, yita_N
    real :: a = -1.0, b = 1.0, yita_1, yita_2, dx, dl = 5.0e-6, i, Efk, En
    real, parameter :: kFPb = 1.57, EFPb = 9.37

    do while(k <= 20.0)
        ! initial N
        N = 1.0     
        yita_1 = 0.
        yita_2 = 0.
        Efk = 0.0

        ! calculate N for finite barrier
        do while(yita_1**2*EFPb/k**2 > Efk .or. Efk >= yita_2**2*EFPb/k**2)
            ! initial
            Efk = 0.0
            i = 1.0

            ! E(N) = yita_1**2*EFPb/k**2
            yita_1 = yita_N(k, N, dl, V0Pb)
            ! print *, "---------------------------------------yita_1 is ", yita_1
        
            ! E(N+1) = yita_2**2*EFPb/k**2
            yita_2 = yita_N(k, (N+1.), dl, V0Pb)
            ! print *, "---------------------------------------yita_2 is ", yita_2

            ! Efk
            do while(i <= N)
                Efk = Efk + yita_N(k, i, dl, V0Pb)**2*EFPb/k**2
                i = i + 1.
                ! print *, "i = ", i
            end do
            Efk = Efk/N + EFPb*2*k/(3*N)
            ! print *, "---------------------------------------Efk is ", Efk

            if(yita_1**2*EFPb/k**2 > Efk .or. Efk >= yita_2**2*EFPb/k**2) then
                N = N + 1.0
                ! print *, "N =  ", N
            end if

            ! print *, "k is ", k
        end do 
        ! end calculate N for finite barrier

        ! EsPb
        i = 1.0
        Espb = 0.
        ! print *, "EsPb!"
        ! print *, "N = ", N
        ! print *, "Efk = ", Efk, "E(N) = ", yita_1**2*EFPb/k**2, "E(N+1) = ", yita_2**2*EFPb/k**2
        do while(i <= N)
            En =  yita_N(k, i, dl, V0Pb)**2*EFPb/k**2
            ! print *, "En = ", En
            EsPb = EsPb + (Efk**2-En**2)/EFPb**2
            ! print *, "EsPb = ", EsPb
            i = i + 1.
        end do
        EsPb = (kFPb**2)*EFPb/(16.*atan(1.)) * (EsPb-4.*k/5.) 
        ! print *, "EsPb is ", EsPb

        open(6,FILE='fermiFinite.dat', position='append')
        write(6,'(5(E12.6,1X))') k, EsPb
        close(6)

        k = k + 0.01
        ! print *, "k = ", k
    end do
end program fermiFinite

! yita_N function for Pb
! only called by yita_N
function fPb(yita, k, N, V0Pb)   
    implicit none
    real :: fPb, EFPb = 9.37
    real, intent(in) :: yita, k, N, V0Pb
    fPb = yita - N + 2.0/(4.0*atan(1.))*asin(yita/(k*sqrt(V0Pb/EFPb)))
end function fPb

! calculate the zero point of yita_N using bisection method
function yita_N(k, N, dl, V0Pb)
    real :: yita_N
    real, intent(in) :: k, N, dl, V0Pb
    real :: a, b, dx
    ! initial
    a = N - 1.0
    b = N + 1.0
    dx = b - a
    do while(abs(dx) > dl)
        yita_N = (a + b)/2.0
        if(fPb(a, k, N, V0Pb)*fPb(yita_N, k, N, V0Pb) >= 0) then
            a = yita_N
            dx = b - a
            ! print *, "dx = ", dx
        else 
            b = yita_N
            dx = b - a
            ! print *, "dx = ", dx
        end if
    end do
    ! print *, "yita_N = ", yita_N
end function yita_N