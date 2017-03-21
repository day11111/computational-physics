MODULE ISING
	implicit none
    integer             :: N_L ! number of lattice points in x and y
    integer             :: N_D ! total number of lattice sites
    integer             :: steps, metro, nt, magnetization
    integer,allocatable :: spin(:,:)
    real*8              :: Pflip(-4:4) ! spin flip probability
    real*8              :: temperature, tmax, dt ! (in Joules)
    real*8              :: average_magnetization

CONTAINS

SUBROUTINE initialize
    integer :: i,j
    real*8  :: r

	open(1,file='ising_2d_input.dat')
	read(1,*) nt
	read(1,*) tmax
	read(1,*) dt
	read(1,*) steps
	close(1)

    N_L = 20
	N_d=N_L**2
	metro = N_d
	Pflip=0.d0
	temperature = tmax

	do i=-4,4,2  ! do i = -4, 4
		Pflip(i)=exp(-i*2.d0/temperature)
        ! print *, "i = ", i
        ! print *, "Pflip = ", Pflip(i)
	end do

	allocate(spin(0:N_L-1,0:N_L-1))

	! Randomly arrange spins
	call random_seed()
	do j=0,N_L-1
		do i=0,N_L-1
			call random_number(r)
			spin(i,j)=2*int(2.d0*r)-1  ! int(3.d0*r)-1
            ! print *, "i, j, spin(i,j) = ", i, j, spin(i,j)
		end do
	end do
	average_magnetization=0.d0
END SUBROUTINE initialize

SUBROUTINE spin_flip(t)
    ! Flip spins using the Metropolis-algorithm
    integer :: i,k,site,x,y,spin_left
    integer :: spin_right,spin_above,spin_below,spin_sum
    real    :: xi
    integer, intent(in) :: t

    ! Choose a random site   
    call random_number(xi)
    site=int(xi*N_D)
    ! print *, "site = ", site

    ! x and y coordinate
    x=mod(site,N_L) 
    y=site/N_L
    ! print *, "x = ",x, "y = ", y, "spin(x,y) = ", spin(x,y)

    ! Four neighbors
    spin_right=spin(mod(x+1,N_L),y)
    spin_left=spin(mod(x-1+N_L,N_L),y)
    spin_above=spin(x,mod(y+1,N_L))
    spin_below=spin(x,mod(y-1+N_L,N_L))
    spin_sum=spin_right+spin_left+spin_above+spin_below
    ! print *, "spin_sum = ", spin_sum

    ! temperature
    temperature = tmax - dt * t
    ! print *, "temperature = ", temperature

    ! flip probability
    ! the parameter k in Pflip(k) stands for dE/2
    do k = -4,4,2  ! do k = -4, 4
        Pflip(k)=exp(-k*2.d0/temperature)
        ! print *, "Pflip = ", Pflip(k)
    end do

    ! Metropolis 
    call random_number(xi)
    ! print *, "cut-off probability xi = ", xi
    ! spin(x,y)*spin_sum = the energy after flip = dE/2
    if(Pflip(spin(x,y)*spin_sum) >= xi) then
        spin(x,y)=-spin(x,y)
        ! spin(x,y) = mod(spin(x,y)+2, 3) - 1
        ! print *, "yes"
    else 
        ! print *, "no"
    end if
END SUBROUTINE spin_flip

SUBROUTINE calculate_properties
    integer :: x, y

    magnetization=0.d0
    do y=0,N_L-1
        do x=0,N_L-1
            magnetization = magnetization + spin(x,y)
        end do
    end do
    ! print *, "magnetization = ", magnetization

END SUBROUTINE calculate_properties

SUBROUTINE output
    average_magnetization = dble(magnetization)/(dble(steps)*dble(metro)*dble(N_d))
    ! print *, "average_magnetization = ", average_magnetization
END SUBROUTINE output

END MODULE ISING

PROGRAM ising2d
    USE ISING
    implicit none
    integer :: i, j, t

    call initialize

    ! Thermalize
    t = 0
    do i=0, metro-1
        call spin_flip(t)
    end do

    ! From tmax to tmin
    do t=0, nt-1
    
        ! steps is the number of Monte Carlo
        do j=0, steps-1
            ! metro is the number of Metropolis
            do i=0, metro-1
                call spin_flip(t)
                call calculate_properties
            end do
        end do

        call output
        open(6, file = 'ising_2d_output.dat', position = 'append')
        write(6,'(5(E12.6,1X))') temperature, average_magnetization
        close(6)
    end do

END PROGRAM ising2d

