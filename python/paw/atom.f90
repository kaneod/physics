
module atom

  implicit none
  
  double precision, parameter :: pi = 3.141592653589793
  
  double precision  :: nuclear_charge, num_electrons
  integer           :: num_states
  double precision, allocatable, dimension(:) :: eigenvalues, occupancies, &
  & density, hartree
  double precision, allocatable, dimension(:,:) :: F, v, u
  integer, allocatable, dimension(:,:)           :: quantum_numbers
  
  ! Grid
  integer           :: grid_size
  double precision  :: rmax, delta, rp
  double precision, allocatable, dimension(:) :: grid
  
  ! parameters
  double precision, parameter :: REALLY_BIG=1.0d6
  double precision, parameter :: REALLY_SMALL=1.0d-10
  
  contains
  
  subroutine init(nuclear_charge_in, charge, num_extra_states, grid_size_in, rmax_in, delta_in)
  
    integer :: i, j, n, l, m, num_extra_states, grid_size_in
    double precision :: nuclear_charge_in, charge, electrons_remaining, rmax_in, delta_in
    
    ! Use inputs
    nuclear_charge = nuclear_charge_in
    num_electrons = nuclear_charge - charge
    num_states = ceiling(num_electrons / 2) + num_extra_states
    grid_size = grid_size_in
    rmax = rmax_in
    delta = delta_in
    
    ! Init grid
    if (allocated(grid)) deallocate(grid)
    allocate(grid(grid_size))
    rp = rmax / (exp(grid_size * delta) - 1)    
    do j=1,grid_size
      grid(j) = rp * (exp(j * delta) - 1)
    end do
    
    ! Quantum numbers
    n = 1
    l = 0
    !m = 0 ! We aren't dealing with m-numbers in this code yet.
    if (allocated(quantum_numbers)) deallocate(quantum_numbers)
    allocate(quantum_numbers(num_states,2))
    do i=1,num_states
      quantum_numbers(i,1) = n
      quantum_numbers(i,2) = l
      !quantum_numbers(i,3) = m
      !if (m .eq. l) then
        if (l .eq. n - 1) then
          n = n + 1
          l = 0
          !m = 0
        else
          l = l + 1
          !m = -l
        end if
      !else
      !  m = m + 1
      !end if
    end do  
    
    ! Allocate and guess eigenvalues and occupancies.
    if (allocated(eigenvalues)) deallocate(eigenvalues)
    allocate(eigenvalues(num_states))
    if (allocated(occupancies)) deallocate(occupancies)
    allocate(occupancies(num_states))    
    electrons_remaining = num_electrons
    do i=1,num_states
      eigenvalues(i) = -0.5d0 * nuclear_charge ** 2 / quantum_numbers(i,1) ** 2
      if (electrons_remaining .gt. 2.0d0) then
        occupancies(i) = 2.0d0
        electrons_remaining = electrons_remaining - 2.0d0
      else
        occupancies(i) = electrons_remaining
        electrons_remaining = 0.0d0
      end if
    end do
    
    ! Allocate F and state vectors
    if (allocated(F)) deallocate(F)
    allocate(F(num_states, grid_size))
    if (allocated(v)) deallocate(v)
    allocate(v(num_states, grid_size))
    if (allocated(u)) deallocate(u)
    allocate(u(num_states, grid_size))
    
    do i=1,num_states
      do j=1,grid_size
        u(i,j) = 0.0d0
        v(i,j) = 0.0d0
        F(i,j) = (rp ** 2) * (delta ** 2) * exp(2 * j * delta) * &
        & schrod_f(quantum_numbers(i,2), grid(j), eigenvalues(i)) +  &
        & 0.25d0 * (delta ** 2)
      end do
    end do
    
    ! Generate the initial states
    do i=1,num_states
      call numerov_schrod(i,"backwards")
    end do
    
    ! Generate the initial density
    call compute_density
    
    
  
  end subroutine
  
  subroutine compute_density
  
    integer :: i, j
    
    if (allocated(density)) deallocate(density)
    allocate(density(grid_size))
    
    ! Ensure density is zero
    do j=1,grid_size
      density(j) = 0.0d0
    end do
    
    do i=1,num_states
      do j=1,grid_size
        density(j) = density(j) + occupancies(i) * u(i,j) * u(i,j)
      end do
    end do
    
  end subroutine
  
  subroutine integrate_poisson
  
    integer :: j, k
    double precision :: integrand(grid_size)
    
    if (allocated(hartree)) deallocate(hartree)
    allocate(hartree(grid_size))
    
    do j=1,grid_size
      do k=1,grid_size
        if (j .ne. k) then
          integrand(k) = density(k) / abs(grid(j) - grid(k))
        else
          ! Eliminate the infinite interactions by brute force
          integrand(k) = 0.0d0
        end if
      end do
      
      hartree(j) = 4.0d0 * pi * integrate(grid, integrand, grid_size)
    
    end do
    
  end subroutine
    
  
  subroutine verlet_poisson
  
    integer :: j
    double precision :: v2, alpha, h0
    
    if (allocated(hartree)) deallocate(hartree)
    allocate(hartree(grid_size))
    
    ! Boundary conditions are U(0) = 0, U(rmax) = num_electrons. 
    ! Can't set the second initially so we guess U(1) = 1 and scale later.
    
    h0 = 0.0d0
    hartree(1) = 1.0d0
    
    do j=2,grid_size
      v2 = -1.0d0 * rp * (delta ** 2) * exp(1.5d0 * (j - 1) * delta) * &
      & density(j-1) / (exp((j - 1) * delta) - 1.0d0) + 0.25d0 * (delta ** 2) &
      & * hartree(j-1)
      hartree(j) = 2.0d0 * hartree(j-1) - h0 + v2
      h0 = hartree(j-1)
    end do
    
    ! Rescale back to r coordinates
    do j=1,grid_size
      hartree(j) = hartree(j) * exp(0.5d0 * j * delta)
    end do
    
    ! Now sort out boundary conditions
    
    alpha = (num_electrons - hartree(grid_size)) / rmax
    
    do j=1,grid_size
      hartree(j) = hartree(j) + alpha * grid(j)
    end do
    
  end subroutine
  
  subroutine numerov_schrod(state_index, direction)
  
    integer :: j, state_index
    character(len=8) :: direction
    double precision :: q0, q1, q2, unorm
    double precision :: u2(grid_size)
    
    select case(trim(direction))
    case('forward')
      v(state_index,1) = grid(1) * exp(-0.5d0 * delta)
      v(state_index,2) = grid(2) * exp(-1.0d0 * delta)
      q0 = 0
      q1 = v(state_index,1) * (1.0d0 - F(state_index,1) / 12.0d0)
      do j=2,grid_size
        q2 = F(state_index,j-1) * v(state_index,j-1) + 2.0d0 * q1 - q0
        q0 = q1
        q1 = q2
        v(state_index, j) = q1 / (1.0d0 - F(state_index, j))
      end do
    case('backward')
      v(state_index,grid_size) = exp(-1.0d0 * sqrt(-2.0d0 * &
      & eigenvalues(state_index)) * grid(grid_size)) * exp(-0.5d0 * & 
      & grid_size * delta)
      v(state_index,grid_size-1) = exp(-1.0d0 * sqrt(-2.0d0 * &
      & eigenvalues(state_index)) * grid(grid_size-1)) * exp(-0.5d0 * & 
      & grid_size * delta)
      q0 = v(state_index,grid_size) * (1.0d0 - F(state_index,grid_size) / &
      & 12.0d0)
      q1 = v(state_index,grid_size-1) * (1.0d0 - F(state_index,grid_size-1) / &
      & 12.0d0)
      do j=grid_size-2,1,-1
        q2 = F(state_index,j+1) * v(state_index,j+1) + 2.0d0 * q1 - q0
        q0 = q1
        q1 = q2
        v(state_index, j) = q1 / (1.0d0 - F(state_index, j))
      end do
    end select
    
    ! Conversion from j-space to r-space.
    
    do j=1,grid_size
      u(state_index,j) = v(state_index,j) * exp(0.5d0 * j * delta)
    end do
    
    ! Normalization
    
    do j=1,grid_size
      u2(j) = u(state_index,j) * u(state_index,j)
    end do
    
    unorm = integrate(grid, u2, grid_size)
    u(state_index, 1:grid_size) = u(state_index, 1:grid_size) / sqrt(unorm)
  
  end subroutine   
      
  double precision function schrod_f(l, r, E)

    integer, intent(in) :: l
    double precision, intent(in) :: r, E
    !double precision :: potential
    
    schrod_f = 2.0d0 * potential(r) + l * (l + 1) / (r ** 2) - 2.0d0 * E
  
  end function
  
  double precision function potential(r)
  
    double precision, intent(in) :: r
    
    ! Just the nuclear potential for now!
    potential = -1.0 * nuclear_charge / r
    
  end function
  
  double precision function integrate(x, y, n)
    
    integer :: n, j
    double precision :: x(n), y(n)
    
    integrate = 0.0d0
      
    do j=1,n-1
      integrate = integrate + 0.5 * (x(j+1) - x(j)) * (y(j+1) + y(j))
    end do
    
  end function
  
end module
  
  
