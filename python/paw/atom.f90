
module atom

  implicit none
  
  double precision  :: nuclear_charge, num_electrons
  integer           :: num_states
  double precision, allocatable, dimension(:) :: eigenvalues, occupancies
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
    if (.not. allocated(grid)) allocate(grid(grid_size))
    rp = rmax / (exp(grid_size * delta) - 1)    
    do j=1,grid_size
      grid(j) = rp * (exp(j * delta) - 1)
    end do
    
    ! Quantum numbers
    n = 1
    l = 0
    !m = 0 ! We aren't dealing with m-numbers in this code yet.
    if (.not. allocated(quantum_numbers)) allocate(quantum_numbers(num_states,2))
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
        F(i,j) = (rp ** 2) * (delta ** 2) * exp(2 * j * delta) * &
        & numerov_f(quantum_numbers(i,2), grid(j), eigenvalues(i)) +  &
        & 0.25d0 * (delta ** 2)
      end do
    end do
  
  end subroutine
  
  subroutine numerov(state_index, direction)
  
    integer :: j, state_index
    character(len=8) :: direction
    double precision, allocatable, dimension(:) :: A, B, C
    
    allocate(A(grid_size))
    allocate(B(grid_size))
    allocate(C(grid_size))
    
    do j=1,grid_size
      A(j) = 2.0d0 + (5.0d0 / 6.0d0) * F(state_index, j)
      B(j) = 1.0d0 - (1.0d0 / 12.0d0) * F(state_index, j)
      C(j) = 1.0d0 / B(j)
    end do
    
    select case(trim(direction))
    case('forward')
      v(state_index,1) = 0.0d0
      v(state_index,2) = 1.0d0
      do j=1,(grid_size-2)
        v(state_index,j+2) = (v(state_index,j+1) * A(j+1) - v(state_index,j) * B(j)) * C(j+2)
      end do
    case('backward')
      v(state_index, grid_size) = 0.0d0
      v(state_index, grid_size-1) = 1.0d0
      do j=grid_size,3,-1
        v(state_index,j-2) = (v(state_index,j-1) * A(j-1) - v(state_index,j) * B(j)) * C(j-2)
      end do
    end select
    
    do j=1,grid_size
      u(state_index,j) = v(state_index,j) * exp(j * delta)
    end do
    
    deallocate(A)
    deallocate(B)
    deallocate(C)
  
  end subroutine
  
  subroutine numerov2(state_index, direction)
  
    integer :: j, state_index
    character(len=8) :: direction
    
    select case(trim(direction))
    case('forward')
      v(state_index,1) = 0.0d0
      v(state_index,2) = 1.0d0
      q0 = v(state_index,
      do j=2,grid_size
        q2 = F(state_index, j-1) * v(state_index, j-1) + 2.0d0 * q1 - q0
        q0 = q1
        q1 = q2
        v(state_index, j) = q1 / (1.0d0 - F(state_index, j))
      end do
    case('backward')
      v(state_index, grid_size) = 0.0d0
      v(state_index, grid_size-1) = 1.0d0
      do j=grid_size,3,-1
        v(state_index,j-2) = (v(state_index,j-1) * A(j-1) - v(state_index,j) * B(j)) * C(j-2)
      end do
    end select
    
    do j=1,grid_size
      u(state_index,j) = v(state_index,j) * exp(j * delta)
    end do
    
    deallocate(A)
    deallocate(B)
    deallocate(C)
  
  end subroutine   
      
  double precision function numerov_f(l, r, E)

    integer, intent(in) :: l
    double precision, intent(in) :: r, E
    !double precision :: potential
    
    numerov_f = 2.0d0 * potential(r) + l * (l + 1) / (r ** 2) - 2.0d0 * E
  
  end function
  
  double precision function potential(r)
  
    double precision, intent(in) :: r
    
    ! Just the nuclear potential for now!
    potential = -1.0 * nuclear_charge / r
    
  end function
  
end module
  
  
