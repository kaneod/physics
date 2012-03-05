
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
  
  
  
  contains
  
  subroutine init(nuclear_charge_in, charge, num_states_in, grid_size_in, rmax_in, delta_in)
  
    integer :: i, j, n, l, m, num_states_in, grid_size_in
    double precision :: nuclear_charge_in, charge, electrons_remaining, rmax_in, delta_in
    
    ! Use inputs
    num_states = num_states_in
    nuclear_charge = nuclear_charge_in
    num_electrons = nuclear_charge - charge
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
    m = 0
    if (.not. allocated(quantum_numbers)) allocate(quantum_numbers(num_states,3))
    do i=1,num_states
      quantum_numbers(i,1) = n
      quantum_numbers(i,2) = l
      quantum_numbers(i,3) = m
      if (l .eq. n - 1) then
        ! We have reached the maximum of the l-cycle: increment n.
        n = n + 1
        l = 0
        m = 0
      else if (m .eq. l) then 
        ! We have reached the maximum of the m-cycle: increment l.
        l = l + 1
        m = -l
      else
        ! Increment m.
        m = m + 1
      end if
    end do  
    
    ! Allocate and guess eigenvalues and occupancies.
    if (.not. allocated(eigenvalues)) allocate(eigenvalues(num_states))
    if (.not. allocated(occupancies)) allocate(occupancies(num_states))    
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
    if (.not. allocated(F)) allocate(F(num_states, grid_size))
    if (.not. allocated(v)) allocate(v(num_states, grid_size))
    if (.not. allocated(u)) allocate(u(num_states, grid_size))
    
    do i=1,num_states
      do j=1,grid_size
        F(i,j) = numerov_f(quantum_numbers(i,2), grid(j), eigenvalues(i))
      end do
    end do
  
  end subroutine
      
  double precision function numerov_f(l, r, E)

    integer, intent(in) :: l
    double precision, intent(in) :: r, E
    double precision :: potential
    
    numerov_f = 2.0d0 * potential(r) + l * (l + 1) / (r ** 2) - 2.0d0 * E
  
  end function
  
  double precision function potential(r)
  
    double precision, intent(in) :: r
    
    ! Just the nuclear potential for now!
    potential = -1.0 * nuclear_charge / r
    
  end function
  
end module
  
  