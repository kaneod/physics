
module multipole
  
  implicit none
  
  !!!!!!!!!!
  !
  ! External functions
  !
  !!!!!!!!!!
  
  double precision, external :: ddot
  
  !!!!!!!!!!
  !
  ! Pre-computed constants
  !
  !!!!!!!!!!
  
  double precision :: pi, y0, y1, y2, y3, y4, y5, y6, y7
  
  !!!!!!!!!!
  !
  ! Counters, tmp variables, etc
  !
  !!!!!!!!!!
  
  logical :: grid_set = .false.  ! If .true., we have already set up everything for a given
                        ! grid size and instead of init_sphere, should use reinit_sphere
  logical :: projector_set = .false. ! Similarly for the projector matrix P (=YtW here).
  
  !!!!!!!!!!
  !
  ! Grid data variables
  !
  !!!!!!!!!!
  
  double precision :: lvec(3,3), h
  
  !!!!!!!!!!
  !
  ! Sphere variables
  !
  !!!!!!!!!!
  
  integer :: ns ! Number of points inside sphere
  double precision :: rs, delta ! test sphere radius
  double precision :: r0(3)
  double precision, allocatable, dimension(:) :: S, rp
  double precision, allocatable, dimension(:,:) :: rpts
  
  !!!!!!!!!!
  !
  ! Basis variables
  !
  !!!!!!!!!!
  
  integer :: nmax, lmax, nb ! Max n and l for basis vectors, number of ind. basis vecs.
  integer, allocatable, dimension(:,:) :: lm
  double precision, allocatable, dimension(:) :: F,flm ! Coefficients matrix
  double precision, allocatable, dimension(:,:) :: G, W, Y ! Metric, weighting, design matrices
  double precision, allocatable, dimension(:,:) :: YtW
  
  contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! clear_grid
  !
  ! Deallocates and unsets all data contained by this module.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine clear_grid
  
    deallocate(F, flm, lm, G, W, Y, YtW, S, rp, rpts)
    grid_set = .false.
    projector_set = .false.
  
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! compute_multipole_coefficients
  !
  ! Main subroutine for this multipole module.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_multipole_coefficients(grid, gridx, gridy, gridz, avec, r_sph, r_c, l_max)
  
    integer, intent(in) :: gridx, gridy, gridz, l_max
    double precision, intent(in) :: grid(gridx, gridy, gridz), avec(3,3), r_sph, r_c(3)
    double precision :: tmp_vec(3)
    
    if (.not. grid_set) then
      print *, "Have not set up for this grid yet: initializing."
      call compute_constants
    
      lvec(1:3,1) = avec(1:3,1) / gridx
      lvec(1:3,2) = avec(1:3,2) / gridy
      lvec(1:3,3) = avec(1:3,3) / gridz
      
      tmp_vec(1) = lvec(2,2) * lvec(3,3) - lvec(3,2) * lvec(2,3)
      tmp_vec(2) = lvec(3,2) * lvec(1,3) - lvec(1,2) * lvec(3,3)
      tmp_vec(3) = lvec(1,2) * lvec(2,3) - lvec(2,2) * lvec(1,3)
      
      h = (ddot(3, lvec(1:3,1), 1, tmp_vec, 1)) ** (1.0d0 / 3)
    
      r0 = r_c
      rs = r_sph
      lmax = l_max
      
      call init_sphere(grid, gridx, gridy, gridz)
      call init_weights
      call init_lms
      
      grid_set = .true.
      
      call compute_coefficients
      
    else
      print *, "Using previously computed grid information."
      r0 = r_c
      rs = r_sph
      call reload_sphere(grid, gridx, gridy, gridz)
      call compute_coefficients
    end if
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! compute_coefficients
  !
  ! Construct Y, YtW and G, then solve GP = YtW and finally apply F = PS.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_coefficients
  
    integer :: i, j, k, p, q
    double precision :: cc ! Coefficient of interest
    
    ! If we have set up the grid before, most of this can be skipped.
    
    if (.not. (projector_set .and. grid_set)) then
      allocate(Y(ns,nb), YtW(nb,ns), G(nb, nb), F(nb))
    
      
      do i=1,ns
        do j=1,nb
          Y(i,j) = bisquare(lm(j,1), rp(i)) * ylm(rpts(i,:), lm(j,2), lm(j,3))
        end do
      end do
    
      call dgemm("T", "N", nb, ns, ns, 1.0d0, Y, ns, W, ns, 0.0d0, YtW, nb)
      call dgemm("N", "N", nb, nb, ns, 1.0d0, YtW, nb, Y, ns, 0.0d0, G, nb)

      ! Compute P as solution of GP = YtW. Note that P is stored in YtW.
      call dpotrf("L", nb, G, nb, i)
      call dpotrs("L", nb, ns, G, nb, YtW, nb, i)
      
      projector_set = .true.
      
    end if
    
    ! Get our derived coefficients F
    call dgemv("N", nb, ns, 1.0d0, YtW, nb, S, 1, 0.0d0, F, 1)

    ! Now take our linear combinations to get the lm coefficients  
    p = 1
    q = 1
    do i=0,lmax
      do j=-i,i
        cc = 0.0d0
        do k=0,nmax
          cc = cc + F(p)
          p = p + 1
        end do
        flm(q) = cc
        q = q + 1
      end (do
    end do

  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! init_lms
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_lms
  
    integer :: n,l,m,p
    
    nmax = 2
    
    ! Count nlm combinations. Note the order (l,m,n) is 
    ! chosen here to make subsequent linear combinations
    ! easier to find because we want to sum over n.
    p=0
    do l=0,lmax
      do m=-l,l
        do n=0,nmax
          p=p+1
        end do
      end do
    end do
    
    nb = p
    allocate(lm(p,3), flm(nb / (nmax+1)))
    
    p=1
    do l=0,lmax
      do m=-l,l
        do n=0,nmax
          lm(p,1) = n
          lm(p,2) = l
          lm(p,3) = m
          p = p + 1
        end do
      end do
    end do
    
  end subroutine
    
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! init_weights
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_weights
  
    integer :: p
    double precision :: u
    
    allocate(W(ns, ns))
    
    W = 0.0

    do p=1,ns
      u = (rp(p) - rs) / delta
      if (abs(u) .lt. 1.0) then
        W(p,p) = (1.0 - u ** 2) ** 2
      else
        W(p,p) = 0.0
      end if
    end do
    
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! init_sphere
  !
  ! Finds points in the sphere and fills S.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_sphere(grid, gridx, gridy, gridz)
    
    integer, intent(in) :: gridx, gridy, gridz
    double precision, intent(in) :: grid(gridx, gridy, gridz)
    integer :: i, j, k, p
    double precision :: r, rv(3), gv(3)
    
    delta = 1.25 * h
    
    p = 0
    do k=1,gridz
      do j=1,gridy
        do i=1,gridx
          gv = (i-1) * lvec(:,1) + (j-1) * lvec(:,2) + (k-1) * lvec(:,3)
          rv = gv - r0
          r = sqrt(rv(1) ** 2 + rv(2) ** 2 + rv(3) ** 2)
          if (abs(r - rs) .lt. delta) then
            p = p + 1
          end if
        end do
      end do
    end do
    
    allocate(S(p), rp(p),rpts(p,3))
    ns = p
    
    p = 1
    do k=1,gridz
      do j=1,gridy
        do i=1,gridx
          gv = (i-1) * lvec(:,1) + (j-1) * lvec(:,2) + (k-1) * lvec(:,3)
          rv = gv - r0
          r = sqrt(rv(1) ** 2 + rv(2) ** 2 + rv(3) ** 2)
          if (abs(r - rs) .lt. delta) then
            S(p) = grid(i,j,k)
            rp(p) = r
            rpts(p,1) = rv(1)
            rpts(p,2) = rv(2)
            rpts(p,3) = rv(3)
            p = p + 1
          end if
        end do
      end do
    end do 
    
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! reload_sphere
  !
  ! If we are just reloading a new function over the same grid we can skip a lot
  ! of the init_sphere function.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine reload_sphere(grid, gridx, gridy, gridz)
    
    integer, intent(in) :: gridx, gridy, gridz
    double precision, intent(in) :: grid(gridx, gridy, gridz)
    integer :: i, j, k, p
    double precision :: r, rv(3), gv(3)
    
    p = 1
    do k=1,gridz
      do j=1,gridy
        do i=1,gridx
          gv = (i-1) * lvec(:,1) + (j-1) * lvec(:,2) + (k-1) * lvec(:,3)
          rv = gv - r0
          r = sqrt(rv(1) ** 2 + rv(2) ** 2 + rv(3) ** 2)
          if (abs(r - rs) .lt. delta) then
            S(p) = grid(i,j,k)
            rp(p) = r
            rpts(p,1) = rv(1)
            rpts(p,2) = rv(2)
            rpts(p,3) = rv(3)
            p = p + 1
          end if
        end do
      end do
    end do 
  
  end subroutine
    
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! init_test
  !
  ! Sets up the grid data for testing.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine init_test(n, a)
!  
!    integer, intent(in) :: n   ! Number of division of the grid sides.
!    double precision, intent(in) :: a ! lattice parameter of the box.
!    integer :: i, j, k
!    double precision :: rv(3), gv(3)
!    
!    gridx = n
!    gridy = n
!    gridz = n
!    
!    r0(1) = a/2
!    r0(2) = a/2
!    r0(3) = a/2
!    
!    h = a / n
!    
!    allocate(grid(n,n,n))
!    
!    avec(:,1) = (/ a, 0.0d0, 0.0d0 /)
!    avec(:,2) = (/ 0.0d0, a, 0.0d0 /)
!    avec(:,3) = (/ 0.0d0, 0.0d0, a /)
!    
!    lvec = avec / n
!    
!    c0 = 1.0
!    c1 = 0.5
!    c2 = 0.33
!    
!    ! Put a test function on the grid
!    do k=1,n
!      do j=1,n
!        do i=1,n
!          gv = (i-1) * lvec(:,1) + (j-1) * lvec(:,2) + (k-1) * lvec(:,3)
!          rv = gv - r0
!          grid(i,j,k) = c0 * ylm(rv, 0, 0) + c1 * ylm(rv, 1, -1) + c2 * ylm(rv, 2,2)
!        end do
!      end do
!    end do
!    
!  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! compute_constants
  !
  ! Calculates the prefactors for the Ylms.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine compute_constants
    
    pi = 4.0 * atan(1.0)
    y0 = sqrt(1.0/pi)
    y1 = sqrt(3.0 / (4.0 * pi))
    y2 = sqrt(5.0 / pi)
    y3 = sqrt(15.0 / pi)
    y4 = sqrt(7.0 / pi)
    y5 = sqrt(17.5 / pi)
    y6 = sqrt(105.0 / pi)
    y7 = sqrt(10.5 / pi)
    
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! bisquare(n, x)
  !
  ! Returns bn(x), polynomials orthonormal
  ! with respect to bisquare weighting integration.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  double precision function bisquare(n, x)
  
    double precision, intent(in) :: x
    integer, intent(in) :: n
    
    select case (n)
      case (0)
        bisquare = sqrt(15.0 / 16)
      case (1)
        bisquare = sqrt(105.0 / 16) * x
      case (2)
        bisquare = sqrt(45.0 / 64) * (1.0 - 7.0 * x ** 2)
      case (3)
        bisquare = sqrt(1155.0 / 64) * (x - 3.0 * x ** 3)
      case (4)
        bisquare = sqrt(1365.0 / 2048) * (1.0 - 18.0 * x ** 2 + 33.0 * x ** 4)
    end select
    return
    
  end function
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! legendre(n, x)
  !
  ! Returns Pn(x), the value of the nth
  ! legendre polynommial at x.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  double precision function legendre(n, x)
  
    double precision, intent(in) :: x
    integer, intent(in) :: n
    
    select case (n)
      case (0)
        legendre = 1.0
      case (1)
        legendre = x
      case (2)
        legendre = 0.5 * (3.0 * x ** 2 - 1)
      case (3)
        legendre = 0.5 * (5.0 * x ** 3 - 3.0 * x)
      case (4)
        legendre = 0.125 * (35.0 * x ** 4 - 30.0 * x ** 2 + 3)
      case (5)
        legendre = 0.125 * (63.0 * x ** 5 - 70.0 * x ** 3 + 15.0 * x)
    end select
    
    return
  
  end function
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! ylm(r, l, m)
  !
  ! Returns the value of Ylm as a function of the vector r.
  !
  ! r should be a 3-component vector relative to some origin.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  double precision function ylm(r, l, m)
    
    double precision, intent(in) :: r(3)
    integer, intent(in) :: l, m
    double precision :: x, y, z, rho
    
    ! Get components
    x = r(1)
    y = r(2)
    z = r(3)
    rho = sqrt(x ** 2 + y ** 2 + z ** 2)
    
    ! If rho == 0, we are undefined. Return 0 if l > 0, 0.5*y0 otherwise.
    if (rho .lt. 1.0d-8) then
      if (l .gt. 0) then
        ylm = 0.0
      else
        ylm = 0.5 * y0
      end if
      return
    end if
    
    select case (l)
      case (0)
        ylm = 0.5 * y0
        return
      case (1)
        select case (m)
          case (-1)
            ylm = y1 * x / rho
          case (0)
            ylm = y1 * y / rho
          case (1)
            ylm = y1 * z / rho
        end select
        return
      case (2)
        select case (m)
          case (-2)
            ylm = 0.25 * y2 * (2.0 * z ** 2 - x ** 2 - y ** 2) / rho ** 2
          case (-1)
            ylm = 0.5 * y3 * y * z / rho ** 2
          case (0)
            ylm = 0.5 * y3 * z * x / rho ** 2
          case (1)
            ylm = 0.5 * y3 * x * y / rho ** 2
          case (2)
            ylm = 0.25 * y3 * (x ** 2 - y ** 2) / rho ** 2
        end select
        return
      case (3)
        select case (m)
          case (-3)
            ylm = 0.25 * y4 * z * (2.0 * z ** 2 - 3.0 * x ** 2 - 3.0 * y ** 2) / rho ** 3
          case (-2)
            ylm = 0.25 * y5 * (3.0 * x ** 2 - y ** 2) * y / rho ** 3
          case (-1)
            ylm = 0.25 * y5 * (x ** 2 - 3.0 * y ** 2) * x / rho ** 3
          case (0)
            ylm = 0.5 * y6 * (x ** 2 - y ** 2) * z / rho ** 3
          case (1)
            ylm = 0.5 * y6 * x * y * z / rho ** 3
          case (2)
            ylm = 0.25 * y7 * y * (4.0 * z ** 2 - x ** 2 - y ** 2) / rho ** 3
          case (3)
            ylm = 0.25 * y7 * x * (4.0 * z ** 2 - x ** 2 - y ** 2) / rho ** 3
        end select
        return
    end select
    
  end function
  
end module