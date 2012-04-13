program multipole

  implicit none
  
  double precision, external :: ddot
  external :: zgemm, dpotrf, zpotrs, zgemv
  
  ! Counters etc
  integer :: i, p, j, k
  double precision :: vec(3)
  double complex :: fsum
  
  ! Timing information
  integer :: clock_start, clock_end, clock_rate
  double precision :: elapsed_time
  
  ! 3D grid
  integer :: gridx=50, gridy=50, gridz=50
  double complex :: grid(50,50,50)
  
  ! Lattice vectors
  double precision :: avec(3,3)
  ! Central position vector
  double precision :: r0(3)
  
  ! Points inside computational shell and the max number of basis functions
  ! M is the number of independent {n,l,m} combinations
  integer :: N, nmax=3, lmax=3, M
  integer, allocatable, dimension(:,:) :: lm ! Set of l,m combinations
  integer, allocatable, dimension(:) :: ipiv
  double complex, allocatable, dimension(:) :: S, F
  ! Matrices. Need work matrices. Gc is a complex version of G, likewise for Wc.
  double complex, allocatable, dimension(:,:) :: Y, Wk1, Wk2, Pr, Gc, W
  ! Real matrices
  double precision, allocatable, dimension(:,:) :: G, WkL1, WrkL2
  double precision, allocatable, dimension(:) :: rx, ry, rz, rp
  integer, allocatable, dimension(:) :: px, py, pz
  double precision :: rs, delta ! Radius of sphere, shell width
  
  ! Known and "unknown" coefficents for the test.
  double complex :: c0, c1, c2, c3, u(0:3)
  
  ! Precomputed values
  double precision :: pi, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, vol
  
  pi = 4.0d0 * atan(1.0d0)
  h1 = 0.5d0 * sqrt(1.0d0 / pi)
  h2 = 0.5d0 * sqrt(3.0d0 / pi)
  h3 = 0.5d0 * sqrt(3.0d0 / (2.0d0 * pi))
  h4 = 0.5d0 * sqrt(5.0d0) * h3
  h5 = 2.0d0 * h4
  h6 = 0.25d0 * sqrt(5.0d0 / pi)
  h7 = 0.125d0 * sqrt(35.0d0 / pi)
  h8 = 0.25d0 * sqrt(105.0d0 / (2.0d0 * pi))
  h9 = 0.125d0 * sqrt(21.0d0 / pi)
  h10 = 0.25d0 * sqrt(7.0d0 / pi)    
  
  c0 = (0.1, 0.9)
  c1 = (0.25, 0.75)
  c2 = (2.26, 0.11)
  c3 = (0.5, 0.5)
  
  grid = 0.0d0
  
  avec(:,1) = (/ 10.0d0, 0.0d0, 0.0d0 /)
  avec(:,2) = (/ 0.0d0, 10.0d0, 0.0d0 /)
  avec(:,3) = (/ 0.0d0, 0.0d0, 10.0d0 /)
  
  ! Grid cell volume
  vec(1) = avec(2,2) * avec(3,3) - avec(3,2) * avec(2,3)
  vec(2) = avec(1,3) * avec(3,2) - avec(1,2) * avec(3,3)
  vec(3) = avec(1,2) * avec(2,3) - avec(2,2) * avec(1,3)  
  vol = ddot(3, avec(:,1), 1, vec, 1) / (gridx * gridy * gridz)
  
  rs = 1.5
  ! Use the suggested value of Fiske (delta = 5/4 * h, or Delta=0.75h)
  delta = 1.25 * vol ** (real(1)/3)
  
  print *, "Cell volume: ", vol * gridx * gridy * gridz
  print *, "Grid cell volume: ", vol
  print *, "delta ", delta
  
  r0 = (/ 5.0d0, 5.0d0, 5.0d0 /)
  
  ! Construct a function with a known set of Ylm coefficients and put on the grid.
  
  call system_clock(count_rate=clock_rate)
  
  call system_clock(count=clock_start)
  call ylm_test
  call system_clock(count=clock_end)
  
  print *, "Elapsed time for ylm_test: ", & 
  &               real(clock_end - clock_start) / clock_rate
  
  ! Print out some test values for checking.
  
  print *, "Grid 25,1,25: ", grid(25,1,25)
  print *, "Grid 25,50,25: ", grid(25,50,25)
  print *, "Grid 1,25,25: ", grid(1,25,25)
  print *, "Grid 50,25,25: ", grid(50,25,25)
  
  ! Next trick is to determine which points lie within the computational shell.
  call system_clock(count_rate=clock_rate)
  
  call system_clock(count=clock_start)
  call shell_init
  call system_clock(count=clock_end)
  
  print *, "Elapsed time for shell_init: ", & 
  &               real(clock_end - clock_start) / clock_rate  
  
  ! Print out the points list
  print *, "Points in total:", gridx*gridy*gridz
  print *, "Found ", N, " points inside sphere."
  !print *, "p S(p) px(p) py(p) pz(p) rx(p) ry(p) rz(p) rp(p)"
  !do i=1,N
  !  print '(1i6, 2f5.2, 3i3, 4f5.2)', i, S(i), px(i), py(i), pz(i), rx(i), & 
  !  &                                 ry(i), rz(i), rp(i)
  !end do
  
  ! Print out linear combinations of entries in F
  print *, F
      
  contains
  
  subroutine ylm_test
    
    integer :: i,j,k
    double precision :: r(3), g(3), rho, theta, phi
    
    ! Note: in future, avoid function calls by using the x,y,z form of the Ylms.
    
    do k=1,gridz
      do j=1,gridy
        do i=1,gridx
          g = real(i-1)/gridx * avec(:,1) + real(j-1)/gridy * avec(:,2) + real(k-1)/gridz * avec(:,3)
          r = g - r0
          rho = sqrt(r(1) ** 2 + r(2) ** 2 + r(3) ** 2)
          if (rho .lt. 1.0d-8) then
            grid(i,j,k) = 0.0d0
          else
            grid(i,j,k) = c0 * ylm(r,0,0) + c1 * ylm(r,1,0) + c2 * ylm(r,2,0) + &
            &             c3 * ylm(r,3,0)
          end if
        end do
      end do
    end do
  
  end subroutine
  
  subroutine shell_init
    
    integer :: i, j, k, p
    double precision :: gv(3), r(3), rho
    
    ! Figure out how many lm components we need, generate and store.
    p = 0
    do i=0,lmax
      do j=-i,i
        do k=0,nmax
          p=p+1
        end do
      end do
    end do
    M = p
    
    if (allocated(lm)) deallocate(lm)
    allocate(lm(p,3))
    
    p = 1
    do i=0,lmax
      do j=-i,i
        do k=0,nmax
          lm(p,1) = k
          lm(p,2) = i
          lm(p,3) = j
          p = p + 1
        end do
      end do
    end do  
    
    p = 0
    do k=1,gridz
      do j=1,gridy
        do i=1,gridx
          gv = real(i-1)/gridx * avec(:,1) + real(j-1)/gridy * avec(:,2) + real(k-1)/gridz * avec(:,3)
          r = gv - r0
          rho = sqrt(r(1) ** 2 + r(2) ** 2 + r(3) ** 2)     
          if ((rho .ge. rs - delta) .and. (rho .le. rs + delta)) then
            p = p + 1
          end if
        end do
      end do
    end do
    
    N = p
    if (allocated(S)) deallocate(S)
    if (allocated(px)) deallocate(px)
    if (allocated(py)) deallocate(py)
    if (allocated(pz)) deallocate(pz)
    if (allocated(rx)) deallocate(rx)
    if (allocated(ry)) deallocate(ry)
    if (allocated(rz)) deallocate(rz)
    if (allocated(rp)) deallocate(rp)
    if (allocated(Y)) deallocate(Y)
    if (allocated(W)) deallocate(W)
    if (allocated(G)) deallocate(G)
    if (allocated(Gc)) deallocate(Gc)
    if (allocated(Wk1)) deallocate(Wk1)
    if (allocated(Wk2)) deallocate(Wk2)
    if (allocated(Pr)) deallocate(Pr)
    if (allocated(F)) deallocate(F)
    allocate(S(p), px(p), py(p), pz(p), rx(p), ry(p), rz(p), rp(p))
    allocate(Y(N,M), W(N,N), G(M,M), Wk1(N, M), F(M), Pr(M,N), Wk2(M,N))
    allocate(Gc(M,M))
    
    p = 1
    do k=1,gridz
      do j=1,gridy
        do i=1,gridx
          gv = real(i-1)/gridx * avec(:,1) + real(j-1)/gridy * avec(:,2) + real(k-1)/gridz * avec(:,3)
          r = gv - r0
          rho = sqrt(r(1) ** 2 + r(2) ** 2 + r(3) ** 2)          
          if ((rho .ge. (rs - delta)) .and. (rho .le. (rs + delta))) then
            S(p) = grid(i,j,k)
            px(p) = i
            py(p) = j
            pz(p) = k
            rx(p) = gv(1)
            ry(p) = gv(2)
            rz(p) = gv(3)
            rp(p) = rho
            p = p + 1
          end if
        end do
      end do
    end do
    
    ! Generate the weights matrix W
    W = 0.0d0
    do p=1,N
      if (abs(rp(p) - rs) .gt. delta) then
        W(p,p) = 0.0d0
      else if (abs(rp(p) - rs) .lt. (delta - (vol ** (real(1)/3)))) then
        W(p,p) = vol
      else
        W(p,p) = (delta - abs(rp(p) - rs)) * (vol ** (real(2)/3))
      end if
      if (real(W(p,p)) .lt. 0) then
        print *, "Warning: Negative weights at p, W(p,p) = ", p, W(p,p)
      end if
    end do
    
    ! Generate the Y matrix
    Y = 0.0d0
    do p=1,N
      do i=1,M
        Y(p,i) = rlegendre(lm(i,1), rp(p)) * ylm(rp(p), lm(i,2), lm(i,3))
        !print *, "p, i, lm(i,:),rp(p), rleg, ylm", p, i, lm(i,:), &
        !&         rp(p), rlegendre(lm(i,1), rp(p)), ylm(rp(p), lm(i,2), lm(i,3))
      end do
    end do
    
    ! Now generate YtW and Gc = YtWY
    Gc = 0.0d0
    Wk1 = 0.0d0
    Wk2 = 0.0d0
    ! This ZGEMM gives us YtW
    call zgemm("C", "N", M, N, N, 1.0d0, Y, N, W, N, 0.0d0, Wk2, M)
    ! This ZGEMM gives Gc = YtWY
    call zgemm("N", "N", M, M, N, 1.0d0, Wk2, M, Y, N, 0.0d0, Gc, M)
    
    deallocate(Wk1)
    
    ! Cholesky factorize G = L*Lt, we get G transformed into L on return.
    !call zpotrf("L", M, Gc, M, i)
    !print *, "ZPOTRF returned ", i
    
    ! Now solve LLtP = YtW
    ! Remember YtW = Wk2. Need a MxN workspace.
    allocate(Wk1(M,N))
    
    !call zpotrs("L", M, N, Gc, M, Wk2, M, i)
    allocate(ipiv(M))
    call zgesv(M, N, Gc, M, ipiv, Wk2, M, i)
    
    print *, "ZGESV returned ", i
    
    ! Finally, compute the coefficients vector F = PS
    ! P is stored in Wk2
    
    call zgemv("N", M, N, 1.0d0, Wk2, M, S, 1, 0.0d0, F, 1) 
    
  end subroutine
  
  double complex function ylm(r, l, m)
  
    double precision, intent(in) :: r(3)
    integer, intent(in) :: l, m
    
    double precision :: rho
  
    ! We use the cartesian ylms here to save function calls
    
    if (l .eq. 0) then
      ylm = h1
      return
    else 
      rho = sqrt(r(1) ** 2 + r(2) ** 2 + r(3) ** 2) 
      if (l .eq. 1) then  
        if (m .eq. -1) then
          ylm = h3 * complex(r(1), -r(2)) / rho
        else if (m .eq. 0) then
          ylm = h2 * r(3) / rho
        else if (m .eq. 1) then
          ylm = -h3 * complex(r(1), r(2)) / rho
        end if
        return
      else if (l .eq. 2) then
        if (m .eq. -2) then
          ylm = h4 * (complex(r(1), -r(2)) / rho) ** 2
        else if (m .eq. -1) then
          ylm = h5 * complex(r(1), -r(2)) * r(3) / (rho ** 2)
        else if (m .eq. 0) then
          ylm = h6 * (2.0d0 * r(3) ** 2 - r(1) ** 2 - r(2) ** 2) / rho ** 2
        else if (m .eq. 1) then
          ylm = -h5 * complex(r(1), r(2)) * r(3) / rho ** 2
        else if (m .eq. 2) then
          ylm = h4 * (complex(r(1), r(2)) / rho) ** 2
        end if
        return
      else if (l .eq. 3) then
        if (m .eq. -3) then
          ylm = h7 * (complex(r(1), -r(2)) / rho) ** 3
        else if (m .eq. -2) then
          ylm = h8 * complex(r(1), -r(2)) ** 2 * r(3) / rho ** 3
        else if (m .eq. -1) then
          ylm = h9 * complex(r(1), -r(2)) * (4.0d0 * r(3) ** 2 - r(1) ** 2 - r(2) ** 2) / &
          &     rho ** 3
        else if (m .eq. 0) then
          ylm = h10 * r(3) * (2.0d0 * r(3) ** 2 - 3.0d0 * r(2) ** 2 - 3.0d0 * r(2) ** 2) &
          &     / rho ** 3
        else if (m .eq. 1) then
          ylm = -h9 * complex(r(1), r(2)) * (4.0d0 * r(3) ** 2 - r(1) ** 2 - r(2) ** 2) / &
          &     rho ** 3
        else if (m .eq. 2) then
          ylm = h8 * complex(r(1), r(2)) ** 2 * r(3) / rho ** 3
        else if (m .eq. 3) then
          ylm = -h7 * (complex(r(1), r(2)) / rho) ** 3
        end if
        return
      end if
    end if
  end function 
  
  double precision function legendre(n, x)
  
    ! There may be a problem with the n=4 and 5 functions here.
    integer, intent(in) :: n
    double precision, intent(in) :: x
    
    if (n .eq. 0) then
      legendre = 1.0d0
      return
    else if (n .eq. 1) then
      legendre = x
      return
    else if (n .eq. 2) then
      legendre = 0.5d0 * (3.0d0 * x ** 2 - 1)
      return   
    else if (n .eq. 3) then
      legendre = 0.5d0 * (5.0d0 * x ** 3 - 3.0d0 * x)
      return
    else if (n .eq. 4) then
      legendre = 0.125d0 * (35.0d0 * x ** 4 - 30.0d0 * x ** 2 + 3)
      return
    else if (n .eq. 5) then
      legendre = 0.125d0 * (63.0d0 * x ** 5 - 70.0 * x ** 3 + 15.0d0 * x)
      return
    end if
  end function
             
  double precision function rlegendre(n, x)
  
    integer, intent(in) :: n
    double precision, intent(in) :: x
    
    rlegendre = (1.0d0 / x) * sqrt((2.0d0 * n + 1) / (2 * delta)) * legendre(n, (x-rs) / delta)
    return
  end function

end program  