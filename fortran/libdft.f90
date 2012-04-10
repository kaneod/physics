!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! libdft.f90
!
! Library of modules for doing things like Fermi's golden rule calculations
! via the generation of optical matrix elements and so on.
!
! Designed to be compilable with f2py (leading to a python-compatible library.)
!
! Written by Kane O'Donnell (Australian Synchrotron), Apr 2012.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright 2012 Kane O'Donnell
!
!     This module is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This module is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this module.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Notes:
!
! 1. We always use atomic units where necessary. We assume this has been taken
! care of before entry into the modules.
!
! 2. f2py currently fails to parse selected_real_kind and (kind=8) seems
! to have a bit of a problem with gfortran46 in some cases so we just
! use double precision and double complex.
! 
! 3. Each module requires an initialize call from Python before it has the data
! to be able to do anything else - this is because it's easier to read files 
! and parse things in Python.
!
! 4. The library has to be linked against fftw3.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module optical

  implicit none
 
  ! FFTW parameters
  
  integer, parameter :: FFTW_FORWARD=-1
  integer, parameter :: FFTW_BACKWARD=1
  integer, parameter :: FFTW_ESTIMATE=64
  integer*8 :: plan
  
  ! Grid dimensions. These must stay constant between initializations!
  
  integer :: nx, ny, nz
  
  ! Initial and final wavefunctions plus a work area for things like gradients.
  
  double complex, allocatable, dimension(:,:,:) :: wfi, wff, work
  
  ! G-vectors and coefficient indices
  
  double precision, allocatable, dimension(:,:) :: G
  integer, allocatable, dimension(:) :: gcx, gcy, gcz
  
  ! Lattice and reciprocal lattice vectors
  
  double precision :: avec(3,3), bvec(3,3) ! Convention is "vectors in columns"
  
  ! The current k-point for wff.
  
  double precision :: kpt(3)
  
  ! Current overlap integral value (if calculated).
  
  double complex :: ovlp
  
  ! Current optical matrix element (if calculated).
  
  double precision :: ome
  
  contains
  
  subroutine test_init(L)
  
    ! Subroutine solely for testing: calls init with test data.
    
    integer :: L
    double complex :: twf(L, L, L)
    double precision :: lvec(3,3)
    
    twf = 0.0d0
    
    lvec(1,1) = 1
    lvec(2,1) = 0
    lvec(3,1) = 0
    lvec(1,2) = 0
    lvec(2,2) = 1
    lvec(3,2) = 0
    lvec(1,3) = 0
    lvec(2,3) = 0
    lvec(3,3) = 1
    
    call init(twf, lvec, lvec, L, L, L)
    
  end subroutine
  
  subroutine init(core_wf, lattice, recip_lattice, L, M, N)
    
    integer :: i,j
    integer, intent(in) :: L, M, N
    double complex, intent(in) :: core_wf(L, M, N)
    double precision, intent(in) :: lattice(3,3), recip_lattice(3,3)
    
    if (allocated(wff)) deallocate(wff)
    if (allocated(wfi)) deallocate(wfi)
    if (allocated(work)) deallocate(work)
    allocate(wfi(L, M, N))
    allocate(wff(L, M, N))
    allocate(work(L, M, N))
    
    ! These are assumed set and unchanged between inits.
    nx = L
    ny = M
    nz = N
    
    do j=1,3
      do i=1,3
        avec(i,j) = lattice(i,j)
        bvec(i,j) = recip_lattice(i,j)
      end do
    end do
    
    call generate_g_vectors
    call wf_copy(core_wf, wfi)
    call wf_copy(core_wf, wff)
    
    ! For now just set the working k-point to gamma.
    kpt(1) = 0.0d0
    kpt(2) = 0.0d0
    kpt(3) = 0.0d0
    
  end subroutine
  
  subroutine set_wff(new_wff, new_kpt, L, M, N)
  
    !! We need to be careful here because we get segfaults when the arrays are
    !! large. So, explicitly declare size of new wf.
    
    integer, intent(in) :: L, M, N
    double precision, intent(in) :: new_kpt(3)
    double complex, intent(in) :: new_wff(L, M, N)
  
    kpt(:) = new_kpt(:)
    
    call wf_copy(new_wff, wff)
    
  end subroutine
  
  subroutine wf_copy(a, b)
  
    ! Copy wf a onto wf b (don't worry about allocation
    ! and assume dimensions are nx,ny, nz
    
    integer :: i,j,k
    double complex, dimension(:,:,:) :: a, b
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
          b(i,j,k) = a(i,j,k)
        end do
      end do
    end do
    
  end subroutine
  
  subroutine generate_g_vectors
  
    integer :: i, j, k, p
    double precision :: tmp_vector(3)
    
    ! We just have a list of g-vectors rather than a 4D array,
    ! as there seem to be MALLOC problems with 4D arrays at the moment
    ! and it's easy enough to transfer between the two.
    
    if (allocated(gcx)) deallocate(gcx)
    if (allocated(gcy)) deallocate(gcy)
    if (allocated(gcz)) deallocate(gcz)
    if (allocated(G)) deallocate(G)
    
    allocate(gcx(nx), gcy(ny), gcz(nz), G(nx*ny*nz,3))
    
    ! G-vector coordinates
    do i=1,(nx/2)+1
      gcx(i) = i - 1
    end do
    
    do i=1,(ny/2)+1
      gcy(i) = i - 1
    end do
    
    do i=1,(nz/2)+1
      gcz(i) = i - 1
    end do
    
    do i=(nx/2)+2,nx
      gcx(i) = i - 1 - nx
    end do
    
    do i=(ny/2)+2,ny
      gcy(i) = i - 1 - ny
    end do
    
    do i=(nz/2)+2,nz
      gcz(i) = i - 1 - nz
    end do
    
    p = 1
    do k=1,nz
      !tmp_vector = gcz(k) * bvec(:,3)
      do j=1,ny
        !tmp_vector = tmp_vector + gcy(j) * bvec(:,2)
        do i=1,nx
          G(p,:) = gcx(i) * bvec(:,1) + gcy(j) * bvec(:,2) + gcz(k) * bvec(:,3)
          p = p + 1
        end do
      end do
    end do
  
  end subroutine
  
  double complex function integrate_grids(a, b)
  
    ! Numerical, extremely naive triple integral over grid.
    ! Note this is a QM-style integral where the integrand
    ! is conjugate(a) * b.
    
    integer :: i, j, k
    double complex, dimension(:,:,:) :: a,b
    
    ! As usual assume a, b have shape (nx,ny,nz)
    
    integrate_grids = (0.0d0, 0.0d0)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          integrate_grids = integrate_grids + conjg(a(i,j,k)) * b(i,j,k)
        end do
      end do
    end do
    
    integrate_grids = integrate_grids / (nx * ny * nz)
    
    return
  
  end function
  
  subroutine wff_gradient(axis)
  
    ! Compute grad(psi)_axis for a wavefunction loaded to wff.
    ! Note: use set_wff to set the wavefunction and the working
    ! kpt. axis can be 1, 2 or 3.
  
    integer :: p, i, j, k
    integer, intent(in) :: axis
    
    ! Copy wff to work and then in place FFT to recip-space.
    
    call wf_copy(wff, work)
    
    ! FFTW plan
    
    call dfftw_plan_dft_3d(plan,nx, ny, nz, work, work, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, work, work)
    call dfftw_destroy_plan(plan)
    
    ! work now contains wff_g (the fft of wff). Now multiply by (K+G).
    
    p = 1
    do k=1,nz
      do j=1,ny
        do i=1,nx
          work(i,j,k) = (kpt(axis) + G(p,axis)) * work(i,j,k)
          p = p + 1
        end do
      end do
    end do
    
    ! Now FFT back (have to normalize by the number of points due to
    ! FFTW not normalizing the FFT) to real space.
    
    call dfftw_plan_dft_3d(plan, nx, ny, nz, work, work, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, work, work)
    call dfftw_destroy_plan(plan)
    
    work = work / (nx * ny * nz)
  
  end subroutine
  
  subroutine fft_work(direction)
  
    integer, intent(in) :: direction ! Can be FFTW_FORWARD (-1) or
                                     ! FFTW_BACKWARD (1). No checks!
    
    call dfftw_plan_dft_3d(plan,nx, ny, nz, work, work, direction, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, work, work)
    call dfftw_destroy_plan(plan)
    
    ! If the direction is backwards, re-normalize.
    
    if (direction .eq. FFTW_BACKWARD) then
      work = work / (nx * ny * nz)
    end if
    
  end subroutine
  
  subroutine optical_matrix_element(axis, out_ome)
  
    integer, intent(in) :: axis
    double precision, intent(out) :: out_ome
  
    ! Calculates overlap integral <wfi|grad|wff>
    
    ! First ensure grad|wff> is in work.
    call wff_gradient(axis)
    
    ! Now integrate <wfi|work>
    ovlp = integrate_grids(wfi, work)
    
    ! The OME itself is the norm squared.
    ome = realpart(conjg(ovlp) * ovlp)
    
    ! Set out_ome to give an output if desired.
    out_ome = ome
  
  end subroutine
      
end module

module paw

  !! Module for calculating properties inside PAW spheres instead of over
  ! linear grids. Not sure if this will work but if it does, should be
  ! much much (much!) faster.

  implicit none
  
  ! These variables store the x,y,z indices and the wf value for linear grid
  ! points only inside a particular PAW sphere. Used to calculate linear grid
  ! integrals over the sphere (things like <phi_i | psi_nk>).
  
  double precision :: paw_sphere_radius
  double precision :: sphere_centre(3)
  integer :: N ! Number of points inside PAW sphere
  integer, allocatable, dimension(:) :: px, py, pz
  double precision, allocatable, dimension(:) :: rx, ry, rz
  double complex, allocatable, dimension(:) :: paw_wf
  
  contains
  
  subroutine sphere_points_locate(grid_size, grid_n, r0, paw_radius)
  
    double precision :: grid_size(3), r0(3), paw_radius
    double precision :: gv(3), cell_size(3), gn
    integer :: grid_n(3), i, j, k, p
    
    ! Use the passed parameters to configure module variables.
    paw_sphere_radius = paw_radius
    sphere_centre(1:3) = r0(1:3)
    cell_size(1) = grid_size(1) / grid_n(1)
    cell_size(2) = grid_size(2) / grid_n(2)
    cell_size(3) = grid_size(3) / grid_n(3)
    
    ! Loop over grid indices, generate grid points and figure out if they are
    ! in the sphere or not.
    p = 0
    do k=1,grid_n(3)
      gv(3) = (k - 1) * cell_size(3)
      do j=1,grid_n(2)
        gv(2) = (j - 1) * cell_size(2)
        do i=1,grid_n(1)
          gv(1) = (i - 1) * cell_size(1)
          gn = sqrt((gv(1)-r0(1)) ** 2 + (gv(2)-r0(2)) ** 2 + &
          &        (gv(3)-r0(3)) ** 2)
          
          if (gn .le. paw_radius) then
            p = p + 1
          end if
        end do
      end do
    end do
          
    ! Now we know how many points we have, we can allocate the arrays and
    ! go through the loop again to add the points.
    
    if (allocated(px)) deallocate(px)
    if (allocated(py)) deallocate(py)
    if (allocated(pz)) deallocate(pz)
    if (allocated(rx)) deallocate(rx)
    if (allocated(ry)) deallocate(ry)
    if (allocated(rz)) deallocate(rz)
    if (allocated(paw_wf)) deallocate(paw_wf)
    
    allocate(px(p), py(p), pz(p), rx(p), ry(p), rz(p), paw_wf(p))
    N = p
    
    p = 1
    do k=1,grid_n(3)
      gv(3) = (k - 1) * cell_size(3)
      do j=1,grid_n(2)
        gv(2) = (j - 1) * cell_size(2)
        do i=1,grid_n(1)
          gv(1) = (i - 1) * cell_size(1)
          gn = sqrt((gv(1)-r0(1)) ** 2 + (gv(2)-r0(2)) ** 2 + &
          &        (gv(3)-r0(3)) ** 2)
          
          if (gn .le. paw_radius) then
            px(p) = i
            py(p) = j
            pz(p) = k
            rx(p) = gv(1)
            ry(p) = gv(2)
            rz(p) = gv(3)
            p = p + 1
          end if
        end do
      end do
    end do
   
    print *, "Number of grid points: ", grid_n(1) * grid_n(2) * grid_n(3)
    print *, "Number of points inside sphere:", N
    print *, "Reduction factor: ", N / (grid_n(1) * grid_n(2) * grid_n(3))
    
  end subroutine
  
end module
          
    
  
