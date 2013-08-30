!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! onespec.f90
!
! Simple NEXAFS spectrum generation code for ONETEP. Uses no tricks in an attempt to make
! the whole thing a bit more reliable! 
!
! 
!
! Written by Kane O'Donnell (Australian Synchrotron), September 2013.
!
! Contact the author on gmail, username kane dot odonnell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright 2013 Kane O'Donnell
!
!     This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Notes:
!
! 1. Usage is just onespec DATFILE PROJECTOR MINBAND. DATFILE is the txt output with
! the momentum matrix elements, PROJECTOR is the core projector to use as an initial 
! state and MINBAND (optional) is the minimum band to be included in the spectrum -  
! effectively allows manual tuning of the Fermi level.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program onespec

  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: DEBUG = 1
  character(len=80),parameter :: ENDIAN = 'big_endian' ! CASTEP is conventionally compiled
  ! as a big-endian code: if you get ridiculous results from nexspec check if changing
  ! this to 'native' fixes it.
  real(kind=dp), parameter :: invpi = 0.3183098861837907d0
  real(kind=dp), parameter :: pi = 3.141592653589793d0
  real(kind=dp), parameter :: invsqrt2pi = 0.3989422804014327d0
  real(kind=dp), parameter :: invfinestruct = 137.035999074d0
  real(kind=dp), parameter :: hart2eV = 27.211396132d0
  real(kind=dp), parameter :: tol = 1.0d-99
  
  real(kind=dp) :: smear_width ! eV
  integer, parameter :: spectrum_points = 2000
  real(kind=dp) :: w_start = -25.0 ! eV
  real(kind=dp) :: w_end = 65.0 ! eV
  real(kind=dp) :: w_step, w(spectrum_points)
  real(kind=dp) :: gwidth = 0.4 ! eV
  real(kind=dp) :: lwidth = 0.2 ! eV
  
  ! File stuff
  logical :: file_exists = .false.
  character(len=80) :: datfile, tmpstr, atstr(4)
  integer :: nargs, orb, nb, ix, iw, icmpt, iminband, cproj
  
  integer :: ncproj, mbands, nkpts, nspins, tmpi, nbands(2)
  real(kind=dp), dimension(:,:,:), allocatable :: optmat
  real(kind=dp), dimension(:,:), allocatable :: spectrum
  real(kind=dp), dimension(:), allocatable :: eigen, transitions
  real(kind=dp) :: linebits(6), cell_volume
  real(kind=dp) :: matrix_cmpt(6), e_nks, f_nks, smear_factor, e_min
  
  print *, "ONESPEC version 1"
  print *, " "
  print *, "Written by Kane O'Donnell, September 2013"
  print *, " "

  ! Get SEED from args
  nargs = command_argument_count()
  if (nargs .eq. 3) then
    call get_command_argument(1,datfile)
    call get_command_argument(2,tmpstr)
    read(tmpstr, *) cproj
    call get_command_argument(3, tmpstr)
    read(tmpstr, *), iminband
  else  
    print *, " "
    print *, "USAGE: onespec DATFILE PROJECTOR MINBAND"
    print *, " "
    print *, "DATFILE is the ONETEP EELS_MAT_ELS file. PROJECTOR is the core initial state index."
    print *, "MINBAND is the minimum band to be included in the spectrum, can be "
    print *, "used to effectively tune the Fermi level."
    call exit(0)
  endif

  ! Inquire if the relevant files exist - if not, die.
  inquire(file=adjustl(trim(datfile)), exist=file_exists)
  if (file_exists) then
    open(unit=100, file=adjustl(trim(datfile)), form='formatted',&
  &     status='old')
  endif

  ! Read everything we need from the EELS_MAT_ELS file.

!  read(100) ncproj
!  read(100) mbands
!  read(100) nkpts
!  read(100) nspins   

  ! First line is a comment
  read(100,*) tmpstr
  ! Number of core orbitals, number of bands
  read(100, '(I5,2x,I5)') ncproj, mbands, tmpi
  ! Cell volume (unused)
  read(100, '(F24.12)') cell_volume
  
  ! Check we specified something sensible for cproj and iminband
  if ((cproj .lt. 1) .or. (cproj .gt. ncproj)) then
    print *, "ERROR: You've specified a projector out of the range 1 - ", ncproj
    call exit(0)
  endif
  
  if (iminband .gt. mbands) then
    print *, "ERROR: You've specified a minimum band outside the range 1 - ", mbands
    call exit(0)
  endif
  
  if (iminband .eq. 1) then
    print *, "WARNING: You've specified iminband = 1 - this is almost certainly a non-physical spectrum!"
  end if

  ! Echo some output.
  
  write(*,*) "Opened file ", adjustl(trim(datfile)), " for NEXAFS spectrum generation."
  write(*,*) "Found ", ncproj, " core projectors."
  write(*,*) "Maximum number of bands per kpt/spin = ", mbands

  ! optmat is the optical transitions matrix, eigen is the eigenvalues
  allocate(optmat(ncproj,mbands,3))
  allocate(eigen(mbands), transitions(mbands))

  ! Loop over cartesian directions
  do ix=1,3
    read(100,*) tmpstr
    do orb=1,ncproj
      do nb=1,mbands
        read(100,*) linebits(:)
        eigen(nb) = linebits(4)
        optmat(orb,nb,ix) = linebits(5)
        transitions(nb) = linebits(6)
      end do
    end do
  end do

  close(100)

  
  ! Right, now we can generate the spectrum. We don't worry *at all* about occupancies
  ! or where the fermi level really is, etc - we just subtract the fermi level from
  ! all eigenvalues and the end user of the spectra can apply a Fermi function afterwards
  ! along with the proper smearing.
  
  ! Establish the spectrum: in order to make this compatible with the CASTEP output,
  ! we simply use -25 eV to 65 eV as with nexspec.
  allocate(spectrum(spectrum_points,6)) ! Separate spins
  
  spectrum = 0.0d0

  w_step = (w_end - w_start) / spectrum_points
  do iw=1,spectrum_points
    w(iw) = w_start + (iw-1)*w_step
  end do
  smear_width = w_step * 0.5
  
  ! Minimum eigenvalue (will shift this to zero and everything else down).
  e_min = eigen(iminband)
  
  ! Loop over bands, kpoints and spins - fix on a particular cproj. Note we can
  ! use OpenMP here if we want to.
  orb = cproj
!$OMP PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(ns, nk, nb, e_nks, f_nks, iw, w, smear_factor, icmpt, tmpstr) &
!$OMP& PRIVATE(matrix_cmpt)
  do nb=iminband,mbands
    
    matrix_cmpt(1) = optmat(orb,nb,1) * optmat(orb,nb,1)
    matrix_cmpt(2) = optmat(orb,nb,2) * optmat(orb,nb,2)
    matrix_cmpt(3) = optmat(orb,nb,3) * optmat(orb,nb,3)
    matrix_cmpt(4) = 2.0d0 * optmat(orb,nb,1) * optmat(orb,nb,2)
    matrix_cmpt(5) = 2.0d0 * optmat(orb,nb,1) * optmat(orb,nb,3)
    matrix_cmpt(6) = 2.0d0 * optmat(orb,nb,2) * optmat(orb,nb,3)

    e_nks = eigen(nb) - e_min
    
    print *, matrix_cmpt(1), optmat(orb,nb,1)
    ! Occupancy weighting: we don't deal with this at all anymore - f_nks is always 1.
    f_nks = 1.0d0
    
    do iw=1,spectrum_points
      ! Tiny lorentzian smear to project onto the nearest spectrum point.
      smear_factor = invpi * smear_width / ((e_nks - w(iw))**2 + (smear_width)**2)
      do icmpt=1,6
        spectrum(iw,icmpt) = spectrum(iw,icmpt) + f_nks * matrix_cmpt(icmpt) * smear_factor
      end do
    end do
  end do
!$OMP END PARALLEL DO
    
  ! Zero tiny values to make sure the output can be read.
  do iw=1,spectrum_points
    do icmpt=1,6
      if (dabs(spectrum(iw,icmpt)) .lt. tol) then
        spectrum(iw,icmpt) = 0.0d0
      end if
    end do
  end do

! Write spectra in the following files:
!   SEED_spinX.raw.nexafs for X=1->ns,
!   SEED.raw.nexafs - both spin channels combined.
!   SEED_spinX.smeared.nexafs for X=1->ns,
!   SEED.smeared.nexafs - smeared with both channels combined.

  write(*,*) 'Writing raw spectrum for orbital ', orb
  write(atstr(1), '(I4)') orb
  tmpstr = '_orb'//trim(adjustl(atstr(1)))
  open(300,file=trim(datfile)//trim(adjustl(tmpstr))//'.raw.nexafs',form='formatted')

  write(300,*) '# NEXAFS core-level spectrum calculated by onespec with ONETEP inputs.'
  write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'

  do iw=1,spectrum_points
    write(300, '(7g16.8)') w(iw), spectrum(iw,1:6)
  end do

  close(300)

  ! Write out a nominally-smeared (0.2 Lorentzian, 0.4 eV Gaussian, no linear) spectrum
  

  do icmpt=1,6
    call lorentzian_convolute(w,spectrum(:,icmpt),spectrum_points, lwidth)
    call gaussian_convolute(w,spectrum(:,icmpt),spectrum_points, gwidth)
  end do

  write(*,*) 'Writing smeared spectrum for orbital ', orb
  write(atstr(1), '(I4)') orb
  tmpstr = '_orb'//trim(adjustl(atstr(1)))
  open(300,file=trim(datfile)//trim(adjustl(tmpstr))//'.smeared.nexafs',form='formatted')

  write(300,*) '# NEXAFS core-level spectrum calculated by onespec with ONETEP inputs.'
  write(300,*) '# Smeared with lorentzian and gaussian broadening,', lwidth, gwidth, 'eV respectively'
  write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'

  do iw=1,spectrum_points
    write(300, '(7g16.8)') w(iw), spectrum(iw,1:6)
  end do

  close(300)
  
  ! Deallocate everything
  deallocate(spectrum)
  deallocate(optmat)
  deallocate(eigen)
  deallocate(transitions)
  
  write(*,*) "Finished onespec - Goodbye!"
  
  contains 
  
  subroutine gaussian_convolute(xdata, ydata, num_points, gwidth)
    !---------------------------------------------------------------------------
    !
    ! gaussian_convolute: a conventional convolution where the gaussian function
    ! of width gwidth is the second factor in the integral.
    !
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: num_points
    real(kind=dp), intent(inout) :: xdata(num_points), ydata(num_points)
    real(kind=dp) :: tmp(num_points), dx
    real(kind=dp), intent(in) :: gwidth
    integer :: i, j
    
    ! Assume the step is even (!!)
    dx = xdata(2) - xdata(1)
    
    tmp = 0.0d0
    
    do i=1,num_points
      do j=1,num_points
        tmp(i) = tmp(i) + dx * ydata(j) * invsqrt2pi * 1.0d0 / gwidth * dexp( &
        &        -0.5d0 * ((xdata(i) - xdata(j)) / gwidth)**2)
      end do
    end do
    
    ydata = tmp
    
    return
  end subroutine
  
  subroutine lorentzian_convolute(xdata, ydata, num_points, lwidth)
    !---------------------------------------------------------------------------
    !
    ! lorentzian_convolute: a conventional convolution where the lorentzian function
    ! of width gwidth is the second factor in the integral.
    !
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: num_points
    real(kind=dp), intent(inout) :: xdata(num_points), ydata(num_points)
    real(kind=dp) :: tmp(num_points), dx
    real(kind=dp), intent(in) :: lwidth
    integer :: i, j
    
    ! Assume the step is even (!!)
    dx = xdata(2) - xdata(1)
    
    tmp = 0.0d0
    
    do i=1,num_points
      do j=1,num_points
        tmp(i) = tmp(i) + dx * ydata(j) * invpi * (lwidth / ((xdata(i) - xdata(j))**2 + &
        & lwidth**2))
      end do
    end do
    
    ydata = tmp
    
    return
  end subroutine
  
  subroutine lorentzian_linear_convolute(xdata, ydata, num_points, lwidth, llin)
    !---------------------------------------------------------------------------
    !
    ! lorentzian_linear_convolute: a conventional convolution where the lorentzian function
    ! of width gwidth is the second factor in the integral. The linear broadening is 
    ! applied as lwidth + llin * | x(j) |, where x(j) is the x value of the spectrum
    ! being convoluted.  
    !
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: num_points
    real(kind=dp), intent(inout) :: xdata(num_points), ydata(num_points)
    real(kind=dp) :: tmp(num_points), dx
    real(kind=dp), intent(in) :: lwidth, llin
    integer :: i, j
    
    ! Assume the step is even (!!)
    dx = xdata(2) - xdata(1)
    
    tmp = 0.0d0
    
    do i=1,num_points
      do j=1,num_points
        tmp(i) = tmp(i) + dx * ydata(j) * invpi * ((lwidth + llin * dabs( &
        & xdata(j)))/ ((xdata(i) - xdata(j))**2 + (lwidth + llin * dabs( &
        & xdata(j)))**2))
      end do
    end do
    
    ydata = tmp
    
    return
  end subroutine
  
end program

