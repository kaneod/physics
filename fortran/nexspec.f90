!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nexspec.f90
!
! Simple NEXAFS spectrum generation code. Uses no tricks in an attempt to make the whole
! thing a bit more reliable!
!
! 
!
! Written by Kane O'Donnell (Australian Synchrotron), August 2013.
!
! Contact the author on gmail, username kane dot odonnell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright 2012 Kane O'Donnell
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
! 1. Usage is just nexspec SEED PROJECTOR [MINBAND]. Assumes the existence of SEED.bands
! and SEED.eels_mat, PROJECTOR is the core projector to use as an initial state and 
! MINBAND (optional) is the minimum band to be included in the spectrum - effectively 
! allows manual tuning of the Fermi level.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program nexspec

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
  character(len=80) :: seed, tmpstr, atstr(4)
  integer :: nargs, nk, orb, ns, nb, iw, icmpt, cproj, iminband
  
  integer :: ncproj, mbands, nkpts, nspins, tmpi, nbands(2)
  integer, dimension(:), allocatable :: core_species, core_ion, core_n, core_lm
  complex(kind=dp), dimension(:,:,:,:,:), allocatable :: optmat
  real(kind=dp), dimension(:), allocatable :: wk
  real(kind=dp), dimension(:,:), allocatable :: kpts, cspectrum
  real(kind=dp), dimension(:,:,:), allocatable :: eigen, spectrum
  real(kind=dp) :: efermi(2), nelectrons(2), lvec(3,3)
  real(kind=dp) :: matrix_cmpt(6), e_nks, f_nks, smear_factor, e_min
  
  print *, "NEXSPEC version 2"
  print *, " "
  print *, "Written by Kane O'Donnell, August 2013"
  print *, " "

  ! Get SEED from args
  nargs = command_argument_count()
  if (nargs .eq. 2) then
    call get_command_argument(1,seed)
    call get_command_argument(2,tmpstr)
    read(tmpstr, *) cproj
    iminband = 0
  elseif (nargs .eq. 3) then
    call get_command_argument(1,seed)
    call get_command_argument(2,tmpstr)
    read(tmpstr, *) cproj
    call get_command_argument(3, tmpstr)
    read(tmpstr, *), iminband
  else  
    print *, " "
    print *, "USAGE: nexspec SEED PROJECTOR [MINBAND]"
    print *, " "
    print *, "Require SEED.bands and SEED.eels_mat. PROJECTOR is the core initial state index."
    print *, "MINBAND (optional) is the minimum band to be included in the spectrum, can be "
    print *, "used to effectively tune the Fermi level."
    call exit(0)
  endif

  ! Inquire if the relevant files exist - if not, die.
  inquire(file=adjustl(trim(seed))//'.eels_mat', exist=file_exists)
  if (file_exists) then
    open(unit=100, file=adjustl(trim(seed))//'.eels_mat', form='unformatted',&
  &     status='old', convert=ENDIAN)
  endif

  inquire(file=adjustl(trim(seed))//'.bands', exist=file_exists)
  if (file_exists) then
    open(unit=200, file=adjustl(trim(seed))//'.bands', form='formatted', &
  &       status='old')
  endif

  ! Read everything we need from the eels_mat file.

  read(100) ncproj
  read(100) mbands
  read(100) nkpts
  read(100) nspins   
  
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
  
  write(*,*) "Opened file ", adjustl(trim(seed))//'.eels_mat', " for NEXAFS spectrum generation."
  write(*,*) "Found ", ncproj, " core projectors."
  write(*,*) "Maximum number of bands per kpt/spin = ", mbands
  write(*,*) "Number of K-points = ", nkpts
  write(*,*) "Number of spins = ", nspins

  ! Allocate the various arrays that store information about the core orbitals.
  allocate(core_species(ncproj), core_ion(ncproj), core_n(ncproj), &
  &        core_lm(ncproj))

  ! optmat is the optical transitions matrix.
  allocate(optmat(ncproj,mbands,3,nkpts,nspins))

  ! Allocate wk (kpoints weighting), eigen (the eigenvalues) and kpts
  allocate(wk(nkpts), eigen(mbands,nkpts,nspins), kpts(nkpts,3))

  ! Read the information about the core orbitals.
  read(100) core_species(1:ncproj)
  read(100) core_ion(1:ncproj)
  read(100) core_n(1:ncproj)
  read(100) core_lm(1:ncproj)
  
  ! Print out a table of available orbitals.
  print *, "Orbital    ", "Species     ", "Ion         ", "n            ", "lm        "
  
  do orb=1,ncproj
    print *, orb, core_species(orb), core_ion(orb), core_n(orb), core_lm(orb)
  end do

  ! Loop over kpts, spins, orbitals, bands and populate optmat.
  do nk=1,nkpts
    do ns=1,nspins
      do orb=1,ncproj
        do nb=1,mbands
          read(100) optmat(orb, nb, 1:3, nk, ns)
        end do
      end do
    end do
  end do

  ! We're now done with the .eels_mat file, start on the .bands file next.
  close(100)
  
  ! Check that we have the same number of kpoints, spins, etc - there is a lot of 
  ! redundant information in the bands file but we might as well check it's the same.
  read(200,'(19x,I5)') tmpi

  if (tmpi .ne. nkpts) then
    write(*,*) "ERROR: Number of k-points in .bands is not equal to the number in .eels_mat."
    write(*,*) "ERROR: The module is now in an uncertain state: re-initialize before use."
    call exit(0)
  end if

  read(200,'(26x,I1)') tmpi
  if (tmpi .ne. nspins) then
    write(*,*) "ERROR: Number of spins in .bands is not equal to the number in .eels_mat."
    write(*,*) "ERROR: The module is now in an uncertain state: re-initialize before use."
    call exit(0)
  end if  

  ! Note: nelectrons here is TOTALLY UNRELIABLE but we read it anyway.
  if (DEBUG .eq. 1) then
    print *, "Made it past integer reads"
  end if
  if (nspins .eq. 1) then
    read(200,'(20x,g10.4)') nelectrons(1)
    read(200,'(22x,I6)') nbands(1)
    read(200,'(31x,F12.6)') efermi(1)
    write(*,*) "Inside the .bands file we get:"
    write(*,*) "Number of electrons (possibly wrong!) = ", nelectrons(1)
    write(*,*) "Number of bands = ", nbands(1)
    write(*,*) "Fermi level = ", efermi(1), " Ha"
  else
    read(200,'(20x,2g10.4)') nelectrons(1:2)
    read(200,'(22x,2I6)') nbands(1:2)
    read(200,'(33x,2F12.6)') efermi(1:2)
    write(*,*) "Inside the .bands file we get:"
    write(*,*) "Number of electrons (wrong!) = ", nelectrons(1:2)
    write(*,*) "Number of bands = ", nbands(1:2)
    write(*,*) "Fermi level = ", efermi(1:2)*hart2eV, " eV"
  end if

  read(200,*)
  read(200,'(3F12.6)') lvec(1:3,1)
  read(200,'(3F12.6)') lvec(1:3,2)
  read(200,'(3F12.6)') lvec(1:3,3)

  do nk=1,nkpts
    read(200,'(8x,I5,4F12.8)') tmpi, kpts(nk,1:3), wk(nk)
    do ns=1,nspins
      read(200,'(15x,I1)') tmpi
        do nb=1,nbands(ns)
          read(200, '(F14.8)') eigen(nb,nk,ns)
        end do
    end do
  end do

  close(200)
  
  ! Right, now we can generate the spectrum. We don't worry *at all* about occupancies
  ! or where the fermi level really is, etc - we just subtract the fermi level from
  ! all eigenvalues and the end user of the spectra can apply a Fermi function afterwards
  ! along with the proper smearing.
  
  efermi = efermi * hart2eV
  eigen = eigen * hart2eV
  
  ! If the user actually specified iminband, use it to find the minimum energy.
  if (iminband .gt. 0) then
    do nk=1,nkpts
      do ns=1,nspins
          do nb=iminband,nbands(ns)
            if (eigen(nb,nk,ns) .lt. e_min) then
              e_min = eigen(nb,nk,ns)
            end if
    !        if (eigen(nb,nk,ns) - efermi(ns) < w_start) then
    !          w_start = eigen(nb,nk,ns) - efermi(ns) - 5.0d0
    !        endif
    !        if (eigen(nb,nk,ns) - efermi(ns) > w_end) then
    !          w_end = eigen(nb,nk,ns) - efermi(ns) + 5.0d0
    !        endif 
          end do
      end do
    end do
  end if 
  
  print *, "Generating a spectrum from ", w_start, "to", w_end, 'eV'
  
  allocate(spectrum(spectrum_points, nspins, 6)) ! Separate spins
  allocate(cspectrum(spectrum_points,6)) ! This is the spectrum with combined spins
  
  spectrum = 0.0d0
  cspectrum = 0.0d0
  
  w_step = (w_end - w_start) / spectrum_points
  do iw=1,spectrum_points
    w(iw) = w_start + (iw-1)*w_step
  end do
  smear_width = w_step * 0.5
  !smear_width = 0.5
  orb = cproj
  
  ! Loop over bands, kpoints and spins - fix on a particular cproj. Note we can
  ! use OpenMP here if we want to.
!$OMP PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(ns, nk, nb, e_nks, f_nks, iw, w, smear_factor, icmpt, tmpstr) &
!$OMP& PRIVATE(matrix_cmpt)
  do ns=1,nspins
    do nk=1,nkpts
      do nb=iminband,mbands
        
        matrix_cmpt(1) = realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,1,nk,ns)))
        matrix_cmpt(2) = realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
        matrix_cmpt(3) = realpart(optmat(orb,nb,3,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
        matrix_cmpt(4) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
        matrix_cmpt(5) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
        matrix_cmpt(6) = 2.0d0 * realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
        
        ! If the user actually set a min band, use it as the zero of energy, otherwise
        ! zero the fermi level.
        if (iminband .gt. 0) then
          e_nks = eigen(nb,nk,ns) - e_min
        else
          e_nks = eigen(nb,nk,ns) - efermi(ns)
        end if
        
        print *, matrix_cmpt(1), optmat(orb,nb,1,nk,ns)
        ! Occupancy weighting: we don't deal with this at all anymore - f_nks is always 1.
        f_nks = 1.0d0
        
        do iw=1,spectrum_points
          ! Tiny lorentzian smear to project onto the nearest spectrum point.
          smear_factor = invpi * smear_width / ((e_nks - w(iw))**2 + (smear_width)**2)
          do icmpt=1,6
            spectrum(iw,ns,icmpt) = spectrum(iw,ns,icmpt) + wk(nk) * f_nks * matrix_cmpt(icmpt) * smear_factor
            cspectrum(iw,icmpt) = cspectrum(iw,icmpt) + wk(nk) * f_nks * matrix_cmpt(icmpt) * smear_factor
          end do
        end do
      end do
    end do
  end do
!$OMP END PARALLEL DO
    
  ! Zero tiny values to make sure the output can be read.
  do iw=1,spectrum_points
    do ns=1,nspins
      do icmpt=1,6
        if (dabs(spectrum(iw,ns,icmpt)) .lt. tol) then
          spectrum(iw,ns,icmpt) = 0.0d0
        end if
      end do
    end do
  end do

! Write spectra in the following files:
!   SEED_spinX.raw.nexafs for X=1->ns,
!   SEED.raw.nexafs - both spin channels combined.
!   SEED_spinX.smeared.nexafs for X=1->ns,
!   SEED.smeared.nexafs - smeared with both channels combined.

  do ns=1,nspins
    write(*,*) 'Writing raw spectrum for orbital ', orb, 'spin ', ns
    write(atstr(1), '(I4)') ns
    !write(atstr(2), '(I4)') core_ion(orb)
    !write(atstr(3), '(I4)') core_n(orb)
    !write(atstr(4), '(I4)') core_lm(orb)
    !tmpstr = trim(adjustl(atstr(1)))//'_'//trim(adjustl(atstr(2)))//'_'//trim(adjustl(atstr(3)))//'_'&
  !&             //trim(adjustl(atstr(4)))
    tmpstr = trim(adjustl(atstr(1)))
    open(300,file=trim(seed)//'_spin'//trim(adjustl(tmpstr))//'.raw.nexafs',form='formatted')
  
    write(300,*) '# NEXAFS core-level spectrum calculated by nexspec with CASTEP inputs.'
    write(300,*) '# Spin channel: ', ns
    write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'
  
    do iw=1,spectrum_points
      write(300, '(7g16.8)') w(iw), spectrum(iw,ns,1:6)
    end do
  
    close(300)
  end do

  write(*,*) 'Writing raw spectrum for orbital ', orb, 'with combined spins'
  open(300,file=trim(seed)//'.raw.nexafs',form='formatted')

  write(300,*) '# NEXAFS core-level spectrum calculated by nexspec with CASTEP inputs.'
  write(300,*) '#'
  write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'

  do iw=1,spectrum_points
    write(300, '(7g16.8)') w(iw), cspectrum(iw,1:6)
  end do

  close(300)
  ! Write out a nominally-smeared (0.2 Lorentzian, 0.4 eV Gaussian, no linear) spectrum
  
  do ns=1,nspins
    do icmpt=1,6
      call lorentzian_convolute(w,spectrum(:,ns,icmpt),spectrum_points, lwidth)
      call gaussian_convolute(w,spectrum(:,ns,icmpt),spectrum_points, gwidth)
    end do
  end do

  do ns=1,nspins  
    write(*,*) 'Writing smeared spectrum for orbital ', orb, 'spin ', ns
    write(atstr(1), '(I4)') ns
    !write(atstr(2), '(I4)') core_ion(orb)
    !write(atstr(3), '(I4)') core_n(orb)
    !write(atstr(4), '(I4)') core_lm(orb)
    !tmpstr = trim(adjustl(atstr(1)))//'_'//trim(adjustl(atstr(2)))//'_'//trim(adjustl(atstr(3)))//'_'&
  !&             //trim(adjustl(atstr(4)))
    tmpstr = trim(adjustl(atstr(1)))
    open(300,file=trim(seed)//'_spin'//trim(adjustl(tmpstr))//'.smeared.nexafs',form='formatted')
  
    write(300,*) '# NEXAFS core-level spectrum calculated by nexspec with CASTEP inputs.'
    write(300,*) '# Smeared with lorentzian and gaussian broadening,', lwidth, gwidth, 'eV respectively'
    write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'
  
    do iw=1,spectrum_points
      write(300, '(7g16.8)') w(iw), spectrum(iw,ns,1:6)
    end do
  
    close(300)
  end do

  do icmpt=1,6
    call lorentzian_convolute(w,cspectrum(:,icmpt),spectrum_points, lwidth)
    call gaussian_convolute(w,cspectrum(:,icmpt),spectrum_points, gwidth)
  end do 
  
  write(*,*) 'Writing raw spectrum for orbital ', orb, 'with combined spins'
  open(300,file=trim(seed)//'.smeared.nexafs',form='formatted')

  write(300,*) '# NEXAFS core-level spectrum calculated by nexspec with CASTEP inputs.'
  write(300,*) '# Smeared with lorentzian and gaussian broadening,', lwidth, gwidth, 'eV respectively'
  write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'

  do iw=1,spectrum_points
    write(300, '(7g16.8)') w(iw), cspectrum(iw,1:6)
  end do

  close(300)
  
  ! Deallocate everything
  deallocate(spectrum)
  deallocate(cspectrum)
  deallocate(wk)
  deallocate(optmat)
  deallocate(eigen)
  deallocate(kpts)
  deallocate(core_species)
  deallocate(core_ion)
  deallocate(core_n)
  deallocate(core_lm)
  
  write(*,*) "Finished nexspec - Goodbye!"
  
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

