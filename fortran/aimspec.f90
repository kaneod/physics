!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! aimspec.f90
!
! Simple NEXAFS spectrum generation code, based on nexspec.
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
! 1. Usage is just aimspec DATFILE CORELEVEL MINBAND. CORELEVEL is the index of the
! initial state, MINBAND is the index of the first unoccupied state. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program aimspec

  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: DEBUG = 1
  real(kind=dp), parameter :: invpi = 0.3183098861837907d0
  real(kind=dp), parameter :: pi = 3.141592653589793d0
  real(kind=dp), parameter :: invsqrt2pi = 0.3989422804014327d0
  real(kind=dp), parameter :: invfinestruct = 137.035999074d0
  real(kind=dp), parameter :: hart2eV = 27.211396132d0
  real(kind=dp), parameter :: tol = 1.0d-99
  
  real(kind=dp) :: smear_width, smear_factor ! eV
  integer, parameter :: spectrum_points = 2000
  real(kind=dp) :: w_start = -25.0 ! eV
  real(kind=dp) :: w_end = 65.0 ! eV
  real(kind=dp) :: w_step, w(spectrum_points)
  real(kind=dp) :: gwidth = 0.4 ! eV
  real(kind=dp) :: lwidth = 0.2 ! eV
  real(kind=dp), dimension(:,:,:,:), allocatable :: spectrum
  
  logical :: file_exists = .false.
  integer :: nargs, ioerr, linecount, tmpi, i
  real(kind=dp) :: tmpf
  character(len=80) :: datfile, tmpstr, atstr(2)
  
  integer :: cproj, iminband, matsize
  integer :: nb, si, sj, nspins, nbands, iw, icmpt
  real(kind=dp) :: linebits(13), e_min, matrix_cmpt(6), eig
  real(kind=dp), dimension(:,:), allocatable :: eigen
  real(kind=dp), dimension(:,:,:,:,:), allocatable :: optmat
  
  print *, "AIMSPEC version 1"
  print *, " "
  print *, "Written by Kane O'Donnell, August 2013"
  print *, " "
  
  ! Get arguments
  nargs = command_argument_count()
  if (nargs .eq. 3) then
    call get_command_argument(1,datfile)
    call get_command_argument(2,tmpstr)
    read(tmpstr, *) cproj
    call get_command_argument(3, tmpstr)
    read(tmpstr, *), iminband
  else  
    print *, " "
    print *, "USAGE: aimspec DATFILE CORELEVEL MINBAND"
    print *, " "
    print *, "Require DATFILE. CORELEVEL is the core initial state index."
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
  
  ! Read dat file twice. First time is to count the lines.
  
  ioerr = 0
  linecount = 0
  do
    read(100,*,IOSTAT=ioerr) tmpstr
    if (ioerr .eq. 0) then
      linecount = linecount + 1
    else
      exit
    end if
  end do
  
  rewind(100)
  
  print *, "Read ", linecount, " lines in the file."
  
  ! Use the linecount to dimension the eigenvalue and matrix element matrices.
  linecount = linecount - 2 ! first two lines are text
  tmpf = real(linecount)
  if (floor(tmpf / 3.0d0) .eq. nint(tmpf / 3.0d0)) then
    ! linecount is divisible by three
    nbands = nint(-0.5d0 + sqrt(1.0d0 + 8.0d0 * tmpf / 3.0d0) / 2.0d0)
    print *, "Size of matrix is, ", nbands
    nspins = 2
  else
    ! linecount is not divisible by three
    nbands = nint(-0.5d0 + sqrt(1.0d0 + 8.0d0 * tmpf) / 2.0d0)
    print *, "Size of matrix is, ", nbands
    nspins = 1
  end if
  
  allocate(eigen(nspins, nbands))
  allocate(optmat(nspins, nbands, nspins, nbands, 6))
  
  eigen = 0.0d0
  optmat = 0.0d0
  
  ! Time to read in the dat file. First two lines are text.
  read(100,*), tmpstr
  read(100,*), tmpstr
  
  ! Read into triangular matrix. We don't worry about the lower half of the triangle
  ! because we assume the final states are always bigger than the initial state.
  do i=1,linecount
    read(100,*), linebits(:)
    eigen(nint(linebits(2)), nint(linebits(1))) = linebits(3)
    optmat(nint(linebits(2)), nint(linebits(1)), nint(linebits(5)), nint(linebits(4)), :) = linebits(7:12)
  end do
  
  ! Establish the spectrum: in order to make this compatible with the CASTEP output,
  ! we simply use -25 eV to 65 eV as with nexspec.
  allocate(spectrum(spectrum_points, nspins, nspins, 6)) ! Separate spins
  
  spectrum = 0.0d0

  w_step = (w_end - w_start) / spectrum_points
  do iw=1,spectrum_points
    w(iw) = w_start + (iw-1)*w_step
  end do
  smear_width = w_step * 0.5
  
  ! Find the minimum eigenvalue to use as a zero.
  e_min = 1d6 ! Something absurdly large
  do si=1,nspins
    if (eigen(si,iminband) .lt. e_min) then
      e_min = eigen(si,iminband)
    end if
  end do
  
  print *, "Prepared for spectrum generation."
  print *, "Spectrum step size:", w_step
  print *, "Base smearing:", smear_width
  print *, "Minimum energy: ", e_min
  print *, "Core energy: ", eigen(:,cproj), " gives transition threshold ", e_min - eigen(:,cproj)
          
  ! Generate the three spectra: spin 1->1, 1->2, 2->2
  if (DEBUG .eq. 1) then
    print *, "Spin 1       Spin 2       Final State Index        Eigenvalue      |<f|dx|i>|^2"
  end if    
  do si=1,nspins
    do sj=si,nspins
      do nb=iminband,nbands
      
        matrix_cmpt(1) = optmat(si, cproj, sj, nb, 1) ** 2 + optmat(si, cproj, sj, nb, 2) ** 2
        matrix_cmpt(2) = optmat(si, cproj, sj, nb, 3) ** 2 + optmat(si, cproj, sj, nb, 4) ** 2
        matrix_cmpt(3) = optmat(si, cproj, sj, nb, 5) ** 2 + optmat(si, cproj, sj, nb, 6) ** 2
        matrix_cmpt(4) = 2.0d0 * (optmat(si, cproj, sj, nb, 1) * optmat(si, cproj, sj, nb, 3) + &
        &                         optmat(si, cproj, sj, nb, 2) * optmat(si, cproj, sj, nb, 4))
        matrix_cmpt(5) = 2.0d0 * (optmat(si, cproj, sj, nb, 1) * optmat(si, cproj, sj, nb, 5) + &
        &                         optmat(si, cproj, sj, nb, 2) * optmat(si, cproj, sj, nb, 6))
        matrix_cmpt(6) = 2.0d0 * (optmat(si, cproj, sj, nb, 3) * optmat(si, cproj, sj, nb, 5) + &
        &                         optmat(si, cproj, sj, nb, 4) * optmat(si, cproj, sj, nb, 6))
        
        eig = eigen(sj, nb) - e_min
        
        if (DEBUG .eq. 1) then
          print *, si, sj, nb, eig, matrix_cmpt(1)
        end if
        
        do iw=1,spectrum_points
          ! Tiny lorentzian smear to project onto the nearest spectrum point.
          smear_factor = invpi * smear_width / ((eig - w(iw))**2 + (smear_width)**2)
          do icmpt=1,6
            spectrum(iw,si,sj,icmpt) = spectrum(iw,si,sj,icmpt) +  matrix_cmpt(icmpt) * smear_factor
          end do
        end do 
      end do
    end do
  end do
  
  do iw=1,spectrum_points
    do si=1,nspins
      do sj=si,nspins
        do icmpt=1,6
          if (dabs(spectrum(iw,si,sj,icmpt)) .lt. tol) then
            spectrum(iw,si,sj,icmpt) = 0.0d0
          end if
        end do
      end do
    end do    
  end do
  
  ! Write out spectra.
  
  do si=1,nspins
    do sj=si,nspins
      write(*,*) 'Writing raw spectrum for orbital ', cproj, 'spin transition: ', si, '->', sj
      write(atstr(1), '(I1)') si
      write(atstr(2), '(I1)') sj
      tmpstr = trim(adjustl(atstr(1)))//trim(adjustl(atstr(2)))
      open(300,file=trim(datfile)//'_spin'//trim(adjustl(tmpstr))//'.raw.nexafs',form='formatted')
  
      write(300,*) '# NEXAFS core-level spectrum calculated by aimspec with FHI-aims inputs.'
      write(300,*) '# Spin transition: ', si, '->', sj
      write(300,*) '# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz'
  
      do iw=1,spectrum_points
        write(300, '(7g16.8)') w(iw), spectrum(iw,si,sj,1:6)
      end do
  
      close(300)
    end do
  end do
        
  close(100)
  
end program