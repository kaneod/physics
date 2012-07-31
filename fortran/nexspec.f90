!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nexspec.f90
!
! Code for generating a NEXAFS spectrum from CASTEP ELNES task output.
!
! Heavily based on Shang-Peng Gao's original EELSplot code. Many thanks to Shang-Peng
! Gao and Chris Pickard for allowing this code to be GPL'd. 
!
!
! Written by Kane O'Donnell (Australian Synchrotron), June 2012.
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
! 1. nexspec doesn't read a parameters file yet - you have to change the options
!     below and recompile.
!
! 2. Compiling with gfortran -o nexspec nexspec.f90 should work - no fancy tricks
!     required in conventional cases. Then call nexspec with nexspec SEED just
!     like CASTEP.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program nexspec

  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: DEBUG = 1
  character(len=80),parameter :: ENDIAN = 'big_endian' ! CASTEP is conventionally compiled
  ! as a big-endian code: if you get ridiculous results from nexspec check if changing
  ! this to 'native' fixes it.
  character(len=80), parameter :: param_file = 'nexspec.param'
  real(kind=dp), parameter :: invpi = 0.3183098861837907d0
  real(kind=dp), parameter :: pi = 3.141592653589793d0
  real(kind=dp), parameter :: invsqrt2pi = 0.3989422804014327d0
  real(kind=dp), parameter :: invfinestruct = 137.035999074d0
  real(kind=dp), parameter :: hart2eV = 27.211396132d0
  real(kind=dp), parameter :: tol = 1.0d-16
  
  ! Spectrum parameters - will take from nexspec.param if it exists.
  integer :: smear_method = 1 ! 0 for Lorentzian, 1 for Gaussian.
  real(kind=dp) :: smear_width = 0.3d0 ! eV
  real(kind=dp) :: w_start = 275.0d0 ! eV
  real(kind=dp) :: w_end = 320.0d0 
  integer :: spectrum_points = 2000
  
  ! Options for more realistic spectra - EDIT THESE FOR YOUR SPECTRUM
  ! What is going on here:
  !
  ! The idea is that we want to be able to apply the proper physical prefactor for 
  ! optical adsorption. However, we can't do this if we're on the wrong energy
  ! scale, as the scaling is energy-dependent and non-linear. So, if scaling_prefactor
  ! is true, we need a non-zero (physical) core level. Can only handle one of these
  ! at the moment - eventually will have an array with e_core for each ion.
  logical :: scaling_prefactor = .true.
  real(kind=dp) :: e_core = -294.39d0 ! eV
  
  ! Parameter file
  logical :: file_exists = .false.
  character(len=80) :: seed, tmpstr, atstr(4)
  integer :: nargs, nk, orb, ns, nb, iw, icmpt
  
  integer :: ncproj, mbands, nkpts, nspins, tmpi, nbands(2)
  integer, dimension(:), allocatable :: core_species, core_ion, core_n, core_lm
  complex(kind=dp), dimension(:,:,:,:,:), allocatable :: optmat
  real(kind=dp), dimension(:), allocatable :: wk
  real(kind=dp), dimension(:,:), allocatable :: kpts
  real(kind=dp), dimension(:,:,:), allocatable :: eigen, spectrum
  real(kind=dp) :: efermi(2), nelectrons(2), lvec(3,3), w_step, w
  real(kind=dp) :: matrix_cmpt(6), e_nks, f_nks, smear_factor
  
  ! Check if nexspec.param is there: if so, set parameters.
  inquire(file=param_file, exist=file_exists)
  
  if (file_exists) then
    open(unit=300, file=param_file, form='formatted', status='old')
    read(300,*) smear_method
    read(300,*) smear_width
    read(300,*) w_start
    read(300,*) w_end
    read(300,*) spectrum_points
    read(300,*) scaling_prefactor
    read(300,*) e_core
    close(300)
  end if
  
  ! Get the seedname from the command line, fail with usage if nargs != 1.
  nargs = command_argument_count()
  if ((nargs .ne. 1) .and. (nargs .ne. 2)) then
    write(*,*) "Usage: nexspec SEED [CORELEVEL_ENERGY]"
    call exit(0)
  end if
  
  call get_command_argument(1,seed)
  if (nargs .eq. 2) then
    call get_command_argument(2,tmpstr)
    read(tmpstr, *) e_core
  end if
  
  ! Open the .eels_mat file and read the header.
  open(unit=100, file=adjustl(trim(seed))//'.eels_mat', form='unformatted',&
&     status='old', convert=ENDIAN)

  read(100) ncproj
  read(100) mbands
  read(100) nkpts
  read(100) nspins
  
  ! Echo some output.
  
  write(*,*) "Opened file ", adjustl(trim(seed))//'.eels_mat', " for NEXAFS spectrum generation."
  write(*,*) "Found ", ncproj, " core projectors."
  write(*,*) "Maximum number of bands per kpt/spin = ", mbands
  write(*,*) "Number of K-points = ", nkpts
  write(*,*) "Number of spins = ", nspins
  
  ! Allocate the various arrays that store information about the core orbitals.
  allocate(core_species(ncproj), core_ion(ncproj), core_n(ncproj), &
&           core_lm(ncproj))

  ! optmat is the optical transitions matrix.
  allocate(optmat(ncproj,mbands,3,nkpts,nspins))
  
  ! Allocate wk (kpoints weighting), eigen (the eigenvalues) and kpts
  allocate(wk(nkpts), eigen(mbands,nkpts,nspins), kpts(nkpts,3))
  
  ! Read the information about the core orbitals.
  read(100) core_species(1:ncproj)
  read(100) core_ion(1:ncproj)
  read(100) core_n(1:ncproj)
  read(100) core_lm(1:ncproj)
  
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
  open(unit=200, file=adjustl(trim(seed))//'.bands', form='formatted', &
&       status='old')

  ! Most of the header is redundant in the bands file but we need to cross-
  ! check the number of kpoints, spins and eigenvalues. We also get the
  ! lattice vectors, number of electrons and fermi-level for free!
  read(200,'(19x,I5)') tmpi
  
  if (tmpi .ne. nkpts) then
    write(*,*) "ERROR: Number of k-points in .bands is not equal to the number in .eels_mat."
    call exit(1)
  end if
  
  read(200,'(26x,I1)') tmpi
  if (tmpi .ne. nspins) then
    write(*,*) "ERROR: Number of spins in .bands is not equal to the number in .eels_mat."
    call exit(1)
  end if  
  
  if (DEBUG .eq. 1) then
    print *, "Made it past integer reads"
  end if
  if (nspins .eq. 1) then
    read(200,'(20x,g10.4)') nelectrons(1)
    read(200,'(22x,I6)') nbands(1)
    read(200,'(31x,F12.6)') efermi(1)
    write(*,*) "Inside the .bands file we get:"
    write(*,*) "Number of electrons = ", nelectrons(1)
    write(*,*) "Number of bands = ", nbands(1)
    write(*,*) "Fermi level = ", efermi(1), " Ha"
  else
    read(200,'(20x,2g10.4)') nelectrons(1:2)
    read(200,'(22x,2I6)') nbands(1:2)
    read(200,'(33x,2F12.6)') efermi(1:2)
    write(*,*) "Inside the .bands file we get:"
    write(*,*) "Number of electrons = ", nelectrons(1:2)
    write(*,*) "Number of bands = ", nbands(1:2)
    write(*,*) "Fermi level = ", efermi(1:2), " Ha"
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
  ! So: we have everything, let's generate the spectrum!
  
  ! First we should really generate an occupancy function.
  ! But, let's just go T=0K and take f_nks = 0 if e_nks < efermi
  
  ! Zero the fermi level.
  if (nspins .eq. 1) then
    eigen = eigen - efermi(1)
    efermi(1) = 0.0d0
  else
    if (efermi(1) .gt. efermi(2)) then
      eigen = eigen - efermi(1)
      efermi(2) = efermi(2) - efermi(1)
      efermi(1) = 0.0d0
    else
      eigen = eigen - efermi(2)
      efermi(1) = efermi(1) - efermi(2)
      efermi(2) = 0.0d0
    end if
  end if
  
  ! 6 spectrum components for arbitrary angle reconstruction
  allocate(spectrum(ncproj, spectrum_points, 6))
  
  spectrum = 0.0d0
  w_step = (w_end - w_start) / spectrum_points
  
  ! For each orbital, sum matrix elements over k, s and w for each b.
  do orb=1,ncproj
    do ns=1,nspins
      do nk=1,nkpts
        do nb=1,mbands
          
          matrix_cmpt(1) = realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,1,nk,ns)))
          matrix_cmpt(2) = realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
          matrix_cmpt(3) = realpart(optmat(orb,nb,3,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          matrix_cmpt(4) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
          matrix_cmpt(5) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          matrix_cmpt(6) = 2.0d0 * realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          
          e_nks = eigen(nb,nk,ns) * hart2eV
          if (eigen(nb,nk,ns) .ge. efermi(ns)) then
            f_nks = 1.0d0
          else
            f_nks = 0.0d0
          end if
          
          do iw=1,spectrum_points
            w = w_start + (iw-1)*w_step
            if (smear_method .eq. 0) then
              ! Lorentzian
              smear_factor = invpi * (smear_width / ((e_nks - e_core - w)**2 + smear_width**2))
            else
              ! Gaussian
              smear_factor = invsqrt2pi * 1.0d0 / smear_width * dexp(-0.5d0 * ((e_nks - e_core - w)/smear_width)**2)
            end if
            do icmpt=1,6
              spectrum(orb,iw,icmpt) = spectrum(orb,iw,icmpt) + wk(nk) * f_nks * matrix_cmpt(icmpt) * smear_factor
            end do
          end do
        end do
      end do
    end do
    
    ! Zero tiny values to make sure the output can be read by Python.
    do iw=1,spectrum_points
      do icmpt=1,6
        if (dabs(spectrum(orb,iw,icmpt)) .lt. tol) then
          spectrum(orb,iw,icmpt) = 0.0d0
        end if
      end do
    end do
    
    ! Optional - apply an inverse-frequency scaling. This makes the spectrum actually
    ! scale like the adsorption spectrum should, but on the other hand, relies on setting 
    ! the core energy correctly, otherwise the scaling is far too strong near
    ! the edge.
    ! Note: A volume factor is missing here. Should use the lattice parameters to generate
    ! the volume and divide that out as well.
    
    if (scaling_prefactor) then
      do iw=1,spectrum_points
        w = w_start + (iw-1)*w_step
        spectrum(1:ncproj,iw,1:6) = (8.0d0 * pi ** 2 * invfinestruct / (3.0d0 * w)) * &
&                                     spectrum(1:ncproj,iw,1:6)
      end do
    end if
    
    ! Output each spectrum to a file.
    !tmpstr = ''
    write(*,*) 'Writing orbital ', orb
    write(atstr(1), '(I4)') core_species(orb)
    write(atstr(2), '(I4)') core_ion(orb)
    write(atstr(3), '(I4)') core_n(orb)
    write(atstr(4), '(I4)') core_lm(orb)
    tmpstr = trim(adjustl(atstr(1)))//'_'//trim(adjustl(atstr(2)))//'_'//trim(adjustl(atstr(3)))//'_'&
&             //trim(adjustl(atstr(4)))
    open(300,file=trim(seed)//'_'//trim(adjustl(tmpstr))//'.nexafs',form='formatted')
    
    write(300,*) '# NEXAFS core-level spectrum calculated by nexspec with CASTEP inputs.'
    write(300,*) '#'
    write(300,*) '# Omega (Ha) Mxx Myy Mzz Mxy Mxz Myz'
    
    do iw=1,spectrum_points
      write(300, '(7g16.8)') w_start + (iw-1)*w_step, spectrum(orb,iw,1:6)
    end do
    close(300)
  end do
  
  ! Deallocate everything
  deallocate(spectrum)
  deallocate(wk)
  deallocate(optmat)
  deallocate(eigen)
  deallocate(kpts)
  deallocate(core_species)
  deallocate(core_ion)
  deallocate(core_n)
  deallocate(core_lm)
  
  write(*,*) "Finished nexspec - Goodbye!"
end program
