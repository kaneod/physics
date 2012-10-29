!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! libpytep.f90
!
! Module code providing a python interface to miscellaneous CASTEP-related 
! outputs.
!
! Written by Kane O'Donnell (Australian Synchrotron), October 2012.
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
! 1. The only interface for the moment is to the eels_mat output (based on 
!    nexspec).
!
! 2. Spectrum generation is NOT THE SAME as nexspec - nexspec has a bug at the
!    moment that means it can't be trusted at the moment and the spectrum
!    stuff here is greatly simplified to try and debug nexspec.
!
! 3. This file should be compiled with f2py to generate the python interface.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants

  implicit none
  
  ! Can't use selected_real_kind with f2py yet.
  integer, parameter :: dp = 8
  integer, parameter :: DEBUG = 1
  
  ! CASTEP is conventionally compiled as a big-endian code.
  character(len=80), parameter :: ENDIAN = 'big_endian' 
  
  real(kind=dp), parameter :: invpi = 0.3183098861837907d0
  real(kind=dp), parameter :: pi = 3.141592653589793d0
  real(kind=dp), parameter :: invsqrt2pi = 0.3989422804014327d0
  real(kind=dp), parameter :: invfinestruct = 137.035999074d0
  real(kind=dp), parameter :: hart2eV = 27.211396132d0
  real(kind=dp), parameter :: tol = 1.0d-40
  
end module

module core_level_spectra

  use constants
  
  implicit none
  
  logical :: file_exists = .false.
  character(len=80) :: seed, tmpstr, atstr(4)
  
  integer :: ncproj, mbands, nkpts, nspins, tmpi, nbands(2)
  integer, dimension(:), allocatable :: core_species, core_ion, core_n, core_lm
  complex(kind=dp), dimension(:,:,:,:,:), allocatable :: optmat
  real(kind=dp), dimension(:), allocatable :: wk
  real(kind=dp), dimension(:,:), allocatable :: kpts
  real(kind=dp), dimension(:,:,:), allocatable :: eigen, spectrum
  real(kind=dp), dimension(:,:), allocatable :: lastspec
  real(kind=dp) :: efermi(2), nelectrons(2), lvec(3,3), w_step, w
  real(kind=dp) :: matrix_cmpt(6), e_nks, f_nks, smear_factor
  
  ! Default settings
  integer :: smear_method = 1 ! 0 for Lorentzian, 1 for Gaussian.
  real(kind=dp) :: smear_width = 0.3d0 ! eV
  real(kind=dp) :: w_start = -5.0d0 ! eV
  real(kind=dp) :: w_end = 50.0d0 
  integer :: spectrum_points = 2000
  !logical :: scaling_prefactor = .false.
  !real(kind=dp) :: e_core = -290.0d0 ! eV
  !logical :: non_uniform_broadening = .false.
  !real(kind=dp) :: smear_offset = 145.0d0 ! eV
  real(kind=dp) :: smear_value = 0.0d0
  !real(kind=dp) :: smear_intercept = 0.0d0
  !real(kind=dp) :: smear_gradient = 0.0d0
  !logical :: zero_fermi = .true.
    
  contains
  
  subroutine deallocate_all
    !---------------------------------------------------------------------------
    !
    ! deallocate_all - deallocates all the allocatable arrays in the eels_mat
    ! module.
    !
    !---------------------------------------------------------------------------
    
    if (allocated(optmat)) deallocate(optmat)
    if (allocated(core_species)) deallocate(core_species)
    if (allocated(core_ion)) deallocate(core_ion)
    if (allocated(core_n)) deallocate(core_n)
    if (allocated(core_lm)) deallocate(core_lm)
    if (allocated(wk)) deallocate(wk)
    if (allocated(kpts)) deallocate(kpts)
    if (allocated(eigen)) deallocate(eigen)
    if (allocated(spectrum)) deallocate(spectrum)
    if (allocated(lastspec)) deallocate(lastspec)
    
  end subroutine
  
  subroutine init(newseed)

    integer :: nk, orb, ns, nb, iw, icmpt
    character(len=80) :: newseed
    
    inquire(file=adjustl(trim(newseed))//'.eels_mat', exist=file_exists)
    
    if (file_exists) then
      call deallocate_all
      open(unit=100, file=adjustl(trim(newseed))//'.eels_mat', &
      &    form='unformatted', status='old', convert=ENDIAN)
      !&    form='unformatted', status='old')
      seed = newseed
    else
      print *, "Error: The required eels_mat file does not exist. Exiting without altering state."
      return
    end if
    
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
    &    status='old')

    ! Most of the header is redundant in the bands file but we need to cross-
    ! check the number of kpoints, spins and eigenvalues. We also get the
    ! lattice vectors, number of electrons and fermi-level for free!
    read(200,'(19x,I5)') tmpi
  
    if (tmpi .ne. nkpts) then
      write(*,*) "ERROR: Number of k-points in .bands is not equal to the number in .eels_mat."
      write(*,*) "ERROR: The module is now in an uncertain state: re-initialize before use."
      return
    end if
  
    read(200,'(26x,I1)') tmpi
    if (tmpi .ne. nspins) then
      write(*,*) "ERROR: Number of spins in .bands is not equal to the number in .eels_mat."
      write(*,*) "ERROR: The module is now in an uncertain state: re-initialize before use."
      return
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
  
  end subroutine
  
  subroutine generate_single_spectrum(orb)
  
    integer :: nk, ns, nb, iw, icmpt
    integer, intent(in) :: orb
    
    ! Subtract the fermi level from all eigenvalues - we calculate the spectra
    ! relative to zero eV - false edge.
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
    
    ! Smearing. For now, just set to a constant value.
    smear_value = smear_width
    
    w_step = (w_end - w_start) / spectrum_points
    
    if (allocated(lastspec)) deallocate(lastspec)
    
    allocate(lastspec(spectrum_points,6))
    
    do ns=1,nspins
      do nk=1,nkpts
        do nb=1,mbands
        
          if (DEBUG .eq. 1) then
            write(*,*) "Spin, Kpt, Band: ", ns, nk, nb
          end if
          
          matrix_cmpt(1) = realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,1,nk,ns)))
          matrix_cmpt(2) = realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
          matrix_cmpt(3) = realpart(optmat(orb,nb,3,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          matrix_cmpt(4) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
          matrix_cmpt(5) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          matrix_cmpt(6) = 2.0d0 * realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          
          if (DEBUG .eq. 1) then
            write(*,*) matrix_cmpt
          end if
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
              smear_factor = invpi * smear_value / ((e_nks - w)**2 + &
              &              (smear_value)**2)
            else
              ! Gaussian
              smear_factor = invsqrt2pi * 1.0d0 / smear_value * dexp(-0.5d0 * &
              &              ((e_nks - w)/smear_value)**2)
            end if 
            do icmpt=1,6
              lastspec(iw,icmpt) = lastspec(iw,icmpt) + wk(nk) * f_nks * &
              &                    matrix_cmpt(icmpt) * smear_factor
            end do
          end do
        end do
      end do
    end do
    
    ! Zero tiny values
    do iw=1,spectrum_points
      do icmpt=1,6
        if (dabs(lastspec(iw,icmpt)) .lt. tol) then
          lastspec(iw,icmpt) = 0.0d0
        end if
      end do
    end do
    
  end subroutine
    
end module    
