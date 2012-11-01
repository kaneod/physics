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

module utilities

  use constants
  
  implicit none
  
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

end module     
 
module core_level_spectra

  use constants
  use utilities
  
  implicit none
  
  logical :: file_exists = .false.
  character(len=80) :: seed
  
  integer :: ncproj, mbands, nkpts, nspins, tmpi, nbands(2)
  integer, dimension(:), allocatable :: core_species, core_ion, core_n, core_lm
  complex(kind=dp), dimension(:,:,:,:,:), allocatable :: optmat
  real(kind=dp), dimension(:), allocatable :: wk, w
  real(kind=dp), dimension(:,:), allocatable :: kpts
  real(kind=dp), dimension(:,:,:), allocatable :: eigen, spectrum
  real(kind=dp), dimension(:,:), allocatable :: lastspec
  real(kind=dp) :: efermi(2), nelectrons(2), lvec(3,3), w_step
  real(kind=dp) :: matrix_cmpt(6), e_nks, f_nks, smear_factor
  integer :: ltstate ! index of the lowest state we can transition to.
  logical :: ltsingle ! account for singly-occupied HOMO states.
  real(kind=dp) :: lteigen ! lowest allowed transition eigenvalue.
  
  ! Default settings
  real(kind=dp) :: lorentzian_width = 0.3d0 ! eV
  real(kind=dp) :: w_start = -5.0d0 ! eV
  real(kind=dp) :: w_end = 50.0d0 ! eV
  integer :: spectrum_points = 2000
  real(kind=dp) :: linear_broadening = 0.0d0
  real(kind=dp) :: gaussian_broadening = 0.3d0 ! eV
  real(kind=dp) :: smear_value = 0.0d0
  logical :: is_core_hole = .true.
  real(kind=dp) :: hole_charge = 1.0d0 ! denotes a full core hole.
    
  contains
  
  subroutine deallocate_all
    !---------------------------------------------------------------------------
    !
    ! deallocate_all - deallocates all the allocatable arrays in the eels_mat
    ! module. Also resets default spectrum values.
    !
    !---------------------------------------------------------------------------
    
    if (allocated(optmat)) deallocate(optmat)
    if (allocated(core_species)) deallocate(core_species)
    if (allocated(core_ion)) deallocate(core_ion)
    if (allocated(core_n)) deallocate(core_n)
    if (allocated(core_lm)) deallocate(core_lm)
    if (allocated(wk)) deallocate(wk)
    if (allocated(w)) deallocate(w)
    if (allocated(kpts)) deallocate(kpts)
    if (allocated(eigen)) deallocate(eigen)
    if (allocated(spectrum)) deallocate(spectrum)
    if (allocated(lastspec)) deallocate(lastspec)
    
    ! Reset defaults
    lorentzian_width = 0.3d0 ! eV
    w_start = -5.0d0 ! eV
    w_end = 50.0d0 ! eV
    spectrum_points = 2000
    linear_broadening = 0.0d0
    gaussian_broadening = 0.3d0 ! eV
    smear_value = 0.0d0
    is_core_hole = .true.
    hole_charge = 1.0d0
    
  end subroutine
  
  subroutine init(newseed)
  
    !---------------------------------------------------------------------------
    !
    ! init - reads in the eels_mat and bands files to populate the optical
    ! matrix and eigenvalues as well as things like the number of atoms etc.
    ! Call this every time you switch SEEDnames or if you change anything.
    !
    !---------------------------------------------------------------------------

    integer :: nk, orb, ns, nb, iw, icmpt, inumelec
    character(len=80) :: newseed
    real(kind=dp) :: numelec
    
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
    
    ! Generate an energy axis based on w_start, w_end and spectrum_points
    allocate(w(spectrum_points))
    
    w_step = (w_end - w_start) / spectrum_points
    
    do iw=1,spectrum_points
      w(iw) = w_start + (iw-1)*w_step
    end do
    
    ! We need to engage in some trickery here. The problem is that typically
    ! in a core hole calculation, we have an electron sitting in what used
    ! to be the LUMO (CBM for crystal) for the neutral pre-excited system. So,
    ! the fermi level calculated by CASTEP is actually above where it would be
    ! before the excitation, and if we aren't careful we exclude the transition
    ! to the LUMO in the calculation of the spectrum. So, we need to find the
    ! subset of bands for which we want to actually calculate the spectrum.
    ! We do this by counting the number of fully occupied states only, and we
    ! set the first single or zero occupancy state to be the zero for the energy
    ! scale. 
    !
    ! This is *not* the conventional way of doing this! But it is a correct
    ! way.
    ! Four possibilities: can be corehole or non-corehole calculation in
    ! conjunction with either an odd or even number of electrons. Here "odd"
    ! includes any fractional values - can use a partial core hole.
    
    if (is_core_hole .eqv. .false.) then
      hole_charge = 0.0d0
    end if
    
    if (nspins .eq. 1) then
      numelec = nelectrons(1)
    else
      numelec = nelectrons(1) + nelectrons(2)
    end if
    
    ! True number of valence electrons in neutral system - must be an integer.
    inumelec = int(numelec - hole_charge)
    ! The lowest transition state is then integer/2 + 1
    ltstate = inumelec / 2 + 1
    
    ! Most common case first: fully-occupied HOMO in the neutral system and
    ! a core hole calculation.
    if ((mod(inumelec, 2) .eq. 0) .and. (is_core_hole)) then
      ltsingle = .false.
    ! Next case: core hole calculation on a single-occupancy HOMO system.
    else if ((mod(inumelec, 2) .eq. 1) .and. (is_core_hole)) then
      ltsingle = .true.
    ! Next case: non-corehole, fully-occupied HOMO.
    else if ((mod(inumelec, 2) .eq. 0) .and. (is_core_hole .eqv. .false.)) then
      ltsingle = .false.
    ! Final case: non-corehole, single-occupancy HOMO.
    else
      ltsingle = .true.
    end if
    
    ! So now we need to iterate over all the kpts and spins to find the lowest
    ! eigenvalue with index ltstate to use as our zero point.
    lteigen = 1e12 ! just some huge number
    do nk=1,nkpts
      do ns=1,nspins
        if (eigen(ltstate, nk, ns) < lteigen) then
          lteigen = eigen(ltstate, nk, ns)
        end if
      end do
    end do
  
  end subroutine
  
  subroutine generate_single_spectrum(orb)
  
    integer :: nk, ns, nb, iw, icmpt
    integer, intent(in) :: orb
    
    ! We want to subtract not the fermi level but the lowest allowed unoccupied
    ! state (determined in init). By doing this, we can add the Pickard/Gao
    ! transition energy and get a spectrum at the 'proper' DFT edge.
    eigen = eigen - lteigen
    !if (nspins .eq. 1) then
    !  eigen = eigen - efermi(1)
    !  efermi(1) = 0.0d0
    !else
    !  if (efermi(1) .gt. efermi(2)) then
    !    eigen = eigen - efermi(1)
    !    efermi(2) = efermi(2) - efermi(1)
    !    efermi(1) = 0.0d0
    !  else
    !    eigen = eigen - efermi(2)
    !    efermi(1) = efermi(1) - efermi(2)
    !    efermi(2) = 0.0d0
    !  end if
    !end if
    
    ! By default the smearing value is just a constant width.
    smear_value = lorentzian_width
    
    if (allocated(lastspec)) deallocate(lastspec)
    
    allocate(lastspec(spectrum_points,6))
    
    lastspec = 0.0d0
    
    do ns=1,nspins
      do nk=1,nkpts
        ! Note: only need to go over bands from ltstate and up.
        do nb=ltstate,mbands
          
          matrix_cmpt(1) = realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,1,nk,ns)))
          matrix_cmpt(2) = realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
          matrix_cmpt(3) = realpart(optmat(orb,nb,3,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          matrix_cmpt(4) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,2,nk,ns)))
          matrix_cmpt(5) = 2.0d0 * realpart(optmat(orb,nb,1,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          matrix_cmpt(6) = 2.0d0 * realpart(optmat(orb,nb,2,nk,ns) * dconjg(optmat(orb,nb,3,nk,ns)))
          
          e_nks = eigen(nb,nk,ns) * hart2eV
          ! Because of the way we've set the eigenvalues, don't need to set
          ! an occupancy here (code commented out, left for comparison).
          !if (eigen(nb,nk,ns) .ge. efermi(ns)) then
          !  f_nks = 1.0d0         
          !else
          !  f_nks = 0.0d0
          !end if
          f_nks = 1.0d0
          if ((nb .eq. ltstate) .and. (ltsingle)) then
            f_nks = 0.5d0 ! This needs to be modified to take into account
                          ! the possibility of a fractional core hole!
          end if
      
          do iw=1,spectrum_points
            smear_value = lorentzian_width + linear_broadening * w(iw)
          
            ! Lifetime broadening.
            smear_factor = invpi * smear_value / ((e_nks - w(iw))**2 + &
            &              (smear_value)**2)

            do icmpt=1,6
              lastspec(iw,icmpt) = lastspec(iw,icmpt) + wk(nk) * f_nks * &
              &                    matrix_cmpt(icmpt) * smear_factor
            end do
          end do
        end do
      end do
    end do
    
    ! Gaussian broadening.
    do icmpt=1,6
      call gaussian_convolute(w, lastspec(:,icmpt), spectrum_points, gaussian_broadening)
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
