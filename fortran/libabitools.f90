!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! libabitools.f90
!
! Library of routines for dealing with the various files the ABINIT code
! produces - _DEN, _WFK, _POT etc.
!
! Designed to be compilable with f2py (leading to a python-compatible library.)
!
! Written by Kane O'Donnell (Australian Synchrotron), Jan 2012.
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
! 1. We only account for abinit file formats from version 5.7 onwards. Earlier
! versions will almost certainly fail at some point, especially very early (2.2)
! files which will fail at the first line.
!
! 2. f2py currently fails to parse selected_real_kind, so for the moment we
! set dp = 8 manually. This is fine for x86_64 and so on, probably not on Cray.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants

  implicit none
  
  integer, parameter :: dp = 8
  integer, parameter :: DEBUG=1
  integer, parameter :: FFTW_FORWARD=-1
  integer, parameter :: FFTW_BACKWARD=1
  integer, parameter :: FFTW_ESTIMATE=64
  
  double precision, parameter :: pi = 4.0d0 * datan(1.0d0)
  
  
end module
    
module io

  use constants
  
  implicit none
    
  !private ! Would like private but f2py can't handle it yet.
  
  ! Header variables (see src/43_abitypes_defs/defs_abitypes.F90)  
  integer :: bantot, date, headform, intxc, ixc, natom, nkpt, npsp, nspden
  integer :: nspinor, nsppol, nsym, ntypat, occopt, pertcase, usepaw, usewvl
  integer :: fform, cplex
  
  character(len=6) :: codvsn
    
  real(kind=dp) :: ecut, ecutdg, ecutsm, ecut_eff, etot, fermie, residm
  real(kind=dp) :: stmbias, tphysel, tsmear  
    
  integer :: ngfft(3), nwvlarr(2)
  
  integer, allocatable, dimension(:) :: istwfk, lmn_size, nband, npwarr, &
  &                                     pspcod, pspdat, pspso, pspxc, so_psp, &
  &                                     symafm, typat
  
  integer, allocatable, dimension(:,:) :: rhoijselect
  
  integer, allocatable, dimension(:,:,:) :: symrel
  
  real(kind=dp) :: qptn(3)
  real(kind=dp) :: rprimd(3,3)
  
  real(kind=dp), allocatable, dimension(:) :: occ, wtk, zionpsp, znuclpsp, &
  &                                           znucltypat
  
  real(kind=dp), allocatable, dimension(:,:) :: kptns, tnons, xred
  real(kind=dp), allocatable, dimension(:,:,:) :: rhoijp
  
  character(len=132), allocatable, dimension(:) :: title
  
  ! Density variables
  real(kind=dp), allocatable, dimension(:) :: rhor
  
  ! Wavefunction variables
  real(kind=dp), allocatable, dimension(:,:) :: cg
  real(kind=dp), allocatable, dimension(:) :: eigen
  integer, allocatable, dimension(:,:,:) :: kg
  integer :: npw
  
    
    
  contains
  
  subroutine deallocate_all
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! deallocate_all
    !
    ! Internal: deallocates all the allocated arrays in the io module.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Header
    if (allocated(istwfk)) deallocate(istwfk)
    if (allocated(lmn_size)) deallocate(lmn_size)
    if (allocated(nband)) deallocate(nband)
    if (allocated(npwarr)) deallocate(npwarr)
    if (allocated(pspcod)) deallocate(pspcod)
    if (allocated(pspdat)) deallocate(pspdat)
    if (allocated(pspso)) deallocate(pspso)
    if (allocated(pspxc)) deallocate(pspxc)
    if (allocated(so_psp)) deallocate(so_psp)
    if (allocated(symafm)) deallocate(symafm)
    if (allocated(typat)) deallocate(typat)
    if (allocated(symrel)) deallocate(symrel)
    if (allocated(occ)) deallocate(occ)
    if (allocated(wtk)) deallocate(wtk)
    if (allocated(zionpsp)) deallocate(zionpsp)
    if (allocated(znuclpsp)) deallocate(znuclpsp)
    if (allocated(znucltypat)) deallocate(znucltypat)
    if (allocated(kptns)) deallocate(kptns)
    if (allocated(tnons)) deallocate(tnons)
    if (allocated(xred)) deallocate(xred)
    if (allocated(title)) deallocate(title)
    if (allocated(eigen)) deallocate(eigen)
    if (allocated(rhoijselect)) deallocate(rhoijselect)
    if (allocated(rhoijp)) deallocate(rhoijp)
    
    ! Density
    if (allocated(rhor)) deallocate(rhor)
    
    ! Wavefunction
    if (allocated(eigen)) deallocate(eigen)
    if (allocated(cg)) deallocate(cg)
    if (allocated(kg)) deallocate(kg)
    
  end subroutine
    
  subroutine read_header(funit)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! read_header(funit)
    !
    ! Internal: reads header of abinit unformatted output and assigns to the
    ! module variables as appropriate.
    ! 
    ! Inputs: 
    !
    ! funit: file unit (file must already be opened).
    !
    ! Notes:
    !
    ! 1. read_header does not deallocate variables: make sure they are
    ! deallocated first!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input variables
        
    integer, intent(in) :: funit
        
    ! Local variables
        
    integer :: ipsp, bsize, lnspden, ii, jj, iatom, nselect, ispden
    integer, allocatable, dimension(:) :: ibuffer, nsel
    real(kind=dp), allocatable, dimension(:) :: buffer
        
    ! First line of the file is the code version, header format and file format.
    
    read(funit) codvsn, headform, fform
    
    if (DEBUG .eq. 1) then
      print *, "codvsn, headform, fform"
      print *, codvsn, headform, fform
    end if
    
    ! Read basic variables
    
    read(funit) bantot, date, intxc, ixc, natom, ngfft, nkpt, nspden, &
    &           nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw, &
    &           ecut, ecutdg, ecutsm, ecut_eff, qptn, rprimd, stmbias, &
    &           tphysel, tsmear, usewvl
    
    if (DEBUG .eq. 1) then
      print *, "bantot, date, intxc, ixc, natom, ngfft(1:3), nkpt, nspden"
      print *, bantot, date, intxc, ixc, natom, ngfft(1:3), nkpt, nspden
      print *, "nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw"
      print *, nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw
      print *, "ecut, ecutdg, ecutsm, ecut_eff, qptn(1:3), rprimd(1:3,1:3)"
      print *, ecut, ecutdg, ecutsm, ecut_eff, qptn(1:3), rprimd(1:3,1:3)
      print *, "stmbias, tphysel, tsmear, usewvl"
      print *, stmbias, tphysel, tsmear, usewvl
    end if
    
    ! We can now allocate our variables
    
    allocate(istwfk(nkpt))
    allocate(kptns(3, nkpt))
    allocate(lmn_size(npsp))
    allocate(nband(nkpt * nsppol))
    allocate(npwarr(nkpt))
    allocate(occ(bantot))
    allocate(pspcod(npsp))
    allocate(pspdat(npsp))
    allocate(pspso(npsp))
    allocate(pspxc(npsp))
    allocate(so_psp(npsp))
    allocate(symafm(nsym))
    allocate(symrel(3,3,nsym))
    allocate(title(npsp))
    allocate(tnons(3,nsym))
    allocate(typat(natom))
    allocate(wtk(nkpt))
    allocate(xred(3,natom))
    allocate(zionpsp(npsp))
    allocate(znuclpsp(npsp))
    allocate(znucltypat(ntypat))
    
    ! Read the non-psp variables
    
    read(funit) istwfk, nband, npwarr, so_psp, symafm, &
    &           symrel, typat, kptns, occ, tnons, &
    &           znucltypat, wtk
    
    ! Read the psp information
    
    do ipsp=1,npsp
      read(funit) title(ipsp), znuclpsp(ipsp), zionpsp(ipsp), pspso(ipsp), &
      &           pspdat(ipsp), pspcod(ipsp), pspxc(ipsp), lmn_size(ipsp)
    end do
    
    ! Final non-PAW record
    
    read(funit) residm, xred, etot, fermie
    
    ! Set cplex in case we need it later and aren't using PAW.
    cplex = 1
    
    ! PAW parts of the header. We don't actually init a pawrhoij data type
    ! for each atom here.
    
    if (usepaw==1) then
      allocate(nsel(natom))
      read(funit) nsel(:), cplex, lnspden
      bsize = sum(nsel)
      allocate(ibuffer(bsize))
      allocate(buffer(bsize * nspden * cplex))
      allocate(rhoijselect(natom, maxval(nsel)))
      allocate(rhoijp(natom, cplex * maxval(nsel), lnspden))
      if (DEBUG .eq. 1) then
        print *, "Allocated rhoijselect as size ", natom, maxval(nsel)
        print *, "Allocated rhoijp as size ", natom, cplex * maxval(nsel), lnspden
      end if
      read(funit) ibuffer, buffer
      ii=0
      jj=0
      do iatom=1,natom
        nselect = nsel(iatom)
        rhoijselect(iatom,1:nselect) = ibuffer(ii+1:ii+nselect)
        ii = ii + nselect
        do ispden=1,lnspden
          rhoijp(iatom, 1:cplex*nselect, ispden) = buffer(jj+1:jj+cplex*nselect)
          jj = jj + cplex*nselect
        end do
      end do
      
      deallocate(ibuffer, buffer, nsel)
    end if
    
  end subroutine
  
  subroutine read_density(funit)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! read_density(funit)
    !
    ! Internal: reads density rhor from an already-opened DEN file.
    !
    ! Input variables:
    !
    ! funit: file unit number for the DEN file.
    !
    ! Notes:
    !
    ! 1. rhor is a single-dimensional. To recast into a three-dimensional array,
    ! take the fast access to be the first primitive vector.
    !
    ! 2. We don't actually handle complex densities anywhere else at the moment
    ! even though we account for the possibilty of them here.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    integer, intent(in) :: funit
    integer :: ispden
    
    allocate(rhor(cplex * ngfft(1) * ngfft(2) * ngfft(3)))
    
    do ispden=1,nspden
      read(funit) rhor
    end do
  
  end subroutine
  
  subroutine read_wavefunction(funit)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! read_wavefunction(funit)
    !
    ! Internal: reads wavefunction coefficients cg and eigenvalues eigen from
    ! an already-opened WFK file.
    !
    ! Input variables:
    !
    ! fnit: file unit number for the WFK file.
    !
    ! Notes:
    !
    ! 1. cg is single dimensional and runs by band first, then kpt, then spin.
    ! 
    ! 2. eigen is also also 1D and runs by band, kpt, spin.
    !
    ! 3. Strictly speaking nband can vary as a function of kpt. However,
    ! in order to be able to allocate cg, we have to assume it is constant. We
    ! use the separate variable nbandk for this constant value.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    integer, intent(in) :: funit
    integer :: isppol, ikpt, ibantot, iband, ii, i, nbandk, ios, ipw, ncg
    
    ! Allocate our plane wave variables
    if (.not. allocated(kg)) allocate(kg(3, maxval(npwarr),nkpt))
    if (.not. allocated(eigen)) allocate(eigen(bantot))
    if (.not. allocated(occ)) allocate(occ(bantot))
    if (.not. allocated(cg)) allocate(cg(2, maxval(npwarr) * nsppol * maxval(nband) * nspinor * nkpt))
    
    ibantot = 0
    i = 0
    ipw = 1
    
    if (DEBUG .eq. 1) then
      print *, "mpw, mband = ", maxval(npwarr), maxval(nband)
    end if
    
    do isppol=1,nsppol
      if (DEBUG .eq. 1) then
        print *, "isppol = ", isppol
      end if
      do ikpt=1,nkpt
        if (DEBUG .eq. 1) then
          print *, "ikpt = ", ikpt
        end if
        read(funit) npw, nspinor, nbandk
        if (DEBUG .eq. 1) then
          print *, "npw, nspinor, nbandk"
          print *, npw, nspinor, nbandk
        end if
        read(funit) kg(1:3,1:npw,ikpt)
        read(funit) eigen(1+ibantot:nbandk+ibantot), &
        & occ(1+ibantot:nbandk+ibantot)
        if (DEBUG .eq. 1) then
          print *, "Got past kg, eigen and occ read."
        end if
        do iband=1,nbandk
          ipw = (iband-1) * npw * nspinor + i
          if (DEBUG .eq. 1) then
            print *, "Trying to read ", npw*nspinor, " values with ipw = ", ipw
          end if
          read(funit) cg(1:2,ipw+1:ipw+npw*nspinor)
        end do
        ibantot = ibantot + nbandk
        i = i + npw * nspinor * nbandk
      end do
    end do
    
    if (DEBUG .eq. 1) then
      print *, "Total size of cg is: ", maxval(npwarr) * nsppol * maxval(nband) * nspinor * nkpt
      print *, "Final i is: ", i
      print *, "Final ipw is: ", ipw
      print *, "Final ibantot is: ", ibantot
      print *, "cg shape is: ", shape(cg)
    end if
    
  end subroutine
  
  subroutine header(filename)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! header(filename)
    !
    ! For external access: opens, reads, and closes, an abinit file, only
    ! paying attention to the header content. Useful if you just want to get
    ! the atomic positions and so on.
    !
    ! Input variables:
    !
    ! filename: Path to the file you want to open.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    character(len=132), intent(in) :: filename
    
    ! Need to deallocate variables before calling read_header
    call deallocate_all    
    open(unit=100, file=filename, status='old', form='unformatted')       
    call read_header(100)        
    close(100)
    
  end subroutine
  
  subroutine density(filename)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! density(filename)
    !
    ! For external access: opens, reads, and closes, an abinit _DEN file.
    !
    ! Input variables:
    !
    ! filename: Path to the file you want to open.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    character(len=132), intent(in) :: filename
    
    call deallocate_all
    open(unit=100, file=filename, status='old', form='unformatted')
    call read_header(100)
    call read_density(100)
    close(100)
    
  end subroutine
  
  subroutine wavefunction(filename)
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! wavefunction(filename)
    ! 
    ! For external access: opens, reads and closes an abinit _WFK file.
    !
    ! Input variables:
    !
    ! filename: Path to the file you want to open.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    character(len=132), intent(in) :: filename
    
    call deallocate_all
    open(unit=100, file=filename, status='old', form='unformatted')
    call read_header(100)
    call read_wavefunction(100)
    close(100)
    
  end subroutine
  
end module

module wave
  
  ! This module is for things like converting a recip-space WF to real space, computing
  ! PW expansions, etc.
  
  use constants
  use io
  
  implicit none
  
  integer :: cband, ckpt
  complex(kind=dp), allocatable, dimension(:,:,:) :: cwf
  real(kind=dp) :: ceig, cocc
  logical :: normalized = .false.
  integer*8 :: plan
  
  contains
  
  subroutine recip_to_real(cn, ck, cs, option)
  
    integer, intent(in) :: cn, ck, cs, option ! option=1: normalize, 0: don't normalize
    
    integer :: isppol, ikpt, npw,ipw, i, j, k, nbandk, iband
    real(kind=dp), allocatable, dimension(:,:) :: cgnk
    real(kind=dp) :: gsum
    
    ! Find the right range of cg to read given cn, ck and cs.
    i = 0
    ipw = 1
    do isppol=1,nsppol
      do ikpt=1,nkpt
        npw = npwarr(ikpt)
        nbandk = nband(ikpt)
        do iband=1,nbandk
          ipw = (iband-1) * npw * nspinor + i
          if ((ikpt .eq. ck) .and. (iband .eq. cn) .and. (isppol .eq. cs)) then
            allocate(cgnk(2,npw*nspinor))
            cgnk = cg(1:2,ipw+1:ipw+npw*nspinor)
            goto 666
          end if
        end do
        ! Need code here to save the eigenvalue and occupancy.
        i = i + 2 * npw * nspinor * nbandk
      end do
    end do

    ! Now, assign the cgnk to the owf grid in the appropriate place.
666 npw = npwarr(ck)
    if (allocated(cwf)) deallocate(cwf)
    allocate(cwf(ngfft(1), ngfft(2), ngfft(3)))
    cwf = complex(0.0d0, 0.0d0)
    do ipw=1,npw
      i = kg(1,ipw,ck) + int(real(ngfft(1) + 1)/2)
      j = kg(2,ipw,ck) + int(real(ngfft(2) + 1)/2)
      k = kg(3,ipw,ck) + int(real(ngfft(3) + 1)/2)
      cwf(i,j,k) = complex(cgnk(1,ipw+(cs-1)*npw), cgnk(2,ipw+(cs-1)*npw))
    end do
    
    ! FFT cwf to real space.
    call dfftw_plan_dft_3d(plan,ngfft(1), ngfft(2), ngfft(3), cwf, cwf, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, cwf, cwf)
    call dfftw_destroy_plan(plan)
    
    ! If option=1, normalize with respect to the grid sum.
    if (option .eq. 1) then
      gsum = 0.0d0
      do k=1,ngfft(3)
        do j=1,ngfft(2)
          do i=1,ngfft(1)
            gsum = gsum + real(conjg(cwf(i,j,k)) * cwf(i,j,k))
          end do
        end do
      end do
      cwf = cwf / sqrt(gsum)
      normalized = .true.
    else
      normalized = .false.
    end if
  
  end subroutine

end module

!module paw

  ! Stores and manipulates PAW projectors and so on.
  
!  use constants
!  use io
!  use wave
  
!  integer :: nmesh, lmax
!  real(kind=dp) :: rmax
!  integer, allocatable, dimension(:) :: mesh_size, mesh_type, orbitals
!  real(kind=dp), allocatable, dimension(:) :: rad_step, log_step
!  
!  implicit none
  
module spectra

  use constants
  
  implicit none
  
  contains
  
  subroutine spectrum_axyz(cmpts, sout, evec, n)
  
    ! Faster version of spectrumAXYZ from esc_lib.
    
    integer, intent(in) :: n
    double precision, intent(in) :: cmpts(n,7), evec(3)
    double precision, intent(out) :: sout(n,2)
    integer :: i
    double precision :: nvec(3)
    
    sout = 0.0d0
    nvec = evec / dsqrt(evec(1)**2 + evec(2)**2 + evec(3)**2)
       
    do i=1,n
      sout(i,1) = cmpts(i,1)
      sout(i,2) = nvec(1)**2 * cmpts(i,2) + nvec(2)**2 * cmpts(i,3) + &
      &           nvec(3)**2 * cmpts(i,4) + nvec(1)*nvec(2) * cmpts(i,5) + &
      &           nvec(1)*nvec(3) * cmpts(i,6) + nvec(2)*nvec(3) * cmpts(i,7)
    end do
    
  end subroutine
  
  subroutine spectrum_atp(cmpts, sout, theta, phi, n)
  
    integer, intent(in) :: n
    double precision, intent(in) :: cmpts(n,7), theta, phi
    double precision, intent(out) :: sout(n,2)
    integer :: i
    double precision :: evec(3), tr, pr
    
    tr = pi * theta / 180.0d0
    pr = pi * phi / 180.0d0
    
    evec(1) = dcos(pr) * dsin(tr)
    evec(2) = dsin(pr) * dsin(tr)
    evec(3) = dcos(tr)
    
    call spectrum_axyz(cmpts, sout, evec, n)
    
  end subroutine
  
  subroutine spectrum_xyz(cmpts, sout, evec, n, m)
  
    ! In this version, cmpts has 3 slots (the first is for the atom).
    
    integer, intent(in) :: n, m
    double precision, intent(in) :: cmpts(m,n,7), evec(3)
    double precision, intent(out) :: sout(n,2)
    integer :: i,j
    double precision :: nvec(3)
    
    sout = 0.0d0
    nvec = evec / dsqrt(evec(1)**2 + evec(2)**2 + evec(3)**2)
    
    do i=1,n
      sout(i,1) = cmpts(1,i,1)
    end do
       
    do i=1,n
      do j=1,m
      sout(i,2) = nvec(1)**2 * cmpts(j,i,2) + nvec(2)**2 * cmpts(j,i,3) + &
      &           nvec(3)**2 * cmpts(j,i,4) + nvec(1)*nvec(2) * cmpts(j,i,5) + &
      &           nvec(1)*nvec(3) * cmpts(j,i,6) + &
      &           nvec(2)*nvec(3) * cmpts(j,i,7)
      end do
    end do
    
    ! Divide by the number of atoms.
    sout(i,2) = sout(i,2) / m  
    
  end subroutine
  
  subroutine spectrum_tp(cmpts, sout, theta, phi, n, m)
  
    integer, intent(in) :: n, m
    double precision, intent(in) :: cmpts(m,n,7), theta, phi
    double precision, intent(out) :: sout(n,2)
    integer :: i
    double precision :: evec(3), tr, pr
    
    tr = pi * theta / 180.0d0
    pr = pi * phi / 180.0d0
    
    evec(1) = dcos(pr) * dsin(tr)
    evec(2) = dsin(pr) * dsin(tr)
    evec(3) = dcos(tr)
    
    call spectrum_xyz(cmpts, sout, evec, n,m)
    
  end subroutine
  
  integer function closest_index_to_energy(energy_array, energy, n)
  
    ! Returns the index i_closest such that energy_array(i_closest) is the 
    ! best match to energy.
    
    integer, intent(in) :: n
    double precision, intent(in) :: energy_array(n), energy
    integer :: i
    double precision :: emin
    
    closest_index_to_energy = 0
    emin = 1.0d10 ! something ridiculously huge to start
    
    do i=1,n
      if (dabs(energy_array(i) - energy).lt.emin) then
        closest_index_to_energy = i
        emin = dabs(energy_array(i) - energy)
      end if
    end do
    
    return
    
  end function
end module
