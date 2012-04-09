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


module io

  implicit none
  
  !integer, parameter :: dp = selected_real_kind(15)
  integer, parameter :: dp = 8 ! see note 2
  integer, parameter :: DEBUG = 1 ! Set to 1 and recompile to echo outputs.
    
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
  
  integer, allocatable, dimension(:,:,:) :: symrel
  
  real(kind=dp) :: qptn(3)
  real(kind=dp) :: rprimd(3,3)
  
  real(kind=dp), allocatable, dimension(:) :: occ, wtk, zionpsp, znuclpsp, &
  &                                           znucltypat
  
  real(kind=dp), allocatable, dimension(:,:) :: kptns, tnons, xred
  
  character(len=132), allocatable, dimension(:) :: title
  
  ! Density variables
  real(kind=dp), allocatable, dimension(:) :: rhor
  
  ! Wavefunction variables
  real(kind=dp), allocatable, dimension(:) :: cg, eigen
  real(kind=dp), allocatable, dimension(:,:) :: kg
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
    if (allocated(cg)) deallocate(cg)
    
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
        
    integer :: ipsp, bsize, lnspden
    integer, allocatable, dimension(:) :: ibuffer, nsel
    real(kind=dp), allocatable, dimension(:) :: buffer
        
    ! First line of the file is the code version, header format and file format.
    
    read(funit) codvsn, headform, fform
    
    if (DEBUG .eq. 1) then
      print *, "codvsn, headform, fform"
      print *, codvsn, headform, fform
    end if
    
    ! Read basic variables
    
    read(funit) bantot, date, intxc, ixc, natom, ngfft(:), nkpt, nspden, &
    &           nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw, &
    &           ecut, ecutdg, ecutsm, ecut_eff, qptn(:), rprimd(:,:), stmbias, &
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
    
    read(funit) istwfk(:), nband(:), npwarr(:), so_psp(:), symafm(:), &
    &           symrel(:,:,:), typat(:), kptns(:,:), occ(:), tnons(:,:), &
    &           znucltypat(:), wtk(:)
    
    ! Read the psp information
    
    do ipsp=1,npsp
      read(funit) title(ipsp), znuclpsp(ipsp), zionpsp(ipsp), pspso(ipsp), &
      &           pspdat(ipsp), pspcod(ipsp), pspxc(ipsp), lmn_size(ipsp)
    end do
    
    ! Final non-PAW record
    
    read(funit) residm, xred(:,:), etot, fermie
    
    ! Set cplex in case we need it later and aren't using PAW.
    cplex = 1
    
    ! PAW parts of the header. We do not save this - just read to put the file
    ! pointer past the header.
    
    if (usepaw==1) then
      allocate(nsel(natom))
      read(funit) nsel(:), cplex, lnspden
      bsize = sum(nsel)
      allocate(ibuffer(bsize))
      allocate(buffer(bsize * nspden * cplex))
      read(funit) ibuffer(:), buffer(:)
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
      read(funit) rhor(:)
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
    integer :: isppol, ikpt, ibantot, iband, ii, i, nbandk
    
    
    ! Have to start the loop before allocating because npw is actually in
    ! the file.
    
    ibantot = 0
    i = 0
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
        if (.not. allocated(kg)) allocate(kg(3,npw))
        read(funit) kg(1:3,1:npw)
        if (DEBUG .eq. 1) then
          print *, "kg(1:3, 1:npw)"
          print *, kg(1:3, 1:npw)
        end if
        read(funit) eigen(1+ibantot:nbandk+ibantot), &
        & occ(1+ibantot:nbandk+ibantot)
        if (DEBUG .eq. 1) then
          print *, "eigen(1+ibantot:nbandk+ibantot), occ(1+ibandtot:etc)"
          print *, eigen(1+ibantot:nbandk+ibantot), & 
          & occ(1+ibantot:nbandk+ibantot)
        end if
        if (.not. allocated(cg)) allocate(cg(2 * npw * nspinor * nbandk))
        do iband=1,nbandk
          read(funit) (cg(ii+i), ii=1,2*npw*nspinor)
        end do
        ibantot = ibantot + nbandk
        i = i + 2 * npw * nspinor * nbandk
      end do
    end do
    
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
