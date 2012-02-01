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
  
  ! Density
  real(kind=dp), allocatable, dimension(:) :: rhor
    
    
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
    
    ! Density
    if (allocated(rhor)) deallocate(rhor)
    
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
    
    ! Read basic variables
    
    read(funit) bantot, date, intxc, ixc, natom, ngfft(:), nkpt, nspden, &
    &           nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw, &
    &           ecut, ecutdg, ecutsm, ecut_eff, qptn(:), rprimd(:,:), stmbias, &
    &           tphysel, tsmear, usewvl
    
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
  
end module
