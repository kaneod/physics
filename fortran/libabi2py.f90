module utilities

    implicit none
    
    integer, parameter :: dp = selected_real_kind(15)

    ! Header variables (v5.7 and up)
    character(6)        :: codvsn
    integer             :: headform, fform
        
    integer             :: bantot, date, intxc, ixc, natom, ngfft(3), nkpt, npsp, &
                    &   nspden, nspinor, nsppol, nsym, ntypat, occopt, pertcase, &
                    &   usepaw, usewvl, cplex
    real(kind=dp)       :: acell(3), ecut, ecutdg, ecutsm, ecut_eff, qptn(3), &
                    &   rprimd(3, 3), stmbias, tphysel, tsmear        
    integer, allocatable, dimension(:) :: istwfk, nband, npwarr, so_psp, symafm, &
                    &   typat, nrhoijsel, lmn_size, pspso, pspdat, pspcod, pspxc
    integer, allocatable, dimension(:,:,:) :: symrel
    real(kind=dp), allocatable, dimension(:) :: occ, znucltypat, wtk, &
                    & znuclpsp, zionpsp
    character(len=132), allocatable, dimension(:) :: title
    real(kind=dp)       :: residm, etotal, fermie
    real(kind=dp), allocatable, dimension(:,:) :: kpt, tnons, xred
    
    ! Fields
    real(kind=dp), allocatable, dimension(:) :: rhor
    
    !private read_header, read_density, deallocate_all
    !public :: header
    
    contains
    
    subroutine deallocate_all
        ! Does exactly what you would think...
        
        if (allocated(istwfk)) deallocate(istwfk)
        if (allocated(nband)) deallocate(nband)
        if (allocated(npwarr)) deallocate(npwarr)
        if (allocated(so_psp)) deallocate(so_psp)
        if (allocated(symafm)) deallocate(symafm)
        if (allocated(typat)) deallocate(typat)
        if (allocated(nrhoijsel)) deallocate(nrhoijsel)
        if (allocated(lmn_size)) deallocate(lmn_size)
        if (allocated(pspso)) deallocate(pspso)
        if (allocated(pspdat)) deallocate(pspdat)
        if (allocated(pspcod)) deallocate(pspcod)
        if (allocated(pspxc)) deallocate(pspxc)
        if (allocated(symrel)) deallocate(symrel)
        if (allocated(occ)) deallocate(occ)
        if (allocated(znucltypat)) deallocate(znucltypat)
        if (allocated(wtk)) deallocate(wtk)
        if (allocated(znuclpsp)) deallocate(znuclpsp)
        if (allocated(zionpsp)) deallocate(zionpsp)
        if (allocated(title)) deallocate(title)
        if (allocated(kpt)) deallocate(kpt)
        if (allocated(tnons)) deallocate(tnons)
        if (allocated(xred)) deallocate(xred)
        if (allocated(rhor)) deallocate(rhor)
        
    end subroutine       

    subroutine read_header(funit)
        ! Reads the header of an ALREADY OPEN file. This is private to prevent
        ! screwups with the file unit - all public functions admit a filename.
        
        ! Input variables
        
        integer, intent(in) :: funit
        
        ! Subroutine variables
        
        integer                                     :: bsize, ipsp, lnspden
        integer, allocatable, dimension(:)          :: ibuffer
        real(kind=dp), allocatable, dimension(:)    :: buffer
       
        read(funit) codvsn, headform, fform
        read(funit) bantot, date, intxc, ixc, natom, ngfft(1:3), nkpt, nspden, &
                & nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw, &
                & ecut, ecutdg, ecutsm, ecut_eff, qptn(1:3), rprimd(1:3, 1:3), &
                & stmbias, tphysel, tsmear, usewvl
                
        allocate(istwfk(nkpt))
        allocate(nband(nkpt * nsppol))
        allocate(npwarr(nkpt))
        allocate(so_psp(npsp))
        allocate(symafm(nsym))
        allocate(symrel(3, 3, nsym))
        allocate(typat(natom))
        allocate(lmn_size(npsp))
        allocate(kpt(3, nkpt))
        allocate(occ(bantot))
        allocate(pspso(npsp))
        allocate(pspcod(npsp))
        allocate(pspdat(npsp))
        allocate(pspxc(npsp))
        allocate(tnons(3, nsym))
        allocate(znucltypat(ntypat))
        allocate(wtk(nkpt))
        allocate(xred(3, natom))
        allocate(title(npsp))
        allocate(znuclpsp(npsp))
        allocate(zionpsp(npsp))
        
        read(funit) istwfk(1:nkpt), nband(1:nkpt * nsppol), npwarr(1:nkpt), &
                & so_psp(1:npsp), symafm(1:nsym), symrel(1:3, 1:3, 1:nsym), &
                & typat(1:natom), kpt(1:3, 1:nkpt), occ(1:bantot), &
                & tnons(1:3, 1:nsym), znucltypat(1:ntypat), wtk(1:nkpt)
        
        do ipsp=1,npsp
            read(funit) title(ipsp), znuclpsp(ipsp), zionpsp(ipsp), pspso(ipsp), &
                    & pspdat(ipsp), pspcod(ipsp), pspxc(ipsp), lmn_size(ipsp)
        end do
        
        read(funit) residm, xred(:,:), etotal, fermie
        
        ! In case usepaw==0, set cplex.
        cplex = 1
        
        if (usepaw==1) then
            allocate(nrhoijsel(natom))
            read(funit) nrhoijsel(:), cplex, lnspden
            bsize = sum(nrhoijsel)
            allocate(ibuffer(bsize))
            allocate(buffer(bsize * nspden * cplex))
            read(funit) ibuffer(:), buffer(:)
            ! Currently we don't do anything with these buffers, we 
            ! just read them to get past the header.
        end if       
        
    end subroutine
    
    subroutine read_density(funit)
        
        ! Inputs
        integer, intent(in) :: funit
        
        ! Subroutine variables
        integer             :: ispden
        
        allocate(rhor(cplex * ngfft(1) * ngfft(2) * ngfft(3)))
        do ispden=1, nspden
            read(funit) rhor(:)
            print *, rhor(:)
        end do
    
    end subroutine
    
    subroutine header(filename)
        
        character(132), intent(in) :: filename
        
        call deallocate_all
        
        open(unit=100, file=filename, status='old', form='unformatted')
        
        call read_header(100)
        
        close(100)
    
    end subroutine
    
    subroutine density(filename)
    
        character(132), intent(in) :: filename
        
        call deallocate_all
        
        open(unit=100, file=filename, status='old', form='unformatted')
        
        call read_header(100)
        call read_density(100)
        
        close(100)
    
    end subroutine 
    
end module        
