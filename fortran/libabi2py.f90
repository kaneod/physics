module abilib

    implicit none
    
    integer, parameter :: dp = selected_real_kind(15)
    character(6)        :: codvsn
    integer             :: headform, fform
        
    integer             :: bantot, date, intxc, ixc, natom, ngfft(3), nkpt, npsp, &
                    &   nspden, nspinor, nsppol, nsym, ntypat, occopt, pertcase, &
                    &   usepaw, usewvl, cplex
    real(kind=dp)       :: acell(3), ecut, ecutdg, ecutsm, ecut_eff, qptn(3), &
                    &   rprimd(3,3), stmbias, tphysel, tsmear        
    integer, allocatable, dimension(:) :: istwfk, nband, npwarr, so_psp, symafm, &
                    &   typat, nrhoijsel
    integer, allocatable, dimension(:,:,:) :: symrel
    integer, allocatable, dimension(:,:) :: rhoijselect
    real(kind=dp), allocatable, dimension(:) :: occ, znucltypat, wtk
    character(132)      :: title
    real(kind=dp)       :: znuclpsp, zionpsp
    integer             :: pspso, pspdat, pspcod, pspxc, lmax, lloc, mmax
    real(kind=dp)       :: residm, etotal, fermie
    real(kind=dp), allocatable, dimension(:,:) :: rhoij, kpt, tnons, xred
    
    contains

    subroutine read_abinit_density(filename)
        
        ! Input variables
        
        character(132), intent(in) :: filename
        
        ! Subroutine variables
        ! (none)

! f2py specifications
!f2py   intent(in)  filename
       
        open(unit=100, file=filename, status='old', form='unformatted')
        
        read(100) codvsn, headform, fform
        
        close(100)
        
    end subroutine
    
end module
