module esclibf90

    implicit none

    subroutine read_abinit_density(filename)
        
        ! Header variables
        character*6         :: codvsn
        integer             :: headform, fform
        integer             :: bantot, date, intxc, ixc, natom, ngfft(3), nkpt, npsp, &
                        &   nspden, nspinor, nsppol, nsym, ntypat, occopt, pertcase, &
                        &   usepaw, usewvl, cplex, nspden
        double precision    :: acell(3), ecut, ecutdg, ecutsm, ecut_eff, qptn(3), &
                        &   rprimd(3,3), stmbias, tphysel, tsmear
        integer             :: istwfk(nkpt), nband(nkpt * nsppol), npwarr(nkpt), &
                        &   so_psp(npsp), symafm(nsym), symrel(3,3,nsym), typat(natom) &
                        &   nrhoijsel(nspden), rhoijselect(*,nspden)
        double precision    :: kpt(3,nkpt), occ(bantot), tnons(3,nsym), &
                        &   znucltypat(ntypat), wtk(nkpt)
        character*132       :: title
        double precision    :: znuclpsp, zionpsp
        integer             :: pspso, pspdat, pspcod, pspxc, lmax, lloc, mmax=integers
        double precision    :: residm, xred(3, natom), etotal, fermie, rhoij(*, nspden)
        
        ! Subroutine variables
        integer             :: status
        
        open(unit=100, file=filename, status='old', form='unformatted')
        
        read(100) codvsn, headform, fform
        
        close(100)
        
    end subroutine
    
end module