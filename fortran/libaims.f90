!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! libaims.f90
!
! Module code providing a python interface to miscellaneous CASTEP-related 
! outputs.
!
! Written by Kane O'Donnell (Australian Synchrotron), July 2013.
!
! Contact the author on gmail, username kane dot odonnell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright 2013 Kane O'Donnell
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
! 1. This is for compilation with f2py and provides functionality to python
! where speed is necessary.
!
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
  
  subroutine gaussian_convolute(xdata, ydata, outdata, num_points, gwidth)
    !---------------------------------------------------------------------------
    !
    ! gaussian_convolute: a conventional convolution where the gaussian function
    ! of width gwidth is the second factor in the integral.
    !
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: num_points
    real(kind=dp), intent(in) :: xdata(num_points), ydata(num_points)
    real(kind=dp), intent(out) :: outdata(num_points)
    real(kind=dp) :: dx
    real(kind=dp), intent(in) :: gwidth
    integer :: i, j
    
    ! Assume the step is even (!!)
    dx = xdata(2) - xdata(1)
    
    outdata = 0.0d0
    
    do i=1,num_points
      do j=1,num_points
        outdata(i) = outdata(i) + dx * ydata(j) * invsqrt2pi * 1.0d0 / gwidth * dexp( &
        &        -0.5d0 * ((xdata(i) - xdata(j)) / gwidth)**2)
      end do
    end do
    
    return
  end subroutine
  
  subroutine lorentzian_convolute(xdata, ydata, outdata, num_points, lwidth)
    !---------------------------------------------------------------------------
    !
    ! lorentzian_convolute: a conventional convolution where the lorentzian function
    ! of width gwidth is the second factor in the integral.
    !
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: num_points
    real(kind=dp), intent(in) :: xdata(num_points), ydata(num_points)
    real(kind=dp), intent(out) :: outdata(num_points)
    real(kind=dp) :: dx
    real(kind=dp), intent(in) :: lwidth
    integer :: i, j
    
    ! Assume the step is even (!!)
    dx = xdata(2) - xdata(1)
    
    outdata = 0.0d0
    
    do i=1,num_points
      do j=1,num_points
        outdata(i) = outdata(i) + dx * ydata(j) * invpi * (lwidth / ((xdata(i) - xdata(j))**2 + &
        & lwidth**2))
      end do
    end do
    
    return
  end subroutine
  
  subroutine lorentzian_linear_convolute(xdata, ydata, outdata, num_points, lwidth, llin)
    !---------------------------------------------------------------------------
    !
    ! lorentzian_linear_convolute: a conventional convolution where the lorentzian function
    ! of width gwidth is the second factor in the integral. The linear broadening is 
    ! applied as lwidth + llin * | x(j) |, where x(j) is the x value of the spectrum
    ! being convoluted.  
    !
    ! We do a subtraction so that the spectrum is 0-based to make the smearing more
    ! convenient (in terms of the linear smearing width being similar to the lorentzian)
    !
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: num_points
    real(kind=dp), intent(in) :: xdata(num_points), ydata(num_points)
    real(kind=dp), intent(out) :: outdata(num_points)
    real(kind=dp) :: dx
    real(kind=dp), intent(in) :: lwidth, llin
    integer :: i, j
    
    ! Assume the step is even (!!)
    dx = xdata(2) - xdata(1)
    
    outdata = 0.0d0
    
    do i=1,num_points
      do j=1,num_points
        outdata(i) = outdata(i) + dx * ydata(j) * invpi * ((lwidth + llin * dabs( &
        & xdata(j)-xdata(1)))/ ((xdata(i) - xdata(j))**2 + (lwidth + llin * dabs( &
        & xdata(j) - xdata(1)))**2))
      end do
    end do
    
    return
  end subroutine

  subroutine integrate(xdata, ydata, integral, num_points)
    !---------------------------------------------------------------------------
    !
    ! integrate: Paired-point Trapezoidal integration of a function. The returned
    ! ydata array contains the cumulative integral so watch out!
    !
    !---------------------------------------------------------------------------
    integer, intent(in) :: num_points
    real(kind=dp), intent(out) :: integral
    real(kind=dp), intent(inout) :: xdata(num_points), ydata(num_points)
    real(kind=dp) :: a, b, fa, fb, tmpsum, tmp(num_points)
    integer :: i
    
    tmpsum = 0.0d0
    integral = 0.0d0
    tmp = 0.0d0
    a = xdata(1)
    fa = ydata(1)
    
    do i=2,num_points
      b = xdata(i)
      fb = ydata(i)
      tmpsum = tmpsum + (b - a) * (fa + fb) * 0.5d0
      tmp(i) = tmpsum
      a = b
      fa = fb
    end do
    
    ydata = tmp
    integral = tmpsum
    
    return
  end subroutine
    
end module 

module spectroscopy

  use constants
  use utilities
  
  implicit none
  
  contains
  
  subroutine generate_spectrum(optmat, opttrans, ecore, w, spectrum, efermi, i_core, &
  &                             nkpts, nstates, ntrans, npts)
    !-------------------------------------------------------------------------------------
    !
    ! generate_spectrum: Does the nested-loop heavy lifting when generating a spectrum
    ! via Fermi's golden rule. Uses a tiny smearing value (half the step width) so the 
    ! spectrum can be broadened for lifetime and instrument effects later on.
    !
    ! This function computes sum_k sum_j M_ijku^2 delta(Ejk - Eik - w) for each of the 6 
    ! spectral components M_ijk1...M_ijk6 for a fixed starting energy level Eik indexed
    ! by i_core. The matrix elements are contained in the input optmat, the transitions
    ! Ejk - Eik are contained in opttrans, ecore is the core level eigenvalue, w is the
    ! array of spectrum points, spectrum is the intended output array, efermi is
    ! the system fermi level used to determine if a transition is included or
    ! not and the other inputs are dimensions.
    !
    !-------------------------------------------------------------------------------------                          

    integer, intent(in) :: i_core, nkpts, nstates, ntrans, npts
    real(kind=dp), intent(in) :: optmat(nkpts, nstates, ntrans, 6)
    real(kind=dp), intent(in) :: opttrans(nkpts, nstates, ntrans)
    real(kind=dp), intent(in) :: ecore
    real(kind=dp), intent(in) :: efermi
    real(kind=dp), intent(in) :: w(npts)
    real(kind=dp), intent(out) :: spectrum(npts, 6)
    
    real(kind=dp) :: smear_value, delta_factor, f_mk, matcmpts(6)
    integer :: ik, iw, icmpt, im
    
    ! Smearing is half the step size so that every eigenvalue is essentially projected 
    ! onto the nearest spectrum point or two.
    smear_value = 0.5d0 * (w(2) - w(1))
    
    do ik=1,nkpts
      !print *, "Inside k = ", ik
      do im=1,ntrans
        !print *, "Inside m = ", im
        ! Simple cutoff test - eventually might add some kind of Fermi function here.
        print *, opttrans(ik, i_core, im), efermi - ecore
        if (opttrans(ik,i_core,im).gt.(efermi - ecore)) then
          f_mk = 1.0d0
          print *, ik, i_core, im
          print *, optmat(ik, i_core, im,:)
          matcmpts(1) = optmat(ik, i_core, im, 1)**2 + &
          &             optmat(ik, i_core, im, 2)**2
          matcmpts(2) = optmat(ik, i_core, im, 3)**2 + & 
          &             optmat(ik, i_core, im, 4)**2
          matcmpts(3) = optmat(ik, i_core, im, 5)**2 + &
          &             optmat(ik, i_core, im, 6)**2
          matcmpts(4) = 2.0d0 * (optmat(ik, i_core, im, 1) * &
          &                      optmat(ik, i_core, im, 3) + &
          &               optmat(ik, i_core, im, 2) * &
          &               optmat(ik, i_core, im, 4))
          matcmpts(5) = 2.0d0 * (optmat(ik, i_core, im, 1) * &
          &                      optmat(ik, i_core, im, 5) + &
          &               optmat(ik, i_core, im, 2) * &
          &               optmat(ik, i_core, im, 6))
          matcmpts(6) = 2.0d0 * (optmat(ik, i_core, im, 3) * &
          &                      optmat(ik, i_core, im, 5) + &
          &               optmat(ik, i_core, im, 4) * &
          &               optmat(ik, i_core, im, 6))
          
          do iw=1,npts
            ! Lorentzian projection factor.
            delta_factor = invpi * smear_value / ((opttrans(ik, i_core, im) - w(iw))**2 + smear_value**2)
            do icmpt=1,6
              spectrum(iw,icmpt) = spectrum(iw,icmpt) + f_mk * matcmpts(icmpt) * delta_factor
            end do
          end do
        else
          f_mk = 0.0d0
        endif
        
      end do
    end do
  end subroutine
end module