!
! Wrappers for functionals and related subroutines 
!


!-----------------------------------------------------------------------
!subroutine py_set_dft_from_name( pydft_ , pyiexch, pyicorr, pyigcx, pyigcc, pyimeta, pyinlc )
!-----------------------------------------------------------------------
!
! translates a string containing the exchange-correlation name
! into internal indices iexch, icorr, igcx, igcc
!
  !USE funct,     ONLY: set_dft_from_name, iexch, icorr, igcx, igcc, imeta, inlc
!  USE funct,     ONLY: set_dft_from_name
!  USE funct,     ONLY: iexch, icorr, igcx, igcc, imeta, inlc
!  implicit none
!  character(len=*), intent(in) :: pydft_

!  integer, intent(out) :: pyiexch, pyicorr, pyigcx, pyigcc, pyimeta, pyinlc


!  CALL set_dft_from_name( pydft_ )
  !pyiexch = iexch
  !
!  return
!end subroutine py_set_dft_from_name

!-----------------------------------------------------------------------
subroutine pyqe_xc (rho, functional, ex, ec, vx, vc)
  !-----------------------------------------------------------------------
  !     lda exchange and correlation functionals - Hartree a.u.
  !
  !     exchange   :  Slater, relativistic Slater
  !     correlation:  Ceperley-Alder (Perdew-Zunger parameters)
  !                   Vosko-Wilk-Nusair
  !                   Lee-Yang-Parr
  !                   Perdew-Wang
  !                   Wigner
  !                   Hedin-Lundqvist
  !                   Ortiz-Ballone (Perdew-Zunger formula)
  !                   Ortiz-Ballone (Perdew-Wang formula)
  !                   Gunnarsson-Lundqvist
  !
  !     input : rho=rho(r)
  !     definitions: E_x = \int E_x(rho) dr, E_x(rho) = rho\epsilon_c(rho)
  !                  same for correlation
  !     output: ex = \epsilon_x(rho) ( NOT E_x(rho) )
  !             vx = dE_x(rho)/drho  ( NOT d\epsilon_x(rho)/drho )
  !             ec, vc as above for correlation
  !
  USE funct,     ONLY: xc, set_dft_from_name, get_iexch, get_icorr, get_igcx
  implicit none

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), intent(in) :: rho
  CHARACTER(len=*), intent(in) :: functional
  REAL(DP), intent(out) :: ec, vc, ex, vx
  integer :: test_iexch = -1
  integer :: test_icorr = -1
  integer :: test_igcx = -1
  !

  CALL set_dft_from_name( functional )
  test_iexch = get_iexch()
  test_icorr = get_icorr()
  test_igcx = get_igcx ()
  !PRINT *, test_iexch, test_icorr, test_igcx
  CALL xc (rho, ex, ec, vx, vc)
  !
  return
end subroutine pyqe_xc
