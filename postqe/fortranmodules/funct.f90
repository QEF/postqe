!
! Copyright (C) 2004-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
!
! MODIFIED VERSION FOR F2PY:
! - only the driver routines are adapted to f2py (other soubroutines may be implemented in Python)
! - USE KINDS only DP substituted with INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
! - USE io_global only stdout substituted with   INTEGER :: stdout = 6    ! unit connected to standard output
! - lines with call to subroutine "errore" removed for now... error handling in calling Python?
! - exx_fraction, exx_started, finite_size_cell_volume, screening_parameter, gau_parameter should be input vars in these functions, they are not for now
! - In gcc_spin routine, variable "zeta" is intent inout and passed to subroutines calls. It must be check that f2py works fine with it.
! - lsd_lyp and LSD_GLYP functionals commented for now. The calls should be restored in this file and also added in an additional fortran file to be compiled with f2py
! Not implemented:
!------- NONLOCAL CORRECTIONS DRIVERS ----------------------------------
!------- META CORRECTIONS DRIVERS ----------------------------------
!------- DRIVERS FOR DERIVATIVES OF XC POTENTIAL -----------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-------  LDA DRIVERS --------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine xc (rho, iexch, icorr, igcx, igcc, imeta, inlc, ex, ec, vx, vc)
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
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  !     exx_fraction, exx_started, finite_size_cell_volume should be input vars
  real(DP):: exx_fraction = 0.0_DP
  logical :: exx_started = .false.
  real(DP):: finite_size_cell_volume = -1 ! not set

  real(DP), intent(in) :: rho
  real(DP), intent(out) :: ec, vc, ex, vx
  real(DP) :: ec__, vc__
  !
  real(DP), parameter :: small = 1.E-10_DP,  third = 1.0_DP / 3.0_DP, &
       pi34 = 0.6203504908994_DP  ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0_DP
     vc = 0.0_DP
     ex = 0.0_DP
     vx = 0.0_DP
     return
  else
     rs = pi34 / rho**third
     ! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  endif
  !..exchange
  if (iexch == 1) THEN             !  'sla'
     call slater (rs, ex, vx)
  ELSEIF (iexch == 2) THEN         !  'sl1'
     call slater1(rs, ex, vx)
  ELSEIF (iexch == 3) THEN         !  'rxc'
     CALL slater_rxc(rs, ex, vx)
  ELSEIF ((iexch == 4).or.(iexch==5)) THEN  ! 'oep','hf'
     IF (exx_started) then
        ex = 0.0_DP
        vx = 0.0_DP
     else
        call slater (rs, ex, vx)
     endif
  ELSEIF (iexch == 6) THEN         !  'pb0x'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = (1.0_DP - exx_fraction) * ex 
        vx = (1.0_DP - exx_fraction) * vx 
     end if
  ELSEIF (iexch == 7) THEN         !  'B3LYP'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = 0.8_DP * ex 
        vx = 0.8_DP * vx 
     end if
  ELSEIF (iexch == 8) THEN         !  'sla+kzk'
!     if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
!          'finite size corrected exchange used w/o initialization',1)
     call slaterKZK (rs, ex, vx, finite_size_cell_volume)
     !
  ELSEIF (iexch == 9) THEN         !  'X3LYP'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = 0.782_DP * ex 
        vx = 0.782_DP * vx 
     end if
  else
     ex = 0.0_DP
     vx = 0.0_DP
  endif
  !..correlation
  if (icorr == 1) then
     call pz (rs, 1, ec, vc)
  elseif (icorr == 2) then
     call vwn (rs, ec, vc)
  elseif (icorr == 3) then
     call lyp (rs, ec, vc)
  elseif (icorr == 4) then
     call pw (rs, 1, ec, vc)
  elseif (icorr == 5) then
     call wigner (rs, ec, vc)
  elseif (icorr == 6) then
     call hl (rs, ec, vc)
  elseif (icorr == 7) then
     call pz (rs, 2, ec, vc)
  elseif (icorr == 8) then
     call pw (rs, 2, ec, vc)
  elseif (icorr == 9) then
     call gl (rs, ec, vc)
  elseif (icorr ==10) then
!     if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
!          'finite size corrected correlation used w/o initialization',1)
     call pzKZK (rs, ec, vc, finite_size_cell_volume)
  elseif (icorr ==11) then
     call vwn1_rpa (rs, ec, vc)
  elseif (icorr ==12) then  ! 'B3LYP'
     call vwn (rs, ec, vc)
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     call lyp( rs, ec__, vc__ )
     ec = ec + 0.81_DP * ec__
     vc = vc + 0.81_DP * vc__
  elseif (icorr ==13) then  ! 'B3LYP-V1R'
     call vwn1_rpa (rs, ec, vc)
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     call lyp( rs, ec__, vc__ )
     ec = ec + 0.81_DP * ec__
     vc = vc + 0.81_DP * vc__
  elseif (icorr ==14) then  ! 'X3LYP'
     call vwn1_rpa (rs, ec, vc)
     ec = 0.129_DP * ec
     vc = 0.129_DP * vc
     call lyp( rs, ec__, vc__ )
     ec = ec + 0.871_DP * ec__
     vc = vc + 0.871_DP * vc__
  else
     ec = 0.0_DP
     vc = 0.0_DP
  endif
  !
  return
end subroutine xc
!!!!!!!!!!!!!!SPIN
!-----------------------------------------------------------------------
subroutine xc_spin (rho, zeta, iexch, icorr, igcx, igcc, imeta, inlc, ex, ec, vxup, vxdw, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  real(DP), intent(in) :: rho, zeta
  real(DP), intent(out) :: ex, ec, vxup, vxdw, vcup, vcdw

  real(DP) :: ec__, vcup__, vcdw__

  !     exx_fraction, exx_started, finite_size_cell_volume should be input vars
  real(DP):: exx_fraction = 0.0_DP
  logical :: exx_started = .false.
  real(DP):: finite_size_cell_volume = -1 ! not set

  !
  real(DP), parameter :: small= 1.E-10_DP, third = 1.0_DP/3.0_DP, &
       pi34= 0.6203504908994_DP ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0_DP
     vcup = 0.0_DP
     vcdw = 0.0_DP
     ex = 0.0_DP
     vxup = 0.0_DP
     vxdw = 0.0_DP
     return
  else
     rs = pi34 / rho**third
  endif
  !..exchange
  IF (iexch == 1) THEN      ! 'sla'
     call slater_spin (rho, zeta, ex, vxup, vxdw)
  ELSEIF (iexch == 2) THEN  ! 'sl1'
     call slater1_spin (rho, zeta, ex, vxup, vxdw)
  ELSEIF (iexch == 3) THEN  ! 'rxc'
     call slater_rxc_spin ( rho, zeta, ex, vxup, vxdw )
  ELSEIF ((iexch == 4).or.(iexch==5)) THEN  ! 'oep','hf'
     IF (exx_started) then
        ex   = 0.0_DP
        vxup = 0.0_DP 
        vxdw = 0.0_DP 
     else
        call slater_spin (rho, zeta, ex, vxup, vxdw)
     endif
  ELSEIF (iexch == 6) THEN  ! 'pb0x'
     call slater_spin (rho, zeta, ex, vxup, vxdw)
     if (exx_started) then
        ex   = (1.0_DP - exx_fraction) * ex
        vxup = (1.0_DP - exx_fraction) * vxup 
        vxdw = (1.0_DP - exx_fraction) * vxdw 
     end if
  ELSEIF (iexch == 7) THEN  ! 'B3LYP'
     call slater_spin (rho, zeta, ex, vxup, vxdw)
     if (exx_started) then
        ex   = 0.8_DP * ex
        vxup = 0.8_DP * vxup 
        vxdw = 0.8_DP * vxdw 
     end if
  ELSE
     ex = 0.0_DP
     vxup = 0.0_DP
     vxdw = 0.0_DP
  ENDIF
  !..correlation
  if (icorr == 0) then
     ec = 0.0_DP
     vcup = 0.0_DP
     vcdw = 0.0_DP
  elseif (icorr == 1) then
     call pz_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 2) then
     call vwn_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 3) then
!     call lsd_lyp (rho, zeta, ec, vcup, vcdw) ! from CP/FPMD (more_functionals)
  elseif (icorr == 4) then
     call pw_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 12) then ! 'B3LYP'
     call vwn_spin (rs, zeta, ec, vcup, vcdw)
     ec = 0.19_DP * ec
     vcup = 0.19_DP * vcup
     vcdw = 0.19_DP * vcdw
!     call lsd_lyp (rho, zeta, ec__, vcup__, vcdw__) ! from CP/FPMD (more_functionals)
     ec = ec + 0.81_DP * ec__
     vcup = vcup + 0.81_DP * vcup__
     vcdw = vcdw + 0.81_DP * vcdw__
  elseif (icorr == 13) then   ! 'B3LYP-V1R'
     call vwn1_rpa_spin (rs, zeta, ec, vcup, vcdw)
     ec = 0.19_DP * ec
     vcup = 0.19_DP * vcup
     vcdw = 0.19_DP * vcdw
!     call lsd_lyp (rho, zeta, ec__, vcup__, vcdw__) ! from CP/FPMD (more_functionals)
     ec = ec + 0.81_DP * ec__
     vcup = vcup + 0.81_DP * vcup__
     vcdw = vcdw + 0.81_DP * vcdw__
!  else
!     call errore ('lsda_functional (xc_spin)', 'not implemented', icorr)
  endif
  !
  return
end subroutine xc_spin
!-----------------------------------------------------------------------
subroutine xc_spin_vec (rho, zeta, length, iexch, icorr, igcx, igcc, imeta, inlc, evx, evc)
  !-----------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  implicit none

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  integer, intent(in)   :: length
  real(DP), intent(in)  :: rho(length), zeta(length)
  real(DP), intent(out) :: evx(length,3), evc(length,3)
  !
  !     exx_fraction, exx_started, finite_size_cell_volume should be input vars
  real(DP):: exx_fraction = 0.0_DP
  logical :: exx_started = .false.
  real(DP):: finite_size_cell_volume = -1 ! not set


  real(DP), parameter :: small= 1.E-10_DP, third = 1.0_DP/3.0_DP, &
       pi34= 0.6203504908994_DP ! pi34=(3/4pi)^(1/3)
  !
  integer  :: i
  logical  :: comp_energy_loc
  real(DP) :: rs(length)
  !
  !..exchange
  select case (iexch)
  case(1)            ! 'sla'
     call slater_spin_vec (rho, zeta, evx, length)
  case(2)            ! 'sl1'
     do i=1,length
        call slater1_spin (rho(i), zeta(i), evx(i,3), evx(i,1), evx(i,2))
     end do
  case(3)            ! 'rxc'
     do i=1,length
        call slater_rxc_spin (rho(i), zeta(i), evx(i,3), evx(i,1), evx(i,2))
     end do
  case(4,5)          ! 'oep','hf'
     if (exx_started) then
        evx = 0.0_DP
     else
        call slater_spin_vec (rho, zeta, evx, length)
     endif
  case(6)            ! 'pb0x'
     call slater_spin_vec (rho, zeta, evx, length)
     if (exx_started) then
        evx = (1.0_DP - exx_fraction) * evx
     end if
  case(7)            ! 'B3LYP'
     call slater_spin_vec (rho, zeta, evx, length)
     if (exx_started) then
        evx = 0.8_DP * evx
     end if
  case default
     evx = 0.0_DP
  end select

  !..correlation
  where (rho.gt.small)
     rs = pi34 / rho**third
  elsewhere
     rs = 1.0_DP ! just a sane default, results are discarded anyway
  end where

  select case(icorr)
  case (0)
     evc = 0.0_DP
  case (1)
     do i=1,length
        call pz_spin (rs(i), zeta(i), evc(i,3), evc(i,1), evc(i,2))
     end do
  case (2)
     do i=1,length
        call vwn_spin (rs(i), zeta(i), evc(i,3), evc(i,1), evc(i,2))
     end do
!  case(3)
!     do i=1,length
!        call lsd_lyp (rho(i), zeta(i), evc(i,3), evc(i,1), evc(i,2)) ! from CP/FPMD (more_functionals)
!     end do
  case(4)
     call pw_spin_vec (rs, zeta, evc, length)
!  case default
!     call errore ('lsda_functional (xc_spin_vec)', 'not implemented', icorr)
  end select
  !
  where (rho.le.small)
     evx(:,1) = 0.0_DP
     evc(:,1) = 0.0_DP

     evx(:,2) = 0.0_DP
     evc(:,2) = 0.0_DP

     evx(:,3) = 0.0_DP
     evc(:,3) = 0.0_DP
  end where
  !
end subroutine xc_spin_vec
!
!-----------------------------------------------------------------------
!------- GRADIENT CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine gcxc (rho, grho, iexch, icorr, igcx, igcc, imeta, inlc, sx, sc, v1x, v2x, v1c, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     See comments at the beginning of module for implemented cases
  !
  !     input:  rho, grho=|\nabla rho|^2
  !     definition:  E_x = \int E_x(rho,grho) dr
  !     output: sx = E_x(rho,grho)
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             sc, v1c, v2c as above for correlation
  !
  implicit none

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  real(DP), intent(in) :: rho, grho
  real(DP), intent(out) :: sx, sc, v1x, v2x, v1c, v2c
  !
  !     exx_fraction, exx_started, finite_size_cell_volume, screening_parameter, gau_parameter should be input vars
  real(DP):: exx_fraction = 0.0_DP
  logical :: exx_started = .false.
  real(DP):: finite_size_cell_volume = -1 ! not set
  real(DP)::  screening_parameter = 0.0_DP
  real(DP)::  gau_parameter = 0.0_DP

  real(DP) :: sx__,v1x__, v2x__
  real(DP) :: sxsr, v1xsr, v2xsr
  real(DP), parameter:: small = 1.E-10_DP

  ! exchange
  if (rho <= small) then
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  elseif (igcx == 1) then
     call becke88 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 2) then
     call ggax (rho, grho, sx, v1x, v2x)
  elseif (igcx == 3) then
     call pbex (rho, grho, 1, sx, v1x, v2x)
  elseif (igcx == 4) then
     call pbex (rho, grho, 2, sx, v1x, v2x)
  elseif (igcx == 5 .and. igcc == 5) then
     call hcth(rho, grho, sx, v1x, v2x)
  elseif (igcx == 6) then
     call optx (rho, grho, sx, v1x, v2x)
  ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
  ! routine: needs kinetic energy density in addition to rho and grad rho
  elseif (igcx == 8) then ! 'PBE0'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if (exx_started) then
        sx  = (1.0_DP - exx_fraction) * sx
        v1x = (1.0_DP - exx_fraction) * v1x
        v2x = (1.0_DP - exx_fraction) * v2x
     end if
  elseif (igcx == 9) then ! 'B3LYP'
     call becke88 (rho, grho, sx, v1x, v2x)
     if (exx_started) then
        sx  = 0.72_DP * sx
        v1x = 0.72_DP * v1x
        v2x = 0.72_DP * v2x
     end if
  elseif (igcx ==10) then ! 'pbesol'
     call pbex (rho, grho, 3, sx, v1x, v2x)
  elseif (igcx ==11) then ! 'wc'
     call wcx (rho, grho, sx, v1x, v2x)
  elseif (igcx ==12) then ! 'pbexsr'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if(exx_started) then
       call pbexsr (rho, grho, sxsr, v1xsr, v2xsr, screening_parameter)
       sx = sx - exx_fraction * sxsr
       v1x = v1x - exx_fraction * v1xsr
       v2x = v2x - exx_fraction * v2xsr
     endif 
  elseif (igcx ==13) then ! 'rPW86'
     call rPW86 (rho, grho, sx, v1x, v2x)
  elseif (igcx ==16) then ! 'C09x'
     call c09x (rho, grho, sx, v1x, v2x)
  elseif (igcx ==17) then ! 'sogga'
     call sogga(rho, grho, sx, v1x, v2x)
  elseif (igcx ==19) then ! 'pbeq2d'
     call pbex (rho, grho, 4, sx, v1x, v2x)
  elseif (igcx ==20) then ! 'gau-pbe'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if(exx_started) then
       call pbexgau (rho, grho, sxsr, v1xsr, v2xsr, gau_parameter)
       sx = sx - exx_fraction * sxsr
       v1x = v1x - exx_fraction * v1xsr
       v2x = v2x - exx_fraction * v2xsr
     endif
  elseif (igcx == 21) then ! 'pw86'
     call pw86 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 22) then ! 'b86b'
     call becke86b (rho, grho, sx, v1x, v2x)
     ! call b86b (rho, grho, 1, sx, v1x, v2x)
  elseif (igcx == 23) then ! 'optB88'
     call pbex (rho, grho, 5, sx, v1x, v2x)
  elseif (igcx == 24) then ! 'optB86b'
     call pbex (rho, grho, 6, sx, v1x, v2x)
     ! call b86b (rho, grho, 2, sx, v1x, v2x)
  elseif (igcx == 25) then ! 'ev93'
     call pbex (rho, grho, 7, sx, v1x, v2x)
  elseif (igcx == 26) then ! 'b86r'
     call b86b (rho, grho, 3, sx, v1x, v2x)
  elseif (igcx == 27) then ! 'cx13'
     call cx13 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 28) then ! 'X3LYP'
     call becke88 (rho, grho, sx, v1x, v2x)
     call pbex (rho, grho, 1, sx__, v1x__, v2x__)
     if (exx_started) then
        sx  = real(0.765*0.709,DP) * sx
        v1x = real(0.765*0.709,DP) * v1x
        v2x = real(0.765*0.709,DP) * v2x
        sx  = sx  + real(0.235*0.709) * sx__
        v1x = v1x + real(0.235*0.709) * v1x__
        v2x = v2x + real(0.235*0.709) * v2x__
     end if
  else
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  endif
  ! correlation
  if (rho.le.small) then
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  elseif (igcc == 1) then
     call perdew86 (rho, grho, sc, v1c, v2c)
  elseif (igcc == 2) then
     call ggac (rho, grho, sc, v1c, v2c)
  elseif (igcc == 3) then
     call glyp (rho, grho, sc, v1c, v2c)
  elseif (igcc == 4) then
     call pbec (rho, grho, 1, sc, v1c, v2c)
  ! igcc == 5 (HCTH) is calculated together with case igcx=5
  ! igcc == 6 (meta-GGA) is treated in a different routine
  elseif (igcc == 7) then !'B3LYP'
     call glyp (rho, grho, sc, v1c, v2c)
     if (exx_started) then
        sc  = 0.81_DP * sc
        v1c = 0.81_DP * v1c
        v2c = 0.81_DP * v2c
     end if
  elseif (igcc == 8) then ! 'PBEsol'
     call pbec (rho, grho, 2, sc, v1c, v2c)
  ! igcc == 9 set to 5, back-compatibility
  ! igcc ==10 set to 6, back-compatibility
  ! igcc ==11 M06L calculated in another routine
  else if (igcc == 12) then ! 'Q2D'
     call pbec (rho, grho, 3, sc, v1c, v2c)
  elseif (igcc == 13) then !'X3LYP'
     call glyp (rho, grho, sc, v1c, v2c)
     if (exx_started) then
        sc  = 0.871_DP * sc
        v1c = 0.871_DP * v1c
        v2c = 0.871_DP * v2c
     end if
  else
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  endif
  !
  return
end subroutine gcxc
!
!
!!!!!!!!!!!!!!SPIN
!-----------------------------------------------------------------------
subroutine gcx_spin (rhoup, rhodw, grhoup2, grhodw2, &
                     iexch, icorr, igcx, igcc, imeta, inlc, &
                     sx, v1xup, v1xdw, v2xup, v2xdw)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !
  implicit none

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  !
  !     dummy arguments
  !
  real(DP), intent(in) :: rhoup, rhodw, grhoup2, grhodw2
  real(DP), intent(out) :: sx, v1xup, v1xdw, v2xup, v2xdw
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  !
  !     exx_fraction, exx_started, finite_size_cell_volume, screening_parameter, gau_parameter should be input vars
  real(DP):: exx_fraction = 0.0_DP
  logical :: exx_started = .false.
  real(DP):: finite_size_cell_volume = -1 ! not set
  real(DP)::  screening_parameter = 0.0_DP
  real(DP)::  gau_parameter = 0.0_DP



  real(DP) :: sxsr, v1xupsr, v2xupsr, v1xdwsr, v2xdwsr
  real(DP), parameter :: small = 1.E-10_DP
  real(DP) :: rho, sxup, sxdw
  integer :: iflag
  !
  !
  ! exchange
  rho = rhoup + rhodw
  if (rho <= small .or. igcx == 0) then
     sx = 0.0_DP
     v1xup = 0.0_DP
     v2xup = 0.0_DP
     v1xdw = 0.0_DP
     v2xdw = 0.0_DP
  elseif (igcx == 1) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = sxup + sxdw
  elseif (igcx == 2) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call ggax (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call ggax (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
  elseif (igcx == 3 .or. igcx == 4 .or. igcx == 8 .or. &
          igcx == 10 .or. igcx == 12 .or. igcx == 20 .or. igcx == 25) then
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8: PBE0, igcx=10: PBEsol
     ! igcx=12: HSE,  igcx=20: gau-pbe, igcx=25: ev93
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
     elseif (igcx == 25) then
        iflag = 7
     else
        iflag = 1
     endif
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pbex (2.0_DP * rhoup, 4.0_DP * grhoup2, iflag, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pbex (2.0_DP * rhodw, 4.0_DP * grhodw2, iflag, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
     if (igcx == 8 .and. exx_started ) then
       sx = (1.0_DP - exx_fraction) * sx
       v1xup = (1.0_DP - exx_fraction) * v1xup
       v1xdw = (1.0_DP - exx_fraction) * v1xdw
       v2xup = (1.0_DP - exx_fraction) * v2xup
       v2xdw = (1.0_DP - exx_fraction) * v2xdw
     end if
     if (igcx == 12 .and. exx_started ) then

        call pbexsr_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
                         v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                         screening_parameter)
        sx  = sx - exx_fraction*sxsr
        v1xup = v1xup - exx_fraction*v1xupsr
        v2xup = v2xup - exx_fraction*v2xupsr
        v1xdw = v1xdw - exx_fraction*v1xdwsr
        v2xdw = v2xdw - exx_fraction*v2xdwsr
     end if

     if (igcx == 20 .and. exx_started ) then
        ! gau-pbe
        call pbexgau_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
                         v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                         gau_parameter)
        sx  = sx - exx_fraction*sxsr
        v1xup = v1xup - exx_fraction*v1xupsr
        v2xup = v2xup - exx_fraction*v2xupsr
        v1xdw = v1xdw - exx_fraction*v1xdwsr
        v2xdw = v2xdw - exx_fraction*v2xdwsr
     end if

  elseif (igcx == 9) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = sxup + sxdw

     if (exx_started ) then
       sx = 0.72_DP * sx
       v1xup = 0.72_DP * v1xup
       v1xdw = 0.72_DP * v1xdw
       v2xup = 0.72_DP * v2xup
       v2xdw = 0.72_DP * v2xdw
     end if

  elseif (igcx == 11) then ! 'Wu-Cohen'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call wcx (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call wcx (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 13) then ! 'revised PW86 for vdw-df2'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call rPW86 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call rPW86 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 16) then ! 'c09x for vdw-df-c09.'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call c09x (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call c09x (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 21) then ! 'PW86'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pw86 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pw86 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 22) then ! 'B86B'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke86b (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke86b (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

   elseif (igcx == 26) then ! 'B86R for rev-vdW-DF2'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call b86b (2.0_DP * rhoup, 4.0_DP * grhoup2, 3, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call b86b (2.0_DP * rhodw, 4.0_DP * grhodw2, 3, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 27) then ! 'cx13 for vdw-df-cx'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call cx13 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call cx13 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw


  ! case igcx == 5 (HCTH) and 6 (OPTX) not implemented
  ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
  ! routine: needs kinetic energy density in addition to rho and grad rho

!  else
!     call errore ('gcx_spin', 'not implemented', igcx)
  endif
  !
  return
end subroutine gcx_spin
!
!-----------------------------------------------------------------------
subroutine gcx_spin_vec(rhoup, rhodw, grhoup2, grhodw2, &
     iexch, icorr, igcx, igcc, imeta, inlc, &
     sx, v1xup, v1xdw, v2xup, v2xdw, length)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !
  implicit none

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  !
  !     dummy arguments
  !
  integer, intent(in) :: length
  real(DP),intent(in) :: rhoup(length), rhodw(length)
  real(DP),intent(in) :: grhoup2(length), grhodw2(length)
  real(DP),intent(out) :: sx(length)
  real(DP),intent(out) :: v1xup(length), v1xdw(length)
  real(DP),intent(out) :: v2xup(length), v2xdw(length)
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !

  !     exx_fraction, exx_started, finite_size_cell_volume, screening_parameter, gau_parameter should be input vars
  real(DP):: exx_fraction = 0.0_DP
  logical :: exx_started = .false.
  real(DP):: finite_size_cell_volume = -1 ! not set
  real(DP)::  screening_parameter = 0.0_DP
  real(DP)::  gau_parameter = 0.0_DP


  real(DP), parameter :: small = 1.E-10_DP
  real(DP) :: rho(length), sxup(length), sxdw(length)
  integer :: iflag
  integer :: i
  ! only used for HSE (igcx == 12):
  real(DP) :: sxsr, v1xupsr, v2xupsr, v1xdwsr, v2xdwsr
  !
  !
  ! exchange
  rho = rhoup + rhodw
  select case(igcx)
  case(0)
     sx = 0.0_DP
     v1xup = 0.0_DP
     v2xup = 0.0_DP
     v1xdw = 0.0_DP
     v2xdw = 0.0_DP
  case(1)
     do i=1,length
        if (rhoup(i) > small .and. sqrt (abs (grhoup2(i)) ) > small) then
           call becke88_spin (rhoup(i), grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt (abs (grhodw2(i)) ) > small) then
           call becke88_spin (rhodw(i), grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = sxup + sxdw
  case(2)
     do i=1,length
        if (rhoup(i) > small .and. sqrt (abs (grhoup2(i)) ) > small) then
           call ggax (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt (abs (grhodw2(i)) ) > small) then
           call ggax (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
  case(3,4,8,10,12,25)
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8 PBE0, igcx=10: PBEsol, 
     ! igcx=25: EV93
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
     elseif (igcx == 25) then
        iflag = 7
     else
        iflag = 1
     endif

     call pbex_vec (2.0_DP * rhoup, 4.0_DP * grhoup2, iflag, sxup, v1xup, v2xup, length, small)
     call pbex_vec (2.0_DP * rhodw, 4.0_DP * grhodw2, iflag, sxdw, v1xdw, v2xdw, length, small)
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
     if (igcx == 8 .and. exx_started ) then
        sx = (1.0_DP - exx_fraction) * sx
        v1xup = (1.0_DP - exx_fraction) * v1xup
        v1xdw = (1.0_DP - exx_fraction) * v1xdw
        v2xup = (1.0_DP - exx_fraction) * v2xup
        v2xdw = (1.0_DP - exx_fraction) * v2xdw
     end if
     if (igcx == 12 .and. exx_started ) then
        ! in this case the subroutine is not really "vector"
        DO i = 1, length
          call pbexsr_lsd (rhoup(i), rhodw(i), grhoup2(i), grhodw2(i), sxsr,  &
                           v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                           screening_parameter)
        sx(i)  = sx(i) - exx_fraction*sxsr
        v1xup(i) = v1xup(i) - exx_fraction*v1xupsr
        v2xup(i) = v2xup(i) - exx_fraction*v2xupsr
        v1xdw(i) = v1xdw(i) - exx_fraction*v1xdwsr
        v2xdw(i) = v2xdw(i) - exx_fraction*v2xdwsr
        ENDDO
     end if
  case(9)
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i)) ) > small) then
           call becke88_spin (rhoup(i), grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call becke88_spin (rhodw(i), grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = sxup + sxdw

     if (exx_started ) then
        sx = 0.72_DP * sx
        v1xup = 0.72_DP * v1xup
        v1xdw = 0.72_DP * v1xdw
        v2xup = 0.72_DP * v2xup
        v2xdw = 0.72_DP * v2xdw
     end if

  case(11) ! 'Wu-Cohen'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call wcx (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call wcx (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(13) ! 'rPW86 for vdw-df2'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call rPW86 (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call rPW86 (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(16) ! 'c09x for vdw-df-c09'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call c09x (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call c09x (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(21) ! 'pw86'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call pw86 (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call pw86 (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(22) ! 'b86b'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call becke86b (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call becke86b (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(26) ! 'B86R for rev-vdW-DF2'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call b86b (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), 3, sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call b86b (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), 3, sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(27) ! 'cx13 for vdw-df-cx'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call cx13 (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call cx13 (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

!  case default
!     call errore ('gcx_spin_vec', 'not implemented', igcx)
  end select
  !
  if (igcx.ne.0) then
     where (rho.le.small)
        sx = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     end where
  end if
  !
end subroutine gcx_spin_vec
!
!
!-----------------------------------------------------------------------
subroutine gcc_spin (rho, zeta, grho, iexch, icorr, igcx, igcc, imeta, inlc, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for correlations - Hartree a.u.
  !     Implemented:  Perdew86, GGA (PW91), PBE
  !
  implicit none

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
  integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

  !
  !     dummy arguments
  !
  real(DP), intent(in) :: rho, grho
  real(DP), intent(inout) :: zeta
!f2py intent(inout) zeta 
  real(DP), intent(out) :: sc, v1cup, v1cdw, v2c
  ! the total charge
  ! the magnetization
  ! the gradient of the charge squared
  ! exchange and correlation energies
  ! derivatives of correlation wr. rho
  ! derivatives of correlation wr. grho

  real(DP), parameter :: small = 1.E-10_DP, epsr=1.E-6_DP
  !
  if ( abs(zeta) > 1.0_DP ) then
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     v2c = 0.0_DP
     return
  else
     !
     ! ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
     zeta = SIGN( MIN( ABS( zeta ), ( 1.0_DP - epsr ) ) , zeta )
  endif

  if (igcc == 0 .or. rho <= small .or. sqrt(abs(grho)) <= small) then
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     v2c = 0.0_DP
  elseif (igcc == 1) then
     call perdew86_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 2) then
     call ggac_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 4) then
     call pbec_spin (rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 8) then
     call pbec_spin (rho, zeta, grho, 2, sc, v1cup, v1cdw, v2c)
!  else
!     call errore ('lsda_functionals (gcc_spin)', 'not implemented', igcc)
  endif
  !
  return
end subroutine gcc_spin
!
!   ==================================================================
    SUBROUTINE gcc_spin_more( RHOA, RHOB, GRHOAA, GRHOBB, GRHOAB, &
                              iexch, icorr, igcx, igcc, imeta, inlc, &
                              SC, V1CA, V1CB, V2CA, V2CB, V2CAB )
!   ==--------------------------------------------------------------==
!   ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
!   ==                                                              ==
!   ==  EXCHANGE  :  BECKE88                                        ==
!   ==               GGAX                                           ==
!   ==  CORRELATION : PERDEW86                                      ==
!   ==                LEE, YANG & PARR                              ==
!   ==                GGAC                                          ==
!   ==--------------------------------------------------------------==

      IMPLICIT NONE
      INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  !
  ! indices for exchange-correlation given in input, after input variables as
  ! in the original routine. In this way, the python code can use a *list to given
  ! them in one list (after all other input variables)
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    inlc:  type of meta-GGA
      integer, intent(in) :: iexch, icorr, igcx, igcc, imeta, inlc

      REAL(DP), intent(in) :: RHOA,RHOB,GRHOAA,GRHOBB,GRHOAB
      REAL(DP), intent(out) :: SC,V1CA,V2CA,V1CB,V2CB,V2CAB

      !     exx_fraction, exx_started
      real(DP):: exx_fraction = 0.0_DP
      logical :: exx_started = .false.

      ! ... Gradient Correction for correlation

      REAL(DP) :: SMALL, RHO
      PARAMETER(SMALL=1.E-20_DP)

      SC=0.0_DP
      V1CA=0.0_DP
      V2CA=0.0_DP
      V1CB=0.0_DP
      V2CB=0.0_DP
      V2CAB=0.0_DP
      IF( igcc == 3 .or. igcc == 7) THEN
        RHO=RHOA+RHOB
        IF(RHO.GT.SMALL) then
!             CALL LSD_GLYP(RHOA,RHOB,GRHOAA,GRHOAB,GRHOBB,SC,&
!                  V1CA,V2CA,V1CB,V2CB,V2CAB)
             if (igcc == 7 .and. exx_started) then
                SC = 0.81d0*SC
                V1CA = 0.81d0*V1CA
                V2CA = 0.81d0*V2CA
                V1CB = 0.81d0*V1CB
                V2CB = 0.81d0*V2CB
                V2CAB = 0.81d0*V2CAB
             endif
         endif
!      ELSE
!        CALL errore( " gcc_spin_more ", " gradiet correction not implemented ", 1 )
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE gcc_spin_more
!
!
