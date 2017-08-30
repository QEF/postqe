!
! Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
! Internazionale Superiore di Studi Avanzati). All rights reserved.
! This file is distributed under the terms of the LGPL-2.1 license. See the
! file 'LICENSE' in the root directory of the present distribution, or
! https://opensource.org/licenses/LGPL-2.1
!
!----------------------------------------------------------------------
subroutine struc_fact2(nat, tau, ngm, g, strf ,check_gg, check_tau)
  !----------------------------------------------------------------------
  !
  !   calculate the structure factors for each type of atoms in the unit
  !   cell
  !
  
  implicit none
  !
  !   Here the dummy variables
  !
  integer, parameter  :: dp = selected_real_kind(14,200) 
  real(dp), parameter :: tpi = 8.d0* atan ( 1.d0) 
  integer, intent(in) :: nat, ngm 
  ! input: the number of identical atom in the unit cell
  ! input: the number of G vectors
  ! input: fft dimension along x
  ! input: fft dimension along y
  ! input: fft dimension along z
  
  real(DP),intent ( in)  :: tau (nat, 3 ), g (ngm, 3 )
  ! input: reciprocal crystal basis vectors
  ! input: the positions of the atoms in the c
  ! input: the coordinates of the g vectors

  complex(DP),intent(out)  :: strf (ngm ) 
  real(DP), intent ( out)  :: check_gg(ngm) 
  real(DP), intent ( out)  :: check_tau(3*nat) 
  ! output: the structure factor
  !
  !    here the local variables
  !
  integer :: nt, na, ng, n1, n2, n3, ipol
  ! counter over atom type
  ! counter over atoms
  ! counter over G vectors
  ! counter over fft dimension along x
  ! counter over fft dimension along y
  ! counter over fft dimension along z
  ! counter over polarizations
  
  real(DP) :: arg
  ! the argument of the exponent
  ! scalar product of bg and tau
  do na = 1, nat 
     do ipol =1 ,3 
        check_tau(3*(nat-1)+ipol ) = tau(na,ipol)
     end do
  end do
  strf(:) = (0.d0,0.d0)
  do na = 1, nat
        do ng = 1, ngm
           check_gg(ng) = sum(g(ng,:)*g(ng,:))
           arg = (g (ng, 1) * tau (na, 1) + g (ng, 2) * tau (na, 2) &
                 + g (ng, 3) * tau (na, 3) )*tpi  
           strf (ng ) = strf (ng ) + CMPLX(cos (arg), -sin (arg),kind=DP)
        enddo
  enddo

  return
end subroutine struc_fact2

