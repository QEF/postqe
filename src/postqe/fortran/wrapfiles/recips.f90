!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !! This routine generates the reciprocal lattice vectors \(b_1,b_2,b_3\)
  !! given the real space vectors \(a_1,a_2,a_3\). The \(b\)'s are units of
  !! \(2 \pi/a\).
  !
  use kinds, ONLY: DP
  implicit none
  real(DP) :: a1 (3)
  !! input: first direct lattice vector
  real(DP) :: a2 (3)
  !! input: second direct lattice vector
  real(DP) :: a3 (3)
  !! input: third direct lattice vector
  real(DP) :: b1 (3)
  !! output: first reciprocal lattice vector
  real(DP) :: b2 (3)
  !! output: second reciprocal lattice vector
  real(DP) :: b3 (3)
  !! output: third reciprocal lattice vector
  !
  ! ... local variables
  end subroutine recips
