!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE latgen_lib(ibrav,celldm,a1,a2,a3,omega, ierr, errormsg)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111)
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ibrav
  real(DP), INTENT(inout) :: celldm(6)
  real(DP), INTENT(inout) :: a1(3), a2(3), a3(3)
  real(DP), INTENT(out) :: omega
  !
  character(len=*),INTENT(out) :: errormsg
  integer,INTENT(out) :: ierr
  !
  real(DP), PARAMETER:: sr2 = 1.414213562373d0, &
                        sr3 = 1.732050807569d0
  INTEGER :: i,j,k,l,iperm,ir
  real(DP) :: term, cbya, s, term1, term2, singam, sen
  !
  !
END SUBROUTINE latgen_lib
!
!-------------------------------------------------------------------------
SUBROUTINE at2celldm (ibrav,alat,a1,a2,a3,celldm)
  !-----------------------------------------------------------------------
  !
  !     Returns celldm parameters computed from lattice vectors a1,a2,a3 
  !     a1, a2, a3 are in "alat" units
  !     If Bravais lattice index ibrav=0, only celldm(1) is set to alat
  !     See latgen for definition of celldm and lattice vectors.
  !     a1, a2, a3, ibrav, alat are not modified
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ibrav
  REAL(DP), INTENT(in) :: alat, a1(3), a2(3), a3(3)
  REAL(DP), INTENT(out) :: celldm(6)
  !
  !
END SUBROUTINE at2celldm
!
FUNCTION at2ibrav (a1, a2, a3) RESULT (ibrav)
  !
  !     Returns ibrav from lattice vectors if recognized, 0 otherwise
  !
  USE kinds, ONLY: dp
  IMPLICIT NONE
  REAL(dp), INTENT (in) :: a1(3), a2(3), a3(3)
  REAL(dp) :: v1, v2, v3, cosab, cosac, cosbc
  !
  INTEGER :: ibrav
  !
  !
END FUNCTION at2ibrav
!
SUBROUTINE abc2celldm ( ibrav, a,b,c,cosab,cosac,cosbc, celldm )
  !
  !  returns internal parameters celldm from crystallographics ones
  !
  USE kinds,     ONLY: dp
  USE constants, ONLY: bohr_radius_angs
  IMPLICIT NONE
  !
  INTEGER,  INTENT (in) :: ibrav
  REAL(DP), INTENT (in) :: a,b,c, cosab, cosac, cosbc
  REAL(DP), INTENT (out) :: celldm(6)
  !
  !
END SUBROUTINE abc2celldm
!
SUBROUTINE celldm2abc ( ibrav, celldm, a,b,c,cosab,cosac,cosbc )
  !
  !  returns crystallographic parameters a,b,c from celldm
  !
  USE kinds,     ONLY: dp
  USE constants, ONLY: bohr_radius_angs
  IMPLICIT NONE
  !
  INTEGER,  INTENT (in) :: ibrav
  REAL(DP), INTENT (in) :: celldm(6)
  REAL(DP), INTENT (out) :: a,b,c, cosab, cosac, cosbc
  !
  !
  !
END SUBROUTINE celldm2abc

SUBROUTINE remake_cell(ibrav, alat, a1,a2,a3, new_alat)
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER,INTENT(in) :: ibrav
  REAL(DP),INTENT(in)  :: alat
  REAL(DP),INTENT(out) :: new_alat
  REAL(DP),INTENT(inout) :: a1(3),a2(3),a3(3)
  REAL(DP) :: e1(3), e2(3), e3(3)
  REAL(DP) :: celldm_internal(6), lat_internal, omega
  ! Better not to do the following, or it may cause problems with ibrav=0 from input
!  ibrav = at2ibrav (a(:,1), a(:,2), a(:,3))
  ! Instead, let's print a warning and do nothing:

END SUBROUTINE

!-------------------------------------------------------------------------
SUBROUTINE latgen(ibrav,celldm,a1,a2,a3,omega)
  !-----------------------------------------------------------------------
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ibrav
  real(DP), INTENT(inout) :: celldm(6)
  real(DP), INTENT(inout) :: a1(3), a2(3), a3(3)
  real(DP), INTENT(out) :: omega
  !
  !
END SUBROUTINE 

