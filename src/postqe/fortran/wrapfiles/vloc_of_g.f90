!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE vloc_of_g( mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, &
                      gl, omega, vloc )
  !----------------------------------------------------------------------
  !! This routine computes the Fourier transform of the local
  !! part of an atomic pseudopotential, given in numerical form.
  !! A term erf(r)/r is subtracted in real space (thus making the
  !! function short-ranged) and added again in G space (for G<>0)
  !! The G=0 term contains \int (V_loc(r)+ Ze^2/r) 4pi r^2 dr.
  !! This is the "alpha" in the so-called "alpha Z" term of the energy.
  !! Atomic Ry units everywhere.
  !
  USE kinds
  USE constants,    ONLY : pi, fpi, e2, eps8
  USE esm,          ONLY : do_comp_esm, esm_bc
  USE Coul_cut_2D,  ONLY : do_cutoff_2D, lz
  !
  IMPLICIT NONE
  !
  integer, intent(in) :: ngl
  !! the number of shells of G vectors
  INTEGER, INTENT(IN) :: mesh
  !! number of grid points in the radial grid
  INTEGER, INTENT(IN) :: msh
  !! as above, used for radial integration
  REAL(DP), INTENT(IN) :: zp
  !! valence pseudocharge
  REAL(DP), INTENT(IN) :: rab(mesh)
  !! the derivative of mesh points
  REAL(DP), INTENT(IN) :: r(mesh)
  !! the mesh points
  REAL(DP), INTENT(IN) :: vloc_at(mesh)
  !! local part of the atomic pseudopotential on the radial mesh
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the moduli of g vectors for each shell 
  REAL(DP), INTENT(OUT) :: vloc(ngl)
  !! the fourier transform of the potential
  !
END SUBROUTINE vloc_of_g
!
!----------------------------------------------------------------------
SUBROUTINE vloc_coul( zp, tpiba2, ngl, gl, omega, vloc )
  !----------------------------------------------------------------------
  !! Fourier transform of the Coulomb potential - For all-electron
  !! calculations, in specific cases only, for testing purposes.
  !
  USE kinds
  USE constants,  ONLY : fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngl
  !! the number of shells of G vectors
  REAL(DP), INTENT(IN) :: zp
  !! valence pseudocharge
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the moduli of g vectors for each shell
  !
END SUBROUTINE vloc_coul

