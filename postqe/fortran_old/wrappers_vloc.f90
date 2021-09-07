!
! Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
! Internazionale Superiore di Studi Avanzati). All rights reserved.
! This file is distributed under the terms of the LGPL-2.1 license. See the
! file 'LICENSE' in the root directory of the present distribution, or
! https://opensource.org/licenses/LGPL-2.1
!
! Wrappers for vloc QE subroutines 
!


subroutine pyqe_vloc_of_g(mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, gl, omega, vloc)
  !----------------------------------------------------------------------
  !    
  !    A wrapper for the vloc_of_g from QE
  !
  !    This routine computes the Fourier transform of the local
  !    part of an atomic pseudopotential, given in numerical form.
  !    A term erf(r)/r is subtracted in real space (thus making the
  !    function short-ramged) and added again in G space (for G<>0)
  !    The G=0 term contains \int (V_loc(r)+ Ze^2/r) 4pi r^2 dr.
  !    This is the "alpha" in the so-called "alpha Z" term of the energy.
  !    Atomic Ry units everywhere.
  !
  implicit none
  integer, intent(in) :: ngl, mesh, msh
  ! ngl : the number of shells of G vectors
  ! mesh: number of grid points in the radial gridlibQ
  ! msh : as above, used for radial integration
  !
  real(8), intent(in) :: zp, rab (mesh), r (mesh), vloc_at (mesh), tpiba2, &
                          omega, gl (ngl)
  ! zp : valence pseudocharge
  ! rab: the derivative of mesh points
  ! r  : the mesh points
  ! vloc_at: local part of the atomic pseudopotential on the radial mesh
  ! tpiba2 : 2 pi / alat
  ! omega  : the volume of the unit cell
  ! gl     : the moduli of g vectors for each shell
  !
  real(8), intent(out):: vloc (ngl)
  !
  call vloc_of_g (mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, &
     gl, omega, vloc) 
end subroutine pyqe_vloc_of_g

subroutine pyqe_struct_fact(nat, tau, ngm, g, strf ,check_gg, check_tau)
  implicit none
  !
  !   calculate the structure factors for each type of atoms in the unit
  !   cell
  !
  !   Here the dummy variables
  !
  integer, intent(in) :: nat, ngm 
  ! input: the number of identical atom in the unit cell
  ! input: the number of G vectors
  ! input: fft dimension along x
  ! input: fft dimension along y
  ! input: fft dimension along z
  
  real(8),intent (in)  :: tau (nat, 3 ), g (ngm, 3 )
  ! input: reciprocal crystal basis vectors
  ! input: the positions of the atoms in the c
  ! input: the coordinates of the g vectors

  complex(8),intent(out)  :: strf (ngm ) 
  real(8), intent (out)  :: check_gg(ngm) 
  real(8), intent (out)  :: check_tau(3*nat) 
  ! output: the structure factor
  call struc_fact2(nat, tau, ngm, g, strf ,check_gg, check_tau)
end subroutine pyqe_struct_fact
