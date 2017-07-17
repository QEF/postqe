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
  ! mesh: number of grid points in the radial grid
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
