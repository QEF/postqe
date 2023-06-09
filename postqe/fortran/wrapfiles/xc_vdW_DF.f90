! Copyright (C) 2001-2009 Quantum ESPRESSO group
! Copyright (C) 2019 Brian Kolb, Timo Thonhauser
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! ----------------------------------------------------------------------


MODULE vdW_DF

!! This module calculates the non-local correlation contribution to the
!! energy and potential according to:
!
!! M. Dion, H. Rydberg, E. Schroeder, D.C. Langreth, and
!! B.I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004).
!
!! henceforth referred to as DION. Further information about the
!! functional and its corresponding potential can be found in:
!
!! T. Thonhauser, V.R. Cooper, S. Li, A. Puzder, P. Hyldgaard,
!! and D.C. Langreth, Phys. Rev. B 76, 125112 (2007).
!
!! The proper spin extension of vdW-DF, i.e. svdW-DF, is derived in:
!
!! T. Thonhauser, S. Zuluaga, C.A. Arter, K. Berland, E. Schroder,
!! and P. Hyldgaard, Phys. Rev. Lett. 115, 136402 (2015).
!
!! henceforth referred to as THONHAUSER.
!
!! Two review articles show many of the vdW-DF applications:
!
!! D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009).
!
!! K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015).
!
!! The method implemented is based on the method of G. Roman-Perez and
!! J.M. Soler described in:
!
!! G. Roman-Perez and J.M. Soler, Phys. Rev. Lett. 103, 096102 (2009).
!
!! henceforth referred to as SOLER.
!
!
! xc_vdW_DF and xc_vdW_DF_spin are the driver routines for vdW-DF
! calculations and are called from Modules/funct.f90. The routines here
! set up the parallel run (if any) and carry out the calls necessary to
! calculate the non-local correlation contributions to the energy and
! potential.


USE kinds,             ONLY : dp
USE constants,         ONLY : pi, fpi, e2
USE mp,                ONLY : mp_sum, mp_barrier, mp_get, mp_size, mp_rank, mp_bcast
USE mp_images,         ONLY : intra_image_comm
USE mp_bands,          ONLY : intra_bgrp_comm
USE io_global,         ONLY : stdout, ionode
USE fft_base,          ONLY : dfftp
USE fft_interfaces,    ONLY : fwfft, invfft
USE control_flags,     ONLY : iverbosity, gamma_only

! ----------------------------------------------------------------------
! No implicit variables

IMPLICIT NONE


! ----------------------------------------------------------------------
! By default everything is private

PRIVATE


! ----------------------------------------------------------------------
! Save all objects in this module

SAVE


! ----------------------------------------------------------------------
! Public functions

PUBLIC  :: xc_vdW_DF, xc_vdW_DF_spin, vdW_DF_stress,                   &
           vdW_DF_energy, vdW_DF_potential,                            &
           generate_kernel, interpolate_kernel,                        &
           initialize_spline_interpolation, spline_interpolation, pw,pw_spin


! ----------------------------------------------------------------------
! Public variables

PUBLIC  :: inlc, vdW_DF_analysis, Nr_points, r_max, q_min, q_cut, Nqs, q_mesh


! ----------------------------------------------------------------------
! General variables

INTEGER                  :: inlc            = 1
! The non-local correlation

INTEGER                  :: vdW_DF_analysis = 0
! vdW-DF analysis tool as described in PRB 97, 085115 (2018)

REAL(DP), PARAMETER      :: epsr            = 1.0D-12
! A small number to cut off densities

INTEGER                  :: idx
! Indexing variable


! ----------------------------------------------------------------------
! Kernel specific parameters and variables

INTEGER, PARAMETER       :: Nr_points = 1024
! The number of radial points (also the number of k points) used in the
! formation of the kernel functions for each pair of q values.
! Increasing this value will help in case you get a run-time error
! saying that you are trying to use a k value that is larger than the
! largest tabulated k point since the largest k point will be 2*pi/r_max
! * Nr_points. Memory usage of the vdW_DF piece of PWSCF will increase
! roughly linearly with this variable.

REAL(DP), PARAMETER      :: r_max     = 100.0D0
! The value of the maximum radius to use for the real-space kernel
! functions for each pair of q values. The larger this value is the
! smaller the smallest k value will be since the smallest k point value
! is 2*pi/r_max. Be careful though, since this will also decrease the
! maximum k point value and the vdW_DF code will crash if it encounters
! a g-vector with a magnitude greater than 2*pi/r_max *Nr_points.

REAL(DP), PARAMETER      :: dr = r_max/Nr_points, dk = 2.0D0*pi/r_max
! Real space and k-space spacing of grid points.

REAL(DP), PARAMETER      :: q_min = 1.0D-5, q_cut = 5.0D0
! The maximum and minimum values of q. During a vdW run, values of q0
! found larger than q_cut will be saturated (SOLER equation 5) to q_cut.

INTEGER,  PARAMETER                 :: Nqs    = 20
REAL(DP), PARAMETER, DIMENSION(Nqs) :: q_mesh = (/  &
   q_min              , 0.0449420825586261D0, 0.0975593700991365D0, 0.159162633466142D0, &
   0.231286496836006D0, 0.315727667369529D0 , 0.414589693721418D0 , 0.530335368404141D0, &
   0.665848079422965D0, 0.824503639537924D0 , 1.010254382520950D0 , 1.227727621364570D0, &
   1.482340921174910D0, 1.780437058359530D0 , 2.129442028133640D0 , 2.538050036534580D0, &
   3.016440085356680D0, 3.576529545442460D0 , 4.232271035198720D0 , q_cut /)

! The above two parameters define the q mesh to be used in the vdW_DF
! code. These are perhaps the most important to have set correctly.
! Increasing the number of q points will DRAMATICALLY increase the
! memory usage of the vdW_DF code because the memory consumption depends
! quadratically on the number of q points in the mesh. Increasing the
! number of q points may increase accuracy of the vdW_DF code, although,
! in testing it was found to have little effect. The largest value of
! the q mesh is q_cut. All values of q0 (DION equation 11) larger than
! this value during a run will be saturated to this value using equation
! 5 of SOLER. In testing, increasing the value of q_cut was found to
! have little impact on the results, although it is possible that in
! some systems it may be more important. Always make sure that the
! variable Nqs is consistent with the number of q points that are
! actually in the variable q_mesh. Also, do not set any q value to 0.
! This will cause an infinity in the Fourier transform.

INTEGER,  PARAMETER      :: Nintegration_points = 256
! Number of integration points for real-space kernel generation (see
! DION equation 14). This is how many a's and b's there will be.

REAL(DP), PARAMETER      :: a_min = 0.0D0, a_max = 64.0D0
! Min/max values for the a and b integration in DION equation 14.

REAL(DP) :: kernel( 0:Nr_points, Nqs, Nqs ), d2phi_dk2( 0:Nr_points, Nqs, Nqs )
! Matrices holding the Fourier transformed kernel function  and its
! second derivative for each pair of q values. The ordering is
! kernel(k_point, q1_value, q2_value).

REAL(DP) :: W_ab( Nintegration_points, Nintegration_points )
! Defined in DION equation 16.

REAL(DP) :: a_points( Nintegration_points ), a_points2( Nintegration_points )
! The values of the "a" points (DION equation 14) and their squares.

CONTAINS








  ! ####################################################################
  !                           |             |
  !                           |  functions  |
  !                           |_____________|
  !
  ! Functions to be used in get_q0_on_grid, get_q0_on_grid_spin, and
  ! phi_value.

  FUNCTION Fs(s)
     IMPLICIT NONE
     REAL(DP) :: s, Fs, Z_ab = 0.0D0
  END FUNCTION Fs




  FUNCTION dFs_ds(s)
     IMPLICIT NONE
     REAL(DP)             :: s, dFs_ds, Z_ab = 0.0D0
  END FUNCTION dFs_ds




  FUNCTION kF(rho)
     IMPLICIT NONE
     REAL(DP)             :: rho, kF
  END FUNCTION kF




  FUNCTION dkF_drho(rho)
     IMPLICIT NONE
     REAL(DP)             :: rho, dkF_drho
  END FUNCTION dkF_drho




  FUNCTION ds_drho(rho, s)
     IMPLICIT NONE
     REAL(DP) :: rho, s, ds_drho
  END FUNCTION ds_drho




  FUNCTION ds_dgradrho(rho)

     IMPLICIT NONE
     REAL(DP) :: rho, ds_dgradrho
  END FUNCTION ds_dgradrho




  FUNCTION dqx_drho(rho, s)
     IMPLICIT NONE
     REAL(DP) :: rho, s, dqx_drho
  END FUNCTION dqx_drho




  FUNCTION h_function(y)

     IMPLICIT NONE
     REAL(DP)             :: y, h_function
  END FUNCTION








  ! ####################################################################
  !                           |             |
  !                           |  XC_VDW_DF  |
  !                           |_____________|

  SUBROUTINE xc_vdW_DF (rho_valence, rho_core, etxc, vtxc, v)
  !
  !! Driver routine for vdW-DF calculations, called from Modules/funct.f90.
  !! The routine here sets up the parallel run (if any) and carry out the 
  !! calls necessary to calculate the non-local correlation contributions
  !! to the energy and potential.
  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba
  IMPLICIT NONE

  ! --------------------------------------------------------------------
  ! Arguments 
  !                                               _
  REAL(DP), INTENT(IN)    :: rho_valence(:,:)    !
  REAL(DP), INTENT(IN)    :: rho_core(:)         !  PWSCF input variables
  REAL(DP), INTENT(INOUT) :: etxc, vtxc, v(:,:)  !_
  END SUBROUTINE xc_vdW_DF








  ! ####################################################################
  !                          |                  |
  !                          |  XC_VDW_DF_spin  |
  !                          |__________________|
  !
  !
  SUBROUTINE xc_vdW_DF_spin (rho_valence, rho_core, etxc, vtxc, v)
  
  !! This subroutine is as similar to \(\texttt{xc_vdW_DF}\) as possible,
  !! but handles the collinear \(\text{nspin}=2\) case.
  
  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba
  IMPLICIT NONE
  ! --------------------------------------------------------------------
  ! Arguments
  !                                              _
  REAL(DP), INTENT(IN) :: rho_valence(:,:)      !
  REAL(DP), INTENT(IN) :: rho_core(:)           ! PWSCF input variables.
  REAL(DP), INTENT(INOUT) :: etxc, vtxc, v(:,:) !_
  END SUBROUTINE xc_vdW_DF_spin








  ! ####################################################################
  !                       |                  |
  !                       |  GET_Q0_ON_GRID  |
  !                       |__________________|
  !
  !
  SUBROUTINE get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)
  !! This routine first calculates the q value defined in (DION equations
  !! 11 and 12), then saturates it according to (SOLER equation 5). More
  !! specifically it calculates the following:
  !
  !! * \(\text{q0(ir)}\) = saturated value of q;
  !! * \(\text{dq0_drho(ir)} = \text{total_rho} * d\text{q0} /d\text{rho}\)
  !! * \(\text{dq0_dgradrho} = \text{total_rho} / |\text{grad_rho}| * d\text{q0}
  !!   / d |\text{grad_rho}|\)
  
  IMPLICIT NONE
  REAL(DP), INTENT(IN)      :: total_rho(:), grad_rho(:,:)         ! Input variables needed.
  REAL(DP), INTENT(OUT)     :: q0(:), dq0_drho(:), dq0_dgradrho(:) ! Output variables that have been allocated
  COMPLEX(DP), INTENT(INOUT):: thetas(:,:)                         ! The thetas from SOLER.
  END SUBROUTINE get_q0_on_grid


  !-----------------------------------------------------------------------
  SUBROUTINE pw( rs, iflag, ec, vc )
    !-----------------------------------------------------------------------
    ! --A provisional copy of the pw routine in XC lib to avoid external calls--
    !! * iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
    !! * iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
    !
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: rs
    INTEGER, INTENT(IN)  :: iflag
    REAL(DP), INTENT(OUT) :: ec, vc
  END SUBROUTINE pw
  !---------------------

  ! ####################################################################
  !                       |                       |
  !                       |  GET_Q0_ON_GRID_spin  |
  !                       |_______________________|

  SUBROUTINE get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, &
             grad_rho_up, grad_rho_down, q0, dq0_drho_up, dq0_drho_down, &
             dq0_dgradrho_up, dq0_dgradrho_down, thetas)
  
  !! Find the value of \(\text{q0}\) for all assigned grid points. q is defined in
  !! equations 11 and 12 of DION and q0 is the saturated version of q
  !! defined in equation 5 of SOLER. In the spin case, q0 is defined by
  !! equation 8 (and text above that equation) of THONHAUSER. This
  !! routine also returns the derivatives of the \(\text{q0}\)s with respect to the
  !! charge-density and the gradient of the charge-density. These are
  !! needed for the potential.
  
  IMPLICIT NONE

  REAL(DP),  INTENT(IN)      :: total_rho(:), grad_rho(:,:)              ! Input variables.
  REAL(DP),  INTENT(IN)      :: rho_up(:), grad_rho_up(:,:)              ! Input variables.
  REAL(DP),  INTENT(IN)      :: rho_down(:), grad_rho_down(:,:)          ! Input variables.

  REAL(DP),  INTENT(OUT)     :: q0(:), dq0_drho_up(:), dq0_drho_down(:)  ! Output variables.
  REAL(DP),  INTENT(OUT)     :: dq0_dgradrho_up(:), dq0_dgradrho_down(:) ! Output variables.
  COMPLEX(DP), INTENT(INOUT) :: thetas(:,:)                              ! The thetas from SOLER.
  END SUBROUTINE get_q0_on_grid_spin

  
  !-----------------------------------------------------------------------
  SUBROUTINE pw_spin( rs, zeta, ec, vc_up, vc_dw )
    !-----------------------------------------------------------------------
    !--A provisional copy of the pw routine in XC lib to avoid external calls--
    !! J.P. Perdew and Y. Wang, PRB 45, 13244 (1992).
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: rs
    !! Wigner-Seitz radius
    REAL(DP), INTENT(IN) :: zeta
    !! zeta = (rho_up - rho_dw)/rho_tot
    REAL(DP), INTENT(OUT) :: ec, vc_up, vc_dw
    !
  END SUBROUTINE pw_spin

  

  ! ####################################################################
  !                            |              |
  !                            |  saturate_q  |
  !                            |______________|

  SUBROUTINE saturate_q (q, q_cutoff, q0, dq0_dq)
  
  !! Here, we calculate \(\text{q0}\) by saturating \(\text{q}\) according to
  !! equation 5 of SOLER. Also, we find the derivative \(\text{dq0/dq}\) needed
  !! for the derivatives \(\text{dq0/drho}\) and \(\text{dq0/dgradrh0}\).
  
  IMPLICIT NONE

  REAL(DP),  INTENT(IN)      :: q             ! Input q.
  REAL(DP),  INTENT(IN)      :: q_cutoff      ! Cutoff q.
  REAL(DP),  INTENT(OUT)     :: q0            ! Output saturated q.
  REAL(DP),  INTENT(OUT)     :: dq0_dq        ! Derivative of dq0/dq.
  END SUBROUTINE saturate_q








  ! ####################################################################
  !                          |               |
  !                          | vdW_DF_energy |
  !                          |_______________|
  !
  !
  SUBROUTINE vdW_DF_energy (thetas, vdW_xc_energy)

  !! This routine carries out the integration of equation 8 of SOLER.  It
  !! returns the non-local exchange-correlation energy and the \(\text{u_alpha(k)}\)
  !! arrays used to find the \(\text{u_alpha(r)}\) arrays via equations 11 and 12 in
  !! SOLER.
  
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  IMPLICIT NONE

  COMPLEX(DP), INTENT(INOUT) :: thetas(:,:)            ! On input this variable holds the theta
                                                       ! functions (equation 8, SOLER) in the format
                                                       ! thetas(grid_point, theta_i). On output
                                                       ! this array holds u_alpha(k) =
                                                       ! Sum_j[theta_beta(k)phi_alpha_beta(k)].

  REAL(DP), INTENT(OUT) :: vdW_xc_energy               ! The non-local correlation energy.
  END SUBROUTINE vdW_DF_energy








  ! ####################################################################
  !                        |                   |
  !                        |  vdW_DF_potential |
  !                        |___________________|
  !
  !
  SUBROUTINE vdW_DF_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, u_vdW, potential)

  !! This routine finds the non-local correlation contribution to the potential (i.e.
  !! the derivative of the non-local piece of the energy with respect to density)
  !! given in SOLER equation 10. The \(\text{u_alpha(k)}\) functions were found 
  !! while calculating the energy. They are passed in as the matrix \(\text{u_vdW}\).
  !! Most of the required derivatives were calculated in the \(\texttt{get_q0_on_grid}\)
  !! routine, but the derivative of the interpolation polynomials, \(P_\alpha(q)\), 
  !! (SOLER equation 3) with respect to q is interpolated here, along with the polynomials
  !! themselves.
  
  USE gvect,               ONLY : g
  USE cell_base,           ONLY : alat, tpiba

  IMPLICIT NONE

  REAL(DP), INTENT(IN) ::  q0(:), grad_rho(:,:)       ! Input arrays holding the value of q0 for
                                                      ! all points assigned to this processor and
                                                      ! the gradient of the charge density for
                                                      ! points assigned to this processor.

  REAL(DP), INTENT(IN) :: dq0_drho(:), dq0_dgradrho(:)! The derivative of q0 with respect to the
                                                      ! charge density and gradient of the charge
                                                      ! density (almost). See comments in the
                                                      ! get_q0_on_grid subroutine above.

  COMPLEX(DP), INTENT(IN) :: u_vdW(:,:)               ! The functions u_alpha(r) obtained by
                                                      ! inverse transforming the functions
                                                      ! u_alph(k). See equations 11 and 12 in SOLER.

  REAL(DP), INTENT(INOUT) :: potential(:)             ! The non-local correlation potential for
                                                      ! points on the grid over the whole cell (not
                                                      ! just those assigned to this processor).

  END SUBROUTINE vdW_DF_potential








  ! ####################################################################
  !                       |                        |
  !                       |  SPLINE_INTERPOLATION  |
  !                       |________________________|
  !
  !
  SUBROUTINE spline_interpolation (x, evaluation_points, values)
  !! This routine is modeled after an algorithm from NUMERICAL_RECIPES It
  !! was adapted for Fortran, of course and for the problem at hand, in
  !! that it finds the bin a particular x value is in and then loops over
  !! all the \(P_i\) functions so we only have to find the bin once.
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x(:), evaluation_points(:)     ! Input variables. The x values used to
                                                         ! form the interpolation (q_mesh in this
                                                         ! case) and the values of q0 for which we
                                                         ! are interpolating the function.

  COMPLEX(DP), INTENT(INOUT) :: values(:,:)              ! An output array (allocated outside this
                                                         ! routine) that stores the interpolated
                                                         ! values of the P_i (SOLER equation 3)
                                                         ! polynomials. The format is
                                                         ! values(grid_point, P_i).

  END SUBROUTINE spline_interpolation








  ! ####################################################################
  !                  |                                   |
  !                  |  INITIALIZE_SPLINE_INTERPOLATION  |
  !                  |___________________________________|
  !
  
  SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)
  
  !! This routine is modeled after an algorithm from NUMERICAL RECIPES It
  !! was adapted for Fortran and for the problem at hand.

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: x(:)               ! The input abscissa values.
  REAL(DP), INTENT(INOUT) :: d2y_dx2(:,:)       ! The output array (allocated outside this routine)
                                                ! that holds the second derivatives required for

  END SUBROUTINE initialize_spline_interpolation








  ! ####################################################################
  !                          |                    |
  !                          | INTERPOLATE_KERNEL |
  !                          |____________________|
  !
  
  SUBROUTINE interpolate_kernel (k, kernel_of_k)

  !! This routine is modeled after an algorithm from NUMERICAL RECIPES
  !! Adapted for Fortran and the problem at hand. This function is used
  !! to find the Phi-alpha-beta needed for equations 8 and 11 of SOLER.
  
  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: k                 ! Input value, the magnitude of the g-vector
                                               ! for the current point.

  REAL(DP), INTENT(INOUT) :: kernel_of_k(:,:)  ! An output array (allocated outside this routine)
                                               ! that holds the interpolated value of the kernel
                                               ! for each pair of q points (i.e. the phi_alpha_beta
                                               ! of the Soler method.

  END SUBROUTINE interpolate_kernel








  ! ####################################################################
  !                         |                 |
  !                         |  VDW_DF_STRESS  |
  !                         |_________________|

  SUBROUTINE vdW_DF_stress (rho_valence, rho_core, nspin, sigma) ! PH adjusted for wrapper spin/nospin
  
  !! vdW-DF stress calculation.
  
  use gvect,           ONLY : ngm, g
  USE cell_base,       ONLY : tpiba

  IMPLICIT NONE

  REAL(DP), INTENT(IN)     :: rho_valence(:,:)       !
  REAL(dp), INTENT(IN)     :: rho_core(:)            ! Input variables.
  INTEGER,  INTENT(IN)     :: nspin                  !
  REAL(dp), INTENT(INOUT)  :: sigma(3,3)             !
  END SUBROUTINE vdW_DF_stress ! PH adjusted for wrapper spin/nospin


! -------------------------------------------------------------------------
 ! Begin Spin vdW-DF_strees_gradient implemented Per Hyldgaard 2019, GPL. No Waranties
 ! Adapted from the original nspin = 1 code (subroutine below) by Thonhauser and coauthors
! -------------------------------------------------------------------------

  ! ####################################################################
  !                     |                               |
  !                     |  VDW_DF_STRESS_GRADIENT_SPIN  |
  !                     |_______________________________|

  SUBROUTINE vdW_DF_stress_gradient_spin (total_rho, grad_rho_up, grad_rho_down, q0, &
                                          dq0_dgradrho_up, dq0_dgradrho_down, &
                                          thetas, sigma)

  USE gvect,                 ONLY : ngm, g, gg, gstart
  USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2

  implicit none

  real(dp), intent(IN)     :: total_rho(:)           !
  real(dp), intent(IN)     :: grad_rho_up (:, :)     ! Input variables.
  real(dp), intent(IN)     :: grad_rho_down(:, :)    !
  real(dp), intent(inout)  :: sigma(:,:)             !
  real(dp), intent(IN)     :: q0(:)                  !
  real(dp), intent(IN)     :: dq0_dgradrho_up(:)     !
  real(dp), intent(IN)     :: dq0_dgradrho_down(:)   !
  complex(dp), intent(IN)  :: thetas(:,:)            !
  END SUBROUTINE vdW_DF_stress_gradient_spin


  ! ####################################################################
  !                     |                          |
  !                     |  VDW_DF_STRESS_GRADIENT  |
  !                     |__________________________|

  SUBROUTINE vdW_DF_stress_gradient (total_rho, grad_rho, q0, &
             dq0_drho, dq0_dgradrho, thetas, sigma)

  USE gvect,                 ONLY : ngm, g, gstart
  USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2

  IMPLICIT NONE

  REAL(DP), INTENT(IN)     :: total_rho(:)           !
  REAL(DP), INTENT(IN)     :: grad_rho(:, :)         !
  REAL(DP), INTENT(INOUT)  :: sigma(:,:)             !
  REAL(DP), INTENT(IN)     :: q0(:)                  ! Input variables.
  REAL(DP), INTENT(IN)     :: dq0_drho(:)            !
  REAL(DP), INTENT(IN)     :: dq0_dgradrho(:)        !
  COMPLEX(DP), INTENT(IN)  :: thetas(:,:)            !
  END SUBROUTINE vdW_DF_stress_gradient








  ! ####################################################################
  !                      |                        |
  !                      |  VDW_DF_STRESS_KERNEL  |
  !                      |________________________|

  SUBROUTINE vdW_DF_stress_kernel (total_rho, q0, thetas, sigma)

  USE gvect,                 ONLY : ngm, g, gg, igtongl, gl, ngl, gstart
  USE cell_base,             ONLY : omega, tpiba, tpiba2

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: q0(:)
  REAL(DP), INTENT(IN)    :: total_rho(:)
  REAL(DP), INTENT(INOUT) :: sigma(3,3)
  COMPLEX(DP), INTENT(IN) :: thetas(:,:)
  END SUBROUTINE vdW_DF_stress_kernel








  ! ####################################################################
  !                        |                        |
  !                        | INTERPOLATE_DKERNEL_DK |
  !                        |________________________|

  SUBROUTINE interpolate_Dkernel_Dk (k, dkernel_of_dk)

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: k                      ! Input value, the magnitude of the g-vector
                                                    ! for the current point.

  REAL(DP), INTENT(INOUT) :: dkernel_of_dk(Nqs,Nqs) ! An output array (allocated outside this
                                                    ! routine) that holds the interpolated value of
                                                    ! the kernel for each pair of q points (i.e. the
                                                    ! phi_alpha_beta of the Soler method.

  END SUBROUTINE interpolate_Dkernel_Dk








  ! ####################################################################
  !                          |              |
  !                          | thetas_to_uk |
  !                          |______________|


  SUBROUTINE thetas_to_uk (thetas, u_vdW)

  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)  :: thetas(:,:)       ! On input this variable holds the theta functions
                                                ! (equation 8, SOLER) in the format
                                                ! thetas(grid_point, theta_i).
  COMPLEX(DP), INTENT(OUT) :: u_vdW(:,:)        ! On output this array holds u_alpha(k) =
                                                ! Sum_j[theta_beta(k)phi_alpha_beta(k)].
  END SUBROUTINE thetas_to_uk








  ! ####################################################################
  !                           |                 |
  !                           | GENERATE_KERNEL |
  !                           |_________________|
  !
  ! The original definition of the kernel function is given in DION
  ! equations 14-16. The Soler method makes the kernel function a
  ! function of only 1 variable (r) by first putting it in the form
  ! phi(q1*r, q2*r). Then, the q-dependence is removed by expanding the
  ! function in a special way (see SOLER equation 3). This yields a
  ! separate function for each pair of q points that is a function of r
  ! alone. There are (Nqs^2+Nqs)/2 unique functions, where Nqs is the
  ! number of q points used. In the Soler method, the kernel is first
  ! made in the form phi(d1, d2) but this is not done here. It was found
  ! that, with q's chosen judiciously ahead of time, the kernel and the
  ! second derivatives required for interpolation could be tabulated
  ! ahead of time for faster use of the vdW-DF functional. Through
  ! testing we found no need to soften the kernel and correct for this
  ! later (see SOLER eqations 6-7).
  !
  ! The algorithm employed here is "embarrassingly parallel," meaning
  ! that it parallelizes very well up to (Nqs^2+Nqs)/2 processors,
  ! where, again, Nqs is the number of q points chosen. However,
  ! parallelization on this scale is unnecessary. In testing the code
  ! runs in under a minute on 16 Intel Xeon processors.
  !
  ! IMPORTANT NOTICE: Results are very sensitive to compilation details.
  ! In particular, the usage of FMA (Fused Multiply-and-Add)
  ! instructions used by modern CPUs such as AMD Interlagos (Bulldozer)
  ! and Intel Ivy Bridge may affect quite heavily some components of the
  ! kernel (communication by Ake Sandberg, Umea University). In practice
  ! this should not be a problem, since most affected elements are the
  ! less relevant ones.
  !
  ! Some of the algorithms here are somewhat modified versions of those
  ! found in:
  !
  !    Numerical Recipes in C; William H. Press, Brian P. Flannery, Saul
  !    A. Teukolsky, and William T. Vetterling. Cambridge University
  !    Press (1988).
  !
  ! hereafter referred to as NUMERICAL_RECIPES. The routines were
  ! translated to Fortran, of course and variable names are generally
  ! different.
  !
  ! For the calculation of the kernel we have benefited from access to
  ! earlier vdW-DF implementation into PWscf and ABINIT, written by Timo
  ! Thonhauser, Valentino Cooper, and David Langreth. These codes, in
  ! turn, benefited from earlier codes written by Maxime Dion and Henrik
  ! Rydberg.

  SUBROUTINE generate_kernel
  IMPLICIT NONE
  ! Indexing variables.
  END SUBROUTINE generate_kernel








  ! ####################################################################
  !                    |                            |
  !                    |  PREP_GAUSSIAN_QUADRATURE  |
  !                    |____________________________|
  !

  SUBROUTINE prep_gaussian_quadrature( weights )

  !! Routine to calculate the points and weights for the
  !! Gaussian-Legendre integration. This routine is modeled after the
  !! routine GAULEG from NUMERICAL RECIPES.
  
  REAL(DP), INTENT(INOUT) :: weights(:)
  ! The points and weights for the Gaussian-Legendre integration.
  END SUBROUTINE prep_gaussian_quadrature








  ! ####################################################################
  !                            |             |
  !                            |  PHI_VALUE  |
  !                            |_____________|
  !
  
  REAL(DP) FUNCTION phi_value(d1, d2)

  !! This function returns the value of the kernel calculated via DION
  !! equation 14.
  
  REAL(DP), INTENT(IN) :: d1, d2
  !! The point at which to evaluate the kernel. d1 = q1*r and d2 = q2*r.
  END FUNCTION phi_value








  ! ####################################################################
  !                            |              |
  !                            |  RADIAL_FFT  |
  !                            |______________|
  !
  
  SUBROUTINE radial_fft(phi)
  
  !! This subroutine performs a radial Fourier transform on the
  !! real-space kernel functions.
  ! Basically, this is just:
  !            int(4*pi*r^2*phi*sin(k*r)/(k*r))dr
  ! integrated from 0 to r_max.
  ! That is, it is the kernel function phi integrated with the 0^th spherical
  ! Bessel function radially, with a 4*pi assumed from angular
  ! integration since we have spherical symmetry. The spherical symmetry
  ! comes in because the kernel function depends only on the magnitude
  ! of the vector between two points. The integration is done using the
  ! trapezoid rule.
  
  REAL(DP), INTENT(INOUT) :: phi(0:Nr_points)
  !! On input holds the real-space function phi_q1_q2(r).
  !! On output hold the reciprocal-space function phi_q1_q2(k).
  END SUBROUTINE radial_fft








  ! ####################################################################
  !                          |                  |
  !                          |  SET UP SPLINES  |
  !                          |__________________|
  !
  
  SUBROUTINE set_up_splines(phi, D2)
  
  !! This subroutine accepts a function (phi) and finds at each point the
  !! second derivative (D2) for use with spline interpolation. This
  !! function assumes we are using the expansion described in SOLER
  !! equation 3.
  ! That is, the derivatives are those needed to interpolate
  ! Kronecker delta functions at each of the q values. Other than some
  ! special modification to speed up the algorithm in our particular
  ! case, this algorithm is taken directly from NUMERICAL_RECIPES.

  REAL(DP), INTENT(IN)    :: phi(0:Nr_points)
  !! The k-space kernel function for a particular q1 and q2.

  REAL(DP), INTENT(INOUT) :: D2(0:Nr_points)
  !! The second derivatives to be used in the interpolation expansion
  !! (SOLER equation 3).
  END SUBROUTINE set_up_splines
  ! ####################################################################
  !                          |            |
  !                          |  VDW_INFO  |
  !                          |____________|

  SUBROUTINE vdW_info( nspin )

  USE xc_lib, ONLY : xclib_dft_is

  IMPLICIT NONE
  INTEGER, INTENT (IN) :: nspin
  END SUBROUTINE
END MODULE vdW_DF
