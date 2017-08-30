!
! Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
! Internazionale Superiore di Studi Avanzati). All rights reserved.
! This file is distributed under the terms of the LGPL-2.1 license. See the
! file 'LICENSE' in the root directory of the present distribution, or
! https://opensource.org/licenses/LGPL-2.1
!
! A collection of utility subroutines better written in Fortran than in Python and wrappers for QE utility routines 
!

subroutine pyqe_getcelldms(alat, at1, at2, at3, ibrav, celldm) 
  ! Returns the celldms parameters from ibrav, alat and the lattice vectors
  implicit  none
  real(8), intent (in)  :: at1(3), at2(3), at3(3), alat
  integer , intent (in)  :: ibrav
  real(8), intent (out) :: celldm(6)
  !
  real(8)     :: at(3,3) 
  !
    celldm = 0.d0
    at(:,1) = at1
    at(:,2) = at2
    at(:,3) = at3
    SELECT CASE  (ibrav ) 
       CASE (0:3,-3) 
          celldm(1) = alat
          celldm(2:6) = 0.d0
       CASE (4) 
          celldm(1) = alat
          celldm(2) = 0.d0
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3)))/alat
          celldm(4:6) = 0.d0
       CASE (5) 
          celldm(1)= alat
          celldm(2:3) = 0.d0
          celldm(4) = DOT_PRODUCT(at(:,1),at(:,2))/(alat**2)
          celldm(5:6) = 0.d0
       CASE (6) 
          celldm(1)= alat 
          celldm(3)= SQRT( DOT_PRODUCT(at(:,3),at(:,3)))/alat
          celldm(2)= 1.d0
          celldm(4:6) = 0.d0
       CASE (7) 
          celldm(1) = alat
          celldm(3) = ABS(at(3,3)/at(1,3)) 
          celldm(2)=0.d0
          celldm(4:6) = 0.d0
       CASE (8)
          celldm(1) = alat
          celldm(2) = SQRT( DOT_PRODUCT (at(:,2),at(:,2)))/alat
          celldm(3) = SQRT( DOT_PRODUCT (at(:,3),at(:,3)))/alat 
          celldm(4:6) = 0.d0
       CASE (9) 
          celldm(1) = alat
          celldm(2) = ABS ( at(2,1)/at(1,1))
          celldm(3) = ABS ( at(3,3)/2.d0/at(1,1))
          celldm(4:6) = 0.d0 
       CASE (10) 
          celldm(1) = alat
          celldm(2) = ABS ( at(2,2)/at(2,1))
          celldm(3) = ABS ( at(3,1)/at(1,1))
          celldm(4:6) = 0.d0
       CASE (11) 
          celldm(1) = alat
          celldm(2) = ABS(at(2,1)/at(1,1))
          celldm(3) = ABS(at(3,1)/at(1,1))
          celldm(4:6) = 0.d0
       CASE (12) 
          celldm(1) = alat 
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(4) = DOT_PRODUCT(at(:,1),at(:,2))/&
                      SQRT(DOT_PRODUCT(at(:,1),at(:,1))*DOT_PRODUCT(at(:,2),at(:,2)))
          celldm(5) =  DOT_PRODUCT(at(:,1),at(:,3))/&
                   SQRT(DOT_PRODUCT(at(:,1),at(:,1))*DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(6) = 0.d0
       CASE (13) 
          celldm(1) = alat
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2)))/(2.d0*at(1,1))
          celldm(3) = ABS (at(3,3)/at(1,3))
          celldm(4) = COS( ATAN2(at(2,2),at(1,2)) )
          celldm(5:6) = 0.d0
       CASE (14) 
          celldm(1) = alat 
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(4) = DOT_PRODUCT(at(:,3),at(:,2))/SQRT(DOT_PRODUCT(at(:,2),at(:,2))*&
                                                   DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(5) = DOT_PRODUCT(at(:,3),at(:,1))/SQRT(DOT_PRODUCT(at(:,1),at(:,1))*&
                                                   DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(6) = DOT_PRODUCT(at(:,1),at(:,2))/SQRT(DOT_PRODUCT(at(:,2),at(:,2))*&
                                                   DOT_PRODUCT(at(:,1),at(:,1)))
       CASE  default  
          celldm(1) = 1.d0
          IF (alat .GT. 0.d0 ) celldm(1) = alat
          celldm (2:6) = 0.d0
    END SELECT 
end subroutine pyqe_getcelldms


subroutine pyqe_recips(a1, a2, a3, b1, b2, b3)
  ! Returns the vectors of the reciprocal cell given those of the direct one
  implicit none 
  real(8), intent(in)  :: a1(3),a2(3),a3(3)
  real(8), intent(out) :: b1(3),b2(3),b3(3) 
  call recips(a1,a2,a3,b1,b2,b3) 
end subroutine pyqe_recips 


subroutine pyqe_latgen(ibrav, celldm, at)
  ! Returns the lattice vectors given ibrav and the celldms parameters
  implicit none
  integer, intent(in)   :: ibrav
  real(8),intent(in)    :: celldm(6) 
  real(8),intent(out)   :: at(3,3) 
  ! 
  real(8)               :: volume
  call latgen(ibrav,celldm, at(:,1), at(:,2), at(:,3), volume) 
end subroutine pyqe_latgen


subroutine pyqe_get_igtongl(ngm, gg, igtongl, ngl)
  implicit none
  real(8), parameter :: eps8 = 1.d-8

  integer , intent (in) :: ngm
  real(8), intent (in) ::   gg(ngm)
  integer, intent (out)  :: igtongl(ngm)
  integer, intent (out) :: ngl

  integer  :: ng, igl
  ngl = 1 
  igtongl(1) = 1 
  do ng = 2, ngm
    if (gg(ng) > gg(ng-1)+eps8) ngl=ngl+1
    igtongl(ng) = ngl
  end do 
end subroutine pyqe_get_igtongl


subroutine pyqe_get_gl(ngm, gg, ngl, gl)  
  implicit none 
  real(8), parameter :: eps8 = 1.d-8
  integer, intent(in)  :: ngm, ngl
  real(8), intent(in) :: gg(ngm)
  real(8), intent(out) :: gl(ngl)
  integer igl, ng
  igl = 1 
  gl(1) = gg(1)
  do ng = 1, ngm 
    if ( gg(ng) > gg( ng-1)+eps8) igl = igl+1 
    if (igl > ngl) stop 'igl > ngl ?' 
    gl(igl) = gg( ng)
  end do  
  !
end subroutine pyqe_get_gl


subroutine pyqe_get_gg_list( nrrr, nr1, nr2, nr3 , bg1, bg2, bg3, g, gg , mill) 
  implicit none
  real(8), parameter :: pi = 4.d0*atan(1.d0) , fpi = 4.d0*pi, eps8 = 1.d-8

  integer, intent(in)   :: nr1, nr2, nr3, nrrr 
  real(8),intent(in)   :: bg1(3), bg2(3), bg3(3)
  real(8), intent(out) :: gg(nrrr) , g(nrrr,3) 
  integer (8), intent(out) :: mill(nrrr,3) 
  real(8)              :: vec(3)
  integer ni, nj, nk, i, ii, j, ij, k, ik,ig
  ni = ( nr1 -1 ) /2  
  nj = ( nr2 -1 ) /2 
  nk = ( nr3 -1 ) /2  
  ig = 0
  do ii = 0, nr1 -1 
    if ( ii .le. nr1/2) then 
        i =  ii 
    else 
        i =  ii - nr1 
    end if
    do ij = 0, nr2 -1
        if ( ij .le. nr2/2 ) then 
          j = ij 
        else
          j = ij - nr2 
        end if  
        do ik = 0, nr3 - 1  
          if ( ik .le. nr3/2 ) then 
              k = ik 
          else 
              k = ik -nr3 
          end if
          ig = ig +1 
          vec    =  i*bg1+j*bg2+k*bg3
          g(ig,:)  = vec(:)
          gg(ig) = sum(vec*vec)
          mill(ig,:) = [ii,ij,ik]
        end do
    end do
  end do  
  !
end subroutine pyqe_get_gg_list

!
!-----------------------------------------------------------------------
subroutine py_w0gauss (x, n, w0gauss)
  !-----------------------------------------------------------------------
  !
  !     the derivative of wgauss:  an approximation to the delta function
  !
  ! --> (n>=0) : derivative of the corresponding Methfessel-Paxton wgauss
  !
  ! --> (n=-1 ): derivative of cold smearing:
  !              1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
  !
  ! --> (n=-99): derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
  !
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: sqrtpi = 1.77245385090551602729_DP 
  REAL(DP), PARAMETER :: sqrtpm1= 1.0_DP / sqrtpi
  real(DP), intent(out) :: w0gauss
  real(DP), intent(in) :: x
  ! output: the value of the function
  ! input: the point where to compute the function

  integer, intent(in) :: n
  ! input: the order of the smearing function
  !
  !    here the local variables
  !
  real(DP) :: a, arg, hp, hd
  ! the coefficients a_n
  ! the argument of the exponential
  ! the hermite function
  ! the hermite function

  integer :: i, ni
  ! counter on n values
  ! counter on 2n values

  ! Fermi-Dirac smearing

  if (n.eq. - 99) then
     if (abs (x) .le.36.0) then
        w0gauss = 1.0d0 / (2.0d0 + exp ( - x) + exp ( + x) )
        ! in order to avoid problems for large values of x in the e
     else
        w0gauss = 0.d0
     endif
     return

  endif
  ! cold smearing  (Marzari-Vanderbilt)
  if (n.eq. - 1) then
     arg = min (200.d0, (x - 1.0d0 / sqrt (2.0d0) ) **2)
     w0gauss = sqrtpm1 * exp ( - arg) * (2.0d0 - sqrt ( 2.0d0) * x)
     return

  endif

  !if (n.gt.10 .or. n.lt.0) call errore('py_w0gauss','higher order smearing is untested and unstable',abs(n))

  ! Methfessel-Paxton
  arg = min (200.d0, x**2)
  w0gauss = exp ( - arg) * sqrtpm1
  if (n.eq.0) return
  hd = 0.0d0
  hp = exp ( - arg)
  ni = 0
  a = sqrtpm1
  do i = 1, n
     hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
     ni = ni + 1
     a = - a / (DBLE (i) * 4.0d0)
     hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
     ni = ni + 1
     w0gauss = w0gauss + a * hp
  enddo
end subroutine  py_w0gauss

