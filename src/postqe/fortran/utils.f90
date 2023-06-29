!
! Copyright (c), 2016-2021, Quantum Espresso Foundation and SISSA (Scuola
! Internazionale Superiore di Studi Avanzati). All rights reserved.
! This file is distributed under the terms of the LGPL-2.1 license. See the
! file 'LICENSE' in the root directory of the present distribution, or
! https://opensource.org/licenses/LGPL-2.1
!
! A collection of utility subroutines better written in Fortran than in Python.
!

subroutine get_igtongl(ngm, gg, igtongl, ngl)
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
end subroutine get_igtongl


subroutine get_gl(ngm, gg, ngl, gl)
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
end subroutine get_gl


subroutine get_gg_list( nrrr, nr1, nr2, nr3 , bg1, bg2, bg3, g, gg , mill)
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
    if ( ii <= nr1/2) then
        i =  ii 
    else 
        i =  ii - nr1 
    end if
    do ij = 0, nr2 -1
        if ( ij <= nr2/2 ) then
          j = ij 
        else
          j = ij - nr2 
        end if  
        do ik = 0, nr3 - 1  
          if ( ik <= nr3/2 ) then
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
end subroutine get_gg_list
