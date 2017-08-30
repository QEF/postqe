!
! Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
! Internazionale Superiore di Studi Avanzati). All rights reserved.
! This file is distributed under the terms of the LGPL-2.1 license. See the
! file 'LICENSE' in the root directory of the present distribution, or
! https://opensource.org/licenses/LGPL-2.1
!
subroutine generate_gg_list(nrrr, nr1, nr2, nr3 , bg1, bg2, bg3, g, gg , mill)
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi = 4.d0*atan(1.d0) , fpi = 4.d0*pi, eps8 = 1.d-8

  integer, intent(in)   :: nr1, nr2, nr3, nrrr 
  real(dp),intent(in)   :: bg1(3), bg2(3), bg3(3)
  real(dp), intent(out) :: gg(nrrr) , g(nrrr,3) 
  integer (dp), intent(out) :: mill(nrrr,3) 
  real(dp)              :: vec(3)
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
end subroutine generate_gg_list


