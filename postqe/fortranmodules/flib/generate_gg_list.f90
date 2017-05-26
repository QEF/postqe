
subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !     first the input variables
  !
  implicit none
  integer, parameter    :: dp = selected_real_kind(14,200) 
  real(DP),intent(in) :: a1 (3), a2 (3), a3 (3) 
  real(dp),intent(out) :: b1 (3), b2 (3), b3 (3)
  ! input: first direct lattice vector
  ! input: second direct lattice vector
  ! input: third direct lattice vector
  ! output: first reciprocal lattice vector
  ! output: second reciprocal lattice vector
  ! output: third reciprocal lattice vector
  !
  !   then the local variables
  !
  real(DP) :: den, s
  ! the denominator
  ! the sign of the permutations
  integer :: iperm, i, j, k, l, ipol
  ! counter on the permutations
  !\
  !  Auxiliary variables
  !/
  !
  ! Counter on the polarizations
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.d0
100 do iperm = 1, 3
     den = den + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s
  if (s.lt.0.d0) goto 100
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
     b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
     b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
     l = i
     i = j
     j = k
     k = l
  enddo
  return
end subroutine recips


subroutine generate_gg_list( nrrr, nr1, nr2, nr3 , bg1, bg2, bg3, g, gg , mill) 
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


