
subroutine get_igtongl(ngm, gg, igtongl, ngl) 
implicit none
integer, parameter  :: dp = selected_real_kind(14,200)
real(dp), parameter :: eps8 = 1.d-8

integer , intent (in) :: ngm
real(dp), intent (in) ::   gg(ngm)
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
integer , parameter :: dp = selected_real_kind(14,200) 
real(dp), parameter :: eps8 = 1.d-8
integer, intent(in)  :: ngm, ngl
real(dp), intent(in) :: gg(ngm)
real(dp), intent(out) :: gl(ngl)
integer igl, ng
igl = 1 
gl(1) = gg(1)
do ng = 1, ngm 
   if ( gg(ng) > gg( ng-1)+eps8) igl = igl+1 
   if (igl > ngl) stop 'igl > ngl ?' 
   gl(igl) = gg( ng)
end do  

end subroutine get_gl

