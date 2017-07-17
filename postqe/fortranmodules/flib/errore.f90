subroutine errore (loc_msg, case_msg, error_code) 
implicit none 
character(len=*),intent(in)    :: loc_msg, case_msg
integer, intent(in)            :: error_code
if (error_code /=0) then 
  print *,loc_msg, case_msg, error_code
  stop
end if

end subroutine errore
