!
! Copyright (c), 2016-2017, Quantum Espresso Foundation and SISSA (Scuola
! Internazionale Superiore di Studi Avanzati). All rights reserved.
! This file is distributed under the terms of the LGPL-2.1 license. See the
! file 'LICENSE' in the root directory of the present distribution, or
! https://opensource.org/licenses/LGPL-2.1
!
!----------------------------------------------------------------------------
SUBROUTINE errore2( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output: 
  ! ... if ierr <= 0 it does nothing, 
  ! ... if ierr  > 0 it stops.
  !
  ! ...          **** Important note for parallel execution ***
  !
  ! ... in parallel execution unit 6 is written only by the first node;
  ! ... all other nodes have unit 6 redirected to nothing (/dev/null).
  ! ... As a consequence an error not occurring on the first node
  ! ... will be invisible. For T3E and ORIGIN machines, this problem
  ! ... is solved by writing an error message to unit * instead of 6.
  ! ... Whenever possible (IBM SP machines), we write to the standard
  ! ... error, unit 0 (the message will appear in the error files 
  ! ... produced by loadleveler).
  !
 !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routine
    ! the output message
  INTEGER,          INTENT(IN) :: ierr
  CHARACTER(LEN=6) :: cerr
  !
  !
  IF ( ierr <= 0 ) RETURN
  !
  ! ... the error message is written un the "*" unit
  !
  WRITE( cerr, FMT = '(I6)' ) ierr
  WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
  WRITE( UNIT = *, FMT = '(5X,"Error in routine ",A," (",A,"):")' ) &
        TRIM(calling_routine), TRIM(ADJUSTL(cerr))
  WRITE( UNIT = *, FMT = '(5X,A)' ) TRIM(message)
  WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
  !
  !
  WRITE( *, '("     stopping ...")' )
  !
  STOP 1
  !
  RETURN
  !
END SUBROUTINE errore2
!

