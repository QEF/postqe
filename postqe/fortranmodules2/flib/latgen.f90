!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE latgen(ibrav,celldm,a1,a2,a3,omega)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111) 
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  use kinds, only: DP
  implicit none
  integer, intent(in) :: ibrav
  real(DP), intent(inout) :: celldm(6)
  real(DP), intent(inout) :: a1(3), a2(3), a3(3)
  real(DP), intent(out) :: omega
  !
  real(DP), parameter:: sr2 = 1.414213562373d0, &
                        sr3 = 1.732050807569d0
  integer :: i,j,k,l,iperm,ir
  real(DP) :: term, cbya, s, term1, term2, singam, sen
  !
  !  user-supplied lattice vectors
  !
  if (ibrav == 0) then
     if (SQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 ) == 0 )  &
         call errore ('latgen', 'wrong at for ibrav=0', 1)
     if (SQRT( a2(1)**2 + a2(2)**2 + a2(3)**2 ) == 0 )  &
         call errore ('latgen', 'wrong at for ibrav=0', 2)
     if (SQRT( a3(1)**2 + a3(2)**2 + a3(3)**2 ) == 0 )  &
         call errore ('latgen', 'wrong at for ibrav=0', 3)

     if ( celldm(1) /= 0.D0 ) then
     !
     ! ... input at are in units of alat => convert them to a.u.
     !
         a1(:) = a1(:) * celldm(1)
         a2(:) = a2(:) * celldm(1)
         a3(:) = a3(:) * celldm(1)
     else
     !
     ! ... input at are in atomic units: define celldm(1) from a1
     !
         celldm(1) = SQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 )
     end if
     !
  else
     a1(:) = 0.d0
     a2(:) = 0.d0
     a3(:) = 0.d0
  end if
  !
  if (celldm (1) <= 0.d0) call errore ('latgen', 'wrong celldm(1)', ABS(ibrav) )
  !
  !  index of bravais lattice supplied
  !
  if (ibrav == 1) then
     !
     !     simple cubic lattice
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
     !
  else if (ibrav == 2) then
     !
     !     fcc lattice
     !
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
     !
  else if (ABS(ibrav) == 3) then
     !
     !     bcc lattice
     !
     term=celldm(1)/2.d0
     do ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     end do
     IF ( ibrav < 0 ) THEN
        a1(1)=-a1(1)
        a2(2)=-a2(2)
        a3(3)=-a3(3)
     ELSE
        a2(1)=-a2(1)
        a3(1)=-a3(1)
        a3(2)=-a3(2)
     END IF
     !
  else if (ibrav == 4) then
     !
     !     hexagonal lattice
     !
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
     !
  else if (ABS(ibrav) == 5) then
     !
     !     trigonal lattice
     !
     if (celldm (4) <= -0.5_dp .or. celldm (4) >= 1.0_dp) &
          call errore ('latgen', 'wrong celldm(4)', ABS(ibrav))
     !
     term1=sqrt(1.0_dp + 2.0_dp*celldm(4))
     term2=sqrt(1.0_dp - celldm(4))
     !
     IF ( ibrav == 5) THEN
        !     threefold axis along c (001)
        a2(2)=sr2*celldm(1)*term2/sr3
        a2(3)=celldm(1)*term1/sr3
        a1(1)=celldm(1)*term2/sr2
        a1(2)=-a1(1)/sr3
        a1(3)= a2(3)
        a3(1)=-a1(1)
        a3(2)= a1(2)
        a3(3)= a2(3)
     ELSE IF ( ibrav == -5) THEN
        !     threefold axis along (111)
        ! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
        ! does not yield the x,y,z axis, but an equivalent rotated triplet:
        !    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
        ! If you prefer the x,y,z axis as cubic limit, you should modify the
        ! definitions of a1(1) and a1(2) as follows:'
        !    a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
        !    a1(2) = celldm(1)*(term1-term2)/3.0_dp
        ! (info by G. Pizzi and A. Cepellotti)
        !
        a1(1) = celldm(1)*(term1-2.0_dp*term2)/3.0_dp
        a1(2) = celldm(1)*(term1+term2)/3.0_dp
        a1(3) = a1(2)
        a2(1) = a1(3)
        a2(2) = a1(1)
        a2(3) = a1(2)
        a3(1) = a1(2)
        a3(2) = a1(3)
        a3(3) = a1(1)
     END IF
  else if (ibrav == 6) then
     !
     !     tetragonal lattice
     !
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
     !
  else if (ibrav == 7) then
     !
     !     body centered tetragonal lattice
     !
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
     !
  else if (ibrav == 8) then
     !
     !     Simple orthorhombic lattice
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
     !
  else if ( ABS(ibrav) == 9) then
     !
     !     One face (base) centered orthorhombic lattice  (C type)
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', &
                                                                 ABS(ibrav))
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', &
                                                                 ABS(ibrav))
     !
     IF ( ibrav == 9 ) THEN
        !   old PWscf description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) = a1(1) * celldm(2)
        a2(1) = - a1(1)
        a2(2) = a1(2)
     ELSE
        !   alternate description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) =-a1(1) * celldm(2)
        a2(1) = a1(1)
        a2(2) =-a1(2)
     END IF
     a3(3) = celldm(1) * celldm(3)
     !
  else if ( ibrav == 91 ) then
     !
     !     One face (base) centered orthorhombic lattice  (A type)
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a1(1) = celldm(1)
     a2(2) = celldm(1) * celldm(2) * 0.5_DP
     a2(3) = - celldm(1) * celldm(3) * 0.5_DP
     a3(2) = a2(2)
     a3(3) = - a2(3)
     !
  else if (ibrav == 10) then
     !
     !     All face centered orthorhombic lattice
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a2(1) = 0.5d0 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
     !
  else if (ibrav == 11) then
     !
     !     Body centered orthorhombic lattice
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
     !
  else if (ibrav == 12) then
     !
     !     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     if (abs(celldm(4))>=1.d0) call errore ('latgen', 'wrong celldm(4)', ibrav)
     !
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
     !
  else if (ibrav ==-12) then
     !
     !     Simple monoclinic lattice, unique axis: b (more common)
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)',-ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)',-ibrav)
     if (abs(celldm(5))>=1.d0) call errore ('latgen', 'wrong celldm(5)',-ibrav)
     !
     sen=sqrt(1.d0-celldm(5)**2)
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(3)=celldm(1)*celldm(3)*sen
     !
  else if (ibrav == 13) then
     !
     !     One face centered monoclinic lattice unique axis c
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     if (abs(celldm(4))>=1.d0) call errore ('latgen', 'wrong celldm(4)', ibrav)
     !
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = 0.5d0 * celldm(1) 
     a1(3) =-a1(1) * celldm(3)
     a2(1) = celldm(1) * celldm(2) * celldm(4)
     a2(2) = celldm(1) * celldm(2) * sen
     a3(1) = a1(1)
     a3(3) =-a1(3)
  else if (ibrav == -13) then
     !
     !     One face centered monoclinic lattice unique axis b
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)',-ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)',-ibrav)
     if (abs(celldm(5))>=1.d0) call errore ('latgen', 'wrong celldm(5)',-ibrav)
     !
     sen = sqrt( 1.d0 - celldm(5) ** 2 )
     a1(1) = 0.5d0 * celldm(1) 
     a1(2) =-a1(1) * celldm(2)
     a2(1) = a1(1)
     a2(2) =-a1(2)
     a3(1) = celldm(1) * celldm(3) * celldm(5)
     a3(3) = celldm(1) * celldm(3) * sen
     !
  else if (ibrav == 14) then
     !
     !     Triclinic lattice
     !
     if (celldm (2) <= 0.d0) call errore ('latgen', 'wrong celldm(2)', ibrav)
     if (celldm (3) <= 0.d0) call errore ('latgen', 'wrong celldm(3)', ibrav)
     if (abs(celldm(4))>=1.d0) call errore ('latgen', 'wrong celldm(4)', ibrav)
     if (abs(celldm(5))>=1.d0) call errore ('latgen', 'wrong celldm(5)', ibrav)
     if (abs(celldm(6))>=1.d0) call errore ('latgen', 'wrong celldm(6)', ibrav)
     !
     singam=sqrt(1.d0-celldm(6)**2)
     term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
          -celldm(4)**2-celldm(5)**2-celldm(6)**2)
     if (term < 0.d0) call errore &
        ('latgen', 'celldm do not make sense, check your data', ibrav)
     term= sqrt(term/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
     !
  else
     !
     call errore('latgen',' nonexistent bravais lattice',ibrav)
     !
  end if
  !
  !  calculate unit-cell volume omega
  !
  CALL volume (1.0_dp, a1, a2, a3, omega)
  !
  RETURN
  !
END SUBROUTINE latgen
!
!-------------------------------------------------------------------------
SUBROUTINE lat2celldm (ibrav,alat,a1,a2,a3,celldm)
  !-----------------------------------------------------------------------
  !
  !     Returns celldm parameters from lattice vectors
  !     See latgen for definition of celldm and lattice vectors
  !
  use kinds, only: DP
  implicit none
  integer, intent(in) :: ibrav
  real(DP), intent(in) :: alat, a1(3), a2(3), a3(3)
  real(DP), intent(out) :: celldm(6)
  !
  celldm = 0.d0
  !
  SELECT CASE  ( ibrav ) 
  CASE (1:3,-3) 
     celldm(1) = alat
     celldm(2:6) = 0.d0
  CASE (4) 
     celldm(1) = alat
     celldm(2) = 0.d0
     celldm(3) = SQRT( DOT_PRODUCT(a3(:),a3(:)))/alat
     celldm(4:6) = 0.d0
  CASE (5, -5 ) 
     celldm(1)= alat
     celldm(2:3) = 0.d0
     celldm(4) = DOT_PRODUCT(a1(:),a2(:))/(alat**2)
     celldm(5:6) = 0.d0
  CASE (6) 
     celldm(1)= alat 
     celldm(3)= SQRT( DOT_PRODUCT(a3(:),a3(:)))/alat
     celldm(2)= 1.d0
     celldm(4:6) = 0.d0
  CASE (7) 
     celldm(1) = alat
     celldm(3) = ABS(a3(3)/a3(1)) 
     celldm(2)=    0.d0
     celldm(4:6) = 0.d0
  CASE (8)
     celldm(1) = alat
     celldm(2) = SQRT( DOT_PRODUCT (a2(:),a2(:)))/alat
     celldm(3) = SQRT( DOT_PRODUCT (a3(:),a3(:)))/alat 
     celldm(4:6) = 0.d0
  CASE (9, -9 ) 
     celldm(1) = alat
     celldm(2) = ABS ( a1(2)/a1(1))
     celldm(3) = ABS ( a3(3)/2.d0/a1(1))
     celldm(4:6) = 0.d0 
  CASE (10) 
     celldm(1) = alat
     celldm(2) = ABS ( a2(2)/a1(2))
     celldm(3) = ABS ( a1(3)/a1(1))
     celldm(4:6) = 0.d0
  CASE (11) 
     celldm(1) = alat
     celldm(2) = ABS(a1(2)/a1(1))
     celldm(3) = ABS(a1(3)/a1(1))
     celldm(4:6) = 0.d0
  CASE (12, -12) 
     celldm(1) = alat 
     celldm(2) = SQRT( DOT_PRODUCT(a2(:),a2(:))/DOT_PRODUCT(a1(:),a1(:)))
     celldm(3) = SQRT( DOT_PRODUCT(a3(:),a3(:))/DOT_PRODUCT(a1(:),a1(:)))
     celldm(4) = DOT_PRODUCT(a1(:),a2(:))/&
          SQRT(DOT_PRODUCT(a1(:),a1(:))*DOT_PRODUCT(a2(:),a2(:)))
     celldm(5) =  DOT_PRODUCT(a1(:),a3(:))/&
          SQRT(DOT_PRODUCT(a1(:),a1(:))*DOT_PRODUCT(a3(:),a3(:)))
     celldm(6) = 0.d0
  CASE (13) 
     celldm(1) = alat
     celldm(2) = SQRT( DOT_PRODUCT(a2(:),a2(:)))/(2.d0*a1(1))
     celldm(3) = ABS (a3(3)/a3(1))
     celldm(4) = COS( ATAN2( a2(2), a2(1) ) )
     celldm(5:6) = 0.d0
  CASE (14) 
     celldm(1) = alat 
     celldm(2) = SQRT( DOT_PRODUCT(a2(:),a2(:))/DOT_PRODUCT(a1(:),a1(:)))
     celldm(3) = SQRT( DOT_PRODUCT(a3(:),a3(:))/DOT_PRODUCT(a1(:),a1(:)))
     celldm(4) = DOT_PRODUCT(a3(:),a2(:))/SQRT(DOT_PRODUCT(a2(:),a2(:))*&
          DOT_PRODUCT(a3(:),a3(:)))
     celldm(5) = DOT_PRODUCT(a3(:),a1(:))/SQRT(DOT_PRODUCT(a1(:),a1(:))*&
          DOT_PRODUCT(a3(:),a3(:)))
     celldm(6) = DOT_PRODUCT(a1(:),a2(:))/SQRT(DOT_PRODUCT(a2(:),a2(:))*&
          DOT_PRODUCT(a1(:),a1(:)))
  CASE  default  
     celldm(1) = 1.d0
     IF (alat > 0.d0 ) celldm(1) = alat
     celldm (2:6) = 0.d0
  END SELECT

END SUBROUTINE lat2celldm
!
INTEGER FUNCTION abc2ibrav ( a,b,c,cosab,cosac,cosbc ) RESULT (ibrav)
  !
  !  guess of bravais lattice from crystallographic parameters
  !  returns valid value for ibrav, or else ibrav=-1 if failed
  !
  USE kinds, ONLY: dp
  IMPLICIT NONE
  !
  REAL(DP), INTENT (IN) :: a,b,c, cosab, cosac, cosbc 
  !
  !
  IF ( (b == 0.0_dp .AND. c == 0.0_dp) .OR. (b == a .AND. c == a) ) THEN
     ! Case: a=b=c
     
     IF ( cosab == 0.0_dp .AND. cosac == 0.0_dp .AND. cosbc == 0.0_dp ) THEN
        ! Case: simple cubic
        ibrav = 1
     ELSE IF ( cosab/= 0.0_dp .AND. cosac== 0.0_dp .AND. cosbc== 0.0_dp) THEN
        ! Case: trigonal R with threefold axis along <111>
        ibrav =-5
     ELSE IF ( cosab==0.5_dp .AND. cosac == cosab .AND. cosbc == cosab ) THEN
        ! Case: fcc
        ibrav = 2
     ELSE IF ( ABS ( cosab + 1.0_dp/sqrt(3.0_dp) ) < 1.0d-6 .AND. &
               cosac == cosab .AND. cosbc == cosab ) THEN
        ! Case: bcc with symmetric basis
        ibrav =-3
     ELSE
        ! Case: unknown (possibly wrong)
        ibrav =-1
     END IF
     
  ELSE IF ( (b == 0.0_dp .AND. c > 0.0_dp) .OR. (b == a .AND. c /= a) ) THEN
     ! Case: a=b/=c

     IF ( cosab == 0.0_dp .AND. cosac == 0.0_dp .AND. cosbc == 0.0_dp ) THEN
        ! Case: simple tetragonal
        ibrav = 6
     ELSE IF ( cosab ==-0.5_dp .AND. cosac == 0.0_dp .AND. cosbc == 0.0_dp ) THEN
        ! Case: simple hexagonal
        ibrav = 4
     ELSE
        ! Case: unknown (body-centered tetragonal or wrong data)
        ibrav =-1
     END IF
     
  ELSE IF ( b > 0.0_dp .AND. c > 0.0_dp .AND. b /= a .AND. c /= a ) THEN
     ! Case: a/=b/=c
     
     IF ( cosab == 0.0_dp .AND. cosac == 0.0_dp .AND. cosbc == 0.0_dp ) THEN
        ! Case: simple orthorhombic
        ibrav = 8
     ELSE IF ( cosab /= 0.0_dp .AND. cosac == 0.0_dp .AND. cosbc == 0.0_dp ) THEN
        ! Case: monoclinic P, unique axis c
        ibrav = 12
     ELSE IF ( cosab == 0.0_dp .AND. cosac /= 0.0_dp .AND. cosbc == 0.0_dp ) THEN
        ! Case: monoclinic P, unique axis b
        ibrav =-12
     ELSE IF ( cosab /= 0.0_dp .AND. cosac /= 0.0_dp .AND. cosbc /= 0.0_dp ) THEN
        ! Case: triclinic
        ibrav = 14
     ELSE
        ! Case: unknown (base-, face-, body-centered orthorombic,
        !               (base-centered monoclinic)
        ibrav = -1
     END IF
  END IF
  !
END FUNCTION abc2ibrav
!
SUBROUTINE abc2celldm ( ibrav, a,b,c,cosab,cosac,cosbc, celldm )
  !
  !  returns internal parameters celldm from crystallographics ones
  !
  USE kinds, ONLY: dp
  USE constants, ONLY: bohr_radius_angs
  IMPLICIT NONE
  !
  INTEGER,  INTENT (IN) :: ibrav 
  REAL(DP), INTENT (IN) :: a,b,c, cosab, cosac, cosbc 
  REAL(DP), INTENT (OUT) :: celldm(6)
  !
  IF (a <= 0.0_dp) CALL errore('abc2celldm','incorrect lattice parameter (a)',1)
  IF (b <  0.0_dp) CALL errore('abc2celldm','incorrect lattice parameter (b)',1)
  IF (c <  0.0_dp) CALL errore('abc2celldm','incorrect lattice parameter (c)',1)
  IF ( ABS (cosab) > 1.0_dp) CALL errore('abc2celldm', &
                   'incorrect lattice parameter (cosab)',1)
  IF ( ABS (cosac) > 1.0_dp) CALL errore('abc2celldm', &
                   'incorrect lattice parameter (cosac)',1)
  IF ( ABS (cosbc) > 1.0_dp) CALL errore('abc2celldm', &
       'incorrect lattice parameter (cosbc)',1)
  !
  celldm(1) = a / bohr_radius_angs
  celldm(2) = b / a
  celldm(3) = c / a
  !
  IF ( ibrav == 14 ) THEN
     !
     ! ... triclinic lattice
     !
     celldm(4) = cosbc
     celldm(5) = cosac
     celldm(6) = cosab
     !
  ELSE IF ( ibrav ==-12 ) THEN
     !
     ! ... monoclinic P lattice, unique axis b
     !
     celldm(5) = cosac
     !
  ELSE
     !
     ! ... trigonal and monoclinic lattices, unique axis c
     !
     celldm(4) = cosab
     !
  ENDIF
  !
END SUBROUTINE abc2celldm
!
SUBROUTINE cell2abc ( alat, at, a,b,c,cosab,cosac,cosbc )
  !
  !  returns crystallographic parameters a,b,c from lattice vectors
  !
  USE kinds, ONLY: dp
  USE constants, ONLY: bohr_radius_angs
  IMPLICIT NONE
  !
  REAL(DP), INTENT (IN) :: alat, at(3,3)
  REAL(DP), INTENT (OUT) :: a,b,c, cosab, cosac, cosbc 
  REAL(DP) :: norm1, norm2, norm3
  !
  !
  norm1 = SQRT ( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 )
  norm2 = SQRT ( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 )
  norm3 = SQRT ( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 )
  !
  a = alat * norm1 / bohr_radius_angs
  b = alat * norm2 / bohr_radius_angs
  c = alat * norm3 / bohr_radius_angs
  !
  cosab = (at(1,1)*at(1,2) + at(2,1)*at(2,2) + at(3,1)*at(3,2))/norm1/norm2
  cosac = (at(1,1)*at(1,3) + at(2,1)*at(2,3) + at(3,1)*at(3,3))/norm1/norm3
  cosbc = (at(1,3)*at(1,2) + at(2,3)*at(2,2) + at(3,3)*at(3,2))/norm3/norm2
  !
END SUBROUTINE cell2abc
