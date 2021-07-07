SUBROUTINE FINDCONST
USE DEFINITION
IMPLICIT NONE

! Find constant !
a_max1 = (gs1*me1**4)/(4.8D1*pi**2*hbar**3)
a_max2 = (gs2*me2**4)/(4.8D1*pi**2*hbar**3)
b_max1 = (mb1/ye1)*((gs1*me1**3)/(6.0E0_DP*pi**2*hbar**3))
b_max2 = (mb2/ye2_old + me2)*((gs2*me2**3)/(6.0E0_DP*pi**2*hbar**3))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the dimensionless fermi momentum used in ideal fermi gas EOS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FERMIMO (dlfmmo, rho_in, type)
USE DEFINITION
IMPLICIT NONE

! Real Variables !
! the inputed parameters !
INTEGER, INTENT (IN) :: type

! the output dimensionless fermi momentum !
REAL (DP), INTENT (IN) :: rho_in
REAL (DP), INTENT (OUT) :: dlfmmo 

! find the fermi momentum !
IF(type == 1) THEN
	dlfmmo = (rho_in/b_max1)**(1.0E0_DP/3.0E0_DP)
ELSEIF(type == 2) THEN
	dlfmmo = (rho_in/b_max2)**(1.0E0_DP/3.0E0_DP)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the derivative of rho wtih respect to x !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDXDRHO (dxdrho, rho_in, type)
USE DEFINITION
IMPLICIT NONE

! Real Variables !
! the inputed parameters !
REAL (DP), INTENT (IN) :: rho_in

! the inputed parameters !
INTEGER, INTENT (IN) :: type

! the output derivative !
REAL (DP), INTENT (OUT) :: dxdrho

! find the fermi momentum !
IF(type == 1) THEN
	dxdrho = 1.0D0/3.0D0/(rho_in**2*b_max1)**(1.0E0_DP/3.0E0_DP)
ELSEIF(type == 2) THEN
	dxdrho = 1.0D0/3.0D0/(rho_in**2*b_max2)**(1.0E0_DP/3.0E0_DP)
END IF

END SUBROUTINE