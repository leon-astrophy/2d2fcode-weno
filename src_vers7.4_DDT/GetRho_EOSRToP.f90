!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the pressure given the density    !
! of a completely degenerate ideal fermi gas. It is used only !
! in solving for the initial star model. Do not confuse it    !
! with the subroutine findpressure, which aims to update the  !
! pressure after hydrodynamic evolution                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_EOSRTOP (ini_p, ini_rho, gs, mb, me, ye, type)
USE DEFINITION
IMPLICIT NONE

! real parameter, including input density ... etc !
REAL (DP), INTENT (IN) :: ini_rho, gs, mb, me, ye

! Integer parameter !
INTEGER, INTENT (IN) :: type

! The output pressure !
REAL (DP), INTENT (OUT) :: ini_p

! dimensionless fermi momentum !
REAL (DP) :: dlfmmo

! we need the dimensionless fermi momentum first !
CALL FERMIMO (dlfmmo, ini_rho, type)

! We find the pressure according to the type, 1 is DM, 2 is NM !
IF (type == 1) THEN
	IF(fermieosdm_flag == 1) THEN 
		IF (dlfmmo <= 1.0E-2_DP) THEN
			ini_p = a_max1*small_pressure(dlfmmo)	
		ELSE
			ini_p = a_max1*large_pressure(dlfmmo)
		END IF
	ELSEIF (polyeosdm_flag == 1) THEN
		ini_p = k_1 * ini_rho ** (gamma1)
	END IF
ELSEIF (type == 2) THEN
	IF (fermieosnm_flag == 1) THEN
		IF (dlfmmo <= 1.0E-2_DP) THEN
			ini_p = a_max2*small_pressure(dlfmmo)	
		ELSE
			ini_p = a_max2*large_pressure(dlfmmo)
		END IF
	ELSEIF (polyeosnm_flag == 1) THEN
		ini_p = k_2 * ini_rho ** (gamma2)
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

contains

	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE