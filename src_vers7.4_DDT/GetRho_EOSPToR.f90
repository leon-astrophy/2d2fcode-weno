!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine aims at finding the density corresponds !
! to a given pressure in the case of ideal degenerate     !
! fermi gas equation of state by interpolation            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_EOSPTOR (ini_rho, ini_p, eostable, eosline, type)
USE DEFINITION
IMPLICIT NONE

! the input pressure !
REAL (DP), INTENT (IN) :: ini_p

! Integer parameter !
INTEGER, INTENT (IN) :: type

! the outputed density !
REAL (DP), INTENT (OUT) :: ini_rho

! related to which matter eos are we interpolating !
REAL (DP), DIMENSION (eoslineno, 2), INTENT (IN) :: eostable

! distinguish between NM and DM !
INTEGER, INTENT (INOUT) :: eosline

! Other essential variables !
INTEGER :: i
REAL (DP) :: x1, x2, x3, x4, y1, y2, y3, y4

! We find the density according to the type, 1 is DM, 2 is NM !
IF (type == 1) THEN
	IF (fermieosdm_flag == 1) THEN

		! We locate the first where the pressure located !
		i = eosline	

		DO
			IF (ini_p > eostable (i, 1)) THEN
				i = i + 1
			ELSE
				IF (ini_p > eostable (i - 1, 1)) THEN
					x1 = eostable (i - 2, 1)
					x2 = eostable (i - 1, 1)
					x3 = eostable (i, 1)
					x4 = eostable (i + 1, 1)
					y1 = eostable (i - 2, 2)
					y2 = eostable (i - 1, 2)
					y3 = eostable (i, 2)
					y4 = eostable (i + 1, 2)
					eosline = i
					EXIT
				ELSE
					i = i - 1
				END IF
			END IF
		END DO

		! We do the interpolation !
		ini_rho = ((ini_p - x2) * (ini_p - x3) * (ini_p - x4)) / ((x1 - x2) * (x1 - x3) * (x1 - x4)) * y1 &
			+ ((ini_p - x1) * (ini_p - x3) * (ini_p - x4)) / ((x2 - x1) * (x2 - x3) * (x2 - x4)) * y2 &
			+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x4)) / ((x3 - x1) * (x3 - x2) * (x3 - x4)) * y3 &
			+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x3)) / ((x4 - x1) * (x4 - x2) * (x4 - x3)) * y4
	ELSEIF (polyeosdm_flag == 1) THEN
		ini_rho = (ini_p / k_1) ** (1.0E0_DP/gamma1)
	END IF
ELSEIF (type == 2) THEN
	IF (fermieosnm_flag == 1) THEN

		! We locate the first where the pressure located !
		i = eosline
		
		DO
			IF (ini_p > eostable (i, 1)) THEN
				i = i + 1
			ELSE
				IF (ini_p > eostable (i - 1, 1)) THEN
					x1 = eostable (i - 2, 1)
					x2 = eostable (i - 1, 1)
					x3 = eostable (i, 1)
					x4 = eostable (i + 1, 1)
					y1 = eostable (i - 2, 2)
					y2 = eostable (i - 1, 2)
					y3 = eostable (i, 2)
					y4 = eostable (i + 1, 2)
					eosline = i
					EXIT
				ELSE
					i = i - 1
				END IF
			END IF
		END DO

		! We do the interpolation !
		ini_rho = ((ini_p - x2) * (ini_p - x3) * (ini_p - x4)) / ((x1 - x2) * (x1 - x3) * (x1 - x4)) * y1 &
			+ ((ini_p - x1) * (ini_p - x3) * (ini_p - x4)) / ((x2 - x1) * (x2 - x3) * (x2 - x4)) * y2 &
			+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x4)) / ((x3 - x1) * (x3 - x2) * (x3 - x4)) * y3 &
			+ ((ini_p - x1) * (ini_p - x2) * (ini_p - x3)) / ((x4 - x1) * (x4 - x2) * (x4 - x3)) * y4

	ELSEIF (polyeosnm_flag == 1) THEN
		ini_rho = (ini_p / k_2) ** (1.0E0_DP/gamma2)
	END IF
END IF

END SUBROUTINE