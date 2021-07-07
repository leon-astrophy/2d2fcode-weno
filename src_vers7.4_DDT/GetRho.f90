!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a two fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-5 method   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_2F
USE DEFINITION
USE HELMEOS_MODULE
IMPLICIT NONE

! The number of hydrostatic equation and the improved accuracy !
INTEGER, PARAMETER :: no_of_eq_ini = 6, more = 10 ** ini_acc

! The extra array arising from the extra accuracy !
INTEGER, PARAMETER :: length_morestep = max(length_step_r_1, length_step_r_2)*more

! Dummy integers
INTEGER :: i, j, k 

! Dummy variables !
REAL (DP) :: dummy, rho1_min, rho2_min

! The grid number of the stellar surface
INTEGER :: r_grid_more1, r_grid_more2

! Dummy arrays for NM and DM density to be interpolated !	
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: den1, den2

! Temporal distance arrays !
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: rtemp

! This is necessary for any finite temperature EOS since the deriative of pressure plays a role !
REAL (DP), DIMENSION (-4: length_morestep + 5) :: dpdrho2_temp

! Center and atmosphere pressure
REAL (DP) :: p1_c, p2_c, p1_a!, p2_a

! Initial density and pressure
REAL (DP) :: ini_rho1, ini_rho2, ini_p1, ini_p2

! Position
REAL (DP) :: x

! The RK5 variables
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: y_zero, y_one, y_two, y_three, y_four, y_five
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: der_one, der_two, der_three, der_four, der_five, der_new
REAL (DP), DIMENSION (1 : no_of_eq_ini, -4 : length_morestep + 5) :: y

! Quantity related to isotope composition if !
! you are using variable compositions of star !
REAL (DP), DIMENSION (totalion) :: xiso_ini
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! special flag for switching composition of envelop !
WRITE(*,*) 'In GetRho2F'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign minimum density !
rho1_min = 1.0D-4*rho1_a
rho2_min = 1.0D-4*rho2_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We read the EOS table for DM !
IF (fermieosdm_flag == 1) THEN
	OPEN (UNIT = 99, FILE = 'EOS_Table1.eos', STATUS = 'OLD')
	DO i = 1, eoslineno
		READ (99, *) eostable1 (i, 1), eostable1 (i, 2)
	END DO
	eosline1 = 1
	CLOSE (99)
END IF

! We read the EOS table for NM !
IF (fermieosnm_flag == 1) THEN
	OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'OLD')
	DO i = 1, eoslineno
		READ (100, *) eostable2 (i, 1), eostable2 (i, 2)
	END DO
	CLOSE (100)
	eosline2 = 1
END IF

! We convert density at center and atmosphere to the corresponding pressure !
CALL GETRHO_EOSRTOP (p1_c, rho1_c, gs1, mb1, me1, ye1, 1)
CALL GETRHO_EOSRTOP (p1_a, rho1_a, gs1, mb1, me1, ye1, 1)

! Finite temperature EOS !
IF(helmeos_flag == 1) THEN

	! We assign initial isotope composition accordingly !
	! If you are using helmeos EOS in construction !
	! It is because the realistic EOS not only depends on !
	! Ye but also Abar and Zbar !	
	xiso_ini = 0.0D0
	xiso_ini(che4) = xhe4_ini1
	xiso_ini(cc12) = xc12_ini1
	xiso_ini(co16) = xo16_ini1
	xiso_ini(cne20) = xne20_ini1
	xiso_a = xiso_ini

	! Now convert the composition into mean atomic and mass number !
	CALL private_helmeos_azbar(xiso_ini, abar_ini, zbar_ini, ye2_ini)

	! assign atmospheric electron fraction !
	ye_a = ye2_ini
ELSE
	ye_a = ye2_old
END IF

! We choose which table is needed to read according to the chosen EOS !
IF(helmeos_flag == 1) THEN
	
	! Read the helmeos eos table for interpolation !
	CALL read_helm_table()

	! We convert density at center and atmosphere to the corresponding pressure !
	! Note that at center the deriative of density should be zero !
	CALL HELMEOS_RtoP(rho2_c, temp_ini, abar_ini, zbar_ini, ye2_ini, p2_c, &
			 dpdrho2_temp(0), dummy, dummy, dummy, dummy)
	CALL HELMEOS_RtoP(rho2_a, temp_ini, abar_ini, zbar_ini, ye2_ini, p2_a, &
			 dummy, dummy, dummy, dummy, dummy)

ELSE

	! We convert density at center and atmosphere to the corresponding pressure !
	CALL GETRHO_EOSRTOP (p2_c, rho2_c, gs2, mb2, me2, ye2_old, 2)
	CALL GETRHO_EOSRTOP (p2_a, rho2_a, gs2, mb2, me2, ye2_old, 2)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Core DM rho and p: ', rho1_c, p1_c
WRITE(*,*) 'Core NM rho and p: ', rho2_c, p2_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign distance !
DO j = 1, length_morestep 
	rtemp(j) = dxmore*DBLE(j - 1/2)
END DO

! We assign the value of y at center !
y (1, 0) = (4.0D0/3.0D0)*pi*rtemp(1)**3*rho1_c
y (2, 0) = p1_c - 0.5D0*y(1,0)*rho1_c/rtemp(1)
y (3, 0) = rho1_c
y (4, 0) = (4.0D0/3.0D0)*pi*rtemp(1)**3*rho2_c
y (5, 0) = p2_c - 0.5D0*y(4,0)*rho2_c/rtemp(1)
y (6, 0) = rho2_c

! initialize the integer parameter !
r_grid_more1 = length_morestep
r_grid_more2 = length_morestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO j = 0, length_morestep - 1	

	! Update the value of x and y !
	DO i = 1, no_of_eq_ini
		y_zero (i) = y (i, j)
	END DO

	x = (DBLE (j) + 0.5D0) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the first step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_one, x, y_zero, no_of_eq_ini)
	
	! We artificially make the value to zero to avoid singularity at center !
	IF (j == 0) THEN
		der_one (2) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	! If the density reach atmospheric values, no changes in all the quantity !
	IF (y_zero (3) == rho1_min) THEN
		der_one (1) = 0.0E0_DP
		der_one (2) = 0.0E0_DP
	END IF
	IF (y_zero (6) == rho2_min) THEN
		der_one (4) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_one (i) = y_zero (i) + (1.0E0_DP / 4.0E0_DP) * dxmore * der_one (i)
	END DO

	! We determine whether the pressure reached atmospheric pressure !
	IF (j < r_grid_more1) THEN
		ini_p1 = y_one (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_one (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	! Assign density according to pressure !
	If (ini_p1 == p1_a) THEN
		y_one (3) = rho1_min
	ELSE
		y_one (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_one (6) = rho2_min
	ELSE
		y_one (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the second step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_two, x, y_one, no_of_eq_ini)

	IF (y_one (3) == rho1_min) THEN
		der_two (1) = 0.0E0_DP
		der_two (2) = 0.0E0_DP
	END IF
	IF (y_one (6) == rho2_min) THEN
		der_two (4) = 0.0E0_DP
		der_two (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_two (i) = y_zero (i) + (1.0E0_DP / 8.0E0_DP) * dxmore * der_one (i) &
				+ (1.0E0_DP / 8.0E0_DP) * dxmore * der_two (i)
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_two (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_two (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_two (3) = rho1_min
	ELSE
		y_two (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_two (6) = rho2_min
	ELSE
		y_two (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the third step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_three, x, y_two, no_of_eq_ini)
	
	IF (y_two (3) == rho1_min) THEN
		der_three (1) = 0.0E0_DP
		der_three (2) = 0.0E0_DP
	END IF
	IF (y_two (6) == rho2_min) THEN
		der_three (4) = 0.0E0_DP
		der_three (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_three (i) = y_zero (i) - (1.0E0_DP / 2.0E0_DP) * dxmore * der_two (i) &
			+ dxmore * der_three (i)
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_three (3) = rho1_min
	ELSE
		y_three (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_three (6) = rho2_min
	ELSE
		y_three (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 2.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the forth step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_four, x, y_three, no_of_eq_ini)
	
	IF (y_three (3) == rho1_min) THEN
		der_four (1) = 0.0E0_DP
		der_four (2) = 0.0E0_DP
	END IF
	IF (y_three (6) == rho2_min) THEN
		der_four (4) = 0.0E0_DP
		der_four (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_four (i) = y_zero (i) + (3.0E0_DP / 1.6E1_DP) * dxmore * der_one (i) &
			+ (9.0E0_DP / 1.6E1_DP) * dxmore * der_four (i) 
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_four (3) = rho1_min
	ELSE
		y_four (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_four(6) = rho2_min
	ELSE
		y_four (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (3.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the fifth step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_five, x, y_four, no_of_eq_ini)
	
	IF (y_four (3) == rho1_min) THEN
		der_five (1) = 0.0E0_DP
		der_five (2) = 0.0E0_DP
	END IF
	IF (y_four (6) == rho2_min) THEN
		der_five (4) = 0.0E0_DP
		der_five (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_five (i) = y_zero (i) - (3.0E0_DP / 7.0E0_DP) * dxmore * der_one (i) &
			+ (2.0E0_DP / 7.0E0_DP) * dxmore * der_two (i) &
			+ (1.2E1_DP / 7.0E0_DP) * dxmore * der_three (i) &
			- (1.2E1_DP / 7.0E0_DP) * dxmore * der_four (i) &
			+ (8.0E0_DP / 7.0E0_DP) * dxmore * der_five (i)
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_five (3) = rho1_min
	ELSE
		y_five (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_five (6) = rho2_min
	ELSE
		y_five (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the final step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_new, x, y_five, no_of_eq_ini)

	IF (y_five (3) == rho1_min) THEN
		der_new (1) = 0.0E0_DP
		der_new (2) = 0.0E0_DP
	END IF
	IF (y_five (6) == rho2_min) THEN
		der_new (4) = 0.0E0_DP
		der_new (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y (i, j + 1) = y_zero (i) + (7.0E0_DP * der_one (i) + 3.2E1_DP * der_three (i) & 
				+ 1.2E1_DP * der_four (i) + 3.2E1_DP * der_five (i) & 
				+ 7.0E0_DP * der_new (i)) * dxmore / 9.0E1_DP
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y (2, j + 1)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y (5, j + 1)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y (3, j + 1) = rho1_min
	ELSE
		y (3, j + 1) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y (6, j + 1) = rho2_min
	ELSE
		y (6, j + 1) = ini_rho2
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_RtoP(y(6, j+1), temp_ini, abar_ini, zbar_ini, ye2_ini, dummy, dpdrho2_temp(j+1), dummy, dummy, dummy, dummy)
	END IF
	
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign atmospheric density !
DO j = r_grid_more1 + 1, length_morestep
	y (3, j) = rho1_min
END DO
DO j = r_grid_more2 + 1, length_morestep
	y (6, j) = rho2_min
END DO

! Assign density !
den1 = y (3, :)
den2 = y (6, :)

! Now we assign the value of density and mass to the original array !
DO j = 1, length_step_r_2
	DO k = 1, length_step_z_2
		IF(rad2(j,k) > rtemp(length_morestep)) THEN
			rho2(j,k) = rho2_min
		ELSE	
			DO i = 1, length_morestep
				IF(rtemp(i) == rad2(j,k)) THEN
					rho2(j,k) = den2(i)
					EXIT
				ELSEIF(rtemp(i) > rad2(j,k)) THEN
					CALL AKIMA(rtemp(i-3), rtemp(i-2), rtemp(i-1), rtemp(i), rtemp(i+1), rtemp(i+2), & 
					den2(i-3), den2(i-2), den2(i-1), den2(i), den2(i+1), den2(i+2), rad2(j,k), rho2(j,k))
					EXIT
				END IF
			END DO
		END IF
	END DO
END DO
DO j = 1, length_step_r_1
	DO k = 1, length_step_z_1
		IF(rad1(j,k) > rtemp(length_morestep)) THEN
			rho1(j,k) = rho1_min
		ELSE	
			DO i = 1, length_morestep
				IF(rtemp(i) == rad1(j,k)) THEN
					rho1(j,k) = den1(i)
					EXIT
				ELSEIF(rtemp(i) > rad1(j,k)) THEN
					CALL AKIMA(rtemp(i-3), rtemp(i-2), rtemp(i-1), rtemp(i), rtemp(i+1), rtemp(i+2), & 
					den1(i-3), den1(i-2), den1(i-1), den1(i), den1(i+1), den1(i+2), rad1(j,k), rho1(j,k))
					EXIT
				END IF
			END DO
		END IF
	END DO
END DO

! Check the density !
DO j = 1, length_step_r_2
	DO k = 1, length_step_z_2
		IF(rho2(j,k) < rho2_a) THEN
			rho2(j,k) = rho2_a
		END IF
	END DO
END DO
DO j = 1, length_step_r_1
	DO k = 1, length_step_z_1
		IF(rho1(j,k) < rho1_a) THEN
			rho1(j,k) = rho1_a
		END IF
	END DO
END DO

! Boundaries !
CALL boundary1D_DMFULL (rho1, even)
CALL boundary1D_NMFULL (rho2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
m_cell = 0.0D0

! Now we assign the value of density to the original array !
DO j = 1, length_step_r_1
	DO k = 1, j
		If (rho1(k, 1) > rho1_a) THEN
			m_cell (j) = m_cell(j) + rho1(k, 1)*(4.0D0/3.0D0)*pi*(rF1(k)**3 - rF1(k-1)**3)
		ELSE
			CYCLE
		END IF
	END DO
ENDDO

! Assign mass centered coordinate by interpolation !
DO j = 1, length_step_r_1
	m_r1(j) = 0.5D0*(m_cell (j) + m_cell (j-1))
ENDDO

! Initialize !
m_cell = 0.0D0

! Now we assign the value of density to the original array !
DO j = 1, length_step_r_2
	DO k = 1, j
		If (rho2(k, 1) > rho2_a) THEN
			m_cell (j) = m_cell(j) + rho2(k, 1)*(4.0D0/3.0D0)*pi*(rF2(k)**3 - rF2(k-1)**3)
		ELSE
			CYCLE
		END IF
	END DO
ENDDO

! Assign mass centered coordinate by interpolation !
DO j = 1, length_step_r_2
	m_r2(j) = 0.5D0*(m_cell (j) + m_cell (j-1))
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now of course, we need boundary condition for composition !
IF (helmeos_flag == 1) THEN

	! We first assign the appropriate value to the arrays !
	DO k = 1, length_step_z_2, 1
		DO j = 1, length_step_r_2, 1
			if (rho2(j,k) > 0.0D0) then	
            			abar2(j,k) = abar_ini
            			zbar2(j,k) = zbar_ini
            			xiso(j,k,:) = xiso_ini(:)
				ye2(j,k) = ye2_ini !zbar_ini / abar_ini
      	    			temp2(j,k) = temp_ini
			END IF
		END DO
	END DO

	! Copy the arrays to the ghost shell !
	CALL BOUNDARY2D_X (xiso)
	CALL BOUNDARY1D_NMFULL (abar2, even)       
	CALL BOUNDARY1D_NMFULL (zbar2, even)
	CALL BOUNDARY1D_NMFULL (temp2, even)
	CALL BOUNDARY1D_NMFULL (ye2, even)  

END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a two fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-5 method   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_1F
USE DEFINITION
USE HELMEOS_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid size for the initial model !
! Hint: One set it to 1% of dx for accuracy !
!REAL (DP) :: dxmore = dx_ini / (10.0D0**ini_acc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The number of hydrostatic equation and the improved accuracy !
INTEGER, PARAMETER :: no_of_eq_ini = 3, more = 10 ** ini_acc

! The extra array arising from the extra accuracy !
INTEGER, PARAMETER :: length_morestep = length_step_r_2 * more

! Dummy integers
INTEGER :: i, j, k 

! Dummy variables !
REAL (DP) :: dummy, rho2_min

! The grid number of the stellar surface
INTEGER :: r_grid_more2

! Center and atmosphere pressure
REAL (DP) :: p2_c!, p2_a

! Initial density and pressure
REAL (DP) :: ini_rho2, ini_p2

! Position
REAL (DP) :: x

! Dummy arrays for NM and DM density to be interpolated !	
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: den2

! Temporal distance arrays !
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: rtemp

! The RK5 variables
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: y_zero, y_one, y_two, y_three, y_four, y_five
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: der_one, der_two, der_three, der_four, der_five, der_new
REAL (DP), DIMENSION (1 : no_of_eq_ini, -4 : length_morestep + 5) :: y		

! This is necessary for any EOS since the deriative of pressure plays a role !
REAL (DP), DIMENSION (-4: length_morestep + 5) :: dpdrho2_temp

! Quantity related to isotope composition if !
! you are using variable compositions of star !
REAL (DP), DIMENSION (totalion) :: xiso_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Special flag for switching composition of envelop
WRITE(*,*) 'In GetRho1F'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign minmium density !
rho2_min = 1.0D-4*rho2_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We read the EOS table for NM !
IF (fermieosnm_flag == 1) THEN
	OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'OLD')
	DO i = 1, eoslineno
		READ (100, *) eostable2 (i, 1), eostable2 (i, 2)
	END DO
	CLOSE (100)
	eosline2 = 1
END IF

IF(helmeos_flag == 1) THEN

	! We assign initial isotope composition accordingly !
	! If you are using helmeos EOS in construction !
	! It is because the realistic EOS not only depends on !
	! Ye but also Abar and Zbar !	
	xiso_ini = 0.0D0
	xiso_ini(che4) = xhe4_ini1
	xiso_ini(cc12) = xc12_ini1
	xiso_ini(co16) = xo16_ini1
	xiso_ini(cne20) = xne20_ini1
	xiso_a = xiso_ini

	! Now convert the composition into mean atomic and mass number !
	CALL private_helmeos_azbar(xiso_ini, abar_ini, zbar_ini, ye2_ini)

	! assign atmospheric electron fraction !
	ye_a = ye2_ini
ELSE
	ye_a = ye2_old
END IF

! We choose which table is needed to read according to the chosen EOS !
IF(helmeos_flag == 1) THEN
	
	! Read the helmeos eos table for interpolation !
	CALL read_helm_table()	

	! We convert density at center and atmosphere to the corresponding pressure !
	! Note that at center the deriative of density should be zero !
	CALL HELMEOS_RtoP(rho2_c, temp_ini, abar_ini, zbar_ini, ye2_ini, p2_c, &
			 dpdrho2_temp(0), dummy, dummy, dummy, dummy)
	CALL HELMEOS_RtoP(rho2_a, temp_ini, abar_ini, zbar_ini, ye2_ini, p2_a, &
			 dummy, dummy, dummy, dummy, dummy)

ELSE

	! We convert density at center and atmosphere to the corresponding pressure !
	CALL GETRHO_EOSRTOP (p2_c, rho2_c, gs2, mb2, me2, ye2_old, 2)
	CALL GETRHO_EOSRTOP (p2_a, rho2_a, gs2, mb2, me2, ye2_old, 2)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Core NM rho and p: ', rho2_c, p2_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign distance !
DO j = 1, length_morestep 
	rtemp(j) = dxmore*DBLE(j - 1/2)
END DO

! We assign the value of y at center !
y (1, 0) = (4.0D0/3.0D0)*pi*rtemp(1)**3*rho2_c
y (2, 0) = p2_c - 0.5D0*y (1, 0)*rho2_c/rtemp(1)
y (3, 0) = rho2_c

! initialize the integer parameter !
r_grid_more2 = length_morestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO j = 0, length_morestep - 1	
	
	! Update the value of x and y !
	DO i = 1, no_of_eq_ini
		y_zero (i) = y (i, j)
	END DO

	x = (DBLE (j) + 0.5D0) * dxmore
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the first step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_one, x, y_zero, no_of_eq_ini)
	
	! We artificially make the value to zero to avoid singularity at center !
	IF (j == 0) THEN
		der_one (2) = 0.0E0_DP
	END IF

	! If the density reach atmospheric values, no changes in all the quantity !
	IF (y_zero (3) == rho2_min) THEN
		der_one (1) = 0.0E0_DP
		der_one (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_one (i) = y_zero (i) + (1.0E0_DP / 4.0E0_DP) * dxmore * der_one (i)
	END DO

	! We determine whether the pressure reached atmospheric pressure !
	IF (j < r_grid_more2) THEN
		ini_p2 = y_one (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF
	
	! Assign density according to pressure !
	If (ini_p2 == p2_a) THEN
		y_one (3) = rho2_min
	ELSE
		y_one (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the second step in RK-5 !	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_two, x, y_one, no_of_eq_ini)

	IF (y_one (3) == rho2_min) THEN
		der_two (1) = 0.0E0_DP
		der_two (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_two (i) = y_zero (i) + (1.0E0_DP / 8.0E0_DP) * dxmore * der_one (i) &
				+ (1.0E0_DP / 8.0E0_DP) * dxmore * der_two (i)
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_two (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_two (3) = rho2_min
	ELSE
		y_two (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the third step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_three, x, y_two, no_of_eq_ini)

	IF (y_two (3) == rho2_min) THEN
		der_three (1) = 0.0E0_DP
		der_three (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_three (i) = y_zero (i) - (1.0E0_DP / 2.0E0_DP) * dxmore * der_two (i) &
			+ dxmore * der_three (i)
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_three (3) = rho2_min
	ELSE
		y_three (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 2.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the forth step in RK-5 !	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_four, x, y_three, no_of_eq_ini)

	IF (y_three (3) == rho2_min) THEN
		der_four (1) = 0.0E0_DP
		der_four (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_four (i) = y_zero (i) + (3.0E0_DP / 1.6E1_DP) * dxmore * der_one (i) &
			+ (9.0E0_DP / 1.6E1_DP) * dxmore * der_four (i) 
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_two (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_four (3) = rho2_min
	ELSE
		y_four (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (3.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the fifth step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_five, x, y_four, no_of_eq_ini)

	IF (y_four (3) == rho2_min) THEN
		der_five (1) = 0.0E0_DP
		der_five (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_five (i) = y_zero (i) - (3.0E0_DP / 7.0E0_DP) * dxmore * der_one (i) &
			+ (2.0E0_DP / 7.0E0_DP) * dxmore * der_two (i) &
			+ (1.2E1_DP / 7.0E0_DP) * dxmore * der_three (i) &
			- (1.2E1_DP / 7.0E0_DP) * dxmore * der_four (i) &
			+ (8.0E0_DP / 7.0E0_DP) * dxmore * der_five (i)
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_five (3) = rho2_min
	ELSE
		y_five (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the final step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_new, x, y_five, no_of_eq_ini)

	IF (y_five (3) == rho2_min) THEN
		der_new (1) = 0.0E0_DP
		der_new (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y (i, j + 1) = y_zero (i) + (7.0E0_DP * der_one (i) + 3.2E1_DP * der_three (i) & 
				+ 1.2E1_DP * der_four (i) + 3.2E1_DP * der_five (i) & 
				+ 7.0E0_DP * der_new (i)) * dxmore / 9.0E1_DP
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y (2, j + 1)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	IF (helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp_ini, abar_ini, zbar_ini, ye2_ini, ini_rho2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y (3, j + 1) = rho2_min
	ELSE
		y (3, j + 1) = ini_rho2
	END IF
	
	if(helmeos_flag == 1) then
		call HELMEOS_RtoP(y(3, j+1), temp_ini, abar_ini, zbar_ini, ye2_ini, dummy, dpdrho2_temp(j+1), dummy, dummy, dummy, dummy)
	END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign atmospheric density !
DO j = r_grid_more2 + 1, length_morestep
	y (3, j) = rho2_min
END DO

! Assign density !
den2 = y (3, :)

! Assign density !
DO j = 1, length_step_r_2
	DO k = 1, length_step_z_2
		IF(rad2(j,k) > rtemp(length_morestep)) THEN
			rho2(j,k) = rho2_min
		ELSE	
			DO i = 1, length_morestep
				IF(rtemp(i) == rad2(j,k)) THEN
					rho2(j,k) = den2(i)
					EXIT
				ELSEIF(rtemp(i) > rad2(j,k)) THEN
					CALL AKIMA(rtemp(i-3), rtemp(i-2), rtemp(i-1), rtemp(i), rtemp(i+1), rtemp(i+2), & 
					den2(i-3), den2(i-2), den2(i-1), den2(i), den2(i+1), den2(i+2), rad2(j,k), rho2(j,k))
					EXIT
				END IF
			END DO
		END IF
	END DO
END DO

! Check the density !
DO j = 1, length_step_r_2
	DO k = 1, length_step_z_2
		IF(rho2(j,k) < rho2_a) THEN
			rho2(j,k) = rho2_a
		END IF
	END DO
END DO

! Boundaries !
call boundary1D_NMFULL (rho2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
m_cell = 0.0D0

! Now we assign the value of density to the original array !
DO j = 1, length_step_r_2
	DO k = 1, j
		If (rho2(k, 1) > rho2_a) THEN
			m_cell (j) = m_cell(j) + rho2(k, 1)*(4.0D0/3.0D0)*pi*(rF2(k)**3 - rF2(k-1)**3)
		ELSE
			CYCLE
		END IF
	END DO
ENDDO

! Assign mass centered coordinate by interpolation !
DO j = 1, length_step_r_2
	m_r2(j) = 0.5D0*(m_cell (j) + m_cell (j-1))
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now of course, we need boundary condition for composition !
! We first assign the appropriate value to the arrays !
! Only if user wants a variable composition of the star !
! Now of course, we need boundary condition for composition !
! We first assign the appropriate value to the arrays !
IF (helmeos_flag == 1) THEN
	DO k = 1, length_step_z_2, 1
		DO j = 1, length_step_r_2, 1
			if (rho2(j,k) > 0.0D0) then	
            			abar2(j,k) = abar_ini
            			zbar2(j,k) = zbar_ini
            			ye2(j,k) = ye2_ini !zbar_ini / abar_ini
            			xiso(j,k,:) = xiso_ini(:)
      	    			temp2(j,k) = temp_ini!
			END IF
		END DO
	END DO

	! Copy the arrays to the ghost shell !
	CALL BOUNDARY2D_X (xiso)
	CALL BOUNDARY1D_NMFULL (abar2, even)       
	CALL BOUNDARY1D_NMFULL (zbar2, even)
	CALL BOUNDARY1D_NMFULL (temp2, even) 
	CALL BOUNDARY1D_NMFULL (ye2, even)  
END IF

END SUBROUTINE