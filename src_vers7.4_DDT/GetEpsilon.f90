!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific energy density epsilon      !
! which is necessary for hydro dynamics evolution. We assume     !
! completely degenerate fermi gas EOS for DM and either          !
! finite temperature EOS or completely degenerate EOS for NM     !
! If you want to use your own EOS, you need to take care of this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETEPSILON
USE HELMEOS_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer variable !
INTEGER :: j, k

! Dummy variable !
REAL (DP) :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
! This part is done only if the users wants DM component
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) THEN
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1
				CALL FERMIMO (dlfmmo1, rho1 (j,k), 1)          
				IF (dlfmmo1<=1.0E-2_DP) THEN
					epsilon1 (j,k) = a_max1*small_energy(dlfmmo1)/rho1(j,k)
				ELSE
					epsilon1 (j,k) = a_max1*large_energy(dlfmmo1)/rho1(j,k)
				END IF
			END DO
		END DO
	ELSEIF (polyeosdm_flag == 1) THEN
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1
				epsilon1 (j,k) = k_1 * rho1(j,k) ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
			END DO
		END DO
	END IF

	! We assign the atmospheric specific energy density !
	epsilon1_a = epsilon1(length_step_r_1, length_step_z_1)

	! Ghost shell !
	CALL BOUNDARY1D_DMFULL (epsilon1, even)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the NM case !
! We see whether finite temperature EOS is used !
IF(helmeos_flag == 1) THEN
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
			CALL HELMEOS_RtoE(rho2 (j,k), temp2(j,k), abar2(j,k), zbar2(j,k), ye2(j,k), epsilon2 (j,k), dummy)
		END DO
	END DO
ELSEIF (fermieosnm_flag == 1) THEN
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
			! We need to get the dimensionless fermi momentum for further calculation !
			CALL FERMIMO (dlfmmo2, rho2 (j,k), 2)      
			IF (dlfmmo2<=1.0E-2_DP) THEN
				epsilon2 (j,k) = a_max2*small_energy(dlfmmo2)/rho2(j,k)
			ELSE
				epsilon2 (j,k) = a_max2*large_energy(dlfmmo2)/rho2(j,k)
			END IF
		END DO
	END DO
ELSEIF (polyeosnm_flag == 1) THEN
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
			epsilon2 (j,k) = k_2 * rho2(j,k) ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
		END DO
	END DO
END IF

! We assign the atmospheric specific energy density !
epsilon2_a = epsilon2(length_step_r_2, length_step_z_2)

! Ghost shell !
CALL BOUNDARY1D_NMFULL(epsilon2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the internal energy for cold EOS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDEPSILON
USE DEFINITION
IMPLICIT NONE

! Integer variable !
INTEGER :: j, k

! Dummy variable !
REAL (DP) :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) THEN
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1

				! We need to get the dimensionless fermi momentum for further calculation !
				CALL FERMIMO (dlfmmo1, rho1 (j,k), 1)      
				IF (dlfmmo1<=1.0E-2_DP) THEN
					epsilon1 (j,k) = a_max1*small_energy(dlfmmo1)/rho1(j,k)
				ELSE
					epsilon1 (j,k) = a_max1*large_energy(dlfmmo1)/rho1(j,k)
				END IF
				
			END DO
		END DO
	ELSE
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1
				epsilon1 (j,k) = k_1 * rho1(j,k) ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
			END DO
		END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the NM case !
IF (NM_epsilon == 0) THEN
	IF (fermieosnm_flag == 1) THEN
		DO k = 1, length_step_z_2
			DO j = 1, length_step_r_2

				! We need to get the dimensionless fermi momentum for further calculation !
				CALL FERMIMO (dlfmmo2, rho2 (j,k), 2)       
				IF (dlfmmo2<=1.0E-2_DP) THEN
					epsilon2 (j,k) = a_max2*small_energy(dlfmmo2)/rho2(j,k)
				ELSE
					epsilon2 (j,k) = a_max2*large_energy(dlfmmo2)/rho2(j,k)
				END IF

			END DO
		END DO
	ELSE
		! The following steps are more or less the same !
		DO k = 1, length_step_z_2
			DO j = 1, length_step_r_2
				epsilon2 (j,k) = k_2 * rho2(j,k) ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
			END DO
		END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update the rhoe for dual energy formalism !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDRHOE
USE DEFINITION 
IMPLICIT NONE

! Integer !
Integer :: j, k

DO k = 1, length_step_z_2
	DO j = 1, length_step_r_2
		rhoe2(j,k) = rho2(j,k)*epsilon2(j,k)
	END DO
END DO

CALL BOUNDARY1D_NMFULL (rhoe2,even)

END SUBROUTINE
