!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure profile and their deriative    !
! once the density profile is being updated through rungekutta time  !
! evolution. It is being used in every time step, do not confused it !
! with subroutine GETRHOEOSRTOP                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE RIEMANN_MODULE
USE HELMEOS_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: j, k

! Dummy variables !
REAL (DP) :: dummy, dxdrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
! This would be done only if users wants DM component !
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) Then
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1
				CALL FERMIMO (dlfmmo1, rho1 (j,k), 1)
				CALL FINDDXDRHO (dxdrho, rho1 (j,k), 1)
				IF (dlfmmo1 <= 1.0E-2_DP) THEN
					p1 (j,k) = a_max1*small_pressure(dlfmmo1)
				ELSE
					p1 (j,k) = a_max1*large_pressure(dlfmmo1)
				END IF

				! Derivative !
				dpdrho1 (j,k) = a_max1*dxdrho*dpdx(dlfmmo1)
				dpdepsilon1 (j,k) = 0.0E0_DP
			END DO
		END DO
	ELSEIF (polyeosdm_flag == 1) THEN
		DO k = 1, length_step_z_1
			DO j = 1, length_step_r_1
				p1 (j,k) = k_1 * rho1(j,k) ** gamma1
				dpdrho1 (j,k) = k_1 * gamma1 * rho1(j,k) ** (gamma1 - 1.0D0)
				dpdepsilon1 (j,k) = 0.0D0
			END DO	
		END DO
	END IF
	CALL BOUNDARY1D_DMFULL (p1, even)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The following steps are more or less similar , so no repeat !
IF (helmeos_flag /= 1) THEN
	IF (fermieosnm_flag == 1) THEN
		DO k = 1, length_step_z_2
			DO j = 1, length_step_r_2
				CALL FERMIMO (dlfmmo2, rho2 (j,k), 2)
				CALL FINDDXDRHO (dxdrho, rho2 (j,k), 2)
				IF (dlfmmo2 <= 1.0E-2_DP) THEN
					p2 (j,k) = a_max2*small_pressure(dlfmmo2)
				ELSE
					p2 (j,k) = a_max2*large_pressure(dlfmmo2)
				END IF

				! Derivative !
				dpdrho2 (j,k) = a_max2*dxdrho*dpdx(dlfmmo2)
			END DO
		END DO
	ELSEIF (polyeosnm_flag == 1) THEN	
		IF(nm_epsilon == 0) THEN
			DO k = 1, length_step_z_2
				DO j = 1, length_step_r_2
					p2 (j,k) = k_2 * rho2(j,k) ** gamma2
					dpdrho2 (j,k) = k_2 * gamma2 * rho2(j,k) ** (gamma2 - 1.0D0)
					dpdepsilon2 (j,k) = 0.0D0
				END DO	
			END DO
		ELSE
			DO k = 1, length_step_z_2
				DO j = 1, length_step_r_2
					p2 (j,k) = rho2(j,k) * epsilon2(j,k) * (gamma2 - 1.0E0_DP) 
					dpdrho2 (j,k) = epsilon2(j,k) * (gamma2 - 1.0E0_DP)
					dpdepsilon2 (j,k) = rho2(j,k) * (gamma2 - 1.0E0_DP)
				END DO
			END DO
		END IF
	END IF
ELSE
	DO k = 1, length_step_z_2
		DO j = 1, length_step_r_2
			IF(rho2(j,k) /= rho2_a) THEN
				CALL HELMEOS_RtoP(rho2 (j,k), temp2(j,k), abar2(j,k), zbar2(j,k), ye2(j,k), p2 (j,k), dpdrho2 (j,k), dpdepsilon2 (j,k))
			ELSE
				p2(j,k) = p2_a
				dpdrho2 (j,k) = dpdrho2_a
				dpdepsilon2 (j,k) = dpdepsilon2_a
			END IF
		END DO
	END DO
END IF

! Copy to boundary !
CALL BOUNDARY1D_NMFULL (p2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the soundspeed of DM or NM !                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SOUNDSPEED
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: j, k

! For DM !
IF (DM_flag == 1) THEN
	DO k = 1, length_step_z_1
		DO j = 1, length_step_r_1
			cs1(j,k) = sqrt(dpdrho1(j,k)+dpdepsilon1(j,k)*p1(j,k)/rho1(j,k)**(2.0E0_DP))
		END DO
	END DO
	CALL BOUNDARY1D_DMFULL (cs1, even)
ENDIF

! For NM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!If(helmeos_flag == 1) THEN
!	DO k = 1, length_step_z_2
!		DO j = 1, length_step_r_2
!			CALL HELMEOS_CS (j, k, cs2(j,k))
!		END DO
!	END DO
!ELSE
!END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO k = 1, length_step_z_2
	DO j = 1, length_step_r_2
		cs2(j,k) = sqrt(dpdrho2(j,k)+dpdepsilon2(j,k)*p2(j,k)/rho2(j,k)**(2.0E0_DP))
	END DO
END DO
CALL BOUNDARY1D_NMFULL (cs2, even)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the pressure gradient !                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDPDR
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j, k

! Loop !
DO k = length_step_z_min_part_2, length_step_z_part_2
	DO j = 1, length_step_r_part_2
		dp2dr (j,k) = 0.5D0*(p2(j+1,k) - p2(j-1,k))/dx2
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the pressure gradient !                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDPDZ
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j, k

! Loop !
DO k = length_step_z_min_part_2, length_step_z_part_2
	DO j = 1, length_step_r_part_2
		dp2dz (j,k) = 0.5D0*(p2(j,k+1) - p2(j,k-1))/dx2
	END DO
END DO

END SUBROUTINE
