!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! This subroutine calculates the central density
! Written by Leung Shing Chi in 2016  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDCENTRALDENSITY
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd order interpolation by Ka Wing
!centralrho1 = ((1.5E0_DP * dx) ** 2 * rho1 (1) - (5.0E-1_DP * dx) ** 2 * rho1 (2)) &
!		/ ((1.5E0_DP * dx) ** 2 - (5.0E-1_DP * dx) ** 2)
!centralrho2 = ((1.5E0_DP * dx) ** 2 * rho2 (1) - (5.0E-1_DP * dx) ** 2 * rho2 (2)) &
!		/ ((1.5E0_DP * dx) ** 2 - (5.0E-1_DP * dx) ** 2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! My version: What you see is what you get
IF(hemisphere_flag == 1) THEN
	centralrho2 = rho2(1,length_step_z_2/2)
ELSE
	centralrho2 = rho2(1,1)
END IF

! For DM !
IF(DM_flag == 1) THEN
	IF(hemisphere_flag == 1) THEN
		centralrho1 = rho1(1,length_step_z_1/2)
	ELSE
		centralrho1 = rho1(1,1)
	END IF
END IF

END SUBROUTINE