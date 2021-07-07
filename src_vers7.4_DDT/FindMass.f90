!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 
! This subroutine calculates the total mass
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDMASS
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k

! threshold density
real (DP) :: rho_min

! Initialization
mass1 = 0.0E0_DP
mass2 = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sum the DM mass
! Do it no for DM_flag == 1

If(dm_flag == 1) THEN
 DO k = length_step_z_min_part_1, length_step_z_part_1
   DO j = 1, length_step_r_part_1

      IF(rho1(j,k) > rho1_a) THEN
         mass1 = mass1 + vol1(j,k) * rho1 (j, k)
      ENDIF

   ENDDO
 ENDDO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sum the normal matter mass

DO k = length_step_z_min_part_2, length_step_z_part_2, 1
   DO j = 1, length_step_r_part_2, 1

      IF(rho2(j,k) > rho2_a) THEN
         mass2 = mass2 + vol2(j,k) * rho2 (j, k)
      ENDIF

   ENDDO
ENDDO

END SUBROUTINE FindMass