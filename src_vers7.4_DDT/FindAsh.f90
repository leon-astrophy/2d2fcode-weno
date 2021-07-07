!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 
! This subroutine calculates the mass of matter 
! which is burnt by deflagration 
! rates and assumes NSE composition 
! Written by Leung Shing Chi in 2016 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE findash
USE definition
USE levelset_module
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k

! Initialization
mass_ash = 0.0D0

DO k = length_step_z_min_part_2, length_step_z_part_2, 1
   DO j = 1, length_step_r_part_2, 1

      ! Just sum the mass up when the grid is burnt
      ! no matter partially or completely
      mass_ash = mass_ash + vol2(j,k) * rho2 (j, k) * flame_ratio(j,k) 

   ENDDO
ENDDO

END SUBROUTINE findash