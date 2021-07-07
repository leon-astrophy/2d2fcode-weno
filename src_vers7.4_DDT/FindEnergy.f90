!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the total energy
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE findenergy
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k

! Component of energy
real (DP) :: kin1_e, int1_e, pot1_e
real (DP) :: kin2_e, int2_e, pot2_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization 1
energy1 = 0.0D0
energy2 = 0.0D0

! Initialization 2
energy1_kin = 0.0D0
energy1_int = 0.0D0
energy1_pot = 0.0D0

! Initialization 3
energy2_kin = 0.0D0
energy2_int = 0.0D0
energy2_pot = 0.0D0

! Initialization 4
kin1_e = 0.0D0
int1_e = 0.0D0
pot1_e = 0.0D0

! Initialization 5
kin2_e = 0.0D0
int2_e = 0.0D0
pot2_e = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the dark matter total energy
if(DM_flag == 1) then

   DO k = length_step_z_min_part_1, length_step_z_part_1
      DO j = 1, length_step_r_part_1 
	
	 IF(rho1(j,k) > rho1_a) THEN  	   
   	   kin1_e = kin1_e + vol1(j,k) * 0.5D0 * rho1(j,k) * (vel1_r(j,k)**2 + vel1_z(j,k)**2 + vel1_p(j,k))
           int1_e = int1_e + vol1(j,k) * rho1(j,k) * epsilon1(j,k)
           pot1_e = pot1_e + vol1(j,k) * 0.5D0 * (rho1(j,k) * phi1(j,k))
   	 ENDIF
              
      ENDDO
   ENDDO
	
endif

! Sum all the components to get the total energy	
energy1_kin = kin1_e
energy1_pot = pot1_e
energy1_int = int1_e
energy1 = energy1_kin + energy1_pot + energy1_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the normal matter total energy
        
DO k = length_step_z_min_part_2, length_step_z_part_2, 1
   DO j = 1, length_step_r_part_2, 1
	      
       IF(rho2(j,k) > rho2_a) THEN
          kin2_e = kin2_e + vol2(j,k) * 0.5D0 * rho2(j,k) * (vel2_r(j,k)**2 + vel2_z(j,k)**2 + vel2_p(j,k))
          int2_e = int2_e + vol2(j,k) * rho2(j,k) * epsilon2(j,k)
          pot2_e = pot2_e + vol2(j,k) * 0.5D0 * (rho2(j,k) * phi2(j,k))
       ENDIF
	      
   ENDDO
ENDDO

! Sum all the components to get the total energy	
energy2_kin = kin2_e
energy2_pot = pot2_e
energy2_int = int2_e
energy2 = energy2_kin + energy2_pot + energy2_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE findenergy