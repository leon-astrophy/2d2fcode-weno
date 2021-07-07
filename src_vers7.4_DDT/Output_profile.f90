!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file contains two subroutines for printing data in files
! 1. output_profile_full 
! 2. output_flame
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUT_PROFILE_FULL (N)
USE DEFINITION
USE HELMEOS_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Input step number
INTEGER, INTENT (IN) :: n

! dummy variables
INTEGER :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM 
If (RUNDM_flag == 1) THEN

	! Potential
	WRITE (101, *) n, global_time
	WRITE (101, *) length_step_r_1, length_step_z_1
	WRITE (101, *) dx1, dt
	DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
		WRITE (101, 701) (phi1(j,k), j=1, length_step_r_1)
	END DO
	WRITE (101, *)

	! Density
	WRITE (102, *) n, global_time
	WRITE (102, *) length_step_r_1, length_step_z_1
	WRITE (102, *) dx1, dt
	DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
		WRITE (102, 701) (rho1(j,k), j=1, length_step_r_1)
	END DO
	WRITE (102, *)

	! Epsilon
	WRITE (103, *) n, global_time
	WRITE (103, *) length_step_r_1, length_step_z_1
	WRITE (103, *) dx1, dt
	DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
		WRITE (103, 701) (epsilon1(j,k), j=1, length_step_r_1)
	END DO
	WRITE (103, *)

	! Pressure
	WRITE (104, *) n, global_time
	WRITE (104, *) length_step_r_1, length_step_z_1
	WRITE (104, *) dx1, dt
	DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
		WRITE (104, 701) (p1(j,k), j=1, length_step_r_1)
	END DO
	WRITE (104, *)

	! Velocity_r
	WRITE (105, *) n, global_time
	WRITE (105, *) length_step_r_1, length_step_z_1
	WRITE (105, *) dx1, dt
	DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
		WRITE (105, 701) (vel1_r(j,k), j=1, length_step_r_1)
	END DO
	WRITE (105, *)

	! Velocity_z
	WRITE (106, *) n, global_time
	WRITE (106, *) length_step_r_1, length_step_z_1
	WRITE (106, *) dx1, dt
	DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
		WRITE (106, 701) (vel1_z(j,k), j=1, length_step_r_1)
	END DO
	WRITE (106, *)

	! Velocity_p
	If (rotationdm_flag == 1) THEN
		WRITE (107, *) n, global_time
		WRITE (107, *) length_step_r_1, length_step_z_1
		WRITE (107, *) dx1, dt
		DO k = 1, length_step_z_1   ! You don't need from 1 to length_step if DM is at the core
			WRITE (107, 701) (vel1_p(j,k), j=1, length_step_r_1)
		END DO
		WRITE (107, *)
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\

! For NM Potential
WRITE (200, *) n, global_time
WRITE (200, *) length_step_r_2, length_step_z_2
WRITE (200, *) dx2, dt 
DO k = 1, length_step_z_2
	WRITE (200, 701) (phi2(j,k), j=1, length_step_r_2)
END DO
WRITE (200, *)

! For NM Density 
WRITE (201, *) n, global_time
WRITE (201, *) length_step_r_2, length_step_z_2
WRITE (201, *) dx2, dt 
DO k = 1, length_step_z_2
	WRITE (201, 701) (rho2(j,k), j=1, length_step_r_2)
END DO
WRITE (201, *)

! Velocity_r
WRITE (202, *) n, global_time
WRITE (202, *) length_step_r_2, length_step_z_2
WRITE (202, *) dx2, dt 
DO k = 1, length_step_z_2
	WRITE (202, 701) (vel2_r(j,k), j=1, length_step_r_2)
END DO
WRITE (202, *)

! Velocity_z
WRITE (203, *) n, global_time
WRITE (203, *) length_step_r_2, length_step_z_2
WRITE (203, *) dx2, dt 
DO k = 1, length_step_z_2
	WRITE (203, 701) (vel2_z(j,k), j=1, length_step_r_2)
END DO
WRITE (203, *)

! Epsilon 
WRITE (204, *) n, global_time
WRITE (204, *) length_step_r_2, length_step_z_2
WRITE (204, *) dx2, dt 
DO k = 1, length_step_z_2
	WRITE (204, 701) (epsilon2(j,k), j=1, length_step_r_2)
END DO
WRITE (204, *)

! Pressure 
WRITE (205, *) n, global_time
WRITE (205, *) length_step_r_2, length_step_z_2
WRITE (205, *) dx2, dt 
DO k = 1, length_step_z_2
	WRITE (205, 701) (p2(j,k), j=1, length_step_r_2)
END DO
WRITE (205, *)

! Velocity_p
If (rotationnm_flag == 1) THEN
	WRITE (206, *) n, global_time
	WRITE (206, *) length_step_r_2, length_step_z_2
	WRITE (206, *) dx2, dt 
	DO k = 1, length_step_z_2
		WRITE (206, 701) (vel2_p(j,k), j=1, length_step_r_2)
	END DO
	WRITE (206, *)
END IF
     
! Temperature !
If (helmeos_flag == 1) THEN
	WRITE (207, *) n, global_time
	WRITE (207, *) length_step_r_2, length_step_z_2
	WRITE (207, *) dx2, dt 
	DO k = 1, length_step_z_2
		WRITE (207, 701) (temp2(j,k), j=1, length_step_r_2)
	END DO
	WRITE (207, *)
END IF

! Electron fractions
If (etran_flag == 1) THEN
	WRITE (208, *) n, global_time
	WRITE (208, *) length_step_r_2, length_step_z_2
	WRITE (208, *) dx2, dt 
	DO k = 1, length_step_z_2
		WRITE (208, 701) (ye2(j,k), j=1, length_step_r_2)
	END DO
	WRITE (208, *)
END IF

701 FORMAT (2000E14.5E3)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!