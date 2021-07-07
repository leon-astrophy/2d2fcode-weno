!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! This subroutine calculates the maximum time step
! which satisfies the Courant condition 
! Written by Leung Shing Chi in 2016  
! 
! If you modify the Euler equation, make sure you change this 
! part to include the new effective sound speed
!
! Limiters are posed based on output time and running time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE finddt
USE definition
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k

! Dummy
REAL (DP) :: dummy

! Local effective speed
REAL (DP) :: lambda2a, lambda2b

! Local maximum effective speed
REAL (DP) :: lambda1, lambda2, lambda3

! Local minimum dt for DM, NM and 1st overlayer
REAL (DP) :: dt_temp1, dt_temp2, dt_temp3

! Local timestep
REAL (DP) :: dt_loc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF(debug_flag == 1) WRITE(*,*) 'In FindDt'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize by setting an arbitrarily large number
dt_temp1 = 1.0D5
dt_temp2 = 1.0D5
dt_temp3 = 1.0D5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we find the minimum time constrained by DM sector
! For non-MHD case
! For dark matter
IF(runDM_flag == 1) THEN
   DO k = length_step_z_min_part_1, length_step_z_part_1, 1
      DO j = 1, length_step_r_part_1, 1

	 ! Only grid with density above threshold density is counted
	 IF(rho1(j,k) > rho1_a) THEN
	    lambda2a = ABS(vel1_r(j,k)) + cs1(j,k)
            lambda2b = ABS(vel1_z(j,k)) + cs1(j,k)
            lambda1 = MAX(lambda2a, lambda2b)
	    dt_temp1 = MIN(cfl * dx1 / lambda1, dt_temp1)
	 ENDIF

      ENDDO
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we find the minimum time constrained by NM sector
! For non-MHD case
! For normal matter
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         ! Only grid with density above threshold density is counted
         IF(rho2(j,k) > rho2_a) THEN

            lambda2a = ABS(vel2_r(j,k)) + cs2(j,k)
            lambda2b = ABS(vel2_z(j,k)) + cs2(j,k)
            lambda2 = MAX(lambda2a, lambda2b)
            dt_temp2 = MIN(dt_temp2, cfl * dx2 / lambda2)

         ENDIF

      ENDDO
   ENDDO

! Only the minimum one is chosen
dt_loc = MIN(dt_temp1, dt_temp2, dt_temp3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!dt = MIN(dt_temp_min, 1.1D0 * dt)
!WRITE(*,*) 'AT time:', global_time
!WRITE(*,*) 'Primitive dt', dt_loc, lambda2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Here we further limit dt by considering output time and running time
! If we reach some of the output time, then lower dt accordingly

!IF(global_time + dt_loc - output_logtime_last >= output_logtime) THEN
!   dt_loc = MIN(output_logtime - global_time + output_logtime_last, dt_loc)
!ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Special flag for CCSN because of the jet heating
!IF(jetexp_flag == 1 .and. global_time < 2.0D0 * jet_time) THEN
!   dt_loc = MIN(10.0D0, dt_loc)
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(output_flag == 1) THEN

   !IF(global_time + dt_loc - output_profiletime_last > output_profiletime) THEN
   !   dt_loc = MIN(output_profiletime - (global_time - output_profiletime_last), dt_loc)
   !   !output_file = .true.
   !ENDIF

   ! For helmeos
   !IF(helmeos_flag == 1) THEN
   !   IF(global_time + dt_loc - output_helmtime_last > output_helmtime) THEN
   !      dt_loc = MIN(output_helmtime - (global_time - output_helmtime_last), dt_loc)
   !	 output_file = .true.
   !   ENDIF
   !ENDIF 

   ! For SGS turbulence
   !IF(turb_flag == 1) THEN
   !   IF(global_time + dt_loc - output_turbtime_last > output_turbtime) THEN 
   !      dt_loc = MIN(output_turbtime - (global_time - output_turbtime_last), dt_loc)
   !	 output_file = .true.
   !   ENDIF
   !ENDIF

   ! For level-set 
   !IF(fusion_flag == 1) THEN
   !   IF(flame_flag == 1) THEN
   !      IF(global_time + dt_loc - output_flametime_last > output_flametime) THEN
   !         dt_loc = MIN(output_flametime - (global_time - output_flametime_last), dt_loc)
   !	    output_file = .true.
   !      ENDIF
   !   ENDIF
   !ENDIF

   ! For PPT
   !IF(tracer_flag == 1) THEN
   !  IF(global_time + dt_loc - output_PPTtime_last > output_PPTtime) THEN 
   !     dt_loc = MIN(output_PPTtime - (global_time - output_PPTtime_last), dt_loc)
   !	output_file = .true.
   !   ENDIF
   !ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check if we reach the maximum time
!IF(global_time + dt_loc > total_time) THEN
!   dt_loc = MIN(total_time - global_time, dt_loc)
!ENDIF

! Set a lower limit to avoid deadlock
dt = MAX(1.0D-10, dt_loc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff
!WRITE(*,*) 'Final dt', dt
!WRITE(*,*)
!IF(Debug_flag == 1) WRITE(*,*) 'Finished FindDt'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE FindDt