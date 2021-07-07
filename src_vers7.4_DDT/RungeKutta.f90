!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine does one single Runge-Kutta full step
! It uses the opeator splitting and separate
! all non-gravitational source term to be done 
! after the hydro step.
!
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RUNGEKUTTA (n)
USE DEFINITION
!USE JETEXPLOSION_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
USE HELMEOS_MODULE
USE PPT_MODULE
IMPLICIT NONE

! The input step number
INTEGER, INTENT (IN) :: n

! Dummy variables
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dummy U array due to Runge-Kutta scheme
!REAL (DP), allocatable, DIMENSION (:,:,:) :: u_temp, u_one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy dx due to Runge Kutta scheme
REAL (DP) :: dx1_old, dx1_two, dx1_three, dxdt1_three
REAL (DP) :: dx2_old, dx2_two, dx2_three, dxdt2_three

! Dummy !
REAL (DP) :: rhoaold, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check timing with or without openmp
!INTEGER :: time_start, time_end
!INTEGER :: time_start2, time_end2
!INTEGER :: cr, cm
!REAL :: rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count time !
!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)
!WRITE(*,*) "system_clock rate ",rate          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First allocate dummy U 
! CALL system_clock(time_start)
!ALLOCATE(u_temp(-4:length_step_r_max+5, -4:length_step_z_max+5,1:no_of_eq))
!ALLOCATE(u_one(-4:length_step_r_max+5, -4:length_step_z_max+5,1:no_of_eq))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update the first half of particle tracer
IF(tracer_flag == 1) THEN
	call evolve_p_1st
END IF

! Find the boundaries !
IF(movinggridnm_flag == 1) THEN
	CALL FINDRADIUS_NM
END IF
IF(movinggriddm_flag == 1) THEN
	CALL FINDRADIUS_DM
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) 'Preparation = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st iteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Movinggrid !
IF(movinggridnm_flag == 1) THEN
	dx2_old = dx2
END IF
IF(movinggriddm_flag == 1) THEN
	dx1_old = dx1
END IF

CALL BACKUPCONS (u_old1, u_temp1, u_old2, u_temp2)
CALL SPATIAL (u_old1, u_old2)

IF (RUNDM_flag == 1) THEN
 DO i = imin1, imax1, 1
   DO k = length_step_z_min_part_1, length_step_z_part_1, 1
      DO j = 1, length_step_r_part_1, 1
         u_new1 (j,k,i) = u_old1 (j,k,i) + 0.391752226571890D0 * dt * l1 (j,k,i)
      ENDDO
   ENDDO
 END DO
 CALL BOUNDARY2D_DM (u_new1)
END IF

DO i = imin2, imax2, 1
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1
         u_new2 (j,k,i) = u_old2 (j,k,i) + 0.391752226571890D0 * dt * l2 (j,k,i)
      ENDDO
   ENDDO
END DO

! Copy the data to ghost cells in U
CALL BOUNDARY2D_NM (u_new2)

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = dx1_old + 0.391752226571890D0 * dt * vel1_max * dx1 / radius1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = dx2_old + 0.391752226571890D0 * dt * vel2_max * dx2 / radius2
   call getgridnm
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL SYSTEM_CLOCK(time_start2)
!CALL SYSTEM_CLOCK(time_end2)           
!WRITE(*,*) '1st FINDAZ =', REAL(time_end2 - time_start2) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert from conservative to primitive
CALL FROMUTORVE (u_new1, u_new2)

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
IF(burn_prog_flag == 1) THEN
	CALL Find_AZBar_fromBPhi
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update the physical quantities
!CALL SYSTEM_CLOCK(time_start2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update 
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU (u_new1, u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL SYSTEM_CLOCK(time_end2)
!WRITE(*,*) '1st update =', REAL(time_end2 - time_start2) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'Done first partial step'
!CALL REPORT_CURRENT_DATA(1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) '1st substep ', REAL(time_end - time_start)/rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd iteration

IF(tracer_flag == 1) THEN
	call evolve_p_2nd
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL BACKUPCONS (u_old1, u_temp1, u_old2, u_temp2)
CALL SPATIAL (u_new1, u_new2)

IF (RUNDM_flag == 1) THEN
 DO i = imin1, imax1, 1
   DO k = length_step_z_min_part_1, length_step_z_part_1
      DO j = 1, length_step_r_part_1
	u_new1 (j, k, i) = 0.444370493651235D0 * u_old1 (j, k, i) + 0.555629506348765D0 * u_new1 (j, k, i) + 0.368410593050371D0 * dt * l1 (j, k, i)
	u2_dm (j, k, i) = u_new1 (j, k, i)
      ENDDO
   END DO
 END DO
 CALL BOUNDARY2D_DM (u_new1)
END IF

DO i = imin2, imax2, 1
   DO k = length_step_z_min_part_2, length_step_z_part_2
      DO j = 1, length_step_r_part_2
	u_new2 (j, k, i) = 0.444370493651235D0 * u_old2 (j, k, i) + 0.555629506348765D0 * u_new2 (j, k, i) + 0.368410593050371D0 * dt * l2 (j, k, i)
	u2_nm (j, k, i) = u_new2 (j, k, i)
      ENDDO
   END DO
END DO

! Copy the data to the ghost cells in U
CALL BOUNDARY2D_NM (u_new2)

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = 0.444370493651235D0 * dx1_old + 0.555629506348765D0 * dx1 + 0.368410593050371D0 * dt * vel1_max * dx1 / radius1
   dx1_two = dx1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = 0.444370493651235D0 * dx2_old + 0.555629506348765D0 * dx2 + 0.368410593050371D0 * dt * vel2_max * dx2 / radius2
   dx2_two = dx2
   call getgridnm
ENDIF

! Convert from conservative to primitive
CALL FROMUTORVE (u_new1,u_new2)

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
IF(burn_prog_flag == 1) THEN
	CALL FIND_AZBAR_fromBPhi()
END IF

! Update physical quantities
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU (u_new1, u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'Done second partial step'
!CALL REPORT_CURRENT_DATA(1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) '2nd substep = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3rd iteration

IF(tracer_flag == 1) THEN
	call evolve_p_3rd
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)
IF (RUNDM_flag == 1) THEN
 DO i = imin1, imax1, 1
   DO k = length_step_z_min_part_1, length_step_z_part_1
      DO j = 1, length_step_r_part_1        
	u_new1 (j, k, i) = 0.620101851488403D0 * u_old1 (j, k, i) + 0.379898148511597D0 * u_new1 (j, k, i) + 0.251891774271694D0 * dt * l1 (j, k, i)
	u3_dm (j, k, i) = u_new1 (j, k, i)
	l3_dm (j, k, i) = l1 (j, k, i)	
      ENDDO
   END DO
 END DO
 CALL BOUNDARY2D_DM (u_new1)
END IF

DO i = imin2, imax2, 1
   DO k = length_step_z_min_part_2, length_step_z_part_2
      DO j = 1, length_step_r_part_2        
	u_new2 (j, k, i) = 0.620101851488403D0 * u_old2 (j, k, i) + 0.379898148511597D0 * u_new2 (j, k, i) + 0.251891774271694D0 * dt * l2 (j, k, i)
	u3_nm (j, k, i) = u_new2 (j, k, i)
	l3_nm (j, k, i) = l2 (j, k, i)
      ENDDO
   END DO
END DO

! Update dx
IF(movinggriddm_flag == 1) then
   dxdt1_three = vel1_max * dx1 / radius1
   dx1 = 0.620101851488403D0 * dx1_old + 0.379898148511597D0 * dx1 + 0.251891774271694D0 * dt * vel1_max * dx1 / radius1
   dx1_three = dx1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dxdt2_three = vel2_max * dx2 / radius2
   dx2 = 0.620101851488403D0 * dx2_old + 0.379898148511597D0 * dx2 + 0.251891774271694D0 * dt * vel2_max * dx2 / radius2
   dx2_three = dx2
   call getgridnm
ENDIF

! Copy the data to the ghost cells in U
CALL BOUNDARY2D_NM (u_new2)

! Convert from conservative to primitive
CALL FROMUTORVE (u_new1, u_new2)

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
IF(burn_prog_flag == 1) THEN
	CALL FIND_AZBAR_fromBPhi()
END IF

! Update physical quantities
CALL UPDATE (0)

! Do conversion again !
CALL FROMRVETOU (u_new1, u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff 
!if(debug_flag == 1) write(*,*) 'Done third partial step' 
!CALL REPORT_CURRENT_DATA(1,1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) '3rd substep = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fourth iteration

IF(tracer_flag == 1) THEN
	call evolve_p_4th
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)
IF (RUNDM_flag == 1) THEN
 DO i = imin1, imax1, 1
   DO k = length_step_z_min_part_1, length_step_z_part_1
      DO j = 1, length_step_r_part_1
	u_new1 (j, k, i) = 0.178079954393132D0 * u_old1 (j, k, i)  + 0.821920045606868D0 * u_new1 (j, k, i)  + 0.544974750228521D0 * dt * l1 (j, k, i) 
      ENDDO
   END DO
 END DO
 CALL BOUNDARY2D_DM (u_new1)
END IF

DO i = imin2, imax2, 1
   DO k = length_step_z_min_part_2, length_step_z_part_2
      DO j = 1, length_step_r_part_2
	u_new2 (j, k, i) = 0.178079954393132D0 * u_old2 (j, k, i) + 0.821920045606868D0 * u_new2 (j, k, i) + 0.544974750228521D0 * dt * l2 (j, k, i)
      ENDDO
   END DO
END DO

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = 0.178079954393132D0 * dx1_old + 0.821920045606868D0 * dx1 + 0.544974750228521D0 * dt * vel1_max * dx1 / radius1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = 0.178079954393132D0 * dx2_old + 0.821920045606868D0 * dx2 + 0.544974750228521D0 * dt * vel2_max * dx2 / radius2
   call getgridnm
ENDIF

! Copy the data to ghost cells in U
CALL BOUNDARY2D_NM (u_new2)

! Convert from conservative to primitive
CALL FROMUTORVE (u_new1, u_new2)

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
IF(burn_prog_flag == 1) THEN
	CALL FIND_AZBAR_fromBPhi()
END IF

! Update physical quantities
CALL UPDATE (0)

! Convert again !
CALL FROMRVETOU (u_new1,u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1)	write(*,*) 'Done fourth partial step'
!	CALL REPORT_CURRENT_DATA(1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) '4th substep = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for next step

IF(tracer_flag == 1) THEN
	call evolve_p_5th
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)
IF (RUNDM_flag == 1) THEN
 DO i = imin1, imax1, 1
   DO k = length_step_z_min_part_1, length_step_z_part_1
      DO j = 1, length_step_r_part_1
	u_new1 (j, k, i) = 0.517231671970585D0 * u2_dm (j, k, i) + 0.096059710526147D0 * u3_dm (j, k, i) + 0.386708617503269D0 * u_new1 (j, k, i) &
		+ 0.063692468666290D0 * dt * l3_dm (j, k, i) + 0.226007483236906D0 * dt * l1 (j, k, i)
      ENDDO
   END DO
 END DO
 CALL BOUNDARY2D_DM (u_new1)
END IF

DO i = imin2, imax2, 1
   DO k = length_step_z_min_part_2, length_step_z_part_2
      DO j = 1, length_step_r_part_2
	u_new2 (j, k, i) = 0.517231671970585D0 * u2_nm (j, k, i) + 0.096059710526147D0 * u3_nm (j, k, i) + 0.386708617503269D0 * u_new2 (j, k, i) &
		+ 0.063692468666290D0 * dt * l3_nm (j, k, i) + 0.226007483236906D0 * dt * l2 (j, k, i)
      ENDDO
   END DO
END DO

! Update dx
IF(movinggriddm_flag == 1) then
   dx1 = 0.517231671970585D0 * dx1_two + 0.096059710526147D0 * dx1_three + 0.386708617503269D0 * dx1 &
		+ 0.063692468666290D0 * dt * dxdt1_three + 0.226007483236906D0 * dt * vel1_max * dx1 / radius1
   call getgriddm
ENDIF
IF(movinggridnm_flag == 1) then
   dx2 = 0.517231671970585D0 * dx2_two + 0.096059710526147D0 * dx2_three + 0.386708617503269D0 * dx2 &
		+ 0.063692468666290D0 * dt * dxdt2_three + 0.226007483236906D0 * dt * vel2_max * dx2 / radius2
   call getgridnm
ENDIF

! Copy the data to ghost cells in U
CALL BOUNDARY2D_NM (u_new2)

! Convert from conservative to primitive
CALL FROMUTORVE (u_new1, u_new2)

! Update Abar and Zbar
IF(xisotran_flag == 1) THEN 
	CALL FIND_AZBAR()
END IF
IF(burn_prog_flag == 1) THEN
	CALL FIND_AZBAR_fromBPhi()
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Adjust the sponge
IF (sponge_flag == 1) THEN
	call getsponge
END IF

! Update temperature
IF (helmeos_flag == 1) THEN
	CALL FindhelmTEMP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(nuceos_flag == 1) call findnuctemp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'Prepare for next step (before fusion)'
!CALL REPORT_CURRENT_DATA(1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) '5th substep = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for burning

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Type Ia supernovae 
IF(fusion_flag == 1) THEN

   ! If there is level-set, update it
   IF(flame_flag == 1) CALL UPDATE_FLAME_RADIUS 

   ! This trigger the burning package proposed by
   ! Reinecke 1999b
   IF(xisotran_flag == 1) THEN

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The change of level-set needs to energy input
      !IF(flame_flag == 1) CALL FLAME
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This does the Carbon burning
      IF(carburn_flag == 1) CALL BURN_PHASE1B

      ! This do the O- and Si- burning
      IF(advburn_flag == 1) CALL BURN_PHASE2B

      ! Update the AZbar and temperature accordingly
      CALL FIND_AZBAR
      CALL FINDHELMTEMP

      ! For completely burnt zone, check if NSE applies
      IF(convert_nse_flag == 1) CALL NSE2

      ! Copy the new Xiso and epsilon to ghost cells
      CALL BOUNDARY2D_X()
      CALL BOUNDARY1D_NM(epsilon2, even)

      ! Check if the change of isotope perserve the sum
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CALL system_clock(time_start2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL CHECKXISOTOPE

      ! Update the burntimw
      last_burntime = global_time

      ! Update Abar and Zbar and temperature again
      CALL Find_AZBAR()
      CALL FindhelmTemp

   ENDIF

   IF(burn_prog_flag == 1) THEN

      ! Call the carbon burning scheme
      CALL burn_phase1

      ! Call the NQSE burning scheme 
      CALL burn_phase2

      ! Call the NSE burning scheme
      CALL NSE3

      ! Get the isotope abundance from the 
      ! bunring-progress variables
      CALL reconstructXIso

      ! Update Abar and Zbar
      CALL FIND_AZBAR_fromBPhi()

      ! Find the new temperature
      CALL FINDhelmTEMP

   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Special subroutine for setting detonation seed
   !call finddetonseed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Plant a delayed He-det seed
   !IF(global_time < 10.0D4 .and. global_time + dt >= 10.0D4) THEN
   !   CALL plantdetonseed         
   !ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For AIC, the electron capture occurs as 
! a function of density
if(ecap_flag == 1) then 
   call findEcap
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!CALL report_current_data(1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!IF(jetexp_flag == 1) THEN
!   CALL findesum
!   IF(global_time <= jet_time .or. (global_time > jet_time .and. erad_sum > 1.0D-10)) then
!      CALL cleaninnerzone 
!   ENDIF
!   IF(global_time < 2.0D0 * jet_time) THEN
!      CALL solveJetDiff
!      CALL solveJetDep
!   ENDIF
!ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!CALL report_current_data(1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fix the atmosphere 
! For AICwDM, the drastic change of typical
! DM density makes the adjustment of definition of
! DM atmosphere become necessary
!if(runDM_flag == 1) then
!   rho1_a = MIN(0.01D0 * rho1(1,1), 5.0D-7)
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Obsolete subroutine for fixing the atmosphere
!call fix_rho2_a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For typical AIC or SNIa, the rapid change
! of density also makes the redefinition
! of atmosphere necessary
!rho2_a = MIN(rho2_a, 1.0D-4 * rho2(1,1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Notice we need to find the thermodyanmics properties 
! of the newly defined atmosphere state
!CALL teos_epsilon(rho2_a, temp_a, abar_ini, zbar_ini, ye_ini1, epsilon2_a)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some artificial limiter of the velocity
!IF(SRHydro_flag == 0 .and. checkvel_flag == 1) THEN
!   Call fix_vel
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Section for adjusting atmospheric density !
IF(fixrhonm_flag == 1) THEN
	rhoaold = rho2_a
	rho2_a = min(min(maxval(rho2), rho2_c)*rhofac_2, rho2_a)

	! Assign atmoshperic epsilon !
	IF(nm_epsilon == 1 .AND. polyeosnm_flag == 1) THEN
		epsilon2_a = k_2 * rho2_a ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
	ELSEIF(helmeos_flag == 1) THEN
		CALL HELMEOS_RtoE (rho2_a, temp_a, abar_ini, zbar_ini, ye2_ini, epsilon2_a, dummy)
	ELSE
		CALL EOSEPSILON (rho2_a, dummy, epsilon2_a, 2)
	END IF

	! Adjust density !
	DO j = -4, length_step_r_2 + 5
		DO k = -4, length_step_z_2 + 5
			IF(rho2(j,k) == rhoaold) THEN
				rho2(j,k) = rho2_a
			END IF
		END DO
	END DO

	! Pressure !
	IF(helmeos_flag == 1) THEN
		CALL HELMEOS_RtoP(rho2_a, temp_a, abar_ini, zbar_ini, ye2_ini, p2_a, dpdrho2_a, dpdepsilon2_a)
	END IF
END IF

! For DM !
IF(fixrhodm_flag == 1) THEN
	rhoaold = rho1_a
	rho1_a = min(min(maxval(rho1), rho1_c)*rhofac_1, rho1_a)

	! Assign atmoshperic epsilon !
	CALL EOSEPSILON (rho1_a, dummy, epsilon1_a, 1)

	! Adjust density !
	DO j = -4, length_step_r_1 + 5
		DO k = -4, length_step_z_1 + 5
			IF(rho1(j,k) == rhoaold) THEN
				rho1(j,k) = rho1_a
			END IF
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check density !
IF (checkrho_flag == 1) THEN
	CALL CHECKRHO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count ending time 
!CALL system_clock(time_end)
!WRITE(*,*) 'Burning = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read to start a new step. Do some preparation
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert from primitive to conservative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end2)
!WRITE(*,*) 'FromRVEToU = ', REAL(time_end2 - time_start2) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update physical quantities
IF (MOD (n, 2) == 0) THEN
	CALL UPDATE (1)
ELSE
	CALL UPDATE (0)
END IF

! Update again !
CALL FROMRVETOU (u_new1, u_new2)
CALL BACKUPCONS (u_new1, u_old1, u_new2, u_old2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug 
!CALL REPORT_CURRENT_DATA(1, 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare data for the output purpose
IF (tracer_flag == 1) THEN
	call evolve_p_final
END IF

! Backup everything 
CALL BACKUP_REF_ARRAY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!
! Deallocate arrays !
!deallocate(u_temp)
!deallocate(u_one)
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) 'Final = ', REAL(time_end - time_start) / rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Format !
100 format (20ES15.7)

END SUBROUTINE