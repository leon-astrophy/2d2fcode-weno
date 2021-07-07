!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
! Welcome to SNIa_Phase IV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This code is a multi-dimensional hydro code which
! aims at modeling Type-Ia supernovae during its
! explosion phase up to homologus expansion
! This is built based on many physics input,
! therefore, the code is designed to be flexible
! to add/delete physics component to fit your need.
!
! The code uses WENO (Weighted essentially non-oscillating)
! fifth-order scheme for spatial discretization and 
! third-order five-step non-strong stability preseving 
! scheme for time discretization. 
!
! Furthermore, the code has applied level-set method
! for front capturing, flame/detonation physics, 
! sub-grid scale turbulence, chemical composition
! tracer scheme and magnetohydrodynamics. They 
! can be turned on/off to fit the purpose of the 
! simulation. 
! 
! Written by Leung Shing Chi in 2016   
! Questions and comments towards the code are welcome.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM STAR_WENO
USE DEFINITION
USE PPT_MODULE
USE HELMEOS_MODULE
USE NUSPEC_MODULE
USE HELMEOS_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k, m

! For atmoshperic density !
REAL (DP) :: dummy, rhoaold

! Initial and local step number
INTEGER :: n, n_ini, n_ini2

! File number 
integer :: file_no, fileno_ini, fileno_ini2

! Signal for starting the run
INTEGER :: flag_run

! Signal for continuing the run
integer :: flag_continue_run
	
! The file number for different output
integer :: filecount_log
integer :: filecount_profile
integer :: filecount_hydro
integer :: filecount_flame
integer :: filecount_helm
integer :: filecount_turb
integer :: filecount_PPT

! Package for file name
INTEGER :: fileno_len
CHARACTER (len = 256) :: fileno

! Flag for finishing simulations
INTEGER :: exit_flag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initial time
global_time = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*)
WRITE(*,*) 'Welcome to WENO version 7.1'
WRITE(*,*) 

! Report all initial parameter setting
CALL ReportParameter

! Do all the initial stuff here
CALL INITIAL (n_backup, fileno_ini, fileno_ini2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now deal with the file I/O matter
! Initial parameters for simulations
filecount_log = 0
filecount_profile = fileno_ini + 1
filecount_hydro = fileno_ini + 1
filecount_flame = fileno_ini + 1
filecount_helm = fileno_ini + 1
filecount_turb = fileno_ini + 1
filecount_PPT = fileno_ini + 1

! Other flag !
flag_error = 0
flag_run = 0

! Initialize !
output_file = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare output for initial profiles
write(*,*) 'Open file for initial storage'
WRITE (fileno, *) filecount_log
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)

! Output !
write(*,*) 'Output initial profile'

WRITE (fileno, *) filecount_profile
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)
CALL OPENFILE_PROFILE(fileno, fileno_len)
CALL OUTPUT_PROFILE_FULL (n_backup)
CALL CLOSEFILE_PROFILE

WRITE (fileno, *) filecount_log
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)
CALL OPENFILE_LOG(fileno, fileno_len)
CALL OUTPUT_LOG (n_backup)

! Output others 
IF (flame_flag == 1) THEN
	call outputlevelset(filecount_flame)
END IF
IF (tracer_flag == 1) THEN
	call outputPPT(filecount_PPT)
END IF
IF (turb_flag == 1) THEN
	call outputturb_profile(filecount_turb)
END IF
IF (helmeos_flag == 1) THEN
	call outputXiso_profile(filecount_helm)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main trunk for the evolution

! Coniditon !
write(*,*) 'Finished preparing inital data, ready for run? (1 for yes, 0 for no,)' 
WRITE(*,*)
!read(*,*) flag_run
flag_run = 1

!!!!!!!!!!!!!!!!!!!!!!!!!
!if(flag_run == 1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!

DO n = 1, time_step

	! Evolve one-step
	CALL RUNGEKUTTA (n)

	! Update the time
	global_time = global_time + dt
	WRITE(*,*) 'iteration time', n
	WRITE(*,*) 'global time', global_time
	IF (RUNDM_flag == 1) THEN
		WRITE(*,*) 'length_step_r_part_1', length_step_r_part_1
  		WRITE(*,*) 'length_step_z_part_1', length_step_z_part_1
		WRITE(*,*) 'DM ATM Density r', rho1(length_step_r_1, 1), rho1_a
		WRITE(*,*) 'DM ATM Density z', rho1(1, length_step_z_1), rho1_a
		WRITE(*,*) 'DM grid size', dx1
	END IF
        WRITE(*,*) 'length_step_r_part_2', length_step_r_part_2
        WRITE(*,*) 'length_step_z_part_2', length_step_z_part_2
	WRITE(*,*) 'NM ATM Density r', rho2(length_step_r_2, 1), rho2_a
	WRITE(*,*) 'NM ATM Density z', rho2(1, length_step_z_2), rho2_a
	WRITE(*,*) 'NM grid size', dx2
        WRITE(*,*) 'time step', dt
	WRITE(*,*) 'cfl', cfl
	WRITE(*,*)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Check if we need to use expanding grid method

	! When DDT is assumed, first do for NM !
	IF(found_movinggridnm_flag == 1 .and. movinggridnm_flag == 0) then
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! When PTD is assumed
 	!IF(found_movinggrid_flag == 1 .and. movinggrid_flag == 0) THEN
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   	  If(length_step_r_2 - length_step_r_part_2 <= 10 .OR. length_step_z_2 - length_step_z_part_2 <= 10) then
	      WRITE (*,*) 'Movinggrid triggered for NM'
	      WRITE (*,*) 'At global time = ', global_time
	      WRITE (*,*)
	      movinggridnm_flag = 1
              checkstepnm_flag = 0
	      length_step_r_part_2 = length_step_r_2
	      length_step_z_part_2 = length_step_z_2
	      output_file = .true.
              !EXIT
	   END IF
        END IF

	! Now do for DM !
	IF(RUNDM_flag == 1) THEN
           IF(found_movinggriddm_flag == 1 .and. movinggriddm_flag == 0) then
	      If(length_step_r_1 - length_step_r_part_1 <= 10 .OR. length_step_z_1 - length_step_z_part_1 <= 10) then
	         WRITE (*,*) 'Movinggrid triggered for DM'
	         WRITE (*,*) 'At global time = ', global_time
	         WRITE (*,*)
	         movinggriddm_flag = 1
                 checkstepdm_flag = 0
	         length_step_r_part_1 = length_step_r_1
	         length_step_z_part_1 = length_step_z_1
                 output_file = .true.
 	      END IF
	   END IF
	END IF
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	! Output log 
	if(ABS(global_time - output_logtime_last) >= output_logtime .or. output_file .eqv. .true.) then
	   output_logtime_last = global_time
	   WRITE (fileno, *) filecount_log
	   fileno = ADJUSTL (fileno)             
	   fileno_len = LEN_TRIM (fileno)                     
	   CALL OUTPUT_LOG (n_backup + n)         
        endif
	
	! Output profile
	if(output_flag == 1) then
   	   if(ABS(global_time - output_profiletime_last) >= output_profiletime .or. output_file .eqv. .true.) then
	      output_profiletime_last = global_time
              filecount_profile = filecount_profile + 1
	      filecount_hydro = filecount_hydro + 1

	      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      !Prepare the filename suffix
	      WRITE (fileno, *) filecount_profile
	      fileno = ADJUSTL (fileno)
	      fileno_len = LEN_TRIM (fileno)
	      CALL OPENFILE_PROFILE(fileno, fileno_len)

	      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      !Output to files
              CALL OUTPUT_PROFILE_FULL (n_backup + n)
	      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      !Close profile files
	      CALL CLOSEFILE_PROFILE
           endif

	   ! Output turbulence
	   IF(turb_flag == 1) THEN
               IF(ABS(global_time - output_turbtime_last) >= output_turbtime .or. output_file .eqv. .true.) then
                 filecount_turb = filecount_turb + 1
                 output_turbtime_last = global_time
                 CALL outputTurb_profile(filecount_turb)
               ENDIF
           ENDIF

	   ! Output flame
 	   IF(flame_flag == 1 .and. fusion_flag == 1) THEN
              IF(ABS(global_time - output_flametime_last) >= output_flametime .or. output_file .eqv. .true.) THEN
                 filecount_flame = filecount_flame + 1
	         output_flametime_last = global_time
                 CALL outputLevelSet(filecount_flame)
              ENDIF
	   ENDIF

	   ! Output chemical composition
   	   IF(helmeos_flag == 1) THEN
              IF(ABS(global_time - output_Helmtime_last) >= output_Helmtime .or. output_file .eqv. .true.) then
	         filecount_helm = filecount_helm + 1
                 output_Helmtime_last = global_time
                 call outputXiso_profile(filecount_helm)
              ENDIF
           ENDIF

	   ! Output tracer particle scheme
   	   IF(tracer_flag == 1) THEN
	      IF(ABS(global_time - output_PPTtime_last) >= output_PPTtime .or. output_file .eqv. .true.) then
	         filecount_PPT = filecount_PPT + 1
                 output_PPTtime_last = global_time
	         CALL outputPPT(filecount_PPT)
	      ENDIF
	   ENDIF

	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! Do the backup	 
	   ! Note: I muted it because the file size is too large to handle
	   ! Besides, the old file reading option can serve somehow 
	   ! similar function.
   	   !if(mod(n,500) == 0) then
 	   !   call backup_rve(n)
	   !endif
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	endif
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Quit if any error arise
	!IF(flag_error == 1) THEN
        !    EXIT
        !ENDIF   
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Check if we need to exit the simulation
        IF(global_time >= total_time) THEN
		WRITE (*,*) 'Reaching the end of simulations...'
 	    	WRITE (*,*)
		EXIT
	END IF
	IF(initmodel_flag /= 0 .AND. dt < 1.0D-30) THEN
		WRITE (*,*) 'Time step too small...'
		WRITE (*,*)
		EXIT
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!if(dx > 15.0D0) exit
	!if(length_step_r_part >= 0.95D0 * length_step_r) exit
        !if(length_step_z_part >= 0.95D0 * length_step_z) exit
        !IF(Exit_flag == 1) EXIT
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Switch off burning after some time to speed up calculation
	IF(fusion_flag == 1) THEN
   	   IF(deton_flag == 0 .and. global_time > 3.0D5) THEN
	      fusion_flag = 0
	      turb_flag = 0
	      WRITE (*,*) 'Burning ends at time =', global_time
	   ELSEIF(deton_flag == 1 .and. global_time > 4.0D5) THEN
	      fusion_flag = 0
              turb_flag = 0
	      WRITE (*,*) 'Burning ends at time =', global_time
	   ENDIF
	ENDIF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Update time step
	IF(Updatedt_flag == 1) THEN
		CALL FindDt
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Use full box if the star has expanded
	!IF(length_step_r_part == length_step_r .and. length_step_z_part == length_step_z) checkstep_flag = 0
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! If found_deton_flag == 1 but after some time, relax print profile frequency
	IF(found_deton_flag == 1) THEN
	   IF(global_time - found_deton_time > 1.0D5) THEN
              output_profiletime = 2.0D4
              output_flametime = 2.0D4       
              output_turbtime = 2.0D4                 
              output_Helmtime = 2.0D4       
	   ENDIF
        ENDIF

	! Use less frequent output after some time for SNIa
	IF(fusion_flag == 0) THEN
	   IF(global_time <= 4.0D5 .and. global_time + dt > 4.0D5) THEN
	      output_logtime = 2.0D3
	      output_profiletime = 1.0D5
              output_flametime = 1.0D5
              output_turbtime = 1.0D5
              output_Helmtime = 1.0D5
	      output_PPTtime = 2.0D4
	   ENDIF
	ENDIF

	! For CCSN, relax to usual profile output after the jet has finished pumping energy
	!IF(jetexp_flag == 1) THEN
	!   IF(global_time > 2.0D0 * jet_time) THEN
	!      output_profiletime = 2.0D4
        !      output_flametime = 2.0D4
        !      output_turbtime = 2.0D4
        !      output_Helmtime = 2.0D4
	!   ENDIF
  	!ENDIF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Automatically switch off output_file flag
	output_file = .false.

END DO

!!!!!!!
!endif
!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Final output for the end of simulations
WRITE (fileno, *) filecount_profile + 1
filecount_hydro = filecount_hydro + 1
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)
CALL OPENFILE_PROFILE(fileno, fileno_len)
CALL OUTPUT_PROFILE_FULL (n_backup + n)
CALL CLOSEFILE_PROFILE

WRITE (fileno, *) filecount_log
fileno = ADJUSTL (fileno)             
fileno_len = LEN_TRIM (fileno)                       
CALL OUTPUT_LOG (n_backup + n)         

! For others
if(turb_flag == 1) then
   filecount_turb = filecount_turb + 1
   call outputturb_profile(filecount_turb)
endif
if(flame_flag == 1) then
   filecount_flame = filecount_flame + 1
   call outputLevelSet(filecount_flame)
endif
if(tracer_flag == 1) then
   filecount_PPT = filecount_PPT + 1
   call outputPPT(filecount_PPT)
endif
if(helmeos_flag == 1) then
   filecount_helm = filecount_helm + 1
   CALL outputXiso_profile (filecount_helm)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output parameter files

! Open 
OPEN (UNIT = 999, FILE = './Outfile//Star_WENO_Parameter.dat', STATUS = 'REPLACE')

! Output !
WRITE (999,*) 'rhoc1', rhomax1
WRITE (999,*) 'rhoc2', rhomax2
WRITE (999,*) 'rhoa1', rho1_a
WRITE (999,*) 'rhoa2', rho2_a
WRITE (999,*) 'filecount_hydro', filecount_hydro
WRITE (999,*) 'filecount_flame', filecount_flame
WRITE (999,*) 'filecount_helm', filecount_helm
WRITE (999,*) 'filecount_turb', filecount_turb
WRITE (999,*) 'filecount_PPT', filecount_PPT

! Close !
CLOSE (999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare to switch off

! Deallocate all previously allocated variables
CALL FINAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the backup
! Note: I muted it again due to the same 
!       problem described above
!call backup_rve(n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the backup
! Note: I muted it because this is inherited
!       from the old version, which is not 
!       sufficient for the current version
!IF(flag_run == 1) THEN
!   CALL BACKUP (n_backup + n - 1)
!   CALL BACKUP_RVE
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print out message 
WRITE (*,*) 'Simulation End!'

! Format !
100 FORMAT (20ES18.8)

END PROGRAM
