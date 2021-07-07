!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the global stellar data
! Written by Leung Shing CHi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUT_LOG (n)
USE DEFINITION
USE HELMEOS_MODULE, ONLY : OUTPUT_XMASS
USE LEVELSET_MODULE
USE TURB_MODULE
USE NUSPEC_MODULE
USE GW_MODULE, ONLY : OUTPUTGRAVQ
IMPLICIT NONE

! input step number
INTEGER, INTENT (IN) :: n

! dummy variable
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the information first
CALL PREPARE_LOGDATA()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output mandatory data

! Common !
WRITE (11, 701) global_time, energy1 + energy2, energy_input

! DM log files 
If (DM_flag == 1) THEN
	WRITE (12, 701) global_time, mass1
END IF
If (RUNDM_flag == 1) THEN
	WRITE (13, 702) global_time, energy1, energy1_kin, energy1_int, energy1_pot
	WRITE (14, 702) global_time, rho1(1,1), maxval(rho1)
END IF

! NM log files 
WRITE (21, 701) global_time, mass2
WRITE (22, 702) global_time, energy2, energy2_kin, energy2_int, energy2_pot	
WRITE (23, 701) global_time, rho2(1,1), MAXVAL(rho2)
WRITE (24, 702) global_time, r2(length_step_r_part_2), z2(length_step_z_part_2)
IF(helmeos_flag == 1) THEN
	WRITE (25, 701) global_time, temp2(1,1), MAXVAL(temp2)
END IF
IF(etran_flag == 1) THEN
	WRITE (26, 701) global_time, ye2(1,1), MINVAL(ye2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output other quantities when needed
if(turb_flag == 1) call outputturb_log
if(nuspec_flag == 1) call output_nuphi
if(gravwave_flag == 1) call outputGravQ
if(helmeos_flag == 1) CALL output_xmass()
!if(nucEOS_flag == 1) call outputNucEOS_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Format !
701 FORMAT (F23.15, 2ES18.10)
702 FORMAT (F23.15, 30ES18.10)
703 FORMAT (F23.15, 30ES18.10)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares all the information which 
! are presented in the log file.
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PREPARE_LOGDATA() 
USE DEFINITION
USE NUSPEC_MODULE, ONLY : FINDNUSPEC
USE GW_MODULE, ONLY : FINDGRAVQ
USE HELMEOS_MODULE, ONLY : FINDCENTRALTEMPERATURE, &
			   FINDNEUTRINOLOSS, &
			   FINDTOTALNEUTRINOLOSS, &
			   FINDLUMINOSITY
USE TURB_MODULE
IMPLICIT NONE

! Find Epsilon !
CALL FINDEPSILON

! Find related to hydro !
CALL FINDMASS
CALL FINDENERGY
CALL FINDCENTRALDENSITY
!!!!!!!!!!!!!!!
!CALL FINDBOX
!!!!!!!!!!!!!!!

! Related to nuclear reaction or finite temperature eos
if(helmeos_flag == 1) CALL FINDNEUTRINOLOSS
if(helmeos_flag == 1) CALL FINDTOTALNEUTRINOLOSS
if(helmeos_flag == 1) CALL FINDLUMINOSITY
if(helmeos_flag == 1) CALL FINDCENTRALTEMPERATURE

! Other stuff !
if(turb_flag == 1) call findturbenergy
if(flame_flag == 1) CALL FINDASH
if(nuspec_flag == 1) call findNuSpec
if(gravwave_flag == 1) call findGravQ

END SUBROUTINE