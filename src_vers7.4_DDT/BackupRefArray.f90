!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine backup all the necessary variables
! which are referred during the Runge-Kutta steps
! Written by Leung Shing Chi in 2016
! If you want your own variables to be backup-ed, 
! simply add your variables in the module and 
! add a command here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BACKUP_REF_ARRAY()
USE DEFINITION
USE HELMEOS_MODULE
IMPLICIT NONE

! Useful hydro variables
rho2_old = rho2
epsilon2_old = epsilon2

! Finite temperature EOS !
IF(helmeos_flag == 1) THEN
	temp2_old = temp2
END IF

! Useful chemical composition
if(helmeos_flag == 1) then
   xiso_old = xiso
endif

END SUBROUTINE BACKUP_REF_ARRAY