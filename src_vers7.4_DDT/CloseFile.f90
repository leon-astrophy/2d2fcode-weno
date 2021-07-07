!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine close all the PROFILE files which 
! are used to record the hydro variables
! Written by Leung Shing Chi in 2016
! You will need to change the file number when 
! more/fewer files are included
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR CLOSING PROFILE FILES

SUBROUTINE CLOSEFILE_PROFILE
USE DEFINITION
IMPLICIT NONE

! Integer !
integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DM profile !
If (RUNDM_flag == 1) THEN
  DO i = 101, 106, 1
    CLOSE (i)
  END DO
  IF(rotationdm_flag == 1) THEN
    CLOSE (107)
  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NM profile !
DO i = 200, 205, 1
   CLOSE (i)
END DO
IF(rotationnm_flag == 1) THEN
   CLOSE (206)
END IF
IF(helmeos_flag == 1) THEN
   CLOSE (207)
END IF
IF(etran_flag == 1) THEN
   CLOSE (208)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CLOSEFILE_PROFILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          
! This subroutine close all the LOG files which
! are used to record the hydro variables
! Written by Leung Shing Chi in 2016
! You will need to change the file number when 
! more/fewer files are included
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For closing log files

SUBROUTINE CLOSEFILE_LOG
USE DEFINITION
IMPLICIT NONE

! Integer 
integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Common log
CLOSE (11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DM log
If (DM_flag == 1) THEN
   CLOSE (12)
END IF
If (RUNDM_flag == 1) THEN
   CLOSE (13)
   CLOSE (14)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NM log
DO i = 21, 24, 1
   CLOSE (i)
END DO
IF(helmeos_flag == 1) THEN
   CLOSE (25)
END IF
IF(etran_flag == 1) THEN
   CLOSE (26)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CLOSEFILE_LOG