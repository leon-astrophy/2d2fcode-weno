!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file stores subroutines for opening files. They include
! 1. openfile_profile
! 2. openfile_profile2 (for more frequent output, such as flame surface)
! 3. openfile_log
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPENFILE_PROFILE(FILENO, FILENO_LEN)
USE DEFINITION
IMPLICIT NONE

! Related to file index !
INTEGER :: fileno_len
CHARACTER (len = 256) :: fileno

! Assign file index !
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM 
If (RUNDM_flag == 1) THEN
	OPEN (UNIT = 101, FILE = './Outfile/Hydro/Star_WENO_Potential_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 102, FILE = './Outfile/Hydro/Star_WENO_Density_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 103, FILE = './Outfile/Hydro/Star_WENO_Epsilon_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 104, FILE = './Outfile/Hydro/Star_WENO_Pressure_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
 	OPEN (UNIT = 105, FILE = './Outfile/Hydro/Star_WENO_Velocity_DM_r_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 106, FILE = './Outfile/Hydro/Star_WENO_Velocity_DM_z_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	IF(rotationdm_flag == 1) THEN
		OPEN (UNIT = 107, FILE = './Outfile/Hydro/Star_WENO_Velocity_DM_p_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM
OPEN (UNIT = 200, FILE = './Outfile/Hydro/Star_WENO_Potential_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 201, FILE = './Outfile/Hydro/Star_WENO_Density_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 202, FILE = './Outfile/Hydro/Star_WENO_Velocity_NM_r_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 203, FILE = './Outfile/Hydro/Star_WENO_Velocity_NM_z_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 204, FILE = './Outfile/Hydro/Star_WENO_Epsilon_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 205, FILE = './Outfile/Hydro/Star_WENO_Pressure_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
IF(rotationdm_flag == 1) THEN
	OPEN (UNIT = 206, FILE = './Outfile/Hydro/Star_WENO_Velocity_NM_p_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF
IF(helmeos_flag == 1) THEN
	OPEN (UNIT = 207, FILE = './Outfile/Hydro/Star_WENO_Temperature_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')   
END IF
IF(etran_flag == 1) THEN
	OPEN (UNIT = 208, FILE = './Outfile/Hydro/Star_WENO_Ye_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! For files that are needed to be opened once such as integrated global quantity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

SUBROUTINE OPENFILE_LOG(FILENO, FILENO_LEN)
USE DEFINITION
IMPLICIT NONE

! File index !
INTEGER :: fileno_len, i
CHARACTER (len = 256) :: fileno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Common 
OPEN (UNIT = 11, FILE ='./Outfile/Hydro/Star_WENO_Energy_Total_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM 
If (DM_flag == 1) THEN
	OPEN (UNIT = 12, FILE ='./Outfile/Hydro/Star_WENO_Mass_DM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
END IF
If (RUNDM_flag == 1) THEN
	OPEN (UNIT = 13, FILE ='./Outfile/Hydro/Star_WENO_Energy_DM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
	OPEN (UNIT = 14, FILE ='./Outfile/Hydro/Star_WENO_CentralDensity_DM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
END IF	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM 
OPEN (UNIT = 21, FILE ='./Outfile/Hydro/Star_WENO_Mass_NM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
OPEN (UNIT = 22, FILE ='./Outfile/Hydro/Star_WENO_Energy_NM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
OPEN (UNIT = 23, FILE ='./Outfile/Hydro/Star_WENO_CentralDensity_NM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
OPEN (UNIT = 24, FILE ='./Outfile/Hydro/Star_WENO_Radius_NM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
IF(helmeos_flag == 1) THEN
	OPEN (UNIT = 25, FILE ='./Outfile/Hydro/Star_WENO_CentralTemperature_NM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
END IF
IF(etran_flag == 1) THEN
	OPEN (UNIT = 26, FILE ='./Outfile/Hydro/Star_WENO_CentralYe_NM_'//fileno (1 : fileno_len)//'.dat', ACTION = 'WRITE', POSITION = 'APPEND')
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine