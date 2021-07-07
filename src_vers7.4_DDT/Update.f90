!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine prepare the data necessary for constructing
! the flux for the spatial discretization.
! It takes input/output of the U array and the 
! mode p which decides whether or not to update
! the gravitational potential
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE UPDATE (p)
USE DEFINITION
USE LEVELSET_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Signal for finding potential
INTEGER, INTENT (IN) :: p

! Dummy variable
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check timing with or without openmp
!INTEGER :: time_start, time_end
!INTEGER :: cr, cm
!REAL :: rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Counting time !
!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'In Update'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) 'Copy u to v = ', REAL(time_end - time_start)/rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find temperature !
IF(helmeos_flag == 1) THEN
	CALL FINDHELMTEMP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(nuceos_flag == 1) call findnuctemp
!CALL system_clock(time_end)
!WRITE(*,*) 'FindHelmtemp = ', REAL(time_end - time_start)/rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find pressure
CALL FINDPRESSURE

! For dual energy
IF(dual_energy == 1) THEN

	! Pressure graident !
	CALL FINDDPDR
	CALL FINDDPDZ
	
	! Find internal energy rhoe !
	CALL FINDRHOE

END IF

! Speed of sound !
CALL SOUNDSPEED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)  
!WRITE(*,*) 'FindPressure = ', REAL(time_end - time_start)/rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update turbulence transport term
IF(turb_flag == 1) THEN
	CALL FINDTURBULENCE
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)  
!WRITE(*,*) 'FindTurb = ', REAL(time_end - time_start)/rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update the potential
IF (p == 1 .OR. (testmodel_flag == 5 .AND. initmodel_flag == 0)) THEN
	CALL FINDPOTENTIAL
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALL system_clock(time_end)
!WRITE(*,*) 'FindPotential = ', REAL(time_end - time_start)/rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine backup the conservative variables for the !
! next temporal evolution				    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BACKUPCONS (v1, u1, v2, u2)
USE DEFINITION
USE WENO_MODULE
IMPLICIT NONE

! The new U-array
REAL (DP), INTENT (IN), DIMENSION (-4:length_step_r_1 + 5 , -4 : length_step_z_1 + 5,imin1:imax1) :: v1
REAL (DP), INTENT (IN), DIMENSION (-4:length_step_r_2 + 5 , -4 : length_step_z_2 + 5,imin2:imax2) :: v2

! The old U array
REAL (DP), INTENT (OUT), DIMENSION (-4:length_step_r_1 + 5 , -4 : length_step_z_1 + 5,imin1:imax1) :: u1
REAL (DP), INTENT (OUT), DIMENSION (-4:length_step_r_2 + 5 , -4 : length_step_z_2 + 5,imin2:imax2) :: u2

! Integer parameter !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! These are the neccesarry variables that need to be updated each time !
! Assign the value of v !
! Copy the new data to old data
DO i = imin2, imax2, 1
   DO k = -4, length_step_z_2 + 5
      DO j = -4, length_step_r_2 + 5
	 u2 (j,k,i) = v2 (j,k,i)
      ENDDO
   END DO
END DO

! For DM 
IF(RUNDM_flag == 1) THEN
 DO i = imin1, imax1, 1
   DO k = -4, length_step_z_1 + 5
      DO j = -4, length_step_r_1 + 5
	 u1 (j,k,i) = v1 (j,k,i)
      ENDDO
   END DO
 END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE