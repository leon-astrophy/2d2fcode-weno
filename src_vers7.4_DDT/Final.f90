!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine close all the allocated
! arrays which is opened at the beginning 
! of the code
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINAL
USE DEFINITION
USE HELMEOS_MODULE
!USE JETEXPLOSION_MODULE
USE PPT_MODULE
USE TURB_MODULE
USE LEVELSET_MODULE
IMPLICIT NONE

! Destroy hydro variables !
CALL DESTROYHYDRO

! Deallocate !
IF (fermieosnm_flag == 1) THEN 
	DEALLOCATE(eostable2)
END IF
IF(DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) THEN
		DEALLOCATE(eostable1)
	END IF
END IF

if(tracer_flag == 1) then
   call destroyPPT
endif

if(flame_flag == 1) then
   call destroyLevelset
endif 

if(turb_flag == 1) then
   call destroyTurb
endif

if(helmeos_flag == 1) then
   call destroyhelm
endif

!if(jetexp_flag == 1) then
!   call destroyjetexp
!endif

! For DM Rungekutta !
DEALLOCATE(l1)
DEALLOCATE(u2_dm)
DEALLOCATE(u3_dm)
DEALLOCATE(l3_dm)
DEALLOCATE(u_old1)
DEALLOCATE(u_new1)
DEALLOCATE(u_temp1)

! For NM Rungekutta !
DEALLOCATE(l2)
DEALLOCATE(u2_nm)
DEALLOCATE(u3_nm)
DEALLOCATE(l3_nm)
DEALLOCATE(u_old2)
DEALLOCATE(u_new2)
DEALLOCATE(u_temp2)

END SUBROUTINE FINAL