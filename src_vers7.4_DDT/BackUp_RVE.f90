!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutines backup arrays such as density, pressure, etc...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BACKUP_RVE(n)
USE DEFINITION
USE HELMEOS_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: n

! Integer !
integer :: j, k, m

! Open files  !
OPEN(unit= 980, file='Backup_data.dat', action='write')

! Output header !
write(980,*) n
write(980,*) global_time

! Do the loop !
do k = 1, length_step_z_2, 1
   do j = 1, length_step_r_2, 1

      write(980, 701) r2(j), z2(k), rho1(j,k), rho2(j,k), vel2_r(j,k), vel2_z(j,k), phi2(j,k), &
		      temp2(j,k), (xiso(j,k,m), m=1, totalion), &
		      scaG(j,k), flame_ratio_old(j,k), turb_q(j,k)

   enddo
enddo

! Close file !
close(980)

! Format !
701 FORMAT (200ES18.10)
702 FORMAT (F33.15, 204(9E13.5))

end subroutine