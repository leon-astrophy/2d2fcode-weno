!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine assigns the initial velocity
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GETVEL
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k

! Assign DM hydro velocity
if(runDM_flag == 1) then

   vel1_r = 0.0D0
   vel1_z = 0.0D0
   vel1_p = 0.0D0

   CALL BOUNDARY1D_DMFULL (vel1_r, oddR)
   CALL BOUNDARY1D_DMFULL (vel1_z, oddZ)
   CALL BOUNDARY1D_DMFULL (vel1_p, oddR)  

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign NM hydro velocity
vel2_r = 0.0D0 
vel2_z = 0.0D0
vel2_p = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use this if homologus expansion/contraction is assumed
!do k = 1, length_step_z, 1
!   do j = 1, length_step_r, 1
!      vel2_r(j,k) = -1.0D-6 * r2(j)
!      vel2_z(j,k) = -1.0D-6 * z2(k)
!   enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use this if rotation is needed
do k = 1, length_step_z_2, 1
   do j = 1, length_step_r_2, 1
      vel2_p(j,k) = r2(j) * omega_ini
   enddo
enddo

! Set up the velocity for the ghost cells
CALL BOUNDARY1D_NMFULL (vel2_r, oddR)
CALL BOUNDARY1D_NMFULL (vel2_z, oddZ)
CALL BOUNDARY1D_NMFULL (vel2_p, oddR)

END SUBROUTINE