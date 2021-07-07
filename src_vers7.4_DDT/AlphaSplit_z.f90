!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates the maximum effective speed along a constant r
! Written by Leung Shing Chi in 2016
! The effective speed is obtained by solving the determinant of the 
! Jacobian of the flux vector written in terms of the primitive variables
! Modification of this subroutine is necessary when we modify the 
! Euler equation to include other physics, such as B-field
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ALPHASPLIT_Z (alpha1, alpha2, j, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type

! Input: The column number
INTEGER, INTENT(IN) :: j

! Output: The effective speed
REAL (DP), INTENT (OUT) :: alpha1(imin1:imax1)
REAL (DP), INTENT (OUT) :: alpha2(imin2:imax2)

! dummy variables
INTEGER :: i, k

! The candidate effective speeds for non-MHD
REAL (DP) :: lambda(3)

! The maximum effective speeds
REAL (DP) :: lambda_max(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For Debug !
!IF(debug_flag == 1) WRITE(*,*) 'In alphasplit_z'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization
lambda_max(:) = 0.0D0

! Include DM motion
if(type == 1) then
   DO k = length_step_z_min_part_1, length_step_z_part_1
      lambda(1) = ABS(vel1_z(j,k)) + DSQRT (cs1(j, k)) !dpdrho1 (j,k) + (p1 (j,k) / rho1 (j,k) ** 2) * dpdepsilon1(j,k))
      lambda_max(1) = MAX(lambda(1), lambda_max(1))
   END DO

   ! Send the DM output
   alpha1(imin1:imax1) = lambda_max(1)       
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For non-MHD 
IF(type == 2) THEN 
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      lambda(2) = ABS(vel2_z(j,k)) + DSQRT (cs2(j, k)) !dpdrho2 (j,k) + (p2 (j,k) / rho2 (j,k) ** 2) * dpdepsilon2(j,k)
      lambda_max(2) = MAX(lambda(2), lambda_max(2))
   enddo 

   ! Send the NM output
   alpha2(imin2:imax2) = lambda_max(2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF(SRhydro_flag == 0) THEN
!ELSEIF(SRhydro_flag == 1) THEN
!   DO k = 1, length_step_z_part_2, 1
!      lambda(2) = ABS(vel2_z(j,k)) + DSQRT (dpdrho2 (j,k) * rho2(j,k) / (rho2(j,k) + rho2(j,k) * epsilon2(j,k) + p2(j,k)) + & 
!                     (p2 (j,k) / rho2 (j,k) / (rho2(j,k) + rho2(j,k) * epsilon2(j,k) + p2(j,k))) * dpdepsilon2(j,k)) 
!      lambda_max(2) = MAX(lambda(2), lambda_max(2))
!   ENDDO
!ELSE
!   STOP 'Check SRhydro_flag (1 = with SR, 0 = no SR)'
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For Debug !
!IF(debug_flag == 1) WRITE(*,*) 'Done alphasplit_z'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE