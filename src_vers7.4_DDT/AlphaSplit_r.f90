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

SUBROUTINE ALPHASPLIT_r (alpha1, alpha2, k, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type

! Input: The row number
INTEGER, INTENT(IN) :: k

! Output: The effective speed
REAL (DP), INTENT (OUT) :: alpha1(imin1:imax1)
REAL (DP), INTENT (OUT) :: alpha2(imin2:imax2)

! dummy variables
INTEGER :: i, j

! The candidate effective speeds for non-MHD
real (DP) :: lambda(3)

! The maximum effective speeds
real (DP) :: lambda_max(3)

! Initialization
lambda_max(:) = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For Debug 
!IF(debug_flag == 1) WRITE(*,*) 'In AlphaSplit_r'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Include DM component
IF(type == 1) THEN
   DO j = 1, length_step_r_part_1, 1
      lambda(1) = ABS(vel1_r(j,k)) + DSQRT (cs1(j, k)) !dpdrho1 (j,k) + (p1 (j,k) / rho1 (j,k) ** 2) * dpdepsilon1(j,k)
      lambda_max(1) = MAX(lambda(1), lambda_max(1))
   END DO

   ! Send the DM output
   alpha1(imin1:imax1) = lambda_max(1)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Include NM component
IF(type == 2) THEN 
   DO j = 1, length_step_r_part_2, 1           
      lambda(2) = ABS(vel2_r(j,k)) + DSQRT (cs2(j, k)) !dpdrho2 (j,k) + (p2 (j,k) / rho2 (j,k) ** 2) * dpdepsilon2(j,k)
      lambda_max(2) = MAX(lambda(2), lambda_max(2))
   ENDDO

   ! Send the NM output 
   alpha2(imin2:imax2) = lambda_max(2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF(SRhydro_flag == 0) THEN
!ELSEIF(SRhydro_flag == 1) THEN
!   DO j = 1, length_step_r_part_2, 1
!      lambda(2) = ABS(vel2_r(j,k)) + DSQRT (dpdrho2 (j,k) * rho2(j,k) / (rho2(j,k) + rho2(j,k) * epsilon2(j,k) + p2(j,k)) + &
!                      (p2 (j,k) / rho2 (j,k) / (rho2(j,k) + rho2(j,k) * epsilon2(j,k) + p2(j,k))) * dpdepsilon2(j,k))
!      lambda_max(2) = MAX(lambda(2), lambda_max(2))
!   ENDDO
!ELSE
!   STOP 'Check SRhydro_flag (1 = with SR, 0 = no SR)'
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For Debug 
!IF(debug_flag == 1) WRITE(*,*) 'Done alphasplit_r'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Alphasplit_r