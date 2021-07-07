!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for the primitive variables
! Written by Leung Shing Chi in 2016
! The subroutine takes the full U arrays as input
! and do the odd/even parity extension
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, switch the checkstep_flag = 0 
! For hybrid boundaries, such as the quadrant star
! Specific modifications are needed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE boundary2D_DM (u)
USE DEFINITION
USE helmeos_module
USE Levelset_module
USE Turb_module
IMPLICIT NONE

! Input conservative flux array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5,imin1 : imax1) :: u

! Dummy variables
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'In Boundary 2D'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the first (r-inner) boundary

IF(boundary_flag(1) == 0) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(1-j,:,k) = u(length_step_r_part_1+1-j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(1-j,:,k) = bfac_r(k) * u(j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(1) == 2) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(1-j,:,k) = u(1,:,k)
      ENDDO
   ENDDO

ENDIF

! Do the second (r-outer) boundary

IF(boundary_flag(2) == 0) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(length_step_r_part_1+j,:,k) = u(j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(length_step_r_part_1+j,:,k) = bfac_r(k) * u(length_step_r_part_1+1-j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(length_step_r_part_1+j,:,k) = u(length_step_r_part_1,:,k)                  
      ENDDO
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(:,length_step_z_min_part_1-j,k) = u(:,length_step_z_part_1+1-j,k)                     
      ENDDO
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(:,length_step_z_min_part_1-j,k) = bfac_z(k) * u(:,length_step_z_min_part_1-1+j,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1              
         u(:,length_step_z_min_part_1-j,k) = u(:,length_step_z_min_part_1,k)
      ENDDO             
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(:,length_step_z_part_1+j,k) = u(:,length_step_z_min_part_1-1+j,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(:,length_step_z_part_1+j,k) = bfac_z(k) * u(:,length_step_z_part_1+1-j,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO k = imin1, imax1, 1
      DO j = 1, 5, 1
         u(:,length_step_z_part_1+j,k) = u(:,length_step_z_part_1,k)
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for Normal matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE boundary2D_NM (u)
USE DEFINITION
USE helmeos_module
USE Levelset_module
USE Turb_module
IMPLICIT NONE

! Input conservative flux array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5,imin2 : imax2) :: u

! Dummy variables
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'In Boundary 2D'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the first (r-inner) boundary

IF(boundary_flag(1) == 0) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(1-j,:,k) = u(length_step_r_part_2+1-j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(1-j,:,k) = bfac_r(k) * u(j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(1) == 2) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(1-j,:,k) = u(1,:,k)
      ENDDO
   ENDDO

ENDIF

! Do the second (r-outer) boundary

IF(boundary_flag(2) == 0) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(length_step_r_part_2+j,:,k) = u(j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(length_step_r_part_2+j,:,k) = bfac_r(k) * u(length_step_r_part_2+1-j,:,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(length_step_r_part_2+j,:,k) = u(length_step_r_part_2,:,k)                  
      ENDDO
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(:,length_step_z_min_part_2-j,k) = u(:,length_step_z_part_2+1-j,k)                     
      ENDDO
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(:,length_step_z_min_part_2-j,k) = bfac_z(k) * u(:,length_step_z_min_part_2-1+j,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1              
         u(:,length_step_z_min_part_2-j,k) = u(:,length_step_z_min_part_2,k)
      ENDDO             
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(:,length_step_z_part_2+j,k) = u(:,length_step_z_min_part_2-1+j,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(:,length_step_z_part_2+j,k) = bfac_z(k) * u(:,length_step_z_part_2+1-j,k)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO k = imin2, imax2, 1
      DO j = 1, 5, 1
         u(:,length_step_z_part_2+j,k) = u(:,length_step_z_part_2,k)
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE