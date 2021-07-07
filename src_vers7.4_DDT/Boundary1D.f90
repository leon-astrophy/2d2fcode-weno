!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for WENO (DM and NM seperately)
! Written by Leung Shing Chi in 2016
! The subroutine takes ARRAY as input/output and SIGN
! for doing odd/even parity extension
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, check Boundary1D_FULL.f90
! For hybrid boundaries, such as the quadrant star 
! Specific modifications are needed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_DM (array, sign)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Input parity
INTEGER, INTENT (IN) :: sign

! Dummy variables
INTEGER :: j, k

! Parity factor
INTEGER :: fac_r, fac_z

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_r = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_r = -1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_r = 1
   fac_z = -1
ELSEIF(sign == 3) THEN
   fac_r = -1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do it case by case
IF(boundary_flag(1) == 0) THEN

   DO j = 1, 5, 1
      array(1-j,:) = array(length_step_r_part_1+1-j,:)
   ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

   DO j = 1, 5, 1
      array(1-j,:) = fac_r * array(j,:)
   ENDDO

ELSEIF(boundary_flag(1) == 2) THEN

   DO j = 1, 5, 1
      array(1-j,:) = array(1,:)
   ENDDO

ENDIF

! Do the second (r-outer) boundary

IF(boundary_flag(2) == 0) THEN

   DO j = 1, 5, 1
      array(length_step_r_part_1+j,:) = array(j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO j = 1, 5, 1
      array(length_step_r_part_1+j,:) = fac_r * array(length_step_r_part_1+1-j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO j = 1, 5, 1
      array(length_step_r_part_1+j,:) = array(length_step_r_part_1,:)                  
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_min_part_1-j) = array(:,length_step_z_part_1+1-j)                     
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO j = 1, 5, 1
      array(:,length_step_z_min_part_1-j) = fac_z * array(:,length_step_z_min_part_1-1+j)
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO j = 1, 5, 1              
      array(:,length_step_z_min_part_1-j) = array(:,length_step_z_min_part_1)
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_part_1+j) = array(:,length_step_z_min_part_1-1+j)
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_part_1+j) = fac_z * array(:,length_step_z_part_1+1-j)
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_part_1+j) = array(:,length_step_z_part_1)
   ENDDO

ENDIF

END SUBROUTINE boundary1D_DM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for dark matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_NM (array, sign)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Input parity
INTEGER, INTENT (IN) :: sign

! Dummy variables
INTEGER :: j, k

! Parity factor
INTEGER :: fac_r, fac_z

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
   fac_r = 1
   fac_z = 1
ELSEIF(sign == 1) THEN
   fac_r = -1
   fac_z = 1
ELSEIF(sign == 2) THEN
   fac_r = 1
   fac_z = -1
ELSEIF(sign == 3) THEN
   fac_r = -1
   fac_z = -1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do it case by case
IF(boundary_flag(1) == 0) THEN

   DO j = 1, 5, 1
      array(1-j,:) = array(length_step_r_part_2+1-j,:)
   ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

   DO j = 1, 5, 1
      array(1-j,:) = fac_r * array(j,:)
   ENDDO

ELSEIF(boundary_flag(1) == 2) THEN

   DO j = 1, 5, 1
      array(1-j,:) = array(1,:)
   ENDDO

ENDIF

! Do the second (r-outer) boundary

IF(boundary_flag(2) == 0) THEN

   DO j = 1, 5, 1
      array(length_step_r_part_2+j,:) = array(j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO j = 1, 5, 1
      array(length_step_r_part_2+j,:) = fac_r * array(length_step_r_part_2+1-j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO j = 1, 5, 1
      array(length_step_r_part_2+j,:) = array(length_step_r_part_2,:)                  
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_min_part_2-j) = array(:,length_step_z_part_2+1-j)                     
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO j = 1, 5, 1
      array(:,length_step_z_min_part_2-j) = fac_z * array(:,length_step_z_min_part_2-1+j)
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO j = 1, 5, 1              
      array(:,length_step_z_min_part_2-j) = array(:,length_step_z_min_part_2)
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_part_2+j) = array(:,length_step_z_min_part_2-1+j)
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_part_2+j) = fac_z * array(:,length_step_z_part_2+1-j)
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_part_2+j) = array(:,length_step_z_part_2)
   ENDDO

ENDIF

END SUBROUTINE boundary1D_NM