!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for WENO
! Written by Leung Shing Chi in 2016
! The subroutine takes ARRAY as input/output and SIGN
! for doing odd/even parity extension
! Notice that this subroutines worked for full size arrays 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_DMFULL (array, sign)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      array(1-j,:) = array(length_step_r_1+1-j,:)
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
      array(length_step_r_1+j,:) = array(j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO j = 1, 5, 1
      array(length_step_r_1+j,:) = fac_r * array(length_step_r_1+1-j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO j = 1, 5, 1
      array(length_step_r_1+j,:) = array(length_step_r_1,:)                  
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO j = 1, 5, 1
      array(:,1-j) = array(:,length_step_z_1+1-j)                     
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO j = 1, 5, 1
      array(:,1-j) = fac_z * array(:,j)
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO j = 1, 5, 1              
      array(:,1-j) = array(:,1)
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_1+j) = array(:,j)
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_1+j) = fac_z * array(:,length_step_z_1+1-j)
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_1+j) = array(:,length_step_z_1)
   ENDDO

ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for Normal matter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D_NMFULL (array, sign)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      array(1-j,:) = array(length_step_r_2+1-j,:)
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
      array(length_step_r_2+j,:) = array(j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO j = 1, 5, 1
      array(length_step_r_2+j,:) = fac_r * array(length_step_r_2+1-j,:)
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO j = 1, 5, 1
      array(length_step_r_2+j,:) = array(length_step_r_2,:)                  
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO j = 1, 5, 1
      array(:,1-j) = array(:,length_step_z_2+1-j)                     
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO j = 1, 5, 1
      array(:,1-j) = fac_z * array(:,j)
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO j = 1, 5, 1              
      array(:,1-j) = array(:,1)
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_2+j) = array(:,j)
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_2+j) = fac_z * array(:,length_step_z_2+1-j)
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO j = 1, 5, 1
      array(:,length_step_z_2+j) = array(:,length_step_z_2)
   ENDDO

ENDIF

END SUBROUTINE