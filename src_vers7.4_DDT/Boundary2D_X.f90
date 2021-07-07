!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for the chemical isotopes
!
! The subroutine automatically extends the isotopes
! assuming all isotopes behave like scalars
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, switch the checkstep_flag = 0
! For hybrid boundaries, such as the quadrant star
! Specific modifications are needed
!
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY2D_X()
USE definition
USE helmeos_module
IMPLICIT NONE

!Dummy variables
integer :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(boundary_flag(1) == 0) THEN

   DO j = 1, 5, 1
      xiso(1-j,:,:) = xiso(length_step_r_part_2+1-j,:,:)
   ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

   DO j = 1, 5, 1
      xiso(1-j,:,:) = xiso(j,:,:)
   ENDDO

ELSEIF(boundary_flag(1) == 2) THEN

   DO j = 1, 5, 1
      xiso(1-j,:,:) = xiso(1,:,:)
   ENDDO

ENDIF

! Do the second (r-outer) boundary

IF(boundary_flag(2) == 0) THEN

   DO j = 1, 5, 1
      xiso(length_step_r_part_2+j,:,:) = xiso(j,:,:)
   ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

   DO j = 1, 5, 1
      xiso(length_step_r_part_2+j,:,:) = xiso(length_step_r_part_2+1-j,:,:)
   ENDDO

ELSEIF(boundary_flag(2) == 2) THEN

   DO j = 1, 5, 1
      xiso(length_step_r_part_2+j,:,:) = xiso(length_step_r_part_2,:,:)                  
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the third (z-inner) boundary

IF(boundary_flag(3) == 0) THEN

   DO j = 1, 5, 1
      xiso(:,length_step_z_min_part_2-j,:) = xiso(:,length_step_z_part_2+1-j,:)                     
   ENDDO

ELSEIF(boundary_flag(3) == 1) THEN                 

   DO k = 1, no_of_eq, 1
      DO j = 1, 5, 1
         xiso(:,length_step_z_min_part_2-j,:) = xiso(:,length_step_z_min_part_2-1+j,:)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(3) == 2) THEN

   DO k = 1, no_of_eq, 1
      DO j = 1, 5, 1              
         xiso(:,length_step_z_min_part_2-j,:) = xiso(:,length_step_z_min_part_2,:)
      ENDDO             
   ENDDO

ENDIF

! Do the fourth (z-outer) boundary

IF(boundary_flag(4) == 0) THEN

   DO k = 1, no_of_eq, 1
      DO j = 1, 5, 1
         xiso(:,length_step_z_part_2+j,:) = xiso(:,length_step_z_min_part_2-1+j,:)
      ENDDO
   ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

   DO j = 1, 5, 1
      xiso(:,length_step_z_part_2+j,:) = xiso(:,length_step_z_part_2+1-j,:)
   ENDDO

ELSEIF(boundary_flag(4) == 2) THEN

   DO j = 1, 5, 1
      xiso(:,length_step_z_part_2+j,:) = xiso(:,length_step_z_part_2,:)
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE boundary2D_X