!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine ensures that the density does not go below rho_atm
! and replace the grids with atmospheric density (tempearture and 
! chemical composition and so on) if found
! Written by Leung Shing Chi in 2016
! The subroutine do all the check automatically
! Notice that this subroutines also check the size of 
! the hydro array, and reduced the simulation grid-number
! to boost the calculation,
! i.e. (1:length_step_r_part, 1:length_step_z_part).
! For full array extension, switch the checkstep_flag = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKRHO
USE DEFINITION
USE HELMEOS_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k

! Flag for finding the outer boundary
INTEGER :: found

! Threshold for atmosphere density
real (DP) :: rho_min1, rho_min2

! Temporal variables !
INTEGER :: temp_step_1, temp_step_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'In Check Rho'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First check the DM
if(DM_flag == 1) then

   ! Set up the threshold
   rho_min1 = 1.1D0 * rho1_a

   ! Initialize
   found = 0

   DO k = -4, length_step_z_1 + 5, 1
      DO j = -4, length_step_r_1 + 5, 1

	 ! IF low-density grid is found
	 ! Replace them with stmospheric
	 ! density and velocity and temp
         if(rho1(j,k) <= rho_min1) then
            rho1 (j,k) = rho1_a
            vel1_r (j,k) = vel1_a
            vel1_z (j,k) = vel1_a  
	    vel1_p (j,k) = vel1_a
            epsilon1 (j,k) = epsilon1_a
         endif

      ENDDO
   ENDDO

   ! Make sure the update also applies to ghost cells
   CALL BOUNDARY1D_DMFULL (rho1, even)              
   CALL BOUNDARY1D_DMFULL (vel1_r, oddR)
   CALL BOUNDARY1D_DMFULL (vel1_z, oddZ)         
   CALL BOUNDARY1D_DMFULL (vel1_p, oddR)
   CALL BOUNDARY1D_DMFULL (epsilon1, even)

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do the normal matter

! Set up the threshold density
rho_min2 = 1.1D0 * rho2_a

do k = -4, length_step_z_2 + 5, 1
   do j = -4, length_step_r_2 + 5, 1          

      ! Check the density of normal matter
      ! Replace the low-density grid by the 
      ! atmospheric one with default density
      ! tempeature, velocity and electron fraction

      if(rho2(j,k) <= rho_min2) then   
         rho2 (j,k) = rho2_a
         vel2_r (j,k) = vel2_a
         vel2_z (j,k) = vel2_a
	 vel2_p (j,k) = vel2_a
         epsilon2 (j,k) = epsilon2_a
      endif

   enddo
enddo

! Make sure the ghost cell knows the udpate
CALL BOUNDARY1D_NMFULL (rho2, even)
CALL BOUNDARY1D_NMFULL (vel2_r, oddR)
CALL BOUNDARY1D_NMFULL (vel2_z, oddZ)
CALL BOUNDARY1D_NMFULL (vel2_p, oddR)
CALL BOUNDARY1D_NMFULL (epsilon2, even)

! Temperauter !
if(dual_energy == 1) then
   do k = -4, length_step_z_2 + 5, 1
      do j = -4, length_step_r_2 + 5, 1
         if(rho2(j,k) <= rho_min2) then
            rhoe2 (j,k) = rho2_a*epsilon2_a
	 endif
      enddo
   enddo
   CALL BOUNDARY1D_NMFULL (rhoe2, even)
endif

! Temperauter !
if(helmeos_flag == 1) then
   do k = -4, length_step_z_2 + 5, 1
      do j = -4, length_step_r_2 + 5, 1
         if(rho2(j,k) <= rho_min2) then
            temp2 (j,k) = temp_a
	 endif
      enddo
   enddo
   CALL BOUNDARY1D_NMFULL (temp2, even)
endif

! Electron fraction !
if(etran_flag == 1) then
   do k = -4, length_step_z_2 + 5, 1
      do j = -4, length_step_r_2 + 5, 1
         if(rho2(j,k) <= rho_min2) then
            ye2 (j,k) = ye_a
	 endif
      enddo
   enddo
   call boundary1D_NMFULL (ye2, even)   
endif

! Also replace the chemical composition
! This part is difficult if the default composition
! has density dependent mixture
if(xisotran_flag == 1) then
   do k = -4, length_step_z_2 + 5, 1
      do j = -4, length_step_r_2 + 5, 1
         if(rho2(j,k) <= rho_min2) then
            xiso(j,k,:) = xiso_a
	 endif
      enddo
   enddo

   ! Check and copy to boundaries !
   call checkxisotope
   CALL BOUNDARY2D_X
endif

! Also replace the turbulence strength
IF(turb_flag == 1) then
   do k = -4, length_step_z_2 + 5, 1
      do j = -4, length_step_r_2 + 5, 1
         if(rho2(j,k) <= rho_min2) then
	    turb_q(j,k) = turb_q_a
         endif
      enddo
   enddo
  
   ! Ghost shell !
   CALL BOUNDARY1D_NMFULL (turb_q, even)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box
! which can contain the whole star, but minimize the calculation

IF(checkstepdm_flag == 1) THEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_r_part_1

   ! Initialized r_grid1
   r_grid1 = 0          

   ! Find the outermost grid which is the star surface    
   DO k = 1, length_step_z_1, 1     
      DO j = 1, length_step_r_1 - 1, 1
         if(rho1(j,k) >= rho_min1 .and. rho1(j+1,k) == rho1_a) then
	    if(j > r_grid1) then
               r_grid1 = j
	    endif
	    !exit
         endif                 
      ENDDO
   ENDDO

   temp_step_1 = length_step_r_part_1
   ! Set the effective length_step_r
   IF(r_grid1 == 0) then
      length_step_r_part_1 = length_step_r_1
   ELSEIF(r_grid1 /= 0 .and. r_grid1 < length_step_r_1 - 15) then
      length_step_r_part_1 = r_grid1 + 15
   ELSE
      length_step_r_part_1 = length_step_r_1
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_r_part_1 = max(temp_step_1, length_step_r_part_1)
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_z_part_1

   IF(hemisphere_flag == 1) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid1
      z_grid1 = length_step_z_1/2-1

      ! Find the outermost grid which is the star surface
      DO j = 1, length_step_r_1, 1
         DO k = length_step_z_1/2-1, 2, -1
            IF(rho1(j,k) >= rho_min1 .and. rho1(j,k-1) == rho1_a) THEN
	       IF(k < z_grid1) THEN
                  z_grid1 = k
   	       ENDIF
	       !exit
            ENDIF
         ENDDO
      ENDDO 

      temp_step_1 = length_step_z_min_part_1
      ! Set the effective length_step_z
      IF(z_grid1 == length_step_z_1/2-1) then
         length_step_z_min_part_1 = 1
      ELSEIF(z_grid1 /= length_step_z_1/2-1 .and. z_grid1 > 15) THEN
         length_step_z_min_part_1 = z_grid1 - 15
      ELSE
         length_step_z_min_part_1 = 1
      ENDIF
      IF(global_time > 0.0D0) THEN
         length_step_z_min_part_1 = min(temp_step_1, length_step_z_min_part_1)
      END IF

   ELSEIF(hemisphere_flag == 0) THEN

      length_step_z_min_part_1 = 1

   ELSE

      STOP 'Check the value of hemisphere_flag'

   ENDIF

   ! Now find the upper limit of length_step_z
   ! Initialized z_grid1
   z_grid1 = 0

   ! Find the outermost grid which is the star surface
   DO j = 1, length_step_r_1, 1
      DO k = 1, length_step_z_1 - 1, 1
         IF(rho1(j,k) >= rho_min1 .and. rho1(j,k+1) == rho1_a) THEN
            IF(k > z_grid1) then
               z_grid1 = k
            ENDIF
            !exit
         ENDIF
      ENDDO
   ENDDO

   temp_step_1 = length_step_z_part_1
   ! Set the effective length_step_z
   IF(z_grid1 == 0) then
      length_step_z_part_1 = length_step_z_1
   ELSEIF(z_grid1 /= 0 .and. z_grid1 < length_step_z_1 - 15) THEN
      length_step_z_part_1 = z_grid1 + 15
   ELSE
      length_step_z_part_1 = length_step_z_1
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_z_part_1 = max(temp_step_1, length_step_z_part_1)
   END IF

ELSE

   ! Nothing changed if you do not need the check
   length_step_z_min_part_1 = 1
   length_step_z_part_1 = length_step_z_1
   length_step_r_part_1 = length_step_r_1

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box
! which can contain the whole star, but minimize the calculation

IF(checkstepnm_flag == 1) THEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_r_part

   ! Initialized r_grid2
   r_grid2 = 0          

   ! Find the outermost grid which is the star surface    
   DO k = 1, length_step_z_2, 1     
      DO j = 1, length_step_r_2 - 1, 1
         if(rho2(j,k) >= rho_min2 .and. rho2(j+1,k) == rho2_a) then
	    if(j > r_grid2) then
               r_grid2 = j
	    endif
	    !exit
         endif                 
      ENDDO
   ENDDO

   temp_step_2 = length_step_r_part_2
   ! Set the effective length_step_r
   IF(r_grid2 == 0) then
      length_step_r_part_2 = length_step_r_2
   ELSEIF(r_grid2 /= 0 .and. r_grid2 < length_step_r_2 - 15) then
      length_step_r_part_2 = r_grid2 + 15
   ELSE
      length_step_r_part_2 = length_step_r_2
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_r_part_2 = max(temp_step_2, length_step_r_part_2)
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Obtain length_step_z_part

   IF(hemisphere_flag == 1) THEN

      ! Now find the lower limit of length_step_z
      ! Initialized z_grid2
      z_grid2 = length_step_z_2/2-1

      ! Find the outermost grid which is the star surface
      DO j = 1, length_step_r_2, 1
         DO k = length_step_z_2/2-1, 2, -1
            IF(rho2(j,k) >= rho_min2 .and. rho2(j,k-1) == rho2_a) THEN
	       IF(k < z_grid2) THEN
                  z_grid2 = k
   	       ENDIF
	       !exit
            ENDIF
         ENDDO
      ENDDO 

      temp_step_2 = length_step_z_min_part_2
      ! Set the effective length_step_z
      IF(z_grid2 == length_step_z_2/2-1) then
         length_step_z_min_part_2 = 1
      ELSEIF(z_grid2 /= length_step_z_2/2-1 .and. z_grid2 > 15) THEN
         length_step_z_min_part_2 = z_grid2 - 15
      ELSE
         length_step_z_min_part_2 = 1
      ENDIF
      IF(global_time > 0.0D0) THEN
         length_step_z_min_part_2 = min(temp_step_2, length_step_z_min_part_2)
      END IF

   ELSEIF(hemisphere_flag == 0) THEN

      length_step_z_min_part_2 = 1

   ELSE

      STOP 'Check the value of hemisphere_flag'

   ENDIF

   ! Now find the upper limit of length_step_z
   ! Initialized z_grid2
   z_grid2 = 0

   ! Find the outermost grid which is the star surface
   DO j = 1, length_step_r_2, 1
      DO k = 1, length_step_z_2 - 1, 1
         IF(rho2(j,k) >= rho_min2 .and. rho2(j,k+1) == rho2_a) THEN
            IF(k > z_grid2) then
               z_grid2 = k
            ENDIF
            !exit
         ENDIF
      ENDDO
   ENDDO

   temp_step_2 = length_step_z_part_2
   ! Set the effective length_step_z
   IF(z_grid2 == 0) then
      length_step_z_part_2 = length_step_z_2
   ELSEIF(z_grid2 /= 0 .and. z_grid2 < length_step_z_2 - 15) THEN
      length_step_z_part_2 = z_grid2 + 15
   ELSE
      length_step_z_part_2 = length_step_z_2
   ENDIF
   IF(global_time > 0.0D0) THEN
      length_step_z_part_2 = max(temp_step_2, length_step_z_part_2)
   END IF

ELSE

   ! Nothing changed if you do not need the check
   length_step_z_min_part_2 = 1
   length_step_z_part_2 = length_step_z_2
   length_step_r_part_2 = length_step_r_2

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE CHECKRHO