!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine solves the Poisson equatiuon by using the 
! Successive Over-relaxation method. Notice if you need
! to solve in other geometries, you need to modify the 
! Poisson equation accordingly. In particular, for 1D
! problem, I do not suggest using over-relaxation method
! since it will reach to a wrong solution. 
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RELAXATION_NM
USE DEFINITION
!USE OPENACC
IMPLICIT NONE

! dummy variables
INTEGER :: j, k, l

! Current Iteration
integer :: n

! Grid parity !
INTEGER :: k_min, k_sub, k_add

! Absoulte error !
REAL (DP) :: abserror, dh

! Weight !
REAL (DP) :: RBSOR_WEIGHT

! Backup array
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: phi_local
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: phiold
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: gravrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check timing with or without openmp     
!INTEGER :: time_start, time_end
!INTEGER :: cr, cm
!REAL :: rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Call the system clock !
!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)
!CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SOR Weight !
dh = (1.0D0/DBLE(max(length_step_r_2,length_step_z_2)))
RBSOR_WEIGHT = 2.0D0/(1.0D0+sin(pi*dh))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Backup
phi_local(:,:) = phi2_nm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fix the boundary first
phi_local(length_step_r_2,:) = 0.0D0
do k = 1, length_step_z_2, 1
   DO l = 0, lmax
   	phi_local(length_step_r_2,k) = phi_local(length_step_r_2,k) + qpole2(2*l)*legendre2(length_step_r_2,k,2*l)/(radbar2(length_step_r_2,k)**(2*l+1))
   END DO
   phi_local(length_step_r_2,k) = phi_local(length_step_r_2,k) * (-1.0D0) * rmax2 **2
enddo              

! Another boundaries !
phi_local(:,length_step_z_2) = 0.0D0
do j = 1, length_step_r_2, 1
   DO l = 0, lmax
   	phi_local(j,length_step_z_2) = phi_local(j,length_step_z_2) + qpole2(2*l)*legendre2(j,length_step_z_2,2*l)/(radbar2(j,length_step_z_2)**(2*l+1))
   END DO
   phi_local(j,length_step_z_2) = phi_local(j,length_step_z_2) * (-1.0D0) * rmax2 **2
enddo          

! Speical care of boundaries for hemisphere 
IF(hemisphere_flag == 1) THEN
   phi_local(:,1) = 0.0D0
   DO j = 1, length_step_r_2, 1 
      DO l = 0, lmax
	phi_local(j,1) = phi_local(j,1) + qpole2(2*l)*legendre2(j,1,2*l)/(radbar2(j,1)**(2*l+1))
      END DO
      phi_local(j,1) = phi_local(j,1) * (-1.0D0) * rmax2 **2
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign gravitational mass !
DO j = 1, length_step_r_2
	DO k = 1, length_step_z_2
		IF(rho2(j,k) > rho2_a) THEN
			gravrho(j,k) = rho2(j,k)
		ELSE
			gravrho(j,k) = 0.0D0
		END IF
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(hemisphere_flag == 0) THEN
	k_min = 1
ELSEIF(hemisphere_flag == 1) THEN
	k_min = 2
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Relaxation loops 
!$acc data copyin(gravrho(:,:), phiold(:,:), r2(:), abserror, rbsor_weight, dx2, k_min) copy(phi_local(:,:), k_add, k_sub)
DO n = 1, relax_max, 1

   ! Backup first to find the residue term
   !$acc parallel loop gang
   DO j = -4, length_step_r_2 + 5
	!$acc loop vector
	DO k = -4, length_step_z_2 + 5
		phiold(j, k) = phi_local(j, k)
	END DO
   END DO
   !$acc end parallel

   ! Set error !
   !$acc serial	
   abserror = 1.0E-100_DP
   !$acc end serial

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Here I try to implement red-black relaxation
   ! Red chessboard
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF(coordinate_flag == 0) THEN

     STOP 'In PotentialRelax: No relaxation when using Cartesian'

   ELSEIF(coordinate_flag == 1) THEN

      !$acc parallel loop gang 
      DO k = k_min, length_step_z_2 - 1, 1
	 k_add = MOD(k,2)
	 k_sub = 1 - k_add
         !$acc loop vector
         DO j = length_step_r_2 - 1 - k_sub, 1, -2
            phi_local(j,k) = (1.0D0 - rbsor_weight) * phi_local(j,k) + (rbsor_weight) * &
	 	 	     (0.25D0 * (phi_local(j+1,k) + phi_local(j-1,k) + phi_local(j,k+1) + phi_local(j,k-1)) + &
			     0.125D0 * (phi_local(j+1,k) - phi_local(j-1,k)) * dx2 / r2(j) - & 
			     pi * dx2 ** 2 * gravrho(j,k))		
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   !cot_theta = 1.0D0 / DTAN(z_temp(j,k))
	   !phinew (j, k) = 0.25D0 * (phi (j+1, k) + phi (j-1, k) + phi(j, k+1) + phi(j, k-1)) + &
   	   !	     0.125D0 * dx2 / ((DBLE (j) - 5.0E-1_DP) * dx2) * & 
 	   !	     (phi (j+1,k) - phi (j-1,k)) - pi * dx2 ** 2 * rho2 (j,k)
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         END DO
      ENDDO
      !$acc end parallel

   ENDIF
	
   ! Copy to boundary 
   !$acc parallel loop gang
   DO j = 1, 5
	phi_local(1-j, :) = phi_local(j, :)
	phi_local(:, 1-j) = phi_local(:, j)
   ENDDO
   !$acc end parallel

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Black chessboard
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF(coordinate_flag == 0) THEN

     STOP 'In PotentialRelax: No relaxation when using Cartesian'

   ELSEIF(coordinate_flag == 1) THEN

      !$acc parallel loop gang 
      DO k = k_min, length_step_z_2 - 1, 1
	 k_add = MOD(k,2)
	 k_sub = 1 - k_add
         !$acc loop vector
         DO j = length_step_r_2 - 1 - k_add, 1, -2
   	    phi_local(j,k) = (1.0D0 - rbsor_weight) * phi_local(j,k) + (rbsor_weight) * &
                             (0.25D0 * (phi_local(j+1,k) + phi_local(j-1,k) + phi_local(j,k+1) + phi_local(j,k-1)) + &
                             0.125D0 * (phi_local(j+1,k) - phi_local(j-1,k)) * dx2 / r2(j) - &
                             pi * dx2 ** 2 * gravrho(j,k))
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   !cot_theta = 1.0D0 / DTAN(z_temp(j,k))        
           !phinew (j, k) = 0.25D0 * (phi (j+1, k) + phi (j-1, k) + phi(j, k+1) + phi(j, k-1)) + &
           !            0.125D0 * dx2 / ((DBLE (j) - 5.0E-1_DP) * dx2) * &
           !            (phi (j+1,k) - phi (j-1,k)) - pi * dx2 ** 2 * rho2 (j,k)
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 END DO
      ENDDO
      !$acc end parallel
	
   ENDIF

   ! Copy to boundary 
   !$acc parallel loop gang
   DO j = 1, 5
	phi_local(1-j, :) = phi_local(j, :)
	phi_local(:, 1-j) = phi_local(:, j)
   ENDDO
   !$acc end parallel

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Here, I check the error
   !$acc parallel loop gang reduction(MAX:abserror)
   DO k = 1, length_step_z_2 - 1, 1
      !$acc loop vector
      DO j = 1, length_step_r_2 - 1, 1
	abserror = max(abserror, abs((phi_local(j,k) - phiold (j,k)) / phiold (j,k)))
      ENDDO
   END DO
   !$acc end parallel
   !$acc update self(abserror) 
	!WRITE (*,*) abserror, n
   IF(abserror <= tolerance) EXIT 

END DO
!$acc end data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now copy the result to the trunk
phi2_nm(:,:) = phi_local(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Speical care for mass cut !
IF(masscut_flag == 1) THEN

   IF(Coordinate_flag /= 2) STOP 'In PotentialRelax: No masscut without using spherical flag'

   DO k = 1, length_step_z_2, 1
      DO j = 1, length_step_r_2, 1
	 phi2_nm(j,k) = phi2_nm(j,k) - mass_cut / r2(j)
      ENDDO
   ENDDO

   CALL BOUNDARY1D_NMFULL(phi2_nm, even)

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Finalize, and check the timing !
!CALL system_clock(time_end)
!WRITE(*,*) 'Total time elapsed is = ', REAL(time_end - time_start) / rate

! Format !
100 format (3I5, ES18.8)

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Relaxation method, but for DM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RELAXATION_DM
USE DEFINITION
!USE OPENACC
IMPLICIT NONE

! dummy variables
INTEGER :: j, k, l

! Current Iteration
integer :: n

! Grid parity !
INTEGER :: k_min, k_sub, k_add

! Absoulte error !
REAL (DP) :: abserror, dh

! Weight !
REAL (DP) :: RBSOR_WEIGHT

! Backup array
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: phi_local
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: phiold
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: gravrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check timing with or without openmp     
!INTEGER :: time_start, time_end
!INTEGER :: cr, cm
!REAL :: rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Call the system clock !
!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)
!CALL system_clock(time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SOR Weight !
dh = (1.0D0/DBLE(max(length_step_r_1,length_step_z_1)))
RBSOR_WEIGHT = 2.0D0/(1.0D0+sin(pi*dh))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Backup
phi_local(:,:) = phi1_dm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fix the boundary first
phi_local(length_step_r_1,:) = 0.0D0
do k = 1, length_step_z_1, 1
   DO l = 0, lmax
   	phi_local(length_step_r_1,k) = phi_local(length_step_r_1,k) + qpole1(2*l)*legendre1(length_step_r_1,k,2*l)/(radbar1(length_step_r_1,k)**(2*l+1))
   END DO
   phi_local(length_step_r_1,k) = phi_local(length_step_r_1,k) * (-1.0D0) * rmax1 **2
enddo              

! Another boundaries !
phi_local(:,length_step_z_1) = 0.0D0
do j = 1, length_step_r_1, 1
   DO l = 0, lmax
   	phi_local(j,length_step_z_1) = phi_local(j,length_step_z_1) + qpole1(2*l)*legendre1(j,length_step_z_1,2*l)/(radbar1(j,length_step_z_1)**(2*l+1))
   END DO
   phi_local(j,length_step_z_1) = phi_local(j,length_step_z_1) * (-1.0D0) * rmax1 **2
enddo          

! Speical care of boundaries for hemisphere 
IF(hemisphere_flag == 1) THEN
   phi_local(:,1) = 0.0D0
   DO j = 1, length_step_r_1, 1 
      DO l = 0, lmax
	phi_local(j,1) = phi_local(j,1) + qpole1(2*l)*legendre1(j,1,2*l)/(radbar1(j,1)**(2*l+1))
      END DO
      phi_local(j,1) = phi_local(j,1) * (-1.0D0) * rmax1 **2
   ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign gravitational mass !
DO j = 1, length_step_r_1
	DO k = 1, length_step_z_1
		IF(rho1(j,k) > rho1_a) THEN
			gravrho(j,k) = rho1(j,k)
		ELSE
			gravrho(j,k) = 0.0D0
		END IF
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(hemisphere_flag == 0) THEN
	k_min = 1
ELSEIF(hemisphere_flag == 1) THEN
	k_min = 2
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Relaxation loops 
!$acc data copyin(gravrho(:,:), phiold(:,:), r1(:), abserror, rbsor_weight, dx1, k_min) copy(phi_local(:,:), k_add, k_sub)
DO n = 1, relax_max, 1

   ! Backup first to find the residue term
   !$acc parallel loop gang
   DO j = -4, length_step_r_1 + 5
	!$acc loop vector
	DO k = -4, length_step_z_1 + 5
		phiold(j, k) = phi_local(j, k)
	END DO
   END DO
   !$acc end parallel

   ! Set error !
   !$acc serial	
   abserror = 1.0E-100_DP
   !$acc end serial

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Here I try to implement red-black relaxation
   ! Red chessboard
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF(coordinate_flag == 0) THEN

     STOP 'In PotentialRelax: No relaxation when using Cartesian'

   ELSEIF(coordinate_flag == 1) THEN

      !$acc parallel loop gang 
      DO k = k_min, length_step_z_1 - 1, 1
	 k_add = MOD(k,2)
	 k_sub = 1 - k_add
         !$acc loop vector
         DO j = length_step_r_1 - 1 - k_sub, 1, -2
            phi_local(j,k) = (1.0D0 - rbsor_weight) * phi_local(j,k) + (rbsor_weight) * &
		      	     (0.25D0 * (phi_local(j+1,k) + phi_local(j-1,k) + phi_local(j,k+1) + phi_local(j,k-1)) + &
			     0.125D0 * (phi_local(j+1,k) - phi_local(j-1,k)) * dx1 / r1(j) - & 
			     pi * dx1** 2 * gravrho(j,k))
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   !cot_theta = 1.0D0 / DTAN(z_temp(j,k))
	   !phinew (j, k) = 0.25D0 * (phi (j+1, k) + phi (j-1, k) + phi(j, k+1) + phi(j, k-1)) + &
   	   !	     0.125D0 * dx1 / ((DBLE (j) - 5.0E-1_DP) * dx1) * & 
 	   !	     (phi (j+1,k) - phi (j-1,k)) - pi * dx1 ** 2 * rho1 (j,k)
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	END DO
      ENDDO
      !$acc end parallel

   ENDIF
	
   ! Copy to boundary 
   !$acc parallel loop gang
   DO j = 1, 5
	phi_local(1-j, :) = phi_local(j, :)
	phi_local(:, 1-j) = phi_local(:, j)
   ENDDO
   !$acc end parallel

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Black chessboard
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF(coordinate_flag == 0) THEN

     STOP 'In PotentialRelax: No relaxation when using Cartesian'

   ELSEIF(coordinate_flag == 1) THEN
 
      !$acc parallel loop gang 
      DO k = k_min, length_step_z_1 - 1, 1
	 k_add = MOD(k,2)
	 k_sub = 1 - k_add
         !$acc loop vector
         DO j = length_step_r_1 - 1 - k_add, 1, -2
   	    phi_local(j,k) = (1.0D0 - rbsor_weight) * phi_local(j,k) + (rbsor_weight) * &
                             (0.25D0 * (phi_local(j+1,k) + phi_local(j-1,k) + phi_local(j,k+1) + phi_local(j,k-1)) + &
                             0.125D0 * (phi_local(j+1,k) - phi_local(j-1,k)) * dx1 / r1(j) - &
                             pi * dx1 ** 2 * gravrho(j,k))
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   !cot_theta = 1.0D0 / DTAN(z_temp(j,k))        
           !phinew (j, k) = 0.25D0 * (phi (j+1, k) + phi (j-1, k) + phi(j, k+1) + phi(j, k-1)) + &
           !            0.125D0 * dx1 / ((DBLE (j) - 5.0E-1_DP) * dx1) * &
           !            (phi (j+1,k) - phi (j-1,k)) - pi * dx1 ** 2 * rho1 (j,k)
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         END DO
      ENDDO
      !$acc end parallel

   ENDIF
   
   ! Copy to boundary 
   !$acc parallel loop gang
   DO j = 1, 5
	phi_local(1-j, :) = phi_local(j, :)
	phi_local(:, 1-j) = phi_local(:, j)
   ENDDO
   !$acc end parallel

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Here, I check the error
   !$acc parallel loop gang reduction(MAX:abserror)
   DO k = 1, length_step_z_2 - 1, 1
      !$acc loop vector
      DO j = 1, length_step_r_2 - 1, 1
	abserror = max(abserror, abs((phi_local(j,k) - phiold (j,k)) / phiold (j,k)))
      ENDDO
   END DO
   !$acc end parallel
   !$acc update self(abserror) 
	!WRITE (*,*) abserror, n
   IF(abserror <= tolerance) EXIT 

END DO
!$acc end data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now copy the result to the trunk
phi1_dm(:,:) = phi_local(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Speical care for mass cut !
IF(masscut_flag == 1) THEN

   IF(coordinate_flag /= 2) STOP 'In PotentialRelax: No masscut without using spherical flag'

   DO k = 1, length_step_z_1, 1
      DO j = 1, length_step_r_1, 1
	 phi1_dm(j,k) = phi1_dm(j,k) - mass_cut / r1(j)
      ENDDO
   ENDDO

   CALL BOUNDARY1D_DMFULL(phi1_dm, even)

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Finalize, and check the timing !
!CALL system_clock(time_end)
!WRITE(*,*) 'Total time elapsed is = ', REAL(time_end - time_start) / rate

! Format !
100 format (3I5, ES18.8)

END SUBROUTINE 
