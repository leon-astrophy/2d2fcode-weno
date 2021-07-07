!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 
! This subroutine provides a quick build of initial profiles 
! based on some well known tests 
! It includes: 
! 1. TORO 1D shocktube
! 2. Sedov Spherical test
! 3. 2D explosion
! 4. Kelvin-helmholtz
! 5. Rayleigh-Taylor
! 6. Implosion
! 7. Gresho problem
!
! Written by Leung Shing Chi in 2015 
! Updated by Leung Shing Chi in 2017 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GETRHO_template
USE DEFINITION
USE helmeos_module
IMPLICIT NONE

INTEGER :: i, j, k
REAL (DP) :: rad_dist
REAL (DP) :: dummy

IF(testmodel_flag == 1) THEN

   ! Test 1, sod shocktube test Use gamma = 1.4
   gamma2 = 1.4E0_DP
   ! Set 3
   DO i = 1, length_step_z_2 * 4 / 10, 1
      rho2(:,i) = 5.99924D0
      p2(:,i) = 460.894D0
      vel2_r(:,i) = 0.0D0
      vel2_z(:,i) = 19.5975D0 
      epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0) 
   ENDDO

   DO i = length_step_z_2 * 4 / 10 + 1, length_step_z_2, 1
      rho2(:,i) = 5.99242D0
      p2(:,i) = 46.0950D0
      vel2_r(:,i) = 0.0D0
      vel2_z(:,i) = -6.19633D0
      epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0) 
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set 4
   !DO i = 1, length_step_z_2 * 8 / 10, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 1000.0D0
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = -19.59745D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !DO i = length_step_z_2 * 8 / 10 + 1, length_step_z_2, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 0.01D0
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = -19.59745D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! Set 1 SR (See Rosswog 2013 SR-SPH) aka weak blast test (gamma = 4/3)
   !DO i = 1, length_step_z_2 / 2, 1
   !   rho2(:,i) = 1.0D0       
   !   p2(:,i) = 1.0D0        
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !DO i = length_step_z_2 / 2 + 1, length_step_z_2, 1
   !   rho2(:,i) = 0.125D0
   !   p2(:,i) = 0.1D0 
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0 
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set 2 SR (See Rosswog 2013 SR-SPH) aka intermediate blast test (gamma = 4/3)
   !DO i = 1, length_step_z_2 / 2, 1
   !   rho2(:,i) = 10.0D0        
   !   p2(:,i) = 13.33D0       
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !DO i = length_step_z_2 / 2 + 1, length_step_z_2, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 0.1D0
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0 
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set 3 SR (See Rosswog 2013 SR-SPH) aka strong blast test (gamma = 5/3)
   ! See also Tominaga 2009 SR Test 1 aka strong blast test (gamma = 5/3)
   !DO i = 1, length_step_z_2 / 2, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 1000.0D0
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !DO i = length_step_z_2 / 2 + 1, length_step_z_2, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 0.01D0
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set 4 SR (See Tominaga2009 Test 2) (gamma = 4/3)
   !DO i = 1, length_step_z_2 / 2, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 1.0D0
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 1.0D0 - 1.0D-6
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !DO i = length_step_z_2 / 2 + 1, length_step_z_2, 1
   !   rho2(:,i) = 1.0D0
   !   p2(:,i) = 1.0D6
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.0D0
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set 5 SR (See Tominaga2009 Test 3) (gamma = 4/3)
   !DO i = 1, length_step_z_2 / 2, 1
   !   rho2(:,i) = 1.0D0      
   !   p2(:,i) = 10.0D0         
   !   vel2_r(:,i) = 0.0D0       
   !   vel2_z(:,i) = -0.9D0   
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !DO i = length_step_z_2 / 2 + 1, length_step_z_2, 1
   !   rho2(:,i) = 10.0D0
   !   p2(:,i) = 100.0D0               
   !   vel2_r(:,i) = 0.0D0
   !   vel2_z(:,i) = 0.9D0    
   !   epsilon2(:,i) = p2(:,i) / rho2(:,i) / (gamma2 - 1.0D0)
   !ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, even)
   CALL boundary1D_NM(vel2_z, even)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 2) THEN

   ! Test 2, 2D shock test with gamma = 1.4
   gamma2 = 1.4E0_DP
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1

         if(DSQRT((r2(i) - 0.5D0) ** 2 + (z2(j) - 0.75D0) ** 2) < 0.1D0) then       
            rho2(i,j) = 1.0D0       
            p2(i,j) = 10.0D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         else
            rho2(i,j) = 1.0D0
            p2(i,j) = 1.0D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         endif
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 3) THEN

   ! Test 3, 2D Explosion (Toro1997) with gamma = 1.4
   gamma2 = 1.4E0_DP
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1

         IF(DSQRT(r2(i) ** 2 + z2(j) ** 2) < 0.4D0) THEN
            rho2(i,j) = 1.0D0       
            p2(i,j) = 1.0D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         ELSE
            rho2(i,j) = 0.125D0
            p2(i,j) = 0.1D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         ENDIF
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

elseif(testmodel_flag == 4) then

   ! Test 4, Kelvin-Helmholtz with gamma = 1.4
   gamma2 = 1.4E0_DP
   do j = 1, length_step_z_2, 1
      do i = 1, length_step_r_2, 1

         if(j <= length_step_z_2 * 3 / 4 .and. j >= length_step_z_2 / 4) then
            rho2(i,j) = 2.0D0 
            p2(i,j) = 2.5D0      
            vel2_r(i,j) = -0.5D0 
            vel2_z(i,j) = 0.05D0 * DSIN(4.0D0 * pi * r2(i))
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         else                  
            rho2(i,j) = 1.0D0
            p2(i,j) = 2.5D0
            vel2_r(i,j) = 0.5D0
            vel2_z(i,j) = 0.05D0 * DSIN(4.0D0 * pi * r2(i))
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         endif
      enddo
   enddo

   call boundary1D_NM(rho2, even)
   call boundary1D_NM(vel2_r, oddR)
   call boundary1D_NM(vel2_z, oddZ)
   call boundary1D_NM(p2, even)     
   call boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 5) THEN

   ! Test 5, Rayleigh-Taylor
   gamma2 = 1.4E0_DP
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1

         IF(j <= length_step_z_2 * 1 / 2) THEN
            rho2(i,j) = 1.0D0
            p2(i,j) = 2.5D0 - 0.1D0 * rho2(i,j) * z2(j)
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.025D0 * (1.0D0 + DCOS(4.0D0 * pi * (r2(i) - 0.25D0))) * &
                                    (1.0D0 + DCOS(3.0D0 * pi * (z2(j) - 0.75D0)))
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         ELSE
            rho2(i,j) = 2.0D0
            p2(i,j) = 2.5D0 - 0.1D0 * rho2(i,j) * z2(j)
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.025D0 * (1.0D0 + DCOS(4.0D0 * pi * (r2(i) - 0.25D0))) * &
                                    (1.0D0 + DCOS(3.0D0 * pi * (z2(j) - 0.75D0)))
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         ENDIF
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)   
   CALL boundary1D_NM(vel2_r, even)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 6) THEN

   ! Test 6, Kelvin-Helmholtz
   ! A unit box is assumed
   gamma2 = 1.4E0_DP
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1

	 IF(j <= length_step_z_2 / 4) THEN
   	    rho2(i,j) = 1.0D0 + 0.5D0 * EXP((z2(j) - 0.25D0) / 0.025D0)
	    p2(i,j) = 2.5D0
	    vel2_r(i,j) = 0.5D0 - 0.5D0 * EXP((z2(j) - 0.25D0) / 0.025D0)
	    vel2_z(i,j) = 0.01D0 * SIN(4.0D0 * pi * r2(i))
	    vel2_p(i,j) = 0.0D0
	    epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
	 ELSEIF(j > length_step_z_2 / 4 .and. j <= length_step_z_2 / 2) THEN
	    rho2(i,j) = 2.0D0 - 0.5D0 * EXP((-z2(j) + 0.25D0) / 0.025D0)
            p2(i,j) = 2.5D0
            vel2_r(i,j) = -0.5D0 + 0.5D0 * EXP((-z2(j) + 0.25D0) / 0.025D0)
            vel2_z(i,j) = 0.01D0 * SIN(4.0D0 * pi * r2(i))
            vel2_p(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
	 ELSEIF(j > length_step_z_2 / 2 .and. j <= 3 * length_step_z_2 / 4) THEN
	    rho2(i,j) = 2.0D0 - 0.5D0 * EXP(-(0.75D0 - z2(j)) / 0.025D0)
            p2(i,j) = 2.5D0
            vel2_r(i,j) = -0.5D0 + 0.5D0 * EXP((-(0.75D0 - z2(j))) / 0.025D0)
            vel2_z(i,j) = 0.01D0 * SIN(4.0D0 * pi * r2(i))
            vel2_p(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
	 ELSE
	    rho2(i,j) = 1.0D0 + 0.5D0 * EXP(-(z2(j) - 0.75D0) / 0.025D0)
            p2(i,j) = 2.5D0
            vel2_r(i,j) = 0.5D0 - 0.5D0 * EXP((-(z2(j) - 0.75D0)) / 0.025D0)
            vel2_z(i,j) = 0.01D0 * SIN(4.0D0 * pi * r2(i))
            vel2_p(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
	 ENDIF
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, even)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even) 

ELSEIF(testmodel_flag == 7) THEN

   ! Test 7, Gresho
   gamma2 = 1.4E0_DP
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1

         dummy = DSQRT((r2(i) - 0.5D0)** 2 + (z2(j) - 0.5D0)** 2)

         IF(dummy <= 0.2D0) THEN

            rho2(i,j) = 1.0D0
            p2(i,j) = 5.0D0 + 12.5D0 * (dummy ** 2)
            vel2_r(i,j) = -5.0D0 * (z2(j) - 0.5D0)
            vel2_z(i,j) = 5.0D0 * (r2(i) - 0.5D0)
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)

         ELSEIF(dummy > 0.2D0 .and. dummy <= 0.4D0) THEN

            rho2(i,j) = 1.0D0
            p2(i,j) = 9.0D0 - 4.0D0 * LOG(0.2D0) + 12.5D0 * (dummy ** 2) - 20.0D0 * dummy + 4.0D0 * LOG(dummy)
            vel2_r(i,j) = -(2.0D0 - 5.0D0 * dummy) * (z2(j) - 0.5D0) / dummy
            vel2_z(i,j) = (2.0D0 - 5.0D0 * dummy) * (r2(i) - 0.5D0) / dummy
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)

         ELSE
            rho2(i,j) = 1.0D0
            p2(i,j) = 3.0D0 + 4.0D0 * LOG(2.0D0)
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
            epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
         ENDIF

      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ) 
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 8) THEN

	! Test 8, Implosion
	gamma2 = 1.4E0_DP
	DO j = 1, length_step_z_2, 1
		DO i = 1, length_step_r_2, 1

			IF ( r2(i) + z2(j) > 0.15 ) THEN
				rho2(i,j) = 1.0E0_DP
				p2(i,j) = 1.0E0_DP
			ELSE
				rho2(i,j) = 0.125E0_DP
				p2(i,j) = 0.14E0_DP
			END IF
			
			epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
			vel2_r(i,j) = 0.0D0
			vel2_z(i,j) = 0.0D0
			vel2_p(i,j) = 0.0D0

		END DO
	END DO

	CALL boundary1D_NM(rho2, even)
	CALL boundary1D_NM(vel2_r, oddR)
	CALL boundary1D_NM(vel2_z, oddZ)
	CALL BOUNDARY1D_NM(vel2_p, oddR)
	CALL boundary1D_NM(p2, even)
	CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 9) THEN

   ! Test 9, 2D Riemann Problem Test 3 with gamma = 1.4
   gamma2 = 1.4D0
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1
	IF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            rho2(i,j) = 0.138D0       
            p2(i,j) = 0.0290D0
            vel2_r(i,j) = 1.2060D0
            vel2_z(i,j) = 1.2060D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            rho2(i,j) = 0.5323D0       
            p2(i,j) = 0.3D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 1.2060D0
	ELSEIF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            rho2(i,j) = 0.5323D0       
            p2(i,j) = 0.3D0
            vel2_r(i,j) = 1.2060D0
            vel2_z(i,j) = 0.0D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            rho2(i,j) = 1.5D0  
            p2(i,j) = 1.5D0 
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
	END IF
	epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 10) THEN

   ! Test 10, 2D Riemann Problem Test 4 with gamma = 1.4
   gamma2 = 1.4D0
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1
	IF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
	    p2(i,j) = 1.1D0
            rho2(i,j) = 1.1D0       
            vel2_r(i,j) = 0.8939D0
            vel2_z(i,j) = 0.8939D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
 	    p2(i,j) = 0.35D0
            rho2(i,j) = 0.5065D0       
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.8939D0
	ELSEIF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
	    p2(i,j) = 0.35D0
            rho2(i,j) = 0.5065D0       
            vel2_r(i,j) = 0.8939D0
            vel2_z(i,j) = 0.0D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.1D0  
            rho2(i,j) = 1.1D0 
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
	END IF
	epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 11) THEN

   ! Test 11, 2D Riemann Problem Test 6 with gamma = 1.4
   gamma2 = 1.4D0
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1
	IF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 1.0D0
            vel2_r(i,j) = -0.75D0
            vel2_z(i,j) = 0.5D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 3.0D0
            vel2_r(i,j) = -0.75D0
            vel2_z(i,j) = -0.5D0
	ELSEIF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 2.0D0
            vel2_r(i,j) = 0.75D0
            vel2_z(i,j) = 0.5D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.0D0  
            rho2(i,j) = 1.0D0 
            vel2_r(i,j) = 0.75D0
            vel2_z(i,j) = -0.5D0
	END IF
	epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 12) THEN

   ! Test 12, 2D Riemann Problem Test 12 with gamma = 1.4
   gamma2 = 1.4D0
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1
	IF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 0.8D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 1.0D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.7276D0
	ELSEIF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 1.0D0
            vel2_r(i,j) = 0.7276D0
            vel2_z(i,j) = 0.5D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 0.4D0  
            rho2(i,j) = 0.5313D0 
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.0D0
	END IF
	epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 13) THEN

   ! Test 13, 2D Riemann Problem Test 15 with gamma = 1.4
   gamma2 = 1.4D0
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1
	IF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 0.4D0       
            rho2(i,j) = 0.8D0
            vel2_r(i,j) = 0.1D0
            vel2_z(i,j) = -0.3D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 0.4D0       
            rho2(i,j) = 0.5313D0
            vel2_r(i,j) = 0.1D0
            vel2_z(i,j) = 0.4276D0
	ELSEIF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 0.4D0       
            rho2(i,j) = 0.5197D0
            vel2_r(i,j) = -0.6259D0
            vel2_z(i,j) = -0.3D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.0D0  
            rho2(i,j) = 1.0D0 
            vel2_r(i,j) = 0.1D0
            vel2_z(i,j) = -0.3D0
	END IF
	epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)

ELSEIF(testmodel_flag == 14) THEN

   ! Test 14, 2D Riemann Problem Test 17 with gamma = 1.4
   gamma2 = 1.4D0
   DO j = 1, length_step_z_2, 1
      DO i = 1, length_step_r_2, 1
	IF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 0.4D0       
            rho2(i,j) = 1.0625D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = 0.2145D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.0D0 .AND. z2(j) < 0.5D0) THEN
            p2(i,j) = 0.4D0       
            rho2(i,j) = 0.5197D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = -1.1259D0
	ELSEIF(r2(i) > 0.0D0 .AND. r2(i) < 0.5D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.0D0       
            rho2(i,j) = 2.0D0
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = -0.3D0
	ELSEIF(r2(i) > 0.5D0 .AND. r2(i) < 1.0D0 .AND. z2(j) > 0.5D0 .AND. z2(j) < 1.0D0) THEN
            p2(i,j) = 1.0D0  
            rho2(i,j) = 1.0D0 
            vel2_r(i,j) = 0.0D0
            vel2_z(i,j) = -0.4D0
	END IF
	epsilon2(i,j) = p2(i,j) / rho2(i,j) / (gamma2 - 1.0D0)
      ENDDO
   ENDDO

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(p2, even)
   CALL boundary1D_NM(epsilon2, even)


ENDIF

! In this phase 
! Do a direct copy from motherboard to overlayer
! to test its function, later add the 
! zoomin function

END SUBROUTINE getrho_template