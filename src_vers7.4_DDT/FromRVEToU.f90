!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine converts the primitive variables to
! conservative variables (or vice versa)
! Written by Leung Shing Chi in 2016
!
! If you add your own variables in the WENO scheme,
! add your conversion step here.
!
! This subroutine takes in the U array and conversion mode
! Mode 0: From primitive to conservative
! Mode 1: From conservative to primitive
!
! Here is a reminder in how to add new physics:
! 1. Add your own module that contains the physics
! 2. Remind BuildWENO to include your quantity
! 3. Add the conversion here
! 4. Write a section in how to calculate the flux in Spatial
! 5. Add a flag to give signal to the program whenever you use the code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FROMRVETOU (u1, u2)
USE DEFINITION
USE turb_module
USE levelset_module
USE helmeos_module
IMPLICIT NONE

! Input U array
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_1 + 5, -4 :length_step_z_1 + 5, imin1 : imax1) :: u1
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 :length_step_z_2 + 5, imin2 : imax2) :: u2

! Dummy variables
INTEGER :: i, j, k

! Varaibles for the SR solver
REAL (DP) :: LorFac
REAL (DP) :: gamfac
REAL (DP) :: rho_out, vel_r_out, vel_z_out, eps_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!if(debug_flag == 1) write(*,*) 'In FromRVEToU', mode_in
!WRITE(*,*) 'Start', length_Step_r_max, length_step_z_max, no_of_eq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM Sector !

! From RVE to U 
IF(runDM_flag == 1) then
	DO k = -4, length_step_z_1 + 5
		DO j = -4, length_step_r_1 + 5
			u1 (j,k,irho1) = rho1(j,k)
	    		u1 (j,k,ivel1_r) = rho1(j,k) * vel1_r(j,k)
	    		u1 (j,k,ivel1_z) = rho1(j,k) * vel1_z(j,k)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    		!u1 (j,k,itau1) = rho1(j,k) * epsilon1(j,k) + 5.0E-1_DP * rho1(j,k) * & 
			!		(vel1_r(j,k) ** 2 + vel1_z(j,k) ** 2 + vel1_p(j,k) ** 2)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         	ENDDO
	ENDDO

	IF (rotationdm_flag == 1) THEN
		DO k = -4, length_step_z_1 + 5
			DO j = -4, length_step_r_1 + 5
				u1 (j,k,ivel1_p) = rho1(j,k) * vel1_p(j,k)
			END DO
		END DO
	END IF
END IF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NM Sector !

! Convert NM hydro
DO k = -4, length_step_z_2 + 5, 1
	DO j = -4, length_step_r_2 + 5, 1
		u2 (j,k,irho2) = rho2(j,k)
		u2 (j,k,ivel2_r) = rho2(j,k) * vel2_r(j,k)
		u2 (j,k,ivel2_z) = rho2(j,k) * vel2_z(j,k)
	END DO
END DO

! Rotation !
IF (rotationnm_flag == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
			u2 (j,k,ivel2_p) = rho2(j,k) * vel2_p(j,k)
		END DO
      	END DO
END IF

! Epsilon !
IF (nm_epsilon == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
            		u2 (j,k,itau2) = rho2(j,k)*epsilon2(j,k) + 5.0E-1_DP*rho2(j,k)*&
				(vel2_r(j,k)**2 + vel2_z(j,k)**2 + vel2_p(j,k)**2) 
		END DO 
	END DO
END IF

! Dual energy !
IF(dual_energy == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
			u2 (j,k,ieps2) = rho2(j,k) * epsilon2(j,k)
		END DO
	END DO
END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !IF(SRhydro_flag == 0) THEN
   !ELSEIF(SRhydro_flag == 1) THEN
   !   DO k = -4, length_step_z_2 + 5, 1
   !      DO j = -4, length_step_r_2 + 5, 1
   !         LorFac = 1.0D0 / SQRT(1.0D0 - (vel2_r(j,k)**2 + vel2_z(j,k)**2))
   !         u2 (j,k,irho2) = rho2(j,k) * LorFac
   !         u2 (j,k,ivel2_r) = rho2(j,k) * (1.0D0 + epsilon2(j,k) + p2(j,k) / rho2(j,k)) * LorFac**2 * vel2_r(j,k)
   !         u2 (j,k,ivel2_z) = rho2(j,k) * (1.0D0 + epsilon2(j,k) + p2(j,k) / rho2(j,k)) * LorFac**2 * vel2_z(j,k)
   !         If (rotationnm_flag == 1) THEN
   !		u2 (j,k,ivel2_p) = rho2(j,k) * (1.0D0 + epsilon2(j,k) + p2(j,k) / rho2(j,k)) * LorFac**2 * vel2_p(j,k)
   !	    END IF
   !         u2 (j,k,itau2) = rho2(j,k) * (1.0D0 + epsilon2(j,k) + p2(j,k) / rho2(j,k)) * LorFac**2 - p2(j,k) - rho2(j,k) * LorFac
   !      ENDDO
   !   ENDDO
   !ELSE
   !   STOP 'Check SRhydro_flag (1 = with SR, 0 = no SR)'
   !ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Other NM scalar equations !

! Convert the electron fraction
IF(etran_flag == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
			u2 (j,k,iye2) = ye2(j,k) * rho2(j,k)
		ENDDO
	ENDDO
ENDIF

   ! Convert the SGS turblence
   IF(turb_flag == 1) THEN
   
       DO k = -4, length_step_z_2 + 5
         DO j = -4, length_step_r_2 + 5
            u2 (j,k,iturbq) = rho2(j,k) * turb_q(j,k)
          ENDDO
      END DO
   ENDIF

   ! Convert the deflagration level set 
   IF(deflevelset_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            u2(j,k,iscaG1) = scaG(j,k) * rho2(j,k)
         ENDDO
      ENDDO
   ENDIF

   ! Convert the detonation level set
   IF(detlevelset_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            u2(j,k,iscaG2) = scaG2(j,k) * rho2(j,k)
         ENDDO
      ENDDO
   ENDIF

   ! Convert the chemical compsoition
   IF(xisotran_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            u2(j,k,ihe4) = rho2(j,k) * xiso(j,k,che4)
   	    u2(j,k,ic12) = rho2(j,k) * xiso(j,k,cc12)
   	    u2(j,k,io16) = rho2(j,k) * xiso(j,k,co16)  
   	    u2(j,k,ine20) = rho2(j,k) * xiso(j,k,cne20)
   	    u2(j,k,img24) = rho2(j,k) * xiso(j,k,cmg24)
   	     u2(j,k,isi28) = rho2(j,k) * xiso(j,k,csi28)
   	    u2(j,k,ini56) = rho2(j,k) * xiso(j,k,cni56)	    
         ENDDO
      ENDDO
   ENDIF

   IF(burn_prog_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1      
         DO j = -4, length_step_r_2 + 5, 1         
            u2(j,k,ibphi1) = rho2(j,k) * burn_phi1(j,k)
            u2(j,k,ibphi2) = rho2(j,k) * burn_phi2(j,k)
            u2(j,k,ibphi3) = rho2(j,k) * burn_phi3(j,k)
   	    u2(j,k,iyiso) = rho2(j,k) * yiso(j,k)
   	    u2(j,k,iqash) = rho2(j,k) * qash(j,k)
         ENDDO
      ENDDO
   ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert back to primitive variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMUTORVE (u1, u2)
USE DEFINITION
USE turb_module
USE levelset_module
USE helmeos_module
IMPLICIT NONE

! Input U array
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_1 + 5, -4 :length_step_z_1 + 5, imin1 : imax1) :: u1
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_2 + 5, -4 :length_step_z_2 + 5, imin2 : imax2) :: u2

! Dummy variables
INTEGER :: i, j, k

! Dummy !
REAL (DP) :: dummy1, dummy2, dummy_max

! Varaibles for the SR solver
REAL (DP) :: LorFac
REAL (DP) :: gamfac
REAL (DP) :: rho_out, vel_r_out, vel_z_out, eps_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From U to RVE for DM sectors !

! Convert the DM hydro
IF(runDM_flag == 1) THEN
	DO k = -4, length_step_z_1 + 5, 1
		DO j = -4, length_step_r_1 + 5, 1
            		rho1(j,k) = u1 (j,k,irho1)
            		vel1_r(j,k) = u1 (j,k,ivel1_r) / rho1(j,k)
            		vel1_z(j,k) = u1 (j,k,ivel1_z) / rho1(j,k)
 	    		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            		!epsilon1(j,k) = u1 (j,k,itau1) / rho1(j,k) - 5.0E-1_DP * & 
	    		!	(vel1_r(j,k) ** 2 + vel1_z(j,k) ** 2 + vel1_p(j,k) ** 2)
            		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ENDDO
	ENDDO
      
	! Rotation for DM !
	IF (rotationdm_flag == 1) THEN
		DO k = -4, length_step_z_1 + 5, 1
			DO j = -4, length_step_r_1 + 5, 1
				vel1_p(j,k) = u1 (j,k,ivel1_p) / rho1(j,k)
         		ENDDO
      		ENDDO
	END IF

	! Copy to boundaries !
      	CALL BOUNDARY1D_DMFULL (rho1,even)
      	CALL BOUNDARY1D_DMFULL (vel1_r,oddR)
      	CALL BOUNDARY1D_DMFULL (vel1_z,oddZ)
      	CALL BOUNDARY1D_DMFULL (vel1_p,oddR)
      	CALL BOUNDARY1D_DMFULL (epsilon1,even)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

! Convert the NM hydro
DO k = -4, length_step_z_2 + 5, 1
	DO j = -4, length_step_r_2 + 5, 1
		rho2(j,k) = u2 (j,k,irho2)
            	vel2_r(j,k) = u2 (j,k,ivel2_r) / rho2(j,k)
            	vel2_z(j,k) = u2 (j,k,ivel2_z) / rho2(j,k)
	ENDDO      
ENDDO

! For NM rotation !
IF (rotationnm_flag == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
			vel2_p(j,k) = u2 (j,k,ivel2_p) / rho2(j,k)
		END DO
	END DO
END IF

! Boundaries !
CALL BOUNDARY1D_NMFULL (rho2,even)
CALL BOUNDARY1D_NMFULL (vel2_r,oddR)
CALL BOUNDARY1D_NMFULL (vel2_z,oddZ)
CALL BOUNDARY1D_NMFULL (vel2_p,oddR)

! Dual energies !
IF (dual_energy == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
			rhoe2 (j,k) = u2 (j, k, ieps2)
			bige2 (j,k) = u2 (j,k,itau2) / rho2 (j,k)
			et2 (j,k) = (u2 (j,k,itau2) - 5.0E-1_DP * rho2(j,k) * (vel2_r(j,k) ** 2 & 
				  + vel2_z(j,k) ** 2 + vel2_p(j,k) ** 2)) / rho2(j,k)
		END DO
	END DO
	CALL BOUNDARY1D_NMFULL (rhoe2,even)
	CALL BOUNDARY1D_NMFULL (bige2,even)
	CALL BOUNDARY1D_NMFULL (et2,even)

	! Determine the epsilon for epsilon equation !
	DO k = 1, length_step_z_2, 1
		DO j = 1, length_step_r_2, 1
			If(et2 (j,k) > 1.0D-1*bige2 (j,k)) THEN
				epsilon2(j,k) = et2 (j,k)
			ELSE
				epsilon2(j,k) = u2 (j,k,ieps2) / rho2 (j,k)
			END IF
		END DO
	END DO

	! Special setting for boundary values !
	CALL BOUNDARY1D_NMFULL (epsilon2, even)
ELSE

	! For NM epsilon equation !
	IF (nm_epsilon == 1) THEN
		DO k = -4, length_step_z_2 + 5, 1
			DO j = -4, length_step_r_2 + 5, 1
				epsilon2(j,k) = (u2 (j,k,itau2) - 5.0E-1_DP * rho2(j,k) * &
					(vel2_r(j,k) ** 2 + vel2_z(j,k) ** 2 + vel2_p(j,k) ** 2)) / rho2(j,k)
			END DO
		END DO
		CALL BOUNDARY1D_NMFULL (epsilon2,even)
	END IF

END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !IF(SRhydro_flag == 0) THEN
   !ELSEIF(SRhydro_flag == 1) THEN
   !   IF(helmeos_flag == 1) THEN
   !      DO k = 1, length_step_z_2
   !         DO j = 1, length_step_r_2 
   !    
   !            ! For general EOS
   !            CALL solveRDTHelm(j, k, u2(j,k,irho2), u2(j,k,ivel2_r), u2(j,k,ivel2_z), u2(j,k,itau2), &
   !                              rho2(j,k), vel2_r(j,k), vel2_z(j,k), epsilon2(j,k), temp2(j,k), &
   !                              abar2(j,k), zbar2(j,k), ye2(j,k), &
   !                              rho_out, vel_r_out, vel_z_out, eps_out)
   !            rho2(j,k) = rho_out
   !            vel2_r(j,k) = vel_r_out
   !            vel2_z(j,k) = vel_z_out
   !            vel2_p(j,k) = 0.0D0
   !            epsilon2(j,k) = eps_out
   !         ENDDO      
   !      ENDDO
   !    
   !      ! Boundaries !
   !   	 CALL BOUNDARY1D_NMFULL (rho2,EVEN)
   !      CALL BOUNDARY1D_NMFULL (vel2_r,oddR)
   !  	 CALL BOUNDARY1D_NMFULL (vel2_z,oddZ)
   !   	 CALL BOUNDARY1D_NMFULL (vel2_p,oddR)
   !   	 CALL BOUNDARY1D_NMFULL (epsilon2,EVEN)
   !    	 CALL BOUNDARY1D_NMFULL (rhoe2,EVEN)
   !   ELSEIF(polyeosnm_flag == 1) THEN
   !	 DO k = 1, length_step_z_2  
   !         DO j = 1, length_step_r_2 
   !
   !            ! For general EOS
   !            CALL solveRDTPoly(j, k, u2(j,k,irho2), u2(j,k,ivel2_r), u2(j,k,ivel2_z), u2(j,k,itau2), &
   !                              rho2(j,k), vel2_r(j,k), vel2_z(j,k), epsilon2(j,k), &
   !                              rho_out, vel_r_out, vel_z_out, eps_out)
   !            rho2(j,k) = rho_out
   !            vel2_r(j,k) = vel_r_out
   !            vel2_z(j,k) = vel_z_out
   !            vel2_p(j,k) = 0.0D0 
   !            epsilon2(j,k) = eps_out
   !
   !         ENDDO
   !      ENDDO
   !
   !      ! Boundaries !
   !	 CALL BOUNDARY1D_NMFULL (rho2,EVEN)
   !      CALL BOUNDARY1D_NMFULL (vel2_r,oddR)
   !	 CALL BOUNDARY1D_NMFULL (vel2_z,oddZ)
   !	 CALL BOUNDARY1D_NMFULL (vel2_p,oddR)
   !	 CALL BOUNDARY1D_NMFULL (epsilon2,EVEN) 
   !	 CALL BOUNDARY1D_NMFULL (rhoe2,EVEN)
   !   ENDIF
   !ELSE
   !   STOP 'Check SRhydro_flag (1 = with SR, 0 = no SR)'
   !ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert the electron fraction
IF(etran_flag == 1) THEN
	DO k = -4, length_step_z_2 + 5, 1
		DO j = -4, length_step_r_2 + 5, 1
			ye2(j,k) = u2 (j,k,iye2) / rho2(j,k) 
	    		! Max/Min check
	    		ye2(j,k) = MIN(MAX(ye2(j,k), ye_min), ye_max)
         	ENDDO
      	ENDDO
ENDIF

   ! Convert the SGS turbulence
   IF(turb_flag == 1) then
      DO k = -4, length_step_z_2 + 5
         DO j = -4, length_step_r_2 + 5
            turb_q(j,k) = u2(j,k,iturbq) / rho2(j,k)
   	    ! Max/Min check
   	    turb_q(j,k) = MIN(MAX(turb_q(j,k), turb_q_a), 1.0D-5)
         ENDDO   
      ENDDO
   
   ENDIF

   ! Convert the deflagration level set
   IF(deflevelset_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            scaG(j,k) = u2(j,k,iscaG1) / rho2(j,k)
         ENDDO
      ENDDO
   ENDIF

   ! Convert the detonation level set
   IF(detlevelset_flag == 1) then
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            scaG2(j,k) = u2(j,k,iscaG2) / rho2(j,k)
        ENDDO
      ENDDO
   ENDIF

   ! Convert the chemical composition
   IF(xisotran_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            xiso(j,k,che4) = u2(j,k,ihe4) / rho2(j,k)
   	    xiso(j,k,cc12) = u2(j,k,ic12) / rho2(j,k)
   	    xiso(j,k,co16) = u2(j,k,io16) / rho2(j,k)
   	    xiso(j,k,cne20) = u2(j,k,ine20) / rho2(j,k)
   	    xiso(j,k,cmg24) = u2(j,k,img24) / rho2(j,k)
   	    xiso(j,k,csi28) = u2(j,k,isi28) / rho2(j,k)
   	    xiso(j,k,cni56) = u2(j,k,ini56) / rho2(j,k)
         ENDDO
      ENDDO
   
      ! Check the chemical composition
      CALL checkxisotope
   ENDIF

   IF(burn_prog_flag == 1) THEN
      DO k = -4, length_step_z_2 + 5, 1
         DO j = -4, length_step_r_2 + 5, 1
            !burn_phi1(j,k) = MAX(u(ibphi1,j,k) / rho2(j,k), burn_phi1(j,k))
            !burn_phi2(j,k) = MAX(u(ibphi2,j,k) / rho2(j,k), burn_phi2(j,k))
            !burn_phi3(j,k) = MAX(u(ibphi3,j,k) / rho2(j,k), burn_phi3(j,k))
   	    burn_phi1(j,k) = u2(j,k,ibphi1) / rho2(j,k)
            burn_phi2(j,k) = u2(j,k,ibphi2) / rho2(j,k)
            burn_phi3(j,k) = u2(j,k,ibphi3) / rho2(j,k)
   	    yiso(j,k) = u2(j,k,iyiso) / rho2(j,k)
   	    qash(j,k) = u2(j,k,iqash) / rho2(j,k)
         ENDDO
      ENDDO 
   
      ! Check the results 
      CALL check_Burnphi
      CALL ReconstructXIso
   END IF

END SUBROUTINE