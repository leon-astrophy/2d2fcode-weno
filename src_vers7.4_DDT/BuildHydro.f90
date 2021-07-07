!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDHYDRO
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Legendre polynominal !
ALLOCATE (legendre2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, 0 : 2*lmax))

! Mass coordinates !
ALLOCATE (m_r2(-4 : length_step_r_2 + 5))
ALLOCATE (m_cell(-4 : length_step_r_2 + 5))

! For Grid variables !
ALLOCATE (r2(-4 : length_step_r_2 + 5))
ALLOCATE (z2(-4 : length_step_z_2 + 5))
ALLOCATE (zF2(-4 : length_step_z_2 + 5))
ALLOCATE (rF2(-4 : length_step_r_2 + 5))
ALLOCATE (COT_z2(-4 : length_step_z_2 + 5))
ALLOCATE (r2f_m(-4 : length_step_r_2 + 5))  
ALLOCATE (sca2_fac1(-4 : length_step_r_2 + 5))  
ALLOCATE (sca2_fac2(-4 : length_step_r_2 + 5))    
ALLOCATE (rad2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (cos2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (vol2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (radbar2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (volbar2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))

! Hydrodynamic variables !
ALLOCATE (rho2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (vel2_r(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (vel2_z(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (vel2_p(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (epsilon2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (p2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (dpdrho2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (dpdepsilon2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (gamma_2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (ye2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (cs2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))

! Potentials !
ALLOCATE (phi2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (phi2_dm(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (phi2_nm(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (phi2_r(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (phi2_z(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (qpole2(0:2*lmax))

! Back up variables !
ALLOCATE (rho2_old(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (temp2_old(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE (epsilon2_old(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))

! Extra Hydrodynamic variables !
IF (helmeos_flag == 1) THEN
	ALLOCATE (temp2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF
IF (dual_energy == 1) THEN
	ALLOCATE (rhoe2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (bige2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (et2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
        ALLOCATE (dp2dr(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
        ALLOCATE (dp2dz(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF
IF(movinggridnm_flag == 1 .OR. found_movinggridnm_flag == 1) THEN
	ALLOCATE (vel_frame_r2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vel_frame_z2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
IF (DM_flag == 1) THEN
	ALLOCATE (m_r1(-4 : length_step_r_1 + 5))
	ALLOCATE (legendre1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, 0 : 2*lmax))
	ALLOCATE (r1(-4 : length_step_r_1 + 5))
	ALLOCATE (z1(-4 : length_step_z_1 + 5))
	ALLOCATE (zF1(-4 : length_step_r_1 + 5))
	ALLOCATE (rF1(-4 : length_step_r_1 + 5))
	ALLOCATE (COT_z1(-4 : length_step_z_1 + 5))
	ALLOCATE (r1f_m(-4 : length_step_r_1 + 5))  
	ALLOCATE (sca1_fac1(-4 : length_step_r_1 + 5))  
	ALLOCATE (sca1_fac2(-4 : length_step_r_1 + 5))    
	ALLOCATE (rad1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (cos1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (vol1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (rho1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (radbar1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (volbar1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (vel1_r(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (vel1_z(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (vel1_p(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (epsilon1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (p1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (dpdrho1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (dpdepsilon1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (gamma_1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (cs1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (phi1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (phi1_dm(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (phi1_nm(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (phi1_r(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (phi1_z(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (qpole1(0:2*lmax))
	IF(movinggriddm_flag == 1 .OR. found_movinggriddm_flag == 1) THEN
		ALLOCATE (vel_frame_r1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vel_frame_z1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocate arrays related to hydro variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYHYDRO
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Legendre function !
DEALLOCATE (legendre2)

! Mass coordinates !
DEALLOCATE (m_r2)
DEALLOCATE (m_cell)

! For Grid variables !
DEALLOCATE (r2)
DEALLOCATE (z2)
DEALLOCATE (rF2)
DEALLOCATE (zF2)
DEALLOCATE (COT_z2)
DEALLOCATE (r2f_m)
DEALLOCATE (sca2_fac1)
DEALLOCATE (sca2_fac2)
DEALLOCATE (rad2)
DEALLOCATE (cos2)
DEALLOCATE (vol2)
DEALLOCATE (radbar2)
DEALLOCATE (volbar2)

! Hydrodynamic variables !
DEALLOCATE (rho2)
DEALLOCATE (vel2_r)
DEALLOCATE (vel2_z)
DEALLOCATE (vel2_p)
DEALLOCATE (epsilon2)
DEALLOCATE (p2)
DEALLOCATE (dpdrho2)
DEALLOCATE (dpdepsilon2)
DEALLOCATE (gamma_2)
DEALLOCATE (ye2)
DEALLOCATE (cs2)

! Potentials !
DEALLOCATE (phi2)
DEALLOCATE (phi2_dm)
DEALLOCATE (phi2_nm)
DEALLOCATE (phi2_r)
DEALLOCATE (phi2_z)
DEALLOCATE (qpole2)

! Back up variables !
DEALLOCATE (rho2_old)
DEALLOCATE (temp2_old)
DEALLOCATE (epsilon2_old)

! Extra Hydrodynamic variables !
IF (helmeos_flag == 1) THEN
	DEALLOCATE (temp2)
END IF
IF (dual_energy == 1) THEN
	DEALLOCATE (rhoe2)
	DEALLOCATE (bige2)
	DEALLOCATE (et2)
        DEALLOCATE (dp2dr)
        DEALLOCATE (dp2dz)
END IF
IF(movinggridnm_flag == 1 .OR. found_movinggridnm_flag == 1) THEN
	DEALLOCATE (vel_frame_r2)
	DEALLOCATE (vel_frame_z2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
IF (DM_flag == 1) THEN
	DEALLOCATE (legendre1)
	DEALLOCATE (m_r1)
	DEALLOCATE (r1)
	DEALLOCATE (z1)
	DEALLOCATE (rF1)
	DEALLOCATE (zF1)
	DEALLOCATE (COT_z1)
	DEALLOCATE (r1f_m)
	DEALLOCATE (sca1_fac1)
	DEALLOCATE (sca1_fac2)
	DEALLOCATE (rad1)
	DEALLOCATE (cos1)
	DEALLOCATE (vol1)
	DEALLOCATE (rho1)
	DEALLOCATE (radbar1)
	DEALLOCATE (volbar1)
	DEALLOCATE (vel1_r)
	DEALLOCATE (vel1_z)
	DEALLOCATE (vel1_p)
	DEALLOCATE (epsilon1)
	DEALLOCATE (p1)
	DEALLOCATE (dpdrho1)
	DEALLOCATE (dpdepsilon1)
	DEALLOCATE (gamma_1)
	DEALLOCATE (cs1)
	DEALLOCATE (phi1)
	DEALLOCATE (phi1_dm)
	DEALLOCATE (phi1_nm)
	DEALLOCATE (phi1_r)
	DEALLOCATE (phi1_z)
	DEALLOCATE (qpole1)
	IF(movinggriddm_flag == 1 .OR. found_movinggriddm_flag == 1) THEN
		DEALLOCATE (vel_frame_r1)
		DEALLOCATE (vel_frame_z1)
	END IF
END IF

END SUBROUTINE