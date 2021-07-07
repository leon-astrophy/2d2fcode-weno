!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Error in small corrections !
REAL (DP), PARAMETER :: corr1 = (1.0D0/(max(length_step_r_1, length_step_z_1)))**5
REAL (DP), PARAMETER :: corr2 = (1.0D0/(max(length_step_r_2, length_step_z_2)))**5

! the alpha in the LF flux !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: alpha1_r, alpha1_z
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: alpha2_r, alpha2_z

! Left and right hydro-states for DM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: rho1R, rho1L, p1R, p1L, eps1R, eps1L 
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vel1rR, vel1rL, vel1zR, vel1zL, vel1pR, vel1pL

! DM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf1rR, vf1rL
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf1zR, vf1zL

! Left and right hydro-states for NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: rho2R, rho2L, p2R, p2L, eps2R, eps2L, rhoe2L, rhoe2R
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vel2rR, vel2rL, vel2zR, vel2zL, vel2pR, vel2pL

! NM moving grid !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf2rR, vf2rL
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: vf2zR, vf2zL

! Speed of sound !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: cs1L, cs1R
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: cs2L, cs2R

! Left and right scalars for NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: sca2L, sca2R

! Left and right fluxes, conserved quantity !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: fluxL1, fluxR1, uL1, uR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: fluxL2, fluxR2, uL2, uR2

! Scalars variables !
REAL (DP), ALLOCATABLE, DIMENSION(:,:,:) :: sca2

! Grid flag !
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: gridflag1
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: gridflag2

! For high order interpolation !
REAL (DP), PARAMETER :: qm2 = -4.6875D-3
REAL (DP), PARAMETER :: qm1 = (29.0D0/480.0D0)
REAL (DP), PARAMETER :: qc = -(107.0D0/960.0D0)
REAL (DP), PARAMETER :: qp1 = (29.0D0/480.0D0)
REAL (DP), PARAMETER :: qp2 = -4.6875D-3

! For high order interpolation !
REAL (DP), PARAMETER :: fm2 = (17.0D0/5760.0D0)
REAL (DP), PARAMETER :: fm1 = -(77.0D0/1440.0D0)
REAL (DP), PARAMETER :: fc = (97.0D0/960.0D0)
REAL (DP), PARAMETER :: fp1 = -(77.0D0/1440.0D0)
REAL (DP), PARAMETER :: fp2 = (17.0D0/5760.0D0)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Arrays for alpha splits !
ALLOCATE (alpha1_r(-4 : length_step_z_1 + 5, imin1:imax1))
ALLOCATE (alpha1_z(-4 : length_step_r_1 + 5, imin1:imax1))
ALLOCATE (alpha2_r(-4 : length_step_z_2 + 5, imin2:imax2))
ALLOCATE (alpha2_z(-4 : length_step_r_2 + 5, imin2:imax2))

! For scalars !
ALLOCATE (sca2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, iminsca2 : imaxsca2))

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	ALLOCATE(rho1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(rho1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(p1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(p1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(eps1L(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(eps1R(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(vel1rL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(vel1rR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(vel1zL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(vel1zR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(vel1pL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(vel1pR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(gridflag1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE(fluxL1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(fluxR1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(uL1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	ALLOCATE(uR1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1))
	IF(movinggriddm_flag == 1) THEN
		ALLOCATE (vf1rL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vf1rR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vf1zL(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
		ALLOCATE (vf1zR(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	END IF
END IF

! NM !
ALLOCATE(rho2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(rho2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(p2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(p2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(eps2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(eps2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(vel2rL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(vel2rR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(vel2zL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(vel2zR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(vel2pL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(vel2pR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(gridflag2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
ALLOCATE(fluxL2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(fluxR2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(uL2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(uR2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2))
ALLOCATE(sca2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, iminsca2 : imaxsca2))
ALLOCATE(sca2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, iminsca2 : imaxsca2))

! Dual energy !
IF(dual_energy == 1) THEN
	ALLOCATE(rhoe2L(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE(rhoe2R(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF

! Moving grid !
IF(movinggridnm_flag == 1) THEN
	ALLOCATE (vf2rL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vf2rR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vf2zL(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
	ALLOCATE (vf2zR(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CLEARRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	rho1L = 0.0D0
	rho1R = 0.0D0
	p1L = 0.0D0
	p1R = 0.0D0
	eps1L = 0.0D0
	eps1R = 0.0D0
	vel1rL = 0.0D0
	vel1rR = 0.0D0
	vel1zL = 0.0D0
	vel1zR = 0.0D0
	gridflag1 = 0.0D0
	fluxL1 = 0.0D0
	fluxR1 = 0.0D0
	uL1 = 0.0D0
	uR1 = 0.0D0
END IF

! NM !
rho2L = 0.0D0
rho2R = 0.0D0
p2L = 0.0D0
p2R = 0.0D0
eps2L = 0.0D0
eps2R = 0.0D0
vel2rL = 0.0D0
vel2rR = 0.0D0
vel2zL = 0.0D0
vel2zR = 0.0D0
gridflag2 = 0.0D0
fluxL2 = 0.0D0
fluxR2 = 0.0D0
uL2 = 0.0D0
uR2 = 0.0D0
sca2L= 0.0D0
sca2R = 0.0D0

! Dual energy !
IF(dual_energy == 1) THEN
	rhoe2L = 0.0D0
	rhoe2R = 0.0D0
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Arrays for alpha splits !
DEALLOCATE (alpha1_r)
DEALLOCATE (alpha1_z)
DEALLOCATE (alpha2_r)
DEALLOCATE (alpha2_z)

! For scalars !
DEALLOCATE (sca2)

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	DEALLOCATE(rho1L)
	DEALLOCATE(rho1R)
	DEALLOCATE(p1L)
	DEALLOCATE(p1R)
	DEALLOCATE(eps1L)
	DEALLOCATE(eps1R)
	DEALLOCATE(vel1rL)
	DEALLOCATE(vel1rR)
	DEALLOCATE(vel1zL)
	DEALLOCATE(vel1zR)
	DEALLOCATE(vel1pL)
	DEALLOCATE(vel1pR)
	DEALLOCATE(gridflag1)
	DEALLOCATE(fluxL1)
	DEALLOCATE(fluxR1)
	DEALLOCATE(uL1)
	DEALLOCATE(uR1)
	IF(movinggriddm_flag == 1) THEN
		DEALLOCATE (vf1rL)
		DEALLOCATE (vf1rR)
		DEALLOCATE (vf1zL)
		DEALLOCATE (vf1zR)
	END IF
END IF

! NM !
DEALLOCATE(rho2L)
DEALLOCATE(rho2R)
DEALLOCATE(p2L)
DEALLOCATE(p2R)
DEALLOCATE(eps2L)
DEALLOCATE(eps2R)
DEALLOCATE(vel2rL)
DEALLOCATE(vel2rR)
DEALLOCATE(vel2zL)
DEALLOCATE(vel2zR)
DEALLOCATE(vel2pL)
DEALLOCATE(vel2pR)
DEALLOCATE(gridflag2)
DEALLOCATE(fluxL2)
DEALLOCATE(fluxR2)
DEALLOCATE(uL2)
DEALLOCATE(uR2)
DEALLOCATE(sca2L)
DEALLOCATE(sca2R)

! Dual energy !
IF(dual_energy == 1) THEN
	DEALLOCATE(rhoe2L)
	DEALLOCATE(rhoe2R)
END IF

! Movinggrid !
IF(movinggridnm_flag == 1) THEN
	DEALLOCATE (vf2rL)
	DEALLOCATE (vf2rR)
	DEALLOCATE (vf2zL)
	DEALLOCATE (vf2zR)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM_R (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO i = imin1, imax1
	DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
		DO j = 0, length_step_r_part_1
			flux_out(j,k,i) = 0.5D0 * (fluxL1 (j,k,i) + fluxR1 (j,k,i) - alpha1_r(k,i) * (uR1 (j,k,i) - uL1 (j,k,i)))
		END DO
	END DO
END DO

! Do high order interpolation !
DO i = imin1, imax1
	temporal(:,:) = flux_out(:,:,i)
	DO j = 0, length_step_r_part_1
		DO k = length_step_z_min_part_1, length_step_z_part_1
			IF(gridflag1(j,k) == 1) THEN
				flux_out(j,k,i) = temporal(j,k)
			ELSE
				flux_out(j,k,i) = temporal(j,k) - (fm2*temporal(j,k-2) + fm1*temporal(j,k-1) + & 
		      	                  fc*temporal(j,k) + fp1*temporal(j,k+1) + fp2*temporal(j,k+2))
			END IF
		END DO
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFDM_Z (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For DM !
DO i = imin1, imax1
	DO j = -4, length_step_r_part_1 + 5
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
			flux_out(j,k,i) = 0.5D0 * (fluxL1 (j,k,i) + fluxR1 (j,k,i) - alpha1_z(j,i) * (uR1 (j,k,i) - uL1 (j,k,i)))
		END DO
	END DO
END DO

! Do high order interpolation !
DO i = imin1, imax1
	temporal(:,:) = flux_out(:,:,i)
	DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
		DO j = 1, length_step_r_part_1
			IF(gridflag1(j,k) == 1) THEN
				flux_out(j,k,i) = temporal(j,k)
			ELSE
				flux_out(j,k,i) = temporal(j,k) - (fm2*temporal(j-2,k) + fm1*temporal(j-1,k) + & 
			      	                  fc*temporal(j,k) + fp1*temporal(j+1,k) + fp2*temporal(j+2,k))
			END IF
		END DO
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM_R (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
DO i = imin2, imax2
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			flux_out(j,k,i) = 0.5D0 * (fluxL2 (j,k,i) + fluxR2 (j,k,i) - alpha2_r(k,i) * (uR2 (j,k,i) - uL2 (j,k,i)))
		END DO
	END DO
END DO

! Do high order interpolation !
DO i = imin2, imax2
	temporal(:,:) = flux_out(:,:,i)
	DO j = 0, length_step_r_part_2
		DO k = length_step_z_min_part_2, length_step_z_part_2
			IF(gridflag2(j,k) == 1) THEN
				flux_out(j,k,i) = temporal(j,k)
			ELSE
				flux_out(j,k,i) = temporal(j,k) - (fm2*temporal(j,k-2) + fm1*temporal(j,k-1) + & 
			      	                  fc*temporal(j,k) + fp1*temporal(j,k+1) + fp2*temporal(j,k+2))
			END IF
		END DO
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LFNM_Z (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! output !
REAL (DP), INTENT (OUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: flux_out

! Temporal arrays !
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: temporal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
DO i = imin2, imax2
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			flux_out(j,k,i) = 0.5D0 * (fluxL2 (j,k,i) + fluxR2 (j,k,i) - alpha2_z(j,i) * (uR2 (j,k,i) - uL2 (j,k,i)))
		END DO
	END DO
END DO

! Do high order interpolation !
DO i = imin2, imax2
	temporal(:,:) = flux_out(:,:,i)
	DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
		DO j = 1, length_step_r_part_2
			IF(gridflag2(j,k) == 1) THEN
				flux_out(j,k,i) = temporal(j,k)
			ELSE
				flux_out(j,k,i) = temporal(j,k) - (fm2*temporal(j-2,k) + fm1*temporal(j-1,k) + & 
			      	                  fc*temporal(j,k) + fp1*temporal(j+1,k) + fp2*temporal(j+2,k))
			END IF
		END DO
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE to find the maximum effective speed !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDALPHA
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Dummy variables !
REAL (DP), DIMENSION(imin1:imax1) :: dummy1
REAL (DP), DIMENSION(imin2:imax2) :: dummy2

! Do for DM !
IF(RUNDM_flag == 1) THEN
	DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
		CALL ALPHASPLIT_R (alpha1_r(k,:), dummy2(:), k, 1)
	END DO	
	DO j = -4, length_step_r_part_1 + 5
		CALL ALPHASPLIT_Z (alpha1_z(j,:), dummy2(:), j, 1)
	END DO
END IF

! Do for NM !
DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
	CALL ALPHASPLIT_R (dummy1(:), alpha2_r(k,:), k, 2)
END DO	
DO j = -4, length_step_r_part_2 + 5
	CALL ALPHASPLIT_Z (dummy1(:), alpha2_z(j,:), j, 2)
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE