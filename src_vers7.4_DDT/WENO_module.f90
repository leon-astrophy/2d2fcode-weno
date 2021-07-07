!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all necessary code for doing the 
! WENO reconstruction. 
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Merged and systematized by Leung Shing Chi in 2016
! More information about WENO, refer Shu (2000)
!
! This module contains the following subroutines
! 1. subroutine GetConst
! 2. subroutine WENO_r
! 3. subroutine WENO_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE WENO_MODULE        
USE DEFINITION
IMPLICIT NONE

! The C constants
REAL (DP), DIMENSION (-1 : 2, 0 : 2) :: c

! The D and tilde-D constants
REAL (DP), DIMENSION (0 : 2) :: d, td

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine reads in the constant for WENO reconstuction
   ! assuming uniform grid anywhere along the row/column
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine GetConst
   implicit none

   ! Dummy variable
   integer :: r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Const C

   c (-1, 0) = 1.1E1_DP / 6.0E0_DP
   c (-1, 1) = - 7.0E0_DP / 6.0E0_DP
   c (-1, 2) = 1.0E0_DP / 3.0E0_DP

   c (0, 0) = 1.0E0_DP / 3.0E0_DP
   c (0, 1) = 5.0E0_DP / 6.0E0_DP
   c (0, 2) = - 1.0E0_DP / 6.0E0_DP

   c (1, 0) = - 1.0E0_DP / 6.0E0_DP
   c (1, 1) = 5.0E0_DP / 6.0E0_DP
   c (1, 2) = 1.0E0_DP / 3.0E0_DP

   c (2, 0) = 1.0E0_DP / 3.0E0_DP
   c (2, 1) = - 7.0E0_DP / 6.0E0_DP
   c (2, 2) = 1.1E1_DP / 6.0E0_DP

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Const D

   d (0) = 3.0E0_DP / 1.0E1_DP
   d (1) = 3.0E0_DP / 5.0E0_DP
   d (2) = 1.0E0_DP / 1.0E1_DP

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Const tilde_D

   DO r = 0, 2
      td (r) = d (2 - r)
   END DO

   end subroutine getconst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_R
USE RIEMANN_MODULE
!USE FLAME_MODULE
!USE LEVELSET_MODULE
!USE NUCLEAR_MODULE
USE DEFINITION  
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j, in, k

! Dummy !
REAL (DP) :: left, right, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do reconstruction steps for DM !
IF(RUNDM_flag == 1) THEN
	
	! Reconstruct hydrodynamics variables at cell boundaries for DM !
	DO k = length_step_z_min_part_1, length_step_z_part_1
		DO j = 0, length_step_r_part_1 + 1
			CALL WENO (j, p1(j-2,k), p1(j-1,k), p1(j,k), p1(j+1,k), p1(j+2,k), p1R(j-1,k), p1L(j,k))
			CALL WENO (j, rho1(j-2,k), rho1(j-1,k), rho1(j,k), rho1(j+1,k), rho1(j+2,k), rho1R(j-1,k), rho1L(j,k))
			CALL WENO (j, vel1_r(j-2,k), vel1_r(j-1,k), vel1_r(j,k), vel1_r(j+1,k), vel1_r(j+2,k), vel1rR(j-1,k), vel1rL(j,k))
			CALL WENO (j, vel1_z(j-2,k), vel1_z(j-1,k), vel1_z(j,k), vel1_z(j+1,k), vel1_z(j+2,k), vel1zR(j-1,k), vel1zL(j,k))
		END DO
	END DO
	
	! For DM rotation !
	If(rotationdm_flag == 1) THEN
		DO k = length_step_z_min_part_1, length_step_z_part_1
			DO j = 0, length_step_r_part_1 + 1
				CALL WENO (j, vel1_p(j-2,k), vel1_p(j-1,k), vel1_p(j,k), vel1_p(j+1,k), vel1_p(j+2,k), vel1pR(j-1,k), vel1pL(j,k))
			END DO
		END DO
	END IF

	! Do high order interpolation !
	CALL STATESDM_Z

	! Get epsilon at boundary !
	DO k = length_step_z_min_part_1, length_step_z_part_1
		DO j = 0, length_step_r_part_1 + 1
			CALL EOSEPSILON(rho1R(j,k), p1R(j,k), eps1R(j,k), 1)
			CALL EOSEPSILON(rho1L(j,k), p1L(j,k), eps1L(j,k), 1)
		END DO
	END DO	

	! Copy to boundaries !
	CALL BOUNDARYDM_Z (eps1R, even)
	CALL BOUNDARYDM_Z (eps1L, even)

	! No reconstruction for frame velocity since we assume analytic continous form !
	If(movinggriddm_flag == 1) THEN
		DO k = length_step_z_min_part_1, length_step_z_part_1
			DO j = 0, length_step_r_part_1 + 1
				vf1rL(j,k) = vel1_max*rF1(j)/radius1
				vf1rR(j,k) = vel1_max*rF1(j)/radius1
			END DO
		END DO
		CALL BOUNDARYDM_Z (vf1rR, oddR)
		CALL BOUNDARYDM_Z (vf1rL, oddR)
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reconstruct hydrodynamics variables at cell boundaries for NM !
DO k = length_step_z_min_part_2, length_step_z_part_2
	DO j = 0, length_step_r_part_2 + 1
		CALL WENO (j, p2(j-2,k), p2(j-1,k), p2(j,k), p2(j+1,k), p2(j+2,k), p2R(j-1,k), p2L(j,k))
		CALL WENO (j, rho2(j-2,k), rho2(j-1,k), rho2(j,k), rho2(j+1,k), rho2(j+2,k), rho2R(j-1,k), rho2L(j,k))
		CALL WENO (j, vel2_r(j-2,k), vel2_r(j-1,k), vel2_r(j,k), vel2_r(j+1,k), vel2_r(j+2,k), vel2rR(j-1,k), vel2rL(j,k))
		CALL WENO (j, vel2_z(j-2,k), vel2_z(j-1,k), vel2_z(j,k), vel2_z(j+1,k), vel2_z(j+2,k), vel2zR(j-1,k), vel2zL(j,k))
	END DO
END DO
! Rotation !
IF(rotationnm_flag == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 0, length_step_r_part_2 + 1
			CALL WENO (j, vel2_p(j-2,k), vel2_p(j-1,k), vel2_p(j,k), vel2_p(j+1,k), vel2_p(j+2,k), vel2pR(j-1,k), vel2pL(j,k))
		END DO
	END DO
END IF

! Do extra reconstuctions for dual energy !
IF (dual_energy == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 0, length_step_r_part_2 + 1
			CALL WENO (j, rhoe2(j-2,k), rhoe2(j-1,k), rhoe2(j,k), rhoe2(j+1,k), rhoe2(j+2,k), rhoe2R(j-1,k), rhoe2L(j,k))
		END DO
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO k = length_step_z_min_part_2, length_step_z_part_2
			DO j = 0, length_step_r_part_2 + 1
				CALL WENO (j, sca2(j-2,k,i), sca2(j-1,k,i), sca2(j,k,i), sca2(j+1,k,i), sca2(j+2,k,i), sca2R(j-1,k,i), sca2L(j,k,i))
			END DO
		END DO
	END DO
END IF

! Do high order interpolation !
CALL STATESNM_Z

! Get epsilon at boundary !
IF (dual_energy == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 0, length_step_r_part_2 + 1
			eps2R(j,k) = rhoe2R(j,k)/rho2R(j,k)
			eps2L(j,k) = rhoe2L(j,k)/rho2L(j,k)
		END DO
	END DO
ELSE
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 0, length_step_r_part_2 + 1
			CALL EOSEPSILON(rho2R(j,k), p2R(j,k), eps2R(j,k), 2)
			CALL EOSEPSILON(rho2L(j,k), p2L(j,k), eps2L(j,k), 2)
		END DO
	END DO
END IF

! Copy to boundaries !
CALL BOUNDARYNM_Z (eps2R, even)
CALL BOUNDARYNM_Z (eps2L, even)

! No reconstruction for frame velocity since we assume analytic continous form !
If(movinggridnm_flag == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 0, length_step_r_part_2 + 1
			vf2rL(j,k) = vel2_max*rF2(j)/radius2
			vf2rR(j,k) = vel2_max*rF2(j)/radius2
		END DO
	END DO
	CALL BOUNDARYNM_Z (vf2rR, oddR)
	CALL BOUNDARYNM_Z (vf2rL, oddR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation, but along the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO_Z
USE RIEMANN_MODULE
!USE FLAME_MODULE
!USE LEVELSET_MODULE
!USE NUCLEAR_MODULE
USE DEFINITION  
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j, in, k

! Dummy !
REAL (DP) :: left, right, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do reconstruction steps for DM !
IF(RUNDM_flag == 1) THEN
	
	! Reconstruct hydrodynamics variables at cell boundaries for DM !
	DO j = 1, length_step_r_part_1
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_1 + 1
			CALL WENO (k, p1(j,k-2), p1(j,k-1), p1(j,k), p1(j,k+1), p1(j,k+2), p1R(j,k-1), p1L(j,k))
			CALL WENO (k, rho1(j,k-2), rho1(j,k-1), rho1(j,k), rho1(j,k+1), rho1(j,k+2), rho1R(j,k-1), rho1L(j,k))
			CALL WENO (k, vel1_r(j,k-2), vel1_r(j,k-1), vel1_r(j,k), vel1_r(j,k+1), vel1_r(j,k+2), vel1rR(j,k-1), vel1rL(j,k))
			CALL WENO (k, vel1_z(j,k-2), vel1_z(j,k-1), vel1_z(j,k), vel1_z(j,k+1), vel1_z(j,k+2), vel1zR(j,k-1), vel1zL(j,k))
		END DO
	END DO
	
	! For DM rotation !
	If(rotationdm_flag == 1) THEN
		DO j = 1, length_step_r_part_1
			DO k = length_step_z_min_part_1 - 1, length_step_z_part_1 + 1
				CALL WENO (k, vel1_p(j,k-2), vel1_p(j,k-1), vel1_p(j,k), vel1_p(j,k+1), vel1_p(j,k+2), vel1pR(j,k-1), vel1pL(j,k))
			END DO
		END DO
	END IF

	! Do high order interpolation !
	CALL STATESDM_R

	! Get epsilon at boundary !
	DO j = 1, length_step_r_part_1
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_1 + 1
			CALL EOSEPSILON(rho1R(j,k), p1R(j,k), eps1R(j,k), 1)
			CALL EOSEPSILON(rho1L(j,k), p1L(j,k), eps1L(j,k), 1)
		END DO
	END DO	

	! Copy to boundaries !
	CALL BOUNDARYDM_R (eps1R, even)
	CALL BOUNDARYDM_R (eps1L, even)

	! No reconstruction for frame velocity since we assume analytic continous form !
	If(movinggriddm_flag == 1) THEN
		DO j = 1, length_step_r_part_1
			DO k = length_step_z_min_part_1 - 1, length_step_z_part_1 + 1
				vf1zL(j,k) = vel1_max*zF1(k)/radius1
				vf1zR(j,k) = vel1_max*zF1(k)/radius1
			END DO
		END DO
		CALL BOUNDARYDM_R (vf1zR, oddZ)
		CALL BOUNDARYDM_R (vf1zL, oddZ)
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reconstruct hydrodynamics variables at cell boundaries for NM !
DO j = 1, length_step_r_part_2
	DO k = length_step_z_min_part_1 - 1, length_step_z_part_2 + 1
		CALL WENO (k, p2(j,k-2), p2(j,k-1), p2(j,k), p2(j,k+1), p2(j,k+2), p2R(j,k-1), p2L(j,k))
		CALL WENO (k, rho2(j,k-2), rho2(j,k-1), rho2(j,k), rho2(j,k+1), rho2(j,k+2), rho2R(j,k-1), rho2L(j,k))
		CALL WENO (k, vel2_r(j,k-2), vel2_r(j,k-1), vel2_r(j,k), vel2_r(j,k+1), vel2_r(j,k+2), vel2rR(j,k-1), vel2rL(j,k))
		CALL WENO (k, vel2_z(j,k-2), vel2_z(j,k-1), vel2_z(j,k), vel2_z(j,k+1), vel2_z(j,k+2), vel2zR(j,k-1), vel2zL(j,k))
	END DO
END DO

! Rotation !
IF(rotationnm_flag == 1) THEN
	DO j = 1, length_step_r_part_2
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_2 + 1
			CALL WENO (k, vel2_p(j,k-2), vel2_p(j,k-1), vel2_p(j,k), vel2_p(j,k+1), vel2_p(j,k+2), vel2pR(j,k-1), vel2pL(j,k))
		END DO
	END DO
END IF

! Do extra reconstuctions for dual energy !
IF (dual_energy == 1) THEN
	DO j = 1, length_step_r_part_2
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_2 + 1
			CALL WENO (k, rhoe2(j,k-2), rhoe2(j,k-1), rhoe2(j,k), rhoe2(j,k+1), rhoe2(j,k+2), rhoe2R(j,k-1), rhoe2L(j,k))
		END DO
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = 1, length_step_r_part_2
			DO k = length_step_z_min_part_1 - 1, length_step_z_part_2 + 1
				CALL WENO (k, sca2(j,k-2,i), sca2(j,k-1,i), sca2(j,k,i), sca2(j,k+1,i), sca2(j,k+2,i), sca2R(j,k-1,i), sca2L(j,k,i))
			END DO
		END DO
	END DO
END IF

! Do high order interpolation !
CALL STATESNM_R

! Get epsilon at boundary !
IF (dual_energy == 1) THEN
	DO j = 1, length_step_r_part_2
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_2 + 1
			eps2R(j,k) = rhoe2R(j,k)/rho2R(j,k)
			eps2L(j,k) = rhoe2L(j,k)/rho2L(j,k)
		END DO
	END DO
ELSE
	DO j = 1, length_step_r_part_2
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_2 + 1
			CALL EOSEPSILON(rho2R(j,k), p2R(j,k), eps2R(j,k), 2)
			CALL EOSEPSILON(rho2L(j,k), p2L(j,k), eps2L(j,k), 2)
		END DO
	END DO
END IF

! Copy to boundaries !
CALL BOUNDARYNM_R (eps2R, even)
CALL BOUNDARYNM_R (eps2L, even)

! No reconstruction for frame velocity since we assume analytic continous form !
If(movinggridnm_flag == 1) THEN
	DO j = 1, length_step_r_part_2
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2 + 1
			vf2zL(j,k) = vel2_max*zF2(k)/radius2
			vf2zR(j,k) = vel2_max*zF2(k)/radius2
		END DO
	END DO
	CALL BOUNDARYNM_R (vf2zR, oddZ)
	CALL BOUNDARYNM_R (vf2zL, oddZ)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the WENO scheme for reconstructing the numerical flux at both the !
! left and right hand side located at the boundary cell. In this version, I !
! provide different WENO scheme that differ by their smoothness indicator   !
! I also include WENO scheme that use combination of high and low order     !
! polynominal as building block. Nonetheless, a monotonicity preserving     !
! limter option is provided so to make the solution to be MPW               !
! For details, please refer to the textbook with ISBN 3-540-65893-9         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO (i, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: i

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vm_out, vp_out

! Temporal arrays !
REAL (DP), DIMENSION (0 : 2) :: vrhs, vlhs

! For assigning weights !
REAL (DP), DIMENSION (0 : 2) :: alpha, talpha, omega, tomega, beta

! Temporal arrays !
REAL (DP), DIMENSION (i - 2 : i + 2) :: v

! Integer !
INTEGER :: j, r, s

! Tempeorary parameter !
REAL (DP) :: tau, temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign temporal arrays !
v(i - 2) = vm2
v(i - 1) = vm1
v(i) = vc
v(i + 1) = vp1
v(i + 2) = vp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! We calculate the value of u at each grid by the following loop !
! Do the right cell boundary !
DO r = 0, 2
	vrhs (r) = 0.0E0_DP
		
	! We calculate the value of u at right boundary !
	DO j = 0, 2
		vrhs (r) = vrhs (r) + c (r, j) * v (i - r + j)
	END DO
END DO

! Do the left cell boundary !
DO r = 0, 2
	vlhs (r) = 0.0E0_DP

	! Do the same for left boundary !
	DO j = 0, 2
		vlhs (r) = vlhs (r) + c (r - 1, j) * v (i - r + j)
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! These are essential parameter for further construction of u !
beta (0) = (1.3E1_DP / 1.2E1_DP) * (v (i) - 2 * v (i + 1) + v (i + 2)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (3 * v (i) - 4 * v (i + 1) + v (i + 2)) ** 2
beta (1) = (1.3E1_DP / 1.2E1_DP) * (v (i - 1) - 2 * v (i) + v (i + 1)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (v (i - 1) - v (i + 1)) ** 2
beta (2) = (1.3E1_DP / 1.2E1_DP) * (v (i - 2) - 2 * v (i - 1) + v (i)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (v (i - 2) - 4 * v (i - 1) + 3 * v (i)) ** 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assigning tau for the WENO-Z corrections !
tau = abs (beta(0) -beta(2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do WENO-Z weight reconstructions !
DO r = 0, 2
	alpha (r) = d (r) * (1.0D0 + (tau/(beta(r) + smallpara)))
END DO

temp = 0.0E0_DP
	
! The denominator in finding omega, a coefficient for the last step of reconstruction  !
DO s = 0, 2
	temp = temp + alpha (s)
END DO
	
! Find the omega !
DO r = 0, 2
	omega (r) = alpha (r) / temp
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the alpha, omega for the value of u at left grid boundary... !
DO r = 0, 2
	talpha (r) = td (r) * (1.0D0 + (tau/(beta(r) + smallpara)))
END DO

temp = 0.0E0_DP

DO s = 0, 2
	temp = temp + talpha (s)
END DO

DO r = 0, 2
	tomega (r) = talpha (r) / temp
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vp_out = 0.0E0_DP
	
! u at the left boundary !
DO r = 0, 2

	! Original WENO !
	vp_out = vp_out + omega (r) * vrhs (r)

END DO

vm_out = 0.0E0_DP
	
! u at the right boundary !
DO r = 0, 2

	! Original WENO !
	vm_out = vm_out + tomega (r) * vlhs (r)	

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE