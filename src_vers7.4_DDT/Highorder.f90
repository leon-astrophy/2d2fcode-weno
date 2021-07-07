!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the z-directions for DM states
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE STATESDM_Z
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Copy to boundaries !
CALL BOUNDARYDM_Z (p1R, even)
CALL BOUNDARYDM_Z (p1L, even)
CALL BOUNDARYDM_Z (rho1R, even)
CALL BOUNDARYDM_Z (rho1L, even)
CALL BOUNDARYDM_Z (vel1rR, oddR)
CALL BOUNDARYDM_Z (vel1rL, oddR)
CALL BOUNDARYDM_Z (vel1zR, oddZ)
CALL BOUNDARYDM_Z (vel1zL, oddZ)
IF(rotationdm_flag == 1) THEN
	CALL BOUNDARYDM_Z (vel1pR, oddR)
	CALL BOUNDARYDM_Z (vel1pL, oddR)
END IF

! Check the grid flag !
CALL CHECKGRIDDM_Z (p1R)
CALL CHECKGRIDDM_Z (p1L)
CALL CHECKGRIDDM_Z (rho1R)
CALL CHECKGRIDDM_Z (rho1L)
CALL CHECKGRIDDM_Z (vel1rR)
CALL CHECKGRIDDM_Z (vel1rL)
CALL CHECKGRIDDM_Z (vel1zR)
CALL CHECKGRIDDM_Z (vel1zL)
IF(rotationdm_flag == 1) THEN
	CALL CHECKGRIDDM_Z (vel1pR)
	CALL CHECKGRIDDM_Z (vel1pL)
END IF

! Do high order interpolation !
CALL HIGHORDERDM_Z (p1R)
CALL HIGHORDERDM_Z (p1L)
CALL HIGHORDERDM_Z (rho1R)
CALL HIGHORDERDM_Z (rho1L)
CALL HIGHORDERDM_Z (vel1rR)
CALL HIGHORDERDM_Z (vel1rL)
CALL HIGHORDERDM_Z (vel1zR)
CALL HIGHORDERDM_Z (vel1zL)
IF(rotationdm_flag == 1) THEN
	CALL HIGHORDERDM_Z (vel1pR)
	CALL HIGHORDERDM_Z (vel1pL)
END IF

! Copy to boundary again !
CALL BOUNDARYDM_Z (p1R, even)
CALL BOUNDARYDM_Z (p1L, even)
CALL BOUNDARYDM_Z (rho1R, even)
CALL BOUNDARYDM_Z (rho1L, even)
CALL BOUNDARYDM_Z (vel1rR, oddR)
CALL BOUNDARYDM_Z (vel1rL, oddR)
CALL BOUNDARYDM_Z (vel1zR, oddZ)
CALL BOUNDARYDM_Z (vel1zL, oddZ)
IF(rotationdm_flag == 1) THEN
	CALL BOUNDARYDM_Z (vel1pR, oddR)
	CALL BOUNDARYDM_Z (vel1pL, oddR)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the z-directions for NM states
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE STATESNM_Z
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Copy to boundaries !
CALL BOUNDARYNM_Z (p2R, even)
CALL BOUNDARYNM_Z (p2L, even)
CALL BOUNDARYNM_Z (rho2R, even)
CALL BOUNDARYNM_Z (rho2L, even)
CALL BOUNDARYNM_Z (vel2rR, oddR)
CALL BOUNDARYNM_Z (vel2rL, oddR)
CALL BOUNDARYNM_Z (vel2zR, oddZ)
CALL BOUNDARYNM_Z (vel2zL, oddZ)
IF(rotationnm_flag == 1) THEN
	CALL BOUNDARYNM_Z (vel2pR, oddR)
	CALL BOUNDARYNM_Z (vel2pL, oddR)
END IF
IF (dual_energy == 1) THEN
	CALL BOUNDARYNM_Z (rhoe2R, even)
	CALL BOUNDARYNM_Z (rhoe2L, even)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL BOUNDARYNM_Z (sca2R(:,:,i), even)
		CALL BOUNDARYNM_Z (sca2L(:,:,i), even)
	END DO
END IF

! Check the grid !
CALL CHECKGRIDNM_Z (p2R)
CALL CHECKGRIDNM_Z (p2L)
CALL CHECKGRIDNM_Z (rho2R)
CALL CHECKGRIDNM_Z (rho2L)
CALL CHECKGRIDNM_Z (vel2rR)
CALL CHECKGRIDNM_Z (vel2rL)
CALL CHECKGRIDNM_Z (vel2zR)
CALL CHECKGRIDNM_Z (vel2zL)
IF(rotationnm_flag == 1) THEN
	CALL CHECKGRIDNM_Z (vel2pR)
	CALL CHECKGRIDNM_Z (vel2pL)
END IF
IF (dual_energy == 1) THEN
	CALL CHECKGRIDNM_Z (rhoe2R)
	CALL CHECKGRIDNM_Z (rhoe2L)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL CHECKGRIDNM_Z (sca2R(:,:,i))
		CALL CHECKGRIDNM_Z (sca2L(:,:,i))
	END DO
END IF

! Do high order interpolation !
CALL HIGHORDERNM_Z (p2R)
CALL HIGHORDERNM_Z (p2L)
CALL HIGHORDERNM_Z (rho2R)
CALL HIGHORDERNM_Z (rho2L)
CALL HIGHORDERNM_Z (vel2rR)
CALL HIGHORDERNM_Z (vel2rL)
CALL HIGHORDERNM_Z (vel2zR)
CALL HIGHORDERNM_Z (vel2zL)
IF(rotationnm_flag == 1) THEN
	CALL HIGHORDERNM_Z (vel2pR)
	CALL HIGHORDERNM_Z (vel2pL)
END IF
IF (dual_energy == 1) THEN
	CALL HIGHORDERNM_Z (rhoe2R)
	CALL HIGHORDERNM_Z (rhoe2L)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL HIGHORDERNM_Z (sca2R(:,:,i))
		CALL HIGHORDERNM_Z (sca2L(:,:,i))
	END DO
END IF

! Copy to boundaries !
CALL BOUNDARYNM_Z (p2R, even)
CALL BOUNDARYNM_Z (p2L, even)
CALL BOUNDARYNM_Z (rho2R, even)
CALL BOUNDARYNM_Z (rho2L, even)
CALL BOUNDARYNM_Z (vel2rR, oddR)
CALL BOUNDARYNM_Z (vel2rL, oddR)
CALL BOUNDARYNM_Z (vel2zR, oddZ)
CALL BOUNDARYNM_Z (vel2zL, oddZ)
IF(rotationnm_flag == 1) THEN
	CALL BOUNDARYNM_Z (vel2pR, oddR)
	CALL BOUNDARYNM_Z (vel2pL, oddR)
END IF
IF (dual_energy == 1) THEN
	CALL BOUNDARYNM_Z (rhoe2R, even)
	CALL BOUNDARYNM_Z (rhoe2L, even)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL BOUNDARYNM_Z (sca2R(:,:,i), even)
		CALL BOUNDARYNM_Z (sca2L(:,:,i), even)
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the z-directions for states
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE STATESDM_R
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Copy to boundaries !
CALL BOUNDARYDM_R (p1R, even)
CALL BOUNDARYDM_R (p1L, even)
CALL BOUNDARYDM_R (rho1R, even)
CALL BOUNDARYDM_R (rho1L, even)
CALL BOUNDARYDM_R (vel1rR, oddR)
CALL BOUNDARYDM_R (vel1rL, oddR)
CALL BOUNDARYDM_R (vel1zR, oddZ)
CALL BOUNDARYDM_R (vel1zL, oddZ)
IF(rotationdm_flag == 1) THEN
	CALL BOUNDARYDM_R (vel1pR, oddR)
	CALL BOUNDARYDM_R (vel1pL, oddR)
END IF

! Check the grid flag !
CALL CHECKGRIDDM_R (p1R)
CALL CHECKGRIDDM_R (p1L)
CALL CHECKGRIDDM_R (rho1R)
CALL CHECKGRIDDM_R (rho1L)
CALL CHECKGRIDDM_R (vel1rR)
CALL CHECKGRIDDM_R (vel1rL)
CALL CHECKGRIDDM_R (vel1zR)
CALL CHECKGRIDDM_R (vel1zL)
IF(rotationdm_flag == 1) THEN
	CALL CHECKGRIDDM_R (vel1pR)
	CALL CHECKGRIDDM_R (vel1pL)
END IF

! Do high order interpolation !
CALL HIGHORDERDM_R (p1R)
CALL HIGHORDERDM_R (p1L)
CALL HIGHORDERDM_R (rho1R)
CALL HIGHORDERDM_R (rho1L)
CALL HIGHORDERDM_R (vel1rR)
CALL HIGHORDERDM_R (vel1rL)
CALL HIGHORDERDM_R (vel1zR)
CALL HIGHORDERDM_R (vel1zL)
IF(rotationdm_flag == 1) THEN
	CALL HIGHORDERDM_R (vel1pR)
	CALL HIGHORDERDM_R (vel1pL)
END IF

! Copy to boundary again !
CALL BOUNDARYDM_R (p1R, even)
CALL BOUNDARYDM_R (p1L, even)
CALL BOUNDARYDM_R (rho1R, even)
CALL BOUNDARYDM_R (rho1L, even)
CALL BOUNDARYDM_R (vel1rR, oddR)
CALL BOUNDARYDM_R (vel1rL, oddR)
CALL BOUNDARYDM_R (vel1zR, oddZ)
CALL BOUNDARYDM_R (vel1zL, oddZ)
IF(rotationdm_flag == 1) THEN
	CALL BOUNDARYDM_R (vel1pR, oddR)
	CALL BOUNDARYDM_R (vel1pL, oddR)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the z-directions for states
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE STATESNM_R
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Copy to boundaries !
CALL BOUNDARYNM_R (p2R, even)
CALL BOUNDARYNM_R (p2L, even)
CALL BOUNDARYNM_R (rho2R, even)
CALL BOUNDARYNM_R (rho2L, even)
CALL BOUNDARYNM_R (vel2rR, oddR)
CALL BOUNDARYNM_R (vel2rL, oddR)
CALL BOUNDARYNM_R (vel2zR, oddZ)
CALL BOUNDARYNM_R (vel2zL, oddZ)
IF(rotationnm_flag == 1) THEN
	CALL BOUNDARYNM_R (vel2pR, oddR)
	CALL BOUNDARYNM_R (vel2pL, oddR)
END IF
IF (dual_energy == 1) THEN
	CALL BOUNDARYNM_R (rhoe2R, even)
	CALL BOUNDARYNM_R (rhoe2L, even)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL BOUNDARYNM_R (sca2R(:,:,i), even)
		CALL BOUNDARYNM_R (sca2L(:,:,i), even)
	END DO
END IF

! Check the grid !
CALL CHECKGRIDNM_R (p2R)
CALL CHECKGRIDNM_R (p2L)
CALL CHECKGRIDNM_R (rho2R)
CALL CHECKGRIDNM_R (rho2L)
CALL CHECKGRIDNM_R (vel2rR)
CALL CHECKGRIDNM_R (vel2rL)
CALL CHECKGRIDNM_R (vel2zR)
CALL CHECKGRIDNM_R (vel2zL)
IF(rotationnm_flag == 1) THEN
	CALL CHECKGRIDNM_R (vel2pR)
	CALL CHECKGRIDNM_R (vel2pL)
END IF
IF (dual_energy == 1) THEN
	CALL CHECKGRIDNM_R (rhoe2R)
	CALL CHECKGRIDNM_R (rhoe2L)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL CHECKGRIDNM_R (sca2R(:,:,i))
		CALL CHECKGRIDNM_R (sca2L(:,:,i))
	END DO
END IF

! Do high order interpolation !
CALL HIGHORDERNM_R (p2R)
CALL HIGHORDERNM_R (p2L)
CALL HIGHORDERNM_R (rho2R)
CALL HIGHORDERNM_R (rho2L)
CALL HIGHORDERNM_R (vel2rR)
CALL HIGHORDERNM_R (vel2rL)
CALL HIGHORDERNM_R (vel2zR)
CALL HIGHORDERNM_R (vel2zL)
IF(rotationnm_flag == 1) THEN
	CALL HIGHORDERNM_R (vel2pR)
	CALL HIGHORDERNM_R (vel2pL)
END IF
IF (dual_energy == 1) THEN
	CALL HIGHORDERNM_R (rhoe2R)
	CALL HIGHORDERNM_R (rhoe2L)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL HIGHORDERNM_R (sca2R(:,:,i))
		CALL HIGHORDERNM_R (sca2L(:,:,i))
	END DO
END IF

! Copy to boundaries !
CALL BOUNDARYNM_R (p2R, even)
CALL BOUNDARYNM_R (p2L, even)
CALL BOUNDARYNM_R (rho2R, even)
CALL BOUNDARYNM_R (rho2L, even)
CALL BOUNDARYNM_R (vel2rR, oddR)
CALL BOUNDARYNM_R (vel2rL, oddR)
CALL BOUNDARYNM_R (vel2zR, oddZ)
CALL BOUNDARYNM_R (vel2zL, oddZ)
IF(rotationnm_flag == 1) THEN
	CALL BOUNDARYNM_R (vel2pR, oddR)
	CALL BOUNDARYNM_R (vel2pL, oddR)
END IF
IF (dual_energy == 1) THEN
	CALL BOUNDARYNM_R (rhoe2R, even)
	CALL BOUNDARYNM_R (rhoe2L, even)
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		CALL BOUNDARYNM_R (sca2R(:,:,i), even)
		CALL BOUNDARYNM_R (sca2L(:,:,i), even)
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check the grid flag conditions for DM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKGRIDDM_R (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Dummy !
REAL (DP) :: dummy

! Dummy variables
INTEGER :: j, k

! Then do the interpolation !
DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
	DO j = 1, length_step_r_part_1
		dummy = qm2*array(j-2,k) + qm1*array(j-1,k) + qc*array(j,k) + qp1*array(j+1,k) + qp2*array(j+2,k)
		IF((abs(array(j,k)) < 2.0D0*abs(dummy)) .AND. (abs(dummy) > 1.0D-20)) THEN
			gridflag1(j,k) = 1
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the r-directions for DM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HIGHORDERDM_R (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Temporal array !
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: temporal

! Dummy variables
INTEGER :: j, k

! First, assign !
temporal = array 

! Then do the interpolation !
DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
	DO j = 1, length_step_r_part_1
		IF(gridflag1(j,k) == 1) THEN
			array(j,k) = temporal(j,k)
		ELSE
			array(j,k) = temporal(j,k) - (qm2*temporal(j-2,k) + qm1*temporal(j-1,k) + & 
			      	     qc*temporal(j,k) + qp1*temporal(j+1,k) + qp2*temporal(j+2,k))
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check the grid flag conditions for DM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKGRIDDM_Z (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Dummy !
REAL (DP) :: dummy

! Dummy variables
INTEGER :: j, k

! Then do the interpolation !
DO j = 0, length_step_r_part_1
	DO k = length_step_z_min_part_1, length_step_z_part_1
		dummy = qm2*array(j,k-2) + qm1*array(j,k-1) + qc*array(j,k) + qp1*array(j,k+1) + qp2*array(j,k+2)
		IF((abs(array(j,k)) < 2.0D0*abs(dummy)) .AND. (abs(dummy) > 1.0D-20)) THEN
			gridflag1(j,k) = 1
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the z-directions for DM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HIGHORDERDM_Z (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: array

! Temporal array !
REAL (DP), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5) :: temporal

! Dummy variables
INTEGER :: j, k

! First, assign !
temporal = array 

! Then do the interpolation !
DO j = 0, length_step_r_part_1
	DO k = length_step_z_min_part_1, length_step_z_part_1
		IF(gridflag1(j,k) == 1) THEN
			array(j,k) = temporal(j,k)
		ELSE
			array(j,k) = temporal(j,k) - (qm2*temporal(j,k-2) + qm1*temporal(j,k-1) + & 
			             qc*temporal(j,k) + qp1*temporal(j,k+1) + qp2*temporal(j,k+2))
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check the grid flag conditions for NM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKGRIDNM_R (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Dummy !
REAL (DP) :: dummy

! Dummy variables
INTEGER :: j, k

! Then do the interpolation !
DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
	DO j = 1, length_step_r_part_2
		dummy = qm2*array(j-2,k) + qm1*array(j-1,k) + qc*array(j,k) + qp1*array(j+1,k) + qp2*array(j+2,k)
		IF((abs(array(j,k)) < 2.0D0*abs(dummy)) .AND. (abs(dummy) > 1.0D-20)) THEN
			gridflag2(j,k) = 1
		END IF
	END DO
END DO

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the r-directions for NM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HIGHORDERNM_R (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Temporal array !
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: temporal

! Dummy variables !
REAL (DP) :: dummy

! Dummy variables
INTEGER :: j, k

! First, assign !
temporal = array 

! Then do the interpolation !
DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
	DO j = 1, length_step_r_part_2
		IF(gridflag2(j,k) == 1) THEN
			array(j,k) = temporal(j,k)
		ELSE
			array(j,k) = temporal(j,k) - (qm2*temporal(j-2,k) + qm1*temporal(j-1,k) + & 
			      	     qc*temporal(j,k) + qp1*temporal(j+1,k) + qp2*temporal(j+2,k))
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check the grid flag conditions for NM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKGRIDNM_Z (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Temporal array !
REAL (DP) :: dummy

! Dummy variables
INTEGER :: j, k

! Then do the interpolation !
DO j = 0, length_step_r_part_2
	DO k = length_step_z_min_part_2, length_step_z_part_2
		dummy = qm2*array(j,k-2) + qm1*array(j,k-1) + qc*array(j,k) + qp1*array(j,k+1) + qp2*array(j,k+2)
		IF((abs(array(j,k)) < 2.0D0*abs(dummy)) .AND. (abs(dummy) > 1.0D-20)) THEN
			gridflag2(j,k) = 1
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do high order interpolation along the z-directions for NM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HIGHORDERNM_Z (array)
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: array

! Temporal array !
REAL (DP), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5) :: temporal

! Dummy variables
INTEGER :: j, k

! Dummy variables !
REAL (DP) :: dummy

! First, assign !
temporal = array 

! Then do the interpolation !
DO j = 0, length_step_r_part_2
	DO k = length_step_z_min_part_2, length_step_z_part_2
		IF(gridflag2(j,k) == 1) THEN
			array(j,k) = temporal(j,k)
		ELSE
			array(j,k) = temporal(j,k) - (qm2*temporal(j,k-2) + qm1*temporal(j,k-1) + & 
			             qc*temporal(j,k) + qp1*temporal(j,k+1) + qp2*temporal(j,k+2))
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Boundary conditions for DM states and fluxes along the r-directions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYDM_R (array, sign)
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

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Boundary conditions for DM states and fluxes along the z-directions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYDM_Z (array, sign)
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

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Boundary conditions for NM states and fluxes along the r-directions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYNM_R (array, sign)
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

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Boundary conditions for NM states and fluxes along the z-directions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYNM_Z (array, sign)
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

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!