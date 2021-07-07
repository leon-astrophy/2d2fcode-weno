!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up the position of each grid
! Notice if dx is changed, this subroutine needs
! to be called every time
!
! This subroutine also finds the local volume to save computing time
!
! Note: Here the Cart. grid style is used, meaning that the 
! r- and z- position is a constant along a constant z and r respectively
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetGridNM
USE definition
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k, x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Get the z-position
IF(coordinate_flag == 0) THEN
   DO k = -4, length_step_z_2 + 5, 1
      z2(k) = (DBLE(k) - 0.5D0) * dx2
      zF2(k) = DBLE(k)* dx2
   ENDDO
ELSEIF(coordinate_flag == 1) THEN
   IF(hemisphere_flag == 1) THEN
      x = length_step_z_2/2
      DO k = -4, length_step_z_2 + 5, 1
         z2(k) = (DBLE(k - x) - 0.5D0) * dx2
     	 zF2(k) = DBLE(k - x)* dx2
      ENDDO
   ELSE
      DO k = -4, length_step_z_2 + 5, 1
         z2(k) = (DBLE(k) - 0.5D0) * dx2
	 zF2(k) = DBLE(k)* dx2
      ENDDO
   ENDIF
ENDIF

! Get the r-position
IF(coordinate_flag == 0) THEN
   DO j = -4, length_step_r_2 + 5, 1
      r2(j) = (DBLE(j) - 0.5D0) * dx2
      rF2(j) = DBLE(j) * dx2
   ENDDO
   DO j = -4, length_step_r_2 + 5, 1
      sca2_fac1(j) = 1.0D0
      sca2_fac2(j) = 1.0D0
   ENDDO
ELSEIF(coordinate_flag == 1) THEN
   DO j = -4, length_step_r_2 + 5, 1
      r2(j) = (DBLE(j) - 0.5D0) * dx2
      rF2(j) = DBLE(j) * dx2
   ENDDO
   DO j = -4, length_step_r_2 + 5, 1
      sca2_fac1(j) = 1.0D0
      sca2_fac2(j) = 1.0D0
   ENDDO
ENDIF

! Radial distances !
DO j = -4, length_step_r_2 + 5, 1
	DO k = -4, length_step_z_2 + 5, 1
		rad2(j,k) = SQRT(r2(j)**2 + z2(k)**2)
	END DO
END DO

! Get maxmimum distance !
rmax2 = max(r2(length_step_r_2), z2(length_step_z_2))

! Dimensionless distances !
DO j = -4, length_step_r_2 + 5, 1
	DO k = -4, length_step_z_2 + 5, 1
		radbar2(j,k) = rad2(j,k)/rmax2
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the volume 
IF(coordinate_flag == 0) THEN
   DO j = -4, length_step_r_2 + 5, 1
      vol2(j,:) = dx2 * dx2
   ENDDO
ELSEIF(coordinate_flag == 1) THEN
   IF(hemisphere_flag == 1) THEN
      ! You need to change this volume if we use
      ! stricktly Cartesian coordinate
      DO j = -4, length_step_r_2 + 5, 1
         vol2(j,:) = pi * (rF2 (j)**2 - rF2 (j-1)**2) * dx2
      ENDDO
   ELSEIF(hemisphere_flag == 0) THEN
      DO j = -4, length_step_r_2 + 5, 1
         vol2(j,:) = 2.0D0 * pi * (rF2 (j)**2 - rF2 (j-1)**2) * dx2
      ENDDO
   ELSE
      STOP 'Check hemisphere flag, stopped at GetGrid'
   ENDIF
ENDIF

! Dimensionless volume !
DO j = -4, length_step_r_2 + 5, 1
	DO k = -4, length_step_z_2+ 5, 1
		volbar2(j,k) = vol2(j,k)/rmax2**3
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do everything the same, but for DM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetGridDM
USE definition
IMPLICIT NONE

! Dummy variables
INTEGER :: j, k, x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For DM motherboard !

! Get the z-position
IF(coordinate_flag == 0) THEN
	DO k = -4, length_step_z_1 + 5, 1
      		z1(k) = (DBLE(k) - 0.5D0) * dx1
		zF1(k) = DBLE(k)* dx1
   	ENDDO
ELSEIF(coordinate_flag == 1) THEN
   	IF(hemisphere_flag == 1) THEN
      		x = length_step_z_1/2
      		DO k = -4, length_step_z_1 + 5, 1
         		z1(k) = (DBLE(k - x) - 0.5D0) * dx1
			zF1(k) = DBLE(k - x)* dx1
      		ENDDO
   	ELSE
		DO k = -4, length_step_z_1 + 5, 1
         		z1(k) = (DBLE(k) - 0.5D0) * dx1
			zF1(k) = DBLE(k)* dx1
      		ENDDO
   	ENDIF
ENDIF

! Get the r-position
IF(coordinate_flag == 0) THEN
   	DO j = -4, length_step_r_1 + 5, 1
      		r1(j) = (DBLE(j) - 0.5D0) * dx1
		rF1(j) = DBLE(j) * dx1
   	ENDDO
   	DO j = -4, length_step_r_1 + 5, 1
      		sca1_fac1(j) = 1.0D0
      		sca1_fac2(j) = 1.0D0
   	ENDDO
ELSEIF(coordinate_flag == 1) THEN
  	 DO j = -4, length_step_r_1 + 5, 1
      		r1(j) = (DBLE(j) - 0.5D0) * dx1
                rF1(j) = DBLE(j) * dx1
   	 ENDDO
   	 DO j = -4, length_step_r_1 + 5, 1
      		sca1_fac1(j) = 1.0D0
      		sca1_fac2(j) = 1.0D0
   	 ENDDO
ENDIF

! Radial distances !
DO j = -4, length_step_r_1 + 5, 1
	DO k = -4, length_step_z_1 + 5, 1
		rad1(j,k) = SQRT(r1(j)**2 + z1(k)**2)
	END DO
END DO

! Get maxmimum distance !
rmax1 = max(r1(length_step_r_1), z1(length_step_z_1))

! Dimensionless distances !
DO j = -4, length_step_r_1 + 5, 1
	DO k = -4, length_step_z_1+ 5, 1
		radbar1(j,k) = rad1(j,k)/rmax1
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For DM !

IF(coordinate_flag == 0) THEN
   	DO j = -4, length_step_r_1 + 5, 1
      		vol1(j,:) = dx1 * dx1
   	ENDDO
ELSEIF(coordinate_flag == 1) THEN
   	IF(hemisphere_flag == 1) THEN
      		! You need to change this volume if we use
      		! stricktly Cartesian coordinate
      		DO j = -4, length_step_r_1 + 5, 1
         		vol1(j,:) = pi * (rF1 (j)**2 - rF1 (j-1)**2) * dx1
      		ENDDO
   	ELSEIF(hemisphere_flag == 0) THEN
      		DO j = -4, length_step_r_1 + 5, 1
         		vol1(j,:) = 2.0D0 * pi * (rF1 (j)**2 - rF1 (j-1)**2) * dx1
      		ENDDO
   	ELSE
      		STOP 'Check hemisphere flag, stopped at GetGrid'
   	ENDIF
ENDIF

! Dimensionless volume !
DO j = -4, length_step_r_1 + 5, 1
	DO k = -4, length_step_z_1 + 5, 1
		volbar1(j,k) = vol1(j,k)/rmax1**3
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE