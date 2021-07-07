!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares the data for spatial discretization,
! due ask the WENO_module to do the reconstruction
! and then combines the results for one Runge-Kutta sub-step
!
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Extended by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Search for the following section for a fast jump
! Section 1: Building f_r, f_z, sa_r, sa_z and sb
! Section 2: WENO reconstruction for DM
! Section 3: WENO reconstruction for NM 
! Section 3: Calculate l
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPATIAL (u1, u2)
USE DEFINITION
USE RIEMANN_MODULE
USE WENO_MODULE
USE HELMEOS_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
IMPLICIT NONE

! dummy variables !
INTEGER :: i, j, k

! Limits of density to be considered
REAL (DP) :: rho_min1, rho_min2

! The input U array
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5, imin1 : imax1) :: u1
REAL (DP), INTENT (IN), DIMENSION (-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5, imin2 : imax2) :: u2

! The auxillary array for the flux term, DM
REAL (DP), allocatable, DIMENSION (:,:,:) :: sa1, sb1
REAL (DP), allocatable, DIMENSION (:,:,:) :: fluxr_1, fluxz_1
REAL (DP), allocatable, DIMENSION (:,:,:) :: dfdr1, dfdz1

! Flux arrays for NM !
REAL (DP), allocatable, DIMENSION (:,:,:) :: sa2, sb2
REAL (DP), allocatable, DIMENSION (:,:,:) :: fluxr_2, fluxz_2
REAL (DP), allocatable, DIMENSION (:,:,:) :: dfdr2, dfdz2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check timing with or without openmp
!INTEGER :: time_start, time_end
!INTEGER :: cr, cm
!REAL :: rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count timing 
!CALL system_clock(count_rate=cr)
!CALL system_clock(count_max=cm)
!rate = REAL(cr)               
!WRITE(*,*) "system_clock rate ",rate
!CALL system_clock(time_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First allocate the arrays for future use
! Debug stuff 
!IF(debug_flag == 1) WRITE(*,*) 'In Spatial'
!IF(debug_flag == 1) WRITE(*,*) 'Allocate variables'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate for NM !
ALLOCATE(dfdr2(-4:length_step_r_2+5, -4:length_step_z_2+5, imin2:imax2))
ALLOCATE(dfdz2(-4:length_step_r_2+5, -4:length_step_z_2+5, imin2:imax2))
ALLOCATE(fluxr_2(-4:length_step_r_2+5, -4:length_step_z_2+5, imin2:imax2))
ALLOCATE(fluxz_2(-4:length_step_r_2+5, -4:length_step_z_2+5, imin2:imax2))
ALLOCATE(sa2(-4:length_step_r_2+5, -4:length_step_z_2+5, imin2:imax2))
ALLOCATE(sb2(-4:length_step_r_2+5, -4:length_step_z_2+5, imin2:imax2))

! For DM !
IF(RUNDM_flag == 1) THEN
	ALLOCATE(dfdr1(-4:length_step_r_1+5, -4:length_step_z_1+5, imin1:imax1))
	ALLOCATE(dfdz1(-4:length_step_r_1+5, -4:length_step_z_1+5, imin1:imax1))
	ALLOCATE(fluxr_1(-4:length_step_r_1+5, -4:length_step_z_1+5, imin1:imax1))
	ALLOCATE(fluxz_1(-4:length_step_r_1+5, -4:length_step_z_1+5, imin1:imax1))
	ALLOCATE(sa1(-4:length_step_r_1+5, -4:length_step_z_1+5, imin1:imax1))
	ALLOCATE(sb1(-4:length_step_r_1+5, -4:length_step_z_1+5, imin1:imax1))
END IF

! Allocate arrays for Riemann Solver !
CALL BUILDRIEMANN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare the flux array
!IF(debug_flag == 1) WRITE(*,*) 'Prepare array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find frame velocity 
IF (movinggriddm_flag == 1) THEN
	CALL FINDFRAMEVEL_DM
END IF
IF (movinggridnm_flag == 1) THEN
	CALL FINDFRAMEVEL_NM
END IF

! Initialize source term and fluxes !
IF(RUNDM_flag == 1) THEN
	sa1 = 0.0D0
	sb1 = 0.0D0
	fluxr_1 = 0.0D0
	fluxz_1 = 0.0D0
	dfdr1 = 0.0D0
	dfdz1 = 0.0D0
END IF

! For NM !
sa2 = 0.0D0
sb2 = 0.0D0
fluxr_2 = 0.0D0
fluxz_2 = 0.0D0
dfdr2 = 0.0D0
dfdz2 = 0.0D0

! Assign scalars !
IF(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -4, length_step_r_2 + 5
			DO k = -4, length_step_z_2 + 5
				sca2(j,k,i) = u2(j,k,i)/u2(j,k,irho2)
			END DO
		END DO
	END DO
END IF

! Find effective speed !
CALL FINDALPHA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1: Reconstructions and build states in the horizontal directions !

! Initialize !
CALL CLEARRIEMANN

! Reconstruct using WENO !
CALL WENO_R

! Build the fluxes and states !
CALL BUILDSTATES_R

! Choose an appropriate riemann solver !
IF(RUNDM_flag == 1) THEN
	CALL LFDM_R(fluxr_1)
END IF

! For NM !
CALL LFNM_R(fluxr_2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2: Reconstructions and build states in the vertical directions !

! Initialize !
CALL CLEARRIEMANN

! Reconstruct using WENO !
CALL WENO_Z

! Build the fluxes and states !
CALL BUILDSTATES_Z

! Choose an appropriate riemann solver !
IF(RUNDM_flag == 1) THEN
	CALL LFDM_Z(fluxz_1)
END IF

! For NM !
CALL LFNM_Z(fluxz_2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 3: Get source terms

! For DM sector
IF(runDM_flag == 1) THEN
	rho_min1 = 1.1D0 * rho1_a

	! Gravitational source term !
	DO k = length_step_z_min_part_1, length_step_z_part_1
		DO j = 1, length_step_r_part_1
			IF(rho1(j,k) > rho_min1) THEN
				sb1(j,k,ivel1_r) = rho1(j,k) * phi1_r(j,k)
				sb1(j,k,ivel1_z) = rho1(j,k) * phi1_z(j,k)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! No DM epsilon equation !
				!sb1(j,k,itau1) = rho1(j,k)*(vel1_r(j,k)*phi1_r(j,k) + vel1_z(j,k)*phi1_z(j,k))
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			END IF
		ENDDO
	ENDDO

	! Add source term to momentum equation according to the coordinate system !
	IF(coordinate_flag == 1) THEN
		DO k = length_step_z_min_part_1, length_step_z_part_1
			DO j = 1, length_step_r_part_1
				sa1(j,k,ivel1_r) = - p1(j,k) - rho1(j,k)*vel1_p(j,k)**2
			ENDDO
		ENDDO
		IF(rotationdm_flag == 1) THEN
			DO k = length_step_z_min_part_1, length_step_z_part_1
				DO j = 1, length_step_r_part_1
					sa1(j,k,ivel1_p) = rho1(j,k)*vel1_p(j,k)*vel1_r(j,k)
				ENDDO
			ENDDO
		END IF
	END IF

	! Moving grid !
	IF(movinggriddm_flag == 1) THEN
		DO k = length_step_z_min_part_1, length_step_z_part_1
			DO j = 1, length_step_r_part_1
				DO i = imin1, imax1
					sb1 (j,k,i) = sb1 (j,k,i) + u1 (j,k,i) * (3.0D0 * vel1_max / radius1)
				ENDDO
			END DO
		ENDDO
	END IF

END IF

! Threshold density !
rho_min2 = 1.1D0 * rho2_a

! For NM sector
DO k = length_step_z_min_part_2, length_step_z_part_2
	DO j = 1, length_step_r_part_2
		IF(rho2(j,k) > rho_min2) THEN
			sb2(j,k,ivel2_r) = rho2(j,k) * phi2_r(j,k)
			sb2(j,k,ivel2_z) = rho2(j,k) * phi2_z(j,k)
		END IF
	END DO
END DO
IF(nm_epsilon == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 1, length_step_r_part_2
			IF(rho2(j,k) > rho_min2) THEN
				sb2(j,k,itau2) = rho2(j,k)*(vel2_r(j,k)*phi2_r(j,k) + vel2_z(j,k)*phi2_z(j,k))
			END IF
		END DO
	END DO
END IF

IF(dual_energy == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 1, length_step_r_part_2
			sb2(j,k,ieps2) = - (vel2_r(j,k)*dp2dr(j,k) + vel2_z(j,k)*dp2dz(j,k))
		END DO
	END DO
END IF

! Choose coordinate system 
IF(coordinate_flag == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 1, length_step_r_part_2
			sa2(j,k,ivel2_r) = - p2(j,k) - rho2(j,k)*vel2_p(j,k)**2
		END DO
	END DO
	IF(rotationnm_flag == 1) THEN
		DO k = length_step_z_min_part_2, length_step_z_part_2
			DO j = 1, length_step_r_part_2
				sa2(j,k,ivel2_p) = rho2(j,k)*vel2_p(j,k)*vel2_r(j,k)
			ENDDO
		ENDDO
	END IF
END IF

! Extra source term for SGS turbulence
IF(turb_flag == 1) then

   DO k = length_step_z_min_part_2, length_step_z_part_2 , 1
      DO j = 1, length_step_r_part_2, 1
         sb2 (j,k,iturbq) = - turb_source(j,k) - turb_divergence(j,k)
      ENDDO
   ENDDO

   DO k = length_step_z_min_part_2, length_step_z_part_2 , 1
      DO j = 1, length_step_r_part_2, 1
         sb2 (j,k,itau2) = sb2 (j,k,itau2) + turb_source(j,k)                      
      ENDDO
   ENDDO

   IF(dual_energy == 1) THEN
      DO k = length_step_z_min_part_2, length_step_z_part_2 , 1
         DO j = 1, length_step_r_part_2, 1
            sb2 (j,k,ieps2) = sb2 (j,k,ieps2) + turb_source(j,k)                     
         ENDDO
      ENDDO
   END IF

ENDIF

! Moving grid !
IF(movinggridnm_flag == 1) THEN
	DO k = length_step_z_min_part_2, length_step_z_part_2
		DO j = 1, length_step_r_part_2
			DO i = imin2, imax2
				sb2 (j,k,i) = sb2 (j,k,i) + u2 (j,k,i) * (3.0D0 * vel2_max / radius2)
			ENDDO
		END DO
	ENDDO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 4: Get flux derivatives

! Now calculate dfdx accordingly to the corrected flux, do it case by case
IF(coordinate_flag == 0) THEN
	IF(RUNDM_flag == 1) THEN
		DO i = imin1, imax1
			DO j = 1, length_step_r_part_1, 1
				DO k = length_step_z_min_part_1, length_step_z_part_1
					dfdr1 (j, k, i) = (fluxr_1 (j, k, i) - fluxr_1 (j - 1, k, i)) / dx1
					dfdz1 (j, k, i) = (fluxz_1 (j, k, i) - fluxz_1 (j, k - 1, i)) / dx1
				END DO
			END DO
		END DO
	END IF

	! Do for NM !
	DO i = imin2, imax2
		DO j = 1, length_step_r_part_2, 1
			DO k = length_step_z_min_part_2, length_step_z_part_2, 1
				dfdr2 (j, k, i) = (fluxr_2 (j, k, i) - fluxr_2 (j - 1, k, i)) / dx2
				dfdz2 (j, k, i) = (fluxz_2 (j, k, i) - fluxz_2 (j, k - 1, i)) / dx2
			END DO
		END DO
	END DO
ELSEIF (coordinate_flag == 1) THEN
	IF(RUNDM_flag == 1) THEN
		DO i = imin1, imax1
			DO j = 1, length_step_r_part_1, 1
				DO k = length_step_z_min_part_1, length_step_z_part_1
					dfdr1 (j, k, i) = (rF1(j)*fluxr_1 (j, k, i) - rF1(j-1)*fluxr_1 (j - 1, k, i)) / (r1(j)*dx1)
					dfdz1 (j, k, i) = (fluxz_1 (j, k, i) - fluxz_1 (j, k - 1, i)) / dx1
				END DO
			END DO
		END DO
	END IF

	! Do for NM !
	DO i = imin2, imax2
		DO j = 1, length_step_r_part_2, 1
			DO k = length_step_z_min_part_2, length_step_z_part_2, 1
				dfdr2 (j, k, i) = (rF2(j)*fluxr_2 (j, k, i) - rF2(j-1)*fluxr_2 (j - 1, k, i)) / (r2(j)*dx2)
				dfdz2 (j, k, i) = (fluxz_2 (j, k, i) - fluxz_2 (j, k - 1, i)) / dx2
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 5: Calculate l (SSPRK54)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF(debug_flag == 1) WRITE(*,*) 'Obtain l'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization
l1 = 0.0D0
l2 = 0.0D0

IF(coordinate_flag == 0) THEN

   IF(runDM_flag == 1) THEN

      DO i = imin1, imax1, 1                
         DO k = length_step_z_min_part_1, length_step_z_part_1, 1
            DO j = 1, length_step_r_part_1, 1
               l1(j,k,i) = - dfdr1 (j,k,i) - dfdz1 (j,k,i) - sb1 (j,k,i) 
            ENDDO
         ENDDO
      ENDDO

   ENDIF

   DO i = imin2, imax2, 1

      DO k = length_step_z_min_part_2, length_step_z_part_2, 1
         DO j = 1, length_step_r_part_2, 1
            l2(j,k,i) = - dfdr2 (j,k,i) - dfdz2 (j,k,i) - sb2 (j,k,i) 
         ENDDO
      ENDDO

   ENDDO

ELSEIF(coordinate_flag == 1) THEN

   IF(runDM_flag == 1) THEN

      DO i = imin1, imax1, 1

         DO k = length_step_z_min_part_1, length_step_z_part_1, 1
            DO j = 1, length_step_r_part_1, 1
               l1(j,k,i) = - dfdr1 (j,k,i) - dfdz1 (j,k,i) - sb1 (j,k,i) - sa1 (j,k,i) / r1(j)
            ENDDO
         ENDDO
      ENDDO

   ENDIF

   DO i = imin2, imax2, 1
  
      DO k = length_step_z_min_part_2, length_step_z_part_2, 1
         DO j = 1, length_step_r_part_2, 1
            l2(j,k,i) = - dfdr2 (j,k,i) - dfdz2 (j,k,i) - sb2 (j,k,i) - sa2 (j,k,i) / r2(j)
         ENDDO 
      ENDDO
      
   ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate the auxillary arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deug stuff !
!IF(debug_flag == 1) WRITE(*,*) 'Deallocate arrays'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate for NM !
DEALLOCATE(dfdr2)
DEALLOCATE(dfdz2)
DEALLOCATE(fluxr_2)
DEALLOCATE(fluxz_2)
DEALLOCATE(sa2)
DEALLOCATE(sb2)

! For DM !
IF(RUNDM_flag == 1) THEN
	DEALLOCATE(dfdr1)
	DEALLOCATE(dfdz1)
	DEALLOCATE(fluxr_1)
	DEALLOCATE(fluxz_1)
	DEALLOCATE(sa1)
	DEALLOCATE(sb1)
END IF

! Deallocate arrays for Riemann Solver !
CALL DESTROYRIEMANN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count time 
!CALL system_clock(time_end)
!WRITE(*,*) REAL(time_end - time_start) / rate
!WRITE(*,*) 'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the horizontal directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES_R
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build conserved variables !
If (RUNDM_flag == 1) then
	DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
		DO j = 0, length_step_r_part_1
			uL1 (j, k, irho1) = rho1L (j, k)
			uL1 (j, k, ivel1_r) = rho1L (j, k) * vel1rL (j, k)
			uL1 (j, k, ivel1_z) = rho1L (j, k) * vel1zL (j, k)
			uR1 (j, k, irho1) = rho1R (j, k)
			uR1 (j, k, ivel1_r) = rho1R (j, k) * vel1rR (j, k)
			uR1 (j, k, ivel1_z) = rho1R (j, k) * vel1zR (j, k)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! No DM internal energy !
			!uL1 (j, k, itau1) = rho1L (j, k) * eps1L (j, k) + 5.0E-1_DP * rho1L (j, k) * (vel1rL (j, k) ** 2 + vel1zL (j, k) ** 2 + vel1pL (j, k) ** 2)
			!uR1 (j, k, itau1) = rho1R (j, k) * eps1R (j, k) + 5.0E-1_DP * rho1R (j, k) * (vel1rR (j, k) ** 2 + vel1zR (j, k) ** 2 + vel1pR (j, k) ** 2)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   	END DO
	END DO

	! Rotation !
	IF(rotationdm_flag == 1) THEN
		DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
			DO j = 0, length_step_r_part_1
				uL1 (j, k, ivel1_p) = rho1L (j, k) * vel1pL (j, k)
				uR1 (j, k, ivel1_p) = rho1R (j, k) * vel1pR (j, k)
		 	END DO
		END DO
	END IF
END IF

! Do the same for NM !
DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
	DO j = 0, length_step_r_part_2
		uL2 (j, k, irho2) = rho2L (j, k)
		uL2 (j, k, ivel2_r) = rho2L (j, k) * vel2rL (j, k)
		uL2 (j, k, ivel2_z) = rho2L (j, k) * vel2zL (j, k)
		uR2 (j, k, irho2) = rho2R (j, k)
		uR2 (j, k, ivel2_r) = rho2R (j, k) * vel2rR (j, k)
		uR2 (j, k, ivel2_z) = rho2R (j, k) * vel2zR (j, k)	
	END DO
END DO

! Rotation !
IF(rotationnm_flag == 1) THEN
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			uL2 (j, k, ivel2_p) = rho2L (j, k) * vel2pL (j, k)
			uR2 (j, k, ivel2_p) = rho2R (j, k) * vel2pR (j, k)
		END DO
	END DO
END IF

! NM energy equation !
IF(nm_epsilon == 1) THEN
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			uL2 (j, k, itau2) = rho2L (j,k) * eps2L (j,k) + 5.0E-1_DP * rho2L (j,k) * (vel2rL (j,k) ** 2 + vel2zL (j,k) ** 2 + vel2pL (j,k) ** 2)
			uR2 (j, k, itau2) = rho2R (j,k) * eps2R (j,k) + 5.0E-1_DP * rho2R (j,k) * (vel2rR (j,k) ** 2 + vel2zR (j,k) ** 2 + vel2pR (j,k) ** 2)
		END DO
	END DO
END IF
IF(dual_energy == 1) THEN 
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			uL2 (j, k, ieps2) = rhoe2L(j,k)
			uR2 (j, k, ieps2) = rhoe2R(j,k)
		END DO
	END DO
END IF

! NM scalar equation !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
			DO j = 0, length_step_r_part_2
				uL2 (j, k, i) = rho2L (j,k) * sca2L(j,k,i)
				uR2 (j, k, i) = rho2R (j,k) * sca2R(j,k,i)
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build fluxes For DM !
If (RUNDM_flag == 1) then
   	DO i = imin1, imax1, 1
		DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
			DO j = 0, length_step_r_part_1
				fluxL1 (j,k,i) = uL1 (j,k,i) * vel1rL (j,k)
				fluxR1 (j,k,i) = uR1 (j,k,i) * vel1rR (j,k)
			ENDDO
		END DO
	ENDDO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
   	! Add the presusre work done term to the energy equation             
	!DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
	!	DO j = 0, length_step_r_part_1
	!		fluxL1 (j,k,itau1) = fluxL1 (j,k,itau1) + p1L(j,k) * vel1rL(j,k)
	!		fluxR1 (j,k,itau1) = fluxR1 (j,k,itau1) + p1R(j,k) * vel1rR(j,k)
	!	END DO
      	!ENDDO   
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   	! Add the pressure term to r-momentum equation !
	DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
		DO j = 0, length_step_r_part_1
         		fluxL1 (j,k,ivel1_r) = fluxL1 (j,k,ivel1_r) + p1L(j,k)
			fluxR1 (j,k,ivel1_r) = fluxR1 (j,k,ivel1_r) + p1R(j,k)
      		ENDDO
	END DO

	! Moving grid !
	If(movinggriddm_flag == 1) THEN
		DO i = imin1, imax1, 1
			DO k = length_step_z_min_part_1 - 5, length_step_z_part_1 + 5
				DO j = 0, length_step_r_part_1
         				fluxL1 (j,k,i) = fluxL1 (j,k,i) - uL1(j,k,i)*vf1rL(j,k)
					fluxR1 (j,k,i) = fluxR1 (j,k,i) - uR1(j,k,i)*vf1rR(j,k)
				END DO
      			ENDDO
		END DO
	END IF

END IF

! For NM !
DO i = imin2, imax2, 1
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			fluxL2 (j,k,i) = uL2 (j,k,i) * vel2rL (j,k)
			fluxR2 (j,k,i) = uR2 (j,k,i) * vel2rR (j,k)
		END DO
	ENDDO
ENDDO

! Add the presusre work done term to the energy equation             
IF(nm_epsilon == 1) THEN   
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			fluxL2 (j,k,itau2) = fluxL2 (j,k,itau2) + p2L(j,k) * vel2rL(j,k)
			fluxR2 (j,k,itau2) = fluxR2 (j,k,itau2) + p2R(j,k) * vel2rR(j,k)
		END DO
	ENDDO 
END IF  
IF(dual_energy == 1) THEN    
	DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
		DO j = 0, length_step_r_part_2
			fluxL2 (j,k,ieps2) = fluxL2 (j,k,ieps2) + p2L(j,k) * vel2rL(j,k)
			fluxR2 (j,k,ieps2) = fluxR2 (j,k,ieps2) + p2R(j,k) * vel2rR(j,k)
		ENDDO 
	END DO
END IF  	

! Add the pressure term to r-momentum equation !
DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
	DO j = 0, length_step_r_part_2
		fluxL2 (j,k,ivel2_r) = fluxL2 (j,k,ivel2_r) + p2L(j,k)
		fluxR2 (j,k,ivel2_r) = fluxR2 (j,k,ivel2_r) + p2R(j,k)
	END DO
ENDDO

! Extra flux term for moving grid !
If(movinggridnm_flag == 1) THEN
	DO i = imin2, imax2, 1
		DO k = length_step_z_min_part_2 - 5, length_step_z_part_2 + 5
			DO j = 0, length_step_r_part_2
         			fluxL2 (j,k,i) = fluxL2 (j,k,i) - uL2(j,k,i)*vf2rL(j,k)
				fluxR2 (j,k,i) = fluxR2 (j,k,i) - uR2(j,k,i)*vf2rR(j,k)
			END DO
      		ENDDO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine built the states for left and right edges for the vertical directions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES_Z
USE DEFINITION 
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build conserved variables !
If (RUNDM_flag == 1) then
	DO j = -4, length_step_r_part_1 + 5
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
			uL1 (j, k, irho1) = rho1L (j, k)
			uL1 (j, k, ivel1_r) = rho1L (j, k) * vel1rL (j, k)
			uL1 (j, k, ivel1_z) = rho1L (j, k) * vel1zL (j, k)
			uR1 (j, k, irho1) = rho1R (j, k)
			uR1 (j, k, ivel1_r) = rho1R (j, k) * vel1rR (j, k)
			uR1 (j, k, ivel1_z) = rho1R (j, k) * vel1zR (j, k)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! No DM internal energy !
			!uL1 (j, k, itau1) = rho1L (j, k) * eps1L (j, k) + 5.0E-1_DP * rho1L (j, k) * (vel1rL (j, k) ** 2 + vel1zL (j, k) ** 2 + vel1pL (j, k) ** 2)
			!uR1 (j, k, itau1) = rho1R (j, k) * eps1R (j, k) + 5.0E-1_DP * rho1R (j, k) * (vel1rR (j, k) ** 2 + vel1zR (j, k) ** 2 + vel1pR (j, k) ** 2)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   	END DO
	END DO

	! Rotation !
	IF(rotationdm_flag == 1) THEN
		DO j = -4, length_step_r_part_1 + 5
			DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
				uL1 (j, k, ivel1_p) = rho1L (j, k) * vel1pL (j, k)
				uR1 (j, k, ivel1_p) = rho1R (j, k) * vel1pR (j, k)
		 	END DO
		END DO
	END IF

END IF

! Do the same for NM !
DO j = -4, length_step_r_part_2 + 5
	DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
		uL2 (j, k, irho2) = rho2L (j, k)
		uL2 (j, k, ivel2_r) = rho2L (j, k) * vel2rL (j, k)
		uL2 (j, k, ivel2_z) = rho2L (j, k) * vel2zL (j, k)
		uR2 (j, k, irho2) = rho2R (j, k)
		uR2 (j, k, ivel2_r) = rho2R (j, k) * vel2rR (j, k)
		uR2 (j, k, ivel2_z) = rho2R (j, k) * vel2zR (j, k)
	END DO
END DO

! Rotation !
IF(rotationnm_flag == 1) THEN
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			uL2 (j, k, ivel2_p) = rho2L (j, k) * vel2pL (j, k)
			uR2 (j, k, ivel2_p) = rho2R (j, k) * vel2pR (j, k)
		END DO
	END DO
END IF

! NM energy equation !
IF(nm_epsilon == 1) THEN
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			uL2 (j, k, itau2) = rho2L (j,k) * eps2L (j,k) + 5.0E-1_DP * rho2L (j,k) * (vel2rL (j,k) ** 2 + vel2zL (j,k) ** 2 + vel2pL (j,k) ** 2)
			uR2 (j, k, itau2) = rho2R (j,k) * eps2R (j,k) + 5.0E-1_DP * rho2R (j,k) * (vel2rR (j,k) ** 2 + vel2zR (j,k) ** 2 + vel2pR (j,k) ** 2)
		END DO
	END DO
END IF
IF(dual_energy == 1) THEN 
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			uL2 (j, k, ieps2) = rhoe2L(j,k)
			uR2 (j, k, ieps2) = rhoe2R(j,k)
		END DO
	END DO
END IF

! NM scalar equation !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -4, length_step_r_part_2 + 5
			DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
				uL2 (j, k, i) = rho2L (j,k) * sca2L(j,k,i)
				uR2 (j, k, i) = rho2R (j,k) * sca2R(j,k,i)
			END DO
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build fluxes For DM !
If (RUNDM_flag == 1) then
   	DO i = imin1, imax1, 1
		DO j = -4, length_step_r_part_1 + 5
			DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
				fluxL1 (j,k,i) = uL1 (j,k,i) * vel1zL (j,k)
				fluxR1 (j,k,i) = uR1 (j,k,i) * vel1zR (j,k)
			ENDDO
		END DO
	ENDDO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
   	! Add the presusre work done term to the energy equation             
	!DO j = -4, length_step_r_part_1 + 5
	!	DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
	!		fluxL1 (j,k,itau1) = fluxL1 (j,k,itau1) + p1L(j,k) * vel1zL(j,k)
	!		fluxR1 (j,k,itau1) = fluxR1 (j,k,itau1) + p1R(j,k) * vel1zR(j,k)
	!	END DO
      	!ENDDO   
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   	! Add the pressure term to z-momentum equation !
	DO j = -4, length_step_r_part_1 + 5
		DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
         		fluxL1 (j,k,ivel1_z) = fluxL1 (j,k,ivel1_z) + p1L(j,k)
			fluxR1 (j,k,ivel1_z) = fluxR1 (j,k,ivel1_z) + p1R(j,k)
      		ENDDO
	END DO

	! Moving grid !
	If(movinggriddm_flag == 1) THEN
		DO i = imin1, imax1, 1
			DO j = -4, length_step_r_part_1 + 5
				DO k = length_step_z_min_part_1 - 1, length_step_z_part_1
         				fluxL1 (j,k,i) = fluxL1 (j,k,i) - uL1(j,k,i)*vf1zL(j,k)
					fluxR1 (j,k,i) = fluxR1 (j,k,i) - uR1(j,k,i)*vf1zR(j,k)
				END DO
      			ENDDO
		END DO
	END IF

END IF

! For NM !
DO i = imin2, imax2, 1
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			fluxL2 (j,k,i) = uL2 (j,k,i) * vel2zL (j,k)
			fluxR2 (j,k,i) = uR2 (j,k,i) * vel2zR (j,k)
		END DO
	ENDDO
ENDDO

! Add the presusre work done term to the energy equation             
IF(nm_epsilon == 1) THEN   
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			fluxL2 (j,k,itau2) = fluxL2 (j,k,itau2) + p2L(j,k) * vel2zL(j,k)
			fluxR2 (j,k,itau2) = fluxR2 (j,k,itau2) + p2R(j,k) * vel2zR(j,k)
		END DO
	ENDDO 
END IF  
IF(dual_energy == 1) THEN    
	DO j = -4, length_step_r_part_2 + 5
		DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
			fluxL2 (j,k,ieps2) = fluxL2 (j,k,ieps2) + p2L(j,k) * vel2zL(j,k)
			fluxR2 (j,k,ieps2) = fluxR2 (j,k,ieps2) + p2R(j,k) * vel2zR(j,k)
		ENDDO 
	END DO
END IF  	

! Add the pressure term to r-momentum equation !
DO j = -4, length_step_r_part_2 + 5
	DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
		fluxL2 (j,k,ivel2_z) = fluxL2 (j,k,ivel2_z) + p2L(j,k)
		fluxR2 (j,k,ivel2_z) = fluxR2 (j,k,ivel2_z) + p2R(j,k)
	END DO
ENDDO

! Moving grid !
If(movinggridnm_flag == 1) THEN
	DO i = imin2, imax2, 1
		DO j = -4, length_step_r_part_2 + 5
			DO k = length_step_z_min_part_2 - 1, length_step_z_part_2
         			fluxL2 (j,k,i) = fluxL2 (j,k,i) - uL2(j,k,i)*vf2zL(j,k)
				fluxR2 (j,k,i) = fluxR2 (j,k,i) - uR2(j,k,i)*vf2zR(j,k)
			END DO
      		ENDDO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSEPSILON (den, pre, eps, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT (IN) :: type

! Input density !
REAL (DP), INTENT (IN) :: den, pre

! Output value !
REAL (DP), INTENT (OUT) :: eps

! Extra patches by Ivan !
REAL (DP) :: p_poly, eps_poly, eps_thermal

! Fermi-momentum !
REAL (DP) :: fermi

! For DM Output !
IF(type == 1) THEN
	IF(fermieosdm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 1)
		IF (fermi<=1.0E-2_DP) THEN
			eps = a_max1*small_energy(fermi)/den
		ELSE
			eps = a_max1*large_energy(fermi)/den
		END IF
	ELSE
		eps = k_1 * den ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
	END IF

! For NM !
ELSEIF(type == 2) THEN
	IF(fermieosnm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 2)
		IF (fermi<=1.0E-2_DP) THEN
			eps = a_max2*small_energy(fermi)/den
		ELSE
			eps = a_max2*large_energy(fermi)/den
		END IF
	ELSE
		IF(nm_epsilon == 0) THEN
			eps = k_2 * den ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
		ELSE
			eps = pre/den/(gamma2 - 1.0D0)
		END IF
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the speed of sound !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSSOUNDSPEED(pre, den, eps, cs_out, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type

! Input density !
REAL (DP), INTENT (IN) :: pre, den, eps

! Output value !
REAL (DP), INTENT (OUT) :: cs_out

! Local real variables !
REAL (DP) :: fermi, dpdden, dpdeps, dxdrho

! Extra patches by Ivan !
REAL (DP) :: p_poly, eps_poly, eps_thermal

! We do the DM case first !
IF (type == 1) Then
	IF(fermieosdm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 1)
		CALL FINDDXDRHO (dxdrho, den, 1)
		dpdden = a_max1*dxdrho*dpdx(fermi)
		dpdeps = 0.0E0_DP
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	ELSE
		dpdden = k_1 * gamma1 * den ** (gamma1 - 1.0D0)
		dpdeps = 0.0D0
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	END IF

! For NM !
ELSEIF(type == 2) THEN
	IF (fermieosnm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 2)
		CALL FINDDXDRHO (dxdrho, den, 2)
		dpdden = a_max2*dxdrho*dpdx(fermi)
		dpdeps = 0.0E0_DP
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	ELSE
		IF(nm_epsilon == 0) THEN
			dpdden = k_2 * gamma2 * den ** (gamma2 - 1.0D0)
			dpdeps = 0.0D0
			cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
		ELSE
			dpdden = eps * (gamma2 - 1.0E0_DP)
			dpdeps = den * (gamma2 - 1.0E0_DP)
			cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
		END IF
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the suitable frame velocity v_f appears in !
! v = v_f * r/R so that we can keep track on the whole star with   !
! co-expanding grid 						   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDFRAMEVEL_NM
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: k, j, count_sum

! Real parameter !
REAL (DP) :: rho_max2, rad_sum, rad_dist, rad_vel

! Here I find the suitable velocity for the !
! boundary to expand such that the whole star !
! is always contained inside the box !

! Define the density range to be considered
rho_max2 = 1.0D2 * rho2_a
	
! Searching NM atmosphere speed ! 
! Initialize !
rad_sum = 0.0D0
vel2_max = 0.0D0
count_sum = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do for NM first !
!DO k = length_step_z_min_part_2, length_step_z_part_2, 1
!	DO j = 1, length_step_r_part_2, 1
!		if(rho2(j,k) < rho_max2) then
!			rad_dist = DSQRT(r2(j) ** 2 + z2(k) ** 2)
!			rad_vel = (vel2_r(j,k) * r2(j) + vel2_z(j,k) * z2(k)) / rad_dist
!			rad_sum = rad_sum + rad_dist
!			vel2_max = vel2_max + rad_vel
!			count_sum = count_sum + 1
!		endif
!	enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If the star surface is too close to the box boundary
! Artificially increase the boundary velocity
!IF(count_sum /= 0) THEN
!	vel2_max = vel2_max / DBLE(count_sum)
!	radius2 = rad_sum / DBLE(count_sum)
! 	IF((radius2/dx2 > 0.9D0 * DBLE(length_step_r_2)) .OR. (radius2/dx2 > 0.9D0 * DBLE(length_step_z_2))) THEN
!		vel2_max = 4.0 * vel2_max
!	END IF
!      	! Create the boundary velocity for NM First!
!	vel_frame_r2 (:,:) = 0.0D0
!	vel_frame_z2 (:,:) = 0.0D0
!	DO k = length_step_z_min_part_2, length_step_z_part_2, 1
!		DO j = 1, length_step_r_part_2, 1
!			vel_frame_r2(j,k) = vel2_max*r2(j)/radius2
!			vel_frame_z2(j,k) = vel2_max*z2(k)/radius2
!		ENDDO
!	ENDDO
!	! Copy the results to ghost cells !
!	CALL boundary1D_NM(vel_frame_r2, oddR)
!	CALL boundary1D_NM(vel_frame_z2, oddZ)
!ELSE
!	! Otherwise, zero all of them !
!	vel_frame_r2 (:,:) = 0.0D0
!	vel_frame_z2 (:,:) = 0.0D0
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! My patch !
vel2_max = SQRT(maxval(vel2_r)**2 + maxval(vel2_z)**2)
radius2 = max(r2(length_step_r_2), z2(length_step_z_2))

! Adjust grid velocity !
IF(boundary2 < 0.90D0*max(r2(length_step_r_2),z2(length_step_z_2))) THEN
        vel2_max = vel2_max*0.5D0
ELSE
        vel2_max = vel2_max
END IF

! Create the boundary velocity for NM First!
vel_frame_r2 (:,:) = 0.0D0
vel_frame_z2 (:,:) = 0.0D0

! Assign !
DO k = length_step_z_min_part_2, length_step_z_part_2, 1
	DO j = 1, length_step_r_part_2, 1
		vel_frame_r2(j,k) = vel2_max*r2(j)/radius2
		vel_frame_z2(j,k) = vel2_max*z2(k)/radius2
	ENDDO
ENDDO

! Copy the results to ghost cells !
CALL boundary1D_NM(vel_frame_r2, oddR)
CALL boundary1D_NM(vel_frame_z2, oddZ)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To repeat everything in moving grid algorithm, but for DM !						   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDFRAMEVEL_DM
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: k, j, count_sum

! Real parameter !
REAL (DP) :: rho_max1, rad_sum, rad_dist, rad_vel

! Here I find the suitable velocity for the !
! boundary to expand such that the whole star !
! is always contained inside the box !

! Define the density range to be considered
rho_max1 = 1.0D2 * rho1_a
	
! Searching NM atmosphere speed ! 
! Initialize !
rad_sum = 0.0D0
vel1_max = 0.0D0
count_sum = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do for NM first !
!DO k = length_step_z_min_part_1, length_step_z_part_1, 1
!	DO j = 1, length_step_r_part_1, 1
!		if(rho1(j,k) < rho_max1) then
!			rad_dist = DSQRT(r1(j) ** 2 + z1(k) ** 2)
!			rad_vel = (vel1_r(j,k) * r1(j) + vel1_z(j,k) * z1(k)) / rad_dist
!			rad_sum = rad_sum + rad_dist
!			vel1_max = vel1_max + rad_vel
!			count_sum = count_sum + 1
!		endif
!	enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If the star surface is too close to the box boundary
! Artificially increase the boundary velocity
!IF(count_sum /= 0) THEN
!	vel1_max = vel1_max / DBLE(count_sum)
!	radius1 = rad_sum / DBLE(count_sum)
! 	IF((radius1/dx1 > 0.9D0 * DBLE(length_step_r_1)) .OR. (radius1/dx1 > 0.9D0 * DBLE(length_step_z_1))) THEN
!		vel1_max = 4.0 * vel1_max
!	END IF
!      	! Create the boundary velocity for DM First!
!	vel_frame_r1 (:,:) = 0.0D0
!	vel_frame_z1 (:,:) = 0.0D0
!	DO k = length_step_z_min_part_1, length_step_z_part_1, 1
!		DO j = 1, length_step_r_part_1, 1
!			vel_frame_r1(j,k) = vel1_max*r1(j)/radius1
!			vel_frame_z1(j,k) = vel1_max*z1(k)/radius1
!		ENDDO
!	ENDDO
!	! Copy the results to ghost cells !
!	CALL boundary1D_DM(vel_frame_r1, oddR)
!	CALL boundary1D_DM(vel_frame_z1, oddZ)
!ELSE
!	! Otherwise, zero all of them !
!	vel_frame_r1 (:,:) = 0.0D0
!	vel_frame_z1 (:,:) = 0.0D0
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! My patch !
vel1_max = SQRT(maxval(vel1_r)**2 + maxval(vel1_z)**2)
radius1 = max(r1(length_step_r_1), z1(length_step_z_1))

! Adjust grid velocity !
IF(boundary1 < 0.90D0*max(r1(length_step_r_1),z1(length_step_z_1))) THEN
	vel1_max = vel1_max*0.5D0
ELSE
	vel1_max = vel1_max
END IF

! Create the boundary velocity for DM First!
vel_frame_r1 (:,:) = 0.0D0
vel_frame_z1 (:,:) = 0.0D0

! Assign !
DO k = length_step_z_min_part_1, length_step_z_part_1, 1
	DO j = 1, length_step_r_part_1, 1
		vel_frame_r1(j,k) = vel1_max*r1(j)/radius1
		vel_frame_z1(j,k) = vel1_max*z1(k)/radius1
	ENDDO
ENDDO

! Copy the results to ghost cells !
CALL boundary1D_DM(vel_frame_r1, oddR)
CALL boundary1D_DM(vel_frame_z1, oddZ)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the outermost radius of the star !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDRADIUS_DM
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Temporal integer !
INTEGER :: jr, kz

! Initialize !
jr = 1
kz = 1

! For DM !
DO k = 1, length_step_z_1, 1     
   DO j = 1, length_step_r_1 - 1, 1
      if(rho1(j,k) > rho1_a .and. rho1(j+1,k) <= rho1_a) then
	jr = max(j, jr)
      endif                 
   ENDDO
ENDDO
DO j = 1, length_step_r_1, 1
   DO k = 1, length_step_z_1 - 1, 1
      IF(rho1(j,k) > rho1_a .and. rho1(j,k+1) <= rho1_a) THEN
	kz = max(k, kz)
      ENDIF
   ENDDO
ENDDO
boundary1 = max(r1(jr), z1(kz))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the outermost radius of the star !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDRADIUS_NM
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Temporal integer !
INTEGER :: jr, kz

! Initialize !
jr = 1
kz = 1

! For NM !
DO k = 1, length_step_z_2, 1     
   DO j = 1, length_step_r_2 - 1, 1
      if(rho2(j,k) > rho2_a .and. rho2(j+1,k) <= rho2_a) then
	jr = max(j, jr)
      endif       
   ENDDO
ENDDO
DO j = 1, length_step_r_2, 1
   DO k = 1, length_step_z_2 - 1, 1
      IF(rho2(j,k) > rho2_a .and. rho2(j,k+1) <= rho2_a) THEN
	kz = max(k, kz)
      ENDIF
   ENDDO
ENDDO
boundary2 = max(r2(jr), z2(kz))

END SUBROUTINE