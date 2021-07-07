!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine estimates the initial gravitational potential
! and pass to the POTENTIALRELAX subroutine 
! Written by Leung Shing Chi in 2016
! The subroutine asumes spherical distribution of matter
! and gives a monopole estimation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDPOTENTIAL(n)
use definition
implicit none

! Input: The find-gravity flag
INTEGER :: n

! Dummy variables
INTEGER :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug stuff !
!IF(debug_flag == 1) WRITE(*,*) 'In Find Potential'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (w_gravity_i == 1) THEN

   ! FInd the mass to match the boundary condition
   CALL FINDMASS

   ! Find potential !
   CALL MULTIPOLE

   ! Set up initial trial potential
   if(potential_flag == 0) then

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Based on the above group, a trial potential (in spherical symmetry) is posed

      ! Guess the 1D corresponding grid-number
      DO k = 1, length_step_z_2, 1
	 DO j = 1, length_step_r_2, 1
	    phi2_nm(j, k) = - (mass1 + mass2) / rad2(j, k)
         ENDDO
      ENDDO

      ! Copy to ghost cells
      CALL BOUNDARY1D_NMFULL(phi2_nm, even)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! For DM 
      IF(DM_flag == 1) THEN
      	DO k = 1, length_step_z_1, 1
	 DO j = 1, length_step_r_1, 1
	    phi1_dm(j, k) = - (mass1 + mass2) / rad1(j, k)
         ENDDO
      	ENDDO
	CALL BOUNDARY1D_DMFULL(phi1_dm, even)
      END IF

      ! The initial guess only do one time!
      potential_flag = 1

   ENDIF
 
   ! Go to successive over-relaxation
   CALL RELAXATION_NM()

   ! For DM !
   IF(DM_flag == 1) THEN
   	CALL RELAXATION_DM()
   END IF

   ! Map the potential !
   CALL INTERPOLATION

ELSEIF (initmodel_flag == 0 .AND. testmodel_flag == 5) THEN

	! A constant gravity = -0.1 is needed !
	! For RT instability test !
	DO k = -4, length_step_z_2 + 5
		DO j = -4, length_step_r_2 + 5
			phi2 (j,k) = 0.1E0_DP * z2(k)
			phi2_dm(j, k) = 0.0D0
			phi2_nm(j, k) = 0.0D0
			phi2_r (j,k) = 0.0E0_DP
			phi2_z (j,k) = 0.1E0_DP
		ENDDO
	END DO
	
ELSE
	
   ! Gravity flag is off, so no need
   ! to find any gravitational potential
   DO k = -4, length_step_z_2 + 5
      DO j = -4, length_step_r_2 + 5
	 phi2 (j,k) = 0.0E0_DP
	 phi2_dm(j, k) = 0.0D0
	 phi2_nm(j, k) = 0.0D0
	 phi2_r (j,k) = 0.0D0
	 phi2_z (j,k) = 0.0D0
      ENDDO
   END DO

   ! For DM !
   IF(DM_flag == 1) THEN
   	DO k = -4, length_step_z_1 + 5
   	   DO j = -4, length_step_r_1 + 5
	     phi1 (j,k) = 0.0E0_DP
	     phi1_dm(j, k) = 0.0D0
	     phi1_nm(j, k) = 0.0D0
	     phi1_r (j,k) = 0.0D0
	     phi1_z (j,k) = 0.0D0
          ENDDO
       END DO
   END IF

ENDIF

END SUBROUTINE findpotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the legendre polynominal at each grid !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDLEGENDRE
USE DEFINITION 
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, l

! Find all the legendre function !
DO j = 1, length_step_r_2
	DO k = 1, length_step_z_2
		DO l = 0, 2*LMAX
			CALL LEGENDRE(legendre2(j,k,l), cos2(j,k), l)
		END DO
	END DO
END DO

! For DM Grid !
IF(DM_flag == 1) THEN
	DO j = 1, length_step_r_1
		DO k = 1, length_step_z_1
			DO l = 0, 2*LMAX
				CALL LEGENDRE(legendre1(j,k,l), cos1(j,k), l)
			END DO
		END DO
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the legendre function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE LEGENDRE(pout, x, in)
USE DEFINITION 
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: in

! Input real number !
REAL (DP), INTENT(IN) :: x

! Output polynominal !
REAL (DP), INTENT(OUT) :: pout

! Dummy integer !
INTEGER :: j

! legender function !
REAL (DP), DIMENSION(0:2*LMAX) :: pnx

! Assign !
pnx(0) = 1.0D0
pnx(1) = x

! For higher legender function !
DO j = 2, in
	pnx(j) = ((2.0D0*DBLE(j) - 1.0D0)*pnx(j-1)*x - (DBLE(j) - 1.0D0)*pnx(j-2))/DBLE(j)
END DO

pout = pnx(in)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate multipole moment of 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MULTIPOLE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Initialize !
qpole1 = 0.0D0
qpole2 = 0.0D0

! Sum over for NM !
DO l = 0, LMAX
	
	! Integrate z-direction !
	DO j = 1, length_step_r_2
		DO k = 1, length_step_z_2
			IF(rho2(j,k) > rho2_a) THEN
				qpole2(2*l) = qpole2(2*l) + legendre2(j,k,2*l)*rho2(j,k)*radbar2(j,k)**(2*l)*volbar2(j,k)
			ELSE
				CYCLE
			END IF
		END DO
	END DO
END DO

! Do for DM !
IF(RUNDM_flag == 1) THEN
	DO l = 0, LMAX
		DO j = 1, length_step_r_1
			DO k = 1, length_step_z_1
				IF(rho1(j,k) > rho1_a) THEN
					qpole1(2*l) = qpole1(2*l) + legendre1(j,k,2*l)*rho1(j,k)*radbar1(j,k)**(2*l)*volbar1(j,k)
				ELSE
					CYCLE
				END IF
			END DO
		END DO
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Bilinear interpolation of potentials between DM and NM grid 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INTERPOLATION
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, m, n, l

! Do for the DM case !
IF(DM_flag == 1) THEN

	! Initialize !
	phi2_dm = 0.0D0

	! First map from DM grid to NM !
	DO j = 1, length_step_r_2
		DO k = 1, length_step_z_2 
			IF(r2(j) > r1(length_step_r_1) .OR. z2(k) > z1(length_step_z_1)) THEN
				DO l = 0, lmax
					phi2_dm(j,k) = phi2_dm(j,k) + (qpole1(2*l)*legendre2(j,k,2*l)/radbar2(j,k)**(2*l+1))*(rmax1/rmax2)**(2*l)
				END DO
    				phi2_dm(j,k) = phi2_dm(j,k) * (-1.0D0) * (rmax1**3/rmax2)
   			ELSE
				CALL BINARY_R (j, m, 2)
				CALL BINARY_Z (k, n, 2)
     				IF(r1(m) == r2(j)) THEN
       					IF(z1(n) == z2(k)) THEN
        					phi2_dm(j,k) = phi1_dm(m,n)
       					ELSE
						CALL LINEAR(z1(n-1), z1(n), phi1_dm(m,n-1), phi1_dm(m,n), z2(k), phi2_dm(j,k))
       					END IF
     				ELSE
       					IF(z1(n) == z2(k)) THEN
						CALL LINEAR(r1(m-1), r1(m), phi1_dm(m-1,n), phi1_dm(m,n), r2(j), phi2_dm(j,k))
       					ELSE
						CALL BILINEAR(r1(m-1), r1(m), z1(n-1), z1(n), phi1_dm(m-1,n-1), phi1_dm(m,n-1), phi1_dm(m-1,n), phi1_dm(m,n), r2(j), z2(k), phi2_dm(j,k))
       					END IF
     				END IF
   			END IF
  		END DO
	END DO

	! Initialize !
	phi1_nm = 0.0D0

	! Then we map for NM grid to DM !
	DO j = 1, length_step_r_1
		DO k = 1, length_step_z_1 
			IF(r1(j) > r2(length_step_r_2) .OR. z1(k) > z2(length_step_z_2)) THEN
				DO l = 0, lmax
					phi1_nm(j,k) = phi1_nm(j,k) + (qpole2(2*l)*legendre1(j,k,2*l)/radbar1(j,k)**(2*l+1))*(rmax2/rmax1)**(2*l)
				END DO
    				phi1_nm(j,k) = phi1_nm(j,k) * (-1.0D0) * (rmax2**3/rmax1)
   			ELSE
				CALL BINARY_R (j, m, 1)
				CALL BINARY_Z (k, n, 1)
     				IF(r2(m) == r1(j)) THEN
       					IF(z2(n) == z1(k)) THEN
        					phi1_nm(j,k) = phi2_nm(m,n)
       					ELSEIF(z2(n) > z1(k)) THEN
						CALL LINEAR(z2(n-1), z2(n), phi2_nm(m,n-1), phi2_nm(m,n), z1(k), phi1_nm(j,k))
       					END IF
     				ELSEIF(r2(m) > r1(j)) THEN
       					IF(z2(n) == z1(k)) THEN
						CALL LINEAR(r2(m-1), r2(m), phi2_nm(m-1,n), phi2_nm(m,n), r1(j), phi1_nm(j,k))
       					ELSEIF(z2(n) > z1(k)) THEN
						CALL BILINEAR(r2(m-1), r2(m), z2(n-1), z2(n), phi2_nm(m-1,n-1), phi2_nm(m,n-1), phi2_nm(m-1,n), phi2_nm(m,n), r1(j), z1(k), phi1_nm(j,k))
       					END IF
     				END IF
   			END IF
  		END DO
	END DO

	! Sum them up !
	phi1 = phi1_nm + phi1_dm
	phi2 = phi2_nm + phi2_dm

	! Copy to boundaries !
	CALL BOUNDARY1D_DMFULL(phi1, even)
	CALL BOUNDARY1D_NMFULL(phi2, even)

ELSE

	! Sum them up !
	phi2 = phi2_nm

	! Copy to boundaries !
	CALL BOUNDARY1D_NMFULL(phi2, even)

END IF

! Then we compute the gravitational force !
DO k = 1, length_step_z_2, 1
   DO j = 1, length_step_r_2, 1
      phi2_r (j, k) = (- phi2(j + 2, k) + 8.0E0_DP * phi2 (j + 1, k) - &
                        8.0E0_DP * phi2 (j - 1, k) + phi2 (j - 2, k)) / (1.2E1_DP * dx2)
   ENDDO
ENDDO
DO k = 1, length_step_z_2, 1
   DO j = 1, length_step_r_2, 1
      phi2_z (j, k) = (- phi2 (j, k + 2) + 8.0E0_DP * phi2 (j, k + 1) - &
                        8.0E0_DP * phi2 (j, k - 1) + phi2 (j, k - 2)) / (1.2E1_DP * dx2)
   ENDDO
ENDDO

! Copy to boundaries !
CALL BOUNDARY1D_NMFULL (phi2_r, oddR)                  
CALL BOUNDARY1D_NMFULL (phi2_z, oddZ) 

! DO for DM !
IF(RUNDM_flag == 1) THEN
	DO k = 1, length_step_z_1, 1
   	 DO j = 1, length_step_r_1, 1
      	  phi1_r (j, k) = (- phi1(j + 2, k) + 8.0E0_DP * phi1 (j + 1, k) - &
                        8.0E0_DP * phi1 (j - 1, k) + phi1 (j - 2, k)) / (1.2E1_DP * dx1)
   	 ENDDO
	ENDDO
	DO k = 1, length_step_z_1, 1
   	 DO j = 1, length_step_r_1, 1
      	  phi1_z (j, k) = (- phi1 (j, k + 2) + 8.0E0_DP * phi1 (j, k + 1) - &
                        8.0E0_DP * phi1 (j, k - 1) + phi1 (j, k - 2)) / (1.2E1_DP * dx1)
   	 ENDDO
	ENDDO

	! Copy to boundaries !
	CALL BOUNDARY1D_DMFULL (phi1_r, oddR)                  
	CALL BOUNDARY1D_DMFULL (phi1_z, oddZ) 
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Binary search for position table in R direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BINARY_R (index, out, type)
USE DEFINITION
IMPLICIT NONE

! Integer input !
INTEGER, INTENT(IN) :: index, type
INTEGER, INTENT(OUT) :: out

! Integer !
INTEGER :: left, right, m

! Binary search, case by case !
IF(type == 2) THEN
	left = 1
	right = length_step_r_1
	DO
		IF(left > right) THEN
			IF(r1(m) > r2(index)) THEN
				m = m
			ELSEIF(r1(m) < r2(index)) THEN
				m = m + 1
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(r1(m) < r2(index)) THEN
			left = m + 1	
		ELSEIF(r1(m) > r2(index)) THEN
			right = m - 1
		ELSEIF(r1(m) == r2(index)) THEN
			EXIT
		END IF
	END DO

	! Output !
	out = m
ELSEIF(type == 1) THEN
	left = 1
	right = length_step_r_2
	DO
		IF(left > right) THEN
			IF(r2(m) > r1(index)) THEN
				m = m
			ELSEIF(r2(m) < r1(index)) THEN
				m = m + 1
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(r2(m) < r1(index)) THEN
			left = m + 1	
		ELSEIF(r2(m) > r1(index)) THEN
			right = m - 1
		ELSEIF(r2(m) == r1(index)) THEN
			EXIT
		END IF
	END DO

	! Output !
	out = m
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Binary search for position table in Z direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BINARY_Z (index, out, type)
USE DEFINITION
IMPLICIT NONE

! Integer input !
INTEGER, INTENT(IN) :: index, type
INTEGER, INTENT(OUT) :: out

! Integer !
INTEGER :: left, right, m

! Binary search, case by case !
IF(type == 2) THEN
	left = 1
	right = length_step_z_1
	DO
		IF(left > right) THEN
			IF(z1(m) > z2(index)) THEN
				m = m
			ELSEIF(z1(m) < z2(index)) THEN
				m = m + 1
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(z1(m) < z2(index)) THEN
			left = m + 1	
		ELSEIF(z1(m) > z2(index)) THEN
			right = m - 1
		ELSEIF(z1(m) == z2(index)) THEN
			EXIT
		END IF
	END DO

	! Output !
	out = m
ELSEIF(type == 1) THEN
	left = 1
	right = length_step_z_2
	DO
		IF(left > right) THEN
			IF(z2(m) > z1(index)) THEN
				m = m
			ELSEIF(z2(m) < z1(index)) THEN
				m = m + 1
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(z2(m) < z1(index)) THEN
			left = m + 1	
		ELSEIF(z2(m) > z1(index)) THEN
			right = m - 1
		ELSEIF(z2(m) == z1(index)) THEN
			EXIT
		END IF
	END DO

	! Output !
	out = m
END IF

END SUBROUTINE