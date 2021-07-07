!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all the tools for calculating the 
! gravitational wave emission and its associated energy loss
! based on the quadruple model
! Written by Leung Shing Chi in 2016
!
! This module contains the following subroutines
! 1. subroutine getGravQ
! 2. subroutine findgravQ
! 3. subroutine findGravLum
! 4. subroutine outputGravQ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE GW_MODULE
USE DEFINITION, ONLY : DP, RUNDM_flag
IMPLICIT NONE

! The quadrupole contribution !
REAL (DP) :: ae220_1, ae220_2

! The plus-direction for the TT mode
REAL (DP) :: hTT_plus1, hTT_plus2

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine find the mass moment of inertia
   ! Notice that the non-diagonal terms are 
   ! automatically zero when there is rotation symmettry
   !
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE FINDGRAVQ
   USE DEFINITION
   IMPLICIT NONE

   ! Definition of 1kpc !
   REAL (DP) :: onekpc = 2.089599257D16

   ! Dummy variables
   integer :: i, j, k

   ! Dummy variables !
   REAL (DP) :: sum1, sum2

   ! Radial velocity !
   REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vr1, vr2

   ! Polar velocity !
   REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vtheta1, vtheta2

   ! Azimuthal velocity !
   REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vphi1, vphi2

   ! Radial potential derivative !
   REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dphidr1, dphidr2

   ! Polar potential derivative !
   REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dphidth1, dphidth2

   ! Integrand for calculating G-wave !
   REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: int1, int2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Allocate !
   ALLOCATE (vr2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   ALLOCATE (vphi2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   ALLOCATE (vtheta2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   ALLOCATE (dphidr2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   ALLOCATE (dphidth2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   ALLOCATE (int2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))

   ! For DM !
   IF(RUNDM_flag == 1) THEN
   	ALLOCATE (vr1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
   	ALLOCATE (vphi1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
   	ALLOCATE (vtheta1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (dphidr1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	ALLOCATE (dphidth1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
        ALLOCATE (int1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize !
   vr2 = 0.0D0
   vtheta2 = 0.0D0
   vphi2 = 0.0D0

   ! First assign the radial and polar velocity !
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

	 ! Atmosphere is not counted
         if(rho2(j,k) > rho2_a) then
            vr2(j,k) = (r2(j)*vel2_r(j,k) + z2(k)*vel2_z(j,k))/rad2(j,k)
	    vtheta2(j,k) = (z2(k)*vel2_r(j,k) - r2(j)*vel2_z(j,k))/rad2(j,k)
	    vphi2(j,k) = vel2_p(j,k)
         endif

      enddo
   enddo

   ! DO for the DM !
   IF(RUNDM_flag == 1) THEN
      ! Initialize !
      vr1 = 0.0D0
      vtheta1 = 0.0D0
      vphi1 = 0.0D0

      do k = length_step_z_min_part_1, length_step_z_part_1, 1
      	do j = 1, length_step_r_part_1, 1

	 ! Atmosphere is not counted
         if(rho1(j,k) > rho1_a) then
            vr1(j,k) = (r1(j)*vel1_r(j,k) + z1(k)*vel1_z(j,k))/rad1(j,k)
	    vtheta1(j,k) = (z1(k)*vel1_r(j,k) - r1(j)*vel1_z(j,k))/rad1(j,k)
	    vphi1(j,k) = vel1_p(j,k)
         endif

        enddo
      enddo
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize !
   dphidr2 = 0.0D0
   dphidth2 = 0.0D0

   ! Then assign the radial and polar potential derivative !
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

	 ! Atmosphere is not counted
         if(rho2(j,k) > rho2_a) then
            dphidr2(j,k) = (r2(j)*phi2_r(j,k) + z2(k)*phi2_z(j,k))/rad2(j,k)
	    dphidth2(j,k) = (z2(k)*phi2_r(j,k) - r2(j)*phi2_z(j,k))
         endif

      enddo
   enddo

   ! DO for the DM !
   IF(RUNDM_flag == 1) THEN

      ! Initialize !
      dphidr1 = 0.0D0
      dphidth1 = 0.0D0

      do k = length_step_z_min_part_1, length_step_z_part_1, 1
      	do j = 1, length_step_r_part_1, 1

	 ! Atmosphere is not counted
         if(rho1(j,k) > rho1_a) then
            dphidr1(j,k) = (r1(j)*phi1_r(j,k) + z1(k)*phi1_z(j,k))/rad1(j,k)
	    dphidth1(j,k) = (z1(k)*phi1_r(j,k) - r1(j)*phi1_z(j,k))
         endif

        enddo
      enddo
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize !
   int2 = 0.0D0

   ! Finally get the integrand !
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

	 ! Atmosphere is not counted
         if(rho2(j,k) > rho2_a) then
            int2(j,k) = rho2(j,k)*((vr2(j,k)**2 - rad2(j,k)*dphidr2(j,k))*(3.0D0*cos2(j,k)**2 - 1.0D0) + &
			(3.0D0*dphidth2(j,k) - 6.0D0*vr2(j,k)*vtheta2(j,k))*cos2(j,k)*sqrt(1.0D0 - cos2(j,k)**2) + &
			vtheta2(j,k)**2*(2.0D0 - 3.0D0*cos2(j,k)**2) - vphi2(j,k)**2)
         endif

      enddo
   enddo

   ! DO for the DM !
   IF(RUNDM_flag == 1) THEN

      ! Initialize !
      int1 = 0.0D0

      do k = length_step_z_min_part_1, length_step_z_part_1, 1
      	do j = 1, length_step_r_part_1, 1

	 ! Atmosphere is not counted
         if(rho1(j,k) > rho1_a) then
            int1(j,k) = rho1(j,k)*((vr1(j,k)**2 - rad1(j,k)*dphidr1(j,k))*(3.0D0*cos1(j,k)**2 - 1.0D0) + &
			(3.0D0*dphidth1(j,k) - 6.0D0*vr1(j,k)*vtheta1(j,k))*cos1(j,k)*sqrt(1.0D0 - cos1(j,k)**2) + &
			vtheta1(j,k)**2*(2.0D0 - 3.0D0*cos1(j,k)**2) - vphi1(j,k)**2)
         endif

        enddo
      enddo
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize !
   sum2 = 0.0D0

   ! Do the integration !
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

	 ! Atmosphere is not counted
         if(rho2(j,k) > rho2_a) then
            sum2 = sum2 + int2(j,k)*vol2(j,k)
         endif

      enddo
   enddo

   ! DO for the DM !
   IF(RUNDM_flag == 1) THEN

      ! Initialize !
      sum1 = 0.0D0

      do k = length_step_z_min_part_1, length_step_z_part_1, 1
      	do j = 1, length_step_r_part_1, 1

	 ! Atmosphere is not counted
         if(rho1(j,k) > rho1_a) then
            sum1 = sum1 + int1(j,k)*vol1(j,k)
         endif

        enddo
      enddo
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign the value !
   ae220_2 = 16.0D0*pi**(1.5D0)/sqrt(15.0D0)*sum2/(2.0D0*pi)
   hTT_plus2 = 0.5D0*(0.125D0*sqrt(15.0D0/pi)*ae220_2/onekpc)

   ! Do for DM !
   IF(RUNDM_flag == 1) THEN
   	ae220_1 = 16.0D0*pi**(1.5D0)/sqrt(15.0D0)*sum1/(2.0D0*pi)
   	hTT_plus1 = 0.5D0*(0.125D0*sqrt(15.0D0/pi)*ae220_1/onekpc)
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\

   ! Allocate !
   DEALLOCATE (vr2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   DEALLOCATE (vphi2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   DEALLOCATE (vtheta2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   DEALLOCATE (dphidr2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   DEALLOCATE (dphidth2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))
   DEALLOCATE (int2(-4 : length_step_r_2 + 5, -4 : length_step_z_2 + 5))

   ! For DM !
   IF(RUNDM_flag == 1) THEN
   	DEALLOCATE (vr1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
   	DEALLOCATE (vphi1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
   	DEALLOCATE (vtheta1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	DEALLOCATE (dphidr1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
	DEALLOCATE (dphidth1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
        DEALLOCATE (int1(-4 : length_step_r_1 + 5, -4 : length_step_z_1 + 5))
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   END SUBROUTINE findgravQ

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 
   ! This subroutine output the results of the above
   ! calculation.
   ! Written by Leung Shing CHi
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE outputGravQ
   USE DEFINITION, only : global_time
   IMPLICIT NONE
 
   ! Open, write and close !
   OPEN(unit=601, file='./Outfile/GWave/Star_WENO_GravQ_0_NM.dat', action='write', position='append')
   WRITE(601, 100) global_time, ae220_2, hTT_plus2
   CLOSE(601)

   ! Open, write and close !
   IF(RUNDM_flag == 1) THEN
   	OPEN(unit=601, file='./Outfile/GWave/Star_WENO_GravQ_0_DM.dat', action='write', position='append')
   	WRITE(601, 100) global_time, ae220_1, hTT_plus1
   	CLOSE(601)
   END IF

   100 format(5ES18.8)
   END SUBROUTINE outputGravQ

END MODULE GW_module