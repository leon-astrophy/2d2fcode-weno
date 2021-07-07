!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine to report all current data at some position for debug and checking 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE REPORT_CURRENT_DATA(J,K)
USE DEFINITION
USE LEVELSET_MODULE
USE HELMEOS_MODULE
IMPLICIT NONE

! Grid number of the grid to be reported
INTEGER :: j, k

! Dummy variable 
INTEGER :: i

WRITE(*,*)
WRITE(*,*) '---------------------------------------------------'
WRITE(*,*) 'Report all data for NM of grid ', j, k
WRITE(*,*) '---------------------------------------------------'
WRITE(*,*)
WRITE(*,*) 'Time step =', dt
WRITE(*,*) 'R position', r2(j) 
WRITE(*,*) 'z position', z2(k)
WRITE(*,*) 'Volume = ', vol2(j,k)
WRITE(*,*) 'sca2_fac1', sca2_fac1(j)
WRITE(*,*) 'sca2_fac2', sca2_fac2(j)
WRITE(*,*)
WRITE(*,*) 'Density = ', rho2(j,k)
WRITE(*,*) 'Velocity_r = ', vel2_r(j,k)
WRITE(*,*) 'Velocity_z = ', vel2_z(j,k)
WRITE(*,*) 'Velocity_p = ', vel2_p(j,k)
WRITE(*,*) 'Pressure = ', p2(j,k)
WRITE(*,*) 'Epsilon = ', epsilon2(j,k)
IF (dual_energy == 1) THEN
	WRITE(*,*) 'Energy per unit volume = ', rhoe2(j,k)
END IF
IF (helmeos_flag == 1) THEN
	WRITE(*,*) 'Temperature = ', temp2(j,k)
END IF
WRITE(*,*) 'Ye = ', ye2(j,k)
IF(RUNDM_flag == 1) THEN
	WRITE(*,*) 
	WRITE(*,*) '---------------------------------------------------'
	WRITE(*,*) 'Report all data for DM of grid ', j, k
	WRITE(*,*) '---------------------------------------------------'
	WRITE(*,*)
	WRITE(*,*) 'Time step =', dt
	WRITE(*,*) 'R position', r1(j) 
	WRITE(*,*) 'z position', z1(k)
	WRITE(*,*) 'Volume = ', vol1(j,k)
	WRITE(*,*) 'sca1_fac1', sca1_fac1(j)
	WRITE(*,*) 'sca1_fac2', sca1_fac2(j)
	WRITE(*,*)
	WRITE(*,*) 'Density = ', rho1(j,k)
	WRITE(*,*) 'Velocity_r = ', vel1_r(j,k)
	WRITE(*,*) 'Velocity_z = ', vel1_z(j,k)
	WRITE(*,*) 'Velocity_p = ', vel1_p(j,k)
	WRITE(*,*) 'Pressure = ', p1(j,k)
	WRITE(*,*) 'Epsilon = ', epsilon1(j,k)
END IF

IF(Helmeos_flag == 1) THEN
   WRITE(*,*) 'Fusion rate = ', flame_qdot(j,k)
   WRITE(*,*) 'nse_flag = ', nse_flag(j,k)
   WRITE(*,*)
ENDIF

IF(Flame_flag == 1) THEN
   WRITE(*,*) 'Flame ratio = ', flame_ratio(j,k)
   WRITE(*,*) 'Deton ratio = ', deton_ratio(j,k)
   WRITE(*,*) 'Burn ratio = ', burn_ratio(j,k)
   WRITE(*,*)
ENDIF

IF(burn_prog_flag == 1) THEN
   WRITE(*,*) 'Burn phi1 = ', burn_phi1(j,k)
   WRITE(*,*) 'Burn phi2 = ', burn_phi2(j,k)
   WRITE(*,*) 'Burn phi3 = ', burn_phi3(j,k)
   WRITE(*,*) 'Qash = ', qash(j,k)
   WRITE(*,*) 'Yiso = ', yiso(j,k)
   WRITE(*,*)
ENDIF

IF(helmeos_flag == 1) THEN
   WRITE(*,*) 'Isotope distribution'

   DO i = 1, totalion, 1
      if(xiso(j,k,i) > 1.0D-4) write(*,*) ionam(i), xiso(j,k,i)
   ENDDO

   WRITE(*,*) 'Abar = ', abar2(j,k)
   WRITE(*,*) 'Zbar = ', zbar2(j,k) 
ENDIF

END SUBROUTINE report_current_data