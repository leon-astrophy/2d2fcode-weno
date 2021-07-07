!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module containing all subroutine that are related to jet explosions due to CCSN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE JetExplosion_module
USE definition
IMPLICIT NONE

! Global arrays and variables !
REAL (DP), DIMENSION(:,:), ALLOCATABLE :: erad
REAL (DP) :: erad_sum
REAL (DP) :: lumino_jet

CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine allocates erad
   !
   ! Written by Leung SHing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE BuildJetExp
   USE DEFINITION
   IMPLICIT NONE

   ! Allocate !
   ALLOCATE(erad(-4:length_step_r+5,-4:length_step_z+5))

   END SUBROUTINE BuildJetExp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine find the total radiation energy
   !
   ! WRitten by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   SUBROUTINE FindEsum
   USE definition
   IMPLICIT NONE

   INTEGER :: j, k

   erad_sum = 0.0D0
   DO k = 1, length_step_z, 1
      DO j = 1, length_step_r, 1

         erad_sum = erad_sum + vol(j,k) * erad(j,k)

      ENDDO
   ENDDO 
   
   END SUBROUTINE FindEsum

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine deallocates erad
   !
   ! Written by Leung SHing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DestroyJetExp 
   USE DEFINITION
   IMPLICIT NONE

   ! Dealloate !
   DEALLOCATE(erad)

   END SUBROUTINE DestroyJetExp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine diffuse the jet energy 
   ! assuming it is from the accretion disk of the BH
   ! inside the massive star
   !
   ! WRITTEN by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SolveJetDiff
   USE DEFINITION
   USE HELMEOS_MODULE
   IMPLICIT NONE

   ! Integer !
   INTEGER :: j, k, kmax 
 
   ! Real variables !
   REAL (DP) :: jet_flux(-4:length_step_r+5)
   REAL (DP) :: total_vol, de
   REAL (DP) :: erad_sum_bf, erad_sum_af

   ! assign !
   kmax = theta_jet / dtheta

   ! Sum to total volume !
   !total_vol = 4.0D0 / 3.0D0 * pi * ((r2(1)+0.5D0*dx)**3 - (r2(1) -0.5D0*dx)**3) * (1.0D0 - COS(theta_jet))    
   total_vol = 0.0D0
   DO k = 1, kmax, 1
      total_vol = total_vol + vol(1,k)
   ENDDO

   !erad_sum_af = 0.0D0
   !DO k = 1, kmax, 1                
   !   DO j = 1, length_step_r, 1
   !      erad_sum_af = erad_sum_af + erad(j,k) * vol(j,k)
   !   ENDDO
   !ENDDO     
   !WRITE(*,*) global_time, erad_sum_af
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Find the jet flux !
   DO k = 1, kmax, 1

      DO j = 1, length_step_r, 1 ! A good approximation
         IF(erad(j,k) > 1.0D-40) then
            jet_flux(j) = erad(j,k) * MIN(dt / dx / sca_fac1(j), 1.0D0)
         ELSE
            jet_flux(j) = 0.0D0
         ENDIF
      ENDDO

      DO j = 1, 1, 1
         IF(global_time < jet_time) then
            de = dedep * dt / total_vol 
         ELSE
            de = 0.0D0
         ENDIF
      ENDDO

      !if(k == 1) then
      !   DO j = 0, 20, 1
      !      write(*,*) erad(j,k), jet_flux(j), de
      !   ENDDO
      !   write(*,*)
      !endif

      DO j = 1, 1, 1
         erad(j,k) = MAX(erad(j,k) - (jet_flux(j) * r2(j)**2 * sca_fac1(j) - jet_flux(j-1) * r2(j-1)**2 * sca_fac1(j-1)) / &
                                      r2(j)**2 / sca_fac1(j) + de, 1.0D-40)
      ENDDO
      DO j = 2, length_step_r, 1
         erad(j,k) = MAX(erad(j,k) - (jet_flux(j) * r2(j)**2 * sca_fac1(j) - jet_flux(j-1) * r2(j-1)**2 * sca_fac1(j-1)) / &
                                      r2(j)**2 / sca_fac1(j), 1.0D-40)
      ENDDO

   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Check the total amount of energy in the radiation field
   !erad_sum_af = 0.0D0              
   !DO k = 1, kmax, 1
   !   DO j = 1, length_step_r, 1     
   !      erad_sum_af = erad_sum_af + erad(j,k) * vol(j,k)
   !   ENDDO        
   !ENDDO        
   !WRITE(*,*) global_time, erad_sum_af

   END SUBROUTINE SolveJetDiff
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This code solves the energy deposition by the jet
   ! This assumes radiation absorption by Compton process
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

   SUBROUTINE SOlveJetDep
   USE definition
   USE helmeos_module
   IMPLICIT NONE

   ! Integer !
   INTEGER :: j, k, kmax

   ! Reals !
   REAL (DP) :: erad_new
   REAL (DP) :: total_vol, de
   REAL (DP) :: erad_sum_bf, erad_sum_af

   ! Assign !
   kmax = theta_jet / dtheta
   lumino_jet = 0.0D0

   !erad_sum_af = 0.0D0
   !DO k = 1, kmax, 1
   !   DO j = 1, length_step_r, 1
   !      erad_sum_af = erad_sum_af + erad(j,k) * vol(j,k)
   !   ENDDO
   !ENDDO
   !WRITE(*,*) global_time, erad_sum_af
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! Sum them all !
   DO k = 1, kmax, 1
      DO j = 1, length_step_r, 1

         IF(rho2(j,k) > 1.0D-14) THEN
            erad_new = MAX(erad(j,k) * EXP(-rho2(j,k) * kap_e * dt), 1.0D-40)
            epsilon2(j,k) = epsilon2(j,k) - (erad_new - erad(j,k)) / rho2(j,k)
            lumino_jet = lumino_jet - (erad_new - erad(j,k)) * VOL(j,k)
            erad(j,k) = erad_new
         ENDIF

      ENDDO
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign, get the rate by diving through dt !
   lumino_jet = lumino_jet / dt   

   CALL findhelmtemp
   CALL boundary1D_NM(erad, EVEN)
   CALL boundary1D_NM(epsilon2, EVEN)

   !erad_sum_af = 0.0D0
   !DO k = 1, kmax, 1
   !   DO j = 1, length_step_r, 1
   !      erad_sum_af = erad_sum_af + erad(j,k) * vol(j,k)
   !   ENDDO
   !ENDDO
   !WRITE(*,*) global_time, erad_sum_af

   END SUBROUTINE SolveJetDep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine clean all the residue grid which
   ! remain non-zero after the jet has passed due
   ! to gravity
   !
   ! Written by Leung Shing Chi in 2017
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   SUBROUTINE cleaninnerzone
   USE definition
   USE helmeos_module
   IMPLICIT NONE

   ! Integer !
   INTEGER :: j, j2, k, kmax, j_in1, j_in2, j_out1, j_out2, j_clean
   INTEGER :: grid_type

   ! Logical flag 
   LOGICAL :: found_flag, found_flag_last, root_flag
   INTEGER :: cont_flag(-4:length_step_r+5, -4:length_step_z+5)

   ! Others !
   REAL (DP) :: rho_min
   INTEGER :: kmax2

   ! Assign !
   kmax = INT(theta_jet / dtheta)
   rho_min = 1.1D0 * rho2_a

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Now try the simple one,
   ! if there is negihbour from the trunk, keep ot

   cont_flag(:,:) = 0
   cont_flag(length_step_r-100:length_step_r,:) = 1

   DO j = length_step_r-101, 1, -1
      DO k = 1, length_step_z, 1
         IF((cont_flag(j+1,k) == 1 .or. cont_flag(j,k-1) == 1) .and. rho2(j,k) > rho_min) then
            cont_flag(j,k) = 1
         endif
      ENDDO
   ENDDO

   DO j = length_step_r-101, 1, -1
      DO k = length_step_z, 1, -1
         if((cont_flag(j+1,k) == 1 .or. cont_flag(j,k+1) == 1) .and. rho2(j,k) > rho_min) then
            cont_flag(j,k) = 1
         endif
      ENDDO
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Clena the grid with cont_flag = 0
   DO k = 1, length_step_z, 1
      DO j = 1, length_step_r, 1

         IF(cont_flag(j,k) == 0) then
            !write(*,*) 'Clean island at ', global_time, j, k
            rho2 (j,k) = rho2_a
            vel2_r (j,k) = vel2_a
            vel2_z (j,k) = vel2_a
            vel2_p (j,k) = vel2_a
            epsilon2 (j,k) = epsilon2_a
            temp2 (j,k) = temp_a  
            ye2 (j,k) = ye_a
         ENDIF

      ENDDO
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL boundary1D_NM(rho2, even)
   CALL boundary1D_NM(vel2_r, oddR)
   CALL boundary1D_NM(vel2_z, oddZ)
   CALL boundary1D_NM(vel2_p, oddR)
   CALL boundary1D_NM(epsilon2, even)
   CALL boundary1D_NM(temp2, even)
   CALL boundary1D_NM(ye2, even)

   END SUBROUTINE cleaninnerzone

END MODULE JetExplosion_module