!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine finds the flame propagation velocity by finding 
! the maximum of the three terms
! 1. laminar flame (Timmes 1992)
! 2. Rayleigh-taylor instability (Schmidty2006a)
! 3. turbulent flame (Niemeyer1995)
!
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_flame_velocity(j_in, k_in, flame_localvel)
   use definition
   use levelset_module
   use helmeos_module
   use turb_module
   implicit none

   ! Critical number for DDT transition
   real (DP) :: Ka_crit = 0.5D0

   ! The Gibson length scale
   real (DP) :: Gibson_length
  
   ! The width of flame
   real (DP) :: flame_length

   ! The various flame speed
   real (DP) :: flame_vel_lam, flame_vel_helm, flame_vel_turb

   ! Output effective flame speed
   real (DP) :: flame_localvel

   ! Check if really the level set is changwed
   integer :: modify_scaG_flag

   ! Input grid number
   integer :: j_in, k_in

   ! Dummy variables
   integer :: j, k

   ! Local dummy variables
   real (DP) :: rho_fuel, dist

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First store in local variables
   rho_fuel = rho2(j_in,k_in)

   ! Then find the laminar flame speed
   flame_vel_lam = 3.06666667D-4 * (rho_fuel / 3.24D-9) ** 0.805D0
   !flame_vel_lam = 4.29D-4 * (rho_fuel / 3.24D-9) ** 0.805D0

   ! When there is turbulence, then map to the effective turblent flame speed
   if(turb_flag == 1) then

      ! The traditional one
      flame_vel_turb = flame_vel_lam * DSQRT(1.0D0 + 8.0D0 / 3.0D0 * turb_q(j_in,k_in) / flame_vel_lam ** 2)  

      ! The updated one (but slower)
      !flame_vel_turb = flame_vel_lam * DSQRT(1.0D0 + 1.228D0 * turb_q(j_in,k_in) / flame_vel_lam ** 2)  !1.228, 200.0

      ! The best fitted one
      !flame_vel_turb = flame_vel_lam * (1.0D0 + 0.654D0 * (turb_q(j_in,k_in) / flame_vel_lam**2) ** 0.5985)

      ! The pre-Schmidt one
      !flame_vel_turb = MAX(flame_vel_lam, DSQRT(2.0D0 * turb_q(j_in,k_in)))

   endif

   ! Only flame above density threshold can propagate
   if(rho_fuel > rho2_flame_min .and. rho_fuel < rho2_flame_max) then

      flame_localvel = max(flame_vel_turb, flame_vel_lam) !flame_vel_turb
      !flame_localvel = flame_vel_lam

   else

      flame_localvel = 0.0D0

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! If deton_flag is on, then check together whether the
   ! flame is thick enough to cause DDT
   if(deton_flag == 1) then

      ! Typing DDT required density
      if(rho_fuel <= 3.24D-11 .and. rho_fuel > rho2_deton_min) then

	 ! Characteristic length scale
	 flame_length = 5.1381D-9 * (rho_fuel/1.62D-9)**(-1.9375D0)
         gibson_length = dx2 * (flame_vel_lam**2 / 2.0D0 / turb_q(j_in,k_in))**1.5D0
	
	 ! Flame width > turbulence eddy-overturn scale is needed
	 ! But the exact value is not yet discovered
	 if(flame_length > Ka_crit * gibson_length) then
         !if(flame_vel_turb >= 0.2268D0 * (rho_fuel/1.62D-9) ** 1.4508) then
         !if(turb_q(j_in,k_in) >= 1.0D0**(2.0D0/3.0D0) * 0.01982D0 * (rho_fuel/1.62D-9) ** 8.705D0/3) then
      
	    ! The grid must be only partially burnt
            if(rho_fuel <= rho2_deton_max .and. flamegrid_flag(j_in,k_in) > 1) then

	       ! For the first flame, then switch to 
	       ! the output timescale and initialize
	       ! the related detonation level set
               IF(found_deton_flag == 0) then

		  write(*,*) 'Occured at = ', global_time
                  write(*,*) 'Found deton at ', j_in, k_in
                  write(*,*) 'Flame length = ', flame_length
                  write(*,*) 'Gibson length = ', gibson_length
                  write(*,*)

                  !cfl = 0.2D0
                  output_file = .true.     
                  output_profiletime = 0.5D4
                  output_flametime = 0.5D4
	          output_turbtime = 0.5D4
	          output_Helmtime = 0.5D4
                  found_deton_flag = 1
                  found_deton_time = global_time

                  deton_ratio_old(:,:) = 0.0D0

		  ! Plant the level set
                  do j = 1, length_step_r_2, 1     
                     do k = 1, length_step_z_2, 1
                        dist = DSQRT((r2(j) - r2(j_in))**2 + (z2(k) - z2(k_in))**2) - 10.0D0
                        scaG2(j,k) = -dist
                     enddo
                  enddo
        
                  CALL boundary1D_NMFULL (scaG2,even)
               
		  ! Make sure the distance property of the level set
	 	  ! needs to be preserved always
                  call reinitialization3()               !You need this for aspherical initial flame
		  !call update_scaG2corn()
                  call identify_detongrid()
                  call compute_detonratio()

                  !deton_ratio_old(:,:) = deton_ratio(:,:)

               ELSEIF(found_deton_flag == 1 .and. scaG2(j_in,k_in) < -100.0D0) THEN

		  ! If this is the 2nd time (or above) for the
		  ! detonation to start, then simply plant
		  ! the detonation seed 
		  modify_scaG_flag = 0

		  write(*,*) 'Occured at = ', global_time
                  write(*,*) 'Found deton at ', j_in, k_in
                  write(*,*) 'Flame length = ', flame_length
                  write(*,*) 'Gibson length = ', gibson_length
                  write(*,*)

                  do j = j_in - 3, j_in + 3, 1
                     do k = k_in - 3, k_in + 3, 1
	
			dist = DSQRT((r2(j) - r2(j_in))**2 + (z2(k) - z2(k_in))**2) - 10.0D0 
                        scaG2(j,k) = -dist
			modify_scaG_flag = 1
                       
		     enddo
                  enddo   

	          if(modify_scaG_flag == 1) then

                     CALL boundary1D_NMFULL (scaG2,even)
               
   		     ! Preserve the distance property
                     call reinitialization3()               !You need this for aspherical initial flame
		     !call update_scaG2corn()
                     call identify_detongrid()
                     call compute_detonratio()

		  endif

               endif
         
            endif
      
         endif

      endif
 
   endif

   end subroutine