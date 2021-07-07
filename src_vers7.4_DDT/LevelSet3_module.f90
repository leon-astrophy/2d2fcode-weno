!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This module contains all the tools for the level-set
! Notice that the level set is designed for general 
! purpose as long as there is some boundary for mixed-phase
! flow. The baryonic physics are encoded in FindFlameVel
! and in FlameTable_module.
! Written by Leung Shing Chi in 2016
! 
! This module contains the following subroutines:
! 1.  subroutine buildLevelSet
! 2.  subroutine destroyLevelset
! 3.  subroutine outputLevelset_log
! 4.  subroutine output_LevelSet
! 5.  subroutine update_flame_radius
! 6.  subroutine GetFlame
! 7.  subroutine GetDeton
! 8.  subroutine identify_flamegrid
! 9.  subroutine identify_detongrid
! 10. subroutine compute_flameratio
! 11. subroutine compute_detonratio
! 12. subroutine update_scaG2
! 13. subroutine update_scaG3
! 14. subroutine reinitialization2
! 15. subroutine reinitialization3
! 16. subroutine locate_flame_grid
! 17. subroutine locate_min_flame_distance
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module levelset_module
use definition
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Call-Code for matching with WENO scheme
integer :: iscaG1, iscaG2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! New variables for LSM used when fusion_flag == 1 and flame_flag == 1
REAL (DP) :: flame_rad
INTEGER :: flame_grid      
REAL (DP) :: flame_vel_coeff = 0.03D0

! Flag for finding detonation in the simulation
INTEGER :: found_deton_flag = 0
REAL (DP) :: found_deton_time = 0.0D0

! The type of grid depending on the geometry
integer, allocatable, dimension(:,:) :: flamegrid_flag, detongrid_flag
integer, allocatable, dimension(:,:) :: flamecorn_flag, detoncorn_flag

! The level sets
real (DP), allocatable, dimension(:,:) :: scaG, scaG2

! The fraction occupied by the level-set (1st)                                        
real (DP), allocatable, dimension(:,:) :: flame_ratio, flame_ratio_old
real (DP), allocatable, dimension(:,:) :: flame_loc_ratio
                               
! The fraction occupied by the level-set (2nd)           
real (DP), allocatable, dimension(:,:) :: deton_ratio, deton_ratio_old
real (DP), allocatable, dimension(:,:) :: deton_loc_ratio

! Sum of fractions by level-set 
real (DP), allocatable, dimension(:,:) :: burn_ratio

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine allocates the neccessary arrays for 
   ! modeling level-sets.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine buildLevelSet
   use definition
   implicit none

   ! First allocate integer-type array
   allocate(flamegrid_flag(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(detongrid_flag(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(flamecorn_flag(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(detoncorn_flag(-4:length_step_r_2+5, -4:length_step_z_2+5))

   ! Then, allocate real-type array
   allocate(scaG(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(scaG2(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(flame_ratio(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(deton_ratio(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(flame_ratio_old(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(deton_ratio_old(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(burn_ratio(-4:length_step_r_2+5, -4:length_step_z_2+5))

   ! Extra array for NSE in HELMEOS_MODULE
   allocate(flame_loc_ratio(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(deton_loc_ratio(-4:length_step_r_2+5, -4:length_step_z_2+5))

   end subroutine buildLevelSet

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine deallocate all array to clean
   ! up memory
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine destroyLevelset
   implicit none

   ! First deallocate integer-type array
   deallocate(flamegrid_flag)
   deallocate(detongrid_flag)
   deallocate(flamecorn_flag)
   deallocate(detoncorn_flag)
  
   ! Second deallocate real-type array
   deallocate(scaG)
   deallocate(scaG2)
   deallocate(flame_ratio)
   deallocate(deton_ratio)
   deallocate(flame_ratio_old)
   deallocate(deton_ratio_old)
   deallocate(burn_ratio)

   ! Deallocate extra array as well
   deallocate(flame_loc_ratio)
   deallocate(deton_loc_ratio)

   end subroutine destroyLevelset

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine output the global quantities
   ! derived from level-set
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine outputLevelset_log
   use definition
   implicit none
 
   ! Open, write and close !
   open(unit=401, file='./Outfile/Star_WENO_Ash_0.dat', action='write', position='append')
   WRITE (401, 701) global_time, mass_ash
   close(401)

   701 format (4ES18.8)

   end subroutine outputLevelset_log

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine opens the files and output the 
   ! profiles of level-set
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine outputLevelSet(n)
   use definition
   implicit none
 
   ! Input file number
   integer :: n

   ! Package for file name
   INTEGER :: fileno_len
   CHARACTER (len = 256) :: fileno

   ! Dummy variables
   INTEGER :: j, k

   ! Set up the files
   WRITE (fileno, *) n
   fileno = ADJUSTL (fileno)
   fileno_len = LEN_TRIM (fileno)

   ! Do the output
   OPEN(UNIT=601, FILE='./Outfile/Flame/Star_WENO_FlameRatioOld_'//fileno(1:fileno_len)//'.dat', action='write')
   write(601, *) n, global_time      
   write(601, *) length_step_r_2, length_step_z_2
   write(601, *) dx2, dt
   DO k = 1, length_step_z_2
       write(601, 701) (flame_ratio_old(j,k), j=1, length_step_r_2)
   enddo
   write(601, *)               
   CLOSE(601)

   OPEN(UNIT=601, FILE='./Outfile/Flame/Star_WENO_FlameRatio_'//fileno(1:fileno_len)//'.dat', action='write')                  
   write(601, *) n, global_time 
   write(601, *) length_step_r_2, length_step_z_2
   write(601, *) dx2, dt
   DO k = 1, length_step_z_2
       write(601, 701) (flame_ratio(j,k), j=1, length_step_r_2)
   enddo
   write(601, *)
   CLOSE(601)

   OPEN(UNIT=601, FILE='./Outfile/Flame/Star_WENO_scaG_'//fileno(1:fileno_len)//'.dat', action='write')        
   write(601, *) n, global_time 
   write(601, *) length_step_r_2, length_step_z_2
   write(601, *) dx2, dt
   DO k = 1, length_step_z_2
       write(601, 701) (scaG(j,k), j=1, length_step_r_2)
   enddo
   write(601, *)
   CLOSE(601)

   OPEN(UNIT=601, FILE='./Outfile/Flame/Star_WENO_DetonRatio_'//fileno(1:fileno_len)//'.dat', action='write')                  
   write(601, *) n, global_time 
   write(601, *) length_step_r_2, length_step_z_2
   write(601, *) dx2, dt
   DO k = 1, length_step_z_2
       write(601, 701) (deton_ratio(j,k), j=1, length_step_r_2)
   enddo
   write(601, *)
   CLOSE(601)

   OPEN(UNIT=601, FILE='./Outfile/Flame/Star_WENO_DetonRatioOld_'//fileno(1:fileno_len)//'.dat', action='write')                  
   write(601, *) n, global_time 
   write(601, *) length_step_r_2, length_step_z_2
   write(601, *) dx2, dt
   DO k = 1, length_step_z_2
       write(601, 701) (deton_ratio_old(j,k), j=1, length_step_r_2)
   enddo
   write(601, *)
   CLOSE(601)

   OPEN(UNIT=601, FILE='./Outfile/Flame/Star_WENO_scaG2_'//fileno(1:fileno_len)//'.dat', action='write')                  
   write(601, *) n, global_time 
   write(601, *) length_step_r_2, length_step_z_2
   write(601, *) dx2, dt
   DO k = 1, length_step_z_2
       write(601, 701) (scaG2(j,k), j=1, length_step_r_2)
   enddo
   write(601, *)
   CLOSE(601)

   701 FORMAT (2000ES13.5)

   end subroutine outputLevelSet

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine collects all the subroutines
   ! for updating the level-set
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_flame_radius()
   use definition
   implicit none

   ! Check timing with or without openmp
   INTEGER :: time_start, time_end
   INTEGER :: cr, cm
   REAL :: rate            

   ! Record the timing !
   !CALL system_clock(count_rate=cr)
   !CALL system_clock(count_max=cm)
   !rate = REAL(cr)
   !WRITE(*,*) "system_clock rate ",rate
   !CALL system_clock(time_start)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! Update the first level-set
   ! The deflagration is neglected as long as detonation has
   ! takes place long enough to surround the whole deflagration
   IF(found_deton_flag == 0 .or. (found_deton_flag == 1 .AND. ABS(global_time - found_deton_time) <= 25000.0D0)) THEN
      CALL identify_flamegrid()
      CALL compute_flameratio()
   ENDIF                

   ! Update the second level-set
   IF(found_deton_flag == 1) THEN
      CALL identify_detongrid()
      CALL compute_detonratio()
   ENDIF

   ! Backup the data
   flame_ratio_old = flame_ratio
   deton_ratio_old = deton_ratio
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Let the first level-set self-Propagate
   IF(found_deton_flag == 0 .OR. (found_deton_flag == 1 .AND. ABS(global_time - found_deton_time) <= 25000.0D0)) THEN
      CALL update_scaG()
   ENDIF

   ! Let the second level-set self-propagate
   IF(found_deton_flag == 1) THEN
      call update_scaG2()
   ENDIF

   ! Do the reinitilization and 
   ! comput the geometry for the 1st level-set
   IF(found_deton_flag == 0 .OR. (found_deton_flag == 1 .and. ABS(global_time - found_deton_time) <= 25000.0D0)) THEN
      CALL reinitialization2()
      CALL identify_flamegrid()
      CALL compute_flameratio()
   ENDIF

   ! Do the reinitilization and 
   ! compute the geometry for the 2nd level-set
   IF(found_deton_flag == 1) THEN
      CALL reinitialization3()
      CALL identify_detongrid()
      CALL compute_detonratio()
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Conmpute the total local fraction
   burn_ratio(:,:) = flame_ratio(:,:) + deton_ratio(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !CALL system_clock(time_end)   
   !WRITE(*,*) 'Level set = ', REAL(time_end - time_start) / rate

   end subroutine update_flame_radius

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine initialize the first level-set
   ! by assuming certain deflagration shape
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GetFlame(temp_flag)
   USE definition
   IMPLICIT NONE

   ! Input signal
   integer :: temp_flag 	! 1 = impose temp. change associated with the flame shape

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialization of the scalar G

   ! Dummy variables
   integer :: j , k

   ! Dist from the front
   real (DP) :: dist

   ! Initilization
   flame_ratio_old = 0.0D0
   flame_ratio = 0.0D0

   ! Assign initial flame
   do j = 1, length_step_r_2, 1   
      do k  = 1, length_step_z_2, 1 

	 ! Leon's Patch b1
	 !scaG(j,k) = -DSQRT((r2(j) - 135.44D0)**2 + (z2(k) - 135.44D0)**2) + 27.09D0

	 ! Leon's Patch c1
         !scaG(j,k) = 81.26D0 - DSQRT(r2(j)**2 + z2(k)**2) + & 
         !            41.31D0 * ABS(DSIN(ASIN(z2(k) / DSQRT(r2(j)**2 + z2(k)**2)) * 2.0D0))

	 ! Leon's Patch c3
         scaG(j,k) = 81.26D0 - DSQRT(r2(j)**2 + z2(k)**2) + & 
                     41.31D0 * ABS(DSIN(ASIN(z2(k) / DSQRT(r2(j)**2 + z2(k)**2)) * 6.0D0))

	 ! big c3 flame
	 !scaG(j,k) = 148.90D0 - DSQRT(r2(j)**2 + z2(k)**2) + &
         !            60.92D0 * ABS(DSIN(ASIN(z2(k) / DSQRT(r2(j)**2 + z2(k)**2)) * 6.0D0))

	 ! c3 flame
         !scaG(j,k) = 74.45D0 - DSQRT(r2(j)**2 + z2(k)**2) + & 
         !            30.46D0 * ABS(DSIN(ASIN(z2(k) / DSQRT(r2(j)**2 + z2(k)**2)) * 6.0D0))

	 ! small c3 flame
	 !scaG(j,k) = 37.23D0 - DSQRT(r2(j)**2 + z2(k)**2) + &
         !            15.23D0 * ABS(DSIN(ASIN(z2(k) / DSQRT(r2(j)**2 + z2(k)**2)) * 6.0D0))

	 ! b1 flame 23.93 -> 50 km, 28.72 -> 60 km
	 !scaG(j,k) = -DSQRT((r2(j) - 47.86D0)**2 + (z2(k) - 47.86D0)**2) + 7.0D0

	 ! b1 flame 23.93 -> 50 km, 28.72 -> 60 km
         !scaG(j,k) = -DSQRT((r2(j))**2 + (z2(k) - 101.52D0)**2) + 15.0D0

	 ! A bubble along z-axis in the helium sphere
	 !scaG(j,k) = -DSQRT((r2(j))**2 + (z2(k) - 1.0D3)**2) + 15.0D0

	 ! Two bubbles along two axis
	 

	 ! b5 flame Type A
         !scaG(j,k) = MAX(10.0D0 - DSQRT((DBLE(j-1)*dx2 - 20.0D0)**2 + (DBLE(k-1)*dx2 - 40.0D0)**2), &      
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 40.0D0)**2 + (DBLE(k-1)*dx2 - 20.0D0)**2), &
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 25.0D0)**2 + (DBLE(k-1)*dx2 - 75.0D0)**2), &
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 50.0D0)**2 + (DBLE(k-1)*dx2 - 50.0D0)**2), &
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 75.0D0)**2 + (DBLE(k-1)*dx2 - 25.0D0)**2))

	 ! b5 flame Type B
 	 !scaG(j,k) = MAX(10.0D0 - DSQRT((DBLE(j-1)*dx2 - 33.84D0)**2 + (DBLE(k-1)*dx2 - 67.68D0)**2), &      
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 67.68D0)**2 + (DBLE(k-1)*dx2 - 33.84D0)**2), &
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 51.80D0)**2 + (DBLE(k-1)*dx2 - 125.07D0)**2), &
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 95.72D0)**2 + (DBLE(k-1)*dx2 - 95.72D0)**2), &
         !               10.0D0 - DSQRT((DBLE(j-1)*dx2 - 125.07D0)**2 + (DBLE(k-1)*dx2 - 51.80D0)**2))

      enddo
   enddo 

   ! Copy the data to ghost cells
   CALL boundary1D_NMFULL(scaG,even)

   ! Readjust the initial data
   call reinitialization2()               !You need this for aspherical initial flame
   call identify_flamegrid()
   call compute_flameratio()
   
   ! Backup initial result
   flame_ratio_old = flame_ratio

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Setup initial temperature due to flame
   ! Modify the temp profile too for a more-consistent model
   !if(temp_flag == 1) then
   !   do j = 1, length_step_r_2, 1
   !      do k = 1, length_step_z_2, 1
   !         if(flame_ratio(j,k) > 0.0D0) temp2 (j,k) = 2.0D0 * flame_ratio(j,k) + temp_a * (1.0D0 - flame_ratio(j,k))
   !      enddo
   !   enddo
   !   temp2_old = temp2
   !   CALL boundary1D_NMFULL(temp2,even)
   !   call getepsilon
   !endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set the deton part
   deton_ratio(:,:) = 0.0D0
   deton_ratio_old(:,:) = 0.0D0
   !call update_scaG2corn()

   ! Sum all level-set ratio up
   burn_ratio(:,:) = flame_ratio(:,:) + deton_ratio(:,:) 

   end subroutine GetFlame

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine is very similar to GetFlame but for 
   ! the second level-set.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine GetDeton(temp_flag)
   use definition
   implicit none

   ! Signal for including temp change
   integer :: temp_flag         ! 1 = impose temp. change associated with the flame shape        

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
   ! Initialization of the scalar G

   ! Dummy variable
   integer :: j , k
  
   ! Distance from the front
   real (DP) :: dist

   ! Initilization
   deton_ratio_old = 0.0D0
   deton_ratio = 0.0D0

   ! This command is needed if we want detonation
   ! to be taken care in the beginning
   !found_deton_flag = .true.

   do j = 1, length_step_r_2, 1
      do k  = 1, length_step_z_2, 1

	 ! Spehrical 
	 !scaG2(j,k) = 50.0D0 - DSQRT(r2(j)**2 + z2(k)**2)

	 ! b1 structure assumed
         scaG2(j,k) = -DSQRT((r2(j) - 478.6D0)**2 + (z2(k) - 478.6D0)**2) + 10.0D0

	 ! A bubble along z-axis in the helium sphere
	 !scaG2(j,k) =-DSQRT((r2(j))**2 + (z2(k) - (radius + 30.0D0))**2) + 15.0D0

	 ! 2 bubbles
         !scaG2(j,k) = MAX(-DSQRT((r2(j))**2 + (z2(k) - (radius + 30.0D0))**2) + 15.0D0, &
	 ! 		  -DSQRT((r2(j))**2 + (z2(k) + (radius + 30.0D0))**2) + 15.0D0)

	 ! Abubble along the r-axis in the helium sphere
	 !scaG2(j,k) = -DSQRT((r2(j))**2 + (z2(k) - (radius + 50.0D0))**2) + 40.0D0 
	 !dist = SQRT(r2(j)**2 + z2(k)**2)
	 !if(dist < radius + 50.0D0) then
   	 !   scaG2(j,k) = -(radius + 30.0D0) + dist
	 !else
	 !   scaG2(j,k) = 20.0D0 - (dist - (radius + 50.0D0))
	 !endif

	 ! 2 bubbles along the r-axis and z-axis
	 !scaG2(j,k) = MAX(-DSQRT((r2(j))**2 + (z2(k) - (radius + 50.0D0))**2) + 40.0D0, &
	 !		  -DSQRT((r2(j) - (radius + 50.0D0))**2 + (z2(k))**2) + 40.0D0)

	 ! 3 bubbles 
         !scaG2(j,k) = MAX(-DSQRT((r2(j))**2 + (z2(k) - (radius + 50.0D0))**2) + 40.0D0, &
	 !		  -DSQRT((r2(j) - (radius + 50.0D0))**2 + (z2(k))**2) + 40.0D0, &
	 !		  -DSQRT((r2(j) - 0.70710678D0*(radius + 50.0D0))**2 + (z2(k) - 0.70710678D0*(radius + 50.0D0))**2) + 40.0D0)

      enddo               
   enddo              

   ! Copy the results to ghost cells
   CALL boundary1D_NMFULL(scaG2,even)

   ! Readjust the initial input
   call reinitialization3()               !You need this for aspherical initial flame
   call identify_detongrid()      
   call compute_detonratio()

   ! Backup the data
   deton_ratio_old = deton_ratio
   burn_ratio(:,:) = flame_ratio(:,:) + deton_ratio(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

   ! Setup initial temperature due to flame

   flame_rad = 0.0D0

   ! Set the flame part
   !flame_ratio(:,:) = 0.0D0
   !flame_ratio_old(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

   ! Sum the fractions from all level-set
   burn_ratio(:,:) = flame_ratio(:,:) + deton_ratio(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

   found_deton_time = 0.0D0
   found_deton_flag = 1
   output_profiletime = 0.5D4
   output_flametime = 0.5D4
   output_turbtime = 0.5D4
   output_Helmtime = 0.5D4   

   end subroutine GetDeton

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine classifies the type of grid
   ! (filled or unfilled?) / (partially or completely?)
   !
   ! The classification is defined as
   ! 
   ! 0 = completely unburnt
   ! 1 = completely burnt
   ! 2 = 
   !
   ! The data transfer includes:
   ! In: scaG_in (the level-set field)
   ! Out: 
   !  
   !
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine identify_flamegrid()
   use definition
   implicit none

   ! Dummy variables 
   integer :: j, k

   ! Initilization
   flamecorn_flag = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

   do j = 0, length_step_r_2, 1
      do k = 0, length_step_z_2, 1

	 ! The technique is as follows
	 ! If all corners are positive/negative, then it is either completely occpied/emptied
	 ! If any one or above corners have different signs, that means the front cut this grid
	 ! which means you need to calculate what fraction this front occupies the grid 

         if(scaG(j,k) > 0.0D0 .or. scaG(j+1,k) > 0.0D0 .or. scaG(j,k+1) > 0.0D0 .or. scaG(j+1,k+1) > 0.0D0) then
	    if(scaG(j,k) > 0.0D0 .and. scaG(j+1,k) > 0.0D0 .and. scaG(j,k+1) > 0.0D0 .and. scaG(j+1,k+1) > 0.0D0) then

	       ! Completely occupied
	       flamecorn_flag(j,k) = 1

	    elseif(scaG(j,k) * scaG(j+1,k) * scaG(j,k+1) * scaG(j+1,k+1) > 0.0D0) then

 	       if((scaG(j,k) * scaG(j,k+1) <= 0.0D0) .and. (scaG(j,k) * scaG(j+1,k) <= 0.0D0)) then
		  ! The ambguous case
	          flamecorn_flag(j,k) = 4
	       else
		  ! The parallelogram case
	          flamecorn_flag(j,k) = 3
	       endif

	    else
	
	       ! The triangle case
	       flamecorn_flag(j,k) = 2

	    endif

         elseif(scaG(j,k) <= 0.0D0 .and. scaG(j+1,k) <= 0.0D0 .and. &
                scaG(j,k+1) <= 0.0D0 .and. scaG(j+1,k+1) <= 0.0D0) then

	    ! Empty case
            flamecorn_flag(j,k) = 0

         endif

      enddo
   enddo

   !CALL BOUNDARY1D_INT(flamegrid_flag, even)

   100 FORMAT (6I5)
   101 FORMAT (6E13.6)

   end subroutine identify_flamegrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine classifies the type of grid 
   ! (filled or unfilled?) / (partially or completely?)
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine identify_detongrid()
   use definition
   implicit none

   ! dummy variables
   integer :: j, k  

   ! initialization
   detoncorn_flag = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k = 0, length_step_z_2, 1   
      do j = 0, length_step_r_2, 1

         if(scaG2(j,k) > 0.0D0 .or. scaG2(j+1,k) > 0.0D0 .or. scaG2(j,k+1) > 0.0D0 .or. scaG2(j+1,k+1) > 0.0D0) then

            if(scaG2(j,k) > 0.0D0 .and. scaG2(j+1,k) > 0.0D0 .and. scaG2(j,k+1) > 0.0D0 .and. scaG2(j+1,k+1) > 0.0D0) then
	      
	       ! Completely filled
               detoncorn_flag(j,k) = 1

            elseif(scaG2(j,k) * scaG2(j+1,k) * scaG2(j,k+1) * scaG2(j+1,k+1) > 0.0D0) then

               if((scaG2(j,k) * scaG2(j,k+1) <= 0.0D0) .and. (scaG2(j,k) * scaG2(j+1,k) <= 0.0D0)) then

		   ! The ambigous case
                  detoncorn_flag(j,k) = 4

               else

		  ! The parallelogram case
                  detoncorn_flag(j,k) = 3

               endif

            else

	       ! The triangle case
               detoncorn_flag(j,k) = 2

            endif

         elseif(scaG2(j,k) <= 0.0D0 .and. scaG2(j+1,k) <= 0.0D0 .and. &
                scaG2(j,k+1) <= 0.0D0 .and. scaG2(j+1,k+1) <= 0.0D0) then

	    ! The empty case
            detoncorn_flag(j,k) = 0

         endif
      enddo
   enddo

   !CALL BOUNDARY1D_INT(detongrid_flag, even)

   100 FORMAT (6I5)
   101 FORMAT (6E13.6)

   end subroutine identify_detongrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the fraction occupied by one of
   ! the phase (ash if in the context of supernova)
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_flameratio()
   use definition
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Local storage
   real (DP) :: loc_flame_ratio 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! I need to include the opposite case as well !

   ! Initialization
   !write(*,*) 'In compute norm snd flame ratio'
   flame_ratio = 0.0D0

   do k = 0, length_step_z_2, 1
      do j = 0, length_step_r_2, 1

         if(flamecorn_flag(j,k) /= 0) then

	    if(flamecorn_flag(j,k) == 1) then

	       ! Completely filled is one	       
	       loc_flame_ratio = 1.0D0

            elseif(flamecorn_flag(j,k) == 2) then

	       ! Incomplete filled needs to be studied case by case
	
	       if(scaG(j,k) >= 0.0D0 .and. scaG(j,k+1) < 0.0D0 .and. &
                  scaG(j+1,k) < 0.0D0 .and. scaG(j+1,k+1) < 0.0D0) then
                  loc_flame_ratio = ((scaG(j,k)/(scaG(j,k)-scaG(j+1,k))) * & 
                                     (scaG(j,k)/(scaG(j,k)-scaG(j,k+1)))) * 0.5D0
               elseif(scaG(j,k) < 0.0D0 .and. scaG(j,k+1) >= 0.0D0 .and. &
                  scaG(j+1,k) >= 0.0D0 .and. scaG(j+1,k+1) >= 0.0D0) then
                  loc_flame_ratio = 1.0D0 - ((scaG(j,k)/(scaG(j,k)-scaG(j+1,k))) * &
                                     (scaG(j,k)/(scaG(j,k)-scaG(j,k+1)))) * 0.5D0

               elseif(scaG(j,k) < 0.0D0 .and. scaG(j,k+1) >= 0.0D0 .and. &
                      scaG(j+1,k) < 0.0D0 .and. scaG(j+1,k+1) < 0.0D0) then
                  loc_flame_ratio = ((scaG(j,k+1)/(scaG(j,k+1)-scaG(j+1,k+1))) * &
                                     (scaG(j,k+1)/(scaG(j,k+1)-scaG(j,k)))) * 0.5D0
               elseif(scaG(j,k) >= 0.0D0 .and. scaG(j,k+1) < 0.0D0 .and. &
                      scaG(j+1,k) >= 0.0D0 .and. scaG(j+1,k+1) >= 0.0D0) then
                  loc_flame_ratio = 1.0D0 - ((scaG(j,k+1)/(scaG(j,k+1)-scaG(j+1,k+1))) * &
                                     (scaG(j,k+1)/(scaG(j,k+1)-scaG(j,k)))) * 0.5D0

               elseif(scaG(j,k) < 0.0D0 .and. scaG(j,k+1) < 0.0D0 .and. &
                      scaG(j+1,k) >= 0.0D0 .and. scaG(j+1,k+1) < 0.0D0) then
                  loc_flame_ratio = ((scaG(j+1,k)/(scaG(j+1,k)-scaG(j+1,k+1))) * &
                                     (scaG(j+1,k)/(scaG(j+1,k)-scaG(j,k)))) * 0.5D0    
               elseif(scaG(j,k) >= 0.0D0 .and. scaG(j,k+1) > 0.0D0 .and. &
                      scaG(j+1,k) < 0.0D0 .and. scaG(j+1,k+1) > 0.0D0) then
                  loc_flame_ratio = 1.0D0 - ((scaG(j+1,k)/(scaG(j+1,k)-scaG(j+1,k+1))) * &
                                     (scaG(j+1,k)/(scaG(j+1,k)-scaG(j,k)))) * 0.5D0

               elseif(scaG(j,k) < 0.0D0 .and. scaG(j,k+1) < 0.0D0 .and. &
                      scaG(j+1,k) < 0.0D0 .and. scaG(j+1,k+1) >= 0.0D0) then
                  loc_flame_ratio = ((scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j+1,k))) * &
                                     (scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j,k+1)))) * 0.5D0
               elseif(scaG(j,k) >= 0.0D0 .and. scaG(j,k+1) >= 0.0D0 .and. &
                      scaG(j+1,k) >= 0.0D0 .and. scaG(j+1,k+1) < 0.0D0) then
                  loc_flame_ratio = 1.0D0 - ((scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j+1,k))) * &
                                     (scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j,k+1)))) * 0.5D0
               endif

   	    elseif(flamecorn_flag(j,k) == 3) then

	       ! Parallelogram case

               if(scaG(j,k) >= 0.0D0 .and. scaG(j,k+1) < 0.0D0) then
                  loc_flame_ratio = ((scaG(j,k)/(scaG(j,k)-scaG(j,k+1))) + &
                                     (scaG(j+1,k)/(scaG(j+1,k)-scaG(j+1,k+1)))) * 0.5D0  
               elseif(scaG(j,k) < 0.0D0 .and. scaG(j,k+1) >= 0.0D0) then      
                  loc_flame_ratio = ((scaG(j,k+1)/(scaG(j,k+1)-scaG(j,k))) + &
                                     (scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j+1,k)))) * 0.5D0

               elseif(scaG(j,k) >= 0.0D0 .and. scaG(j+1,k) < 0.0D0) then
                  loc_flame_ratio = ((scaG(j,k)/(scaG(j,k)-scaG(j+1,k))) + &
                                     (scaG(j,k+1)/(scaG(j,k+1)-scaG(j+1,k+1)))) * 0.5D0
               elseif(scaG(j,k) < 0.0D0 .and. scaG(j+1,k) >= 0.0D0) then
                  loc_flame_ratio = ((scaG(j+1,k)/(scaG(j+1,k)-scaG(j,k))) + &  
                                     (scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j,k+1)))) * 0.5D0

               endif
	       
       	    elseif(flamecorn_flag(j,k) == 4) then

	       ! The ambigorous case averages all possibilities

	       if(scaG(j,k) > 0.0D0) then
                  loc_flame_ratio = 0.5D0 * &
                                     ((((scaG(j,k)/(scaG(j,k) - scaG(j,k+1))) * &
                                     (scaG(j,k)/(scaG(j,k) - scaG(j+1,k)))) * 0.5D0 + &
                                     ((scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j+1,k))) * &
                                     (scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j,k+1)))) * 0.5D0) + &
                                     (1.0D0 - ((scaG(j,k+1)/(scaG(j,k+1) - scaG(j,k))) * &
                                     (scaG(j,k+1)/(scaG(j,k+1) - scaG(j+1,k+1)))) * 0.5D0 - &  
                                     ((scaG(j+1,k)/(scaG(j+1,k)-scaG(j,k))) * &
                                     (scaG(j+1,k)/(scaG(j+1,k)-scaG(j+1,k+1)))) * 0.5D0))
               else
                  loc_flame_ratio = 0.5D0 * &
                                     ((1.0D0 - ((scaG(j,k)/(scaG(j,k) - scaG(j,k+1))) * &
                                     (scaG(j,k)/(scaG(j,k) - scaG(j+1,k)))) * 0.5D0 - &
                                     ((scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j+1,k))) * &
                                     (scaG(j+1,k+1)/(scaG(j+1,k+1)-scaG(j,k+1)))) * 0.5D0) + &
                                     (((scaG(j,k+1)/(scaG(j,k+1) - scaG(j,k))) * &  
                                     (scaG(j,k+1)/(scaG(j,k+1) - scaG(j+1,k+1)))) * 0.5D0 + & 
                                     ((scaG(j+1,k)/(scaG(j+1,k)-scaG(j,k))) * &             
                                     (scaG(j+1,k)/(scaG(j+1,k)-scaG(j+1,k+1)))) * 0.5D0))
               endif
            
   	    endif

         else

	    ! Completely empty is zero
            loc_flame_ratio = 0.0D0

         endif

	 ! Notice that the ratio is defined with grid corner
	 ! as center, so the contribution of the fraction is 
	 ! applied to all the overlapped grids
	 loc_flame_ratio = 0.25D0 * loc_flame_ratio

	 flame_ratio(j,k) = flame_ratio(j,k) + loc_flame_ratio
	 flame_ratio(j+1,k) = flame_ratio(j+1,k) + loc_flame_ratio
	 flame_ratio(j,k+1) = flame_ratio(j,k+1) + loc_flame_ratio
   	 flame_ratio(j+1,k+1) = flame_ratio(j+1,k+1) + loc_flame_ratio

      enddo
   enddo

   ! Copy the result to the current fraction
   flame_loc_ratio(:,:) = flame_ratio(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! According to the summed fraction, find the 
   ! grid type according to the geometry
   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1

	 ! Note: For this case, only three types
	 ! of grid are needed because we have 
	 ! already found the fraction	 

	 ! Get the effective flame ratio
         !flame_ratio(j,k) = MIN(MAX(flame_ratio(j,k), flame_ratio_old(j,k)), 1.0D0 - deton_ratio(j,k))
	 flame_ratio(j,k) = MAX(flame_ratio(j,k), flame_ratio_old(j,k))
	 !flame_ratio(j,k) = MIN(flame_ratio(j,k), 1.0D0 - deton_ratio(j,k))

	 if(flame_ratio(j,k) == 0.0D0) then

	    ! Completely empty
            flamegrid_flag(j,k) = 0

         elseif(flame_ratio(j,k) == 1.0D0) then

	    ! Completely filled
            flamegrid_flag(j,k) = 1

         else
	
	    ! Partially filled	
            flamegrid_flag(j,k) = 2

         endif

      enddo
   enddo

   ! Copy the results to ghost cells
   CALL BOUNDARY1D_NMFULL(flame_ratio,even)
   !CALL BOUNDARY1D_INT(flamegrid_flag,even)

   100 FORMAT (6E13.6)

   end subroutine compute_flameratio

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the fraction occupied by one of
   ! the phase (ash if in the context of supernova)
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_detonratio()
   use definition
   implicit none

   ! Dummy variables         
   integer :: j, k

   ! Local stroage
   real (DP) :: loc_deton_ratio

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
   ! I need to include the opposite case as well !
   !write(*,*) 'In compute norm snd flame ratio'
   deton_ratio(:,:) = 0.0D0
       
   do k = 0, length_step_z_2, 1
      do j = 0, length_step_r_2, 1

         if(detoncorn_flag(j,k) /= 0) then

            if(detoncorn_flag(j,k) == 1) then

	       ! Completely filled is one
               loc_deton_ratio = 1.0D0

            elseif(detoncorn_flag(j,k) == 2) then

	       ! For triangle case, need to study case by case

	       if(scaG2(j,k) >= 0.0D0 .and. scaG2(j,k+1) < 0.0D0 .and. &
                  scaG2(j+1,k) < 0.0D0 .and. scaG2(j+1,k+1) < 0.0D0) then
                  loc_deton_ratio = ((scaG2(j,k)/(scaG2(j,k)-scaG2(j+1,k))) * &
                                     (scaG2(j,k)/(scaG2(j,k)-scaG2(j,k+1)))) * 0.5D0
               elseif(scaG2(j,k) < 0.0D0 .and. scaG2(j,k+1) >= 0.0D0 .and. &     
                  scaG2(j+1,k) >= 0.0D0 .and. scaG2(j+1,k+1) >= 0.0D0) then
                  loc_deton_ratio = 1.0D0 - ((scaG2(j,k)/(scaG2(j,k)-scaG2(j+1,k))) * &
                                     (scaG2(j,k)/(scaG2(j,k)-scaG2(j,k+1)))) * 0.5D0
                                  
               elseif(scaG2(j,k) < 0.0D0 .and. scaG2(j,k+1) >= 0.0D0 .and. &
                      scaG2(j+1,k) < 0.0D0 .and. scaG2(j+1,k+1) < 0.0D0) then
                  loc_deton_ratio = ((scaG2(j,k+1)/(scaG2(j,k+1)-scaG2(j+1,k+1))) * &
                                     (scaG2(j,k+1)/(scaG2(j,k+1)-scaG2(j,k)))) * 0.5D0
               elseif(scaG2(j,k) >= 0.0D0 .and. scaG2(j,k+1) < 0.0D0 .and. &   
                      scaG2(j+1,k) >= 0.0D0 .and. scaG2(j+1,k+1) >= 0.0D0) then
                  loc_deton_ratio = 1.0D0 - ((scaG2(j,k+1)/(scaG2(j,k+1)-scaG2(j+1,k+1))) * &
                                     (scaG2(j,k+1)/(scaG2(j,k+1)-scaG2(j,k)))) * 0.5D0
                                  
               elseif(scaG2(j,k) < 0.0D0 .and. scaG2(j,k+1) < 0.0D0 .and. &
                      scaG2(j+1,k) >= 0.0D0 .and. scaG2(j+1,k+1) < 0.0D0) then
                  loc_deton_ratio = ((scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j+1,k+1))) * &
                                     (scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j,k)))) * 0.5D0
               elseif(scaG2(j,k) >= 0.0D0 .and. scaG2(j,k+1) > 0.0D0 .and. &
                      scaG2(j+1,k) < 0.0D0 .and. scaG2(j+1,k+1) > 0.0D0) then
                  loc_deton_ratio = 1.0D0 - ((scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j+1,k+1))) * &
                                     (scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j,k)))) * 0.5D0

               elseif(scaG2(j,k) < 0.0D0 .and. scaG2(j,k+1) < 0.0D0 .and. &
                      scaG2(j+1,k) < 0.0D0 .and. scaG2(j+1,k+1) >= 0.0D0) then
                  loc_deton_ratio = ((scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j+1,k))) * &
                                     (scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j,k+1)))) * 0.5D0
               elseif(scaG2(j,k) >= 0.0D0 .and. scaG2(j,k+1) >= 0.0D0 .and. &
                      scaG2(j+1,k) >= 0.0D0 .and. scaG2(j+1,k+1) < 0.0D0) then
                  loc_deton_ratio = 1.0D0 - ((scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j+1,k))) * &
                                     (scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j,k+1)))) * 0.5D0
               endif
                                               
            elseif(detoncorn_flag(j,k) == 3) then

	       ! The parallelogram case

	       if(scaG2(j,k) >= 0.0D0 .and. scaG2(j,k+1) < 0.0D0) then
                  loc_deton_ratio = ((scaG2(j,k)/(scaG2(j,k)-scaG2(j,k+1))) + &
                                     (scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j+1,k+1)))) * 0.5D0
               elseif(scaG2(j,k) < 0.0D0 .and. scaG2(j,k+1) >= 0.0D0) then
                  loc_deton_ratio = ((scaG2(j,k+1)/(scaG2(j,k+1)-scaG2(j,k))) + &
                                  (scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j+1,k)))) * 0.5D0
            
               elseif(scaG2(j,k) >= 0.0D0 .and. scaG2(j+1,k) < 0.0D0) then
                  loc_deton_ratio = ((scaG2(j,k)/(scaG2(j,k)-scaG2(j+1,k))) + &
                                     (scaG2(j,k+1)/(scaG2(j,k+1)-scaG2(j+1,k+1)))) * 0.5D0
               elseif(scaG2(j,k) < 0.0D0 .and. scaG2(j+1,k) >= 0.0D0) then  
                  loc_deton_ratio = ((scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j,k))) + &
                                     (scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j,k+1)))) * 0.5D0

               endif
               
            elseif(detoncorn_flag(j,k) == 4) then

	       ! The ambigous case, average of all possible configuration

	       if(scaG2(j,k) > 0.0D0) then
                  loc_deton_ratio = 0.5D0 * &
                                     ((((scaG2(j,k)/(scaG2(j,k) - scaG2(j,k+1))) * &
                                     (scaG2(j,k)/(scaG2(j,k) - scaG2(j+1,k)))) * 0.5D0 + &  
                                     ((scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j+1,k))) * &
                                     (scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j,k+1)))) * 0.5D0) + &
                                     (1.0D0 - ((scaG2(j,k+1)/(scaG2(j,k+1) - scaG2(j,k))) * &
                                     (scaG2(j,k+1)/(scaG2(j,k+1) - scaG2(j+1,k+1)))) * 0.5D0 - &
                                     ((scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j,k))) * &
                                     (scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j+1,k+1)))) * 0.5D0))
               else
                  loc_deton_ratio = 0.5D0 * &
                                     ((1.0D0 - ((scaG2(j,k)/(scaG2(j,k) - scaG2(j,k+1))) * &
                                     (scaG2(j,k)/(scaG2(j,k) - scaG2(j+1,k)))) * 0.5D0 - &
                                     ((scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j+1,k))) * &
                                     (scaG2(j+1,k+1)/(scaG2(j+1,k+1)-scaG2(j,k+1)))) * 0.5D0) + &
                                     (((scaG2(j,k+1)/(scaG2(j,k+1) - scaG2(j,k))) * &
                                     (scaG2(j,k+1)/(scaG2(j,k+1) - scaG2(j+1,k+1)))) * 0.5D0 + &
                                     ((scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j,k))) * &
                                     (scaG2(j+1,k)/(scaG2(j+1,k)-scaG2(j+1,k+1)))) * 0.5D0))
               endif
                                  
            endif

         else   

	    ! Completely empty mean zero
            loc_deton_ratio = 0.0D0

         endif

	! Notice that the ratio is defined with grid corner
        ! as center, so the contribution of the fraction is
        ! applied to all the overlapped grids		
	loc_deton_ratio = 0.25D0 * loc_deton_ratio
        deton_ratio(j,k) = deton_ratio(j,k) + loc_deton_ratio
        deton_ratio(j+1,k) = deton_ratio(j+1,k) + loc_deton_ratio
        deton_ratio(j,k+1) = deton_ratio(j,k+1) + loc_deton_ratio
        deton_ratio(j+1,k+1) = deton_ratio(j+1,k+1) + loc_deton_ratio

     enddo                          
   enddo

   ! Copy the result to the current fraction
   deton_loc_ratio(:,:) = deton_ratio(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
   ! Now classify the grid according to the summed fraction    
   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1

	 ! Find out the effective ratio
         deton_ratio(j,k) = MIN(MAX(deton_ratio(j,k), deton_ratio_old(j,k)), 1.0D0 - flame_ratio(j,k))
	 !deton_ratio(j,k) = MIN(deton_ratio(j,k), 1.0D0 - flame_ratio(j,k))

	 if(deton_ratio(j,k) == 0.0D0) then

	    ! Completely empty
	    detongrid_flag(j,k) = 0

	 elseif(deton_ratio(j,k) == 1.0D0) then

	    ! Completely filled
	    detongrid_flag(j,k) = 1

	 else

	    ! Partially filled
	    detongrid_flag(j,k) = 2

	 endif

      enddo
   enddo

   ! Copy the results to ghost cells
   CALL BOUNDARY1D_NMFULL(deton_ratio,even)
   !CALL BOUNDARY1D_INT(detongrid_flag, even)
             
   100 FORMAT (6E13.6)
                                  
   end subroutine compute_detonratio

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine do the propagation of level-set due 
   ! to its normal-flow (propagation of flame in the 
   ! supernovae context)
   ! Written by Leung SHing CHi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_scaG()
   use definition
   use flametable_module
   implicit none

   ! dummy variables
   integer :: j, k  

   ! local geometry information
   real (DP) :: norm_r, norm_z

   ! local sound speed
   real (DP) :: cs_grid

   ! Local propagation
   real (DP), dimension(-4:length_step_r_2+5,-4:length_step_z_2+5) :: scaG_flux1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! initilization
   cs_grid = 0.0D0
   scaG_flux1 = 0.0D0

   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1

	 if(flamegrid_flag(j,k) > 1) then
            call update_flame_velocity(j, k, cs_grid)   
            scaG_flux1(j,k) = cs_grid
         endif

      enddo
   enddo       

   ! Copy the results to ghost cells
   call boundary1D_NMfull(scaG_flux1, even)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! And then, I calculate the flux by fluid motion

   ! Now update the level set
   do j = 1, length_step_r_2, 1
      do k = 1, length_step_z_2, 1
         scaG(j,k) = scaG(j,k) + dt * scaG_flux1(j,k)
      enddo
   enddo

   ! Fill the level set ghost grid
   CALL BOUNDARY1D_NMFULL(scaG, even)

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine do the propagation of level-set due
   ! to its normal-flow (propagation of flame in the
   ! supernovae context)
   ! Written by Leung SHing CHi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_scaG2()
   use definition
   use flametable_module
   implicit none
      
   ! dummy variables
   integer :: j, k

   ! Geometry information
   real (DP) :: norm_r, norm_z

   ! local sound speed
   real (DP) :: cs_grid

   ! local propagation
   real (DP), dimension(-4:length_step_r_2+5,-4:length_step_z_2+5) :: scaG_flux1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
   ! initialization
   cs_grid = 0.0D0
   scaG_flux1 = 0.0D0

   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1

	 ! Only the front propagate
         if(detongrid_flag(j,k) > 1) then

            if(rho2(j,k) > 3.2D-11) then        
                
               ! Pathological substained detonation
               call readtable_detonvel(rho2(j,k), cs_grid)

            elseif(rho2(j,k) > 1.6D-12 .and. rho2(j,k) < 3.2D-11) then  ! between 1x10^6 and 2x10^7

               ! Chapman-Jouguet detonation
               call findsoundspeed(j, k, cs_grid)

            else

               ! No detonation otherwise
               cs_grid = 0.0D0

            endif

            scaG_flux1(j,k) = cs_grid       

         endif

      enddo
   enddo
       
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! First, I calculate the flux by the flame
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! And then, I calculate the flux by fluid motion

   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1
         scaG2(j,k) = scaG2(j,k) + dt * scaG_flux1(j,k)
      enddo
   enddo
         
   CALL BOUNDARY1D_NMFULL(scaG2, even)

   end subroutine update_scaG2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine do the reinitialization in order to 
   ! maintain the distance property of the 2nd level set
   ! Based on the level set paper Reinecke (1999a)
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!
 
   subroutine reinitialization2
   use definition
   implicit none

   ! Dummy variables
   integer :: j, k

   ! The number of point intersected by the grid and the front
   integer :: flame_count

   ! The scaled distance
   real (DP) :: H_d
 
   ! The signed distance
   real (DP) :: sgn_scaG

   ! Array storing the minimal distance of the grid to the nearest front
   real (DP), dimension(-4:length_step_r_2+5, -4:length_step_z_2+5) :: flame_distance

   ! Array storing all intersection points
   real (DP), dimension(length_step_r_2 * length_step_z_2, 2):: flame_position

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First find all intersections
   call locate_flame_grid(1, flame_count, flame_position)

   ! Then find the minimal distance
   call locate_min_flame_distance(flame_count, flame_position, flame_distance)

   !$OMP PARALLEL DO PRIVATE(j, k, H_d, sgn_scaG)
   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1

	 ! Correct the level set
         !if(scaG(j,k) > 10.0D0 * dx2 .or. scaG(j,k) < -10.0D0 * dx2) cycle

         H_d = (1.0D0 - TANH(3.0D0 * (flame_distance(j,k) - 3.0D0 * dx2) / dx2)) / & 
               (1.0D0 - TANH(-9.0D0))
         if(scaG(j,k) >= 0.0D0) then
            sgn_scaG = 1.0D0
         else
            sgn_scaG = -1.0D0
         endif
         scaG(j,k) = H_d * scaG(j,k) + (1.0D0 - H_d) * sgn_scaG * flame_distance(j,k)

      enddo
   enddo
   !$OMP END PARALLEL DO
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Copy the results to ghost cells
   CALL BOUNDARY1D_NMFULL(scaG, even)  

   100 FORMAT (6E13.6)

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine do the reinitialization in order to 
   ! maintain the distance property of the 2nd level set
   ! Based on the level set paper Reinecke (1999a)
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This subroutine set the reinitialization, using grid as counter
   
   subroutine reinitialization3
   use definition
   implicit none

   ! Dummy variables                 
   integer :: j, k          

   ! The number of point intersected by the grid and the front
   integer :: flame_count

   ! The scaled distance
   real (DP) :: H_d

   ! The signed distance
   real (DP) :: sgn_scaG            

   ! Array storing the minimal distance of the grid to the nearest front
   real (DP), dimension(-4:length_step_r_2+5, -4:length_step_z_2+5) :: flame_distance

   ! Array storing all intersection points
   real (DP), dimension(length_step_r_2 * length_step_z_2, 2):: flame_position

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
   ! First find all intersections         
   call locate_flame_grid(2, flame_count, flame_position)

   ! Then find the minimal distance
   call locate_min_flame_distance(flame_count, flame_position, flame_distance)

   !$OMP PARALLEL DO PRIVATE(j,k,H_d,sgn_scaG)   
   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1  
        
	 ! Correct the level set
         !if(scaG(j,k) > 10.0D0 * dx2 .or. scaG(j,k) < -10.0D0 * dx2) cycle

         H_d = (1.0D0 - TANH(3.0D0 * (flame_distance(j,k) - 3.0D0 * dx2) / dx2)) / &
               (1.0D0 - TANH(-9.0D0))
         if(scaG2(j,k) >= 0.0D0) then
            sgn_scaG = 1.0D0
         else
            sgn_scaG = -1.0D0
         endif
         scaG2(j,k) = H_d * scaG2(j,k) + (1.0D0 - H_d) * sgn_scaG * flame_distance(j,k)

      enddo
   enddo
   !$OMP END PARALLEL DO
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Copy the results to the ghost cells
   CALL BOUNDARY1D_NMFULL(scaG2, even)

   100 FORMAT (6E13.6)

   end subroutine reinitialization3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine find all points where the 
   ! deflagration/detonation front cuts the grid
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutines find all the intersections points of the 
   ! front and the grid. The subroutine takes input of the mode
   ! Mode 1 = first level set
   ! Mode 2 = second level set
   ! and gives the number of intersection point and their
   ! corresponding positions. 
   !
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   subroutine locate_flame_grid(mode, flame_count, flame_position)
   use definition
   implicit none

   !dummy variables
   integer :: j, k
 
   ! Output intersection point number
   integer :: flame_count

   ! Input mode
   integer :: mode

   ! distance of grid to the intersection
   real (DP) :: flame_dx

   ! Output array of intersection point
   real (DP), dimension(length_step_r_2 * length_step_z_2, 2):: flame_position

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialization
   flame_count = 0 
   flame_rad = 0.0D0
 
   ! Debug !
   !write(*,*) 'In locate flame grid'

   if(mode == 1) then

      ! For the first level set
      do k = 1, length_step_z_2, 1
         do j = 1, length_step_r_2, 1
 
            if(scaG(j,k) * scaG(j+1,k) < 0.0D0) then
              flame_count = flame_count + 1
	      flame_dx = ABS(scaG(j,k)) / ABS(scaG(j+1,k) - scaG(j,k))
              flame_position(flame_count, 1) = r2(j) + flame_dx * dx2
	      flame_position(flame_count, 2) = z2(k)
              if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1)
            endif

            if(scaG(j,k) * scaG(j,k+1) < 0.0D0) then             
              flame_count = flame_count + 1
	      flame_dx = ABS(scaG(j,k)) / ABS(scaG(j,k+1) - scaG(j,k))
              flame_position(flame_count, 1) = r2(j)
              flame_position(flame_count, 2) = z2(k) + flame_dx * dx2
              if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1) 
            endif

         enddo
      enddo

   elseif(mode == 2) then

      ! For the second level set
      do k = 1,length_step_z_2, 1
	 do j = 1, length_step_r_2, 1

   	   if(scaG2(j,k) * scaG2(j+1,k) < 0.0D0) then
              flame_count = flame_count + 1
              flame_dx = ABS(scaG2(j,k)) / ABS(scaG2(j+1,k) - scaG2(j,k))
              flame_position(flame_count, 1) = r2(j) + flame_dx * dx2
              flame_position(flame_count, 2) = z2(k)
              if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1)
           endif

           if(scaG2(j,k) * scaG2(j,k+1) < 0.0D0) then
              flame_count = flame_count + 1
              flame_dx = ABS(scaG2(j,k)) / ABS(scaG2(j,k+1) - scaG2(j,k))
              flame_position(flame_count, 1) = r2(j) 
              flame_position(flame_count, 2) = z2(k) + flame_dx * dx2
              if(flame_position(flame_count, 1) > flame_rad) flame_rad = flame_position(flame_count, 1)
            endif

         enddo
      enddo

   endif

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the minimum distance from that
   ! particular grid point to the flame surface.
   ! The subroutines take input of the intersection
   ! point number, their positions. Then it gives
   ! the minimum distances from all grid point to 
   ! the surface
   !
   ! Written by Leung Shing Chi
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine locate_min_flame_distance(flame_count, flame_position, flame_distance)
   use definition
   implicit none

   ! dummy variables
   integer :: j, k, k2

   ! Input intersection point number
   integer :: flame_count

   ! The local minimal distance
   real (selected_real_kind(15,307)):: distance, last_distance

   ! Output array for the minimal distance
   real (selected_real_kind(15,307)), dimension(-4:length_step_r_2+5, -4:length_step_z_2+5):: flame_distance

   ! Input array of intersection point positions
   real (selected_real_kind(15,307)), dimension(length_step_r_2 * length_step_z_2, 2) :: flame_position

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !$OMP PARALLEL DO SHARED(flame_distance) PRIVATE(j,k,k2,last_distance,distance)
   do j = 1, length_step_r_2, 1
      do k = 1, length_step_z_2, 1

         !if(scaG(j,k) > 10.0D0 * dx2 .or. scaG(j,k) < -10.0D0 * dx2) cycle

         last_distance = 10.0D0 * DBLE(length_step_z_2) * dx2

	 ! Search for the minimal distance
         do k2 = 1, flame_count, 1

   	    distance = DSQRT((flame_position(k2,1) - r2(j)) ** 2 + & 
		   	     (flame_position(k2,2) - z2(k)) ** 2) 

  	    if((flame_count > 1 .and. distance <= last_distance) .or. flame_count == 1) then
	       flame_distance(j,k) = distance 
   	       last_distance = distance
            endif

	    if(k2 == 1) last_distance = distance
  
         enddo

      enddo
   enddo
   !$OMP END PARALLEL DO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Debug !
   !WRITE(*,*) flame_distance(1,1), flame_count

   100 FORMAT (6E13.6)

   end subroutine

end module levelset_module
