!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is the tracer particle module that records the		
! thermodynamics history of each particule trajectory		
! The information will be used to input as data post processing 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module PPT_module
use definition, only : DP
implicit none

! Number of particles !
integer, parameter :: nop = 40000 !25600

! Grid variables !
integer :: grid_p(1:nop,2)
real (DP), allocatable, dimension(:,:) :: dgrid_p

! Hydro dynamic variables for particles
real (DP), allocatable, dimension(:,:) :: posold_p
real (DP), allocatable, dimension(:,:) :: pos_p
real (DP), allocatable, dimension(:,:) :: pos2_p
real (DP), allocatable, dimension(:,:) :: pos3_p
real (DP), allocatable, dimension(:,:) :: vel_p
real (DP), allocatable, dimension(:,:) :: vel3_p
real (DP), allocatable, dimension(:) :: temp_p
real (DP), allocatable, dimension(:) :: rho_p
real (DP), allocatable, dimension(:) :: phi_p

! Whether particle leave box
logical, allocatable, dimension(:) :: leavebox_p

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine build arrays for PPT scheme
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine buildPPT
   implicit none

   ! Allocate
   allocate(dgrid_p(1:nop,2))
   allocate(pos_p(1:nop, 2))
   allocate(pos2_p(1:nop, 2))
   allocate(pos3_p(1:nop, 2))
   allocate(posold_p(1:nop, 2))
   allocate(vel_p(1:nop, 2))        
   allocate(vel3_p(1:nop, 2)) 
   allocate(temp_p(1:nop))  
   allocate(rho_p(1:nop))
   allocate(leavebox_p(1:nop))
   allocate(phi_p(1:nop))

   end subroutine buildPPT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine destroy arrays for PPT scheme
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine destroyPPT
   implicit none

   ! Deallocate !
   deallocate(dgrid_p)
   deallocate(pos_p)
   deallocate(pos2_p)
   deallocate(pos3_p)
   deallocate(posold_p)
   deallocate(vel_p)        
   deallocate(vel3_p) 
   deallocate(temp_p)  
   deallocate(rho_p)
   deallocate(leavebox_p)
   deallocate(phi_p)

   end subroutine destroyPPT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine output ppt profile for post-processing 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   subroutine outputPPT(n)
   use definition, only : global_time
   implicit none
  
   ! Integer !
   integer :: i, n
   integer :: sqrt_nop

   ! For profile output !
   character(len=256) :: fileno
   integer :: fileno_len

   write(fileno,*) n
   fileno = ADJUSTL(fileno)
   fileno_len = LEN_TRIM(fileno)

   ! Square of number of particle 
   sqrt_nop = DSQRT(DBLE(nop))

   ! Open file
   open(unit=401,file='./Tracer/Star_WENO_PPT_'//fileno(1:fileno_len)//'.dat',action='write',position='append')

   !!!!!!!!!!!!!!!!!!!!!!!!!
   !do i = 1, sqrt_nop, 1
   !   write(401,100) 
   !enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!

   ! Write it out and close !
   write(401,*) n, global_time
   write(401,*) nop
   do i = 1, nop, 1
      write(401,100) i, rho_p(i), temp_p(i), pos_p(i,1), pos_p(i,2), vel_p(i,1), vel_p(i,2), leavebox_p(i), phi_p(i)
   enddo 
   write(401,*)
   close(401)

   100 format (I6, 6ES16.8, L3, ES16.8)

   end subroutine outputPPT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine initialize all tracer particle variables  
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getppt
   implicit none

   write(*,*) 'Now initialize tracer particles'
   call getpos_p
   call getgrid_p
   call getnse_p
   call getdata_p

   end subroutine getppt

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine initialize the status of tracer particles 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getnse_p
   implicit none

   ! Integer !
   integer :: i

   do i = 1, nop, 1
      leavebox_p(i) = .false.
   enddo

   end subroutine getnse_p

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position before RK time step
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine evolve_p_1st
   use definition
   implicit none 

   ! Integer !
   INTEGER :: i     
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In Evolve_p_1st'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   posold_p = pos_p
   call getgrid_p
   call getvel_p
   !!!!!!!!!!!!!!!!!!!!
   !call update_p_half
   !!!!!!!!!!!!!!!!!!!!

   ! Do update !
   DO i = 1, nop, 1
      pos_p(i,1) = posold_p(i,1) + vel_p(i,1) * 0.391752226571890D0 * dt
      pos_p(i,2) = posold_p(i,2) + vel_p(i,2) * 0.391752226571890D0 * dt
   ENDDO

   ! Special care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      IF(hemisphere_flag == 0) THEN

         DO i = 1, nop, 1
   	    IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
            IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
         ENDDO

      ELSEIF(hemisphere_flag == 1) THEN

	 DO i = 1, nop, 1
            IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
         ENDDO

      ELSE
	
	 STOP 'In Update_p_half: Check hemisphere_flag'

      ENDIF
   ENDIF
 
   100 format(100ES15.7)

   end subroutine evolve_p_1st

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position before RK time step
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine evolve_p_2nd
   use definition
   implicit none     

   ! Integer !
   INTEGER :: i  
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In Evolve_p_1st'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   call getgrid_p
   call getvel_p
   !!!!!!!!!!!!!!!!!!!!
   !call update_p_half
   !!!!!!!!!!!!!!!!!!!!

   ! Do update !
   DO i = 1, nop, 1
      pos_p(i,1) = 0.444370493651235D0*posold_p(i,1)+0.555629506348765D0*pos_p(i,1)+0.368410593050371D0*dt*vel_p(i,1)
      pos_p(i,2) = 0.444370493651235D0*posold_p(i,2)+0.555629506348765D0*pos_p(i,2)+0.368410593050371D0*dt*vel_p(i,2)
      pos2_p(i,1) = pos_p(i,1)
      pos2_p(i,2) = pos_p(i,2)
   ENDDO

   ! Special care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      IF(hemisphere_flag == 0) THEN

         DO i = 1, nop, 1
   	    IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
            IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
         ENDDO

      ELSEIF(hemisphere_flag == 1) THEN

	 DO i = 1, nop, 1
            IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
         ENDDO

      ELSE
	
	 STOP 'In Update_p_half: Check hemisphere_flag'

      ENDIF
   ENDIF
 
   100 format(100ES15.7)

   end subroutine evolve_p_2nd

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position before RK time step
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine evolve_p_3rd
   use definition
   implicit none      

   ! Integer !
   INTEGER :: i  

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In Evolve_p_1st'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   call getgrid_p
   call getvel_p
   !!!!!!!!!!!!!!!!!!!!
   !call update_p_half
   !!!!!!!!!!!!!!!!!!!!

   ! Do update !
   DO i = 1, nop, 1
      pos_p(i,1) = 0.620101851488403D0*posold_p(i,1)+0.379898148511597D0*pos_p(i,1)+0.251891774271694D0*dt*vel_p(i,1)
      pos_p(i,2) = 0.620101851488403D0*posold_p(i,2)+0.379898148511597D0*pos_p(i,2)+0.251891774271694D0*dt*vel_p(i,2)
      pos3_p(i,1) = pos_p(i,1)
      pos3_p(i,2) = pos_p(i,2)
      vel3_p(i,1) = vel_p(i,1)
      vel3_p(i,2) = vel_p(i,2)
   ENDDO

   ! Special care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      IF(hemisphere_flag == 0) THEN

         DO i = 1, nop, 1
   	    IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
            IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
         ENDDO

      ELSEIF(hemisphere_flag == 1) THEN

	 DO i = 1, nop, 1
            IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
         ENDDO

      ELSE
	
	 STOP 'In Update_p_half: Check hemisphere_flag'

      ENDIF
   ENDIF
 
   100 format(100ES15.7)

   end subroutine evolve_p_3rd

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position before RK time step
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine evolve_p_4th
   use definition
   implicit none     

   ! Integer !
   INTEGER :: i  
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In Evolve_p_1st'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   call getgrid_p
   call getvel_p
   !!!!!!!!!!!!!!!!!!!!
   !call update_p_half
   !!!!!!!!!!!!!!!!!!!!

   ! Do update !
   DO i = 1, nop, 1
      pos_p(i,1) = 0.178079954393132D0*posold_p(i,1)+0.821920045606868D0*pos_p(i,1)+0.544974750228521D0*dt*vel_p(i,1)
      pos_p(i,2) = 0.178079954393132D0*posold_p(i,2)+0.821920045606868D0*pos_p(i,2)+0.544974750228521D0*dt*vel_p(i,2)
   ENDDO

   ! Special care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      IF(hemisphere_flag == 0) THEN

         DO i = 1, nop, 1
   	    IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
            IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
         ENDDO

      ELSEIF(hemisphere_flag == 1) THEN

	 DO i = 1, nop, 1
            IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
         ENDDO

      ELSE
	
	 STOP 'In Update_p_half: Check hemisphere_flag'

      ENDIF
   ENDIF
 
   100 format(100ES15.7)

   end subroutine evolve_p_4th

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position before RK time step
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine evolve_p_5th
   use definition
   implicit none      

   ! Integer !
   INTEGER :: i 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In Evolve_p_1st'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   call getgrid_p
   call getvel_p
   !!!!!!!!!!!!!!!!!!!!
   !call update_p_half
   !!!!!!!!!!!!!!!!!!!!

   ! Do update !
   DO i = 1, nop, 1
      pos_p(i,1) = 0.517231671970585D0*pos2_p(i,1)+0.096059710526147D0*pos3_p(i,1)+0.386708617503269D0*pos_p(i,1) &
		+ 0.063692468666290D0*dt*vel3_p(i,1) + 0.226007483236906D0*dt*vel_p(i,1)
      pos_p(i,2) = 0.517231671970585D0*pos2_p(i,2)+0.096059710526147D0*pos3_p(i,2)+0.386708617503269D0*pos_p(i,2) &
		+ 0.063692468666290D0*dt*vel3_p(i,2) + 0.226007483236906D0*dt*vel_p(i,2)
   ENDDO

   ! Special care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      IF(hemisphere_flag == 0) THEN

         DO i = 1, nop, 1
   	    IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
            IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
         ENDDO

      ELSEIF(hemisphere_flag == 1) THEN

	 DO i = 1, nop, 1
            IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
         ENDDO

      ELSE
	
	 STOP 'In Update_p_half: Check hemisphere_flag'

      ENDIF
   ENDIF
 
   100 format(100ES15.7)

   end subroutine evolve_p_5th

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position after RK time step
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine evolve_p_final
   use definition
   implicit none

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In Evolve_p_2nd'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!
   !call update_p_half
   !!!!!!!!!!!!!!!!!!!!!
   call getgrid_p
   call getdata_p
 
   100 format(100ES15.7)

   end subroutine evolve_p_final

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine update the PPT position using velocity 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_p_half
   USE definition
   IMPLICIT NONE

   ! Variables and integers !
   INTEGER :: i
   REAL (DP) :: xmax 
   REAL (DP) :: half_pi 
   REAL (DP) :: half_dt 

   ! half dt and pi !
   half_dt = 0.5D0 * dt
   half_pi = 0.5D0 * pi

   ! assign !
   xmax = (DBLE(length_step_r_2) - 0.5D0) * dx2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !IF(debug_flag == 1) write(*,*) 'In Update p half'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO i = 1, nop, 1
      pos_p(i,1) = pos_p(i,1) + vel_p(i,1) * half_dt
      pos_p(i,2) = pos_p(i,2) + vel_p(i,2) * half_dt
   ENDDO

   ! Special care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      IF(hemisphere_flag == 0) THEN

         DO i = 1, nop, 1
   	    IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
            IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
         ENDDO

      ELSEIF(hemisphere_flag == 1) THEN

	 DO i = 1, nop, 1
            IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)
         ENDDO

      ELSE
	
	 STOP 'In Update_p_half: Check hemisphere_flag'

      ENDIF
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ELSEIF(coordinate_flag == 2) THEN
   !   DO i = 1, nop, 1
   !      IF(pos_p(i,1) < 0.0D0) pos_p(i,1) = -pos_p(i,1)   
   !      IF(pos_p(i,2) < 0.0D0) pos_p(i,2) = -pos_p(i,2)
   !      IF(pos_p(i,2) > half_pi) pos_p(i,2) = half_pi - (pos_p(i,2) - half_pi)
   !   ENDDO
   !ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine update_p_half

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine get the initial position of PPT particles 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE getpos_p
   USE definition
   IMPLICIT NONE

   ! Integer !
   INTEGER :: i, j
   INTEGER :: start_j, grid_j
   INTEGER :: sqrt_nop

   ! Related to enclosed mass shell !
   REAL (DP) :: pos_r
   REAL (DP) :: target_mass
   REAL (DP) :: mass_shell

   ! Assign 
   mass_shell = mass2 / DSQRT(DBLE(nop))
   sqrt_nop = SQRT(DBLE(nop))

   ! Input the position 
   start_j = 1

   ! Speical assign for different initial profile !
   !IF(read_oldfile_flag == 0) THEN

      do i = 1, sqrt_nop, 1
         target_mass = (DBLE(i) - 0.5D0) * mass_shell
         do j = 1, length_step_r_2, 1
            if(m_r2(j) >= target_mass) then
                !!!!!!!!!!!!!
	        !grid_j = j
                !!!!!!!!!!!!!
		! Use linear interpolation !
		pos_r = r2(j-1) + (r2(j) - r2(j-1))* & 
		(target_mass - m_r2(j-1))/(m_r2(j) - m_r2(j-1))
                exit
            endif
         enddo

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !pos_r = (DBLE(grid_j) - 0.5D0) * dx
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 IF(hemisphere_flag == 1) THEN
            do j = sqrt_nop*(i-1)+1, sqrt_nop*i, 1 
       	       pos_p(j,1) = pos_r * DSIN((DBLE(MOD(j,sqrt_nop)) + 0.5D0) * pi / DBLE(sqrt_nop))
	       pos_p(j,2) = pos_r * DCOS((DBLE(MOD(j,sqrt_nop)) + 0.5D0) * pi / DBLE(sqrt_nop))
            enddo
	 ELSEIF(hemisphere_flag == 0) THEN
	    do j = sqrt_nop*(i-1)+1, sqrt_nop*i, 1 
               pos_p(j,1) = pos_r * DSIN((DBLE(MOD(j,sqrt_nop)) + 0.5D0) * 0.5D0 * pi / DBLE(sqrt_nop))
               pos_p(j,2) = pos_r * DCOS((DBLE(MOD(j,sqrt_nop)) + 0.5D0) * 0.5D0 * pi / DBLE(sqrt_nop))
            enddo
	 ENDIF
      enddo
   !END IF
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ELSEIF(read_oldfile_flag == 2) THEN
   !   IF(coordinate_flag /= 2) STOP 'In GetPos_P: CAnnot read stellar profile without spherical coordinate'
   !   IF(masscut_flag == 0) STOP 'In GetPos_P: Cannot import stellar model without masscut'
   !   IF(fornax_flag == 0) STOP 'In GetPos_P:CAnnot import stellat model without fornax'
   !   do i = 1, sqrt_nop, 1
   !      target_mass = (DBLE(i) - 0.5D0) * mass_shell
   !      do j = start_j, length_step_r_2 * 2, 1
   !         if(m_r(j) >= target_mass .and. m_r(j-1) <= target_mass) then
   !            grid_j = j
   !            exit
   !         endif
   !      enddo
   !      !pos_r = (DBLE(grid_j) - 0.5D0) * dx
   !      pos_r = r2(grid_j) + (target_mass - m_r(j-1)) / (m_r(j) - m_r(j-1)) * (r2(grid_j) - r2(grid_j-1))
   !      do j = sqrt_nop*(i-1)+1, sqrt_nop*i, 1 
   !         pos_p(j,1) = pos_r
   !         pos_p(j,2) = (DBLE(MOD(j,sqrt_nop)) + 0.5D0) * 0.5D0 * pi / DBLE(sqrt_nop)
   !      enddo
   !   ENDDO
   !ENDIF  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   END SUBROUTINE getpos_p

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine to get the PPT particle grid point 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE getgrid_p
   USE definition
   IMPLICIT NONE

   ! Integer !
   integer :: i

   ! Real variables !
   REAL (DP) :: half_dx
   REAL (DP) :: ratio, grid_p_raw

   ! For estimating homologous expansion
   REAL (DP) :: rad

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug !
   !if(debug_flag == 1) write(*,*) 'In Getgrid_p'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign  !
   half_dx = dx2 * 0.5D0

   ! Speical care for different coordinate system 
   IF(coordinate_flag == 1) THEN

      ! In Cylindrical coordinate
      IF(hemisphere_flag == 1) THEN

         DO i = 1, nop, 1
            
            ! Assume the hydro quantity is defined on the grid center
            grid_p(i,1) = INT((pos_p(i,1) + half_dx) / dx2) + 1
            grid_p(i,2) = INT((pos_p(i,2) + half_dx) / dx2) + 1
            dgrid_p(i,1) = (pos_p(i,1) - (DBLE(grid_p(i,1)) - 1.5D0) * dx2) / dx2
            dgrid_p(i,2) = (pos_p(i,2) - (DBLE(grid_p(i,2)) - 1.5D0) * dx2) / dx2	    

            ! Assume the hydro quantity satisfies the whole gird
            !grid_p(i,1) = INT(pos_p(i,1) / dx2) + 1
            !grid_p(i,2) = INT(pos_p(i,2) / dx2) + 1 + length_step_z_2 / 2
            !dgrid_p(i,1) = 0.0D0
            !dgrid_p(i,2) = 0.0D0

         ENDDO
 
      ELSEIF(hemisphere_flag == 0) THEN

	 DO i = 1, nop, 1
 
            ! Assume the hydro quantity is defined on the grid center
            grid_p(i,1) = INT((pos_p(i,1) + half_dx) / dx2) + 1              
            grid_p(i,2) = INT((pos_p(i,2) + half_dx) / dx2) + 1
            dgrid_p(i,1) = (pos_p(i,1) - (DBLE(grid_p(i,1)) - 1.5D0) * dx2) / dx2
            dgrid_p(i,2) = (pos_p(i,2) - (DBLE(grid_p(i,2)) - 1.5D0) * dx2) / dx2
            
            ! Assume the hydro quantity satisfies the whole gird
            !grid_p(i,1) = INT(pos_p(i,1) / dx2) + 1
            !grid_p(i,2) = INT(pos_p(i,2) / dx2) + 1 
            !dgrid_p(i,1) = 0.0D0
            !dgrid_p(i,2) = 0.0D0

         ENDDO

      ENDIF

   ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ELSEIF(coordinate_flag == 2) THEN
   !   IF(masscut_flag == 0) THEN
   !	 STOP 'In GetGrid_p, not implemented at this moment'
   !   ELSEIF(masscut_flag == 1) THEN
   !	 IF(fornax_flag == 1) THEN
   !	    DO i = 1, nop, 1
   !	       ratio = (pos_p(i,1) - rad_cut) / A_fornax
   !            IF(ratio > 0.0D0) THEN
   !               grid_p_raw = A_fornax * LOG(ratio + SQRT(ratio**2 + 1.0D0))
   !               grid_p(i,1) = INT(grid_p_raw / dx) + 1
   !               grid_p(i,2) = INT(pos_p(i,2) / dtheta) + 1
   !            ELSE
   !               grid_p(i,1) = 0
   !               grid_p(i,2) = 0
   !            ENDIF
   !	    ENDDO
   !	 ENDIF
   !   ENDIF
   !ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'Done getgrid_p!'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   END SUBROUTINE getgrid_p

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine to get the PPT particle hydro status 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE getdata_p
   USE definition
   USE helmeos_module
   IMPLICIT NONE

   ! Integers
   INTEGER :: i

   ! Related to grid variables
   INTEGER :: grid_x, grid_y
   REAL (DP) :: dgrid_x, dgrid_y

   ! Related to hydro variables !
   REAL (DP) :: rad, vel_mag
   REAL (DP) :: eps_old, eps_new, temp_new

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !IF(debug_flag == 1) WRITE(*,*) 'In Getdata_p'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign !
   DO i = 1, nop, 1

      grid_x = grid_p(i,1)
      grid_y = grid_p(i,2)

      dgrid_x = dgrid_p(i,1)
      dgrid_y = dgrid_p(i,2)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !IF(i == 1) WRITE(*,*) i, grid_x, grid_y
      !if(grid_x < 0 .or. grid_x > length_step_r) write(*,*) 'Found x', i, grid_x, pos_p(i,1), vel_p(i,1)
      !if(grid_y < 0 .or. grid_y > length_step_z) write(*,*) 'Found y', i, grid_y, pos_p(i,2), vel_p(i,2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(leavebox_p(i) == .false.) THEN

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !IF(i == 1) WRITE(*,*) rho2(grid_x,grid_y), temp2(grid_x,grid_y)
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if((grid_x > 0 .and. grid_x < length_step_r_2 .and. &
   	     grid_y > 0 .and. grid_y < length_step_z_2)) then

	    ! Grid value !
            temp_p(i) = temp2(grid_x,grid_y)
            rho_p(i) = rho2(grid_x,grid_y)
	    phi_p(i) = phi2(grid_x,grid_y)

	    ! Interpolation !
      	    temp_p(i) = temp2(grid_x-1,grid_y-1) * (1.0D0 - dgrid_x) * (1.0D0 - dgrid_y) + &
      		  	temp2(grid_x,grid_y-1) * dgrid_x * (1.0D0 - dgrid_y) + &
      		  	temp2(grid_x-1,grid_y) * (1.0D0 - dgrid_x) * dgrid_y + &
      		  	temp2(grid_x,grid_y) * dgrid_x * dgrid_y
      	    rho_p(i) = rho2(grid_x-1,grid_y-1) * (1.0D0 - dgrid_x) * (1.0D0 - dgrid_y) + &
                       rho2(grid_x,grid_y-1) * dgrid_x * (1.0D0 - dgrid_y) + &
                       rho2(grid_x-1,grid_y) * (1.0D0 - dgrid_x) * dgrid_y + &
                       rho2(grid_x,grid_y) * dgrid_x * dgrid_y
      	    phi_p(i) = phi2(grid_x-1,grid_y-1) * (1.0D0 - dgrid_x) * (1.0D0 - dgrid_y) + &
                       phi2(grid_x,grid_y-1) * dgrid_x * (1.0D0 - dgrid_y) + &
                       phi2(grid_x-1,grid_y) * (1.0D0 - dgrid_x) * dgrid_y + &
                       phi2(grid_x,grid_y) * dgrid_x * dgrid_y
         else

	    leavebox_p(i) = .true.

      	    if(grid_x > length_step_r_2) then
      	 	grid_x = length_step_r_2
      	 	dgrid_x = 0.0D0
      	    endif
            if(grid_y > length_step_z_2) then
      	 	grid_y = length_step_z_2
      	  	dgrid_y = 0.0D0
            endif

            temp_p(i) = temp2(grid_x,grid_y)
            rho_p(i) = rho2(grid_x,grid_y)
	    phi_p(i) = phi2(grid_x,grid_y)
            vel_p(i,1) = vel2_r(grid_x,grid_y)
            vel_p(i,2) = vel2_z(grid_x,grid_y)
	 
	 endif

      ELSE

	 ! Assume homologous expansion

      ENDIF

   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !IF(debug_flag == 1) WRITE(*,*) 'Done getdata_p!'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   END SUBROUTINE getdata_p

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine get the velocity of each particles 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getvel_p
   use definition
   implicit none

   ! Grid variables !
   integer :: i
   integer :: grid_x, grid_y
   real (DP) :: dgrid_x, dgrid_y

   ! Related to homologus expansion
   REAL (DP) :: rad

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In GetVel_p'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign 
   do i = 1, nop, 1

      grid_x = grid_p(i,1)
      grid_y = grid_p(i,2)

      dgrid_x = dgrid_p(i,1)
      dgrid_y = dgrid_p(i,2)

      IF(leavebox_p(i) == .false.) THEN
  
         IF((grid_x > 0 .and. grid_x < length_step_r_2) .and. &
             (grid_y > 0 .and. grid_y < length_step_z_2)) THEN

	    ! Grid value !
            !vel_p(i,1) = vel2_r(grid_x,grid_y)
            !vel_p(i,2) = vel2_z(grid_x,grid_y)

            ! Interpolation !
      	    vel_p(i,1) = vel2_r(grid_x-1,grid_y-1) * (1.0D0 - dgrid_x) * (1.0D0 - dgrid_y) + &
                   	 vel2_r(grid_x,grid_y-1) * dgrid_x * (1.0D0 - dgrid_y) + &
                   	 vel2_r(grid_x-1,grid_y) * (1.0D0 - dgrid_x) * dgrid_y + &
                   	 vel2_r(grid_x,grid_y) * dgrid_x * dgrid_y
      	    vel_p(i,2) = vel2_z(grid_x-1,grid_y-1) * (1.0D0 - dgrid_x) * (1.0D0 - dgrid_y) + &
                   	 vel2_z(grid_x,grid_y-1) * dgrid_x * (1.0D0 - dgrid_y) + & 
                   	 vel2_z(grid_x-1,grid_y) * (1.0D0 - dgrid_x) * dgrid_y + & 
                   	 vel2_z(grid_x,grid_y) * dgrid_x * dgrid_y

         ELSE

            leavebox_p(i) = .true.

	    IF(coordinate_flag == 1) THEN
               rad = DSQRT(pos_p(i,1)**2 + pos_p(i,2)**2)
               vel_p(i,1) = vel_p(i,1) * pos_p(i,1) / rad
               vel_p(i,2) = vel_p(i,2) * pos_p(i,2) / rad
 	    ENDIF
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    !ELSEIF(coordinate_flag == 2) THEN
	    !   vel_p(i,1) = vel_p(i,1)	
	    !   vel_p(i,2) = 0.0D0
	    !ENDIF	
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ENDIF

      ENDIF

   ENDDO

   END SUBROUTINE getvel_p

end module PPT_module
