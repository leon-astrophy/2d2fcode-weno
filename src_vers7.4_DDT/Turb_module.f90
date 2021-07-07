!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! New variables for turbulent flame used when turbulence_flag == 1
!
! This module contains the components for calculating
! sub-grid scale (SGS) turbulence using the 
! one-equation form (See Niemeyer1995b). 
! This assumes the turbulence can be characterized by the 
! local velocity flucation v', and the associated 
! turbulent kinetic energy q.
! By solving the evolution of q, one can approximate
! the local eddy motion, which can be used for 
! 
! Written by Leung Shing Chi in 2015
! Documented by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE turb_module
USE definition
IMPLICIT NONE

! Equation index !
INTEGER :: iturbq

! Atmospheric values !
REAL (DP) :: turb_q_a = 1.0D-10
REAL (DP) :: turb_qtotal
REAL (DP) :: turb_q_max

! Output arrays
REAL (DP), allocatable, dimension(:,:) :: turb_source
REAL (DP), allocatable, dimension(:,:) :: turb_divergence
REAL (DP), allocatable, dimension(:,:) :: turb_diff_r
REAL (DP), allocatable, dimension(:,:) :: turb_diff_z

! The k-eps component
REAL (DP), allocatable, dimension(:,:) :: turb_q
REAL (DP), allocatable, dimension(:,:) :: turb_eps

CONTAINS

   SUBROUTINE buildTurb
   USE definition
   IMPLICIT NONE
   
   ! Source terms !
   ALLOCATE(turb_source(-4:length_step_r_2+5, -4:length_step_z_2+5))
   ALLOCATE(turb_divergence(-4:length_step_r_2+5, -4:length_step_z_2+5))
   ALLOCATE(turb_diff_r(-4:length_step_r_2+5, -4:length_step_z_2+5))
   ALLOCATE(turb_diff_z(-4:length_step_r_2+5, -4:length_step_z_2+5))

   ! Turbulence energy !
   ALLOCATE(turb_q(-4:length_step_r_2+5, -4:length_step_z_2+5))
   ALLOCATE(turb_eps(-4:length_step_r_2+5, -4:length_step_z_2+5))

   END SUBROUTINE buildTurb

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE destroyTurb
   USE definition
   IMPLICIT NONE

   ! Source terms !
   DEALLOCATE(turb_source)
   DEALLOCATE(turb_divergence)
   DEALLOCATE(turb_diff_r)
   DEALLOCATE(turb_diff_z)

   ! Turbulence energy !
   DEALLOCATE(turb_q)
   DEALLOCATE(turb_eps)

   END SUBROUTINE destroyTurb

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GetTurb
   USE definition
   IMPLICIT NONE

   ! Integer !
   INTEGER :: j, k

   ! This is to initialize the sub-grid turbulence-dissipation field
   ! First give the turbulent kinetic energy field a very small number
   turb_q = turb_q_a

   END SUBROUTINE GetTurb

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE outputturb_log
   USE definition
   IMPLICIT NONE

   OPEN(unit=401, FILE='./Outfile/Turbulence/Star_WENO_TotalQ_0.dat', ACTION='write', POSITION='append')
   WRITE (401, 701) global_time, turb_qtotal
   CLOSE(401)

   OPEN(unit=401, FILE='./Outfile/Turbulence/Star_WENO_MaxQ_0.dat', ACTION='write', POSITION='append')
   WRITE (401, 701) global_time, turb_q_max
   CLOSE(401)

   ! Format !
   701 FORMAT (F23.15, 2ES18.10)

   END SUBROUTINE outputturb_log

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE outputturb_profile(n)
   USE definition
   IMPLICIT NONE
 
   ! Integer !
   INTEGER :: j, k, n

   ! File index !
   INTEGER :: fileno_len
   CHARACTER (len = 256) :: fileno
   
   ! File name !
   WRITE (fileno, *) n
   fileno = ADJUSTL (fileno)
   fileno_len = LEN_TRIM (fileno)

   open(unit=401, file='./Outfile/Turbulence/Star_WENO_TurbQ_'//fileno(1:fileno_len)//'.dat', action='write', position='append')
   WRITE (401, *) n, global_time              
   WRITE (401, *) length_step_r_2, length_step_z_2
   WRITE (401, *) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (401, 701) (turb_q(j,k), j=1, length_step_r_2)
   END DO
   WRITE (401, *)
   CLOSE(401)

   ! Format !
   701 FORMAT (2000ES13.5)

   END SUBROUTINE outputturb_profile

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE findturbenergy
   USE definition
   IMPLICIT NONE

   ! Integer !
   INTEGER :: j, k

   ! Initialize !
   turb_qtotal = 0.0D0
   turb_q_max = 0.0D0

   ! Do the loop !
   DO k = length_step_z_min_part_2, length_step_z_part_2
      DO j = 1, length_step_r_part_2
         IF(rho2(j,k) > rho2_a) THEN
            turb_qtotal = turb_qtotal +  vol2(j,k) * rho2(j,k) * turb_q(j,k)
            turb_q_max = MAX(turb_q(j,k), turb_q_max)
         ENDIF
      ENDDO
   ENDDO

   END SUBROUTINE findturbenergy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This subroutine finds all component for the sub-grid turbulence

   subroutine findturbulence
   use definition
   use levelset_module
   implicit none

   ! Integer 
   integer :: j, k

   ! Real number 
   real (DP) :: turb_C, turb_D, turb_F
   real (DP) :: at_no, g_eff
   real (DP) :: dx_eff
   real (DP) :: inv_sqrt2 = 1.0D0 / DSQRT(2.0D0)
   real (DP) :: c_lambda = 1.0D0

   ! For Archimedis Production
   real (DP), allocatable, dimension(:,:) :: turb_RT

   ! For turbulence compression
   real (DP), allocatable, dimension(:,:) :: turb_comp

   ! For turbulence production
   real (DP), allocatable, dimension(:,:) :: turb_div_v
   real (DP), allocatable, dimension(:,:) :: turb_str
   real (DP), allocatable, dimension(:,:) :: turb_eta
   real (DP), allocatable, dimension(:,:,:,:) :: turb_velgrad

   ! For turbulence diffusion
   real (DP), allocatable, dimension(:,:,:) :: turb_qdev

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !if(debug_flag == 1) write(*,*) 'In Find turbulence'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Allocate arrays !
   allocate(turb_RT(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(turb_comp(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(turb_div_v(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(turb_str(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(turb_eta(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(turb_velgrad(-4:length_step_r_2+5, -4:length_step_z_2+5,2,2))
   allocate(turb_qdev(-4:length_step_r_2+5, -4:length_step_z_2+5,2))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize !
   turb_source = 0.0D0
   turb_divergence = 0.0D0
   turb_diff_r = 0.0D0
   turb_diff_z = 0.0D0
   turb_eps = 0.0D0

   ! Initialize !
   turb_RT = 0.0D0
   turb_comp = 0.0D0
   turb_div_v = 0.0D0
   turb_str = 0.0D0
   turb_eta = 0.0D0
   turb_velgrad = 0.0D0
   turb_qdev = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   dx_eff = dx2 * 1.6D0 

   if(flame_flag == 1) then

      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1

            turb_F = MIN(100.0D0, MAX(0.1D0, 1.0D-4 * epsilon2(j,k) / turb_q(j,k)))
            turb_C = 0.1D0 * turb_F
            turb_D = 0.5D0 / turb_F

            !turb_eta(j,k) = rho2(j,k) * turb_C * dx_eff * DSQRT(turb_q(j,k))
            !turb_eps(j,k) = rho2(j,k) * turb_D * (turb_q(j,k) ** 1.5D0) / dx_eff
	    turb_eta(j,k) = rho2(j,k) * turb_C * dx_eff * DSQRT(DSQRT(turb_q(j,k)**2))
            turb_eps(j,k) = rho2(j,k) * turb_D * DSQRT(DSQRT(turb_q(j,k)**6)) / dx_eff

            if(flamegrid_flag(j,k) /= 0 .and. flamegrid_flag(j,k) /= 1) then
               at_no = MAX(0.5D0 * (0.0522D0 + 0.145D0 / DSQRT(rho2(j,k)/1.62D-9) - 0.01D0 / (rho2(j,k)/1.62D-9)), 0.0D0)
               g_eff = at_no * DSQRT(phi2_r(j,k)**2 + phi2_z(j,k)**2)
               !turb_RT(j,k) = 0.625 * rho2(j,k) * DSQRT(turb_q(j,k)) * g_eff
	       turb_RT(j,k) = 0.625D0 * rho2(j,k) * DSQRT(DSQRT(turb_q(j,k)**2)) * g_eff
            else
               turb_RT(j,k) = 0.0D0
            endif

         enddo
      enddo

   else

      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1

            turb_F = MIN(100.0D0, MAX(0.1D0, 1.0D-4 * epsilon2(j,k) / turb_q(j,k)))
            turb_C = 0.1D0 * turb_F
            turb_D = 0.5D0 / turb_F

            !turb_eta(j,k) = rho2(j,k) * turb_C * dx_eff * DSQRT(turb_q(j,k))
            !turb_eps(j,k) = rho2(j,k) * turb_D * (turb_q(j,k) ** 1.5D0) / dx_eff
	    turb_eta(j,k) = rho2(j,k) * turb_C * dx_eff * DSQRT(DSQRT(turb_q(j,k)**2))
            turb_eps(j,k) = rho2(j,k) * turb_D * DSQRT(DSQRT(turb_q(j,k)**6)) / dx_eff

            turb_RT(j,k) = 0.0D0

         enddo
      enddo

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get the component using the PPM scheme 
   !call splitppm(vel2_r(:,:), turb_velgrad(:,:,1,1), 1, 0)
   !call splitppm(vel2_r(:,:), turb_velgrad(:,:,1,2), 2, 3)
   !call splitppm(vel2_z(:,:), turb_velgrad(:,:,2,1), 1, 3)
   !call splitppm(vel2_z(:,:), turb_velgrad(:,:,2,2), 2, 0)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k = length_step_z_min_part_2 - 4, length_step_z_part_2 + 4, 1
      do j = -3, length_step_r_part_2 + 4, 1
   	 turb_velgrad(j,k,1,1) = 0.5D0 * (vel2_r(j+1,k) - vel2_r(j-1,k)) / dx2
   	 turb_velgrad(j,k,1,2) = 0.5D0 * (vel2_r(j,k+1) - vel2_r(j,k-1)) / dx2
   	 turb_velgrad(j,k,2,1) = 0.5D0 * (vel2_z(j+1,k) - vel2_z(j-1,k)) / dx2
   	 turb_velgrad(j,k,2,2) = 0.5D0 * (vel2_z(j,k+1) - vel2_z(j,k-1)) / dx2
      enddo
   enddo

   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1
         turb_div_v(j,k) = turb_velgrad(j,k,1,1) + turb_velgrad(j,k,2,2) + vel2_r(j,k) / r2(j)
         turb_comp(j,k) = c_lambda * rho2(j,k) * turb_q(j,k) * turb_div_v(j,k)
      enddo
   enddo

   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1
         turb_str(j,k) = turb_eta(j,k) * &
 			 (2.0D0*(turb_velgrad(j,k,1,1)**2 + turb_velgrad(j,k,2,2)**2 + (vel2_r(j,k)/r2(j))**2) + &
			 (turb_velgrad(j,k,1,2) + turb_velgrad(j,k,2,1))**2 - c_lambda*turb_div_v(j,k)**2)
			 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         !((turb_velgrad(j,k,1,1) - turb_velgrad(j,k,2,2)) ** 2 + &
                         !(turb_velgrad(j,k,1,2) + turb_velgrad(j,k,2,1)) ** 2 - & 
			 ! vel2_r(j,k) / r2(j) * (turb_velgrad(j,k,1,1) + turb_velgrad(j,k,2,2)))
			 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
   enddo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Find diffusion
   !call splitppm(turb_q(:,:), turb_qdev(:,:,1), 1, 1)
   !call splitppm(turb_q(:,:), turb_qdev(:,:,2), 2, 2)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Old scheme
   do k = length_step_z_min_part_2-1, length_step_z_part_2+1, 1
      do j = 0, length_step_r_part_2+1, 1
    	 turb_qdev(j,k,1) = (turb_q(j+1,k) - turb_q(j-1,k)) * 0.5D0 / dx2
    	 turb_qdev(j,k,2) = (turb_q(j,k+1) - turb_q(j,k-1)) * 0.5D0 / dx2
      enddo
   enddo

   ! Old scheme
   do k = length_step_z_min_part_2-1, length_step_z_part_2+1, 1
      do j = 0, length_step_r_part_2+1, 1
         turb_diff_r(j,k) = turb_eta(j,k) * turb_qdev(j,k,1)
         turb_diff_z(j,k) = turb_eta(j,k) * turb_qdev(j,k,2)
      enddo
   enddo
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !call boundary1D(turb_diff_r, oddR)
   !call boundary1D(turb_diff_z, oddZ)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Find source

   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1
   	 turb_divergence(j,k) = (turb_diff_r(j+1,k) - turb_diff_r(j-1,k)) / 2.0D0 / dx2 + &
	    		    (turb_diff_z(j,k+1) - turb_diff_z(j,k-1)) / 2.0D0 / dx2 + turb_diff_r(j,k)/r2(j)
      enddo
   enddo

   if(flame_flag == 0) then

      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1
   	    turb_source(j,k) = -turb_comp(j,k) + turb_str(j,k) - turb_eps(j,k)
         enddo
      enddo

   else

      DO k = length_step_z_min_part_2, length_step_z_part_2, 1
         DO j = 1, length_step_r_part_2, 1
	    IF(flamegrid_flag(j,k) == 1 .or. flamegrid_flag(j,k) == 0) THEN
               turb_source(j,k) = -turb_comp(j,k) + turb_str(j,k) - turb_eps(j,k)
   	    ELSE
	       turb_source(j,k) = turb_RT(j,k) - turb_eps(j,k)
   	    ENDIF
         ENDDO
      ENDDO

   ENDIF

   CALL boundary1D_NM(turb_source, even)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Deallocate !
   DEALLOCATE(turb_RT)
   DEALLOCATE(turb_comp)
   DEALLOCATE(turb_div_v)
   DEALLOCATE(turb_str)
   DEALLOCATE(turb_eta)
   DEALLOCATE(turb_velgrad)
   DEALLOCATE(turb_qdev)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   END SUBROUTINE findturbulence

END MODULE turb_module