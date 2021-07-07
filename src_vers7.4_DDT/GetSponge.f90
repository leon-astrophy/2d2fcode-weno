!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine supply a sponge on the stellar !
! surface to avoid spurious oscillation due to   !
! the imperfect initial model			 !
! Written by Leung Shing Chi in 2016		 !
! based on the Zingale 2002			 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getsponge
use definition
implicit none

! Dummy variables
integer :: i, j, k

! Position factor defined in the paper
real (DP) :: r_tp, r_md, r_sp, rad_dist
real (DP) :: f_damp

! Threshold density and damping factor
real (DP) :: rho1_sponge = 1.0D-13
real (DP) :: rho2_sponge = 1.0D-13
real (DP) :: kap_sponge = 4.9282D-3

If(rundm_flag == 1) THEN

	! Give r_tp and r_sp an arbitrarily large number 
	rho1_sponge = rho1(1,1) * 0.01D0

	r_md = 1.0D10
	r_sp = 1.0D10

	! First find r_sp according to the prescription in paper
	do k = 1, length_step_z_1, 1
		do j = 1, length_step_r_1, 1
			rad_dist = DSQRT(r2(j)**2 + z2(k)**2)
			if(rho1(j,k) <= rho1_sponge) then
				r_sp = MIN(rad_dist, r_sp)
			endif
			if(rho1(j,k) <= rho1_a) then
				r_md = MIN(rad_dist, r_md)
			endif
		enddo
	enddo

	r_tp = 2.0D0 * r_md - r_sp
	! Then assign a relaxed veloicty according
	! to the damping strength
	do k = 1, length_step_z_1, 1
		do j = 1, length_step_r_1, 1
			rad_dist = DSQRT(r2(j)**2 + z2(k)**2)
			if(rad_dist > r_tp) then
				f_damp = 1.0D0
			elseif(rad_dist <= r_tp .and. rad_dist > r_sp) then
				f_damp = 0.5D0 * (1.0D0 - DCOS(pi * (rad_dist - r_sp)/(r_tp - r_sp)))
			else
				f_damp = 0.0D0
			endif
			vel1_r(j,k) = vel1_r(j,k) / (1.0D0 + kap_sponge * f_damp * dt)
			vel1_z(j,k) = vel1_z(j,k) / (1.0D0 + kap_sponge * f_damp * dt) 
		enddo
	enddo

END IF

! Give r_tp and r_sp an arbitrarily large number 
rho2_sponge = rho2(1,1) * 0.01D0

r_md = 1.0D10
r_sp = 1.0D10

! First find r_sp according to the prescription in paper
do k = 1, length_step_z_2, 1
	do j = 1, length_step_r_2, 1
		rad_dist = DSQRT(r2(j)**2 + z2(k)**2)
		if(rho2(j,k) <= rho2_sponge) then
			r_sp = MIN(rad_dist, r_sp)
		endif

		if(rho2(j,k) <= rho2_a) then
			r_md = MIN(rad_dist, r_md)
		endif
	enddo
enddo

r_tp = 2.0D0 * r_md - r_sp
! Then assign a relaxed veloicty according
! to the damping strength
do k = 1, length_step_z_2, 1
	do j = 1, length_step_r_2, 1

		rad_dist = DSQRT(r2(j)**2 + z2(k)**2)

		if(rad_dist > r_tp) then
			f_damp = 1.0D0
		elseif(rad_dist <= r_tp .and. rad_dist > r_sp) then
			f_damp = 0.5D0 * (1.0D0 - DCOS(pi * (rad_dist - r_sp)/(r_tp - r_sp)))
		else
			f_damp = 0.0D0
		endif

		vel2_r(j,k) = vel2_r(j,k) / (1.0D0 + kap_sponge * f_damp * dt)
		vel2_z(j,k) = vel2_z(j,k) / (1.0D0 + kap_sponge * f_damp * dt) 

	enddo
enddo

end subroutine getsponge