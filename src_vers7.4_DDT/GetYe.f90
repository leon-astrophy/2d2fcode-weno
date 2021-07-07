!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine assign the initial Ye assuming the matter is 
! in its equilibrium (density-dependent) Ye, based
! on the model by Liebendoerfer (2005a)
! Written by Leung Shing Chi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getye
use definition
implicit none

! Dummy variables
integer :: i, j, k

! The scaled density
real (DP) :: xrho

! The expected Ye
real (DP) :: ye_new

! Assign Ye for normal matter
do k = length_step_z_min_part_2, length_step_z_part_2, 1
   do j = 1, length_step_r_part_2, 1

      xrho = MAX(-1.0D0, MIN(1.0D0, LOG(rho2(j,k)**2 / 1.5756D-15) / 13.410D0))
      ye_new = 0.389D0 - 0.111D0 * xrho + 0.035D0 * (1.0D0 - ABS(xrho) + &
               4.0D0 * ABS(xrho) * (ABS(xrho) - 0.5D0) * (ABS(xrho) - 1.0D0))

      ye2(j,k) = ye_new

   enddo
enddo

! Boundaries 
call boundary1D_NM(ye2, even)

end subroutine getye