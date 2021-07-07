!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine finds the adiabatic index for both DM and NM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDGAMMA
USE DEFINITION
IMPLICIT NONE

! Integer and dummy variables !
integer :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do for DM !
if(runDM_flag == 1) then
   do k = length_step_z_min_part_1, length_step_z_part_1, 1
      do j = 1, length_step_r_part_1, 1
         gamma_1(j,k) = rho1(j,k)*cs1(j,k)**2/p1(j,k)
      enddo
   enddo

   ! Boundaries !
   CALL BOUNDARY1D_DM(gamma_1)
endif

! For NM !
do k = length_step_z_min_part_2, length_step_z_part_2, 1
   do j = 1, length_step_r_part_2, 1
         gamma_2(j,k) = rho2(j,k)*cs2(j,k)**2/p2(j,k)
   enddo
enddo

! Boundaries !
CALL BOUNDARY1D_NM (gamma_2)

end subroutine findgamma