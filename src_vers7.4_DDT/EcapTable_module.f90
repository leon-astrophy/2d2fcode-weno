!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all the tools for calculating
! the electron capture rate assuming NSE composition
! Written by Leung Shing Chi in 2016
!
! This module contains the following subroutines
! You can directly go there by searching
!
! subroutine readEcapRate
! subroutine GetECapRate 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ecaptable_module
use definition, only : DP
implicit none

! The size of the E-cap table
integer, parameter :: rho_row = 23
integer, parameter :: temp_row = 48
integer, parameter :: ye_row = 122

! The electron capture rate table
real (DP), dimension(0:ye_row, 0:temp_row, 0:rho_row) :: ecaptable

! The neutrino energy loss rate table
real (DP), dimension(0:ye_row, 0:temp_row, 0:rho_row) :: eneutable

contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine reads the E-cap table reported in 
   ! Seitenzahl(2010). The rate used the newest weak
   ! rates and assumes NSE composition
   ! Written by Leung Shing Chi in 2016  
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine readEcapRate
   use definition, only : DP
   implicit none

   ! Dummy variables
   integer :: i, j, k

   ! Dummy
   real (DP) :: dummy

   ! Read-in data from file
   real (DP) :: input1, input2, input3, input4, input5

   ! Read in the low-Ye patch
   open(unit=441, file='../lib/ecaptable_Ye040_wNabi_full.dat', action='read')
   do i = 0, ye_row, 1  !0 - 15
      do j = 0, temp_row, 1
         do k = 0, rho_row, 1       
            read(441,*) dummy, dummy, input1, input2, input3, dummy, input4, input5
            eneutable(i,j,k) = input5 / 9.0D20 * 4.9282D-6
            ecaptable(i,j,k) = input4 * 4.9282D-6
            !bindetable(i,j,k) = 1.07227D-3 * input2 * input4 / input1
         enddo
      enddo
   enddo
   close(441)

   !open(unit=441, file='ecaptable_Ye040.dat', action='read')
   !do i = 0, 15, 1  !0 - 15
   !   do j = 0, temp_row, 1
   !      do k = 0, rho_row, 1
   !         read(441,*) dummy, dummy, input1, input2, input3, dummy, input4, input5
   !         eneutable(i,j,k) = input5 / 9.0D20 * 4.9282D-6
   !         ecaptable(i,j,k) = input4 * 4.9282D-6
   !         !bindetable(i,j,k) = 1.07227D-3 * input2 * input4 / input1
   !      enddo
   !   enddo
   !enddo
   !close(441)

   ! Read in the high-Ye Table
   !open(unit=441, file='mmc2.txt', action='read')
   !do i = 16, ye_row, 1
   !   do j = 0, temp_row, 1
   !      do k = 0, rho_row, 1
   !         read(441,*) dummy, dummy, input1, input2, input3, dummy, input4, input5
   !	    eneutable(i,j,k) = input4 / 9.0D20 * 4.9282D-6
   !         ecaptable(i,j,k) = input5 * 4.9282D-6
   !         !bindetable(i,j,k) = 1.07227D-3 * input2 * input4 / input1
   !      enddo
   !   enddo
   !enddo
   !close(441)

   end subroutine ReadEcapRate

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine takes the density, temperature 
   ! and Ye as input and give the corresponding 
   ! electron capture rate and neutrino energy loss rate
   ! using the trilinear interpolation
   ! Written by Leung Shing Chi in 2016             
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   subroutine GetECapRate(rho_in, temp_in, ye_in, ecaprate_out, eneurate_out)
   use definition, only : DP
   implicit none

   integer :: printflag 

   ! Input quantities
   real (DP) :: rho_in, temp_in, ye_in
 
   ! Output quantities
   real (DP) :: ecaprate_out
   real (DP) :: eneurate_out

   ! LOG-10 density for table entry
   real (DP) :: log10rho

   ! The corresponding grid number
   integer :: rho_grid, temp_grid, ye_grid

   ! The corresponding grid derivation
   real (DP) :: rho_dgrid, temp_dgrid, ye_dgrid

   log10rho = LOG10(rho_in * 6.171D17)                

   ! Make sure the input quantities is within the table limit
   if(log10rho > 7.0D0) then		! Minimum = 5.3

      if(temp_in > 5.0D0) then		! Minimum = 3.0

         if(ye_in > 0.20D0) then

	    ! First locate the grid in the table
            rho_grid = INT((log10rho - 5.3D0) / 0.25D0)
            temp_grid = INT((temp_in - 3.0D0) / 0.25D0)
            ye_grid = INT((ye_in - 0.200D0) / 0.0025D0)

	    ! Second locate how far the requested point from the grid boundary
            rho_dgrid = (log10rho - (5.3D0 + (DBLE(rho_grid) * 0.25D0))) / 0.25D0
            temp_dgrid = (temp_in - (3.0D0 + (DBLE(temp_grid) * 0.25D0))) / 0.25D0
            ye_dgrid = (ye_in - (0.20D0 + (DBLE(ye_grid) * 0.0025D0))) / 0.0025D0

	    if(rho_grid >= rho_row) then
	       rho_grid = rho_row
	       rho_dgrid = 0.0D0
            endif
	    if(temp_grid >= temp_row) then
	       temp_grid = temp_row
	       temp_dgrid = 0.0D0
	    endif
	    if(ye_grid >= ye_row) then
	       ye_grid = ye_row
	       ye_dgrid = 0.0D0
            endif

	    ! Some debug stuff
	    !if(printflag == 1) then
            ! write(*,*) rho_in, temp_in, ye_in
            ! write(*,*) rho_grid, temp_grid, ye_grid
            ! write(*,*) rho_dgrid, temp_dgrid, ye_dgrid
	    !endif

	    !Do the interpolation for the electron capture rate
            ecaprate_out =  ecaptable(ye_grid, temp_grid, rho_grid) + &
            (1.0D0 - rho_dgrid) * (1.0D0 - temp_dgrid) * ye_dgrid * &
            (ecaptable(ye_grid+1, temp_grid, rho_grid) - ecaptable(ye_grid, temp_grid, rho_grid)) + & 
             rho_dgrid * (1.0D0 - temp_dgrid) * (1.0D0 - ye_dgrid) * &
            (ecaptable(ye_grid, temp_grid, rho_grid+1) - ecaptable(ye_grid, temp_grid, rho_grid)) + &
            (1.0D0 - rho_dgrid) * temp_dgrid * (1.0D0 - ye_dgrid) * &
            (ecaptable(ye_grid, temp_grid+1, rho_grid) - ecaptable(ye_grid, temp_grid, rho_grid)) + &
            rho_dgrid * temp_dgrid * (1.0D0 - ye_dgrid) * &
            (ecaptable(ye_grid, temp_grid+1, rho_grid+1) - ecaptable(ye_grid, temp_grid, rho_grid)) + &
            rho_dgrid * (1.0D0 - temp_dgrid) * ye_dgrid * &
            (ecaptable(ye_grid+1, temp_grid, rho_grid+1) - ecaptable(ye_grid, temp_grid, rho_grid)) + &
            (1.0D0 - rho_dgrid) * temp_dgrid * ye_dgrid * &
            (ecaptable(ye_grid+1, temp_grid+1, rho_grid) - ecaptable(ye_grid, temp_grid, rho_grid)) + &
            rho_dgrid * temp_dgrid * ye_dgrid * &
            (ecaptable(ye_grid+1, temp_grid+1, rho_grid+1) - ecaptable(ye_grid, temp_grid, rho_grid))


	    ! Do the interpolation for the neutrino energy loss rate
	    eneurate_out =  eneutable(ye_grid, temp_grid, rho_grid) + &
            (1.0D0 - rho_dgrid) * (1.0D0 - temp_dgrid) * ye_dgrid * &
            (eneutable(ye_grid+1, temp_grid, rho_grid) - eneutable(ye_grid, temp_grid, rho_grid)) + & 
            rho_dgrid * (1.0D0 - temp_dgrid) * (1.0D0 - ye_dgrid) * &
            (eneutable(ye_grid, temp_grid, rho_grid+1) - eneutable(ye_grid, temp_grid, rho_grid)) + &
            (1.0D0 - rho_dgrid) * temp_dgrid * (1.0D0 - ye_dgrid) * &
            (eneutable(ye_grid, temp_grid+1, rho_grid) - eneutable(ye_grid, temp_grid, rho_grid)) + &
            rho_dgrid * temp_dgrid * (1.0D0 - ye_dgrid) * & 
            (eneutable(ye_grid, temp_grid+1, rho_grid+1) - eneutable(ye_grid, temp_grid, rho_grid)) + &
            rho_dgrid * (1.0D0 - temp_dgrid) * ye_dgrid * & 
            (eneutable(ye_grid+1, temp_grid, rho_grid+1) - eneutable(ye_grid, temp_grid, rho_grid)) + &
            (1.0D0 - rho_dgrid) * temp_dgrid * ye_dgrid * & 
            (eneutable(ye_grid+1, temp_grid+1, rho_grid) - eneutable(ye_grid, temp_grid, rho_grid)) + &
            rho_dgrid * temp_dgrid * ye_dgrid * & 
            (eneutable(ye_grid+1, temp_grid+1, rho_grid+1) - eneutable(ye_grid, temp_grid, rho_grid))

            ! To map it back in code unit
            !ecaprate_out = ecaprate_out

            !write(*,*) ecaprate_out
            !read(*,*)

         else

	   ! Zero output for invalid Ye
           ecaprate_out = 0.0D0
	   eneurate_out = 0.0D0

         endif

      else

	 ! Zero output for invalid temp
         ecaprate_out = 0.0D0
	 eneurate_out = 0.0D0

      endif

   else

      ! Zero output for invalid rho
      ecaprate_out = 0.0D0
      eneurate_out = 0.0D0

   endif

   end subroutine GetEcapRate

end module ecaptable_module