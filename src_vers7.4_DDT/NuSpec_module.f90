!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all the tools for calculating the 
! neutrino spectra for some given temperature and density
! It includes two major emission channels:
! The Plasmaneutrino and pair annihilation.
! Refer to Odrzywolek (2007) and Misiaszek (2006) 
! for more details
!
! Written by Leung Shing Chi in 2016
!
! The module contains the following subroutines:
! 1. subroutine FindNuSpec
! 2. subroutine read_nutable
! 3. subroutine Output_NuPhi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module NuSpec_module
use definition, only : DP
implicit none

! Size of table
INTEGER, PARAMETER :: temp_rowno3 = 10
INTEGER, PARAMETER :: den_rowno3 = 30    

! The effective electron mass
REAL (DP), DIMENSION(temp_rowno3, den_rowno3):: nutable_mass

! The calculated emissivities
REAL (DP), DIMENSION(temp_rowno3, den_rowno3):: nutable_emiss

! The neutrino spectra
real (DP) :: nu_phi(1:10)

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine interpolates the neutrino table to 
   ! Get the effective mass and the emissivities
   !
   ! Written by Leung Shing Chi in 2016 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine findnuspec
   use definition
   use helmeos_module
   implicit none

   ! dummy variables
   integer :: i, j, k

   ! Target grid in the table
   integer :: itemp, irho

   ! Derivation of the actual point from the grid boundary
   real (DP) :: dtemp, drho

   ! transverse mass
   real (DP) :: mass_t
   
   ! Neutrino emissivities
   real (DP) :: emiss

   ! Neutrino frequence
   real (DP) :: nu_freq

   ! Initialization
   nu_phi = 0.0D0

   do k = length_step_z_min_part_2, length_step_z_part_2
      do j = 1, length_Step_r_part_2

         if(rho2(j,k) > 1.62D-11 .and. temp2(j,k) > 1.0D0) then

	    ! Do the standard bilinear interpolation 
            irho = INT((LOG10(rho2(j,k) * 6.171D17) - 7.0D0) / 0.1D0) + 1
            itemp = INT(temp2(j,k)) 

            drho = ((LOG10(rho2(j,k) * 6.171D17) - 7.0D0) / 0.1D0 + 1.0D0) - DBLE(irho)
            dtemp = temp2(j,k) - DBLE(itemp) 

     	    mass_t = nutable_mass(itemp, irho) + &
                     dtemp * (1.0D0 - drho) * (nutable_mass(itemp + 1, irho) - nutable_mass(itemp, irho)) + &
                     drho * (1.0D0 - dtemp) * (nutable_mass(itemp, irho + 1) - nutable_mass(itemp, irho)) + &
                     dtemp * drho * (nutable_mass(itemp + 1, irho + 1) - nutable_mass(itemp, irho))

   	    emiss = nutable_emiss(itemp, irho) + &
                    dtemp * (1.0D0 - drho) * (nutable_emiss(itemp + 1, irho) - nutable_emiss(itemp, irho)) + &
                    drho * (1.0D0 - dtemp) * (nutable_emiss(itemp, irho + 1) - nutable_emiss(itemp, irho)) + &
                    dtemp * drho * (nutable_emiss(itemp + 1, irho + 1) - nutable_emiss(itemp, irho))

            do i = 1, 10, 1

	       ! Fix the frequence
               !nu_freq = 10.0D0 ** (-1.0D0 + 0.1D0 * DBLE(i))
   	       nu_freq = DBLE(i)

	       ! Use the analytic form given in the article
	       ! 3.22474D15 = length unit^3 = (1.4774 km)^3
	       nu_phi(i) = nu_phi(i) + (4.5215D25 * emiss * 0.1425776426D0 / (0.08614D0 * temp2(j,k)) * &
		   	   (nu_freq / 0.08614D0 / temp2(j,k)) ** 3.180657028D0 * &
		           EXP(-1.108192299D0 * nu_freq / 0.08614D0 / temp2(j,k)) + &
			   2.115D30 * (0.08614D0 * temp2(j,k)) * (mass_t * 0.511D0 ** 2) ** 3 * &
			   EXP(-nu_freq / 0.08614D0 / temp2(j,k))) * vol2(j,k) * 3.22474D15  
		
            enddo

         endif

      enddo
   enddo

   end subroutine findnuspec

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 
   ! This subroutine reads in the nu_mass.dat and 
   ! nu_emiss.dat table. 
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_nutable
   implicit none

   ! dummy variables
   integer :: i, j

   ! open the files
   open(unit=898, file='../lib/nu_mass.dat', action='read')
   open(unit=899, file='../lib/nu_emiss.dat', action='read')

   ! read the table
   do i = 1, temp_rowno3, 1
      read(898,*) (nutable_mass(i,j), j=1,den_rowno3,1)
      read(899,*) (nutable_emiss(i,j), j=1,den_rowno3,1)
   enddo

   close(898)
   close(899)

   end subroutine read_nutable

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! The subroutine output the calculated spectra
   ! Written by Leung Shing Chi
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_nuphi
   use definition, only : global_time
   implicit none

   ! dummy variables
   integer :: j

   ! output the data
   open(unit=601, file='./Outfile/Neutrino/Star_WENO_NuSpec_0.dat', action='write', position='append')
   write(601, 702) global_time, (nu_phi(j), j = 1, 10, 1)
   close(601)

   702 FORMAT (F23.15, 30ES18.10)

   end subroutine output_nuphi

end module NuSpec_module