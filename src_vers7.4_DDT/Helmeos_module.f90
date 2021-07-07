!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains all the subroutine related
! to chemical reaction and flame/detonation  physics
! Written by Leung Shing Chi in 2016
! 
! This module contains the following subroutine
! You can search the subroutine from here
! 1.  subroutine buildHelm
! 2.  subroutine destroyHelm
! 3.  subroutine initialize_network
! 4.  subroutine init_my_iso7
! 5.  subroutine read_nse_table
! 6.  subroutine GetNSEState
! 7.  SUBROUTINE FINDCENTRALTEMPERATURE
! 8.  subroutine findluminosity
! 9.  subroutine FINDTOTALNEUTRINOLOSS
! 10. subroutine OutputXiso_profile19
! 11. subroutine output_xmass
! 12. subroutine FLAME
! 13. subroutine FLAME_INI
! 14. subroutine DETON_INI
! 15. subroutine BURN
! 16. subroutine NSE
! 17. subroutine NSE2
! 18. subroutine compute_binde
! 19. subroutine FIND_AZBAR
! 20. subroutine private_helmeos_azbar
! 21. subroutine checkXisotope
! 22. subroutine findneutrinoloss
! 23. subroutine updateXIsotope (muted)
! 24. subroutine AZBAR
! 25. subroutine invert_helm_pt_quiet
! 26. subroutine constructXISO
! 27. subroutine burn_Phase1
! 28. subroutine burn_Phase2
! 29. subroutine checkBurnPhi
! 30. subroutine init_burnPhi
! 31. subroutine read_NSE_table2
! 32. subroutine getNSEState2
! 33. subroutine Find_AZBar2
! 34. subroutine burn_Phase1_ini
! 35. subroutine burn_Phase1b
! 36. subroutine burn_Phase1b_ini
! 37. subroutine burn_Phase2b
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! All variables and subroutines related to open-source Helmholtz EOS
! New variables designed for finite temperature EOS used only when helmeos_flag == 1 .or. fusion_flag == 1
module helmeos_module
use definition
implicit none

! Number of isotope
INTEGER, PARAMETER :: totalion = 7 !19 !204

! The last global time for burning
REAL (DP) :: last_burntime

! Density limit for nuclear burning
real (DP), parameter :: rho2_burn_max = 1.62D-11
real (DP), parameter :: rho2_burn_min = 8.10D-12

! Density limit for deflagration
real (DP), parameter :: rho2_flame_max = 1.62D-8
real (DP), parameter :: rho2_flame_min = 1.62D-11

! Density limit for detonation
real (DP), parameter :: rho2_deton_max = 1.62D-9
REAL (DP), parameter :: rho2_deton_min = 1.62D-13

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Some global quantities
! 1. Central temperature
! 2. Total thermal neutrino loss
! 3. Total non-thermal neutrino loss
! 4. Luminsoity
! 5. Luminosity by burning, def. and det.
REAL (DP) :: centraltemperature
REAL (DP) :: total_nu_qdot
REAL (DP) :: total_ecap_nu_qdot
REAL (DP) :: lumino
REAL (DP) :: lumino_burn, lumino_flame, lumino_deton

! MAss burned by deflgration/detonation
REAL (DP) :: burn_mass

! Initial chemical composition
REAL (DP) :: abar_ini, abar_ini1, abar_ini2
REAL (DP) :: zbar_ini, zbar_ini1, zbar_ini2
REAL (DP) :: ye2_ini, ye2_ini1, ye2_ini2

! Atmospheric composition !
REAL (DP), dimension(totalion) :: xiso_a, xiso_a1, xiso_a2
REAL (DP) :: yiso_a, qash_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! My iso7 prescription

integer :: ihe4, ic12, io16, ine20, img24, isi28, ini56
integer :: che4, cc12, co16, cne20, cmg24, csi28, cni56
character(len=5) :: ionam(1:totalion)

real (selected_real_kind(15,307)), dimension(1:totalion) :: aion, zion, bion, nion, mion, wion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! My iso13 extension

integer :: is32, iar36, ica40, iti44, icr48, ife52 
integer :: cs32, car36, cca40, cti44, ccr48, cfe52

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! My iso19 extension   

integer :: ih1, ihe3, in14, ife54, iprot, ineut                                             
integer :: ch1, che3, cn14, cfe54, cprot, cneut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! My burn progress variables

integer :: ibphi1, ibphi2, ibphi3, iyiso, iqash

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 1
! Table size
INTEGER, PARAMETER :: temp_rowno_nse = 70
INTEGER, PARAMETER :: den_rowno_nse = 30

! Binding energy and composition of NSE
REAL (DP), DIMENSION(0:den_rowno_nse, 0:temp_rowno_nse):: nsetable_binde
REAL (DP), DIMENSION(0:den_rowno_nse, 0:temp_rowno_nse, 1:totalion):: nsetable_xiso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 2
! Table size         
INTEGER, PARAMETER :: ent_rowno_nse2 = 50
INTEGER, PARAMETER :: ye_rowno_nse2 = 26
INTEGER, PARAMETER :: den_rowno_nse2 = 50

! New binding energy and composition of NSE with a larger network
REAL (DP), DIMENSION(0:den_rowno_nse2, 0:ye_rowno_nse2+1, 3):: nsetable2_head
REAL (DP), DIMENSION(0:den_rowno_nse2, 0:ye_rowno_nse2+1, 0:ent_rowno_nse2+1, 6):: nsetable2_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Section of NSE table 3                       
! Table size
INTEGER, PARAMETER :: temp_rowno_nse3 = 48     
INTEGER, PARAMETER :: ye_rowno_nse3 = 122      
INTEGER, PARAMETER :: den_rowno_nse3 = 23

! New binding energy and composition of NSE with a larger network
REAL (DP), DIMENSION(0:den_rowno_nse3+1, 0:ye_rowno_nse3+1, 0:temp_rowno_nse3+1, 6):: nsetable3_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocatable objects

! Flag for being in NSE state
! nse_flag = 0 means in C-burning 
! nse_flag = 1 means in O- and Mg- burning
! nse_flag = 2 means in NSE burning
INTEGER  , allocatable, DIMENSION (:,:) :: nse_flag

! Energy loss by neutrino
REAL (DP), allocatable, DIMENSION (:,:) :: nu_qdot

! Energy input by burning
real (DP), allocatable, dimension (:,:) :: burn_qdot

! Energy input by deflagration
real (DP), allocatable, dimension (:,:) :: flame_qdot

! Energy input by detonation
real (DP), allocatable, dimension (:,:) :: deton_qdot

! Mean atomic mass
real (DP), allocatable, dimension (:,:) :: abar2

! Mean atomic number
real (DP), allocatable, dimension (:,:) :: zbar2

! Mean atomic mass (for 1st over-layer)
real (DP), allocatable, dimension (:,:) :: abar3

! Mean atomic mass (for 1st over-layer)
real (DP), allocatable, dimension (:,:) :: zbar3

! Chemical composition (for 1st over-layer)
REAL (DP), allocatable, DIMENSION (:,:,:) :: xiso3

! Chemical composition and its backup
REAL (DP), allocatable, DIMENSION (:,:,:) :: xiso, xiso_old

! Dummy variables for time-evolution
REAL (DP), allocatable, DIMENSION (:,:,:) :: xiso1, delta_xiso

! The progress variables
REAL (DP), allocatable, dimension (:,:) :: burn_phi1
REAL (DP), allocatable, dimension (:,:) :: burn_phi2
REAL (DP), allocatable, dimension (:,:) :: burn_phi3
REAL (DP), allocatable, dimension (:,:) :: yiso
REAL (DP), allocatable, dimension (:,:) :: qash

CONTAINS
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine allocates the necessary variables 
   ! that will be used in this module.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine buildHelm
   use definition
   implicit none

   ! Allocate all the arrays 
   allocate(nse_flag(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(nu_qdot(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(burn_qdot(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(flame_qdot(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(deton_qdot(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(abar2(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(zbar2(-4:length_step_r_2+5, -4:length_step_z_2+5))
   allocate(xiso(-4:length_step_r_2+5, -4:length_step_z_2+5, totalion))

   if(xisotran_flag == 1) then
      allocate(xiso_old(-4:length_step_r_2+5, -4:length_step_z_2+5, totalion))
      allocate(xiso1(-4:length_step_r_2+5, -4:length_step_z_2+5, totalion))
      allocate(delta_xiso(-4:length_step_r_2+5, -4:length_step_z_2+5, totalion))
   endif

   if(burn_prog_flag == 1) then 
      allocate(burn_phi1(-4:length_step_r_2+5, -4:length_step_z_2+5))
      allocate(burn_phi2(-4:length_step_r_2+5, -4:length_step_z_2+5))
      allocate(burn_phi3(-4:length_step_r_2+5, -4:length_step_z_2+5))
      allocate(yiso(-4:length_step_r_2+5, -4:length_step_z_2+5))
      allocate(qash(-4:length_step_r_2+5, -4:length_step_z_2+5))
   endif

   end subroutine buildHelm

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine deallocates all assigned array to 
   ! free the memory.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE destroyHelm
   USE definition
   IMPLICIT NONE

   ! Deallocate them
   DEALLOCATE(nse_flag)
   DEALLOCATE(nu_qdot)
   DEALLOCATE(burn_qdot)
   DEALLOCATE(flame_qdot)
   DEALLOCATE(deton_qdot)
   DEALLOCATE(abar2)
   DEALLOCATE(zbar2)
   DEALLOCATE(xiso)

   IF(xisotran_flag == 1) THEN
      DEALLOCATE(xiso_old)
      DEALLOCATE(xiso1)
      DEALLOCATE(delta_xiso)
   ENDIF

   IF(burn_prog_flag == 1) THEN
      DEALLOCATE(burn_phi1)  
      DEALLOCATE(burn_phi2)
      DEALLOCATE(burn_phi3)
      DEALLOCATE(yiso)
      DEALLOCATE(qash)
   ENDIF

   END SUBROUTINE destroyHelm

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calls all network related subroutines.
   ! It includes the setting up of isotopes, NSE table,
   ! deflagration and detonation table and electron capture rate
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE initialize_network()
   USE flametable_module
   USE ecaptable_module
   IMPLICIT NONE

   CALL init_my_iso7
   CALL read_nse_table
   CALL read_nse_table2
   CALL read_flame_table
   CALL read_deton_table
   CALL readEcapRate

   END SUBROUTINE initialize_network

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine initialize all variables related to isotope network
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init_my_iso7
   USE definition, only : DP
   IMPLICIT NONE

   INTEGER :: i
   REAL (DP) :: mev2erg,mev2gr
   REAL (DP), PARAMETER :: ev2erg = 1.60217648740d-12
   REAL (DP), PARAMETER :: clight  = 2.99792458d10
   REAL (DP), PARAMETER :: avo     = 6.0221367d23
   PARAMETER        (mev2erg = ev2erg*1.0d6, &
                     mev2gr  = mev2erg/clight**2)

   burn_mass = 0.0D0
 
   ! set the id numbers of the elements
   che4  = 1
   cc12  = 2
   co16  = 3
   cne20 = 4
   cmg24 = 5
   csi28 = 6
   cni56 = 7

   ! set the names of the elements
   ionam(che4)  = 'he4 '
   ionam(cc12)  = 'c12 '
   ionam(co16)  = 'o16 '
   ionam(cne20) = 'ne20'
   ionam(cmg24) = 'mg24'
   ionam(csi28) = 'si28'
   ionam(cni56) = 'ni56'

   ! set the number of nucleons in the element
   aion(che4)  = 4.0d0
   aion(cc12)  = 12.0d0
   aion(co16)  = 16.0d0
   aion(cne20) = 20.0d0
   aion(cmg24) = 24.0d0
   aion(csi28) = 28.0d0
   aion(cni56) = 56.0d0

   ! set the number of protons in the element
   zion(che4)  = 2.0d0
   zion(cc12)  = 6.0d0
   zion(co16)  = 8.0d0
   zion(cne20) = 10.0d0
   zion(cmg24) = 12.0d0
   zion(csi28) = 14.0d0
   zion(cni56) = 28.0d0

   ! set the binding energy of the element
   bion(che4)  =  28.29603d0
   bion(cc12)  =  92.16294d0
   bion(co16)  = 127.62093d0
   bion(cne20) = 160.64788d0
   bion(cmg24) = 198.25790d0
   bion(csi28) = 236.53790d0
   bion(cni56) = 484.00300d0

   ! set the number of neutrons and mass
   DO i = 1, totalion, 1
      nion(i) = aion(i) - zion(i)
   ENDDO

   ! mass of each isotope
   DO i = 1, totalion, 1
      mion(i) = nion(i)*1.67492721184d-24 + zion(i)*1.67262163783d-24 - bion(i)*mev2gr
   ENDDO

   ! molar mass      
   DO i = 1, totalion, 1
    wion(i) = avo * mion(i)             
   ENDDO

   ! a common approximation   
   DO i = 1, totalion, 1
      wion(i) = aion(i)       
   ENDDO

   END SUBROUTINE init_my_iso7

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine reads the pre created NSE table 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE read_nse_table
   USE definition, only : DP
   IMPLICIT NONE

   INTEGER :: i, j, k
   REAL (DP) :: dummy

   ! Open file !
   OPEN(unit=500, file='../lib/nse_table_7iso.dat',action='read')

   ! Read it !
   DO i = 0, den_rowno_nse, 1
      DO j = 0, temp_rowno_nse, 1
         READ(500,*) dummy, dummy, nsetable_binde(i,j), (nsetable_xiso(i,j,k), k = 1, totalion)
         nsetable_binde(i,j) = nsetable_binde(i,j) / 9.0D20
      ENDDO
   ENDDO

   ! Close !
   CLOSE(500)

   END SUBROUTINE read_nse_table

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subrouinte reads and gets the nse state - isotope mass tractions
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GetNSEState(rho_in, temp_in, xiso_nse_out)
   USE definition, only : DP
   IMPLICIT NONE

   REAL (DP) :: rho_in, temp_in
   REAL (DP) :: xiso_nse_out(1:totalion)

   REAL (DP) :: log10rho

   INTEGER :: rho_grid, temp_grid, ye_grid
   REAL (DP) :: rho_dgrid, temp_dgrid, ye_dgrid

   log10rho = LOG10(rho_in * 6.171D17)                

   IF(log10rho > 7.0D0) THEN            ! Minimum = 5.3
      if(temp_in > 5.0D0) THEN          ! Minimum = 3.0 and maximum 11.0

         rho_grid = INT((log10rho - 7.0D0) / 0.1D0)
         temp_grid = INT((temp_in - 4.0D0) / 0.1D0)

         rho_dgrid = (log10rho - (7.0D0 + (DBLE(rho_grid) * 0.1D0))) / 0.1D0
         temp_dgrid = (temp_in - (4.0D0 + (DBLE(temp_grid) * 0.1D0))) / 0.1D0

	 IF(rho_grid >= 30) THEN
            rho_grid = 30         
            rho_dgrid = 0.0D0
         ENDIF

         IF(temp_grid >= 70) THEN
            temp_grid = 70
            temp_dgrid = 0.0D0                    
         ENDIF

         xiso_nse_out(:) =  nsetable_xiso(rho_grid, temp_grid, :) + &
                         rho_dgrid * temp_dgrid * &
			(nsetable_xiso(rho_grid+1, temp_grid+1, :) - nsetable_xiso(rho_grid, temp_grid, :)) + & 
                         rho_dgrid * (1.0D0 - temp_dgrid) * &
			(nsetable_xiso(rho_grid+1, temp_grid, :) - nsetable_xiso(rho_grid, temp_grid, :)) + &
                         (1.0D0 - rho_dgrid) * temp_dgrid * &
			(nsetable_xiso(rho_grid, temp_grid+1, :) - nsetable_xiso(rho_grid, temp_grid, :)) 

      ENDIF
   ENDIF

   END SUBROUTINE GetNSEState

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the central temperature
   ! Written by Wong Ka Wing in 2010 (or before?)
   !
   ! Merged by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE FINDCENTRALTEMPERATURE
   USE DEFINITION
   IMPLICIT NONE

   centraltemperature = temp2(1,1)
  
   END SUBROUTINE findcentraltemperature

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the total luminosity
   ! by summing up all energy-generating source
   !
   ! Merged by Leung Shing Chi in 2016
   ! Updated by Leung Shing CHi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE findluminosity
   USE definition
   IMPLICIT NONE

   ! dummy variables
   INTEGER :: j, k

   ! local mass
   REAL (DP) :: dmass

   ! Initilization
   lumino_flame = 0.0D0
   lumino_deton = 0.0D0
   lumino_burn = 0.0D0
   
   ! Do the integration
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1
	 lumino_flame = lumino_flame + vol2(j,k) * rho2(j,k) * flame_qdot(j,k)
      ENDDO
   ENDDO

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1
	 lumino_deton = lumino_deton + vol2(j,k) * rho2(j,k) * deton_qdot(j,k)
      ENDDO
   ENDDO

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1
	 lumino_burn = lumino_burn + vol2(j,k) * rho2(j,k) * burn_qdot(j,k)
      ENDDO
   ENDDO

   ! Divide dt to get the time-rate
   lumino_flame = lumino_flame / dt
   lumino_deton = lumino_deton / dt
   lumino_burn = lumino_burn / dt

   ! Sum them up
   lumino = (lumino_flame + lumino_deton + lumino_burn) 

   END SUBROUTINE findluminosity

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the total neutrino loss
   ! This one is supplemented to GetNuQdot.f90 
   ! Merged by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE FindTotalNeutrinoLoss
   USE definition
   IMPLICIT NONE

   ! Dummy variables
   INTEGER :: j, k

   ! Initilization
   total_nu_qdot = 0.0D0

   ! Do the integration
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1
         total_nu_qdot = total_nu_qdot + vol2(j,k) * nu_qdot (j, k) * rho2 (j, k)
      ENDDO
   ENDDO

   END SUBROUTINE FindTotalNeutrinoLoss

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine outputs the chemical composition
   ! profiles. Additional elements can be insert.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE outputXiso_profile(n)
   USE definition
   IMPLICIT NONE

   ! Input file no
   INTEGER :: n

   !dummy variables
   INTEGER :: i, j, k

   ! File name package
   INTEGER :: fileno_len          
   CHARACTER (len = 256) :: fileno

   ! First assign the file name
   WRITE (fileno, *) n
   fileno = ADJUSTL (fileno)
   fileno_len = LEN_TRIM (fileno)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Output all the elements
   open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XHe4_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,che4), j=1,length_step_r_2)
   enddo         
   write(601,*)
   close(601)

   open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XC12_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cc12), j=1,length_step_r_2)
   enddo
   write(601,*)

      open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XO16_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time    
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,co16), j=1,length_step_r_2)
   enddo     
   write(601,*)
   close(601)

   open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XNe20_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cne20), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XMg24_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt      
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cmg24), j=1,length_step_r_2)
   enddo                          
   write(601,*)
   close(601)

   open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XSi28_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time   
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cSi28), j=1,length_step_r_2)
   enddo          
   write(601,*)
   close(601)

   open(unit=601, FILE = './Outfile/Isotope/Star_WENO_XNi56_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time     
   write(601,*) length_step_r_2, length_step_z_2       
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cni56), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(burn_prog_flag == 1) then

      open(unit=601, FILE = './Outfile/Flame/Star_WENO_BurnPhi1_'//fileno (1 : fileno_len)//'.dat', action='write')
      write(601,*) n, global_time
      write(601,*) length_step_r_2, length_step_z_2
      write(601,*) dx2, dt
      DO k = 1, length_step_z_2  
         WRITE (601, 702) (burn_phi1(j,k), j=1,length_step_r_2)
      enddo
      write(601,*)
      close(601)

      open(unit=601, FILE = './Outfile/Flame/Star_WENO_BurnPhi2_'//fileno (1 : fileno_len)//'.dat', action='write')
      write(601,*) n, global_time
      write(601,*) length_step_r_2, length_step_z_2
      write(601,*) dx2, dt
      DO k = 1, length_step_z_2
         WRITE (601, 702) (burn_phi2(j,k), j=1,length_step_r_2)
      enddo
      write(601,*)
      close(601)

      open(unit=601, FILE = './Outfile/Flame/Star_WENO_BurnPhi3_'//fileno (1 : fileno_len)//'.dat', action='write')
      write(601,*) n, global_time
      write(601,*) length_step_r_2, length_step_z_2
      write(601,*) dx2, dt
      DO k = 1, length_step_z_2    
         WRITE (601, 702) (burn_phi3(j,k), j=1,length_step_r_2)
      enddo
      write(601,*)
      close(601)

      open(unit=601, FILE = './Outfile/Flame/Star_WENO_Yiso_'//fileno (1 : fileno_len)//'.dat', action='write')
      write(601,*) n, global_time
      write(601,*) length_step_r_2, length_step_z_2
      write(601,*) dx2, dt
      DO k = 1, length_step_z_2    
         WRITE (601, 702) (yiso(j,k), j=1,length_step_r_2)
      enddo
      write(601,*)
      close(601)

      open(unit=601, FILE = './Outfile/Flame/Star_WENO_Qash_'//fileno (1 : fileno_len)//'.dat', action='write')
      write(601,*) n, global_time
      write(601,*) length_step_r_2, length_step_z_2
      write(601,*) dx2, dt
      DO k = 1, length_step_z_2    
         WRITE (601, 702) (qash(j,k), j=1,length_step_r_2)
      enddo
      write(601,*)
      close(601)

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(advburn_flag == 1) then

      open(unit=601, FILE = './Outfile/Flame/Star_WENO_BurnStage_'//fileno (1 : fileno_len)//'.dat', action='write')
      write(601,*) n, global_time
      write(601,*) length_step_r_2, length_step_z_2
      write(601,*) dx2, dt                 
      DO k = 1, length_step_z_2
         WRITE (601, 703) (nse_flag(j,k), j=1,length_step_r_2)
      enddo                 
      write(601,*)               
      close(601)

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   702 FORMAT (2000ES18.8)
   703 FORMAT (2000I3)

   end subroutine outputXiso_profile

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine outputs all the elements for 19 isotopes
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung SHing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   subroutine outputXiso_profile19(n)
   use definition
   implicit none

   ! Input file no
   integer :: n

   ! Dummy variables
   integer :: i, j, k

   ! File name package
   INTEGER :: fileno_len
   CHARACTER (len = 256) :: fileno

   ! First assign the file name
   WRITE (fileno, *) n
   fileno = ADJUSTL (fileno)
   fileno_len = LEN_TRIM (fileno)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Output all the isotopes
   open(unit=601, FILE = 'Star_WENO_XHe4_'//fileno (1 : fileno_len)//'.dat', action='write')
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,che4), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XC12_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2 
   write(601,*) dx2, dt
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cc12), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XO16_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,co16), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XNe20_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cne20), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XMg24_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2  
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cmg24), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XSi28_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cSi28), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XS32_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cs32), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XAr36_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,car36), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XCa40_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2  
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cca40), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XTi44_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cti44), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XCr48_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2  
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,ccr48), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XFe52_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2 
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cfe52), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XFe54_'//fileno (1 : fileno_len)//'.dat', action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cfe54), j=1,length_step_r_2)
   enddo
   write(601,*)
   close(601)

   open(unit=601, FILE = 'Star_WENO_XNi56_'//fileno (1 : fileno_len)//'.dat',action='write')   
   write(601,*) n, global_time
   write(601,*) length_step_r_2, length_step_z_2
   write(601,*) dx2, dt  
   DO k = 1, length_step_z_2
      WRITE (601, 702) (xiso(j,k,cni56), j=1,length_step_r_2)
   enddo
   write(601,*) 
   close(601)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   702 FORMAT (2000ES18.8)

   end subroutine OutputXiso_profile19

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds and outputs the total mass of 
   ! all isotopes.
   ! Merged by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE output_xmass()
   USE definition
   IMPLICIT NONE

   ! dummy variables
   INTEGER :: i, j, k, n

   ! Local sum
   REAL (DP), DIMENSION(totalion) :: xmass_sum

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialization
   xmass_sum = 0.0D0
  
   ! Do the integration
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         xmass_sum(:) = xmass_sum(:) + vol2(j,k) * rho2 (j,k) * xiso(j,k,:)

      ENDDO
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Now output the result   
   OPEN(unit=401, FILE='./Outfile/Isotope/Star_WENO_XMass_Selected_0.dat', action='write', position='append')
   WRITE(401,702) global_time, (xmass_sum(i), i = 1, totalion)
   CLOSE(401)

   ! Output other related results
   ! Total neutrino energy loss
   OPEN (UNIT = 401, FILE = './Outfile/Neutrino/Star_WENO_Energy_Loss_0.dat', action='write', position='append')
   WRITE (401, 702) global_time, total_nu_qdot, total_ecap_nu_qdot !/ dt
   close(401)

   ! Central mean atomic mass
   OPEN (UNIT = 401, FILE = './Outfile/Isotope/Star_WENO_CentralABar_0.dat', action='write', position='append')
   WRITE (401, 702) global_time, abar2(1,1)
   close(401)

   ! Central mean atomic number
   OPEN (UNIT = 401, FILE = './Outfile/Isotope/Star_WENO_CentralZBar_0.dat', action='write', position='append')
   WRITE (401, 702) global_time, zbar2(1,1)                      
   CLOSE(401)

   ! Total luminosity and its components
   OPEN (UNIT = 401, FILE = './Outfile/Isotope/Star_WENO_Luminosity_0.dat', action='write', position='append')
   WRITE (401, 702) global_time, lumino, lumino_flame, lumino_deton, lumino_burn, burn_mass
   CLOSE(401)

   IF(burn_prog_flag == 1) THEN
  
      OPEN (UNIT = 401, FILE = './Outfile/Isotope/Star_WENO_CentralYiso_0.dat', action='write', position='append')
      WRITE (401, 702) global_time, yiso(1,1)
      CLOSE(401)

      OPEN (UNIT = 401, FILE = './Outfile/Isotope/Star_WENO_CentralBPhi_0.dat', action='write', position='append')
      WRITE (401, 702) global_time, burn_phi1(1,1), burn_phi2(1,1), burn_phi3(1,1)
      CLOSE(401)

      OPEN (UNIT = 401, FILE = './Outfile/Isotope/Star_WENO_CentralQash_0.dat', action='write', position='append')
      WRITE (401, 702) global_time, qash(1,1)
      CLOSE(401)

   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   702 FORMAT (F23.10, 20ES23.10)
   703 FORMAT (F23.10, 300ES23.10)

   end subroutine output_xmass

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the energy released by 
   ! both deflagration and detonation
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE FLAME()
   USE definition
   USE levelset_module
   USE flametable_module, only : readtable_flameenergy, readtable_detonenergy
   IMPLICIT NONE

   ! Dummy variables
   INTEGER :: j, k

   ! Flag for finding temp
   INTEGER :: flag_notfindtemp

   ! Local variables
   REAL (DP) :: rho_mid, temp_mid, vol_mid

   ! Change of fraction for both level sets
   REAL (DP) :: x1, x2

   ! Local Atwood Number
   REAL (DP) :: at_no

   ! local flame energy
   real (DP) :: flame_ene

   ! local deton energy
   real (DP) :: deton_ene

   ! Ash composition
   real (DP), dimension(totalion) :: x_burn

   ! Other variables
   real (DP) :: esum

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug !
   !if(debug_flag == 1) write(*,*) 'In Flame'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First, initialize the grid, to erase all past data
   flame_qdot(:,:) = 0.0D0
   deton_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !call findenergy
   !write(*,*) 'Before: ', energy2
   !esum = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Then, find the energy release due to nuclear fusion
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

         rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)
         vol_mid = vol2(j,k)

	 ! Just compute the change of fraction occupies
	 ! by deflagration or detonation
         x1 = flame_ratio(j,k) - flame_ratio_old(j,k)
         x2 = deton_ratio(j,k) - deton_ratio_old(j,k)

	 ! Remember to switch flame_ini as well
	 ! When there is a change in area fraction
	 ! Change the chemical composition accordingly
	 ! And inject temperature
         !if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min .and. rho_mid < rho2_flame_max) then
         if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min) then
            call readtable_flameenergy(temp_mid, rho_mid, flame_ene, x_burn(:))
	    flame_qdot(j,k) = flame_ene * x1
            epsilon2(j,k) = epsilon2(j,k) + flame_qdot(j,k) * x1
            xiso(j,k,:) = xiso(j,k,:) + x1 * x_burn(:) - x1 * xiso_a(:)	
	    burn_mass = burn_mass + x1 * rho_mid * vol_mid
         endif
	
	 !Repeat the same procedure for detonation
         !if(x2 > 0.0D0 .and. rho_mid > rho2_deton_min .and. rho_mid < rho2_deton_max) then
	 if(x2 > 0.0D0 .and. rho_mid > rho2_deton_min) then
	    !call readtable_flameenergy(temp_mid, rho_mid, deton_qdot(j,k), x_burn(:))
            call readtable_detonenergy(temp_mid, rho_mid, deton_ene, x_burn(:))
	    deton_qdot(j,k) = deton_ene * x2 
            epsilon2(j,k) = epsilon2(j,k) + deton_qdot(j,k) * x2
            xiso(j,k,:) = xiso(j,k,:) + x2 *  x_burn(:) - x2 * xiso_a(:)
 	    burn_mass = burn_mass + x2 * rho_mid * vol_mid
         endif

      enddo
   enddo

   100 format(10ES15.6)

   end subroutine flame

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the initial input
   ! of energy due to deflagration/detonation
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE flame_ini()
   USE definition
   USE levelset_module
   USE flametable_module, only : readtable_flameenergy, readtable_flamestate
   IMPLICIT NONE

   ! dummy variables
   INTEGER :: j, k

   ! Local variables
   REAL (DP):: epsilon_temp, rho_mid, temp_mid, rho_ash, eps_ash

   ! Local dummy
   REAL (DP) :: dummy
   REAL (DP) :: flame_mid

   ! Ash composition
   REAL (DP), dimension(totalion) :: x_burn

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initilzation
   WRITE(*,*) 'Set initial flame_qdot'
   flame_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Find the energy input
   !WRITE(*,*) 'Give the initial def/det energy'
   !WRITE(*,*) length_step_z_min_part_2, length_step_r_part_2, length_step_z_part_2
   !WRITE(*,"(10ES15.7)") rho2(1,1), temp2(1,1), flame_ratiO(1,1), epsilon2(1,1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

	 ! Store to local variables
         rho_mid = rho2(j,k) 
         temp_mid = temp2(j,k) 
   
	 ! Energy injection and change in chemical composition are done
	 ! to grids which are first assumed to be burnt
         IF(rho_mid > rho2_flame_min .and. flame_ratio(j,k) > 0.0D0) THEN

            CALL readtable_flameenergy(temp_mid, rho_mid, flame_mid, x_burn(:))
	    CALL readtable_flamestate(temp_mid, rho_mid, rho_ash, dummy, eps_ash)

            !epsilon2(j,k) = epsilon2(j,k) + flame_mid * flame_ratio(j,k)
	    epsilon2(j,k) = (1.0D0 - flame_ratio(j,k)) * epsilon2(j,k) + flame_ratio(j,k) * eps_ash
	    rho2(j,k) = (1.0D0 - flame_ratio(j,k)) * rho_mid + flame_ratio(j,k) * rho_ash
            xiso(j,k,:) = (1.0D0 - flame_ratio(j,k)) * xiso(j,k,:) + flame_ratio(j,k) * x_burn(:) 
	    burn_mass = burn_mass + rho_mid * vol2(j,k) !2.0D0 * pi * r2(j) * dx2 * dx2
         endif
  
      ENDDO
   ENDDO

   ! Copy to new results to ghost grids
   CALL BOUNDARY1D_NM(rho2, even)
   CALL BOUNDARY2D_X(xiso)
   CALL BOUNDARY1D_NM(epsilon2, even)

   END SUBROUTINE flame_ini

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the initial input
   ! of energy due to deflagration/detonation
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE flame_run()
   USE definition
   USE levelset_module        
   USE flametable_module, only : readtable_flameenergy, readtable_flamestate
   IMPLICIT NONE

   ! dummy variables
   INTEGER :: j, k

   ! Local variables
   REAL (DP):: epsilon_temp, rho_mid, temp_mid, rho_ash, eps_ash

   ! Local dummy
   REAL (DP) :: dummy
   REAL (DP) :: flame_mid

   ! Ash composition            
   REAL (DP), dimension(totalion) :: x_burn

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initilzation
   WRITE(*,*) 'Set initial flame_qdot'
   flame_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Find the energy input    
   !WRITE(*,*) 'Give the initial def/det energy'
   !WRITE(*,*) length_step_z_min_part_2, length_step_r_part_2, length_step_z_part_2
   !WRITE(*,"(10ES15.7)") rho2(1,1), temp2(1,1), flame_ratiO(1,1), epsilon2(1,1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         ! Store to local variables
         rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)

         ! Energy injection and change in chemical composition are done
         ! to grids which are first assumed to be burnt
         IF(rho_mid > rho2_flame_min .and. flame_ratio(j,k) > 0.0D0) THEN

            CALL readtable_flameenergy(temp_mid, rho_mid, flame_mid, x_burn(:))
            CALL readtable_flamestate(temp_mid, rho_mid, rho_ash, dummy, eps_ash)

            !epsilon2(j,k) = epsilon2(j,k) + flame_mid * flame_ratio(j,k)
            epsilon2(j,k) = (1.0D0 - flame_ratio(j,k)) * epsilon2(j,k) + flame_ratio(j,k) * eps_ash
            xiso(j,k,:) = (1.0D0 - flame_ratio(j,k)) * xiso(j,k,:) + flame_ratio(j,k) * x_burn(:)
            burn_mass = burn_mass + rho_mid * vol2(j,k) !2.0D0 * pi * r2(j) * dx2 * dx2
         endif                  

      ENDDO         
   ENDDO

   ! Copy to new results to ghost grids      
   CALL BOUNDARY2D_X(xiso)
   CALL BOUNDARY1D_NM(epsilon2, even)

   END SUBROUTINE flame_run

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the initial energy
   ! input and composition change due to detonation
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine DETON_ini()
   use definition
   use levelset_module
   use flametable_module
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Local variables
   real (DP):: rho_mid, temp_mid

   ! Ash composition
   real (DP), dimension(totalion) :: x_burn

   ! Local variables
   REAL (DP):: epsilon_temp, rho_ash, eps_ash

   ! Local dummy
   REAL (DP) :: dummy
   REAL (DP) :: flame_mid       

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initilziation
   write(*,*) 'Set initial flame_qdot'
   deton_qdot(:,:) = 0.0D0 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

   ! Find the energy input
   write(*,*) 'Read the initial flametable'            
   do k = length_step_z_min_part_2, length_step_z_part_2, 1         
      do j = 1, length_step_r_part_2, 1

	 ! Store the data to local variables
         rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)

	 ! Energy injection and change in chemical composition are done
         ! to grids which are first assumed to be burnt
         if(rho_mid > 1.62D-13 .and. deton_ratio(j,k) > 0.0D0) then

	    CALL readtable_flameenergy(temp_mid, rho_mid, flame_mid, x_burn(:))
            CALL readtable_flamestate(temp_mid, rho_mid, rho_ash, dummy, eps_ash)
	    epsilon2(j,k) = (1.0D0 - deton_ratio(j,k)) * epsilon2(j,k) + deton_ratio(j,k) * eps_ash
            rho2(j,k) = (1.0D0 - deton_ratio(j,k)) * rho_mid + deton_ratio(j,k) * rho_ash
            xiso(j,k,:) = (1.0D0 - deton_ratio(j,k)) * xiso(j,k,:) + deton_ratio(j,k) * x_burn(:)
            burn_mass = burn_mass + rho_mid * vol2(j,k)

            !call readtable_detonenergy(temp_mid, rho_mid, deton_qdot(j,k), x_burn(:))
	    !x_burn(:) = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0/)
            !epsilon2(j,k) = epsilon2(j,k) + deton_qdot(j,k) * deton_ratio(j,k)
            !xiso(j,k,:) = (1.0D0 - deton_ratio(j,k)) * xiso(j,k,:) + deton_ratio(j,k) * x_burn(:)
         endif

      enddo
   enddo

   ! Copy the result to the ghost grids
   CALL BOUNDARY2D_X(xiso)
   CALL BOUNDARY1D_NM(epsilon2, even)

   end subroutine DETON_INI

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the initial energy
   ! input and composition change due to nuclear reaction
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine BURN()
   use definition
   use levelset_module
   implicit none

   ! Dummy variables 
   integer :: j, k

   ! Local variables
   real (DP) :: temp_mid, rho_mid

   ! Ash composition
   real (DP), dimension(totalion) :: x_burn

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

   ! Initialzation
   burn_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

   if(flame_flag == 1) then

      ! When there is flame, only completely
      ! burnt grids are allowed for nuclear reaction
   
      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1

	    ! First store the data in local variables
            rho_mid = rho2(j,k)
            temp_mid = temp2(j,k)

	    ! Do all the check
            if(temp_mid < 2.0D0) cycle 
            if(rho_mid > rho2_burn_max) cycle 
            if(rho_mid < rho2_burn_min) cycle 
            if(burn_ratio(j,k) /= 1.0D0) cycle
            if(xiso(j,k,cc12) < 1.0D-2) cycle 
            if(xiso(j,k,co16) < 1.0D-2) cycle

	    ! Call subroutine for the energy input
            !call drive_aprox19(dt, temp_mid, rho_mid, xiso(j,k,:), burn_qdot(j,k), x_burn(:))  !For old one
            !call drive_aprox13(dt, temp_mid, rho_mid, xiso(j,k,:), burn_qdot(j,k), x_burn(:))
            epsilon2(j,k) = epsilon2(j,k) + burn_qdot(j,k)
            xiso(j,k,:) = x_burn(:)

         enddo
      enddo

   else

      ! For usual case, all grids with sufficient
      ! temperature and density are burnt

      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1

	    ! Store the data in local variables
            rho_mid = rho2(j,k)
            temp_mid = temp2(j,k)

	    ! Do the density and temperature check
            if(temp_mid < 2.0D0) cycle
            if(rho_mid > rho2_burn_max) cycle
            if(rho_mid < rho2_burn_min) cycle
            if(xiso(j,k,cc12) < 1.0D-2) cycle
            if(xiso(j,k,co16) < 1.0D-2) cycle

	    ! Call subroutine for the energy input
            !call drive_aprox19(dt, temp_mid, rho_mid, xiso(j,k,:), burn_qdot(j,k), x_burn(:))  !For old one
            !call drive_aprox13(dt, temp_mid, rho_mid, xiso(j,k,:), burn_qdot(j,k), x_burn(:))
            epsilon2(j,k) = epsilon2(j,k) + burn_qdot(j,k)
            xiso(j,k,:) = x_burn(:)

         enddo
      enddo

   endif

   ! Copy the result to ghost cells

   end subroutine burn

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the chemical composition
   ! due to nuclear statistical equilibrium (NSE)    
   ! Written by Leung Shing Chi in 2016
   ! Detailed procedure refer to Reinecke (2002a)
   ! and Timme's website
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine NSE()
   use definition
   use levelset_module
   use ecaptable_module   
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Local variables
   real (DP):: temp_mid, rho_mid, ecaprate_mid

   ! Ash composition
   real (DP), dimension(totalion) :: x_burn

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug stuff !
   !if(debug_flag == 1) write(*,*) 'In NSE'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

   if(flame_flag == 1) then

      ! When there is flame, then only completely
      ! burnt grids can have NSE

      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1

	    ! Store in local variables
            rho_mid = rho2(j,k)
            temp_mid = temp2(j,k)
         
	    ! Do the check
            if(temp_mid < 5.0D0) cycle 
            if(burn_ratio(j,k) /= 1.0D0) cycle
      
	    ! Get the composition
            !!call nse_interface(temp_mid ,rho_mid, zbar2(j,k)/abar2(j,k), x_burn)  
	    !call getecaprate(rho_mid, temp_mid, ye2(j,k), ecaprate_mid)
            !ye2(j,k) = ye2(j,k) + ecaprate_mid * dt
	    call getnsestate(rho_mid, temp_mid, x_burn)

	    ! Replace them 
            xiso(j,k,:) = x_burn(:)
            nse_flag(j,k) = 2
   
         enddo
      enddo

   else

      ! Othewise, all grids with high enough temperature
      ! are assumed to be in NSE

      do k = length_step_z_min_part_2, length_step_z_part_2, 1
         do j = 1, length_step_r_part_2, 1

	    ! Stored in local variables
            rho_mid = rho2(j,k)              
            temp_mid = temp2(j,k)            
        
	    ! Do the check
            if(temp_mid < 5.0D0) cycle

	    ! Get the composition
            !!call nse_interface(temp_mid ,rho_mid, zbar2(j,k)/abar2(j,k), x_burn)
	    !call getecaprate(rho_mid, temp_mid, ye2(j,k), ecaprate_mid)
            !ye2(j,k) = ye2(j,k) + ecaprate_mid * dt
	    call getnsestate(rho_mid, temp_mid, x_burn)

	    ! Replace them
            xiso(j,k,:) = x_burn(:)
            nse_flag(j,k) = 2
   
         enddo
      enddo

   endif

   ! Copy the data to ghost cells
   call boundary2D_X()
   call boundary1D_NM(ye2)
   
   end subroutine nse

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the chemical composition
   ! due to nuclear statistical equilibrium (NSE)
   ! Written by Leung Shing Chi in 2016, but 
   ! assuming the binding energy change due to 
   ! the composition difference causes gain/loss in 
   ! local energy density
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine NSE2
   use definition                               
   use levelset_module
   use ecaptable_module
   implicit none                               

   ! Flag for finding temperature
   integer :: flag_notfindtemp             

   ! Dummy variables
   integer :: i, j, k, k2

   ! Dummy variables
   real (DP) :: ye_sample = 0.5D0, dummy

   ! Local variables
   real (DP) :: abar_mid, zbar_mid, ye_mid

   ! Input variables
   real (DP) :: temp_beg    
   real (DP) :: eps_beg

   ! Trial variables
   real (DP) :: temp_mid                        
   real (DP) :: eps_mid         
   real (DP) :: rho_mid

   ! Change in temperature
   real (DP) :: dtemp_mid 

   ! Number of successful trial                
   integer  :: count_digit

   ! Ecap rate and neutrino energy loss
   real (DP) :: ecaprate_mid
   real (DP) :: eneurate_mid

   ! Initial and Expected chemical composition
   real (DP), dimension(totalion) :: x_mid
   real (DP), dimension(totalion) :: x_burn

   ! Binding energy
   real (DP) :: binde_bf, binde_af, deps_nuc

   ! Energy balance equation
   real (DP) :: check_e, check_e_last

   ! Timescale for burning to NSE
   real (DP) :: nse_burntime, temp_nse, fnse

   ! Other !
   REAL (DP) :: generate, eps_old

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   total_ecap_nu_qdot = 0.0D0
   generate = 0.0D0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !if(debug_flag == 1) write(*,*) 'In NSE2'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialization
   ! Already done in burn_phase2b
   !burn_qdot(:,:) = 0.0D0
   nu_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k = length_step_z_min_part_2, length_step_z_part_2, 1              
      do j = 1, length_step_r_part_2, 1                    

         ! Do the checking first               
         if(temp2(j,k) < 5.0D0) cycle
	 if(temp2(j,k) > 11.00) cycle
	 if(ye2(j,k) < 0.40D0) cycle
	 if(ye2(j,k) > 0.50D0) cycle
	 if(rho2(j,k) < 1.62D-11) cycle
         if(burn_ratio(j,k) /= 1.0D0) cycle

	 ! Special test if multi-stage burning is used
	 !if(advburn_flag == 1) then
	 if(nse_flag(j,k) == 1) cycle
	 !endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! If they pass then start the iteration
         nse_flag(j,k) = 2                

	 ! The density is not changed, so 
	 ! the trial value is also the 
	 ! final value
         rho_mid = rho2(j,k)

         ! Give a trial temperature
         temp_beg = temp2(j,k)                  
         eps_beg = epsilon2(j,k)

	 ! Also give the trial composition
         abar_mid = abar2(j,k)
         zbar_mid = zbar2(j,k)
         ye_mid = ye2(j,k)                
	 x_mid(:) = xiso(j,k,:)
 	 eps_old = epsilon2(j,k)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !if(flame_loc_ratio(j,k) == 1.0D0 .or. deton_loc_ratio(j,k) == 1.0D0) then
	 if(burn_ratio(j,k) == 1.0D0) then

            ! Scheme for electron capture
	    ! Get the Ecap rate and energy loss
            call getecaprate(rho_mid, temp_beg, ye2(j,k), ecaprate_mid, eneurate_mid)
	    !if(ye_mid < 0.46D0) then
	    !   ecaprate_mid = 0.0D0; eneurate_mid = 0.0D0
	    !endif


	    ! Update the Ye rate
	    ! Note: No iteration is done here because
	    ! the temperature sensitivity of Ecap
	    ! rate is much smaller than NSE composition
            ye_mid = ye2(j,k) + ecaprate_mid * dt
            !ye_mid = ye2(j,k)                 

	    ! Compute the binding energy of the 
	    ! initial composition
            call compute_binde(x_mid, binde_bf)

            !write(*,*) 'Temp Before: ', temp_beg
            !write(*,*) 'Abar before: ', abar2(j,k)
            !write(*,*) 'Zbar before: ', zbar2(j,k)         
            !write(*,*) 'Eps before: ', eps_beg
            !write(*,*) 'Binde before: ', binde_bf
            !write(*,*)                 

	    ! Prepare for the search of temperature
            count_digit = 0
            temp_mid = temp_beg
            dtemp_mid = 0.01D0 * temp_beg          

	    ! Some debug stuff
	    !if(j == 1 .and. k == 1) write(*,*) 'Now in grid', j, k
	    !if(j == 1 .and. k == 1) write(*,*) 'Input density = ', rho_mid
	    !if(j == 1 .and. k == 1) write(*,*) 'Old temp = ', temp_beg
	    !if(j == 1 .and. k == 1) write(*,101) (xiso(j,k,k2),k2=1,totalion)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do i = 0, 400, 1   

	       ! Now do the search

               !call nse_interface(temp_mid ,rho_mid, ye_mid, x_burn)

	       ! Get the trial NSE state by the trial temp
	       call getnsestate(rho_mid, temp_mid, x_burn)

	       ! Compute the trial binding energy
               call compute_binde(x_burn, binde_af)

	       ! Calculate the trial abar and zbar
               call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)

	       ! Get the trial epsilon
               !call teos_epsilon(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid)
	       CALL HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)	       

	       ! Calculate the binding energy change
               deps_nuc = binde_af - binde_bf

	       ! Check if the energy is balanced
               check_e = eps_mid - eps_beg - deps_nuc - &  !Original +
			 8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt

	       ! Make sure you go to the right direction of dtemp
               if(i == 0) then
                  if(check_e > 0.0D0) dtemp_mid = -dtemp_mid
               endif

               ! Use shooting method    
               if(check_e_last * check_e < 0.0D0 .and. i /= 0) then
                  temp_mid = temp_mid - dtemp_mid
                  dtemp_mid = dtemp_mid * 0.1D0
                  temp_mid = temp_mid + dtemp_mid  
                  count_digit = count_digit + 1
               else
                  temp_mid = temp_mid + dtemp_mid
                  check_e_last = check_e
               endif

               if(count_digit == 4) exit
               if(i == 400 .or. temp_mid < 5.0D0) then
                  !stop 'Check NSE solver'
                  temp_mid = temp_beg        
                  eps_mid = eps_beg                  
                  !call nse_interface(temp_mid ,rho_mid, ye_mid, x_burn)
	          call getnsestate(rho_mid, temp_mid, x_burn)
	 	  exit
               endif            

	       !if(j == 1 .and. k == 1) write(*,100) i, temp_mid, check_e, ecaprate_mid, eneurate_mid, (x_burn(k2), k2 = 1, totalion)
	       !if(mod(i,100) == 0 .and. j == 1 .and. k == 1) read(*,*)

            enddo

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    !nse_burntime = EXP(196.02D0 / temp_mid - 41.646D0) / 4.9282D-6
	    temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
            nse_burntime = EXP(196.02D0 / temp_nse - 41.646D0) / 4.9282D-6

	    if(dt > nse_burntime) then

  	       ! When things are done, move the refined
	       ! trial results to the outpuit
               temp2(j,k) = temp_mid       
               epsilon2(j,k) = eps_mid
               ye2(j,k) = ye_mid                
               xiso(j,k,:) = x_burn(:) 
               burn_qdot(j,k) = burn_qdot(j,k) + binde_af - binde_bf
	       !nu_qdot(j,k) = eneurate_mid 

	    else

	       ! When things are partially done, use
	       ! linear interpolation
		fnse =  -1.0604D0 * (dt / nse_burntime)**2 + 2.0604D0 * (dt / nse_burntime)**2
	       temp_mid = temp_beg + (temp_mid - temp_beg) * fnse !dt / nse_burntime
	       x_burn(:) = x_mid(:) + (x_burn(:) - x_mid(:)) * fnse !dt / nse_burntime
	       call compute_binde(x_burn, binde_af)
	       call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)	
	       !call teos_epsilon(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid)	   
	       CALL HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)

	       ! Now patch the result
	       temp2(j,k) = temp_mid
               epsilon2(j,k) = eps_mid
               ye2(j,k) = ye_mid    
               xiso(j,k,:) = x_burn(:)
               burn_qdot(j,k) = burn_qdot(j,k) + binde_af - binde_bf
	       !nu_qdot(j,k) = eneurate_mid 

	    endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         else        

	    ! For detonation, no NSE energy is accounted
	    temp_mid = temp2(j,k)

	    ! Electron capture for deflagration zone
	    !if(flame_loc_ratio(j,k) == 1.0D0) then
	    !   call getecaprate(rho_mid, temp_mid, ye2(j,k), ecaprate_mid, eneurate_mid)
            !   ye2(j,k) = ye2(j,k) + ecaprate_mid * dt
	    !endif

	    ! Get the NSE state
	    call getnsestate(rho_mid, temp_mid, x_burn)

	    ! Copy the result to output
            xiso(j,k,:) = x_burn(:)
            burn_qdot(j,k) = 0.0D0  

         endif
	 generate = generate + rho2(j,k)*(epsilon2(j,k) - eps_old)*vol2(j,k)
 	 total_ecap_nu_qdot = total_ecap_nu_qdot + vol2(j,k) * rho_mid * (eneurate_mid * dt)
      enddo         
   enddo
   energy_input = energy_input + generate

   ! Copy the results to ghost cells
   call boundary1D_NM(temp2, even)                 
   call boundary1D_NM(epsilon2, even)
   call boundary1D_NM(ye2, even)
   call boundary2D_X()   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              

   100 format(I5, 20ES15.7)
   101 format(10ES15.7)

   end subroutine nse2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine compute the binding energy for 
   ! a given composition. The subroutine takes 
   ! chemical composition as input and return the
   ! binding energy.
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_binde(x_in, binde_out)
   use definition, only : DP
   implicit none

   ! Input chemical composition
   real (DP) :: x_in(1:totalion)

   ! Output binding energy 
   real (DP) :: binde_out

   ! Dummy variable
   integer :: i

   ! Binding energy
   real (DP) :: binde_sum

   ! Initialiation
   binde_sum = 0.0D0
  
   ! Sum all the isotopes
   do i = 1, totalion, 1
      binde_sum = binde_sum + x_in(i) / mion(i) * bion(i)
   enddo

   ! Change to code unit
   binde_out = 1.78D-27 * binde_sum

   end subroutine compute_binde

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the mean atomic mass and mean atomic
   ! number of all grids. 
   ! Written by Leung Shing Chi
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This subroutines finds abar and zbar in all grid

   subroutine find_AZbar()
   use definition
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Dummy variable
   real (DP) :: dummy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !if(debug_flag == 1) write(*,*) 'In FindAZBar'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Simply call the subroutine to help you
   ! find the abar and zbar directly
   do k = 1, length_step_z_2, 1
      do j = 1, length_step_r_2, 1
         CALL PRIVATE_HELMEOS_AZBAR(xiso(j,k,:), abar2(j,k), zbar2(j,k), dummy)
      enddo
   enddo

   ! Copy the result to ghost cells
   CALL BOUNDARY1D_NMFULL(abar2, even)
   CALL BOUNDARY1D_NMFULL(zbar2, even)
   call boundary1D_NMFULL(ye2, even)

   end subroutine find_AZbar

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine bridges the NUCLEAR_NET subroutine which 
   ! aims at finding the abar and zbar.
   ! Written by Leung Shin Chi in 2016.
   ! For more details about nuclear net, refer Timmes (2000b)
   ! or his website.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine private_helmeos_azbar(xiso_in, abar_out, zbar_out, electron_frac_out)
   implicit none

   !For 7 isotopes only, expand the function when more isotopes are included
   !The sequence is as follows: He4, C12, O16, Ne20, Mg24, Si28, Ni56

   ! Input quantities
   real (selected_real_kind(15,307)), dimension(totalion):: xiso_in

   ! Output quantitites
   real (selected_real_kind(15,307)) :: abar_out, zbar_out, electron_frac_out

   ! Dummy variables
   real (selected_real_kind(15,307)) :: dummy
   real (selected_real_kind(15,307)), dimension(totalion) :: dummy_x

   ! The traditional way
   !ymass(1:ionmax) = xmass(1:ionmax)/wion(1:ionmax)      
   !wbar  = 1.0d0/sum(ymass(1:ionmax))
   !sum1  = sum(aion(1:ionmax)*ymass(1:ionmax))
   !abar2  = wbar * sum1
   !electron_frac2  = sum(zion(1:ionmax)*ymass(1:ionmax))
   !zbar2  = wbar * ye

   ! Call the preset subroutine
   call azbar(xiso_in,aion,zion,wion,totalion,dummy_x,abar_out,zbar_out,dummy,electron_frac_out,dummy)

   end subroutine private_helmeos_azbar

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine checks if the local isotope mass fraction sum
   ! to one. If not, normalize it. 
   ! Written by Leung Shing Chi in 2016.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine checkxisotope()
   use definition
   implicit none

   ! dummy variables
   integer :: i, j, k
 
   ! local variables
   real (DP):: xiso_sum

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Debug !
   !if(debug_flag == 1) write(*,*) 'In CheckXIsotope'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First check if any grid has any unusual mass fraction
   do k = 1, totalion, 1

      do j = length_step_z_min_part_2, length_step_z_part_2, 1
         do i = 1, length_step_r_part_2, 1
            if(xiso(i,j,k) < 1.0D-30) xiso(i,j,k) = 1.0D-30
         enddo
      enddo

   enddo

   ! Then do the sum and normalization
   do j = length_step_z_min_part_2, length_step_z_part_2, 1
      do i = 1, length_step_r_part_2, 1

         xiso_sum = 0.0D0
         do k = 1, totalion, 1
            xiso_sum = xiso_sum + xiso(i,j,k)
         enddo
         xiso(i,j,:) = xiso(i,j,:) / xiso_sum

      enddo
   enddo

   CALL BOUNDARY2D_X()

   end subroutine checkxisotope

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           
   !
   ! The subroutine bridges the GetNuQdot and finds the energy
   ! loss rate by neutrino of the hydro grids. 
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine findneutrinoloss()
   use definition
   implicit none

   ! dummy variables
   integer :: j, k

   ! local variables
   real (DP) :: temp_in, rho_in, abar_in, zbar_in   

   ! initialization
   nu_qdot = 0.0D0

   DO k = length_step_z_min_part_2, length_step_z_part_2
      DO j = 1, length_step_r_part_2

	 ! Pass the data to local storage
         rho_in = rho2(j,k) * 6.144D17
         temp_in = temp2(j,k) * 1.0D9
         abar_in = abar2(j,k)
         zbar_in = zbar2(j,k)

	 ! Call GetNuQdot if the matter is sufficiently hot and dense
         IF(temp2(j,k)>0.5D0 .and. rho2(j,k) > rho2_a) THEN
            call private_sneut5(temp_in, rho_in, abar_in, zbar_in, nu_qdot(j,k))
            nu_qdot(j,k) = nu_qdot(j,k) / 1.8262D26
         ENDIF

      ENDDO
   ENDDO

   ! Copy the result to ghost grids
   CALL BOUNDARY1D_NM(nu_qdot, even)

   end subroutine findneutrinoloss

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine advect the isotope composition by PPM.
   ! For PPM, refer Collela (1984). 
   ! Written by Leung Shing Chi in 2016
   !
   ! Note: In this version, all PPM is muted to ensure accuracy	
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !subroutine updatexisotopePPM()
   !use definition
   !implicit none
   !
   !! Dummy variables
   !integer :: j, k, k2
   !
   ! Local variables for spatial derivative
   !real (DP), dimension(-4:length_step_r_2+5, -4:length_step_z_2+5) :: dxdr, dxdz
   !
   !do k2 = 1, totalion, 1
   !
   !   call splitPPM(xiso(:,:,k2), dxdr, 1, 1)
   !   call splitPPM(xiso(:,:,k2), dxdz, 2, 2)
   !
   !   do k = 1, length_step_z_part_2, 1
   !      do j = 1, length_step_r_part_2, 1
   !         delta_xiso(j,k,k2) = ((vel2_r(j,k) - vel_frame_r(j,k)) * dxdr(j,k) + &
   !                               (vel2_z(j,k) - vel_frame_z(j,k)) * dxdz(j,k))
   !      enddo
   !   enddo
   !
   !enddo
   !
   !end subroutine updatexisotopePPM
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                             
   !
   ! This subroutine is obtained from Timmes' nuclear network
   ! program to calculate the abar and zbar for a given 
   ! chemical composition in a flexible structure
   ! Merged by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine azbar(xmass,aion,zion,wion,ionmax, &
                    ymass,abar,zbar,wbar,ye,nxcess)
   use definition, only : DP
   implicit none

   ! this routine calculates composition variables
   ! input:
   ! mass fractions               = xmass(1:ionmax)  dimensionless
   ! number of nucleons           = aion(1:ionmax)   dimensionless
   ! charge of nucleus            = zion(1:ionmax)   dimensionless
   ! atomic weight or molar mass  = wion(1:ionmax)    g/mole
   ! number of isotopes           = ionmax
   !
   ! output:
   ! molar abundances        = ymass(1:ionmax)   mole/g
   ! mean number of nucleons = abar              dimensionless
   ! mean nucleon charge     = zbar              dimensionless
   ! mean weight             = wbar              g/mole
   ! electron fraction       = ye                mole/g
   ! neutron excess          = xcess


  ! declare the pass
   integer          ionmax
   real (selected_real_kind(15,307)), dimension(1:totalion) :: xmass,aion,zion,wion,ymass
   real (selected_real_kind(15,307)) :: abar,zbar,wbar,ye,nxcess


   ! local variables
   real (selected_real_kind(15,307)) asum,sum1

   ! molar abundances
   ymass(1:totalion) = xmass(1:totalion)/wion(1:totalion)

   ! mean molar mass
   wbar  = 1.0d0/sum(ymass(1:totalion))

   ! mean number of nucleons
   sum1  = sum(aion(1:totalion)*ymass(1:totalion))
   abar  = wbar * sum1

   ! mean charge
   ye  = sum(zion(1:totalion)*ymass(1:totalion))
   zbar  = wbar * ye

   ! neutron excess
   nxcess = sum1 - 2.0d0 * ye

   return
   end subroutine

   !---------------------------------------------------------------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine find the density when given a pressure 
   ! and temperature, using the subroutine from nuclear 
   ! network developed by Timmes (See Timmes, Apj 1999)
   ! Merged by Leung Shing Chi
   ! Not needed in the new version 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine private_invert_helm_pt_quiet
   include 'implno.dek'
   include 'const.dek'
   include 'vector_eos.dek'


   ! given the pressure, temperature, and composition
   ! find everything else

   ! it is assumed that ptot_row(j), temp_row(j), abar_row(j),   
   ! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
   ! been set before calling this routine.

   ! on input den_row(j) conatins a guess for the density,   
   ! on output den_row(j) contains the converged density.

   ! To get the greatest speed advantage, the eos should be fed a
   ! large pipe of data to work on.

   ! this version is quiet on all errors

   ! local variables
   integer          i,j,jlo_save,jhi_save
   double precision den,f,df,dennew,eostol,fpmin
   parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)

   ! initialize
   jlo_save = jlo_eos
   jhi_save = jhi_eos
   do j=jlo_eos, jhi_eos
      eoswrk01(j) = 0.0d0
      eoswrk02(j) = 0.0d0
      eoswrk03(j) = ptot_row(j)
      eoswrk04(j) = den_row(j)
   end do

   ! do the first newton loop with all elements in the pipe
   call helmeos

   do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

   ! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

   ! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

   ! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

   enddo

   ! now loop over each element of the pipe individually
   do j = jlo_save, jhi_save

       do i=2,40

          if (eoswrk01(j) .lt. eostol .or. &
              abs(eoswrk02(j)) .le. fpmin) goto 20

          jlo_eos = j
          jhi_eos = j

          call helmeos

          f     = ptot_row(j)/eoswrk03(j) - 1.0d0
          df    = dpd_row(j)/eoswrk03(j)
          eoswrk02(j) = f/df

          ! limit excursions to factor of two changes
          den    = den_row(j)
          dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

          ! compute the error
          eoswrk01(j)  = abs((dennew - den)/den)

          ! store the new density, keep it within the table limits
          den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

	  ! end of netwon loop
      end do

   ! we did not converge if we land here
   !      write(6,*)
   !      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
   !      write(6,*) 'pipeline element',j
   !      write(6,01) 'pres  =',eoswrk03(j)
   ! 01   format(1x,5(a,1pe16.8))
   !      write(6,01) 'error =',eoswrk01(j),
   !     1            '  eostol=',eostol,'  fpmin =',fpmin
   !      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
   !      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
   !      write(6,*)
   !      stop 'could not find a density in routine invert_helm_pt'

   ! land here if newton loop converged, back for another pipe element
   20    continue
   end do

   ! call eos one more time with the converged value of the density

   jlo_eos = jlo_save
   jhi_eos = jhi_save

   call helmeos

   return
   end subroutine
   !---------------------------------------------------------------------

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine reconstruct the chemical composition
   ! based on the burning-progress variables, for detilas
   ! see Townsley et al, ApJ 2007 and Calder et al, ApJ 2007
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE reconstructXiso
   USE definition
   IMPLICIT NONE

   INTEGER :: i, j, k

   ! Local variables
   REAL (DP) :: burn_phi1_mid, burn_phi2_mid, burn_phi3_mid, yiso_mid

   ! Dummy variables
   REAL (DP) :: xiso_sum, yiso_sum, yiso_nse

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

	 burn_phi1_mid = burn_phi1(j,k)
	 burn_phi2_mid = burn_phi2(j,k)
	 burn_phi3_mid = burn_phi3(j,k)
	 yiso_mid = yiso(j,k)

	 ! For CO burning
         !xiso(j,k,cc12) = (1.0D0 - burn_phi1_mid) * xiso_a(cc12)
	 !xiso(j,k,co16) = (1.0D0 - burn_phi2_mid) * (1.0D0 - xiso_a(cc12))
	 !xiso(j,k,cne20) = 0.0D0
	 !xiso(j,k,cmg24) = (burn_phi1_mid - burn_phi2_mid) * xiso_a(cc12)
	 !xiso(j,k,csi28) = (burn_phi2_mid - burn_phi3_mid)

	 ! For ONe burning
	 xiso(j,k,cc12) = 0.0D0
         xiso(j,k,co16) = (1.0D0 - burn_phi2_mid) * (1.0D0 - xiso_a(cne20))
         xiso(j,k,cne20) = (1.0D0 - burn_phi1_mid) * xiso_a(cne20)
         xiso(j,k,cmg24) = (burn_phi1_mid - burn_phi2_mid) * xiso_a(cne20)
         xiso(j,k,csi28) = (burn_phi2_mid - burn_phi3_mid)

	 ! Solve for X_He and X_Ni based on yiso 
	 xiso_sum = xiso(j,k,cc12) + xiso(j,k,co16) + xiso(j,k,cne20) + xiso(j,k,cmg24) + xiso(j,k,csi28)
         yiso_sum = xiso(j,k,cc12) / aion(cc12) + xiso(j,k,co16) / aion(co16) + &
		    xiso(j,k,cne20) / aion(cne20) + xiso(j,k,cmg24) / aion(cmg24) + xiso(j,k,csi28) / aion(csi28)
	 yiso_nse = yiso_mid - yiso_sum
         xiso(j,k,che4) = (yiso_nse - (1.0D0 - xiso_sum) / aion(cni56)) / (1.0D0 / aion(che4) - 1.0D0 / aion(cni56))
	 xiso(j,k,cni56) = 1.0D0 - xiso_sum - xiso(j,k,che4)

	 !if(j == 1 .and. k == 1) then 
         !   write(*,*) burn_phi1_mid, burn_phi2_mid, burn_phi3_mid
	 !   write(*,*) xiso(j,k,:)
	 !   read(*,*)
	 !endif

      ENDDO
   ENDDO

   CALL boundary2D_X

   END SUBROUTINE reconstructXIso

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine do the first phase of burning
   ! from 12C or 20Ne to Mg24 using the burning varaible algorithm
   !
   ! 
   !
   ! Written by Leung Shing Chi in 2-16
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   SUBROUTINE burn_phase1
   use definition
   use levelset_module
   implicit none

   integer :: i, j, k

   ! local variables
   REAL (DP) :: rho_mid

   ! Change of fraction for both level sets
   REAL (DP) :: x1, x2

   REAL (DP) :: def_energy, det_energy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First initialize flame_qdot and deton_qdot
   flame_qdot(:,:) = 0.0D0
   deton_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !IF(debug_flag == 1) write(*,*) 'In Burn Phase 1'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

	 rho_mid = rho2(j,k) 

	 ! Calculate the local change of flame_ratio and deton_ratio
         x1 = MIN(MAX(flame_ratio(j,k) - flame_ratio_old(j,k), 0.0D0), xiso(j,k,cne20) / xiso_a(cne20))
      	 x2 = MIN(MAX(deton_ratio(j,k) - deton_ratio_old(j,k), 0.0D0), xiso(j,k,cne20) / xiso_a(cne20))

	 ! If there is change in flame_ratio, add in def. energy
         IF(x1 > 0.0D0) THEN

	    ! For CO matter
            !def_energy = 6.2768D-4 * xiso_a(cc12) * x1

            ! For ONe matter
            def_energy = 2.4718D-4 * xiso_a(cne20) * x1

	    flame_qdot(j,k) = def_energy
	    epsilon2(j,k) = epsilon2(j,k) + def_energy

	    burn_phi1(j,k) = burn_ratio(j,k)
	    burn_mass = burn_mass + x1 * rho_mid * vol2(j,k)

	    ! For CO matter

	    ! For ONe matter
            yiso(j,k) = yiso(j,k) - x1 * xiso_a(cne20) / aion(cne20) + x1 * xiso_a(cne20) / aion(cmg24)

 	 ENDIF

	 ! If there is change in deton ratio, add in det. energy
	 IF(x2 > 0.0D0) THEN
	
	    ! For CO matter
            !def_energy = 6.2768D-4 * xiso_a(cc12) * x2

            ! For ONe matter
            def_energy = 2.4718D-4 * xiso_a(cne20) * x2

            deton_qdot(j,k) = det_energy          
            epsilon2(j,k) = epsilon2(j,k) + det_energy       

	    burn_phi1(j,k) = burn_ratio(j,k)
	    burn_mass = burn_mass + x2 * rho_mid * vol2(j,k)

 	    ! For ONe matter
            yiso(j,k) = yiso(j,k) - x2 * xiso_a(cne20) / aion(cne20) + x2 * xiso_a(cne20) / aion(cmg24)

	 ENDIF

      ENDDO
   ENDDO

   ! Copy the results to ghost cells
   CALL boundary1D_NM(epsilon2, even)
   CALL boundary1D_NM(burn_phi1, even)

   end subroutine burn_phase1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculate the energy released by 
   ! the 2nd phase of nuclear burning, from 16O and 24Mg to 28Si
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE burn_phase2
   USE definition
   USE levelset_module
   IMPLICIT NONE

   INTEGER :: j, k

   ! Local burn_phi           
   REAL (DP) :: burn_phi1_new 
   REAL (DP) :: burn_phi2_ini
   REAL (DP) :: burn_phi2_new

   REAL (DP) :: dburn_phi2

   ! Local hydro variables
   REAL (DP) :: rho_mid, temp_mid

   ! NQSE burning time
   REAL (DP) :: nqse_burntime, nqse_burntime_rho, nqse_burntime_temp
   REAL (DP) :: temp_nse
   REAL (DP) :: eps_ini

   ! Local flame energy   
   REAL (DP) :: burn_ene1 = 5.1004D-4   ! Energy release form 16O  --> 28Si
   REAL (DP) :: burn_ene2 = 2.0237D-4   ! Energy release from 24Mg --> 28Si

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1  

         ! Pass to local variables
         eps_ini = epsilon2(j,k)              
         rho_mid = rho2(j,k)            
         temp_mid = temp2(j,k)    

	 ! Pass the burning variables
         burn_phi1_new = burn_phi1(j,k)
         burn_phi2_ini = burn_phi2(j,k)   

	 IF(burn_ratio(j,k) == 1.0D0) THEN
        
	    nqse_burntime_temp = EXP(182.06D0 / temp_mid - 46.054D0) / 4.9282D-6
	    temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
            nqse_burntime_rho = EXP(182.06D0 / temp_nse - 46.054D0) / 4.9282D-6
	    nqse_burntime = MIN(nqse_burntime_temp, nqse_burntime_rho)

            burn_phi2_new = burn_phi1_new * (1.0D0 - EXP(-dt / nqse_burntime)) + burn_phi2_ini * EXP(-dt / nqse_burntime)

	    ! This model assumes NQSE burning burn O16 and Mg24 into Si28
            dburn_phi2 = MAX(burn_phi2_new - burn_phi2_ini, 0.0D0)

            epsilon2(j,k) = eps_ini + (burn_ene1 * xiso(j,k,co16) + burn_ene2 * xiso(j,k,cmg24)) * dburn_phi2
            burn_phi2(j,k) = burn_phi2_new

            yiso(j,k) = yiso(j,k) - dburn_phi2 * xiso_a(co16) / aion(co16) - dburn_phi2 * xiso_a(cne20) / aion(cmg24) + &
                        dburn_phi2 / aion(csi28)

            IF(burn_phi2_new == 1.0D0) THEN
               nse_flag(j,k) = 2
            ELSE
               nse_flag(j,k) = 1
            ENDIF

         ENDIF

      ENDDO
   ENDDO

   CALL boundary1D_NM(burn_phi2, even)
   CALL boundary1D_NM(epsilon2, even)

   END SUBROUTINE burn_phase2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine caluclates the energy released by 
   ! the 2nd phase of nuclear burning, from 16O and 24Mg to 28Si
   ! the 3rd phase of nuclear burning, from 28Si to Fe-peak
   ! And then solve consistenly the final temp
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine burn_phase3
   use definition
   use levelset_module
   use ecaptable_module
   implicit none

   integer :: j, k

   real (DP) :: rho_mid, temp_mid
   real (DP) :: nqse_burntime
   real (DP) :: nse_burntime

   real (DP) :: dburn_phi2
   real (DP) :: dburn_phi3

   real (DP) :: burn_energy2
   real (DP) :: burn_energy3

   ! Local burn_phi
   real (DP) :: burn_phi1_new
   real (DP) :: burn_phi2_ini
   real (DP) :: burn_phi2_new
   real (DP) :: burn_phi3_ini
   real (DP) :: burn_phi3_new

   ! The old and new NSE binding energy
   real (DP) :: qash_old, qash_new

   ! The old and new averaged binding energy
   real (DP) :: qbar_old, qbar_new

   ! The old and new inverse atomic mass
   real (DP) :: yiso_old, yiso_new

   ! The expected NSE state
   real (DP) :: temp_fin, qash_fin, qbar_fin, eps_fin, yiso_fin

   ! Change of eps
   real (DP) :: eps_ini
   real (DP) :: deps

   ! Old and new electron fraction
   real (DP) :: ye_old, ye_new
   
   ! limited timestep
   real (DP) :: dt_step

   ! My patch
   ! New temperature and azbar
   real (DP) :: temp_new, eps_new
   real (DP) :: abar_new, zbar_new

   ! Energy from electron capture and neutrino loss
   real (DP) :: ecaprate_mid, eneurate_mid

   real (DP), parameter :: qbar_C12 = 8.23266D-3
   real (DP), parameter :: qbar_O16 = 8.55261D-3
   real (DP), parameter :: qbar_Mg24 = 8.86028D-3
   real (DP), parameter :: qbar_Si28 = 9.06265D-3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First initialize burn_qdot to avoid overlap
   burn_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

	 ! Pass to local variables

	 eps_ini = epsilon2(j,k)
         rho_mid = rho2(j,k)
	 temp_mid = temp2(j,k)

	 !if(temp_mid < 5.0D0) cycle
	
	 burn_phi1_new = burn_phi1(j,k)         
	 burn_phi2_ini = burn_phi2(j,k)
	 burn_phi3_ini = burn_phi3(j,k)

	 yiso_old = yiso(j,k)
	 qash_old = qash(j,k)
	 ye_old = ye2(j,k)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 if(flamegrid_flag(j,k) == 1 .or. detongrid_flag(j,k) == 1) then

   	    ! Step 1: Solve for the NSE equilibrium state

	    ! Get the explicit electron capture rate
            call getecaprate(rho_mid, temp_mid, ye_old, ecaprate_mid, eneurate_mid)
            ye_new = ye_old + ecaprate_mid * dt

       	    qbar_old = burn_phi2_ini * qash_old + (1.0D0 - burn_phi1_new) * xiso_a(cc12) * qbar_C12 + &
		       (1.0D0 - burn_phi2_ini) * (1.0D0 - xiso_a(cc12)) * qbar_o16 + &
		       (burn_phi1_new - burn_phi2_ini) * qbar_mg24 
	    !write(*,100) j, k, rho_mid, ye_old, eps_ini - qash_old
	    ! Include all energy change from weak interaction
	    call GetNSEstate2(rho_mid, ye_old, eps_ini - qash_old - burn_phi3_new * &
			      (8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt), & 
			      temp_fin, qash_fin, yiso_fin)
	    eps_fin = eps_ini - qash_old + qash_fin

	    !if(j == 1 .and. k == 1) then
	    !   write(*,*) eps_ini, qash_old, qbar_old
	    !   write(*,*) rho_mid, ye_old, eps_ini - qash_old
	    !   write(*,*) temp_fin, qash_fin, yiso_fin
	    !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    ! When the grid is completely burn, allow advanced burning
	    nqse_burntime = EXP(182.06D0 / temp_fin - 46.054D0) / 4.9282D-6
	    nse_burntime = EXP(196.02D0 / temp_fin - 41.646D0) / 4.9282D-6    

	    !if(j == 1 .and. k == 1) then
            !   write(*,*) nqse_burntime, nse_burntime
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  	    ! If there is burn phase 2, then mute the burn phase 2 here
	    !burn_phi2_new = burn_phi1_new * (1.0D0 - EXP(-dt / nqse_burntime)) + burn_phi2_ini * EXP(-dt / nqse_burntime)
	    burn_phi2_new = burn_phi2_ini
	    burn_phi3_new = burn_phi1_new * (1.0D0 - &
			    (1.0D0 / (1.0D0 - nse_burntime / nqse_burntime)) * EXP(-dt / nqse_burntime) - &
                            (1.0D0 / (1.0D0 - nse_burntime / nqse_burntime)) * EXP(-dt / nse_burntime)) + &
                            burn_phi2_ini / (1.0D0 - nse_burntime / nqse_burntime) * &
			    (EXP(-dt / nqse_burntime) - EXP(-dt / nse_burntime)) + &
                            burn_phi3_ini * EXP(-dt / nse_burntime)

	    ! If there is burn phase 2, then mute the burn phase 2 here
 	    !dburn_phi2 = MAX(burn_phi2_new - burn_phi2(j,k), 0.0D0)
	    dburn_phi2 = 0.0D0
	    dburn_phi3 = MAX(burn_phi3_new - burn_phi3(j,k), 0.0D0)

	    !if(j == 1 .and. k == 1) then
	    !   write(*,*) burn_phi1_new
            !   write(*,*) burn_phi2_ini, burn_phi2_new, dburn_phi2
            !   write(*,*) burn_phi3_ini, burn_phi3_new, dburn_phi3
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    ! Update yiso and qbar by the conservative scheme (Townsley, ApJ 2007)
 	    dt_step = MIN(dt, nqse_burntime)
	    qash_new = (burn_phi2_ini * (qash_old + &
		       (qash_fin - qash_old) * dt_step / nqse_burntime) + dburn_phi2 * qash_fin) / burn_phi2_new
	    yiso_new = (burn_phi2_ini * (yiso_old + &
		       (yiso_fin - yiso_old) * dt_step / nqse_burntime) + dburn_phi2 * yiso_fin) / burn_phi2_new
	    !temp_new = (burn_phi2_ini * (temp_mid + (temp_fin - temp_mid) * dt_step / nqse_burntime) + dburn_phi2 * temp_fin) / burn_phi2_new
	    !eps_new = (burn_phi2_ini * (eps_ini + (eps_fin - eps_ini) * dt_step / nqse_burntime) + dburn_phi2 * eps_ini) / burn_phi2_new
	      
	    !if(j == 1 .and. k == 1) then
            !   write(*,*) dt_step, qash_new, yiso_new
            !endif

  	    ! Get the implicit electron capture rate
	    !call getecaprate(rho_mid, temp_fin, ye_old, ecaprate_mid, eneurate_mid)
	    !ye_new = ye_old + ecaprate_mid * dt
	    
	    !if(j == 1 .and. k == 1) then
            !   write(*,*) ecaprate_mid, eneurate_mid, ye_new
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  	    ! Now we are ready to update the energy
	    qbar_new = burn_phi2_new * qash_new + (1.0D0 - burn_phi1_new) * xiso_a(cc12) * qbar_C12 + &
                       (1.0D0 - burn_phi2_new) * (1.0D0 - xiso_a(cc12)) * qbar_o16 + &
                       (burn_phi1_new - burn_phi2_new) * qbar_mg24

	    ! For implicit electron capture rate
	    !deps = (qbar_new - qbar_old) - burn_phi3_new * (8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt)

	    ! For explicit electron capture rate in the dQ formalism
	    deps = (qbar_new - qbar_old)	

	    ! For explicit electron capture rate in the deps formalism
	    !deps = eps_new - eps_ini

	    !if(j == 1 .and. k == 1) then
            !   write(*,*) qbar_new, deps
            !endif

	    ! Now paste all the new results to the hydro variables
	    burn_qdot(j,k) = deps

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    ! Get the expected azbar and then the expected epsilon
	    !abar_new = 1.0D0 / yiso_new
            !zbar_new = abar_new * ye_new
	    !call teos_epsilon(rho_mid, temp_new, abar_new, zbar_new, ye_new, eps_new)

	    ! Copy the data back to the hydro
	    ye2(j,k) = ye_new
	    !burn_qdot(j,k) = eps_new - eps_ini
	    !epsilon2(j,k) = eps_new		! For temp version and eps version
	    epsilon2(j,k) = eps_ini + deps	! For eps version
	    burn_phi2(j,k) = burn_phi2_new
	    burn_phi3(j,k) = burn_phi3_new
	    qash(j,k) = qash_new
	    yiso(j,k) = yiso_new    	  
	
	    !if(j == 1 .and. k == 1) then
	    !   write(*,*) 'Azbar = ', abar2(j,k), zbar2(j,k)
	    !   write(*,*) 'temp', temp_mid, temp_fin
	    !   write(*,*) 'yiso', yiso_old, yiso_new, yiso_fin
	    !   write(*,*) 'eps', eps_ini, eps_new, eps_fin
	    !   write(*,*) 'phi3', burn_phi3_ini, burn_phi3_new
	    !   write(*,*) 'Qash', qash_old, qash_new
	    !   write(*,*) 'Qbar', qbar_old, qbar_new
	    !   write(*,*) 'Eps', eps_ini, eps_ini + deps
	    !endif

	    !if(j == 1 .and. k == 1) then
            !   write(*,*) ye_new, eps_ini + deps
            !   write(*,*) burn_phi2_new, burn_phi3_new
	    !   write(*,*) qash_new, yiso_new
	    !   read(*,*)
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 endif

      enddo
   enddo

   !write(*,*)

   call boundary1D_NM(ye2, even)
   call boundary1D_NM(temp2, even)
   call boundary1D_NM(burn_phi2, even)
   call boundary1D_NM(burn_phi3, even)
   call boundary1D_NM(qash, even)
   call boundary1D_NM(yiso, even)
   call boundary1D_NM(epsilon2, even)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   100 format (2I4, 5ES13.5)

   end subroutine burn_phase3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine check if the burn_phi1, burn_phi2, burn_phi3
   ! are well defined
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE check_burnPhi
   USE definition
   IMPLICIT NONE

   ! Integer 
   INTEGER :: j, k

   DO k = 1, length_step_z_2, 1
      DO j = 1, length_step_r_2, 1

	 IF(burn_phi1(j,k) > 1.0D0) then
	    burn_phi1(j,k) = 1.0D0
	 ELSEIF(burn_phi1(j,k) < 0.0D0) then
	    burn_phi1(j,k) = 0.0D0
	 ENDIF

	 IF(burn_phi2(j,k) > 1.0D0) then
            burn_phi2(j,k) = 1.0D0
         ELSEIF(burn_phi2(j,k) < 0.0D0) then
            burn_phi2(j,k) = 0.0D0
         ENDIF

         IF(burn_phi3(j,k) > 1.0D0) then
            burn_phi3(j,k) = 1.0D0
         ELSEIF(burn_phi3(j,k) < 0.0D0) then
            burn_phi3(j,k) = 0.0D0
         ENDIF

      ENDDO
   ENDDO

   END SUBROUTINE check_burnPhi

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine initialize the burn phase
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init_burnPhi
   USE definition
   IMPLICIT NONE

   ! Integer 
   INTEGER :: i, j, k

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         burn_phi1(j,k) = 0.0D0
	 burn_phi2(j,k) = 0.0D0
	 burn_phi3(j,k) = 0.0D0
	 yiso(j,k) = 1.0D0 / abar_ini
	 qash(j,k) = 9.06265D-3

      ENDDO
   ENDDO

   CALL Boundary1D_NM(burn_phi1, even)
   CALL Boundary1D_NM(burn_phi2, even)
   CALL Boundary1D_NM(burn_phi3, even)
   CALL boundary1D_NM(yiso, even)
   CALL boundary1D_NM(qash, even)

   END SUBROUTINE init_burnphi

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !
   ! This subroutine read the NSE table specific for 
   ! binding energy and equilibrium temperature
   !
   ! written by Leung Shing Chi in 2016
   ! written by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_nse_table2
   USE definition
   IMPLICIT NONE
 
   ! Integer and dummy variables !
   INTEGER :: i, j, k, k2
   REAL (DP) :: dummy

   ! Open file !
   OPEN(unit=500, file='../lib/nse_table_495iso.dat',action='read')

   ! Read it !
   DO i = 0, den_rowno_nse2, 1
      DO j = 0, ye_rowno_nse2, 1
	 READ(500,*) (nsetable2_head(i,j,k2), k2 = 1, 3)
	 DO k = 0, ent_rowno_nse2, 1
            READ(500,*) (nsetable2_binde(i,j,k,k2), k2 = 1, 6)
	    nsetable2_binde(i,j,k,4) = nsetable2_binde(i,j,k,4) / 1.0D9
	 ENDDO
	 READ(500,*)
      ENDDO
   ENDDO

   !write(*,*) 'Header'
   !write(*,*) (nsetable2_head(0,0,k2), k2 = 1, 3)
   !write(*,*) (nsetable2_binde(0,0,0,k2), k2 = 1, 6)
   !write(*,*) (nsetable2_head(0,1,k2), k2 = 1, 3)
   !write(*,*) (nsetable2_binde(0,1,0,k2), k2 = 1, 6)
   !write(*,*) (nsetable2_head(1,0,k2), k2 = 1, 3)
   !write(*,*) (nsetable2_binde(1,0,0,k2), k2 = 1, 6)

   ! Close it !
   CLOSE(500)

   END SUBROUTINE read_nse_table2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine search from the NSE table to
   ! get the equilibrium temperature for a given
   ! material enthalpy
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GetNSEstate2(rho_in, ye_in, neteps_in, temp_out, qash_out, yiso_out)
   USE definition
   IMPLICIT NONE

   REAL (DP) :: rho_in, ye_in, neteps_in
   REAL (DP) :: temp_out, qash_out, yiso_out

   INTEGER :: igrid, jgrid, kgrid1, kgrid2, kgrid3, kgrid4
   REAL (DP) :: digrid, djgrid, dkgrid1, dkgrid2, dkgrid3, dkgrid4

   igrid = INT((LOG10(rho_in * 6.171D17) - 6.0D0) / 0.1D0)
   digrid = ((LOG10(rho_in * 6.171D17)- 6.0D0) / 0.1D0 - DBLE(igrid))

   jgrid = INT((ye_in - 0.44D0) / 0.0025D0)
   djgrid = ((ye_in - 0.44D0) / 0.0025D0 - DBLE(jgrid)) 

   IF(igrid < 0) THEN
      igrid = 0
      digrid = 0.0D0
   ELSEIF(igrid > den_rowno_nse2) THEN
      igrid = den_rowno_nse2
      digrid = 0.0D0
   ENDIF

   IF(jgrid < 0) THEN
      jgrid = 0
      djgrid = 0.0D0
   ELSEIF(jgrid > ye_rowno_nse2) THEN
      jgrid = ye_rowno_nse2
      djgrid = 0.0D0
   ENDIF

   kgrid1 = INT((neteps_in - nsetable2_head(igrid,jgrid,1)) / nsetable2_head(igrid,jgrid,3))
   dkgrid1 = ((neteps_in - nsetable2_head(igrid,jgrid,1)) / nsetable2_head(igrid,jgrid,3) - DBLE(kgrid1)) 

   IF(kgrid1 < 0) THEN
      kgrid1 = 0
      dkgrid1 = 0.0D0             
   ELSEIF(kgrid1 > ent_rowno_nse2) THEN
      kgrid1 = ent_rowno_nse2
      dkgrid1 = 0.0D0
   ENDIF

   kgrid2 = INT((neteps_in - nsetable2_head(igrid+1,jgrid,1)) / nsetable2_head(igrid+1,jgrid,3))               
   dkgrid2 = ((neteps_in - nsetable2_head(igrid+1,jgrid,1)) / nsetable2_head(igrid+1,jgrid,3) - DBLE(kgrid2))

   IF(kgrid2 < 0) THEN
      kgrid2 = 0
      dkgrid2 = 0.0D0
   ELSEIF(kgrid2 > ent_rowno_nse2) THEN
      kgrid2 = ent_rowno_nse2
      dkgrid2 = 0.0D0
   ENDIF

   kgrid3 = INT((neteps_in - nsetable2_head(igrid,jgrid+1,1)) / nsetable2_head(igrid,jgrid+1,3))               
   dkgrid3 = ((neteps_in - nsetable2_head(igrid,jgrid+1,1)) / nsetable2_head(igrid,jgrid+1,3) - DBLE(kgrid3))

   IF(kgrid3 < 0) THEN
      kgrid3 = 0
      dkgrid3 = 0.0D0
   ELSEIF(kgrid3 > ent_rowno_nse2) THEN
      kgrid3 = ent_rowno_nse2
      dkgrid3 = 0.0D0
   ENDIF

   kgrid4 = INT((neteps_in - nsetable2_head(igrid+1,jgrid+1,1)) / nsetable2_head(igrid+1,jgrid+1,3))               
   dkgrid4 = ((neteps_in - nsetable2_head(igrid+1,jgrid+1,1)) / nsetable2_head(igrid+1,jgrid+1,3) - DBLE(kgrid4)) 

   IF(kgrid4 < 0) THEN
      kgrid4 = 0
      dkgrid4 = 0.0D0
   ELSEIF(kgrid4 > ent_rowno_nse2) THEN
      kgrid4 = ent_rowno_nse2
      dkgrid4 = 0.0D0
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The guess temp_out at lower Ye
   temp_out = nsetable2_binde(igrid,jgrid,kgrid1,4) * (1.0D0 - digrid) * (1.0D0 - djgrid) * (1.0D0 - dkgrid1) + &
	         nsetable2_binde(igrid+1,jgrid,kgrid2,4) * digrid * (1.0D0 - djgrid) * (1.0D0 - dkgrid2) + &
	         nsetable2_binde(igrid,jgrid+1,kgrid3,4) * (1.0D0 - digrid) * djgrid * (1.0D0 - dkgrid3) + &
                 nsetable2_binde(igrid+1,jgrid+1,kgrid4,4) * digrid * djgrid * (1.0D0 - dkgrid4) + &
	         nsetable2_binde(igrid,jgrid,kgrid1+1,4) * (1.0D0 - digrid) * (1.0D0 - djgrid) * dkgrid1 + &
                 nsetable2_binde(igrid+1,jgrid,kgrid2+1,4) * digrid * (1.0D0 - djgrid) * dkgrid2 + &
	 	 nsetable2_binde(igrid,jgrid+1,kgrid3+1,4) * (1.0D0 - digrid) * djgrid * dkgrid3 + &
                 nsetable2_binde(igrid+1,jgrid+1,kgrid4+1,4) * digrid * djgrid * dkgrid4

   qash_out = nsetable2_binde(igrid,jgrid,kgrid1,5) * (1.0D0 - digrid) * (1.0D0 - djgrid) * (1.0D0 - dkgrid1) + &   
                 nsetable2_binde(igrid+1,jgrid,kgrid2,5) * digrid * (1.0D0 - djgrid) * (1.0D0 - dkgrid2) + &
                 nsetable2_binde(igrid,jgrid+1,kgrid3,5) * (1.0D0 - digrid) * djgrid * (1.0D0 - dkgrid3) + &
                 nsetable2_binde(igrid+1,jgrid+1,kgrid4,5) * digrid * djgrid * (1.0D0 - dkgrid4) + &
                 nsetable2_binde(igrid,jgrid,kgrid1+1,5) * (1.0D0 - digrid) * (1.0D0 - djgrid) * dkgrid1 + &
                 nsetable2_binde(igrid+1,jgrid,kgrid2+1,5) * digrid * (1.0D0 - djgrid) * dkgrid2 + &
                 nsetable2_binde(igrid,jgrid+1,kgrid3+1,5) * (1.0D0 - digrid) * djgrid * dkgrid3 + &
                 nsetable2_binde(igrid+1,jgrid+1,kgrid4+1,5) * digrid * djgrid * dkgrid4

   yiso_out = nsetable2_binde(igrid,jgrid,kgrid1,6) * (1.0D0 - digrid) * (1.0D0 - djgrid) * (1.0D0 - dkgrid1) + &   
                 nsetable2_binde(igrid+1,jgrid,kgrid2,6) * digrid * (1.0D0 - djgrid) * (1.0D0 - dkgrid2) + &
                 nsetable2_binde(igrid,jgrid+1,kgrid3,6) * (1.0D0 - digrid) * djgrid * (1.0D0 - dkgrid3) + &
                 nsetable2_binde(igrid+1,jgrid+1,kgrid4,6) * digrid * djgrid * (1.0D0 - dkgrid4) + &
                 nsetable2_binde(igrid,jgrid,kgrid1+1,6) * (1.0D0 - digrid) * (1.0D0 - djgrid) * dkgrid1 + &
                 nsetable2_binde(igrid+1,jgrid,kgrid2+1,6) * digrid * (1.0D0 - djgrid) * dkgrid2 + &
                 nsetable2_binde(igrid,jgrid+1,kgrid3+1,6) * (1.0D0 - digrid) * djgrid * dkgrid3 + &
                 nsetable2_binde(igrid+1,jgrid+1,kgrid4+1,6) * digrid * djgrid * dkgrid4

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Debug !
   !write(*,*) igrid, jgrid, kgrid1
   !write(*,*) kgrid2, kgrid3, kgrid4
   !write(*,*) digrid, djgrid, dkgrid1
   !write(*,*) dkgrid2, dkgrid3, dkgrid4
   !write(*,*) nsetable2_binde(igrid,jgrid,kgrid1,6), nsetable2_binde(igrid,jgrid,kgrid1,6), nsetable2_binde(igrid,jgrid,kgrid1,6)
   !write(*,*) nsetable2_binde(igrid+1,jgrid+1,kgrid4,6), nsetable2_binde(igrid,jgrid,kgrid1+1,6), nsetable2_binde(igrid+1,jgrid,kgrid2+1,6)
   !write(*,*) nsetable2_binde(igrid,jgrid+1,kgrid3+1,6), nsetable2_binde(igrid+1,jgrid+1,kgrid4+1,6), yiso_out

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The guess temp_out at higher Ye

   END SUBROUTINE GetNSEstate2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine get the abar and zbar from
   ! yiso (See Townlsey et al. 2007)
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE Find_AZBar_fromBPhi
   USE definition
   IMPLICIT NONE

   ! Integer !
   INTEGER :: j, k

   ! Do the loop !
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         abar2(j,k) = 1.0D0 / yiso(j,k)
	 zbar2(j,k) = abar2(j,k) * ye2(j,k)

      ENDDO
   ENDDO

   CALL boundary1D_NM(abar2, even)
   CALL boundary1D_NM(zbar2, even)

   END SUBROUTINE Find_AZbar_fromBPhi

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! The subroutine tries to give the initial energy release
   ! due to carbon burning for the inner part of the star
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculate the burning of Variable 1
   ! using the burning variable algorithm (SEe Townsley et al., 2017)
   !
   ! Written by Leung SHing CHi in 2016
   ! Updated by Leung SHing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE burn_phase1_ini
   USE definition
   USE levelset_module
   IMPLICIT NONE

   ! Integer !
   INTEGER :: i, j, k

   ! Change of fraction for both level sets
   REAL (DP) :: x1, x2

   ! Energy input !
   REAL (DP) :: def_energy, det_energy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First initialize flame_qdot and deton_qdot
   flame_qdot(:,:) = 0.0D0
   deton_qdot(:,:) = 0.0D0
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1             
      DO j = 1, length_step_r_part_2, 1

         ! Calculate the local change of flame_ratio and deton_ratio
         x1 = flame_ratio(j,k)

         ! If there is change in flame_ratio, add in def. energy
         IF(x1 > 0.0D0) THEN

	    ! For CO matter
            !def_energy = 6.2762D-4 * xiso_a(cc12) * x1

	    ! For ONe matter
            def_energy = 2.4712D-4 * xiso_a(cne20) * x1

            flame_qdot(j,k) = def_energy
            epsilon2(j,k) = epsilon2(j,k) + def_energy
            burn_phi1(j,k) = MIN(burn_phi1(j,k) + x1, 1.0D0)

         ENDIF

      ENDDO
   ENDDO
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Copy the results to ghost cells 
   CALL boundary1D_NM(epsilon2, even)
   CALL boundary1D_NM(burn_phi1, even)

   END SUBROUTINE burn_phase1_ini

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the energy release 
   ! By assuming 12C --> 24Mg
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE burn_phase1b
   USE definition
   USE levelset_module
   IMPLICIT NONE

   ! Dummy variables
   INTEGER :: j, k

   ! Flag for finding temp
   integer :: flag_notfindtemp

   ! Local variables
   REAL (DP) :: rho_mid, temp_mid, vol_mid

   ! Change of fraction for both level sets
   REAL (DP) :: x1, x2

   ! Local Atwood Number
   REAL (DP) :: at_no

   ! Ne-22 binding energy                       
   !real (DP), parameter :: qbar_He4 = 7.57791D-3
   !real (DP), parameter :: qbar_C12 = 8.23266D-3
   !real (DP), parameter :: qbar_O16 = 8.55261D-3
   !real (DP), parameter :: qbar_Ne20 = 8.61316D-3
   !real (DP), parameter :: qbar_Mg24 = 8.86028D-3
   !real (DP), parameter :: qbar_Si28 = 9.06265D-3
   !real (DP), parameter :: qbar_Ne22 = 8.66501D-3
   !real (DP), parameter :: qbar_Ni56 =  9.27328D-3

   ! local flame energy, defined by Dqbar_C12 * xiso_ini(c12) + Dqbar_Ne22 * metallicity
   ! i.e. 6.2762E-4 * (xiso_ini(c12) - Z/2) + 1.9527E-4 * Z
   ! Example: Z = 0.00, flame_ene = 3.13810E-4
   ! Example: Z = 0.02, flame_ene = 3.11439E-4
   ! Example: Z = 0.04, flame_ene = 3.09068E-4
   ! Example: Z = 0.06, flame_ene = 3.06698E-4
   ! Example: Z = 0.08, flame_ene = 3.04327E-4
   ! Example: Z = 0.10, flame_ene = 3.01956E-4
   real (DP), parameter :: flame_ene = 3.11439D-4

   ! Ash composition
   real (DP), dimension(totalion) :: x_fuel1 = (/0.0D0, 0.51D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
   real (DP), dimension(totalion) :: x_ash1 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.51D0, 0.0D0, 0.0D0/)

   ! local deton energy, defined by Dqbar_C12 * xiso_ini(c12) + Dqbar_Ne22 * metallicity
   ! i.e. 6.2762E-4 * (xiso_ini(c12) - Z/2) + 1.9527E-4 * Z
   real (DP), parameter :: deton_ene = 3.11439D-4  ! He to Ne

   ! Ash composition
   real (DP), dimension(totalion) :: x_fuel2 = (/1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
   real (DP), dimension(totalion) :: x_ash2 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/)

   ! Other variables
   real (DP) :: esum, nqse_burntime2, ratio, generate, eps_old

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !if(debug_flag == 1) write(*,*) 'In Flame'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First, initialize the grid, to erase all past data
   flame_qdot(:,:) = 0.0D0
   deton_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !call findenergy
   !write(*,*) 'Before: ', energy2
   !esum = 0.0D0
   generate = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Then, find the energy release due to nuclear fusion
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

         rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)
	 vol_mid = vol2(j,k)
	 eps_old = epsilon2(j,k)

         ! Just compute the change of fraction occupies
         ! by deflagration or detonation
         !x1 = flame_ratio(j,k) - flame_ratio_old(j,k)
	 x1 = MIN(flame_ratio(j,k) - flame_ratio_old(j,k), xiso(j,k,cc12) / 0.51D0)
	 x2 = MIN(deton_ratio(j,k) - deton_ratio_old(j,k), xiso(j,k,cc12) / 0.51D0)
         !x2 = MIN(deton_ratio(j,k) - deton_ratio_old(j,k), xiso(j,k,che4))

         ! Remember to switch flame_ini as well
         ! When there is a change in area fraction
         ! Change the chemical composition accordingly
         ! And inject temperature
         if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min .and. rho_mid < rho2_flame_max) then
         !if(x1 > 0.0D0 .and. rho_mid > rho2_deton_min) then
            flame_qdot(j,k) = flame_ene * x1
            epsilon2(j,k) = epsilon2(j,k) + flame_qdot(j,k) 
            xiso(j,k,:) = xiso(j,k,:) - x1 * x_fuel1(:) + x1 * x_ash1(:)     
            burn_mass = burn_mass + x1 * rho_mid * vol_mid
         endif
        
         !Repeat the same procedure for detonation
         if(x2 > 0.0D0 .and. rho_mid > rho2_deton_min .and. rho_mid < rho2_deton_max) then
         !if(x2 > 0.0D0 .and. rho_mid > 1.62D-13) then

	    !nqse_burntime2 = 3.4838D0 * (rho_mid / 1.62D-10)**(-2)
	    !ratio = MIN(1.0D0, dt / nqse_burntime2)

	    deton_qdot(j,k) = deton_ene * x2 
            epsilon2(j,k) = epsilon2(j,k) + deton_qdot(j,k) 
            xiso(j,k,:) = xiso(j,k,:) - x2 *  x_fuel1(:) + x2 * x_ash1(:)
            burn_mass = burn_mass + x2 * rho_mid * vol_mid

         endif
 	 generate = generate + rho2(j,k)*(epsilon2(j,k) - eps_old)*vol2(j,k)
      enddo
   enddo
   energy_input = energy_input + generate

   ! Copy the results to ghost cells
   !call boundary1D_NM(epsilon2, even)
   !CALL BOUNDARY2D_X

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   100 format(10ES15.6)   
   101 format(A6, 2I5, 10ES15.6)

   end subroutine burn_phase1b   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the energy release
   ! By assuming 12C --> 24Mg at the beginning of 
   ! the simulation
   ! Written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE burn_phase1b_ini
   USE definition
   USE levelset_module
   IMPLICIT NONE

   ! Dummy variables
   INTEGER :: j, k                

   ! Flag for finding temp
   INTEGER :: flag_notfindtemp

   ! Local variables              
   REAL (DP) :: rho_mid, temp_mid    

   ! Change of fraction for both level sets
   real (DP) :: x1, x2

   ! Local Atwood Number
   real (DP) :: at_no

   ! local flame energy
   real (DP), parameter :: flame_ene = 6.7002D-4

   ! Ash composition
   REAL (DP), dimension(totalion) :: x_fuel = (/0.0D0, xc12_ini1, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
   REAL (DP), dimension(totalion) :: x_ash = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, xc12_ini1, 0.0D0, 0.0D0/)

   ! Other variables
   REAL (DP) :: esum

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !IF(debug_flag == 1) wRITE(*,*) 'In Flame'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First, initialize the grid, to erase all past data
   flame_qdot(:,:) = 0.0D0
   deton_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !call findenergy
   !write(*,*) 'Before: ', energy2
   !esum = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   ! Debug !
   WRITE(*,*) 'Check'
   WRITE(*,*) length_step_z_min_part_2, length_step_r_part_2, length_step_z_part_2
   WRITE(*,*) rho2(1,1), temp2(1,1), flame_ratio(1,1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Then, find the energy release due to nuclear fusion
   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)

         ! Just compute the change of fraction occupies
         ! by deflagration or detonation
         x1 = flame_ratio(j,k) 

         ! Remember to switch flame_ini as well
         ! When there is a change in area fraction
         ! Change the chemical composition accordingly
         ! And inject temperature
         !if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min .and. rho_mid < rho2_flame_max) then
         if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min) then
	    WRITE(*,*) 'Input energy', flame_ene*x1, j, k
            epsilon2(j,k) = epsilon2(j,k) + flame_ene * x1
            xiso(j,k,:) = xiso(j,k,:) - x1 * x_fuel(:) + x1 * x_ash(:)
            burn_mass = burn_mass + x1 * rho_mid * vol2(j,k)
         endif

     enddo
   enddo
 
   ! Copy the results to ghost cells   
   CALL BOUNDARY1D_NM(rho2, even)
   CALL BOUNDARY2D_X(xiso)
   CALL BOUNDARY1D_NM(epsilon2, even)

   end subroutine burn_phase1b_ini

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine calculates the energy release
   ! By assuming 16O + 24Mg --> 28Si
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by LEung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine burn_phase2b
   use definition
   use levelset_module
   implicit none

   ! Dummy variables
   integer :: j, k

   ! Flag for finding temp
   integer :: flag_notfindtemp               

   ! Local variables
   real (DP) :: rho_mid, temp_mid       

   ! Change of fraction for both level sets
   real (DP) :: x1, x2

   ! Local Atwood Number
   real (DP) :: at_no

   ! Ne-22 binding energy
   !real (DP), parameter :: qbar_C12 = 8.23266D-3
   !real (DP), parameter :: qbar_O16 = 8.55261D-3
   !real (DP), parameter :: qbar_Ne20 = 8.61316D-3
   !real (DP), parameter :: qbar_Mg24 = 8.86028D-3
   !real (DP), parameter :: qbar_Si28 = 9.06265D-3
   !real (DP), parameter :: qbar_Ne22 = 8.66501D-3
   !real (DP), parameter :: qbar_Ni56 = 9.27328D-3

   ! local flame energy, defined by 
   ! Energy release form 16O  --> 28Si: Dqbar_O16 * (xiso_ini(o16) - Z/2))
   ! Energy release from 24Mg --> 28Si: Dqbar_Mg24 * (xiso_ini(c12) + Z/2)
   ! i.e. 5.1004E-4 * (xiso_ini(o16) - Z/2)
   ! i.e. 2.2037E-4 * (xiso_ini(c12) + Z/2)
   ! Example: 
   real (DP) :: burn_ene1 = 5.1004D-4	! Energy release form 16O  --> 28Si
   real (DP) :: burn_ene2 = 2.0237D-4   ! Energy release from 24Mg --> 28Si  
   real (DP) :: burn_ene3 = 1.69537D-3   ! Energy release from 4He  --> 56Ni
   !real (DP) :: burn_ene3 = 6.6012D-4   ! Energy release from 20Ne  --> 56Ni
   !real (DP) :: burn_ene4 = 2.0163D-4   ! Energy release from 28Si  --> 56Ni

   ! local deton energy            
   real (DP) :: deton_ene 

   ! Ash composition
   real (DP), dimension(totalion) :: x_mid
   real (DP), dimension(totalion) :: x_fuel1 = (/0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
   real (DP), dimension(totalion) :: x_fuel2 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0/)
   real (DP), dimension(totalion) :: x_ash = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0/)

   real (DP), dimension(totalion) :: x_fuelb = (/1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
   real (DP), dimension(totalion) :: x_ashb = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/)

   ! NQSE timescale
   real (DP) :: nqse_burntime, nqse_burntime2, temp_nse, fnse
 
   ! Other variables
   real (DP) :: esum, generate, eps_old

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !if(debug_flag == 1) write(*,*) 'In Flame'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! First, initialize the grid, to erase all past data
   burn_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !call findenergy
   !write(*,*) 'Before: ', energy2
   !esum = 0.0D0 
   generate = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

   ! Then, find the energy release due to nuclear fusion
   do k = length_step_z_min_part_2, length_step_z_part_2, 1
      do j = 1, length_step_r_part_2, 1

	 ! Pass the data to local variables
	 rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)
         x_mid(:) = xiso(j,k,:)
  	 eps_old = epsilon2(j,k)

	 !if(rho_mid < 2.43D-11) cycle  !~3.6 standardized by the rho3, Z2, Ka1 case
         !if(rho_mid < 1.62D-11) cycle  !~3.6 standardized by the rho3, Z2, Ka3.2 (SQRT(10)) case
	 !if(ABS(p2(j+1,k) - p2(j-1,k)) / p2(j,k) >= 0.5D0 .or. &
	 !   ABS(p2(j,k+1) - p2(j,k-1)) / p2(j,k) >= 0.5D0) cycle


	 if(nse_flag(j,k) == 2) cycle

	 ! Only the flame_ratio or deton_ratio = 1 can burn
	 ! Meaning that completely C-burnt cells can start
	 ! Notice the MAX nature of flame_ratio or deton_ratio
	 !if(flame_loc_ratio(j,k) == 1.0D0 .and. rho_mid > rho2_flame_min) then
     	 if(burn_ratio(j,k) == 1.0D0 .and. rho_mid > 2.43D-11) then

            ! When there is a change in area fraction
            ! Change the chemical composition accordingly
            ! And inject temperature

	     ! Calculate the energy timescale
             !nqse_burntime = EXP(182.06D0 / temp_mid - 46.054D0) / 4.9282D-6
	     temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
             nqse_burntime = EXP(182.06D0 / temp_nse - 46.054D0) / 4.9282D-6

             ! Just compute the change of fraction by considering
             ! the maximum energy release equals to the amount of fuel remained
             x1 = x_mid(co16)
             x2 = x_mid(cmg24)

  	     if(dt > nqse_burntime) then
	   
	        ! Complete burning occurs
                burn_qdot(j,k) = burn_ene1 * x1 + burn_ene2 * x2
                epsilon2(j,k) = epsilon2(j,k) + burn_qdot(j,k)
                xiso(j,k,:) = x_mid(:) - x1 * x_fuel1(:) - x2 * x_fuel2(:) + (x1 + x2) * x_ash(:)

	        ! If it is completely burnt, then go to the next burning phase
	        nse_flag(j,k) = 2

	     else
		fnse =  -1.0604D0 * (dt / nqse_burntime)**2 + 2.0604D0 * (dt / nqse_burntime)**2
	        ! InComplete burning occurs
                burn_qdot(j,k) = (burn_ene1 * x1 + burn_ene2 * x2) * fnse !(dt / nqse_burntime)
                epsilon2(j,k) = epsilon2(j,k) + burn_qdot(j,k)
                xiso(j,k,:) = x_mid(:) + (-x1 * x_fuel1(:) - x2 * x_fuel2(:) + &
			      (x1 + x2) * x_ash(:)) * fnse !(dt / nqse_burntime)

	  	nse_flag(j,k) = 1

            endif

         endif
	 generate = generate + rho2(j,k)*(epsilon2(j,k) - eps_old)*vol2(j,k)
      enddo
   enddo                  
   energy_input = energy_input + generate   

   end subroutine burn_phase2b

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This suborutine did the last stage of burning from IME to NSE with IGEs
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE burn_phase3_ini
   USE definition
   USE levelset_module
   USE ecaptable_module
   IMPLICIT NONE

   ! Integer
   INTEGER :: j, k

   ! density and temperature
   REAL (DP) :: rho_mid, temp_mid

   ! nse related 
   REAL (DP) :: temP_nse
   REAL (DP) :: nqse_burntime_temp, nqse_burntime_rho
   REAL (DP) :: nqse_burntime
   REAL (DP) :: nse_burntime

   ! progress variables
   REAL (DP) :: dburn_phi2
   REAL (DP) :: dburn_phi3

   ! energy 
   REAL (DP) :: burn_energy2
   REAL (DP) :: burn_energy3     

   ! Local burn_phi
   REAL (DP) :: burn_phi1_new      
   REAL (DP) :: burn_phi2_ini             
   REAL (DP) :: burn_phi2_new
   REAL (DP) :: burn_phi3_ini
   REAL (DP) :: burn_phi3_new

   ! The old and new NSE binding energy
   REAL (DP) :: qash_old, qash_new

   ! The old and new averaged binding energy            
   REAL (DP) :: qbar_old, qbar_new            

   ! The old and new inverse atomic mass
   REAL (DP) :: yiso_old, yiso_new

   ! The expected NSE state
   REAL (DP) :: temp_fin, qash_fin, qbar_fin, eps_fin, yiso_fin

   ! Change of eps
   REAL (DP) :: eps_ini
   REAL (DP) :: deps

   ! Old and new electron fraction
   REAL (DP) :: ye_old, ye_new

   ! limited timestep
   REAL (DP) :: dt_step   

   ! My patch
   ! New temperature and azbar
   REAL (DP) :: temp_new, eps_new
   REAL (DP) :: abar_new, zbar_new

   ! Energy from electron capture and neutrino loss
   REAL (DP) :: ecaprate_mid, eneurate_mid

   ! change in binding energy q
   REAL (DP), parameter :: qbar_C12 = 8.23266D-3
   REAL (DP), parameter :: qbar_O16 = 8.55261D-3
   REAL (DP), parameter :: qbar_Ne20 = 8.61316D-3
   REAL (DP), parameter :: qbar_Mg24 = 8.86028D-3
   REAL (DP), parameter :: qbar_Si28 = 9.06265D-3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

   ! First initialize burn_qdot to avoid overlap           
   burn_qdot(:,:) = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1
      DO j = 1, length_step_r_part_2, 1

         ! Pass to local variables
         eps_ini = epsilon2(j,k)
         rho_mid = rho2(j,k)
         temp_mid = temp2(j,k)
         !if(temp_mid < 5.0D0) cycle

         burn_phi1_new = burn_phi1(j,k)
         burn_phi2_ini = burn_phi2(j,k)
         burn_phi3_ini = burn_phi3(j,k)

         yiso_old = yiso(j,k)
         qash_old = qash(j,k)
         ye_old = ye2(j,k) 

         ! Do the check 
         IF(temp2(j,k) < 5.0D0) CYCLE
         IF(ye2(j,k) < 0.20D0) CYCLE               
         IF(ye2(j,k) > 0.50D0) CYCLE   
         IF(rho2(j,k) < 1.62D-11) CYCLE
         IF(nse_flag(j,k) == 1) CYCLE     
         IF(burn_ratio(j,k) /= 1.0D0) CYCLE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
          
         ! If pass, procede
         IF(burn_ratio(j,k) == 1.0D0) THEN

            nse_flag(j,k) = 2

            ! Step 1: Solve for the NSE equilibrium state

            ! Get the explicit electron capture rate
            call getecaprate(rho_mid, temp_mid, ye_old, ecaprate_mid, eneurate_mid)
            ye_new = ye_old          

            qbar_old = burn_phi2_ini * qash_old + (1.0D0 - burn_phi1_new) * xiso_a(cc12) * qbar_C12 + &
                       (1.0D0 - burn_phi2_ini) * (1.0D0 - xiso_a(cc12)) * qbar_o16 + &
                       (burn_phi1_new - burn_phi2_ini) * qbar_mg24
            !write(*,100) j, k, rho_mid, ye_old, eps_ini - qash_old

            ! Include all energy change from weak interaction
            call GetNSEstate2(rho_mid, ye_old, eps_ini - qash_old, &
                              temp_fin, qash_fin, yiso_fin)
            eps_fin = eps_ini - qash_old + qash_fin

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

            ! When the grid is completely burn, allow advanced burning
	    
    	    temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
	    nqse_burntime_rho = EXP(182.06D0 / temp_nse - 46.054D0) / 4.9282D-6
            nqse_burntime_temp = EXP(182.06D0 / temp_fin - 46.054D0) / 4.9282D-6
   	    nqse_burntime = MIN(nqse_burntime_rho, nqse_burntime_temp)
            nse_burntime = EXP(196.02D0 / temp_fin - 41.646D0) / 4.9282D-6

            !if(j == 1 .and. k == 1) then        
            !   write(*,*) nqse_burntime, nse_burntime
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

            ! If there is burn phase 2, then mute the burn phase 2 here
            !burn_phi2_new = burn_phi1_new * (1.0D0 - EXP(-dt / nqse_burntime)) + burn_phi2_ini * EXP(-dt / nqse_burntime)
            burn_phi2_new = burn_phi2_ini                
            burn_phi3_new = burn_phi1_new

            ! If there is burn phase 2, then mute the burn phase 2 here
            !dburn_phi2 = MAX(burn_phi2_new - burn_phi2(j,k), 0.0D0)
            dburn_phi2 = 0.0D0
            dburn_phi3 = MAX(burn_phi3_new - burn_phi3(j,k), 0.0D0)

            !if(j == 1 .and. k == 1) then
            !   write(*,*) burn_phi1_new
            !   write(*,*) burn_phi2_ini, burn_phi2_new, dburn_phi2
            !   write(*,*) burn_phi3_ini, burn_phi3_new, dburn_phi3 
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

            ! Update yiso and qbar by the conservative scheme (Townsley, ApJ 2007)
            dt_step = MIN(dt, nqse_burntime)
            qash_new = qash_fin
            yiso_new = yiso_fin

            if(j == 1 .and. k == 1) then
               write(*,*) dt_step, qash_new, yiso_new                 
            endif

            ! Get the implicit electron capture rate
            !call getecaprate(rho_mid, temp_fin, ye_old, ecaprate_mid, eneurate_mid)
            !ye_new = ye_old + ecaprate_mid * dt      

            if(j == 1 .and. k == 1) then
               write(*,*) ecaprate_mid, eneurate_mid, ye_new           
            endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

            ! Now we are ready to update the energy
            qbar_new = burn_phi2_new * qash_new + (1.0D0 - burn_phi1_new) * xiso_a(cc12) * qbar_C12 + &
                       (1.0D0 - burn_phi2_new) * (1.0D0 - xiso_a(cc12)) * qbar_o16 + &
                       (burn_phi1_new - burn_phi2_new) * qbar_mg24     

            ! For implicit electron capture rate
            !deps = (qbar_new - qbar_old) - burn_phi3_new * (8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt)

            ! For explicit electron capture rate in the dQ formalism
            deps = (qbar_new - qbar_old)

            ! For explicit electron capture rate in the deps formalism
            !deps = eps_new - eps_ini

            if(j == 1 .and. k == 1) then
               write(*,*) qbar_new, deps    
            endif           

            ! Now paste all the new results to the hydro variables
            burn_qdot(j,k) = deps

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

            ! Get the expected azbar and then the expected epsilon
            !abar_new = 1.0D0 / yiso_new
            !zbar_new = abar_new * ye_new                             
            !call teos_epsilon(rho_mid, temp_new, abar_new, zbar_new, ye_new, eps_new)

            ! Copy the data back to the hydro       
            ye2(j,k) = ye_old
            !burn_qdot(j,k) = eps_new - eps_ini       
            !epsilon2(j,k) = eps_new            ! For temp version and eps version
            epsilon2(j,k) = eps_ini + deps      ! For eps version
            burn_phi2(j,k) = burn_phi2_new                             
            burn_phi3(j,k) = burn_phi3_new
            qash(j,k) = qash_new
            yiso(j,k) = yiso_new

        ENDIF

      ENDDO
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

   !write(*,*)                  

   CALL boundary1D_NM(ye2, even)            
   CALL boundary1D_NM(temp2, even)
   CALL boundary1D_NM(burn_phi2, even)        
   CALL boundary1D_NM(burn_phi3, even)             
   CALL boundary1D_NM(qash, even)
   CALL boundary1D_NM(yiso, even)
   CALL boundary1D_NM(epsilon2, even)

   100 format (2I4, 5ES13.5)

   END SUBROUTINE burn_phase3_ini

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine read the pre calculated NSE table of 495 isotope network
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_nse_table3
   use definition
   implicit none
 
   ! Dummy integers and variables !
   integer :: i, j, k, k2
   real (DP) :: dummy

   ! Open file !
   open(unit=500, file='../lib/nse_table_495iso_full.dat',action='read')

   ! Read it !
   do i = 0, den_rowno_nse3, 1
      do j = 0, ye_rowno_nse3, 1
         do k = 0, temp_rowno_nse3, 1
            read(500,*) (nsetable3_binde(i,j,k,k2), k2 = 1, 6)
         enddo
         read(500,*)
      enddo
   enddo
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Debug !
   !write(*,*) 'Header'
   !write(*,*) (nsetable2_head(0,0,k2), k2 = 1, 3)
   !write(*,*) (nsetable2_binde(0,0,0,k2), k2 = 1, 6)
   !write(*,*) (nsetable2_head(0,1,k2), k2 = 1, 3)
   !write(*,*) (nsetable2_binde(0,1,0,k2), k2 = 1, 6)
   !write(*,*) (nsetable2_head(1,0,k2), k2 = 1, 3)
   !write(*,*) (nsetable2_binde(1,0,0,k2), k2 = 1, 6)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Close file !
   close(500)

   end subroutine read_nse_table3

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine search from the NSE table to
   ! get the equilibrium temperature for a given
   ! material enthalpy
   ! written by Leung Shing Chi in 2016
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine GetNSEstate3(rho_in, ye_in, temp_in, binde_out, yiso_out, qash_out)
   use definition
   implicit none

   ! Input and output !
   real (DP), intent(in) :: rho_in, ye_in, temp_in
   real (DP), intent(Out) :: binde_out, yiso_out, qash_out

   ! Grid in table !
   integer :: igrid, jgrid, kgrid
   real (DP) :: digrid, djgrid, dkgrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign 
   igrid = INT((LOG10(rho_in * 6.171D17) - 5.3D0) / 0.25D0)
   digrid = ((LOG10(rho_in * 6.171D17)- 5.3D0) / 0.25D0 - DBLE(igrid))

   jgrid = INT((ye_in - 0.20D0) / 0.0025D0)
   djgrid = ((ye_in - 0.20D0) / 0.0025D0 - DBLE(jgrid)) 

   kgrid = INT((temp_in - 3.0D0) / 0.25D0)
   dkgrid = ((temp_in - 3.0D0) / 0.25D0 - DBLE(kgrid))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(igrid < 0) then
      igrid = 0
      digrid = 0.0D0
   elseif(igrid >= den_rowno_nse3) then
      igrid = den_rowno_nse3
      digrid = 0.0D0
   endif

   if(jgrid < 0) then
      jgrid = 0
      djgrid = 0.0D0
   elseif(jgrid >= ye_rowno_nse3) then
      jgrid = ye_rowno_nse3
      djgrid = 0.0D0
   endif

   if(kgrid < 0) then
      kgrid = 0
      dkgrid = 0.0D0
   elseif(kgrid >= temp_rowno_nse3) then
      kgrid = temp_rowno_nse3
      dkgrid = 0.0D0
   endif 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The guess temp_out at lower Ye
   binde_out = nsetable3_binde(igrid,jgrid,kgrid,4) * (1.0D0 - digrid) * (1.0D0 - djgrid) * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid+1,jgrid,kgrid,4) * digrid * (1.0D0 - djgrid) * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid,jgrid+1,kgrid,4) * (1.0D0 - digrid) * djgrid * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid+1,jgrid+1,kgrid,4) * digrid * djgrid * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid,jgrid,kgrid+1,4) * (1.0D0 - digrid) * (1.0D0 - djgrid) * dkgrid + &
                 nsetable3_binde(igrid+1,jgrid,kgrid+1,4) * digrid * (1.0D0 - djgrid) * dkgrid + &
                 nsetable3_binde(igrid,jgrid+1,kgrid+1,4) * (1.0D0 - digrid) * djgrid * dkgrid + &
                 nsetable3_binde(igrid+1,jgrid+1,kgrid+1,4) * digrid * djgrid * dkgrid

   yiso_out = nsetable3_binde(igrid,jgrid,kgrid,5) * (1.0D0 - digrid) * (1.0D0 - djgrid) * (1.0D0 - dkgrid) + &   
                 nsetable3_binde(igrid+1,jgrid,kgrid,5) * digrid * (1.0D0 - djgrid) * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid,jgrid+1,kgrid,5) * (1.0D0 - digrid) * djgrid * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid+1,jgrid+1,kgrid,5) * digrid * djgrid * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid,jgrid,kgrid+1,5) * (1.0D0 - digrid) * (1.0D0 - djgrid) * dkgrid + &
                 nsetable3_binde(igrid+1,jgrid,kgrid+1,5) * digrid * (1.0D0 - djgrid) * dkgrid + &
                 nsetable3_binde(igrid,jgrid+1,kgrid+1,5) * (1.0D0 - digrid) * djgrid * dkgrid + &
                 nsetable3_binde(igrid+1,jgrid+1,kgrid+1,5) * digrid * djgrid * dkgrid

   qash_out = nsetable3_binde(igrid,jgrid,kgrid,6) * (1.0D0 - digrid) * (1.0D0 - djgrid) * (1.0D0 - dkgrid) + &   
                 nsetable3_binde(igrid+1,jgrid,kgrid,6) * digrid * (1.0D0 - djgrid) * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid,jgrid+1,kgrid,6) * (1.0D0 - digrid) * djgrid * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid+1,jgrid+1,kgrid,6) * digrid * djgrid * (1.0D0 - dkgrid) + &
                 nsetable3_binde(igrid,jgrid,kgrid+1,6) * (1.0D0 - digrid) * (1.0D0 - djgrid) * dkgrid + &
                 nsetable3_binde(igrid+1,jgrid,kgrid+1,6) * digrid * (1.0D0 - djgrid) * dkgrid + &
                 nsetable3_binde(igrid,jgrid+1,kgrid+1,6) * (1.0D0 - digrid) * djgrid * dkgrid + &
                 nsetable3_binde(igrid+1,jgrid+1,kgrid+1,6) * digrid * djgrid * dkgrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !write(*,100) igrid, jgrid, kgrid, digrid, djgrid, dkgrid, binde_out, yiso_out, qash_out

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   100 format(3I4, 9ES15.7)

   end subroutine GetNSEstate3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! This subroutine finds the NSE state and correct the temperature
   ! according to energy conserving scheme
   ! This subroutine works in the burning variable algorithm (Townslet et al., 2007)
   !
   ! Written by Leung Shing Chi in 2016
   ! Updated by Leung Shing Chi in 2017
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE NSE3
   USE definition                               
   USE levelset_module
   USE ecaptable_module
   IMPLICIT NONE

   ! Flag for finding temperature
   INTEGER :: flag_notfindtemp             

   ! Dummy variables
   INTEGER :: i, j, k, k2

   ! Dummy variables
   REAL (DP) :: ye_sample = 0.5D0, dummy

   ! Local variables
   REAL (DP) :: abar_mid, zbar_mid, ye_mid

   ! Input variables
   REAL (DP) :: temp_beg    
   REAL (DP) :: eps_beg

   ! Trial variables
   REAL (DP) :: temp_mid                        
   REAL (DP) :: eps_mid         
   REAL (DP) :: rho_mid

   ! Change in temperature
   REAL (DP) :: dtemp_mid 

   ! Number of successful trial                
   INTEGER  :: count_digit

   ! Ecap rate and neutrino energy loss
   REAL (DP) :: ecaprate_mid
   REAL (DP) :: eneurate_mid

    ! Initial and Expected chemical composition
   REAL (DP), dimension(totalion) :: x_mid
   REAL (DP), dimension(totalion) :: x_burn

   ! Binding energy
   REAL (DP) :: binde_bf, binde_af, deps_nuc

   ! Prog variables
   REAL (DP) :: qash_mid, qash_af
   REAL (DP) :: yiso_mid, yiso_af

   ! Energy balance equation
   REAL (DP) :: check_e, check_e_last

   ! Timescale for burning to NSE
   REAL (DP) :: nse_burntime
   REAL (DP) :: nse_burntime_temp
   REAL (DP) :: nse_burntime_rho
   REAL (DP) :: temp_nse

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !IF(debug_flag == 1) WRITE(*,*) 'In NSE3'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialization
   ! Already done in burn_phase2b
   !burn_qdot(:,:) = 0.0D0
   total_ecap_nu_qdot = 0.0D0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO k = length_step_z_min_part_2, length_step_z_part_2, 1              
      DO j = 1, length_step_r_part_2, 1                    

         ! Do the checking first               
         IF(temp2(j,k) < 5.0D0) CYCLE
         IF(ye2(j,k) < 0.20D0) CYCLE
         IF(ye2(j,k) > 0.50D0) CYCLE
         IF(rho2(j,k) < 1.62D-11) CYCLE
         IF(burn_ratio(j,k) /= 1.0D0) CYCLE

         ! Special test if multi-stage burning is used
         IF(nse_flag(j,k) == 1) CYCLE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! If they pass then start the iteration
         nse_flag(j,k) = 2                
         burn_phi3(j,k) = 1.0D0

         ! The density is not changed, so 
         ! the trial value is also the 
         ! final value
         rho_mid = rho2(j,k)

         ! Give a trial temperature
         temp_beg = temp2(j,k)                  
         eps_beg = epsilon2(j,k)

         ! Also give the trial composition
         abar_mid = abar2(j,k)
         zbar_mid = zbar2(j,k)
         ye_mid = ye2(j,k)                
         binde_bf = qash(j,k)

         yiso_mid = yiso(j,k)
         qash_mid = qash(j,k)
     
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !if(flame_loc_ratio(j,k) == 1.0D0 .or. deton_loc_ratio(j,k) == 1.0D0) then
         if(temp2(j,k) < 15.00) then

            ! Scheme for electron capture
            ! Get the Ecap rate and energy loss
            call getecaprate(rho_mid, temp_beg, ye2(j,k), ecaprate_mid, eneurate_mid)
            !ecaprate_mid = 0.0D0; eneurate_mid = 0.0D0


            ! Update the Ye rate
            ! Note: No iteration is done here because
            ! the temperature sensitivity of Ecap
            ! rate is much smaller than NSE composition
            ye_mid = ye2(j,k) + ecaprate_mid * dt
            !ye_mid = ye2(j,k)                 

            ! Compute the binding energy of the 
            ! initial composition

            !write(*,*) 'Temp Before: ', temp_beg
            !write(*,*) 'Abar before: ', abar2(j,k)
            !write(*,*) 'Zbar before: ', zbar2(j,k)         
            !write(*,*) 'Eps before: ', eps_beg
            !write(*,*) 'Binde before: ', binde_bf
            !write(*,*)                 

            ! Prepare for the search of temperature
            count_digit = 0
            temp_mid = temp_beg
            dtemp_mid = 0.01D0 * temp_beg          

            ! Some debug stuff
            !if(j == 1 .and. k == 1) write(*,*) 'Now in grid', j, k
            !if(j == 1 .and. k == 1) write(*,*) 'Input density = ', rho_mid
            !if(j == 1 .and. k == 1) write(*,*) 'Old temp = ', temp_beg
            !if(j == 1 .and. k == 1) write(*,*) 'Current ye = ', ye_mid
            !if(j == 1 .and. k == 1) write(*,*) 'Current qasg = ', qash_mid
            !if(j == 1 .and. k == 1) write(*,101) (xiso(j,k,k2),k2=1,totalion)

            !if(global_time > 7.0D4) then
            !write(*,*) 'Now in grid', j, k
            !write(*,*) 'Input density = ', rho_mid
            !write(*,*) 'Old temp = ', temp_beg            
            !write(*,*) 'Current ye = ', ye_mid
            !write(*,*) 'Current qash = ', qash_mid
            !write(*,*) 'Current eps = ', eps_beg
            !write(*,*) 'Current yiso = ', yiso_mid
            !write(*,101) (xiso(j,k,k2),k2=1,totalion)
            !endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do i = 0, 400, 1   

               ! Now do the search

               !call nse_interface(temp_mid ,rho_mid, ye_mid, x_burn)

               ! Get the trial NSE state by the trial temp
               CALL getnsestate3(rho_mid,ye_mid,  temp_mid, qash_af, yiso_af, dummy)
               binde_af = qash_af

               ! Compute the trial binding energy
               !call compute_binde(x_burn, binde_af)

               ! Calculate the trial abar and zbar
               !call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)
               abar_mid = 1.0D0 / yiso_af
               zbar_mid = ye_mid * abar_mid

               ! Get the trial epsilon
               !call teos_epsilon(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid)
	       CALL HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)

               ! Calculate the binding energy change
               deps_nuc = binde_af - binde_bf

               ! Check if the energy is balanced
               check_e = eps_mid - eps_beg - deps_nuc - &
                         8.398670D-4 * ecaprate_mid * dt + eneurate_mid * dt

               ! Make sure you go to the right direction of dtemp
               if(i == 0) then
                  if(check_e > 0.0D0) dtemp_mid = -dtemp_mid
               endif

               ! Use shooting method    
               if(check_e_last * check_e < 0.0D0 .and. i /= 0) then
                  temp_mid = temp_mid - dtemp_mid
                  dtemp_mid = dtemp_mid * 0.1D0
                  temp_mid = temp_mid + dtemp_mid  
                  count_digit = count_digit + 1
               else
                  temp_mid = temp_mid + dtemp_mid
                  check_e_last = check_e
               endif

               IF(count_digit == 4) EXIT
               IF(i == 400 .or. temp_mid < 5.0D0) THEN
                  !stop 'Check NSE solver'
                  temp_mid = temp_beg        
                  eps_mid = eps_beg                  
                  !call nse_interface(temp_mid ,rho_mid, ye_mid, x_burn)
                  call getnsestate(rho_mid, temp_mid, x_burn)
                  exit
               ENDIF

               !if(j==1 .and. k==1) then
               !write(*,100) i, j, k, temp_mid, check_e, binde_af, yiso_af, eps_mid, 8.398670D-4 * ecaprate_mid * dt, eneurate_mid * dt
               !endif

            ENDDO

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid * 6.171D17)
            nse_burntime_rho = EXP(196.02D0 / temp_nse - 41.646D0) / 4.9282D-6
	    nse_burntime_temp = EXP(196.02D0 / temp_mid - 41.646D0) / 4.9282D-6
	    nse_burntime = MIN(nse_burntime_temp, nse_burntime_rho)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
            IF(dt > nse_burntime) THEN

               ! When things are done, move the refined
               ! trial results to the outpuit
               temp2(j,k) = temp_mid       
               epsilon2(j,k) = eps_mid
               ye2(j,k) = ye_mid   
               yiso(j,k) = yiso_af             
               qash(j,k) = binde_af 
               burn_qdot(j,k) = burn_qdot(j,k) + binde_af - binde_bf

               !if(global_time > 8.63D4) write(*,100) 1, j, k, eps_mid, ye_mid, yiso_af, binde_af, binde_af, binde_bf

            ELSE

               ! When things are partially done, use
               ! linear interpolation

               temp_mid = temp_beg + (temp_mid - temp_beg) * dt / nse_burntime
               yiso_af = yiso_mid + (yiso_af - yiso_mid) * dt / nse_burntime
               qash_af = qash_mid + (qash_af - qash_mid) * dt / nse_burntime

               binde_af = qash_af
               abar_mid = 1.0D0 / yiso_af
               zbar_mid = ye_mid * abar_mid
               !call teos_epsilon(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid)           
	       CALL HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)

               ! Now patch the result
               temp2(j,k) = temp_mid
               epsilon2(j,k) = eps_mid
               ye2(j,k) = ye_mid    
               yiso(j,k) = yiso_af
               qash(j,k) = qash_af
               burn_qdot(j,k) = burn_qdot(j,k) + binde_af - binde_bf

               !if(global_time > 8.63D4) write(*,100) 2, j, k, nse_burntime, temp_mid, eps_mid, ye_mid, yiso_af, qash_af, binde_af, binde_bf

            ENDIF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ELSE    

            ! If temperature is too high, ignore NSE
            ! composition change and focus on E-capture

            ! For detonation, no NSE energy is accounted
            temp_mid = temp2(j,k)

            call getecaprate(rho_mid, temp_beg, ye_mid, ecaprate_mid, eneurate_mid)

            ye2(j,k) = ye_mid + ecaprate_mid
            epsilon2(j,k) = eps_beg + 8.32696D-4 * ecaprate_mid * dt - eneurate_mid * dt
            call invert_helm_ed(epsilon2(j,k), rho2(j,k), abar2(j,k), zbar2(j,k), ye2(j,k), &
                                        temp2_old(j,k), temp2(j,k), dummy)


         ENDIF

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         total_ecap_nu_qdot = total_ecap_nu_qdot + vol2(j,k) * rho_mid * (eneurate_mid * dt)

      ENDDO
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Debug !
   !if(global_time > 8.63D4) write(*,*)
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Copy the results to ghost cells
   CALL boundary1D_NM(temp2, even)                 
   CALL boundary1D_NM(epsilon2, even)
   CALL boundary1D_NM(ye2, even)
   CALL boundary1D_NM(qash,even)
   CALL boundary1D_NM(yiso,even)
   CALL boundary1D_NM(burn_phi3,even)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   100 FORMAT(3I5, 20ES15.7)
   101 FORMAT(10ES15.7)

   END SUBROUTINE nse3

END MODULE Helmeos_module