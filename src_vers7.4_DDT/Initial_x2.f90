!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up every thing you need for running a star
! Based on the prototype by Wong Ka Wing
! Largely modified by Leung SHing Chi 
!
! To add new items, simply plug them here to make them
! appear when the simulation starts
!
! Written by Leung Shing Chi in 2015
! Commented by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIAL(n, fileno_ini, fileno_ini2)
USE DEFINITION 
USE HELMEOS_MODULE
!USE JETEXPLOSION_MODULE
USE PPT_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
USE NUSPEC_MODULE
USE GW_MODULE
USE WENO_MODULE
USE TURB_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!
!use NucEOS_module
!!!!!!!!!!!!!!!!!!!!

! The tool pack for naming files
INTEGER :: n, fileno_ini, fileno_ini2
INTEGER :: fileno_len
CHARACTER(LEN=256) :: fileno

! Signal for the end of PPT file
INTEGER :: flag_end

! Dummy variables
INTEGER :: i, j, k, m, o

! Signal for approving a given density
! Note: Obsolete
INTEGER :: flag_ans

! Dummy variables
REAL (DP) :: temp, dummy
REAL (DP) :: drho1c, check_m, check_m_last, log10rho1_c

! For reading the stellar profile
! for read_oldfile_flag = 2
REAL (DP) :: data_in(15)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Reminder if you read old files:
! If I want to read an old file, do the following
!
! Step 1: Change the starting FILENO in Main.f90
!         This should be a number bigger than the
!	  current file you want to read by one
! Step 2: Set flag_read_newfile = 1
! Step 3: Set the fileno_ini to the file I want to read
!	  Special attention is needed for PPT file
! Step 4: Set the found_deton_flag and found_deton_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First build up all database, EOS table and arrays !
WRITE(*,*) 'In initial'
WRITE(*,*) 'We shall now build the initial model for simulation'
WRITE(*,*)

! Usually new star is set
! The initial step is 0
n_backup = 0

! Give an initial timestep
dt = cfl * min(dx1, dx2) / 10000.0D0

! Atmospheric density !
rho1_a = rho1_c*rhofac_1
rho2_a = rho2_c*rhofac_2

!!!!!!!!!!!!!
!dt = 1.0D-3
!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get polytropic index !
CALL GETPOLY
   
! Read my EOS table
CALL EOSTABLE

! Allocate EOSTABLE !
IF (fermieosdm_flag == 1) THEN
	ALLOCATE(eostable1(eoslineno,2))
END IF
IF (fermieosnm_flag == 1) THEN
	ALLOCATE(eostable2(eoslineno,2))
END IF   

!!!!!!!!!!!!!!!!!!!!!
!CALL readEOSTable()
!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Build hydro variables !
WRITE(*,*) 'Build hydro variables'
CALL BUILDHYDRO
WRITE(*,*) 'Done building hydro variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) 'Build position table'
CALL GETGRIDNM
IF(DM_flag == 1) THEN
	CALL GETGRIDDM
END IF
WRITE(*,*) 'Done building position table'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) 'Get Legendre Function'
CALL FINDLEGENDRE
WRITE(*,*) 'Done getting Legendre function'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read helmholtz EOS
WRITE(*,*) 'Read Helmeos table'
IF(helmeos_flag == 1) THEN
   CALL buildhelm
   CALL read_helm_table()
   CALL initialize_network()
ENDIF
WRITE(*,*) 'Done reading Helmeos table'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(flame_flag == 1) THEN
   WRITE(*,*) 'Build level set variables'
   CALL buildLevelset
   WRITE(*,*) 'Done building level set'
   WRITE(*,*) 
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WRITE(*,*) 'Read Nuclear table'
!if(nuceos_flag == 1) then
!   call readtable(eos_table_name)
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(nuspec_flag == 1) THEN
   WRITE(*,*) 'Read neutirno table'
   CALL read_nutable
   WRITE(*,*) 'Done reading neutrino table'
   WRITE(*,*)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(tracer_flag == 1) THEN
   WRITE(*,*) 'Build PPT variables'
   CALL buildPPT
   WRITE(*,*) 'Done building PPT variables'
   WRITE(*,*)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(turb_flag == 1) THEN
   WRITE(*,*) 'Build sub-grid turbulence variables'
   CALL buildturb
   WRITE(*,*) 'Done building sub-grid turbulence variables'
   WRITE(*,*)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF(jetexp_flag == 1) THEN
!   WRITE(*,*) 'Build jet explosion variables'
!   CALL buildjetexp
!   WRITE(*,*) 'Done building jet explosion variables'
!   WRITE(*,*)
!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) 'Build WENO variables'
CALL BUILDWENO
WRITE(*,*) 'Done building WENO variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*) 'Read WENO constant'
CALL GETCONST
WRITE(*,*) 'Done reading WENO constant'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build initial models !
!IF(read_oldfile_flag == 0) THEN

   WRITE(*,*) 'Now build the initial model'

   ! Now set up a new star
   n = 1
   fileno_ini = -1
   fileno_ini2 = -1

   IF(initmodel_flag == 1) THEN

      ! Set up a real star
      WRITE(*,*) 'Build a real star'

      ! Solve the hydrostatic structure
      If(DM_flag == 1) THEN
	
	! Print out 
	WRITE (*,*) 'Bisection method'

	! Set step size !
	log10rho1_c = log10(rho1_c)
	drho1c = -0.01D0*log10rho1_c

	! Bisection method !
	DO o = 0, 100
		rho1_c = 10.0D0**(log10rho1_c)
		rho1_a = rho1_c*rhofac_1
		CALL GETRHO_2F
		CALL FINDMASS
		check_m = mass1 - mass_dm
		if(o == 0) then
			if(check_m > 0.0D0) drho1c = -drho1c
		endif
		if(check_m_last * check_m < 0.0D0 .and. o /= 0) then
			drho1c = -drho1c * 0.01D0
		end if
		log10rho1_c = log10rho1_c + drho1c
		check_m_last = check_m
		IF(ABS(check_m/mass_dm) < 1.0D-4) THEN
			EXIT
		END IF
		WRITE (*,*) 'DM Mass', mass_dm
	END DO
      ELSE
	CALL GETRHO_1F
      END IF

      ! Assign initial epsilon
      CALL GETEPSILON

      ! Assign initial velocity
      CALL GETVEL

      ! Assign initial equilibrium Ye
      IF(ecap_flag == 1) THEN
      	CALL GETYE
      END IF

   ELSE

      WRITE(*,*) 'Build according to default tests'
           
      ! Set up a quick density field for code test
      CALL getrho_template

      ! No atmosphere required
      rho1_a = 0.0D0
      rho2_a = 0.0D0

   ENDIF

   ! Atmospheric pressure ! 
   IF(helmeos_flag == 1) THEN
	CALL HELMEOS_RtoP(rho2_a, temp_a, abar_ini, zbar_ini, ye_a, p2_a, dpdrho2_a, dpdepsilon2_a)
   END IF

   ! Assign maximum density !
   rhomax1 = maxval(rho1)
   rhomax2 = maxval(rho2)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set up initial turbulence
   IF(turb_flag == 1) THEN
   	CALL GETTURB
   END IF

   ! Set up initial burnphi
   IF(burn_prog_flag == 1) THEN
   	call init_burnPhi
   END IF

   ! Set up the initial deflagration level-set
   IF(flame_flag == 1) THEN
   	CALL GETFLAME(1)
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Set up the initial detonation level-set
   !IF(deton_flag == 1) THEN
   !	CALL GETDETON(1)
   !END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Check if the density is okay
   IF(checkrho_flag == 1) THEN
	CALL CHECKRHO
   END IF

   ! Find the initial mass (for calculating potential)
   CALL FINDMASS
          
   ! Report some stellar parameters
   WRITE(*,*) 
   WRITE(*,*) 'Report initial model parameters'
   WRITE(*,*)
   WRITE(*,*) '---------------------------------'
   WRITE(*,*) 'Stellar parameters (rho_c, M, R):'
   WRITE(*,*) '---------------------------------'
   WRITE(*,*) rho1_c, mass1
   WRITE(*,*) rho2_c, mass2
   WRITE(*,*)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign initial Abar and Zbar 
   IF(xisotran_flag == 1) THEN
   	CALL FIND_AZBAR()
   END IF
   IF(burn_prog_flag == 1) THEN
   	CALL Find_AZBAR_fromBPhi()
   END IF

   ! Almost done
   WRITE(*,*) 'Finish construct initial model'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Assign the initial energy input from deflagration/detonation
   IF(burn_prog_flag == 1) THEN

      ! Set a large timestep at the beginning
      ! because we need the dt to reach 
      ! some equilibrium state
      dt = 1.0D4

      ! Here we burn C+Ne
      CALL burn_Phase1_ini
      CALL reconstructXIso              
      CALL Find_AZbar_fromBPhi()
      CALL Findhelmtemp

      ! Call the NQSE burning
      CALL burn_phase2
      CALL reconstructXIso
      CALL find_AZbar_fromBPhi()
      CALL findhelmtemp

      CALL burn_phase3_ini
      CALL reconstructXIso
      CALL find_AZbar_fromBPhi()
      CALL findhelmtemp

      ! Switch it back because the dt
      ! is just to search for equilibrium state
      dt = 1.0D0

   END IF

   ! For Flame !
   IF(xisotran_flag == 1) THEN
      IF(flame_flag == 1) THEN
	CALL FLAME_INI
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !IF(deton_flag == 1) CALL DETON_INI !DIRECT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ENDIF

   ! Print out !
   WRITE(*,*) 'Finished building the initial model'

!ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is for testing expanding grid, no use in normal case !
!CALL UPDATE (1) 
!DO j = -4, length_step_r_2 + 5
! DO k = -4, length_step_z_2 + 5
!	epsilon2(j,k) = epsilon2(j,k) + abs(phi2(j,k))
!	vel2_r(j,k) = 1.0D-4*rad2(j,k)/rmax2*2.0D0*cos2(j,k)*SQRT(1.0D0 - cos2(j,k)**2)
! END DO
!END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check density again
WRITE(*,*) 'Do the density check'
IF(checkrho_flag == 1) THEN
	CALL CHECKRHO	
END IF
WRITE(*,*) 'Done the initial density check'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute AZbar for both fuel and ash region
WRITE(*,*) 'Find AZbar and then Find temperature'
IF(flame_flag == 1) then
   IF(xisotran_flag == 1) THEN
    	CALL FIND_AZBAR()	
   END IF
   IF(burn_prog_flag == 1) THEN
   	CALL Find_AZBAR_fromBPhi()
   END IF
   IF(helmeos_flag == 1) THEN
   	CALL FINDHELMTEMP
   END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !if(nuceos_flag == 1) call findnuctemp
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDIF
WRITE(*,*) 'Done finding initial azbar and temperature'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the first primitive-conservative conversion
WRITE(*,*) 'Build conservative variables'
CALL FROMRVETOU (u_new1, u_new2)
CALL BACKUPCONS (u_new1, u_old1, u_new2, u_old2)
WRITE(*,*) 'Done building initial conservative variables'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate all the essential stellar quantities
WRITE(*,*) 'Find mass and energy'
CALL FINDMASS
CALL FINDENERGY
WRITE(*,*) 'Done finding initial mass and energy'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the items needed for SPATIAL
WRITE(*,*) 'Do Update'
CALL UPDATE (1) 
WRITE(*,*) 'Done initial update'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the initial dt
WRITE(*,*) 'Find initial dt'
CALL FindDt
Write(*,*) 'Done finding initial dt'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Miscellaneous setting
nse_flag = 0
energy_input = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Back up due to the use of reference data
CALL BACKUP_REF_ARRAY

! Print out !
WRITE(*,*) 'Final decided rho_a (DM, NM)'
WRITE(*,"(10ES15.7)") rho1_a, rho2_a

! Report what you get at center
IF(hemisphere_flag == 1) THEN
   CALL REPORT_CURRENT_DATA(1, length_step_z_2/2)
ELSEIF(hemisphere_flag == 0) THEN
   CALL REPORT_CURRENT_DATA(1, 1)
ELSE
   STOP 'Check the value of hemisphere_flag'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare gravitational wve 
! Not needed in this version
!IF(gravwave_flag == 1) THEN
!   WRITE(*,*) 'Do GetGravQ'
!   CALL getGravQ
!   WRITE(*,*) 'Done GetGravQ'
!   WRITE(*,*)
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare tracer particle scheme
IF(tracer_flag == 1) THEN !read_oldfile_flag /= 1 .and. 
   WRITE(*,*) 'Do getppt'
   CALL getppt
   WRITE(*,*) 'Done getppt'
   WRITE(*,*)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Finish initial...'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialize the polytropic index and exponent for !
! Polytropic equation of state 					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETPOLY
USE DEFINITION 
IMPLICIT NONE

! assign newtonian polytropic index and exponent for NM !
If (polyeosnm_flag == 1) THEN
	!k_2 = (3.0E0_DP ** (2.0E0_DP/3.0E0_DP) * pi ** (4.0E0_DP/3.0E0_DP) * hbar ** (2.0E0_DP) * ye2_old ** (5.0E0_DP/3.0E0_DP) ) &
	!	/ (5.0E0_DP * me2 * mb2 ** (5.0E0_DP/3.0E0_DP) )
	!gamma2 = (5.0E0_DP/3.0E0_DP)
	k_2 = (3.0E0_DP ** (1.0E0_DP/3.0E0_DP) * pi ** (2.0E0_DP/3.0E0_DP) * hbar * ye2_old ** (4.0E0_DP/3.0E0_DP) ) &
		/ (4.0E0_DP * mb2 ** (4.0E0_DP/3.0E0_DP) )
	gamma2 = 1.34D0
END IF

WRITE (*,*) 'report polytropic index for NM'
WRITE (*,*) 'k_2', k_2
WRITE (*,*) 'gamma2', gamma2
WRITE (*,*) 

! assign the polytopic index and exponent for DM if it exist !
If (DM_flag == 1) THEN
	k_1 = (3.0E0_DP ** (2.0E0_DP/3.0E0_DP) * pi ** (4.0E0_DP/3.0E0_DP) * hbar ** (2.0E0_DP) * ye1 ** (5.0E0_DP/3.0E0_DP) ) &
		/ (5.0E0_DP * me1 * mb1 ** (5.0E0_DP/3.0E0_DP) )
	gamma1 = (5.0E0_DP/3.0E0_DP)

	WRITE (*,*) 'report polytropic index for DM'
	WRITE (*,*) 'k_1', k_1
	WRITE (*,*) 'gamma1', gamma1
	WRITE (*,*) 
END IF

END SUBROUTINE