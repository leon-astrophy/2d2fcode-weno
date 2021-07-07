!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine initialized the WENO module
! Written by Leung Shing Chi in 2016
! The subroutine set up the primitive and conservative
! flux array by using the input from the Parameter.h
! When more flags is used, this subroutine will count
! all these addition quantities altogether
! Notice that if you want your own quantities to be 
! transport by the WENO scheme, you need to add your 
! own primitive and their declaration here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BUILDWENO
USE DEFINITION
USE HELMEOS_MODULE
USE LEVELSET_MODULE
USE TURB_MODULE
IMPLICIT NONE

! Dummy integer !
INTEGER :: dummy

! Initialize no_of_eq
no_of_eq = 0

! Initialize imax imin !
imin1 = 0
imax1 = 0
imin2 = 0
imax2 = 0

! Initialize max/min scalar equation !
iminsca1 = 0
iminsca2 = 0
imaxsca1 = 0
imaxsca2 = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find out how many equations are needed to solve

WRITE(*,*)
WRITE(*,*) 'Now we arrange no_of_eq according to the setting in Parameter.h'
WRITE(*,*) 'For each advectable scalar, no_of_eq increases by 1'
WRITE(*,*)

IF(runDM_flag == 1) THEN

   imin1 = 1

   ! DM density
   no_of_eq = no_of_eq + 1
   imax1 = no_of_eq
   irho1 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make irho1 = ', no_of_eq
   
   ! DM vel-r
   no_of_eq = no_of_eq + 1
   imax1 = no_of_eq
   ivel1_r = no_of_eq 
   bfac_r(no_of_eq) = -1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make ivel1_r = ', no_of_eq

   ! DM vel-z
   no_of_eq = no_of_eq + 1
   imax1 = no_of_eq
   ivel1_z = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = -1
   write(*,*) 'Make ivel1_z = ', no_of_eq

   IF(rotationdm_flag == 1) THEN
   	! DM vel-p
   	no_of_eq = no_of_eq + 1
   	imax1 = no_of_eq
   	ivel1_p = no_of_eq
   	bfac_r(no_of_eq) = -1
   	bfac_z(no_of_eq) = 1
   	write(*,*) 'Make ivel1_p = ', no_of_eq
   END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! DM epsilon, not needed
   !no_of_eq = no_of_eq + 1
   !imax1 = no_of_eq
   !itau1 = no_of_eq
   !bfac_r(no_of_eq) = 1
   !bfac_z(no_of_eq) = 1
   !write(*,*) 'Make itau1 = ', no_of_eq
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we do the NM sector

! Set up the code for NM hydro
! NM density
imin2 = no_of_eq + 1

no_of_eq = no_of_eq + 1
imax2 = no_of_eq
irho2 = no_of_eq
bfac_r(no_of_eq) = 1
bfac_z(no_of_eq) = 1
write(*,*) 'Make irho2 = ', no_of_eq

! NM vel-r
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_r = no_of_eq
bfac_r(no_of_eq) = -1
bfac_z(no_of_eq) = 1
write(*,*) 'Make ivel2_r = ', no_of_eq

! NM vel-z
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2_z = no_of_eq
bfac_r(no_of_eq) = 1
bfac_z(no_of_eq) = -1
write(*,*) 'Make ivel2_z = ', no_of_eq

! NM vel-p
If(rotationnm_flag == 1) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	ivel2_p = no_of_eq
	bfac_r(no_of_eq) = -1
	bfac_z(no_of_eq) = 1
	write(*,*) 'Make ivel2_p = ', no_of_eq
END IF

! NM epsilon
IF(nm_epsilon == 1) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	itau2 = no_of_eq
	bfac_r(no_of_eq) = 1
	bfac_z(no_of_eq) = 1
	write(*,*) 'Make itau2 = ', no_of_eq
END IF

! Dual energies !
IF(dual_energy == 1) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	ieps2 = no_of_eq
	bfac_r(no_of_eq) = 1
	bfac_z(no_of_eq) = 1
	write(*,*) 'Make ieps2 = ', no_of_eq
END IF

! Dummy integer !
dummy = no_of_eq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now check if we need to include extra component
! Set up the code for SGS
! turb_q
IF(turb_flag == 1) THEN
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   iturbq = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make iturbq = ', no_of_eq
ENDIF

! Set up the code for electron fraction
IF(Etran_flag == 1) then 

   ! Electron fraction
   no_of_eq = no_of_eq + 1   
   imaxsca2 = no_of_eq             
   imax2 = no_of_eq
   iye2 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make iye2 = ', no_of_eq

ENDIF

! Set up the code for deflagration
IF(deflevelset_flag == 1) THEN

   ! Deflagration level-set
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   iscaG1 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make iscaG1 = ', no_of_eq

ENDIF

! Set up the code for detonation
IF(detlevelset_flag == 1) THEN

   ! Detonation level-set
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   iscaG2 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make iscaG2 = ', no_of_eq

ENDIF

! Set up the code for chemical composition
IF(xisotran_flag == 1) THEN

   ! Helium-4
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ihe4 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make ihe4 = ', no_of_eq

   ! Carbon-12
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ic12 = no_of_eq   
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make ic12 = ', no_of_eq
   
   ! Oxygen-16
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   io16 = no_of_eq   
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   write(*,*) 'Make io16 = ', no_of_eq

   ! Neon-20
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ine20 = no_of_eq   
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make ine20 = ', no_of_eq

   ! Magnesium-24
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   img24 = no_of_eq   
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make img24 = ', no_of_eq

   ! Silicon-28
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   isi28 = no_of_eq   
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make isi28 = ', no_of_eq

   ! Nickel-56
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ini56 = no_of_eq   
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make ini56 = ', no_of_eq

ENDIF

IF(burn_prog_flag == 1) THEN

   ! Burn phi1 
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ibphi1 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make ibphi1 = ', no_of_eq

   ! Burn phi2
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ibphi2 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make ibphi2 = ', no_of_eq

   ! Burn phi3 
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   ibphi3 = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make ibphi3 = ', no_of_eq

   ! Yiso
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   iyiso = no_of_eq
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make Yiso = ', no_of_eq

   ! qash
   no_of_eq = no_of_eq + 1
   imaxsca2 = no_of_eq 
   imax2 = no_of_eq
   iqash = no_of_eq 
   bfac_r(no_of_eq) = 1
   bfac_z(no_of_eq) = 1
   WRITE(*,*) 'Make qash = ', no_of_eq

ENDIF

! Minimum scalar equation !
IF(imaxsca2 > 0) THEN
	iminsca2 = dummy + 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we know how many equations to be solved
! So set up the array accordingly
! Allocate the U array
! Allocate the source-term array

! For DM !
ALLOCATE(l1(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))
ALLOCATE(u2_dm(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))
ALLOCATE(u3_dm(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))
ALLOCATE(l3_dm(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))
ALLOCATE(u_old1(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))
ALLOCATE(u_new1(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))
ALLOCATE(u_temp1(-4:length_step_r_1+5, -4:length_step_z_1+5,imin1:imax1))

! For NM !
ALLOCATE(l2(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))
ALLOCATE(u2_nm(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))
ALLOCATE(u3_nm(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))
ALLOCATE(l3_nm(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))
ALLOCATE(u_old2(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))
ALLOCATE(u_new2(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))
ALLOCATE(u_temp2(-4:length_step_r_2+5, -4:length_step_z_2+5,imin2:imax2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now report the equation structure      
WRITE(*,*)
WRITE(*,*) 'Report equation number' 
WRITE(*,*) 'Dark matter ------------------------------------'
WRITE(*,*) 'imin1 = ', imin1
WRITE(*,*) 'imax1 = ', imax1
WRITE(*,*) 'iminsca1 = ', iminsca1
WRITE(*,*) 'imaxsca1 = ', imaxsca1 
WRITE(*,*)             
WRITE(*,*) 'Data structure'            
WRITE(*,*) 'Motherboard r1:', length_step_r_1
WRITE(*,*) 'Motherboard z1:', length_step_z_1
WRITE(*,*)    
WRITE(*,*) 'Normal matter ----------------------------------'
WRITE(*,*) 'imin2 = ', imin2
WRITE(*,*) 'imax2 = ', imax2
WRITE(*,*) 'iminsca2 = ', iminsca2
WRITE(*,*) 'imaxsca2 = ', imaxsca2
WRITE(*,*)             
WRITE(*,*) 'Data structure'            
WRITE(*,*) 'Motherboard r2:', length_step_r_2
WRITE(*,*) 'Motherboard z2:', length_step_z_2
WRITE(*,*)
WRITE(*,*) 'Finished BuildWENO'
WRITE(*,*)

END SUBROUTINE buildWENO