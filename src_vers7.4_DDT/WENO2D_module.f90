!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This module contains all the arrays and variables that are neccessary to run the hydro
! simulations of either pure hydro, or self gravitating fluid objects like stars
! In the 2D extension, all positions are defined as (r,z) where r is the distance from the z-axis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DEFINITION
IMPLICIT NONE
SAVE
INCLUDE "Parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Supplementary data for eos table

! Number of line in the polytropic EOS table
INTEGER :: eoslineno = 27644

! Number of lines in the other EOS tables
INTEGER :: eosline1, eosline2

! THe corresponding EOS tables for DM and NM
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: eostable1
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: eostable2

! Read error flag
INTEGER :: flag_error

! Back up 
INTEGER :: n_backup

! For grid identifications !
INTEGER :: r_grid1, z_grid1
INTEGER :: r_grid2, z_grid2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Boundary flag for u
! A large array is prepared
! so that no need to change even if we model a complex system
INTEGER :: bfac_r(100)
INTEGER :: bfac_z(100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables for WENO call !

! Minimum Eq. number for DM
INTEGER :: imin1
INTEGER :: iminsca1

! Maximum Eq. number for DM
INTEGER :: imax1
INTEGER :: imaxsca1

! Call code for the common variables of DM
! irho1 = DM density
! ivel1_r = DM r-velocity
! ivel1_z = DM z-velocity
! ivel1_p = DM p-velocity
! itau1 = DM energy density
INTEGER :: irho1
INTEGER :: ivel1_r
INTEGER :: ivel1_z
INTEGER :: ivel1_p
INTEGER :: itau1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Minimum Eq. number for NM
INTEGER :: imin2
INTEGER :: iminsca2

! Maximum Eq. number for NM
INTEGER :: imax2
INTEGER :: imaxsca2

! Call code for the common variables of NM
! irho2 = NM density
! ivel2_r = NM r-velocity
! ivel2_z = NM z-velocity
! ivel2_p = NM p-velocity
! itau2 = NM energy density
INTEGER :: irho2
INTEGER :: ivel2_r
INTEGER :: ivel2_z
INTEGER :: ivel2_p
INTEGER :: itau2
INTEGER :: ieps2
INTEGER :: iye2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Geometrical variables

! R and Z coordinate of the grid for DM !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: r1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: z1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: zF1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: rF1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: rad1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: cos1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vol1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: radbar1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: volbar1

! R and Z coordinate of the grid for NM !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: r2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: z2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: zF2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: rF2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: rad2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: cos2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vol2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: radbar2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: volbar2

! For spherical coordinate
REAL (DP), ALLOCATABLE, DIMENSION (:) :: COT_z1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: COT_z2

! Patch for fornax grid
! Derivative from the mapping function r = A SINH(x/A) + B
REAL (DP), ALLOCATABLE, DIMENSION (:) :: r1f_m
REAL (DP), ALLOCATABLE, DIMENSION (:) :: sca1_fac1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: sca1_fac2

! For NM !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: r2f_m
REAL (DP), ALLOCATABLE, DIMENSION (:) :: sca2_fac1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: sca2_fac2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Hydrodynamical variables
! The following are the hydro set for DM sector
! They include
! 1. density
! 2. R-velocity
! 3. Z-velocity
! 4. P-velocity
! 5. specific internal energy density
! 6. Thermo pressure
! 7. Pressure derivative w.r.t. density
! 8. Pressure derivative w.r.t. epsilon
! 9. Adiabatic index
! 10. Speed of sound
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: rho1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel1_r
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel1_z
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel1_p
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: epsilon1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: p1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dpdrho1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dpdepsilon1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: gamma_1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: cs1

! A supplement for the gravity
! 11. Gravitational potential experienced by DM 
! 12. Gravitational potential contribution from DM
! 13. Gravitational potential contribution from NM
! 14. R-dev of gravitational potential
! 15. Z-dev of gravitational potential
! 16. DM multipole moment
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi1_dm
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi1_nm
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi1_r
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi1_z
REAL (DP), ALLOCATABLE, DIMENSION (:) :: qpole1

! A second patch for the moving grid algorithm
! 17. R-velocity of the frame
! 18. Z-velocity of the frame
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel_frame_r1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel_frame_z1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The following are the hydro set for NM sector

! They include
! 1. density
! 2. R-velocity
! 3. Z-velocity
! 4. P-velocity
! 5. specific internal energy density
! 6. Thermo pressure
! 7. Pressure derivative w.r.t. density
! 8. Pressure derivative w.r.t. epsilon
! 9. Adiabatic index 
! 10. Electron fraction
! 11. Temperature
! 12. Speed of sound
! 13. Internal energy per unit volume
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: rho2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel2_r
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel2_z
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel2_p
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: epsilon2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: p2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dpdrho2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dpdepsilon2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: gamma_2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: ye2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: temp2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: cs2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: rhoe2

! A supplement for the gravity
! 14. Gravitational potential experienced by NM 
! 15. Gravitational potential contribution from DM
! 16. Gravitational potential contribution from NM
! 17. R-dev of gravitational potential
! 18. Z-dev of gravitational potential
! 19. NM multipole moment
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi2_dm
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi2_nm
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi2_r
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: phi2_z
REAL (DP), ALLOCATABLE, DIMENSION (:) :: qpole2

! Supplement for dual energy !
! 20. Total energy density for NM 
! 21. Derived internal energy density for NM 
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: bige2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: et2

! A third patch for backup variables
! 20. density from last step
! 21 temp from last step
! 22. epsilon from last step
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: rho2_old
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: temp2_old
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: epsilon2_old

! For pressure gradient 
! 23. dpdr
! 24. dpdz
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dp2dr
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: dp2dz

! A second patch for the moving grid algorithm
! 25. R-velocity of the frame
! 26. Z-velocity of the frame
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel_frame_r2
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: vel_frame_z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Legendre polynominal for DM and NM
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: legendre1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: legendre2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Mass coordinate !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: m_r1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: m_r2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: m_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For RK-Time evolution 
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: u_temp1, u_old1, u_new1, u2_dm, u3_dm, l3_dm, l1
REAL (DP), ALLOCATABLE, DIMENSION (:,:,:) :: u_temp2, u_old2, u_new2, u2_nm, u3_nm, l3_nm, l2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For output !
REAL (DP) :: rmax1, rmax2
REAL (DP) :: rhomax1, rhomax2

! Atmospheric values !
REAL (DP) :: epsilon1_a, dlfmmo1
REAL (DP) :: epsilon2_a, dlfmmo2
REAL (DP) :: p2_a, dpdrho2_a, dpdepsilon2_a

! For Fermi gas !
REAL (DP) :: a_max1, a_max2
REAL (DP) :: b_max1, b_max2

! Masses ! 
REAL (DP) :: mass1, mass2, mass_ash

! Total Energy !
REAL (DP) :: energy1, energy2, energy_input

! Individual energyies !
REAL (DP) :: energy1_kin, energy1_int, energy1_pot
REAL (DP) :: energy2_kin, energy2_int, energy2_pot

! Central Density !
REAL (DP) :: centralrho1, centralrho2

! Some logical count flag !
INTEGER :: count, potential_flag

! Global simulation time !
REAL (DP) :: global_time

! Variables of expanding grid !
real (DP) :: vel1_max, radius1, radius1_ini
real (DP) :: vel2_max, radius2, radius2_ini

! Boundaries !
real (DP) :: boundary1, boundary2

END MODULE DEFINITION