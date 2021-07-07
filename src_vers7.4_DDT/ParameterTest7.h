!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SNIa_PhaseVII version 1.0 (Last Modified: Leon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Foreword:
! Developed by Leung Shing Chi based on the WENO prototype developed by Wong Ka Wing in 2010
! In this version, I tried to summarize all the setting in this file
! So that minimum changes is required for modifying the code to 
! fit different geometry or EOS or physics.
! The code is designed to deal with supernovae simulations.
! But they are also fit enough to deal with ordinary hydrodynamical simulations
! For further information, refer to my old paper 
! Leung et al., MNRAS 454, 1238 (2015).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update in PhaseVII
! In this version some major changes are done
! 1. The use of boundary flag: I make it simpler for variables
! 2. I removed the temporaily set OpenMP function 
! 3. I optimize the coding in levelset_module
! 4. I optimized the coding in turb_module
! 5. I add the 2nd layer (pseudo-AMR) function (Now removed)
! 6. Documentation of this code
! 7. I updated the helmholtz table
! 8. I set the read/save function for a quick restart
! 9. I renamed the variables for easier understanding
! 10. I optimized the data structure in spatial.f90 
! 11. I patched the SR solver which I used for CC
! 12. I provided sample parameter.h for a quick start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here are the summary of each section
! Section 1: Constant and universal values
! Section 2: Core part of the simulation box
! Section 3: Properties of the star
! Section 4: Advanced property of star configuration
! Section 5: Central part of the code admin
! Section 6: Output setting
! Section 7: Mass cut function
! Section 8: Other Physcial Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1: Constant and universal values 
! This means you don'y change it until the universe changes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Double precision declarer, Define double precision

INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary flag notation for boundary conditions
! Range (0,1,2,3) stands for the scalar, R-type vector, Z-type scalar and RZ-type tensor

INTEGER, PARAMETER :: even = 0, oddR = 1, oddZ = 2, oddB = 3	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pi

REAL (DP), PARAMETER :: pi = 3.1415926535897932384626433832795E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Planck constant

REAL (DP), PARAMETER :: hbar = 1.1965E-76_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2: Core part of the simulation box
! This is where the system variables are controled
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag for 
! 0 = Cartesian coordinate
! 1 = Cylindrical coordinate (quadrant or hemisphere)

INTEGER, PARAMETER :: coordinate_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag for using hemisphere 
! 1 = Geometry in hemisphere
! Caution, be aware of the boundary conditions when using hemisphere

INTEGER, PARAMETER :: hemisphere_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The number of grid in the r-(x-) direction

INTEGER, PARAMETER :: length_step_r_1 = 100
INTEGER, PARAMETER :: length_step_r_2 = 200

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The number of grid in the z-(y-) direction
! Note: For an octant simulation, set length_step_z = length_step_r
! Note: For an quadrant simulation, set length_step_z = 2 * length_step_r + 1

INTEGER, PARAMETER :: length_step_z_1 = 100
INTEGER, PARAMETER :: length_step_z_2 = 200

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The number of highest grid in the r-(x-) direction 
! This is the minimum r-number for containing the whole star in the box
! This does not need to be changed unless you want a fixed smaller box

! For DM and NM !
INTEGER :: length_step_r_part_1 = length_step_r_1
INTEGER :: length_step_r_part_2 = length_step_r_2
								
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The number of highest grid in the z-(y-) direction  
! This is the minimum r-number for containing the whole star in the box
! This does not need to be changed unless you want a fixed smaller box

! For DM and NM !
INTEGER :: length_step_z_part_1 = length_step_z_1
INTEGER :: length_step_z_part_2 = length_step_z_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The number of lowest grid in the z-(y-) direction  
! This is the minimum r-number for containing the whole star in the box
! This does not need to be changed unless you want a fixed smaller box

! For DM and NM !
INTEGER	:: length_step_z_min_part_1 = 1
INTEGER	:: length_step_z_min_part_2 = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fundamental Grid size 
! 
! Notice that it can change when the grid expands                   

REAL (DP), PARAMETER :: dx1_ini = 0.01D0
REAL (DP), PARAMETER :: dx2_ini = (1.0D0/200.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid size 
!
! Notice that it can change when the grid expands
! Hint: The standard grid size (Note: 1 unit = 1.4774 km)

REAL (DP) :: dx1 = dx1_ini			
REAL (DP) :: dx2 = dx2_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control on the accuracy of the initial model
! The initial model resolution is 10^(ini_acc) finer than 
! the hydro one for a quiet initial model 
! 
! Note: Also change dx_more accordingly for the grid size
! Defined as dx_more = dx / 10^{ini_acc)

INTEGER, PARAMETER :: ini_acc = 2                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hardware configuratios of WENO
! This is an important parameters for bridging smooth flow to discontinuous flow
! without leading to singularity
!
! Note: Do NOT change this unless you know what you are doing

REAL (DP), PARAMETER :: smallpara = 1.0D-40    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid size for the initial model
!
! Hint: One set it to 1% of dx for accuracy

REAL (DP), PARAMETER :: dxmore = min(dx1_ini, dx2_ini) / (10.0D0**ini_acc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Number of iteration
! 
! Usually I pick simply a very large number
! and truncate the code by some conditions
! Hints: Time unit (1 unit = 4.9282 micro s)

INTEGER, PARAMETER :: time_step = 10000000			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cournat-Friedrich-Levy constant
!
! Definedas dt = cfl * dx / MAX(vel + cs)
! Also can vary depending if there is some very strong shock
! which one needs to process it carefully to avoid code crash

REAL (DP) :: cfl = 0.20E0_DP				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maximum time to be simulated in the model
!
! This avoids producing too much data and focus
! the simulation to the specific durection 

REAL (DP), PARAMETER :: total_time = 3.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time step for one iteration
!
! This will be changed immediated by FindDt during iteration

REAL (DP) :: dt = 0.001D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 3: Properties of the star
! This is where the star properties are modified
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial central density of dark matter part

REAL (DP), PARAMETER :: rho1_c = 3.0E0_DP * 1.6190E-10_DP !1.0E-13_DP				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial central density of normal matter
!
! Hint1: One sets rho2_c = 7.11D-11 for 1.05 solar mass WD (Ye = 0.5)
! Hint2L One sets rho2_c = 4.83D-9 for 1.37 solar mass WD (Ye = 0.5)

REAL (DP), PARAMETER :: rho2_c = 3.0D0 * 1.6190E-9_DP !4.83E-9_DP		  			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Angular veloicty of the initial model
!
! Note: Setting this one does not mean setting a rotating star density profile, it just tells
! the initial velocity distribution for rigid rotation.
! The density profiles need to be input by hand due to the complexity of the initial model
! Refer to hachisu (1986) for how to constructure 

REAL (DP), PARAMETER :: omega_ini = 0.0D-5			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! An isothermal star is constructed
! Here the temperature is declared
! Initial NM temperature
!
! Note: 1 unit = 1 billion kelvin

REAL (DP) :: temp_ini = 1.0D-1                                  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial composition of the star
! I notice that sometimes a star can have multiple layers
! So one can set at most of 3 layers of different compositions
!
! This includes
! 1. Initial He-4 mass fraction
! 2. Initial C-12 mass fraction
! 3. Initial O-16 mass fraction
! 4. Initial electron fraction

REAL (DP), PARAMETER :: xhe4_ini1  = 0.00D0                        
REAL (DP), PARAMETER :: xc12_ini1  = 0.50D0                       
REAL (DP), PARAMETER :: xo16_ini1  = 0.50D0                       
REAL (DP), PARAMETER :: xne20_ini1 = 0.00D0
REAL (DP), PARAMETER :: ye_ini1    = 0.5D0                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag for the second layer
! 1 = Allow the composition to change 
! Condition needs to be put in by hand in GetRho.f90

INTEGER, PARAMETER :: initmodel_layer2_flag = 0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2nd layer composition
!
! This includes
! 1. Initial He-4 mass fraction
! 2. Initial C-12 mass fraction           
! 3. Initial O-16 mass fraction
! 4. Initial electron fraction

REAL (DP), PARAMETER :: xhe4_ini2  = 0.00D0			
REAL (DP), PARAMETER :: xc12_ini2  = 0.00D0                       
REAL (DP), PARAMETER :: xo16_ini2  = 0.55D0                       
REAL (DP), PARAMETER :: xne20_ini2 = 0.45D0
REAL (DP), PARAMETER :: ye_ini2    = 0.50D0                          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag for the third layer
! 1 = Allow the initial model to contain multiple layers
! Condition needs to be put in by hand in GetRho.f90

INTEGER, PARAMETER :: initmodel_layer3_flag = 0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3rd layer composition
!
! This includes
! 1. Initial He-4 mass fraction
! 2. Initial C-12 mass fraction           
! 3. Initial O-16 mass fraction
! 4. Initial electron fraction            

REAL (DP), PARAMETER :: xhe4_ini3  = 0.0D0                                                
REAL (DP), PARAMETER :: xc12_ini3  = 0.0D0   
REAL (DP), PARAMETER :: xo16_ini3  = 0.55D0   
REAL (DP), PARAMETER :: xne20_ini3 = 0.45D0 
REAL (DP), PARAMETER :: ye_ini3    = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here we set the atmosphere of the star
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Density of the dark matter atmosphere
!
! Note: One can still modify it in RungeKutta.f90 and GetRho.f90

REAL (DP), PARAMETER :: rhofac_1 = 1.0E-10_DP
REAL (DP) :: rho1_a = rho1_c * rhofac_1 !1.0E-13_DP			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Density of the normal matter atmosphere
!
! Note: One can still modify it in RungeKutta.f90 and GetRho.f90

REAL (DP), PARAMETER :: rhofac_2 = 1.0E-10_DP
REAL (DP) :: rho2_a = rho2_c * rhofac_2 !1.0E-14_DP			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Velocity of the dark matter atmosphere

REAL (DP), PARAMETER :: vel1_a = 0.0E0_DP                                  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Velocity of the normal matter atmosphere

REAL (DP), PARAMETER :: vel2_a = 0.0E0_DP                                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atmosphere temperature

REAL (DP), PARAMETER :: temp_a = 1.0D-1		

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atmosphere electron fraction

REAL (DP), PARAMETER :: ye_a = 0.5D0				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Limiter of temperature
! Note: Setting them according to your EOS choice
!
! temp_max = Maximum temperature allowed
! temp_min = Minimum temperature allowed

REAL (DP), PARAMETER :: temp_max = 7.0D1		
REAL (DP), PARAMETER :: temp_min = 1.0D-4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Limiter of electron fraction
! Note: Setting them according to your EOS choice
!
! temp_max = Maximum Ye allowed
! temp_min = Minimum Ye allowed

REAL (DP), PARAMETER :: ye_max = 0.5D0
REAL (DP), PARAMETER :: ye_min = 0.2D0				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 4: Advanced property of star configuration
! This gives more fine-tune control of the star model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing gravity
! 
! Notice that the Poisson solver give gravity in 
! Cylindrical coordinates. For other coordinates
! the corresponding Poisson equation is needed

INTEGER, PARAMETER :: w_gravity_i = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Variables for controling the (successive over-)relaxation method
! They Can be tuned to fit the accuracy
! But note that the time required for solving the system
! also grow exponentially when the tolerance becomes stricter
!
! relax_max = Maximum run time in relaxation of the potential
! tolerance = Maximum residue allowed for approval in relaxation
! sor_weight = Weight for successive-over relaxation in gauss-seidel method

INTEGER, PARAMETER :: lmax = 16
INTEGER , PARAMETER :: relax_max = 300000                        
REAL (DP), PARAMETER :: tolerance = 3.0E-8_DP                   
REAL (DP), PARAMETER :: sor_weight = 0.99D0    
	               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 5: Central part of the code admin
! This controls how the simulation should be run
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Number of equations to be solved by WENO

INTEGER :: no_of_eq = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing DM to exist
!
! 1 = run with DM

INTEGER, PARAMETER :: DM_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing DM to advect
!
! 1 = run with DM dynamics

INTEGER, PARAMETER :: runDM_flag = 0	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for including vel_p as a dynamical variables
!
! 1 = Run with angular motion allowed

! For DM !
INTEGER, PARAMETER :: rotationdm_flag = 0

! For NM !
INTEGER, PARAMETER :: rotationnm_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for solving energy equation

! Solve epsilon equation? !
INTEGER, PARAMETER :: nm_epsilon = 1
INTEGER, PARAMETER :: dual_energy = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Map the conservative quantities to coarser grid

INTEGER, PARAMETER :: mapping_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for star models
! 1 = Set the initial condition to hydrostatic equilibrium
! 0 = Put in handmade initial profiles 

INTEGER, PARAMETER :: initmodel_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for output profiles
! 1 = output profile 

INTEGER, PARAMETER :: output_flag = 1                   	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for boundary condition
! The boundary flag is defined by four scalar
! 1st one for R-inner boundary
! 2nd one for R-outer boundary
! 3rd one for Z-inner boundary
! 4th one for Z-outer boundary
!
! 0 = periodic
! 1 = reflecting boundary (depend on scalar/vector)
! 2 = outgoing (1st derivative = 0)

INTEGER, PARAMETER :: boundary_flag(4) = (/2,2,2,2/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for checking the density
! 1 = Reset the meshes with density below rho_a

INTEGER, PARAMETER :: checkrho_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for checking the simulation box
! 1 = Only simulate the box area with matter, other meshes without matter is ignored

! For DM
INTEGER, PARAMETER :: checkstepdm_flag = 0

! For NM
INTEGER, PARAMETER :: checkstepnm_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for checking the velocity
! 1 = Reset the meshes with velocity close to 1

INTEGER, PARAMETER :: checkvel_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for update the dt at each step
! 1 = Update dt based on the local velocity and sound speed

INTEGER, PARAMETER :: updatedt_flag = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for the initial profiles
! It works if you set initmoddel_flag = 0
! Each number stand for the corresponding test suite
!
! List of providied tests
! 1 = TORO 1D Shocktube
! 2 = Sedov spherical blast
! 3 = 2D explosion
! 4 = Kelvin-Helmholtz
! 5 = Rayleight-Taylor
! 6 = Kelvin-Helmholtz-New
! 7 = Gresho problem
! 8 = Implosion 

INTEGER, PARAMETER :: testmodel_flag = 7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! EOS flag, Noticd that you can set multiple EOS 
! simultaneously to deal with physics with a wide 
! range of density temperature and compositions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For DM

! Flags for the Polytropic EOS
! 1 = use Polytropic EOS
INTEGER, PARAMETER :: polyeosdm_flag = 0

! Flags for the Fermi Gas EOS
! 1 = use Fermi gas EOS
INTEGER, PARAMETER :: fermieosdm_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM

! Flags for the Polytropic EOS
! 1 = use polytropic EOS
INTEGER, PARAMETER :: polyeosnm_flag = 1

! Flags for the Fermi Gas EOS
! 1 = use Fermi gas EOS
INTEGER, PARAMETER :: fermieosnm_flag = 0

! Flags for the HELMHOLTZ EOS, refer Timmes (1999)
! 1 = use helmholtz EOS
INTEGER, PARAMETER :: helmeos_flag = 0                    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing isotope transport
! 1 = Allow isotope mass fraction to be advected
! 1 = with transport of isotopes

INTEGER, PARAMETER :: xisotran_flag = 0                     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing the nuclear burning using progressive variables
! See Townley et al., 2007 for further details
! 1 = Use progressive variables rather than isotopes
!
! Notice that xisotran_flag and fusion_flag are needed
! in order to switch on this function

INTEGER, PARAMETER :: burn_prog_flag = 0			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing the electron fraction to be transported
! 1 = Allow Ye to be advected

INTEGER , PARAMETER :: etran_flag = 0              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flags for allowing electron capture for white dwarf collapse
! A density-Ye relation is provided in Abdikamalov et al. (2010)
! This makes the local Ye to be changed according to its density
! during the white dwarf collapses into a neutron star
! 1 = Use the analytic fit for the density-Ye relation

INTEGER, PARAMETER :: Ecap_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flags for the polytropic prefactor
! One takes pressure = k * rho^(gamma)
! Or in general pressure = (gamma - 1) * rho * epsilon

! For DM
REAL (DP) :: k_1 = 5.0D0 / 3.0D0

! For NM 
REAL (DP) :: k_2 = 5.0D0 / 3.0D0        	                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flags for the polytropic index
! One takes pressure = constant * rho^(gamma)
! Or in general pressure = (gamma - 1) * rho * epsilon

! For DM
REAL (DP) :: gamma1 = 5.0D0 / 3.0D0

! For NM 
REAL (DP) :: gamma2 = 5.0D0 / 3.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here are the series of flags which controls the nuclear reactions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for switch on the nuclear reaction option
! 1 = Run with the nuclear burning package

INTEGER, PARAMETER :: fusion_flag = 0			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for using level set scalar 1 for tracking burning
! 1 = Set up the level set field for deflagration

INTEGER, PARAMETER ::  deflevelset_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for using level set scalar 2 for tracking burning
! 1 = Set up the level set field for detonation
! Note: The physics can be changed by altering the burning scheme

INTEGER, PARAMETER :: detlevelset_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing the interaction between level set 1 and energy input
! 1 = Allow energy input by deflagration

INTEGER, PARAMETER :: flame_flag = 0	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for allowing the interaction between level set 1 and energy input
! 1 = Allow energy input by detonation + finding detonation

INTEGER, PARAMETER :: deton_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Flag for allowing 1st step burning input for level set 1 & 2
! 1 = Allow energy input by carbon burning

INTEGER, PARAMETER :: carburn_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Flag for allowing 2nd step burning input for level set 1 & 2
! 1 = Allow energy input by advanced burning

INTEGER, PARAMETER :: advburn_flag = 0                           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Flag for allowing final burning input for level set 1 & 2
! 1 = Allow energy input by NSE evolution

INTEGER, PARAMETER :: convert_nse_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flags for sub-grid turbulence, See Niemeyer et al. (1995a)
! This subroutine is physical consistent if the grid
! sizein all directions are uniformed. Switch off this
! flag if you have curveilinear coordinates, such as 
! 3D cylindrical or 2D spherical coordinates
!
! 1 = Run with sub-grid scale turbulence

INTEGER, PARAMETER :: turb_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for neutrino spectra
! This calculates the thermal neutrino based on
! the analytic formula by Odrzywolek et al. (2005) and Mischizek et al. (2005)
! Pair- and plasmon- neutrino are included
!
! 1 = Switch on the neutrino spectra calculator

INTEGER, PARAMETER :: nuspec_flag = 0			

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for post-processing by tracer particles
!
! 1 = With tracer particle to record the thermodynamics history

INTEGER, PARAMETER :: tracer_flag = 0
								
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for calculating the gravitational wave by the quadruple model
!
! 1 = Calculate the energy loss due to gravitational wave emission

INTEGER, PARAMETER:: gravwave_flag = 0				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for the option sponge near the atmosphere
! This model is good when you need a quiet star
! But tuning the initial model appears to be 
! important as well
!
! 1 = Put a sponge on the star surface

INTEGER, PARAMETER :: sponge_flag = 0				

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Flag for a global output
!
! This can be manipulated if one needs an output
! in all profiles outside the regular output time

LOGICAL :: output_file = .true.		

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 6: Output setting
! This sets how frequent each type of profile is output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output configurations
! Notice they are not parameter because
! there are cases you want more frequent
! output at some dynamical scenarios
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Physical time interval for all log file

REAL (DP) :: output_logtime = 0.1D0                     
REAL (DP) :: output_logtime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Physical time interval for each hydro profile

REAL (DP) :: output_profiletime = 0.5D0
REAL (DP) :: output_profiletime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Physical time interval for each level set profile

REAL (DP), PARAMETER :: output_flametime = 0.5D0                  	
REAL (DP), PARAMETER :: output_flametime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Physical time interval for each SGS turb profile

REAL (DP), PARAMETER :: output_turbtime = 0.5D0
REAL (DP), PARAMETER :: output_turbtime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Physical time interval for each isotope profiles

REAL (DP), PARAMETER :: output_Helmtime = 0.5D0		
REAL (DP), PARAMETER :: output_Helmtime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Physical time interval for each tracer profiles

REAL (DP), PARAMETER :: output_PPTtime = 0.5D0		
REAL (DP), PARAMETER :: output_PPTtime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Physical time interval for each backup
! But I muted it because it can be done by the new read/save option

REAL (DP), PARAMETER :: output_backuptime = 0.5D0
REAL (DP), PARAMETER :: output_backuptime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 7: Mass cut function (SS9)
! This part controls the star with a mass cut
! and the jet formation physics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The flag for mass cut and enclosed mass

INTEGER, PARAMETER :: masscut_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The radius of the mass cut 

REAL (DP), PARAMETER :: rad_cut = 609.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The enclosed mass inside cut
! hint: 1.42 solar mass

REAL (DP), PARAMETER :: mass_cut = 1.42D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The flag controls the energy injection by jet

INTEGER, PARAMETER :: jetexp_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This controls the energy deposition rate of the jet
! Hint: the rate of energy deposition (3.2855D-7 --> 120x10^51 erg/s)

REAL (DP), PARAMETER :: Dedep = 3.2855D-7                      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This controls the energy deposition angle
! Hint: the angle spanned by the jet  (0.2818D0 = 15 degrees)

REAL (DP), PARAMETER :: theta_jet = 0.2818D0                    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This controls the ratio between internal energy vs kinetic energy

REAL (DP), PARAMETER :: fth = 1.0D-3                      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This controls the  energy deposition time 
! Hint: jet to pump energyu (2.5364D4 --> 0.125s)

REAL (DP), PARAMETER :: jet_time = 2.5364D4                     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This controls the adiabatic index of the jet matter

REAL (DP), PARAMETER :: adia_jet = 4.0D0/3.0D0        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This controls the Lorentz factor of the jet

REAL (DP), PARAMETER :: gamma_jet = 100.0D0                     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This controls the electron scattering opacitiy (basically a constant)
! Now it is kap_e = 0.1 cm^2 / g

REAL (DP), PARAMETER :: kap_e = 9.1629D21

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 8: Further constants
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Baryonic mass for dark matter

REAL (DP), PARAMETER :: mb1 = 8.4158E-59_DP                     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Baryonic mass for normal matter

REAL (DP), PARAMETER :: mb2 = 8.4158E-58_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fermionic mass (dark matter particles) for dark matter

REAL (DP), PARAMETER :: me1 = 8.4158E-59_DP                     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fermionic mass (electrons) for normal matter

REAL (DP), PARAMETER :: me2 = 4.5771E-61_DP                     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Degeneracy parameter for DM

REAL (DP), PARAMETER :: ye1 = 1.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Degeneracy parameter for DM

REAL (DP), PARAMETER :: ye2_old = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiplicity factor for dark matter

REAL (DP), PARAMETER :: gs1 = 2.0E0_DP                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiplicity factor for normal matter

REAL (DP), PARAMETER :: gs2 = 2.0E0_DP                          
