!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine is for the initial check of the 
! Parameter.h input
! This helps you figure out if you have correctly
! input all the parameters in this file
!
! Written by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ReportParameter
USE DEFINITION
IMPLICIT NONE

WRITE(*,*) 'In Report Parameter'
WRITE(*,*)
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'The parameters in parameter.h are listed below:'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 1:'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,I10)") 'DP', DP
WRITE(*,"(A30,I10)") 'Even', even
WRITE(*,"(A30,I10)") 'oddR', oddR
WRITE(*,"(A30,I10)") 'oddZ', oddZ
WRITE(*,"(A30,I10)") 'oddB', oddB
WRITE(*,"(A30,ES13.5)") 'pi', pi
WRITE(*,"(A30,ES13.5)") 'hbar', hbar
WRITE(*,*)

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 2:'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,I10)") 'coordinate flag', coordinate_flag
WRITE(*,"(A30,I10)") 'hemisphere flag', hemisphere_flag
WRITE(*,"(A30,I10)") 'length_step_r_1', length_step_r_1
WRITE(*,"(A30,I10)") 'length_step_z_1', length_step_z_1
WRITE(*,"(A30,I10)") 'length_step_r_2', length_step_r_2
WRITE(*,"(A30,I10)") 'length_step_z_2', length_step_z_2
WRITE(*,"(A30,I10)") 'length_step_r_part_1', length_step_r_part_1
WRITE(*,"(A30,I10)") 'length_step_r_part_2', length_step_r_part_2
WRITE(*,"(A30,I10)") 'length_step_z_part_1', length_step_z_part_1
WRITE(*,"(A30,I10)") 'length_step_z_part_2', length_step_z_part_2
WRITE(*,"(A30,I10)") 'length_step_z_min_part_1', length_step_z_min_part_1
WRITE(*,"(A30,I10)") 'length_step_z_min_part_2', length_step_z_min_part_2
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'dx1_ini1', dx1_ini
WRITE(*,"(A30,ES13.5)") 'dx1', dx1
WRITE(*,"(A30,ES13.5)") 'dx2_ini', dx2_ini
WRITE(*,"(A30,ES13.5)") 'dx2', dx2
WRITE(*,*)
WRITE(*,"(A30,I10)") 'ini_acc', ini_acc
WRITE(*,"(A30,ES13.5)") 'smallpara', smallpara
WRITE(*,"(A30,ES13.5)") 'dxmore', dxmore
WRITE(*,"(A30,I10)") 'time_step', time_step
WRITE(*,"(A30,ES13.5)") 'cfl', cfl
WRITE(*,"(A30,ES13.5)") 'total_time', total_time
WRITE(*,"(A30,ES13.5)") 'dt', dt
WRITE(*,*)

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 3:'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'rho1_c', rho1_c
WRITE(*,"(A30,ES13.5)") 'rho2_c', rho2_c
WRITE(*,"(A30,ES13.5)") 'omega_ini', omega_ini
WRITE(*,"(A30,ES13.5)") 'temp_ini', temp_ini
WRITE(*,"(A30,ES13.5)") 'xhe4_ini1', xhe4_ini1
WRITE(*,"(A30,ES13.5)") 'xc12_ini1', xc12_ini1
WRITE(*,"(A30,ES13.5)") 'xo16_ini1', xo16_ini1
WRITE(*,"(A30,ES13.5)") 'xne20_ini1', xne20_ini1
WRITE(*,"(A30,ES13.5)") 'ye_ini1', ye_ini1
WRITE(*,*)
WRITE(*,"(A30,I10)") 'initmodel_layer2_flag', initmodel_layer2_flag
WRITE(*,"(A30,ES13.5)") 'xhe4_ini2', xhe4_ini2
WRITE(*,"(A30,ES13.5)") 'xc12_ini2', xc12_ini2
WRITE(*,"(A30,ES13.5)") 'xo16_ini2', xo16_ini2
WRITE(*,"(A30,ES13.5)") 'xne20_ini2', xne20_ini2
WRITE(*,"(A30,ES13.5)") 'ye_ini2', ye_ini2
WRITE(*,*)
WRITE(*,"(A30,I10)") 'initmodel_layer3_flag', initmodel_layer3_flag
WRITE(*,"(A30,ES13.5)") 'xhe4_ini3', xhe4_ini3
WRITE(*,"(A30,ES13.5)") 'xc12_ini3', xc12_ini3
WRITE(*,"(A30,ES13.5)") 'xo16_ini3', xo16_ini3 
WRITE(*,"(A30,ES13.5)") 'xne20_ini3', xne20_ini3
WRITE(*,"(A30,ES13.5)") 'ye_ini3', ye_ini3
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'rho1_a', rho1_a
WRITE(*,"(A30,ES13.5)") 'rho2_a', rho2_a
WRITE(*,"(A30,ES13.5)") 'vel1_a', vel1_a
WRITE(*,"(A30,ES13.5)") 'vel2_a', vel2_a
WRITE(*,"(A30,ES13.5)") 'temp_a', temp_a
WRITE(*,"(A30,ES13.5)") 'ye_a', ye_a
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'temp_max', temp_max
WRITE(*,"(A30,ES13.5)") 'temp_min', temp_min
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'ye_max', ye_max
WRITE(*,"(A30,ES13.5)") 'ye_min', ye_min
WRITE(*,*)

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 4'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,I10)") 'w_gravity_i', w_gravity_i
WRITE(*,"(A30,I10)") 'relax_max', relax_max
WRITE(*,"(A30,ES13.5)") 'tolerance', tolerance
WRITE(*,"(A30,ES13.5)") 'sor_weight', sor_weight
WRITE(*,*)

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 5'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,I10)") 'no_of_eq (tent.)', no_of_eq
WRITE(*,*)
WRITE(*,"(A30,I10)") 'DM_flag', DM_flag
WRITE(*,"(A30,I10)") 'runDM_flag', runDM_flag
WRITE(*,"(A30,I10)") 'rotationdm_flag', rotationdm_flag
WRITE(*,"(A30,I10)") 'rotationnm_flag', rotationnm_flag
WRITE(*,"(A30,I10)") 'fixrhodm_flag', fixrhodm_flag
WRITE(*,"(A30,I10)") 'fixrhonm_flag', fixrhonm_flag
WRITE(*,"(A30,I10)") 'found_movinggriddm_flag', found_movinggriddm_flag
WRITE(*,"(A30,I10)") 'found_movinggridnm_flag', found_movinggridnm_flag
WRITE(*,"(A30,I10)") 'movinggriddm_flag', movinggriddm_flag
WRITE(*,"(A30,I10)") 'movinggridnm_flag', movinggridnm_flag
WRITE(*,"(A30,I10)") 'initmodel_flag', initmodel_flag
WRITE(*,"(A30,I10)") 'output_flag', output_flag
WRITE(*,*)
WRITE(*,"(A30,I10)") 'Boundary flag 1', boundary_flag(1)
WRITE(*,"(A30,I10)") 'Boundary flag 2', boundary_flag(2)
WRITE(*,"(A30,I10)") 'Boundary flag 3', boundary_flag(3)
WRITE(*,"(A30,I10)") 'Boundary flag 4', boundary_flag(4)
WRITE(*,*)
WRITE(*,"(A30,I10)") 'Checkrho_flag', checkrho_flag
WRITE(*,"(A30,I10)") 'checkstepdm_flag', checkstepdm_flag
WRITE(*,"(A30,I10)") 'checkstepnm_flag', checkstepnm_flag
WRITE(*,"(A30,I10)") 'Checkvel_flag', checkvel_flag
WRITE(*,"(A30,I10)") 'Updatedt_flag', updatedt_flag
WRITE(*,"(A30,I10)") 'testmodel_flag', testmodel_flag
WRITE(*,*)
WRITE(*,"(A30,I10)") 'polyeosdm_flag', polyeosdm_flag
WRITE(*,"(A30,I10)") 'fermieosdm_flag', fermieosdm_flag
WRITE(*,"(A30,I10)") 'polyeosnm_flag', polyeosnm_flag
WRITE(*,"(A30,I10)") 'fermieosnm_flag', fermieosnm_flag
WRITE(*,"(A30,I10)") 'helmeos_flag', helmeos_flag
WRITE(*,"(A30,I10)") 'xisotran_flag', xisotran_flag
WRITE(*,"(A30,I10)") 'burn_prog_flag', burn_prog_flag
WRITE(*,"(A30,I10)") 'etran_flag', etran_flag
WRITE(*,"(A30,I10)") 'ecap_flag', ecap_flag
WRITE(*,"(A30,ES13.5)") 'k_1', k_1
WRITE(*,"(A30,ES13.5)") 'k_2', k_2
WRITE(*,"(A30,ES13.5)") 'gamma1', gamma1
WRITE(*,"(A30,ES13.5)") 'gamma2', gamma2
WRITE(*,*)
WRITE(*,"(A30,I10)") 'fusion_flag', fusion_flag
WRITE(*,"(A30,I10)") 'deflevelset_flag', deflevelset_flag
WRITE(*,"(A30,I10)") 'detlevelset_flag', detlevelset_flag
WRITE(*,"(A30,I10)") 'flame_flag', flame_flag
WRITE(*,"(A30,I10)") 'deton_flag', deton_flag
WRITE(*,"(A30,I10)") 'carburn_flag', carburn_flag
WRITE(*,"(A30,I10)") 'advburn_flag', advburn_flag
WRITE(*,"(A30,I10)") 'convert_nse_flag', convert_nse_flag
WRITE(*,"(A30,I10)") 'turb_flag', turb_flag
WRITE(*,"(A30,I10)") 'nuspec_flag', nuspec_flag
WRITE(*,"(A30,I10)") 'tracer_flag', tracer_flag
WRITE(*,"(A30,I10)") 'gravwave_flag', gravwave_flag
WRITE(*,"(A30,I10)") 'sponge_flag', sponge_flag
WRITE(*,*)

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 6'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'Output_logtime', output_logtime
WRITE(*,"(A30,ES13.5)") 'Output_logtime_last', output_logtime_last
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'Output_profiletime', output_profiletime
WRITE(*,"(A30,ES13.5)") 'Output_profiletime_last', output_profiletime_last
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'Output_flametime', output_flametime  
WRITE(*,"(A30,ES13.5)") 'Output_flametime_last', output_flametime_last
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'Output_turbtime', output_turbtime  
WRITE(*,"(A30,ES13.5)") 'Output_turbtime_last', output_turbtime_last
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'Output_helmtime', output_helmtime  
WRITE(*,"(A30,ES13.5)") 'Output_helmtime_last', output_helmtime_last
WRITE(*,*)
WRITE(*,"(A30,ES13.5)") 'Output_PPTtime', output_PPTtime  
WRITE(*,"(A30,ES13.5)") 'Output_PPTtime_last', output_PPTtime_last
WRITE(*,*) 

WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*) 'In Section 8'
WRITE(*,*) '-----------------------------------------------------------------'
WRITE(*,*)
WRITE(*,"(A30,I10)") 'masscut_flag', masscut_flag
WRITE(*,"(A30,ES13.5)") 'rad_cut', rad_cut
WRITE(*,"(A30,ES13.5)") 'mass_cut', mass_cut
WRITE(*,*)

WRITE(*,"(A30,I10)") 'jetexp_flag', jetexp_flag
WRITE(*,"(A30,ES13.5)") 'Dedep', dedep
WRITE(*,"(A30,ES13.5)") 'theta_jet', theta_jet
WRITE(*,"(A30,ES13.5)") 'fth', fth
WRITE(*,"(A30,ES13.5)") 'jet_time', jet_time
WRITE(*,"(A30,ES13.5)") 'adia_jet', adia_jet
WRITE(*,"(A30,ES13.5)") 'gamma_jet', gamma_jet
WRITE(*,"(A30,ES13.5)") 'kap_e', kap_e
WRITE(*,*)

WRITE(*,*) 'Please check if all parameters are correctly input'
WRITE(*,*) 'Press enter to continue (The read command can muted for automatic run'
!READ(*,*)

WRITE(*,*) 'Finisied report parameter'

END SUBROUTINE ReportParameter