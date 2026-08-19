import subprocess
import sys
import torch
import random
import numpy as np
import os
import multiprocessing


class CFG:
    
    have_sig_events_1 = 2.162835 * 1e6  # test
    have_bg_events_1 = 39.54 * 1e6
    have_sig_events_2 = 2.16 * 1e6      # train
    have_bg_events_2 = 39.909958 * 1e6
    mass_interval = (2.24763, 2.32497)
    
    # Physical constants
    sigma_MB = 38 * 1e-3 * 1e-24                # Minimum bias cross-section in cm2
    cross_section_Lc_plus = 4 * 1e-6 * 1e-24    # Lc+(4122) cross-section in cm2
    L = 1e32                                    # Luminosity in cm-2s-1
    t = 31536000                                # Time period in seconds (1 year)
    branching = 0.035                           # Lc decay branching ratio
    
    # Calculate expected events number
    real_sig_events = L * t * cross_section_Lc_plus * branching
    real_bg_events = L * t * sigma_MB
    
    m_proton = 0.938272
    m_kaon = 0.493677
    m_pion = 0.139570
    
    # Test data
    sig_eff_presel = 0.9471810669404919
    bg_eff_presel = 0.6200400630780377
    
    seed = 58
    
    #######################################################################################
    # GPU
    gpu_available = torch.cuda.is_available()

    print(f"CUDA available: {gpu_available}")
    if gpu_available:
        print(f"GPU name: {torch.cuda.get_device_name(0)}")
        print(f"Number of GPUs: {torch.cuda.device_count()}")
    else:
        print(f'Use CPU, \nNumber of CPUs: {multiprocessing.cpu_count()}')
    
    #######################################################################################
    # Seed
    
    @staticmethod
    def seed_all(seed=42):
        random.seed(seed)
        np.random.seed(seed)
        os.environ['PYTHONHASHSEED'] = str(seed)

        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        
    #######################################################################################
    # Features sets
    
    target_name='tag'
    
    bad_features = [
        'mass_Lc', 'mass_significance', 'mass_pt_correlation', 'Mt_Lc',
        'multiplicity', 'is_mass_window', 'mass_pt_product'
    ]
    
    continual_features = [
        # Initial kinematic
        'mass_Lc',
        'P_p', 'P_pip', 'P_K', 'P_Lc',
        'eta_p', 'eta_pip', 'eta_K', 'eta_Lc',
        'Pt_p', 'Pt_K', 'Pt_pip', 'Pt_Lc',
        'OA_p', 'OA_K', 'OA_pip',
        'xF',
        'phi_angle',
        
        'probab_p', 'probab_k', 'probab_pip',
        
        # Initial topology
        'lengthXY_Lc', 'dlengthXY_Lc', 'ctau_Lc', 'dctau_Lc',
        'chi2_Lc_PV_xy', 'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'chi2_p_PV_xy', 'dist_p_PV_xy', 'dist_p_PV_xy_custom',
        'chi2_K_PV_xy', 'dist_K_PV_xy', 'dist_K_PV_xy_custom',
        'chi2_pip_PV_xy', 'dist_pip_PV_xy', 'dist_pip_PV_xy_custom',
        'chi2_p_Lc_xy', 'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'chi2_K_Lc_xy', 'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'chi2_pip_Lc_xy', 'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'chi2_K_pip_xy', 'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'chi2_p_K_xy', 'dist_p_K_xy', 'dist_p_K_xy_custom',
        'chi2_p_pip_xy', 'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'cosAngle_r_Lc_momentum_Lc_xy',
        'cosAngle_r_Lc_sum_momentum_xy',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
        
        # Added kinematic
        'mass_significance', 'mass_pt_correlation', 'Mt_Lc',
        'p_pv_significance', 'K_pv_significance', 'pip_pv_significance',
        'p_momentum_ratio', 'pip_momentum_ratio', 'K_momentum_ratio',
        'p_Pt_fraction', 'pip_Pt_fraction', 'K_Pt_fraction',
        'delta_eta_p', 'delta_eta_pip', 'delta_eta_K',
        'min_pt_daughter', 'max_pt_daughter', 'pt_range',
        'min_OA_daughter', 'max_OA_daughter', 'OA_range',
        'E_p', 'E_pip', 'E_K', 'E_total', 'pt_over_E',
        'Mt_p', 'Mt_K', 'Mt_pip',
        'opening_angle_sum', 'opening_angle_product', 'lifetime_opening_interaction',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',
        'Pt_p_over_Pt_K', 'Pt_p_over_Pt_pip', 'Pt_K_over_Pt_pip',
        'P_p_over_P_K', 'P_p_over_P_pip', 'P_K_over_P_pip',
        'sum_squares_pt', 'sum_squares_p', 'sum_squares_OA',
        'sum_abs_pt', 'sum_abs_p', 'sum_abs_OA',
        'pt_norm_ratio', 'p_norm_ratio',
        'momentum_balance',
        'OA_p_ratio', 'OA_K_ratio', 'OA_pip_ratio',
        'pt_over_ctau', 'pt_over_length',
        'triple_product',
        'energy_asymmetry_pK', 'energy_asymmetry_ppip', 'energy_asymmetry_Kpip', 'energy_asymmetry_p_total',
        'pt_asymmetry_p_K', 'pt_asymmetry_p_pip', 'pt_asymmetry_K_pip',
        'asymmetry_p_K', 'asymmetry_p_pip', 'asymmetry_K_pip',
        'OA_asymmetry_p_K', 'OA_asymmetry_p_pip', 'OA_asymmetry_K_pip',
        'pt_cv', 'OA_cv',
        'pt_eta_correlation_p', 'pt_eta_correlation_K', 'pt_eta_correlation_pip',
        'pt_p_times_OA_p', 'pt_K_times_OA_K', 'pt_pip_times_OA_pip',
        'ctau_pt_product', 'mass_pt_product', 'xF_pt_product',
        
        # Added topology
        'l_over_dl_XY', 'l_over_dl_pt_product',
        'p_displacement_ratio', 'K_displacement_ratio', 'pip_displacement_ratio',
        'mean_chi2_to_PV', 'mean_chi2_to_Lc',
        'topology_score',
        'Lc_dca_significance',
        'p_PV_dca_significance', 'K_PV_dca_significance', 'pip_PV_dca_significance',
        'p_Lc_dca_significance', 'K_Lc_dca_significance', 'pip_Lc_dca_significance',
        'max_dca_Lc_xy', 'min_dca_Lc_xy',
        'max_dca_PV_xy', 'min_dca_PV_xy',
        'max_dca_Lc_xy_custom', 'min_dca_Lc_xy_custom',
        'max_dca_PV_xy_custom', 'min_dca_PV_xy_custom',
        'dca_sum_consistency', 'dca_sum_consistency_custom',
        'mean_dca', 'mean_dca_custom',
        'min_chi2_daughter', 'max_chi2_daughter', 'chi2_range',
        'alignment_score',
        'sum_squares_chi2_Lc_xy', 'sum_squares_dca_Lc_xy', 'sum_squares_dca_custom_Lc_xy',
        'sum_squares_chi2_PV_xy', 'sum_squares_dca_PV_xy', 'sum_squares_dca_custom_PV_xy',
        'sum_abs_chi2_Lc_xy', 'sum_abs_dca_Lc_xy', 'sum_abs_dca_custom_Lc_xy',
        'sum_abs_chi2_PV_xy', 'sum_abs_dca_PV_xy', 'sum_abs_dca_custom_PV_xy',
        'vertex_consistency', 'chi2_consistency',
    ]

    nominal_features = [
        'highest_pt_daughter', 'highest_momentum_daughter', 'highest_OA_daughter',
        'highest_eta_daughter', 'highest_dca_Lc_daughter', 'highest_dca_PV_daughter',
    ]

    bunary_features = [
        'is_high_pt', 'is_long_lived', 'is_well_separated',
    ]
    
    cols_to_log_transform = [
        'P_p', 'P_K', 'P_pip', 'Pt_p', 'Pt_K', 'Pt_pip',
        'dist_p_PV_xy_custom', 'dist_K_PV_xy_custom', 'dist_pip_PV_xy_custom',
        'dist_p_Lc_xy_custom', 'dist_K_Lc_xy_custom', 'dist_pip_Lc_xy_custom',
        'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'dist_p_K_xy', 'dist_p_K_xy_custom',
        'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'P_Lc', 'Pt_Lc',
        'dist_Lc_PV_xy_custom',
        'dlengthXY_Lc', 'chi2_Lc',
        'p_pv_significance', 'K_pv_significance', 'pip_pv_significance',
        'p_momentum_ratio', 'K_momentum_ratio', 'pip_momentum_ratio',
        'p_Pt_fraction', 'K_Pt_fraction', 'pip_Pt_fraction',
        'delta_eta_p', 'delta_eta_K', 'delta_eta_pip',
        'E_p', 'E_K', 'E_pip',
        'E_total', 'Mt_p', 'Mt_K', 'Mt_pip',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',
        'P_p_over_P_K', 'P_p_over_P_pip', 'P_K_over_P_pip',
        'Pt_p_over_Pt_K', 'Pt_p_over_Pt_pip', 'Pt_K_over_Pt_pip',
        'sum_squares_p', 'sum_squares_pt', 'sum_squares_OA',
        'momentum_balance', 'OA_p_ratio', 'OA_K_ratio', 'OA_pip_ratio',
        'pt_p_times_OA_p', 'pt_K_times_OA_K', 'pt_pip_times_OA_pip',
        'triple_product',
        'min_pt_daughter', 'max_pt_daughter', 'pt_range',
        'max_chi2_daughter', 'chi2_range', 'pt_norm_ratio',
        'sum_squares_chi2_Lc_xy', 'sum_squares_dca_Lc_xy', 'sum_squares_dca_custom_Lc_xy',
        'sum_squares_chi2_PV_xy', 'sum_squares_dca_PV_xy', 'sum_squares_dca_custom_PV_xy',
        'p_PV_dca_significance', 'K_PV_dca_significance', 'pip_PV_dca_significance',
        'p_Lc_dca_significance', 'K_Lc_dca_significance', 'pip_Lc_dca_significance',
        'max_dca_Lc_xy', 'max_dca_PV_xy',
        'max_dca_Lc_xy_custom', 'min_dca_Lc_xy_custom',
        'max_dca_PV_xy_custom', 'min_dca_PV_xy_custom',
        'mean_dca', 'mean_dca_custom',
    ]
    
    feature_set_0 = [
        # Initial kinematic
        'P_p', 'P_pip', 'P_K', 'P_Lc',
        'eta_p', 'eta_pip', 'eta_K', 'eta_Lc',
        'Pt_p', 'Pt_K', 'Pt_pip', 'Pt_Lc',
        'OA_p', 'OA_K', 'OA_pip',
        'xF',
        
        'probab_p', 'probab_k', 'probab_pip',
        
        # Initial topology
        'lengthXY_Lc', 'dlengthXY_Lc', 'ctau_Lc', 'dctau_Lc',
        'chi2_Lc_PV_xy', 'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'chi2_p_PV_xy', 'dist_p_PV_xy', 'dist_p_PV_xy_custom',
        'chi2_K_PV_xy', 'dist_K_PV_xy', 'dist_K_PV_xy_custom',
        'chi2_pip_PV_xy', 'dist_pip_PV_xy', 'dist_pip_PV_xy_custom',
        'chi2_p_Lc_xy', 'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'chi2_K_Lc_xy', 'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'chi2_pip_Lc_xy', 'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'chi2_K_pip_xy', 'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'chi2_p_K_xy', 'dist_p_K_xy', 'dist_p_K_xy_custom',
        'chi2_p_pip_xy', 'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'cosAngle_r_Lc_momentum_Lc_xy',
        'cosAngle_r_Lc_sum_momentum_xy',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
    ]
    
    feature_set_1 = [
        
        'probab_p', 'probab_k', 'probab_pip',
        
        'is_high_pt', 'is_long_lived', 'is_well_separated',
        
        'highest_pt_daughter', 'highest_momentum_daughter', 'highest_OA_daughter',
        'highest_eta_daughter', 'highest_dca_Lc_daughter', 'highest_dca_PV_daughter',
        
        # Initial kinematic
        'P_p', 'P_pip', 'P_K', 'P_Lc',
        'eta_p', 'eta_pip', 'eta_K', 'eta_Lc',
        'Pt_p', 'Pt_K', 'Pt_pip', 'Pt_Lc',
        'OA_p', 'OA_K', 'OA_pip',
        'xF',
        
        # Initial topology
        'lengthXY_Lc', 'dlengthXY_Lc', 'ctau_Lc', 'dctau_Lc',
        'chi2_Lc_PV_xy', 'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'chi2_p_PV_xy', 'dist_p_PV_xy', 'dist_p_PV_xy_custom',
        'chi2_K_PV_xy', 'dist_K_PV_xy', 'dist_K_PV_xy_custom',
        'chi2_pip_PV_xy', 'dist_pip_PV_xy', 'dist_pip_PV_xy_custom',
        'chi2_p_Lc_xy', 'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'chi2_K_Lc_xy', 'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'chi2_pip_Lc_xy', 'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'chi2_K_pip_xy', 'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'chi2_p_K_xy', 'dist_p_K_xy', 'dist_p_K_xy_custom',
        'chi2_p_pip_xy', 'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'cosAngle_r_Lc_momentum_Lc_xy',
        'cosAngle_r_Lc_sum_momentum_xy',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
        
        # Added kinematic
        'p_pv_significance', 'K_pv_significance', 'pip_pv_significance',
        'p_momentum_ratio', 'pip_momentum_ratio', 'K_momentum_ratio',
        'p_Pt_fraction', 'pip_Pt_fraction', 'K_Pt_fraction',
        'delta_eta_p', 'delta_eta_pip', 'delta_eta_K',
        'min_pt_daughter', 'max_pt_daughter', 'pt_range',
        'min_OA_daughter', 'max_OA_daughter', 'OA_range',
        'E_p', 'E_pip', 'E_K', 'E_total', 'pt_over_E',
        'Mt_p', 'Mt_K', 'Mt_pip',
        'opening_angle_sum', 'opening_angle_product', 'lifetime_opening_interaction',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',
        'Pt_p_over_Pt_K', 'Pt_p_over_Pt_pip', 'Pt_K_over_Pt_pip',
        'P_p_over_P_K', 'P_p_over_P_pip', 'P_K_over_P_pip',
        'sum_squares_pt', 'sum_squares_p', 'sum_squares_OA',
        'sum_abs_pt', 'sum_abs_p', 'sum_abs_OA',
        'pt_norm_ratio', 'p_norm_ratio',
        'momentum_balance',
        'OA_p_ratio', 'OA_K_ratio', 'OA_pip_ratio',
        'pt_over_ctau', 'pt_over_length',
        'triple_product',
        'energy_asymmetry_pK', 'energy_asymmetry_ppip', 'energy_asymmetry_Kpip', 'energy_asymmetry_p_total',
        'pt_asymmetry_p_K', 'pt_asymmetry_p_pip', 'pt_asymmetry_K_pip',
        'asymmetry_p_K', 'asymmetry_p_pip', 'asymmetry_K_pip',
        'OA_asymmetry_p_K', 'OA_asymmetry_p_pip', 'OA_asymmetry_K_pip',
        'pt_cv', 'OA_cv',
        'pt_eta_correlation_p', 'pt_eta_correlation_K', 'pt_eta_correlation_pip',
        'pt_p_times_OA_p', 'pt_K_times_OA_K', 'pt_pip_times_OA_pip',
        'ctau_pt_product', 'xF_pt_product',
        
        # Added topology
        'l_over_dl_XY', 'l_over_dl_pt_product',
        'p_displacement_ratio', 'K_displacement_ratio', 'pip_displacement_ratio',
        'mean_chi2_to_PV', 'mean_chi2_to_Lc',
        'topology_score',
        'Lc_dca_significance',
        'p_PV_dca_significance', 'K_PV_dca_significance', 'pip_PV_dca_significance',
        'p_Lc_dca_significance', 'K_Lc_dca_significance', 'pip_Lc_dca_significance',
        'max_dca_Lc_xy', 'min_dca_Lc_xy',
        'max_dca_PV_xy', 'min_dca_PV_xy',
        'max_dca_Lc_xy_custom', 'min_dca_Lc_xy_custom',
        'max_dca_PV_xy_custom', 'min_dca_PV_xy_custom',
        'dca_sum_consistency', 'dca_sum_consistency_custom',
        'mean_dca', 'mean_dca_custom',
        'min_chi2_daughter', 'max_chi2_daughter', 'chi2_range',
        'alignment_score',
        'sum_squares_chi2_Lc_xy', 'sum_squares_dca_Lc_xy', 'sum_squares_dca_custom_Lc_xy',
        'sum_squares_chi2_PV_xy', 'sum_squares_dca_PV_xy', 'sum_squares_dca_custom_PV_xy',
        'sum_abs_chi2_Lc_xy', 'sum_abs_dca_Lc_xy', 'sum_abs_dca_custom_Lc_xy',
        'sum_abs_chi2_PV_xy', 'sum_abs_dca_PV_xy', 'sum_abs_dca_custom_PV_xy',
        'vertex_consistency', 'chi2_consistency',
    ]
        
    # Corr < 0.85, no bad features, ranking by L2 weights + f_regressor
    feature_set_2 = [
        'probab_p', 'probab_k', 'probab_pip',
        
        # Initial kinematic
        'P_p', 'P_pip', 'P_K', 'P_Lc',
        'eta_p', 'eta_pip', 'eta_K', 'eta_Lc',
        'Pt_p', 'Pt_K', 'Pt_pip', 'Pt_Lc',
        'OA_p', 'OA_K', 'OA_pip',
        'xF',
        
        # Initial topology
        'lengthXY_Lc', 'dlengthXY_Lc', 'ctau_Lc', 'dctau_Lc',
        'chi2_Lc_PV_xy', 'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'chi2_p_PV_xy', 'dist_p_PV_xy', 'dist_p_PV_xy_custom',
        'chi2_K_PV_xy', 'dist_K_PV_xy', 'dist_K_PV_xy_custom',
        'chi2_pip_PV_xy', 'dist_pip_PV_xy', 'dist_pip_PV_xy_custom',
        'chi2_p_Lc_xy', 'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'chi2_K_Lc_xy', 'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'chi2_pip_Lc_xy', 'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'chi2_K_pip_xy', 'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'chi2_p_K_xy', 'dist_p_K_xy', 'dist_p_K_xy_custom',
        'chi2_p_pip_xy', 'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'cosAngle_r_Lc_momentum_Lc_xy',
        'cosAngle_r_Lc_sum_momentum_xy',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
        
        # Added kinematic
        'p_pv_significance', 'K_pv_significance', 'pip_pv_significance',
        'p_momentum_ratio', 'pip_momentum_ratio', 'K_momentum_ratio',
        'p_Pt_fraction', 'pip_Pt_fraction', 'K_Pt_fraction',
        'delta_eta_p', 'delta_eta_pip', 'delta_eta_K',
        'min_pt_daughter', 'max_pt_daughter', 'pt_range',
        'min_OA_daughter', 'max_OA_daughter', 'OA_range',
        'E_p', 'E_pip', 'E_K', 'E_total', 'pt_over_E',
        'Mt_p', 'Mt_K', 'Mt_pip',
        'opening_angle_sum', 'opening_angle_product', 'lifetime_opening_interaction',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',
        'Pt_p_over_Pt_K', 'Pt_p_over_Pt_pip', 'Pt_K_over_Pt_pip',
        'P_p_over_P_K', 'P_p_over_P_pip', 'P_K_over_P_pip',
        'sum_squares_pt', 'sum_squares_p', 'sum_squares_OA',
        'sum_abs_pt', 'sum_abs_p', 'sum_abs_OA',
        'pt_norm_ratio', 'p_norm_ratio',
        'momentum_balance',
        'OA_p_ratio', 'OA_K_ratio', 'OA_pip_ratio',
        'pt_over_ctau', 'pt_over_length',
        'triple_product',
        'energy_asymmetry_pK', 'energy_asymmetry_ppip', 'energy_asymmetry_Kpip', 'energy_asymmetry_p_total',
        'pt_asymmetry_p_K', 'pt_asymmetry_p_pip', 'pt_asymmetry_K_pip',
        'asymmetry_p_K', 'asymmetry_p_pip', 'asymmetry_K_pip',
        'OA_asymmetry_p_K', 'OA_asymmetry_p_pip', 'OA_asymmetry_K_pip',
        'pt_cv', 'OA_cv',
        'pt_eta_correlation_p', 'pt_eta_correlation_K', 'pt_eta_correlation_pip',
        'pt_p_times_OA_p', 'pt_K_times_OA_K', 'pt_pip_times_OA_pip',
        'ctau_pt_product', 'xF_pt_product',
        
        # Added topology
        'l_over_dl_XY', 'l_over_dl_pt_product',
        'p_displacement_ratio', 'K_displacement_ratio', 'pip_displacement_ratio',
        'mean_chi2_to_PV', 'mean_chi2_to_Lc',
        'topology_score',
        'Lc_dca_significance',
        'p_PV_dca_significance', 'K_PV_dca_significance', 'pip_PV_dca_significance',
        'p_Lc_dca_significance', 'K_Lc_dca_significance', 'pip_Lc_dca_significance',
        'max_dca_Lc_xy', 'min_dca_Lc_xy',
        'max_dca_PV_xy', 'min_dca_PV_xy',
        'max_dca_Lc_xy_custom', 'min_dca_Lc_xy_custom',
        'max_dca_PV_xy_custom', 'min_dca_PV_xy_custom',
        'dca_sum_consistency', 'dca_sum_consistency_custom',
        'mean_dca', 'mean_dca_custom',
        'min_chi2_daughter', 'max_chi2_daughter', 'chi2_range',
        'alignment_score',
        'sum_squares_chi2_Lc_xy', 'sum_squares_dca_Lc_xy', 'sum_squares_dca_custom_Lc_xy',
        'sum_squares_chi2_PV_xy', 'sum_squares_dca_PV_xy', 'sum_squares_dca_custom_PV_xy',
        'sum_abs_chi2_Lc_xy', 'sum_abs_dca_Lc_xy', 'sum_abs_dca_custom_Lc_xy',
        'sum_abs_chi2_PV_xy', 'sum_abs_dca_PV_xy', 'sum_abs_dca_custom_PV_xy',
        'vertex_consistency', 'chi2_consistency',
    ]    
    
    # RFE XGBoost
    feature_set_3 = [
        'highest_momentum_daughter',
        'probab_p', 'probab_k', 'probab_pip',
        'P_p', 'P_pip', 'P_K',
        'eta_p',
        'Pt_p', 'Pt_K', 'Pt_pip',
        'OA_K',
        'cosAngle_r_Lc_sum_momentum_xy',
        'pip_momentum_ratio', 'K_momentum_ratio',
        'min_pt_daughter',
        'opening_angle_sum',
        'lifetime_opening_interaction',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',
        'Pt_p_over_Pt_pip', 'P_p_over_P_K', 'P_K_over_P_pip',
        'sum_squares_OA',
        'sum_abs_pt',
        'OA_K_ratio', 'OA_pip_ratio',
        'energy_asymmetry_pK', 'energy_asymmetry_Kpip',
        'OA_asymmetry_p_K', 'OA_asymmetry_K_pip',
        'l_over_dl_XY',
        'l_over_dl_pt_product',
    ]
    
    # Check also:
    # 'eta_p', 'eta_pip', 'eta_K', 'eta_Lc',
    # 'is_high_pt', 'is_long_lived', 'is_well_separated',
    feature_cut_based_sel = [
        
        # Initial kinematic
        'Pt_p', 'Pt_K', 'Pt_pip', 'Pt_Lc',
        'OA_p', 'OA_K', 'OA_pip',
        
        # Initial topology
        'ctau_Lc',
        'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'dist_p_PV_xy', 'dist_p_PV_xy_custom',
        'dist_K_PV_xy', 'dist_K_PV_xy_custom',
        'dist_pip_PV_xy', 'dist_pip_PV_xy_custom',
        'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'dist_p_K_xy', 'dist_p_K_xy_custom',
        'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'cosAngle_r_Lc_momentum_Lc_xy',
        'cosAngle_r_Lc_sum_momentum_xy',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
        
        # Added kinematic
        'sum_abs_pt',
        'p_momentum_ratio', 'pip_momentum_ratio', 'K_momentum_ratio',
        'p_Pt_fraction', 'pip_Pt_fraction', 'K_Pt_fraction',
        'delta_eta_p', 'delta_eta_pip', 'delta_eta_K',
        'min_pt_daughter', 'max_pt_daughter', 'pt_range',
        'min_OA_daughter', 'max_OA_daughter', 'OA_range',
        'E_p', 'E_pip', 'E_K', 'E_total', 'pt_over_E',
        'Mt_p', 'Mt_K', 'Mt_pip',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',
        'momentum_balance',
        
        # Added topology
        'p_pv_significance', 'K_pv_significance', 'pip_pv_significance',
        'l_over_dl_XY',
        'mean_chi2_to_PV', 'mean_chi2_to_Lc',
        'topology_score',
        'Lc_dca_significance',
        'p_PV_dca_significance', 'K_PV_dca_significance', 'pip_PV_dca_significance',
        'p_Lc_dca_significance', 'K_Lc_dca_significance', 'pip_Lc_dca_significance',
        'max_dca_Lc_xy', 'min_dca_Lc_xy',
        'max_dca_PV_xy', 'min_dca_PV_xy',
        'max_dca_Lc_xy_custom', 'min_dca_Lc_xy_custom',
        'max_dca_PV_xy_custom', 'min_dca_PV_xy_custom',
        'mean_dca', 'mean_dca_custom',
        'min_chi2_daughter', 'max_chi2_daughter',
        'alignment_score',
    ]
    
    cut_based_sel_restriction = {
        
        'P_p': 'right', 'P_pip': 'right', 'P_K': 'right', 'P_Lc': 'right',
        'Pt_p': 'right', 'Pt_K': 'right', 'Pt_pip': 'right', 'Pt_Lc': 'right',
        
        'lengthXY_Lc': 'right', 'dlengthXY_Lc': 'right', 'ctau_Lc': 'right',
        'dist_Lc_PV_xy': 'left', 'dist_Lc_PV_xy_custom': 'left',
        'dist_p_PV_xy_custom': 'left',
        'dist_K_PV_xy_custom': 'left',
        'dist_p_Lc_xy_custom': 'left',
        'dist_K_Lc_xy_custom': 'left',
        'dist_pip_Lc_xy_custom': 'left',
        'dist_K_pip_xy': 'left', 'dist_K_pip_xy_custom': 'left',
        'dist_p_K_xy': 'left', 'dist_p_K_xy_custom': 'left',
        'dist_p_pip_xy': 'left', 'dist_p_pip_xy_custom': 'left',
        'cosAngle_r_Lc_momentum_Lc_xy': 'right',
        'cosAngle_r_Lc_sum_momentum_xy': 'right',
        'cosAngle_momentum_Lc_sum_momentum_xy': 'right',
        'chi2_Lc': 'left',
        
        'sum_abs_pt': 'right',
        
        'p_pv_significance': 'left', 'K_pv_significance': 'left', 'pip_pv_significance': 'left',
        'delta_eta_p': 'left', 'delta_eta_pip': 'left', 'delta_eta_K': 'left',
        'min_pt_daughter': 'right', 'max_pt_daughter': 'right', 'pt_range': 'right',
        'pt_over_E': 'right',
        'Mt_p': 'right', 'Mt_K': 'right', 'Mt_pip': 'right',
        'sum_squares_pt': 'right',
        'triple_product': 'right',
        'pt_p_times_OA_p': 'right', 'pt_K_times_OA_K': 'right', 'pt_pip_times_OA_pip': 'right',
                
        'l_over_dl_XY': 'right',
        'p_PV_dca_significance': 'left', 'K_PV_dca_significance': 'left', 'pip_PV_dca_significance': 'left',
        'p_Lc_dca_significance': 'left', 'K_Lc_dca_significance': 'left', 'pip_Lc_dca_significance': 'left',
        'max_dca_Lc_xy_custom': 'left', 'min_dca_Lc_xy_custom': 'left',
        'max_dca_PV_xy_custom': 'left', 'min_dca_PV_xy_custom': 'left',
        'mean_dca': 'left', 'mean_dca_custom': 'left',
        'alignment_score': 'right',
    }
    
    # Custom selection
    feature_cut_based_sel_custom = [
        
        # Initial kinematic
        'Pt_Lc',
        'OA_p', 'OA_K', 'OA_pip',
        
        # Initial topology
        'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'dist_p_PV_xy', 'dist_p_PV_xy_custom',
        'dist_K_PV_xy', 'dist_K_PV_xy_custom',
        'dist_pip_PV_xy', 'dist_pip_PV_xy_custom',
        'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'dist_K_pip_xy', 'dist_K_pip_xy_custom',
        'dist_p_K_xy', 'dist_p_K_xy_custom',
        'dist_p_pip_xy', 'dist_p_pip_xy_custom',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
        
        # Added kinematic
        'sum_abs_pt',
        'delta_eta_p', 'delta_eta_pip', 'delta_eta_K',
        'rapidity_p', 'rapidity_K', 'rapidity_pip',

        
        # Added topology
        'l_over_dl_XY',
        'mean_chi2_to_Lc',
        'topology_score',
        'p_Lc_dca_significance',
        'max_dca_Lc_xy',
        'max_dca_PV_xy', 'min_dca_PV_xy',
        'min_dca_Lc_xy_custom',
        'mean_dca', 'mean_dca_custom',
        'min_chi2_daughter',
    ]
    
    # abs(corr) < 0.7,  shap val ranking
    feature_cut_based_sel_low_corr = [
        
        # Initial kinematic
        'Pt_K', 'Pt_pip',
        'OA_p', 'OA_K', 'OA_pip',
        
        # Initial topology
        'dist_Lc_PV_xy', 'dist_Lc_PV_xy_custom',
        'dist_p_PV_xy',
        'dist_K_PV_xy',
        'dist_pip_PV_xy',
        'dist_p_Lc_xy', 'dist_p_Lc_xy_custom',
        'dist_K_Lc_xy', 'dist_K_Lc_xy_custom',
        'dist_pip_Lc_xy', 'dist_pip_Lc_xy_custom',
        'cosAngle_momentum_Lc_sum_momentum_xy',
        'chi2_Lc',
        
        # Added kinematic
        'sum_abs_pt',
        'p_momentum_ratio', 'pip_momentum_ratio', 'K_momentum_ratio',
        'p_Pt_fraction', 'pip_Pt_fraction', 'K_Pt_fraction',
        'delta_eta_p', 'delta_eta_pip',
        'min_pt_daughter',
        'min_OA_daughter', 'max_OA_daughter',
        'E_pip', 'E_K',
        'rapidity_p', 'rapidity_pip',
        
        # Added topology
        'l_over_dl_XY',
        'mean_chi2_to_Lc',
        'topology_score',
        'p_Lc_dca_significance',
        'max_dca_Lc_xy',
        'max_dca_PV_xy', 'min_dca_PV_xy',
        'min_dca_Lc_xy_custom',
        'mean_dca', 'mean_dca_custom',
        'min_chi2_daughter',
    ]
    
    #######################################################################################
    # models params
    
    # XGBoost, features_set_3, Optuna
    params_xgboost_1 = {
        'early_stopping_rounds': 80,
        'n_estimators': 900,
        'max_depth': 9,
        'eta': 0.055739548129716016,
        'lambda': 1.023827972679177e-05,
        'alpha': 0.04117995219623043,
        'gamma': 0.7146754017014212,
        'min_child_weight': 1.8765916208440991,
        'max_delta_step': 7.076090336641161,
        'subsample': 0.8538320684738833,
        'colsample_bytree': 0.9494974450800385,
        'colsample_bylevel': 0.9770967962524991,
        'colsample_bynode': 0.9648182200894126
    }
    #######################################################################################