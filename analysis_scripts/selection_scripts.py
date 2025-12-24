import numpy as np
import pandas as pd
from sklearn.metrics import auc
from typing import Tuple, Union, List, Optional, Dict
import sys


sys.path.append('../analysis_scripts')

from estimate_scripts import signal_estimates


###########################################################################################################
# auto_preselection
###########################################################################################################

def auto_preselection(df: pd.DataFrame,
                      features: List,
                      safety_interval: float = 0.90,
                      indent: float = 2,
                      quantile_left: float = 1e-3,
                      quantile_right: float = 1-1e-3,
                      recursive: bool = False) -> Tuple[pd.DataFrame, str]:
    
    """Clean Data from outliers using quantile-based bounds.

    Args:
        - df (pandas.DataFrame): Input dataset containing raw data with columns for features and event tags.
        - features (list): List of feature column names to check for outliers.
        - safety_interval (float, optional): Central proportion of data used to define the "normal" range.
            The safe zone is calculated as mean +- safety_interval * indent. Must be between 0 and 1.Defaults to 0.90.
        - indent (float, optional): Multiplier that controls the width of the safe zone boundary. Defaults to 2.
        - quantile_left (float, optional): Lower quantile threshold for outlier removal.
            Data below this quantile will be considered outliers. Defaults to 1e-3.
        - quantile_right (float, optional): Upper quantile threshold for outlier removal.
            Data above this quantile will be considered outliers. Defaults to 1-1e-3.
        - recursive (bool, optional): If True, applies outlier removal sequentially feature-by-feature,
            where each subsequent feature is cleaned on the already-filtered dataset.
            If False, calculates all cuts on the original dataset. Defaults to False.

    Returns:
        tuple: A tuple containing two elements:
        - selection_df (pandas.DataFrame): Summary table of cleaning operations with columns:
            * cut_feature: Feature name where cuts were applied
            * left_cut: Lower bound value (NaN if no lower cut)
            * right_cut: Upper bound value (NaN if no upper cut)
            * sig_efficiency: Fraction of signal events passing the cut
            * MB_efficiency: Fraction of minimum bias events passing the cut  
            * OC_bg_efficiency: Fraction of open charm background events passing the cut
        - mask (str): Boolean expression string combining all individual cuts with 'and' operators,
    """
    
    # Input validation
    if df.empty:
        raise ValueError("Input DataFrame is empty")
    
    if not features:
        raise ValueError("No features provided")
    
    missing_features = [f for f in features if f not in df.columns]
    if missing_features:
        raise ValueError(f"Features not in DataFrame: {missing_features}")
    
    # Data preparation
    df_MB = df[df.tag == 'Bg'].reset_index(drop=True).copy()
    # df_OC_bg = df.query("(tag == 'OC' and true_decay == False)").reset_index(drop=True).copy()
    df = df[df.tag == 'Sig'].reset_index(drop=True).copy()
    
    cut_features = []
    left_cuts = []
    right_cuts = []
    
    mask_parts = []
    sig_efficiency_list = []
    # OC_bg_efficiency_list = []
    MB_efficiency_list = []
    
    for feature in features:
        # Skip if feature has no data
        if df[feature].empty:
            continue
        
        # Calculate safety range
        up_quantile = safety_interval + (1 - safety_interval) / 2
        down_quantile = (1 - safety_interval) / 2
        check_range = df[feature].quantile(up_quantile) - df[feature].quantile(down_quantile)
        
        # Calculate cuts
        left_cut = df[feature].quantile(quantile_left)
        right_cut = df[feature].quantile(quantile_right)
        
        # Check if outliers exist
        cond1 = df[feature].min() < df[feature].mean() - indent * check_range
        cond2 = df[feature].max() > df[feature].mean() + indent * check_range
        
        selector = ''
        
        if cond1 and cond2:
            selector = f'{left_cut} <= {feature} and {feature} <= {right_cut}'
            cut_features.append(feature)
            left_cuts.append(left_cut)
            right_cuts.append(right_cut)
        elif cond1:
            selector = f'{left_cut} <= {feature}'
            cut_features.append(feature)
            left_cuts.append(left_cut)
            right_cuts.append(np.nan)
        elif cond2:
            selector = f'{feature} <= {right_cut}'
            cut_features.append(feature)
            left_cuts.append(np.nan)
            right_cuts.append(right_cut)
        else:
            # No outliers detected
            continue
        
        # Handle special angular features
        if feature in (
            'cosAngle_r_Lc_momentum_Lc_xy_left',
            'cosAngle_r_Lc_momentum_Lc_xy_right',
            'cosAngle_r_Lc_sum_momentum_xy_left',
            'cosAngle_r_Lc_sum_momentum_xy_right'
        ):
            suff = feature.split('_')[-1]
            base_feature = feature.rstrip(f'_{suff}')
            
            selector = selector.replace(f"{feature}", f"{base_feature}")
            
            if suff == 'left':    
                selector = f'({selector} or {base_feature} >= 0)'
            elif suff == 'right':
                selector = f'({selector} or {base_feature} <= 0)'
                    
        # Calculate efficiencies
        sig_efficiency = df.query(selector)[feature].count() / df[feature].count() if not df.empty else 0
        # OC_bg_efficiency = df_OC_bg.query(selector)[feature].count() / df_OC_bg[feature].count() if not df_OC_bg.empty else 0
        MB_efficiency = df_MB.query(selector)[feature].count() / df_MB[feature].count() if not df_MB.empty else 0
        
        sig_efficiency_list.append(sig_efficiency)
        # OC_bg_efficiency_list.append(OC_bg_efficiency)
        MB_efficiency_list.append(MB_efficiency)
        
        mask_parts.append(selector)
      
        # Apply recursive filtering if requested
        if recursive:
            df = df.query(selector)
            # df_OC_bg = df_OC_bg.query(selector)
            df_MB = df_MB.query(selector)
            
    # Create final mask
    mask = ' and '.join(mask_parts)
    
    # Create selection summary            
    selection_df = pd.DataFrame({
        'cut_feature': cut_features,
        'left_cut': left_cuts,
        'right_cut': right_cuts,
        'sig_efficiency': sig_efficiency_list,
        'MB_efficiency': MB_efficiency_list,
        # 'OC_bg_efficiency': OC_bg_efficiency_list
    })
    
    return selection_df, mask


###########################################################################################################
# calculate_rates_at_cut_point
###########################################################################################################

def calculate_rates_at_cut_point(
    df: pd.DataFrame,   # mass, feature, tag
    feature: str,
    x_cut: float,
    selection_direction: str,
    bounds: Optional[Tuple[float, float]] = None,
    mass_interval: Tuple[float, float] = (2.24763, 2.32497)
) -> Union[Tuple[float, float], float]:
    
    # Filter by mass interval
    mass_mask = (df['mass_Lc'] >= mass_interval[0]) & (df['mass_Lc'] <= mass_interval[1])
    
    df = df[mass_mask].copy()
    
    # Input validation
    if 'tag' not in df.columns:
        raise KeyError('Input Data Frame should contain "tag" feature.')
    
    unique_tags = np.sort(df['tag'].unique())
    if len(unique_tags) != 2 or unique_tags[0] != 'Bg' or unique_tags[1] != 'Sig':
        raise KeyError('Dataset "tag" should contain both "Sig" and "Bg" values')
    
    if selection_direction not in ['left', 'right']:
        raise ValueError("selection_direction must be 'left' or 'right'")
    
    if feature not in df.columns:
        raise ValueError(f"Feature '{feature}' not found in DataFrame columns")
    
    if bounds is not None:
        # Filter by bounds
        bounds_mask = (df[feature] >= bounds[0]) & (df[feature] <= bounds[1])
        
        df = df[bounds_mask]
    
    if not df.shape[0]:
        raise ValueError('Check bounds coundition. There are no events in dataset after bounds restriction.')
    
    
    if selection_direction == 'right':
            selector = f'{feature} > {x_cut}'
    else:
        selector = f'{feature} < {x_cut}'
    
    selected_df = df.query(selector)
    rejected_df = df[~df.index.isin(selected_df.index)]
    
    sig_selected = selected_df[selected_df.tag == 'Sig']
    bg_selected = selected_df[selected_df.tag == 'Bg']
    
    sig_rejected = rejected_df[rejected_df.tag == 'Sig']
    bg_rejected = rejected_df[rejected_df.tag == 'Bg']
    
    if len(sig_selected) == 0:
        # print('There are no signal left! Return (0, 0)')
        return 0, 0, 0, 0
    elif len(bg_selected) == 0:
        # print('There are no backgrounds left! Return (0, 0)')
        return 0, 0, 0, 0

    tp = len(sig_selected)
    fn = len(sig_rejected)
    fp = len(bg_selected)
    tn = len(bg_rejected)
    
    return tp, fn, fp, tn
    

###########################################################################################################
# find_optimal_cut_point
###########################################################################################################

def find_optimal_cut_point(
    df: pd.DataFrame, # mass, feature, tag
    feature: str,
    nbins: int = 1000,
    select_direction: str = 'right',
    bounds: Tuple = None, 
    min_sig_sel: float = 0.3,
    metric_type: str = 'tpr_fpr',
    mass_interval: Tuple[float, float] = (2.24763, 2.32497)
):
    
    # Filter by mass interval
    mass_mask = (df['mass_Lc'] >= mass_interval[0]) & (df['mass_Lc'] <= mass_interval[1])
    
    df = df[mass_mask].copy()
    
    # Input validation
    if 'tag' not in df.columns:
        raise KeyError('Input Data Frame should contain "tag" feature.')
    
    unique_tags = np.sort(df['tag'].unique())
    if len(unique_tags) != 2 or unique_tags[0] != 'Bg' or unique_tags[1] != 'Sig':
        raise KeyError('Dataset "tag" should be either "Sig" or "Bg"')
    
    if select_direction not in ['left', 'right']:
        raise ValueError("selection_direction must be 'left' or 'right'")
    
    if metric_type not in ['tpr_fpr', 'f1']:
        raise ValueError(f"Invalid metric_type: {metric_type}. \n Should be 'tpr_fpr' or 'f1'.")
    
    if feature not in df.columns:
        raise ValueError(f"Feature '{feature}' not found in DataFrame columns")
    
    if not (0 <= min_sig_sel <= 1):
        raise ValueError("min_sig_sel must be between 0 and 1")
    
    if bounds is not None:
        bounds_mask = (df[feature] >= bounds[0]) & (df[feature] <= bounds[1])
        df = df[bounds_mask]
    
    if not df.shape[0]:
        raise ValueError('Check bounds coundition. There are no events in dataset after bounds restriction.')
    
    # Convert inputs to numpy arrays
    distr_sig = np.asarray(df.loc[df.tag == 'Sig', feature])
    distr_bg = np.asarray(df.loc[df.tag == 'Bg', feature])
        
    # Determine threshold range covering both distributions
    start_x = min(np.min(distr_sig[~np.isnan(distr_sig)]), np.min(distr_bg[~np.isnan(distr_bg)]))
    stop_x = max(np.max(distr_sig[~np.isnan(distr_sig)]), np.max(distr_bg[~np.isnan(distr_bg)]))

    # Create evenly spaced thresholds
    thresholds = np.linspace(start_x, stop_x, nbins)
    
    if len(thresholds) == 0:
        raise ValueError("No thresholds generated - check input distributions")
    
    if select_direction == 'right':
        # Calculate fraction of signal above each threshold
        sig_efficiencies = np.mean(distr_sig > thresholds.reshape(-1, 1), axis=1)
    else:
        # Calculate fraction of signal below each threshold
        sig_efficiencies = np.mean(distr_sig < thresholds.reshape(-1, 1), axis=1)
    
    # Create mask for thresholds meeting minimum signal efficiency requirement
    min_sig_sel_mask = sig_efficiencies >= min_sig_sel

    # Check if any thresholds meet the minimum requirement
    if not np.any(min_sig_sel_mask):
        raise ValueError(f'No selection above minimum signal efficiency {min_sig_sel}!')

    tpr = []
    fpr = []
    metric = []

    eps = 1e-10

    for num, cut_x in enumerate(thresholds):
        
        if not num % 200:
            print(f'Cut point search progress: {num}/{len(thresholds)}')
        
        tp, fn, fp, tn = calculate_rates_at_cut_point(
            df=df,  # mass , feature, tag
            feature=feature,
            bounds=bounds,
            x_cut=cut_x,
            selection_direction=select_direction,
            mass_interval=(2.24763, 2.32497)
        )
        
        tpr_val = tp / (tp + fn + eps)
        fpr_val = fp / (fp + tn + eps)
        
        if metric_type == 'tpr_fpr':
            metric_val = tpr_val - fpr_val
        elif metric_type == 'f1':
            rec = tp / (tp + fn + eps)
            pre = tp / (tp + fp + eps)
            metric_val = 2 * pre * rec / (pre + rec + eps)
         
        metric.append(metric_val)
        tpr.append(tpr_val)
        fpr.append(fpr_val)
        

    metric = np.array(metric)
    tpr = np.array(tpr)
    fpr = np.array(fpr)

    # Apply mask to consider only thresholds meeting minimum signal efficiency
    selected_metric = metric[min_sig_sel_mask]
    temp_thresholds = thresholds[min_sig_sel_mask]

    best_arg = np.argmax(selected_metric)
    best_threshold = temp_thresholds[best_arg]
        
    # Calculate final efficiencies at best threshold
    if select_direction == 'right':
        sig_efficiency = np.sum(distr_sig > best_threshold) / np.sum(~np.isnan(distr_sig))
        bg_efficiency = np.sum(distr_bg > best_threshold) / np.sum(~np.isnan(distr_bg))
    else:
        sig_efficiency = np.sum(distr_sig < best_threshold) / np.sum(~np.isnan(distr_sig))
        bg_efficiency = np.sum(distr_bg < best_threshold) / np.sum(~np.isnan(distr_bg))

    return min_sig_sel_mask, best_arg, metric, thresholds, tpr, fpr, sig_efficiency, bg_efficiency


###########################################################################################################
# calculate_feature_importance
###########################################################################################################

def calculate_feature_importance(
    df: pd.DataFrame,
    features: List[str],
    direction_restrictions: Optional[Dict[str, str]] = None,
    metric_type: str = 'tpr_fpr',       # Respond only for cut point selection, feature importance ranking controls by AUC
    mass_interval: Tuple[float, float] = (2.24763, 2.32497)
) -> pd.DataFrame:
    
    # Validate input parameters
    if metric_type not in ['tpr_fpr', 'f1']:
        raise ValueError(f"Invalid metric_type: {metric_type}. \n Must be one of: 'tpr_fpr', 'f1'")
    
    # Check required columns
    required_columns = ['mass_Lc', 'tag'] + features
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise KeyError(f"Missing required columns: {missing_columns}")
    
    metric_value_list = []
    feature_list = []
    featuree_sel_direction_list = []
    signal_efficiency_list = []
    background_efficiency_list = []
    cut_list = []
    
    if direction_restrictions is None:
        direction_restrictions = {}
    
    for num, feature in enumerate(features):
        
        print('*' * 40)
        print(f'Progress: {num + 1}/{len(features)}')
        print(f'Process feature: {feature}')
        
        # Initialize temp values
        top_metric_value = 0
        top_feature = None
        top_feature_direction = None
        correspond_signal_efficiency = None
        correspond_bg_efficiency = None
        best_cut = None
                
        # Selection direction optimization   
        for direction in ['left', 'right']:
            
            # Apply direction restrictions if specified
            if feature in direction_restrictions and direction != direction_restrictions[feature]:
                print(f'Direction restriction encounter. Skip direction: {direction}')
                continue
            
            bounds = None
            min_sig_sel =  0.3
            
            min_sig_sel_mask, best_arg, metric, thresholds, tpr, fpr, sig_efficiency, bg_efficiency = find_optimal_cut_point(
                df=df[['mass_Lc', feature, 'tag']],
                feature=feature,
                bounds=bounds,
                select_direction=direction,
                metric_type=metric_type,
                min_sig_sel=min_sig_sel,
                mass_interval=mass_interval
            )
        
            # Check if we have valid results
            if (not min_sig_sel_mask.any() or 
                best_arg >= len(thresholds[min_sig_sel_mask]) or
                len(metric[min_sig_sel_mask]) == 0):
                raise ValueError(f"No valid cuts found for {feature} with {direction} direction")
        
            best_cut_x = thresholds[min_sig_sel_mask][best_arg]
            # best_metric = metric[min_sig_sel_mask][best_arg]
            
            # Feature importance evaluates by AUC independently from min_sig_sel_mask
            metric_value = auc(fpr, tpr)
            
            print(f'AUC value: {metric_value:.4f}, direction: {direction}')
                
            # Update best AUC results
            if metric_value > top_metric_value:
                top_metric_value = metric_value
                top_feature = feature
                top_feature_direction = direction
                correspond_signal_efficiency = sig_efficiency
                correspond_bg_efficiency = bg_efficiency
                best_cut = best_cut_x

                
        if top_feature is None:
            raise ValueError('No selected feature after direction optimization')
        
        metric_value_list.append(top_metric_value)   
        feature_list.append(top_feature)
        featuree_sel_direction_list.append(top_feature_direction)
        signal_efficiency_list.append(correspond_signal_efficiency)
        background_efficiency_list.append(correspond_bg_efficiency)
        cut_list.append(best_cut)
        
        
    # Check if we have any results
    if not feature_list:
        raise ValueError("No valid features found after processing all features")
    
    # Create results dataframes
    importance_df = pd.DataFrame({
        f'AUC': metric_value_list,
        f'feature': feature_list,
        f'feature_sel_direction': featuree_sel_direction_list,
        f'cut': cut_list,
        f'signal_efficiency': signal_efficiency_list,
        f'background_efficiency': background_efficiency_list
    }) 
    
    # Sort by importance metrics (descending order)
    importance_df = importance_df.sort_values(by=f'AUC', ascending=False).reset_index(drop=True)
    
    return importance_df


###########################################################################################################
# create_best_selection_path
###########################################################################################################

def create_best_selection_path(
    df: pd.DataFrame, 
    features: List[str], 
    metric_type: str, 
    have_sig_events: int,
    have_bg_events: int, 
    n_features_to_use: int = None, 
    mass_interval=(2.24763, 2.32497),
    direction_restrictions: Dict = None # Format: {'feature_name': 'left'/'right'}
) -> pd.DataFrame:
    
    # Input validation
    if metric_type not in ['tpr_fpr', 'f1']:
        raise ValueError(
            f"Invalid metric_type: {metric_type}. \n Must be one of: 'tpr_fpr', 'f1'")
    
    # Validate required columns
    required_columns = ['tag', 'mass_Lc'] + features
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    # Create copies to avoid modifying original data
    df = df.copy()
    features_to_use = features.copy()
    
    # Number of iterations
    if n_features_to_use is not None:
        num_iter = min(len(features), n_features_to_use)
    else:
        num_iter = len(features)
        
    selected_feature_list = []
    cut_point_list = []
    select_direction_list = []
    signal_efficiency_list = []
    background_efficiency_list = []
    ratio_list = []
    significance_list = [] 
    total_sig_unscaled_list = []
    total_bg_unscaled_list = []
        
    # Sequential feature selection
    for i in range(num_iter):
    
        print('=' * 70)
        print(f'Selection Progress: {i + 1}/{num_iter}')
        print(f'Remaining features: {len(features_to_use)}')
        print(f'Remaining events - Signal: {len(df[df["tag"] == "Sig"])}, '
              f'Background: {len(df[df["tag"] == "Bg"])}')
        
        # Calculate signal and background distributions after selection
        sig_mass_distr = df.loc[(df['tag'] == 'Sig'), 'mass_Lc'].copy()
        bg_mass_distr = df.loc[(df['tag'] == 'Bg'), 'mass_Lc'].copy()
        
        # Early stopping if no events remain
        if len(sig_mass_distr) == 0:
            print('Warning: No Signal events remaining after selection. Stop Path Search!')
            break
        elif len(bg_mass_distr) == 0:
            print('Warning: No Background events remaining after selection. Stop Path Search!')
            break
        
        # Calculate feature importance
        try:  
            
            importance_df = calculate_feature_importance(
                df=df[['mass_Lc', 'tag']+features_to_use],
                features=features_to_use,
                direction_restrictions=direction_restrictions,
                metric_type=metric_type,
                mass_interval=mass_interval
            )
            
        except (TypeError, ValueError) as e:
            print(f'Feature importance calculation failed: {e}. Stop path search.')
            # raise e
            break
        
        # Check if we have valid results
        if importance_df.empty:
            print('No valid features found. Stop path search.')
            break
        
        # Get the best feature and its properties
        top_feature_row = importance_df.iloc[0]
        best_cut_point = top_feature_row['cut']
        top_feature = top_feature_row['feature']
        top_feature_direction = top_feature_row['feature_sel_direction']
        top_signal_efficiency = top_feature_row['signal_efficiency']
        top_background_efficiency = top_feature_row['background_efficiency']
        
        print(f'Selected feature: {top_feature}, Direction: {top_feature_direction}, Cut: {best_cut_point:.4f}')    
        
        # Remove selected feature from future consideration
        features_to_use.remove(top_feature)
        
        # Apply selection to data
        if top_feature_direction == 'right':
            df = df.loc[df[top_feature] > best_cut_point].reset_index(drop=True)
        else:
            df = df.loc[df[top_feature] < best_cut_point].reset_index(drop=True)
        
        df = df.drop(top_feature, axis=1).reset_index(drop=True)
        
        # Calculate signal and background distributions after selection
        sig_mass_distr = df.loc[(df['tag'] == 'Sig'), 'mass_Lc'].copy()
        bg_mass_distr = df.loc[(df['tag'] == 'Bg'), 'mass_Lc'].copy()
        
        # Early stopping if no events remain
        if len(sig_mass_distr) == 0:
            print('Warning: No Signal events remaining after selection. Stop Path Search!')
            break
        elif len(bg_mass_distr) == 0:
            print('Warning: No Background events remaining after selection. Stop Path Search!')
            break
        
        # Calculate performance metrics using signal estimation
        try:
            overall_s_b, overall_s_sqrt_s_b, total_sig, total_bg, total_sig_unscaled, total_bg_unscaled = signal_estimates(
                sig_mass_distr=sig_mass_distr,
                bg_mass_distr=bg_mass_distr,
                have_sig_events=have_sig_events,
                have_bg_events=have_bg_events, 
                mass_interval=mass_interval,
                visualization=False,
                verbose=False
            )
        except Exception as e:
            print(f'Error in signal estimation: {e}. Stop path search.')
            break
        
        # Store results
        selected_feature_list.append(top_feature)
        cut_point_list.append(best_cut_point)
        select_direction_list.append(top_feature_direction)
        signal_efficiency_list.append(top_signal_efficiency)
        background_efficiency_list.append(top_background_efficiency)
        ratio_list.append(overall_s_b)
        significance_list.append(overall_s_sqrt_s_b) 
        total_sig_unscaled_list.append(total_sig_unscaled)
        total_bg_unscaled_list.append(total_bg_unscaled)
        
            
        # Check if we have any features left to select
        if not features_to_use:
            print('No more features available for selection.')
            break
        
    else:
        print(f'Selection Path Search Finished Successfully! Total features processed: {i + 1}')
        
    # Create results dataframe
    resuls_df = pd.DataFrame({
        'selected_feature': selected_feature_list,
        'cut_point': cut_point_list,
        'select_direction': select_direction_list,
        'signal_efficiency': signal_efficiency_list,
        'background_efficiency': background_efficiency_list,
        'ratio': ratio_list,
        'overall_s_sqrt_s_b': significance_list,
        'total_sig_unscaled': total_sig_unscaled_list,
        'total_bg_unscaled': total_bg_unscaled_list
    })
    
    return resuls_df


