import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from scipy.optimize import minimize
from typing import Tuple, Union, List, Optional


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
        if feature in ('cosAngle_r_Lc_momentum_Lc_left','cosAngle_r_Lc_sum_momentum_left',
                       'cosAngle_r_Lc_momentum_Lc_right', 'cosAngle_r_Lc_sum_momentum_right'):
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


def draw_distribution_pair(df: pd.DataFrame,
                       distr_name: str,
                       bins: int = 40) -> Tuple:
    """Draw combined boxplot and histogram for feature distribution analysis."""
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

    # Boxplot
    sns.boxplot(
        data=df,
        x=distr_name,
        color='#2E86AB',
        whis=1.5,
        ax=axes[0]
    )

    axes[0].set_title(f'{distr_name} Distribution', fontsize=14, fontweight='bold', pad=20)
    axes[0].set_xlabel(f'{distr_name}', fontweight='bold')

    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    axes[0].grid(True, alpha=0.8, linestyle='--')

    # Histplot
    hist_ax = sns.histplot(
        data=df,
        x=distr_name,
        color='#2E86AB',
        alpha=1,
        element='step',
        fill=False,
        linewidth=2,
        bins=bins,
        ax=axes[1]
    )

    axes[1].set_title(f'{distr_name} Distribution', fontsize=14, fontweight='bold', pad=20)
    axes[1].set_xlabel(f'{distr_name}', fontweight='bold')
    axes[1].set_ylabel('Counts', fontweight='bold')

    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    axes[1].grid(True, alpha=0.8, linestyle='--')
    
    return fig, axes
    

def draw_feature_distribution(df: pd.DataFrame,
                              distr_name: str,
                              tag: str = None,
                              hue: str = None,
                              bins: int = 40,
                              norma: bool = False,
                              # include_OC_bg: bool = False,
                              cut_point: float = None,
                              select_direction: str = None) -> Tuple:
    """Draw combined boxplot and histogram for feature distribution analysis.

    Args:
        - df (pd.DataFrame): Input DataFrame
        - distr_name (str): Name of the feature/column to visualize
        - tag (str, optional): Filter data to include only specific tag values. 
            If None, uses all available data. Defaults to None.
        - hue (str, optional): Column name to use for color encoding/splitting distributions.
            Creates separate distributions for each unique value. Defaults to None.
        - bins (int, optional): Number of bins to use in histogram. Defaults to 40.
        - norma (bool, optional): If True, normalizes histograms to display density.
            If False, shows raw counts. Defaults to False.
        - include_OC_bg (bool, optional): If False, filters OC events to only include true decays
            and renames 'OC' tag to 'Sig'. If True, uses all OC events including background.
            Defaults to False.
        - cut_point (float, optional): Threshold value to display as vertical line on histogram.
            Requires select_direction to be specified. Defaults to None.
        - select_direction (str, optional): Selection direction for threshold annotation.
            Must be 'left' or 'right' if cut_point is provided. Determines arrow direction.
            Defaults to None.
            
    Returns:
        - fig: created figure
        - axes: figure axes
            
        Notes:
        - When include_OC_bg=False (default): Filters data to MB events and true OC decays,
          renaming 'OC' tag to 'Sig' for clarity
        - Creates a 2-panel figure with boxplot (left) and histogram (right)
        - Includes optional threshold line with directional arrow when cut_point and 
          select_direction are both specified
    """

    df = df.copy()

    if tag is not None:
        df = df[df.tag == tag]

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

    # Boxplot
    sns.boxplot(
        data=df,
        x=distr_name,
        hue=hue,
        palette='Set2' if hue is not None else None,
        color=None if hue is not None else '#2E86AB',
        whis=1.5,
        ax=axes[0]
    )

    axes[0].set_title(f'{distr_name} Distribution', fontsize=14, fontweight='bold', pad=20)
    axes[0].set_xlabel(f'{distr_name}', fontweight='bold')

    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    axes[0].grid(True, alpha=0.8, linestyle='--')

    # Histplot
    hist_ax = sns.histplot(
        data=df,
        x=distr_name,
        hue=hue,
        palette='Set2' if hue is not None else None,
        color=None if hue is not None else '#2E86AB',
        stat="density" if norma else 'count',
        common_norm=False,
        alpha=1,
        element='step',
        fill=False,
        linewidth=2,
        bins=bins
    )

    axes[1].set_title(f'{distr_name} Distribution', fontsize=14, fontweight='bold', pad=20)
    axes[1].set_xlabel(f'{distr_name}', fontweight='bold')
    y_label = 'Density' if norma else 'Counts'
    axes[1].set_ylabel(y_label, fontweight='bold')

    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    axes[1].grid(True, alpha=0.8, linestyle='--')

    # Thrashold Line
    if cut_point is not None and select_direction is not None:
        y_max = hist_ax.get_ylim()[1]
        vline_height = y_max * 0.95
        axes[1].axvline(cut_point, color='black', linestyle='dashed', linewidth=2, ymin=0, ymax=vline_height)

        data = df[distr_name].dropna()
        counts, bin_edges = np.histogram(data, bins=bins)
        bin_width = bin_edges[1] - bin_edges[0]

        if select_direction == 'right':
            arrow_direction = cut_point + bin_width * 3
        else:
            arrow_direction = cut_point - bin_width * 3
            
        axes[1].annotate('', xy=(cut_point, vline_height), xytext=(arrow_direction, vline_height),
                    arrowprops=dict(arrowstyle='<-', color='black', lw=2))

    plt.tight_layout()
    
    return fig, axes


def draw_roc(tpr: List,
             fpr: List,
             best_cut_x: float,
             best_tpr: float,
             best_fpr: float,
             select_direction: str,
             mask: np.ndarray = None,
             min_sig_sel: float = None) -> Tuple:
    """
    Plot a Receiver Operating Characteristic (ROC) curve with optimal threshold point.
    
    Args:
        - tpr : array-like
            True Positive Rate values for ROC curve
        - fpr : array-like 
            False Positive Rate values for ROC curve
        - best_cut_x : float
            Optimal threshold value
        - best_tpr : float
            True Positive Rate at optimal threshold
        - best_fpr : float
            False Positive Rate at optimal threshold  
        - select_direction : str
            Direction of selection criterion ('right' for > threshold, 'left' for < threshold)
        - mask : array-like, optional
            Boolean mask to separate main ROC curve from dashed section
        - min_sig_sel : float, optional
            Minimum signal efficiency value, used for labeling dashed section
    
    Returns:
        fig, ax : matplotlib figure and axes objects
            The generated plot figure and axes
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Prepare main ROC curve data
    if mask is not None:
        roc_data = pd.DataFrame({
            'False Positive Rate': fpr[mask],
            'True Positive Rate': tpr[mask]     
        })
    else:
        roc_data = pd.DataFrame({
            'False Positive Rate': fpr,
            'True Positive Rate': tpr     
        })

    # Calculate AUC score
    roc_auc = auc(fpr, tpr)

    # Plot main ROC curve
    sns.lineplot(
        data=roc_data,
        x='False Positive Rate',
        y='True Positive Rate',
        color='darkorange',
        linewidth=2,
        label=f'ROC curve (AUC = {roc_auc:.3f})',
        ax=ax
    )
    
    # Plot dashed section if mask and min_sig_sel are provided
    if mask is not None and min_sig_sel is not None:
        
        # Create dataframe for dashed section (where mask is False)
        roc_data_dash = pd.DataFrame({
            'False Positive Rate Dash': fpr[~ (mask)],
            'True Positive Rate Dash': tpr[~ (mask)]   
        })
        
        # Plot dashed ROC curve section
        sns.lineplot(
            data=roc_data_dash,
            x='False Positive Rate Dash',
            y='True Positive Rate Dash',
            color='darkorange',
            linewidth=2,
            linestyle='--',
            label=f'Zone where Sig Eff < {min_sig_sel}',
            ax=ax
        )

    # Plot optimal threshold point
    sns.scatterplot(
        x=[best_fpr],
        y=[best_tpr],
        color='red',
        s=100,
        label='Optimal threshold',
        zorder=5,
        ax=ax
    )

    # Add diagonal line for random classifier reference
    ax.plot([0, 1], [0, 1], color='navy', linestyle='--', 
            label='Random classifier (AUC = 0.5)', alpha=0.8)

    # Define selection direction symbols
    sel_dir = {
        'right': '>',
        'left': '<'
    }

    # Validate select_direction parameter
    if select_direction not in sel_dir:
        raise ValueError(f"select_direction must be 'right' or 'left', got '{select_direction}'")

    # Add threshold annotation text
    ax.text(best_fpr + 0.02, best_tpr - 0.05,
            f'Threshold: {sel_dir[select_direction]} {best_cut_x:.3f}\n'
            f'TPR: {best_tpr:.3f}, FPR: {best_fpr:.3f}',
            fontsize=10, 
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Set plot limits and labels
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('Receiver Operating Characteristic (ROC) Curve', fontsize=14, pad=20)
    
    # Add grid and legend
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower right', framealpha=0.9)
    
    # Improve tick formatting
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Tight layout for better spacing
    plt.tight_layout()
    
    return fig, ax


def signal_estimates(sig_mass_distr: Union[pd.Series, np.ndarray], 
                     bg_mass_distr: Union[pd.Series, np.ndarray], 
                     have_sig_events: int, 
                     have_bg_events: int, 
                     mass_interval: Tuple[float, float] = (2.24763, 2.32497), 
                     visualization: bool = False,
                     verbose: bool = True) -> Tuple[float, float, int, int, int, int]:
    """
    Estimate signal significance metrics.
    
    This function calculates signal-to-background ratios and significance metrics
    for a given mass distribution, scaling Monte Carlo samples to expected real-world
    event counts based on cross-sections and luminosity.
    
    Args:
        - sig_mass_distr : pd.Series or np.ndarray
            Mass distribution of signal events from simulation
        - bg_mass_distr : pd.Series or np.ndarray  
            Mass distribution of background events from simulation
        - have_sig_events : int
            Number of signal events in the input simulation
        - have_bg_events : int
            Number of background events in the input simulation
        - mass_interval : tuple of float, optional
            Mass range of interest (low, high) in GeV. Default: (2.24763, 2.32497)
        - visualization : bool, optional
            Whether to generate visualization plots. Default: False
    
    Returns:
        tuple : (overall_s_b, overall_s_sqrt_s_b, total_sig, total_bg, total_sig_unscaled, total_bg_unscaled)
            - overall_s_b : float
                Overall signal-to-background ratio
            - overall_s_sqrt_s_b : float
                Overall significance S/√(S+B)
            - total_sig : int
                Total scaled signal events
            - total_bg : int
                Total scaled background events
            - total_sig_unscaled: int
                Total unscaled signal events
            - total_bg_unscaled: int
                Total unscaled background events

    """
     # Physical constants and parameters
    sigma_MB = 38 * 1e-3 * 1e-24                # Minimum bias cross-section in cm2
    cross_section_Lc_plus = 4 * 1e-6 * 1e-24    # Lc+(4122) cross-section in cm2
    L = 1e32                                    # Luminosity in cm-2s-1
    t = 31536000                                # Time period in seconds (1 year)
    branching = 0.035                           # Lc decay branching ratio
    
    # Calculate expected real events from cross-sections and luminosity
    real_sig_events = L * t * cross_section_Lc_plus * branching
    real_bg_events = L * t * sigma_MB
    
    # Calculate scaling weights to match real event expectations
    sig_weight = real_sig_events / have_sig_events
    bg_weight = real_bg_events / have_bg_events
    
    # Convert inputs to numpy arrays for processing
    sig_mass_distr = np.asarray(sig_mass_distr)
    bg_mass_distr = np.asarray(bg_mass_distr)
    
    # Create weight arrays for signal and background
    weights_sig = np.full_like(sig_mass_distr, sig_weight, dtype=float)
    weights_bg = np.full_like(bg_mass_distr, bg_weight, dtype=float)
    
    # Combine signal and background data
    mass_combined = np.concatenate([sig_mass_distr, bg_mass_distr])
    weights_combined = np.concatenate([weights_sig, weights_bg])
    labels_combined = np.concatenate([np.full(len(sig_mass_distr), 'Signal'), 
                                      np.full(len(bg_mass_distr), 'Background')])

    mass_df = pd.DataFrame({
        'mass': mass_combined,
        'weights': weights_combined,
        'type': labels_combined,
    })

    # Filter by mass interval
    mass_mask = (mass_df['mass'] >= mass_interval[0]) & (mass_df['mass'] <= mass_interval[1])
    mass_df = mass_df[mass_mask]

     # Histogram configuration
    bins = 100
    binrange = (mass_df['mass'].min(), mass_df['mass'].max())

    # Calculate the same bins parameters as in histogram
    bin_edges = np.linspace(binrange[0], binrange[1], bins + 1)  # 100 bins + edges
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    
    # Calculate weighted histograms
    signal_data = mass_df[mass_df['type'] == 'Signal']
    background_data = mass_df[mass_df['type'] == 'Background']

    sig_heights, _ = np.histogram(
        signal_data['mass'],
        weights=signal_data['weights'],
        bins=bin_edges
    )

    bg_heights, _ = np.histogram(
        background_data['mass'],
        weights=background_data['weights'],
        bins=bin_edges
    )

    # Calculate significance metrics per bin
    s_b_ratio = np.divide(
        sig_heights,
        bg_heights,
        out=np.zeros_like(sig_heights),
        where=bg_heights > 0
    )
    
    s_over_sqrt_s_b = np.divide(
        sig_heights, np.sqrt(sig_heights + bg_heights),
        out=np.zeros_like(sig_heights),
        where=(sig_heights + bg_heights) > 0
    )

    # Calculate overall statistics
    total_sig = np.sum(sig_heights)
    total_bg = np.sum(bg_heights)
    total_sig_unscaled = np.sum(~np.isnan(sig_mass_distr))
    total_bg_unscaled = np.sum(~np.isnan(bg_mass_distr))
    overall_s_b = total_sig / total_bg if total_bg > 0 else 0
    overall_s_sqrt_s_b = total_sig / np.sqrt(total_sig + total_bg) if (total_sig + total_bg) > 0 else 0

    # Visualization
    if visualization:
        fig, axes = plt.subplots(3, 1, figsize=(12, 14), sharex=False, gridspec_kw={'height_ratios': [3, 1, 1]})

        colors = {'Signal': '#1f77b4', 'Background': 'black'}

        # Plot 1: Mass Distribution
        hist = sns.histplot(
            data=mass_df,
            x='mass',
            weights='weights',
            hue='type',
            hue_order=['Signal', 'Background'],
            bins=bins,
            binrange=binrange,
            palette=colors,
            log_scale=(False, True),
            alpha=0.8,
            element='step',
            fill=False,
            linewidth=2,
            ax=axes[0]
        )

        axes[0].set_title('Yield Mass Distribution', fontsize=16, fontweight='bold', pad=20)
        axes[0].set_xlabel('Mass (GeV)', fontweight='bold', fontsize=12)
        axes[0].set_ylabel('Counts (log scale)', fontweight='bold', fontsize=12)

        sns.move_legend(
            axes[0],
            loc='best',
            title=f'Basic Statistics:\nS/B: {overall_s_b:.6f}\nS/√(S+B): {overall_s_sqrt_s_b:.2f}',
            labels=['Signal (Open Charm)', 'Background (MB)'],
            frameon=True,
            fancybox=True,
            shadow=True
        )

        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)
        axes[0].grid(True, alpha=0.3, linestyle='--')

        # Plot 2: S/B Ratio
        sns.lineplot(
            x=bin_centers,
            y=s_b_ratio,
            color='#d62728',
            linewidth=2.5,
            label='S/B Ratio',
            ax=axes[1]
        )

        axes[1].set_ylabel('S/B Ratio', fontweight='bold', fontsize=12)
        axes[1].set_xlabel('Mass (GeV)', fontweight='bold', fontsize=12)
        axes[1].grid(True, alpha=0.3, linestyle='--')
        axes[1].legend(loc='upper right', frameon=True, fancybox=True, shadow=True)

        # Add mean S/B ratio line
        mean_s_b = np.mean(s_b_ratio[s_b_ratio > 0])
        axes[1].axhline(y=mean_s_b, color='black', linestyle='--', alpha=0.7, label=f'Mean S/B: {mean_s_b:.6f}')
        axes[1].legend(loc='upper right', frameon=True, fancybox=True, shadow=True)

        # Plot 3: S/sqrt(S+B) Ratio
        sns.lineplot(
            x=bin_centers,
            y=s_over_sqrt_s_b,
            color='#9467bd',
            linewidth=2.5,
            label='S/√(S+B) Ratio',
            ax=axes[2]
        )

        axes[2].set_ylabel('S/√(S+B)', fontweight='bold', fontsize=12)
        axes[2].set_xlabel('Mass (GeV)', fontweight='bold', fontsize=12)
        axes[2].grid(True, alpha=0.3, linestyle='--')

        # Add mean S/sqrt(S+B) ratio line
        mean_s_sqrt = np.mean(s_over_sqrt_s_b[s_over_sqrt_s_b > 0])
        axes[2].axhline(y=mean_s_sqrt, color='black', linestyle='--', alpha=0.7, label=f'Mean S/√(S+B): {mean_s_sqrt:.2f}')
        axes[2].legend(loc='upper right', frameon=True, fancybox=True, shadow=True)

        # Add peak value annotations for the ratios
        max_s_b_idx = np.argmax(s_b_ratio)
        max_s_sqrt_idx = np.argmax(s_over_sqrt_s_b)

        axes[1].annotate(f'Max: {s_b_ratio[max_s_b_idx]:.6f}',
                        xy=(bin_centers[max_s_b_idx], s_b_ratio[max_s_b_idx]),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

        axes[2].annotate(f'Max: {s_over_sqrt_s_b[max_s_sqrt_idx]:.2f}',
                        xy=(bin_centers[max_s_sqrt_idx], s_over_sqrt_s_b[max_s_sqrt_idx]),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

        plt.tight_layout()
        plt.show()

    if verbose:
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)
        print(f"Total Signal Events: {total_sig:.0f}")
        print(f"Total Background Events: {total_bg:.0f}")
        print(f"Total Signal Events Unscaled: {total_sig_unscaled:.0f}")
        print(f"Total Background Events Unscaled: {total_bg_unscaled:.0f}")
        print(f"Overall S/B Ratio: {overall_s_b:.6f}")
        print(f"Overall S/sqrt(S+B) Significance: {overall_s_sqrt_s_b:.2f}")
        if visualization:
            print(f"Maximum S/B Ratio: {s_b_ratio.max():.6f} at {bin_centers[max_s_b_idx]:.3f} GeV")
            print(f"Maximum S/sqrt(S+B): {s_over_sqrt_s_b.max():.2f} at {bin_centers[max_s_sqrt_idx]:.3f} GeV")

    return overall_s_b, overall_s_sqrt_s_b, total_sig, total_bg, total_sig_unscaled, total_bg_unscaled


def calculate_rates_at_cut_point(
    df: pd.DataFrame,
    feature: str,
    x_cut: float,
    selection_direction: str,
    bounds: Optional[Tuple[float, float]] = None,
    metric_type: str = 'tpr_fpr',
    have_sig_events: Optional[int] = None,
    have_bg_events: Optional[int] = None,
    mass_interval: Tuple[float, float] = (2.24763, 2.32497)
) -> Union[Tuple[float, float], float]:
    
    """
    Calculate performance metrics at a specified cut point for a given feature.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing at least 'mass_Lc', the specified feature, and 'tag' columns.
        The 'tag' column should indicate signal ('Sig') and background ('Bg') events.
    feature : str
        Name of the feature column to apply the cut on.
    x_cut : float
        The cut value to apply to the feature.
    selection_direction : str
        Direction of selection: 'left' for values less than x_cut, 'right' for values greater than x_cut.
    bounds : Tuple[float, float], optional
        Optional bounds to filter the feature values before applying the cut.
        If provided, only feature values within [bounds[0], bounds[1]] are considered.
    metric_type : str, default='tpr_fpr'
        Type of metric to calculate:
        - 'tpr_fpr': Returns (true positive rate, false positive rate)
        - 'efficiency': Returns (signal efficiency, background efficiency)
        - 'significance': Returns significance measures (requires have_sig_events and have_bg_events)
        - 'ratio': Returns significance measures (requires have_sig_events and have_bg_events)
    have_sig_events : int, optional
        Total number of signal events for significance calculation.
        Required for 'significance' and 'ratio' metric types.
    have_bg_events : int, optional
        Total number of background events for significance calculation.
        Required for 'significance' and 'ratio' metric types.
    mass_interval : Tuple[float, float], default=(2.24763, 2.32497)
        Mass interval to filter events before analysis.
        
    Returns
    -------
    Union[Tuple[float, float], float]
        Depending on metric_type:
        - 'tpr_fpr': Tuple of (true positive rate, false positive rate)
        - 'efficiency': Tuple of (signal efficiency, background efficiency)
        - 'significance' or 'ratio': Tuple of (significance measure, signal-to-background ratio)
    """
    
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
    
    if metric_type not in ['tpr_fpr', 'efficiency', 'significance', 'ratio']:
        raise ValueError(f"Invalid metric_type: {metric_type}")
    
    if feature not in df.columns:
        raise ValueError(f"Feature '{feature}' not found in DataFrame columns")
    
    # Validate required arguments for significance metrics
    if metric_type in ['significance', 'ratio']:
        if have_sig_events is None or have_bg_events is None:
            raise ValueError('"have_sig_events" and "have_bg_events" are required for "significance" and "ratio" metrics')
    
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
        return 0, 0
    elif len(bg_selected) == 0:
        # print('There are no backgrounds left! Return (0, 0)')
        return 0, 0
    
    # Get total signal and background counts
    total_sig = len(df[df.tag == 'Sig'])
    total_bg = len(df[df.tag == 'Bg'])
    
    # Calculate metrics
    if metric_type == 'tpr_fpr':
        tp = len(sig_selected)
        fn = len(sig_rejected)
        fp = len(bg_selected)
        tn = len(bg_rejected)
        
        epsilon = 1e-10
        tpr = tp / (tp + fn + epsilon)  # True Positive Rate (Sensitivity)
        fpr = fp / (fp + tn + epsilon)  # False Positive Rate
        
        return tpr, fpr
    
    elif metric_type == 'efficiency':
        
        sig_eff = len(sig_selected) / (total_sig + 1e-10)
        bg_eff = len(bg_selected) / (total_bg + 1e-10)
        
        return sig_eff, bg_eff
    
    elif metric_type in ['significance', 'ratio']:
        
        sig_mass_selected = sig_selected['mass_Lc']
        bg_mass_selected = bg_selected['mass_Lc']
        
        overall_s_b, overall_s_sqrt_s_b, total_sig, total_bg, total_sig_unscaled, total_bg_unscaled = signal_estimates(
            sig_mass_distr=sig_mass_selected, 
            bg_mass_distr=bg_mass_selected, 
            have_sig_events=have_sig_events, 
            have_bg_events=have_bg_events, 
            mass_interval=(2.24763, 2.32497), 
            visualization=False,
            verbose=False
        )
        
        return overall_s_sqrt_s_b, overall_s_b
    
        

def find_optimal_cut_point(
    df: pd.DataFrame, # only mass_Lc and faeture and tag
    feature: str,
    nbins: int = 1000,
    select_direction: str = 'right',
    bounds: Tuple = None, 
    min_sig_sel: float = 0.3,
    metric_type: str = 'tpr_fpr',
    have_sig_events: int = None,
    have_bg_events: int = None,
    mass_interval: Tuple[float, float] = (2.24763, 2.32497)
):
    
    # Filter by mass interval
    mass_mask = (df['mass_Lc'] >= mass_interval[0]) & (df['mass_Lc'] <= mass_interval[1])
    
    df = df[mass_mask].copy()
    
    # Input validation
    if 'tag' not in df.columns:
        raise KeyError('Input Data Frame should contain "tag" feature.')
    
    if (np.sort(df['tag'].unique())[0] != 'Bg') | (np.sort(df['tag'].unique())[1] != 'Sig'):
        raise KeyError('Dataset "tag" should be either "Sig" or "Bg"')
    
    if select_direction not in ['left', 'right']:
        raise ValueError("selection_direction must be 'left' or 'right'")
    
    if metric_type not in ['tpr_fpr', 'efficiency', 'significance', 'ratio']:
        raise ValueError(f"Invalid metric_type: {metric_type}")
    
    if feature not in df.columns:
        raise ValueError(f"Feature '{feature}' not found in DataFrame columns")
    
    if not (0 <= min_sig_sel <= 1):
        raise ValueError("min_sig_sel must be between 0 and 1")
    
    if bounds is not None:
        # Filter by bounds
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

    
    metric_1 = []
    metric_2 = []

    for num, cut_x in enumerate(thresholds):
        
        if not num % 200:
            print(f'Cut point search progress: {num}/{len(thresholds)}')
        
        new_metric_1, new_metric_2 = calculate_rates_at_cut_point(
            df=df,
            feature=feature,
            bounds=bounds,
            x_cut=cut_x,
            selection_direction=select_direction,
            metric_type=metric_type,
            have_sig_events=have_sig_events,
            have_bg_events=have_bg_events,
            mass_interval=(2.24763, 2.32497)
        )
                
        metric_1.append(new_metric_1)
        metric_2.append(new_metric_2)

    metric_1 = np.array(metric_1)
    metric_2 = np.array(metric_2)

    # Apply mask to consider only thresholds meeting minimum signal efficiency
    temp_metric_1 = metric_1[min_sig_sel_mask]
    temp_metric_2 = metric_2[min_sig_sel_mask] 
    temp_thresholds = thresholds[min_sig_sel_mask]


    if metric_type in ['tpr_fpr', 'efficiency']:
        best_arg = np.argmax(temp_metric_1 - temp_metric_2)
    elif metric_type in ['significance']:
        best_arg = np.argmax(temp_metric_1)
    elif metric_type in ['ratio']:
        best_arg = np.argmax(temp_metric_2)
        
    best_threshold = temp_thresholds[best_arg]
        
    # Calculate final efficiencies at best threshold
    if select_direction == 'right':
        sig_efficiency = np.sum(distr_sig > best_threshold) / np.sum(~np.isnan(distr_sig))
        bg_efficiency = np.sum(distr_bg > best_threshold) / np.sum(~np.isnan(distr_bg))
    else:
        sig_efficiency = np.sum(distr_sig < best_threshold) / np.sum(~np.isnan(distr_sig))
        bg_efficiency = np.sum(distr_bg < best_threshold) / np.sum(~np.isnan(distr_bg))

    return min_sig_sel_mask, best_arg, metric_1, metric_2, thresholds, sig_efficiency, bg_efficiency



from typing import List, Dict, Optional
from sklearn.metrics import auc


def calculate_feature_importance(
    df: pd.DataFrame,
    features: List[str],
    direction_restrictions: Optional[Dict[str, str]] = None,
    metric_type: str = 'tpr_fpr',
    have_sig_events: Optional[int] = None,
    have_bg_events: Optional[int] = None,
    mass_interval: Tuple[float, float] = (2.24763, 2.32497)
) -> pd.DataFrame:
    
    """
    Calculate feature importance by finding optimal cut points for each feature.
    It considers both selection directions ('left' and 'right') and can incorporate direction restrictions.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing at least 'mass_Lc', the specified features, and 'tag' columns.
        The 'tag' column should indicate signal ('Sig') and background ('Bg') events.
    features : List[str]
        List of feature names to evaluate for importance.
    direction_restrictions : Dict[str, str], optional
        Dictionary specifying direction restrictions for specific features.
        Keys are feature names, values are either 'left' or 'right'.
        If provided, only the specified direction will be considered for that feature.
        Format: {'feature_name': 'left'/'right'}
    metric_type : str, default='tpr_fpr'
        Type of metric to optimize:
        - 'tpr_fpr': Uses AUC of TPR-FPR curve (ROC AUC)
        - 'efficiency': Uses AUC of signal efficiency vs background efficiency curve
        - 'significance': Uses significance measure (requires have_sig_events and have_bg_events)
        - 'ratio': Uses signal-to-background ratio (requires have_sig_events and have_bg_events)
    have_sig_events : int, optional
        Total number of signal events for significance calculation.
        Required for 'significance' and 'ratio' metric types.
    have_bg_events : int, optional
        Total number of background events for significance calculation.
        Required for 'significance' and 'ratio' metric types.
    mass_interval : Tuple[float, float], default=(2.24763, 2.32497)
        Mass interval to filter events before analysis.
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing feature importance results sorted by metric value (descending).
        Columns include:
        - 'metric': The optimized metric value for each feature
        - 'feature': The feature name
        - 'feature_sel_direction': Optimal selection direction ('left' or 'right')
        - 'cut': Optimal cut value for the feature
        - 'signal_efficiency': Signal efficiency at optimal cut
        - 'background_efficiency': Background efficiency at optimal cut
        
    """
    
    # Validate input parameters
    if metric_type not in ['tpr_fpr', 'efficiency', 'significance', 'ratio']:
        raise ValueError(
            f"Invalid metric_type: {metric_type}. "
            "Must be one of: 'tpr_fpr', 'efficiency', 'significance', 'ratio'"
        )
    
    # Validate required arguments for significance metrics
    if metric_type in ['significance', 'ratio']:
        if have_sig_events is None or have_bg_events is None:
            raise ValueError(
                '"have_sig_events" and "have_bg_events" are required for '
                '"significance" and "ratio" metric types'
            )
            
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
            
            min_sig_sel_mask, best_arg, metric_1, metric_2, thresholds, sig_efficiency, bg_efficiency = find_optimal_cut_point(
                df=df,
                feature=feature,
                bounds=bounds,
                select_direction=direction,
                metric_type=metric_type,
                min_sig_sel=min_sig_sel,
                have_sig_events=have_sig_events,
                have_bg_events=have_bg_events,
                mass_interval=mass_interval
            )
        
            # Check if we have valid results
            if (not min_sig_sel_mask.any() or 
                best_arg >= len(thresholds[min_sig_sel_mask]) or
                len(metric_1[min_sig_sel_mask]) == 0):
                raise ValueError(f"No valid cuts found for {feature} with {direction} direction")
                continue
        
            best_cut_x = thresholds[min_sig_sel_mask][best_arg]
            best_metric_1 = metric_1[min_sig_sel_mask][best_arg]
            best_metric_2 = metric_2[min_sig_sel_mask][best_arg]
            
            if metric_type in ('tpr_fpr', 'efficiency'):
                metric_value = auc(metric_2[min_sig_sel_mask], metric_1[min_sig_sel_mask])
            elif metric_type == 'significance':
                metric_value = best_metric_1
            elif metric_type == 'ratio':
                metric_value = best_metric_2
                
            print(f'{metric_type} value: {metric_value:.4f}, direction: {direction}')
                
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
        f'metric': metric_value_list,
        f'feature': feature_list,
        f'feature_sel_direction': featuree_sel_direction_list,
        f'cut': cut_list,
        f'signal_efficiency': signal_efficiency_list,
        f'background_efficiency': background_efficiency_list
    }) 
    
    # Sort by importance metrics (descending order)
    importance_df = importance_df.sort_values(by=f'metric', ascending=False).reset_index(drop=True)
    
    return importance_df
    
    
def create_best_selection_path(
    df: pd.DataFrame, 
    features: List[str], 
    n_features_to_use: int, 
    metric_type: str, 
    have_sig_events: int, 
    have_bg_events: int, 
    mass_interval=(2.24763, 2.32497),
    direction_restrictions: Dict = None # Format: {'feature_name': 'left'/'right'}
) -> pd.DataFrame:
    
    """
    Create an optimal feature selection path by sequentially selecting the most important features.
    
    This function performs forward feature selection by iteratively:
    1. Calculating feature importance for all remaining features
    2. Selecting the most important feature based on the specified metric
    3. Applying the optimal cut to the data
    4. Repeating until the desired number of features is selected or no events remain
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing signal and background events. Must include:
        - 'tag': Event type ('Sig' or 'Bg')
        - 'mass_Lc': Mass values for signal estimation
        - All specified feature columns
    features : List[str]
        List of feature names to consider for selection
    n_features_to_use : int
        Maximum number of features to select. If None, uses all available features
    metric_type : str
        Metric to optimize for feature selection:
        - 'tpr_fpr': ROC AUC (Area Under Receiver Operating Characteristic curve)
        - 'efficiency': AUC of signal vs background efficiency curve
        - 'significance': Statistical significance measure
        - 'ratio': Signal-to-background ratio
    have_sig_events : int
        Total number of signal events in the dataset for proper normalization
    have_bg_events : int
        Total number of background events in the dataset for proper normalization
    mass_interval : Tuple[float, float], default=(2.24763, 2.32497)
        Mass interval used for signal estimation
    direction_restrictions : Dict[str, str], optional
        Dictionary specifying selection direction restrictions for specific features.
        Format: {'feature_name': 'left'/'right'} where:
        - 'left': select events where feature < cut
        - 'right': select events where feature > cut
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the sequential selection path with columns:
        - 'selected_feature': Feature selected at each step
        - 'cut_point': Optimal cut value applied
        - 'select_direction': Direction of selection ('left' or 'right')
        - 'signal_efficiency': Signal efficiency after selection
        - 'background_efficiency': Background efficiency after selection
        - 'ratio': Signal-to-background ratio
        - 'overall_s_sqrt_s_b': Statistical significance (S/√(S+B))
        - 'total_sig_unscaled': Unscaled signal count
        - 'total_bg_unscaled': Unscaled background count
        
    """
    
    # Input validation
    if metric_type not in ['tpr_fpr', 'efficiency', 'significance', 'ratio']:
        raise ValueError(
            f"Invalid metric_type: {metric_type}. "
            "Must be one of: 'tpr_fpr', 'efficiency', 'significance', 'ratio'"
        )
    
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
                have_sig_events=have_sig_events,
                have_bg_events=have_bg_events,
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


def fit_distr_triple_gauss(distr: Union[pd.Series, np.ndarray],
              distr_name: str,
              title: str,
              x_label: str,
              optimize_method: str = 'L-BFGS-B',
              plot: bool = False,
              bins: int = 40,
              max_iter: int=1000) -> Tuple:

    data = np.asarray(distr)

    counts, bin_edges = np.histogram(data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    errors = np.sqrt(counts)
    
    # Filter division by zero
    mask = counts > 0
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    errors = errors[mask]

    #=======================================================================================
    # Fit Model
    #=======================================================================================
    
    def triple_gaussian(x, *params):
        """
        Composite model of 3 Gaussian functions
        
        Args:
            x: array-like, input values
            params: 9 parameters [A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3]
        
        Returns:
            y: sum of 3 Gaussian functions
        """
        A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = params
        
        g1 = A1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
        g2 = A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
        g3 = A3 * np.exp(-0.5 * ((x - mu3) / sigma3) ** 2)
        
        return g1 + g2 + g3

    # Better Initial Guess
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])  # Scale by bin width
    
    # Use percentiles to ensure means are within bin range
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    # Spread initial means across the data range
    initial_guess = [
        initial_amplitude, 0, 0.1 * data_range,  # Gaussian 1
        initial_amplitude, 0, 0.1 * data_range,  # Gaussian 2  
        initial_amplitude, 0, 0.1 * data_range   # Gaussian 3
    ]
    
    # Set bounds for parameters
    bounds = [
        (0, None), (-0.01, 0.01), (1e-6, 0.5 * data_range),  # Gaussian 1
        (0, None), (-0.01, 0.01), (1e-6, 0.5 * data_range),  # Gaussian 2
        (0, None), (-0.01, 0.01), (1e-6, 0.5 * data_range)   # Gaussian 3
    ]
    
    def calculate_effective_sigma_triple_gaussian(fit_params):
        """
        Calculate effective sigma for triple Gaussian mixture
        
        For a mixture distribution: f(x) = w1*g1(x) + w2*g2(x) + w3*g3(x)
        where w_i = A_i / (A1 + A2 + A3) are the weights
        
        The variance of the mixture is:
        σ²_eff = Σ w_i * (σ_i² + μ_i²) - (Σ w_i * μ_i)²
        """
        A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = fit_params
        
        # Calculate weights (normalized by total amplitude/area)
        total_amplitude = A1 + A2 + A3
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
        w3 = A3 / total_amplitude
        
        # Calculate effective mean
        effective_mean = w1 * mu1 + w2 * mu2 + w3 * mu3
        
        # Calculate effective variance
        # effective_variance = (w1 * (sigma1**2 + mu1**2) + 
        #                      w2 * (sigma2**2 + mu2**2) + 
        #                      w3 * (sigma3**2 + mu3**2) - 
        #                      effective_mean**2)
        
        effective_variance = w1 * sigma1**2 + w2 * sigma2**2 + w3 * sigma3**2
        
        effective_sigma = np.sqrt(effective_variance)
        
        return effective_sigma, effective_mean
    
    #=======================================================================================
    #
    #=======================================================================================
    
    # Loss Function
    def weighted_loss(params, x, y, errors):
        """Objective function to minimize (sum of weighted squared residuals)"""
        predictions = triple_gaussian(x, *params)
        residuals = y - predictions
        weights = 1.0 / (errors**2 + 1e-6)
        return np.sum(weights * residuals**2)
    
    # Fit
    result = minimize(
        fun=weighted_loss,
        x0=initial_guess, 
        args=(bin_centers, counts, errors), 
        method=optimize_method,
        bounds=bounds,
        options={'maxiter': max_iter}
    )

    if not result.success:
        print('Fit problems!')
        # raise ValueError('Fit problems!')

    fit_params = result.x

    effective_sigma, effective_mean = calculate_effective_sigma_triple_gaussian(fit_params)

    print("Optimization result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of iterations: {result.nit}")
    print(f"Final objective value: {result.fun:.6f}")
    print(f"Gauss 1 parameters: A={fit_params[0]}, mean={fit_params[1]}, sigma={fit_params[2]}")
    print(f"Gauss 2 parameters: A={fit_params[3]}, mean={fit_params[4]}, sigma={fit_params[5]}")
    print(f"Gauss 3 parameters: A={fit_params[6]}, mean={fit_params[7]}, sigma={fit_params[8]}")
    print(f'Effective Mean: {effective_mean}')
    print(f'Effective Sigma: {effective_sigma}')

    # Draw Data and Fit Curve
    if plot:
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [2, 1]})
        
        if result.success:
            
            # Calculate fit values and pulls
            fit_values = triple_gaussian(bin_centers, *fit_params)
            residuals = counts - fit_values
            pulls = residuals / errors  # (data - fit) / error
            
            # Enable LaTeX
            plt.rcParams.update({
                "text.usetex": True,
                "font.family": "serif",
                "font.size": 12
            })
            
            # Upper plot: Data and fit
            axes[0].errorbar(
                x=bin_centers,
                y=counts,
                yerr=errors, 
                fmt='o',
                markersize=5,
                capsize=3,
                capthick=1, 
                color='black',
                ecolor='gray',
                alpha=0.7,
                linewidth=1,
                label='Data'
            )
            
            # Plot the fitted curve
            x_fit = np.linspace(bin_edges[0], bin_edges[-1], 200)
            y_fit = triple_gaussian(x_fit, *fit_params)
            
            axes[0].plot(
                x_fit,
                y_fit,
                color='red',
                linestyle='-',
                linewidth=2, 
                label='Triple Gaussian Fit'
            )
            
            # Plot individual Gaussians
            g1 = fit_params[0] * np.exp(-0.5 * ((x_fit - fit_params[1]) / fit_params[2]) ** 2)
            g2 = fit_params[3] * np.exp(-0.5 * ((x_fit - fit_params[4]) / fit_params[5]) ** 2)
            g3 = fit_params[6] * np.exp(-0.5 * ((x_fit - fit_params[7]) / fit_params[8]) ** 2)
            
            axes[0].plot(x_fit, g1, 'g:', label=fr'Gaussian 1: $\mu$={fit_params[1] * 1e4:.2f} $\mu m$', alpha=0.7)
            axes[0].plot(x_fit, g2, 'b:', label=fr'Gaussian 2: $\mu$={fit_params[4] * 1e4:.2f} $\mu m$', alpha=0.7)
            axes[0].plot(x_fit, g3, 'm:', label=fr'Gaussian 3: $\mu$={fit_params[7] * 1e4:.2f} $\mu m$', alpha=0.7)
            
            # Add fit parameters as text
            fit_text = (fr'Fit Parameters:'
                        fr'\\G1: A={fit_params[0]:.1f}, $\mu$={fit_params[1] * 1e4:.3f} $\mu m$, $\sigma$={fit_params[2] * 1e4:.3f} $\mu m$'
                        fr'\\G2: A={fit_params[3]:.1f}, $\mu$={fit_params[4] * 1e4:.3f} $\mu m$, $\sigma$={fit_params[5] * 1e4:.3f} $\mu m$'
                        fr'\\G3: A={fit_params[6]:.1f}, $\mu$={fit_params[7] * 1e4:.3f} $\mu m$, $\sigma$={fit_params[8] * 1e4:.3f} $\mu m$'
                        fr'\\Effective $\sigma$: {effective_sigma * 1e4: .4f} $\mu m$')
            
            axes[0].text(0.02, 0.98, fit_text, transform=axes[0].transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)
            
            axes[0].set_title(f'{title}', fontsize=14, fontweight='bold', pad=20)
            axes[0].set_xlabel(x_label, fontweight='bold')
            axes[0].set_ylabel('Counts', fontweight='bold')
            axes[0].legend(loc='best')
            axes[0].grid(True, alpha=0.3, linestyle='--')
            
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            
            # Lower plot: Pulls
            axes[1].axhline(y=0, color='red', linestyle='-', linewidth=1, alpha=0.7)
            axes[1].axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=-1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            axes[1].axhline(y=-2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            
            bin_width = bin_edges[1] - bin_edges[0]
            # axes[1].scatter(bin_centers, pulls, s=12, c='blue', marker='o', alpha=0.7, label='Pulls')
            # axes[1].bar(bin_centers, pulls, width=bin_width, color='blue', alpha=0.7, label='Pulls')
            axes[1].bar(bin_centers, pulls, 
            width=bin_width*0.8,  # 80% of bin width for spacing
            color='skyblue', 
            edgecolor='darkblue',
            alpha=0.7, 
            label='Pulls')
            
            
            # Calculate pull statistics
            mean_pull = np.mean(pulls)
            std_pull = np.std(pulls)
            chi2 = np.sum(pulls**2)
            ndf = len(pulls) - 9  # number of bins - number of parameters (9 for triple Gaussian)
            
            pull_stats = (fr'Pull Statistics:'
                          fr'\\$\mu_{{pull}}$ = {mean_pull:.3f}'
                          fr'\\$\sigma_{{pull}}$ = {std_pull:.3f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{ndf}')
            
            axes[1].text(0.02, 0.98, pull_stats, transform=axes[1].transAxes,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8), fontsize=9)
            
            axes[1].set_ylabel('(data - fit) / sqrt(data) at bin', fontweight='bold')
            axes[1].legend(loc='upper right')
            
            axes[1].set_xlabel(f'{x_label}', fontweight='bold')
            axes[1].grid(True, alpha=0.3, linestyle='--')
            
            axes[1].spines['top'].set_visible(False)
            axes[1].spines['right'].set_visible(False)
        
        else:
            axes[0].text(0.5, 0.5, 'Fit failed!', transform=axes[0].transAxes, 
                    ha='center', va='center', fontsize=12, color='red')    
        
        plt.tight_layout()
        
        # Disable LaTeX
        plt.rcParams.update({
            "text.usetex": False
        })
        
        return result, fig, axes
    
    return result


def fit_distr_double_gauss(distr: Union[pd.Series, np.ndarray],
              distr_name: str,
              title: str,
              x_label: str,
              optimize_method: str = 'L-BFGS-B',
              plot: bool = False,
              bins: int = 40,
              max_iter: int = 1000) -> Tuple:

    data = np.asarray(distr)

    counts, bin_edges = np.histogram(data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    errors = np.sqrt(counts)
    
    # Filter division by zero
    mask = counts > 0
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    errors = errors[mask]

    #=======================================================================================
    # Fit Model - Double Gaussian with Constant
    #=======================================================================================
    
    def double_gaussian(x, *params):
        """
        Composite model of 2 Gaussian functions with constant background
        
        Args:
            x: array-like, input values
            params: 7 parameters [A1, mu1, sigma1, A2, mu2, sigma2, C]
        
        Returns:
            y: sum of 2 Gaussian functions and Constant
        """
        A1, mu1, sigma1, A2, mu2, sigma2, C = params
        
        g1 = A1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
        g2 = A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
        
        return g1 + g2 + C
    
    # Better Initial Guess
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])  # Scale by bin width
    
    # Use percentiles to ensure means are within bin range
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    # Spread initial means across the data range
    initial_guess = [
        initial_amplitude, data_min + 0.25 * data_range, 0.1 * data_range,  # Gaussian 1
        initial_amplitude, data_min + 0.75 * data_range, 0.1 * data_range,  # Gaussian 2  
        0   # Constant
    ]
    
    # Set bounds for parameters
    bounds = [
        (0, None), (data_min, data_max), (1e-6, 0.5 * data_range),  # Gaussian 1
        (0, None), (data_min, data_max), (1e-6, 0.5 * data_range),  # Gaussian 2
        (0, None)   # Constant
    ]
    
    
    def calculate_effective_sigma_double_gaussian(fit_params):
        """
        Calculate effective sigma for triple Gaussian mixture
        
        For a mixture distribution: f(x) = w1*g1(x) + w2*g2(x) + w3*g3(x)
        where w_i = A_i / (A1 + A2 + A3) are the weights
        
        The variance of the mixture is:
        σ²_eff = Σ w_i * (σ_i² + μ_i²) - (Σ w_i * μ_i)²
        """
        A1, mu1, sigma1, A2, mu2, sigma2, C = fit_params
        
        # Calculate weights (normalized by total amplitude/area)
        total_amplitude = A1 + A2
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
        
        # Calculate effective mean
        effective_mean = w1 * mu1 + w2 * mu2
        
        # Calculate effective variance
        # effective_variance = (w1 * (sigma1**2 + mu1**2) + 
        #                      w2 * (sigma2**2 + mu2**2) + 
        #                      w3 * (sigma3**2 + mu3**2) - 
        #                      effective_mean**2)
        
        effective_variance = w1 * sigma1**2 + w2 * sigma2**2
        
        effective_sigma = np.sqrt(effective_variance)
        
        return effective_sigma, effective_mean
    
    #=======================================================================================
    # Loss Function and Optimization
    #=======================================================================================
    
    def weighted_loss(params, x, y, errors):
        """Objective function to minimize (sum of weighted squared residuals)"""
        predictions = double_gaussian(x, *params)
        residuals = y - predictions
        weights = 1.0 / (errors**2 + 1e-6)
        return np.sum(weights * residuals**2)
    
    # Fit
    result = minimize(
        fun=weighted_loss,
        x0=initial_guess, 
        args=(bin_centers, counts, errors), 
        method=optimize_method,
        bounds=bounds,
        options={'maxiter': max_iter}
    )

    if not result.success:
        print('Fit problems!')
        # raise ValueError('Fit problems!')

    fit_params = result.x

    effective_sigma, effective_mean = calculate_effective_sigma_double_gaussian(fit_params)

    print("Optimization result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of iterations: {result.nit}")
    print(f"Final objective value: {result.fun:.6f}")
    print(f"Gauss 1 parameters: A={fit_params[0]:.4f}, mean={fit_params[1]:.6f}, sigma={fit_params[2]:.6f}")
    print(f"Gauss 2 parameters: A={fit_params[3]:.4f}, mean={fit_params[4]:.6f}, sigma={fit_params[5]:.6f}")
    print(f"Constant: C={fit_params[6]:.4f}")
    print(f'Effective Mean: {effective_mean:.6f}')
    print(f'Effective Sigma: {effective_sigma:.6f}')

    # Draw Data and Fit Curve
    if plot:
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [2, 1]})
        
        if result.success:
            
            # Calculate fit values and pulls
            fit_values = double_gaussian(bin_centers, *fit_params)
            residuals = counts - fit_values
            pulls = residuals / errors  # (data - fit) / error
            
            # Enable LaTeX
            plt.rcParams.update({
                "text.usetex": True,
                "font.family": "serif",
                "font.size": 12
            })
            
            # Upper plot: Data and fit
            axes[0].errorbar(
                x=bin_centers,
                y=counts,
                yerr=errors, 
                fmt='o',
                markersize=5,
                capsize=3,
                capthick=1, 
                color='black',
                ecolor='gray',
                alpha=0.7,
                linewidth=1,
                label='Data'
            )
            
            # Plot the fitted curve
            x_fit = np.linspace(bin_edges[0], bin_edges[-1], 200)
            y_fit = double_gaussian(x_fit, *fit_params)
            
            axes[0].plot(
                x_fit,
                y_fit,
                color='red',
                linestyle='-',
                linewidth=2, 
                label='Triple Gaussian Fit'
            )
            
            # Plot individual Gaussians
            g1 = fit_params[0] * np.exp(-0.5 * ((x_fit - fit_params[1]) / fit_params[2]) ** 2)
            g2 = fit_params[3] * np.exp(-0.5 * ((x_fit - fit_params[4]) / fit_params[5]) ** 2)
            C = fit_params[6]
            
            axes[0].plot(x_fit, g1, 'g:', label=fr'Gaussian 1: $\mu$={fit_params[1] * 1e4:.2f} $\mu m$', alpha=0.7)
            axes[0].plot(x_fit, g2, 'b:', label=fr'Gaussian 2: $\mu$={fit_params[4] * 1e4:.2f} $\mu m$', alpha=0.7)
            axes[0].plot(x_fit, C * np.ones(shape=g1.size), 'm:', label=fr'Const: C={fit_params[6]:.2f}', alpha=0.7)
            
            # Add fit parameters as text
            fit_text = (fr'Fit Parameters:'
                        fr'\\G1: A={fit_params[0]:.1f}, $\mu$={fit_params[1] * 1e4:.3f} $\mu m$, $\sigma$={fit_params[2] * 1e4:.3f} $\mu m$'
                        fr'\\G2: A={fit_params[3]:.1f}, $\mu$={fit_params[4] * 1e4:.3f} $\mu m$, $\sigma$={fit_params[5] * 1e4:.3f} $\mu m$'
                        fr'\\Const: C={fit_params[6]:.1f}'
                        fr'\\Effective $\sigma$: {effective_sigma * 1e4: .4f} $\mu m$')
            
            axes[0].text(0.02, 0.98, fit_text, transform=axes[0].transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)
            
            axes[0].set_title(f'{title}', fontsize=14, fontweight='bold', pad=20)
            axes[0].set_xlabel(x_label, fontweight='bold')
            axes[0].set_ylabel('Counts', fontweight='bold')
            axes[0].legend(loc='best')
            axes[0].grid(True, alpha=0.3, linestyle='--')
            
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            
            # Lower plot: Pulls
            axes[1].axhline(y=0, color='red', linestyle='-', linewidth=1, alpha=0.7)
            axes[1].axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=-1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            axes[1].axhline(y=-2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            
            bin_width = bin_edges[1] - bin_edges[0]
            # axes[1].scatter(bin_centers, pulls, s=12, c='blue', marker='o', alpha=0.7, label='Pulls')
            # axes[1].bar(bin_centers, pulls, width=bin_width, color='blue', alpha=0.7, label='Pulls')
            axes[1].bar(bin_centers, pulls, 
            width=bin_width*0.8,  # 80% of bin width for spacing
            color='skyblue', 
            edgecolor='darkblue',
            alpha=0.7, 
            label='Pulls')
            
            
            # Calculate pull statistics
            mean_pull = np.mean(pulls)
            std_pull = np.std(pulls)
            chi2 = np.sum(pulls**2)
            ndf = len(pulls) - 9  # number of bins - number of parameters (9 for triple Gaussian)
            
            pull_stats = (fr'Pull Statistics:'
                          fr'\\$\mu_{{pull}}$ = {mean_pull:.3f}'
                          fr'\\$\sigma_{{pull}}$ = {std_pull:.3f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{ndf}')
            
            axes[1].text(0.02, 0.98, pull_stats, transform=axes[1].transAxes,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8), fontsize=9)
            
            axes[1].set_ylabel('(data - fit) / sqrt(data) at bin', fontweight='bold')
            axes[1].legend(loc='upper right')
            
            axes[1].set_xlabel(f'{x_label}', fontweight='bold')
            axes[1].grid(True, alpha=0.3, linestyle='--')
            
            axes[1].spines['top'].set_visible(False)
            axes[1].spines['right'].set_visible(False)
        
        else:
            axes[0].text(0.5, 0.5, 'Fit failed!', transform=axes[0].transAxes, 
                    ha='center', va='center', fontsize=12, color='red')    
        
        plt.tight_layout()
        
        # Disable LaTeX
        plt.rcParams.update({
            "text.usetex": False
        })
        
        return result, fig, axes
    
    return result


def classification_quality(
    df: pd.DataFrame,
    visualize: bool = False,
    title: str = None
) -> pd.DataFrame:
    """
    Calculate and optionally visualize classification rates for particle identification.
    
    This function computes the fraction of correctly and incorrectly classified particles
    by comparing Monte Carlo truth PID (mc_pid) with predicted PID (pred_pid). It returns
    classification rates and can generate both bar plots and confusion matrix heatmaps.

    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe containing 'mc_pid' and 'pred_pid' columns
    visualize : bool, default=False
        If True, generates bar plot and confusion matrix visualization
    title : str, optional
        Title to be added to the plots when visualize=True
    particle_order : list, default=[2212, 321, 211, 0]
        Order of particles for visualization (proton, kaon, pion, other)
    particle_labels : dict, optional
        Mapping from particle IDs to labels. Default: {2212: 'proton', 321: 'kaon', 211: 'pion', 0: 'other'}
    color_palette : list, optional
        Color palette for the bar plot. Default: ['#E6AB02', '#1E66A6', '#29E786', '#FFCDA4']
    figsize : tuple, default=(12, 5)
        Figure size when visualize=True

    Returns:
    --------
    pandas.DataFrame
        DataFrame with columns:
        - mc_pid: Monte Carlo truth particle ID
        - pred_pid: Predicted particle ID  
        - fraction: Classification rate (normalized count)
        - count: Raw count of each classification

    """
    
    # Input validation
    required_columns = ['mc_pid', 'pred_pid']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    df = df.copy()

    particle_labels = {2212: 'proton', 321: 'kaon', 211: 'pion', 0: 'other'}
    color_palette = ['#E6AB02', '#1E66A6', '#29E786', '#FFCDA4']

    # Calculate classification rates with both normalized and raw counts
    counts_df = df.groupby(['mc_pid', 'pred_pid']).size().reset_index(name='count')
    # fractions_df = df.groupby('mc_pid')['pred_pid'].value_counts(normalize=True).reset_index(name='fraction')
    fractions_df = df.groupby('mc_pid').value_counts(normalize=True).reset_index(name='fraction')
    
    # Merge counts and fractions
    res = pd.merge(fractions_df, counts_df, on=['mc_pid', 'pred_pid'])

    # Create comparison column for correct/incorrect classification
    df['comparison'] = (df['mc_pid'] == df['pred_pid']).astype(np.int32)

    # Create confusion matrix as pivot table
    piv_table = pd.pivot_table(
        data=df,
        values='comparison',
        columns='pred_pid',
        index='mc_pid',
        aggfunc='count',
        fill_value=0
    )

    if visualize:
        
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
        
        # Bar plot showing classification rates
        sns.barplot(
            data=res,
            x='pred_pid',
            y='fraction',
            hue='mc_pid',
            hue_order=[2212, 321, 211, 0],
            palette=color_palette,
            ax=axes[0],
            width=0.5,
            alpha=0.8
        )

        max_value = res['fraction'].max()
        axes[0].set_ylim(0, max_value * 1.2)

        for p in axes[0].patches:
            height = p.get_height()
            if height > 0.01:
                axes[0].text(
                    p.get_x() + p.get_width()/2., 
                    height + height * 0.01,
                    f'{height:.2f}', 
                    ha='center', 
                    va='bottom', 
                    fontsize=6
                )

        pred_labels = [particle_labels.get(pid, str(pid)) for pid in sorted(res['pred_pid'].unique())]
        axes[0].set_xticklabels([f'as {label}' for label in pred_labels])

        plot_title = f'Classification Rate {title}' if title else 'Classification Rate'
        axes[0].set_title(plot_title, fontsize=14, fontweight='bold', pad=20)
        axes[0].set_xlabel('PID Predictions', fontweight='bold')
        axes[0].set_ylabel('Classification Rate', fontweight='bold')

        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)
        axes[0].grid(True, alpha=0.3, linestyle='--')

        legend_labels = [particle_labels.get(pid, str(pid)) for pid in [2212, 321, 211, 0]]
        sns.move_legend(
            axes[0],
            loc='best',
            title='MC PID',
            labels=legend_labels,
            frameon=True,
            fancybox=True,
            shadow=True
        )

        # Heatmap visualization of confusion matrix
        if not piv_table.empty:
            
            sns.heatmap(
                data=piv_table,
                cmap='Blues',
                annot=True,
                fmt='.0f',
                annot_kws={'fontsize': 10},
                linewidths=1,
                linecolor='gray',
                cbar_kws={'label': 'Number of predictions'},
                center=piv_table.values.mean() if piv_table.size > 0 else 0,
                square=True,
                ax=axes[1]
            )

            y_labels = [particle_labels.get(pid, str(pid)) for pid in piv_table.index]
            x_labels = [particle_labels.get(pid, str(pid)) for pid in piv_table.columns]
            axes[1].set_yticklabels(y_labels)
            axes[1].set_xticklabels(x_labels)

            axes[1].set_title(f'Confusion Matrix {title}', fontsize=14, pad=25, fontweight='bold')
            axes[1].set_xlabel('PID Predictions', fontsize=12, fontweight='bold')
            axes[1].set_ylabel('MC Truth PID', fontsize=12, fontweight='bold')
        else:
            axes[1].text(0.5, 0.5, 'No data available\nfor confusion matrix', 
                        ha='center', va='center', transform=axes[1].transAxes)
            axes[1].set_title('Confusion Matrix', fontsize=14, pad=25, fontweight='bold')

        plt.tight_layout()
        plt.show()
        
        return res, fig

    return res, None


def choice(
    df: pd.DataFrame,
    thresholds: list = None
):
    
    """
    Classifies particle types based on probability thresholds and maximum probability selection.

    This function processes a DataFrame containing particle identification probabilities,
    optionally applies threshold filtering, and selects the most probable particle type
    for each entry based on the highest probability among pion (211), kaon (321), and proton (2212).

    Args:
        ll_df (pd.DataFrame): DataFrame containing particle probabilities with columns:
                            [mc_pid, pip_prob, K_prob, p_prob] where:
                            - mc_pid: Monte Carlo truth particle ID
                            - pip_prob: Pion probability
                            - K_prob: Kaon probability  
                            - p_prob: Proton probability
        thresholds (list, optional): List of three threshold values [pip_thresh, K_thresh, p_thresh].
                                   If provided, filters rows where any probability exceeds its threshold
                                   and subtracts thresholds from probabilities. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame with two columns:
                     - mc_pid: Original Monte Carlo truth particle ID
                     - pred_pid: Predicted particle type (211, 321, or 2212) based on maximum probability

    """
    
    df = df.copy()
    
    cols_name = df.columns

    if thresholds is not None:
        
        # Create threshold mapping for each particle type
        mask_thresh = {
            0: thresholds[0],  # Pion threshold
            1: thresholds[1],  # Kaon threshold  
            2: thresholds[2]   # Proton threshold
        }

        # Filter rows where any probability exceeds its threshold
        thresh_mask = ((df[cols_name[1]] > mask_thresh[0]) | 
                      (df[cols_name[2]] > mask_thresh[1]) | 
                      (df[cols_name[3]] > mask_thresh[2]))

        df = df[thresh_mask].reset_index(drop=True)

        # Subtract thresholds from probabilities (normalization step)
        df.iloc[:, 1] -= mask_thresh[0]
        df.iloc[:, 2] -= mask_thresh[1]
        df.iloc[:, 3] -= mask_thresh[2]

    df.columns = ['mc_pid', 211, 321, 2212]  # 211: pion, 321: kaon, 2212: proton

    # Select particle type with highest probability for each row
    res_series = df[[211, 321, 2212]].idxmax(axis=1).rename('pred_pid')

    # Combine original MC truth with predictions
    res_df = pd.concat([df.mc_pid, res_series], axis=1)

    return res_df


def check_thresholds(
    df: pd.DataFrame,
    thresholds: list
):
    """
    Evaluate classification performance across different probability thresholds.
    
    This function analyzes how different probability thresholds affect:
    - Track selection efficiency (what percentage of tracks pass the thresholds)
    - Classification accuracy for each particle type (pi, K, p)
    
    Args:
        df (pd.DataFrame): Input dataframe containing particle classification data.
                          Expected columns: 'mc_pid' (true particle ID), 
                          'pred_pid' (predicted particle ID), and probability columns
                          for each particle type.
        thresholds (list): List of threshold tuples, where each tuple contains 
                          three thresholds for [pi, K, p] probabilities.
    
    Returns:
        None
    """
    
    track_part = []  # Track selection efficiency
    pi_as_pi = []    # pi classification accuracy
    k_as_k = []      # kaon classification accuracy  
    p_as_p = []      # proton classification accuracy
    
    for thresh in thresholds:
        
        input_size = df.shape[0]
        
        selected_df = choice(df, [thresh[0], thresh[1], thresh[2]])
                
        output_size = selected_df.shape[0]
        
        # Calculate classification rates by true particle type
        rate_df = selected_df.groupby('mc_pid').value_counts(normalize=True).reset_index(name='fraction')
        
        # Extract correctly classified samples (where true ID matches predicted ID)
        corrected_class = rate_df[rate_df.mc_pid == rate_df.pred_pid]
        
        track_part.append(output_size / input_size)  # Tracks Selection Efficiency
        
        # Assuming fixed order: [pi, K, p]
        if len(corrected_class) >= 3:
            pi_as_pi.append(corrected_class.fraction.iloc[0])
            k_as_k.append(corrected_class.fraction.iloc[1]) 
            p_as_p.append(corrected_class.fraction.iloc[2])
        else:
            pi_as_pi.append(0)
            k_as_k.append(0)
            p_as_p.append(0)
    
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    sns.lineplot(x=np.arange(len(thresholds)), y=track_part, label='Tracks efficiency', marker='o', ax=ax)
    sns.lineplot(x=np.arange(len(thresholds)), y=pi_as_pi, label='pi as pi', marker='o', ax=ax)
    sns.lineplot(x=np.arange(len(thresholds)), y=k_as_k, label='k as k', marker='o', ax=ax) 
    sns.lineplot(x=np.arange(len(thresholds)), y=p_as_p, label='p as p', marker='o', ax=ax)
    
    ax.grid(visible=True, alpha=0.3, linestyle='--')
    ax.set_xlabel('Threshold Set')
    ax.set_ylabel('Classification Rate')
    ax.set_title('Classification Performance vs. Probability Thresholds')
    
    ax.set_xticks(np.arange(len(thresholds)))
    ax.set_xticklabels([str(thresh) for thresh in thresholds], rotation=45)
    
    ax.legend()
    plt.tight_layout()
    
    
def bar_plot(
    data,
    x_data: str,
    y_data: str,
    xticklabels: list,
    title,
    xlabel,
    ylabel,
    color='#4C72B0',
    width=0.5,
    alpha=0.8,
):
    
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0, 0, 1, 1])

    sns.barplot(
        data=data,
        x=x_data,
        y=y_data,
        color='#4C72B0',
        width=0.5,
        alpha=0.8,
        ax=ax
    )

    max_value = data['proportion'].max()
    ax.set_ylim(0, max_value * 1.2)

    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.text(p.get_x() + p.get_width()/2., height + height * 0.001,
                    f'{height: 0.2f}%', ha='center', va='bottom', fontsize=10)

    ax.set_xticklabels(xticklabels)

    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(True, alpha=0.8, linestyle='--')

    plt.show()
    
    
def create_subplots_dynamic(ncols, df, vars: list, bins: list = None, plot_type='hist', hue=False, norm=False):

    temp_df = df[['mc_pid'] + vars].copy()

    nrows = (temp_df.shape[1] + ncols - 1) // ncols  # Calculate rows needed

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5*ncols, 4*nrows))
    fig.suptitle('Initial track parameters distributions', fontsize=16, fontweight='bold')

    # Flatten axes for easy iteration if it's a 2D array
    if nrows > 1 or ncols > 1:
        axes_flat = axes.flatten()
    else:
        axes_flat = [axes]

    # Plot data
    for i, var in enumerate(vars):
        if i < len(axes_flat):

            labels = {
                'p': 'p (GeV)',
                'pt': 'pt (GeV)',
                'px': 'px (GeV)',
                'py': 'py (GeV)',
                'pz': 'pz (GeV)',
                'eta': 'eta',
                'theta': 'theta (rad)',
                'phi': 'phi (rad)',
                'chi2overndf': 'Chi2 ver NDF',
                'nhitsits': 'N Hits Its',
                'nhitstsb': 'N Hits TsB',
                'nhitstsec': 'N Hits TsEc',
                'nhitstsb + nhitstsec': 'N Hits TsB + N Hits TsEc',
                'nhitsits + nhitstsb + nhitstsec': 'N Hits Its + N Hits TsB + N Hits TsEc',
                'isfitted': 'IsFitted',
                'isgood': 'IsGood',
                'convergency': 'Convergency'
            }

            if plot_type == 'hist':
                sns.histplot(
                    data=temp_df,
                    x=var,
                    stat='density' if norm else 'count',
                    common_norm=False,
                    hue='mc_pid' if hue else None,
                    hue_order=[0, 2212, 321, 211] if hue and not norm else None,
                    palette='Set2' if hue else None,
                    alpha=0.8,
                    bins=bins[i] if bins is not None else 'auto',
                    ax=axes_flat[i]
                )
            elif plot_type == 'count':
                sns.countplot(
                    data=temp_df,
                    x=var,
                    alpha=0.8,
                    hue='mc_pid' if hue else None,
                    palette='Set2' if hue else None,
                    ax=axes_flat[i]
                )
            elif plot_type == 'box':
                sns.boxplot(
                    data=temp_df,
                    x=var,
                    hue='mc_pid' if hue else None,
                    palette='Set2' if hue else None,
                    width=0.6,
                    linewidth=1.5,
                    fliersize=4,
                    ax=axes_flat[i]
                )


            axes_flat[i].set_xlabel(f'{labels[var]}', fontweight='bold')
            axes_flat[i].set_ylabel(f'Density' if norm else f'Counts' , fontweight='bold')

            axes_flat[i].spines['top'].set_visible(False)
            axes_flat[i].spines['right'].set_visible(False)

            axes_flat[i].grid(True, alpha=0.8, linestyle='--')


    # Hide unused subplots
    for i in range(len(vars), len(axes_flat)):
        axes_flat[i].set_visible(False)


    plt.tight_layout()
    return fig, axes