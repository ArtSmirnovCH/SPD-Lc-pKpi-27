import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import erf
from typing import Tuple, Union
import sys
import os
from math import ceil
from scipy.stats import beta


# Path setup
try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.getcwd()
    
parent_dir = os.path.join(current_dir, '..', '..')
sys.path.insert(0, parent_dir)

from analysis_scripts.fit_scripts import fit_distr_flat_bg
from analysis_scripts.config import CFG

###########################################################################################################
# signal_estimates
###########################################################################################################

def signal_estimates(sig_mass_distr: Union[pd.Series, np.ndarray], 
                     bg_mass_distr: Union[pd.Series, np.ndarray], 
                     have_sig_events: int,
                     have_bg_events: int,
                     mass_interval: Tuple[float, float] = (2.24763, 2.32497),
                     total_spectrum: bool = False,
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
     # Physical constants
    sigma_MB = 38 * 1e-3 * 1e-24                # Minimum bias cross-section in cm2
    cross_section_Lc_plus = 4 * 1e-6 * 1e-24    # Lc+(4122) cross-section in cm2
    L = 1e32                                    # Luminosity in cm-2s-1
    t = 31536000                                # Time period in seconds (1 year)
    branching = 0.035                           # Lc decay branching ratio
    
    # Calculate expected events number
    real_sig_events = L * t * cross_section_Lc_plus * branching
    real_bg_events = L * t * sigma_MB
    
    # Calculate scaling weights to match real event expectations
    sig_weight = real_sig_events / have_sig_events
    bg_weight = real_bg_events / have_bg_events
    
    sig_mass_distr = np.asarray(sig_mass_distr)
    bg_mass_distr = np.asarray(bg_mass_distr)
    
    # Create weight arrays for signal and background
    weights_sig = np.full_like(sig_mass_distr, sig_weight, dtype=float)
    weights_bg = np.full_like(bg_mass_distr, bg_weight, dtype=float)
    
    mass_combined = np.concatenate([sig_mass_distr, bg_mass_distr])
    weights_combined = np.concatenate([weights_sig, weights_bg])
    labels_combined = np.concatenate([np.full(len(sig_mass_distr), 'Signal'), 
                                      np.full(len(bg_mass_distr), 'Background')])

    mass_df = pd.DataFrame({
        'mass': mass_combined,
        'weights': weights_combined,
        'type': labels_combined,
    })

    mass_mask = (mass_df['mass'] >= mass_interval[0]) & (mass_df['mass'] <= mass_interval[1])
    mass_df = mass_df[mass_mask]

    bins = 100
    binrange = (mass_df['mass'].min(), mass_df['mass'].max())

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

    total_heights = sig_heights + bg_heights

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

    sigma_sig = sig_weight * np.sqrt(np.sum(sig_mass_distr))
    sigma_bg = bg_weight * np.sqrt(np.sum(bg_mass_distr))
    sigma_part_sig = sigma_sig / total_bg
    sigma_part_bg = sigma_bg * total_sig / total_bg / total_bg
    
    s_b_uncertainty = np.sqrt(sigma_part_sig**2 + sigma_part_bg**2)

    if visualization:
        fig, axes = plt.subplots(3, 1, figsize=(6, 8), sharex=False, gridspec_kw={'height_ratios': [3, 1, 1]})

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

        if total_spectrum:
            axes[0].step(
                bin_edges[:-1], 
                total_heights,
                where='post',
                color='#2ca02c',
                linewidth=2.5,
                linestyle='--',
                label='Total (S+B)'
            )
               

        axes[0].set_title('Yield Mass Distribution', fontsize=16, fontweight='bold', pad=20)
        axes[0].set_xlabel('Mass (GeV)', fontweight='bold', fontsize=12)
        axes[0].set_ylabel('Counts (log scale)', fontweight='bold', fontsize=12)

        sns.move_legend(
            axes[0],
            loc='best',
            title=f'S/B: {overall_s_b:.6f}',
            # title=f'S/B: {overall_s_b:.6f}\nS/√(S+B): {overall_s_sqrt_s_b:.2f}',
            labels=['Signal', 'Background'],
            frameon=True,
            fancybox=True,
            shadow=True
        )

        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)
        axes[0].grid(True, alpha=0.8, linestyle='--')

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
        axes[1].grid(True, alpha=0.8, linestyle='--')
        axes[1].legend(loc='upper right', frameon=True, fancybox=True, shadow=True)

        # Plot 3: S/sqrt(S+B) Ratio
        sns.lineplot(
            x=bin_centers,
            y=s_over_sqrt_s_b,
            color='#9467bd',
            linewidth=2.5,
            label='S/sqrt(S+B) Ratio',
            ax=axes[2]
        )

        axes[2].set_ylabel('S/√(S+B)', fontweight='bold', fontsize=12)
        axes[2].set_xlabel('Mass (GeV)', fontweight='bold', fontsize=12)
        axes[2].grid(True, alpha=0.8, linestyle='--')
        
        plt.tight_layout()

    if verbose:
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)
        print(f"Total Signal Events: {total_sig:.0f}")
        print(f"Total Background Events: {total_bg:.0f}")
        print(f"Total Signal Events Unscaled: {total_sig_unscaled:.0f}")
        print(f"Total Background Events Unscaled: {total_bg_unscaled:.0f}")
        print(f"Overall S/B Ratio: {overall_s_b:.6f}")
        print(f's_b_uncertainty: {s_b_uncertainty}')
        print(f"Overall S/sqrt(S+B) Significance: {overall_s_sqrt_s_b:.2f}")

    return overall_s_b, overall_s_sqrt_s_b, total_sig, total_bg, total_sig_unscaled, total_bg_unscaled


###########################################################################################################
# classification_quality
###########################################################################################################

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

    Args:
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

    # Create confusion matrix
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
        
        # Plot classification rates
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
        axes[0].grid(True, alpha=0.8, linestyle='--')

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

        # Visualization of confusion matrix
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


###########################################################################################################
# choice
###########################################################################################################

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
        
        # CThreshold mapping for each particle type
        mask_thresh = {
            0: thresholds[0],  # Pion threshold
            1: thresholds[1],  # Kaon threshold  
            2: thresholds[2]   # Proton threshold
        }

        # Filter rows where any probability exceeds its threshold
        thresh_mask = ((df[cols_name[1]] > mask_thresh[0]) | 
                      (df[cols_name[2]] > mask_thresh[1]) | 
                      (df[cols_name[3]] > mask_thresh[2]))

        # Unclassified track are set to 211
        # df.iloc[~thresh_mask, 1] = 1
        # df.iloc[~thresh_mask, 2] = 0
        # df.iloc[~thresh_mask, 3] = 0
        df = df[thresh_mask].reset_index(drop=True)     # comment this and uncomment 3 lines above for 
                                                        # Unclassified track are set to 211
        
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


###########################################################################################################
# check_thresholds
###########################################################################################################

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
    
    ax.grid(visible=True, alpha=0.8, linestyle='--')
    ax.set_xlabel('Threshold Set')
    ax.set_ylabel('Classification Rate')
    ax.set_title('Classification Performance vs. Probability Thresholds')
    
    ax.set_xticks(np.arange(len(thresholds)))
    ax.set_xticklabels([str(thresh) for thresh in thresholds], rotation=45)
    
    ax.legend()
    plt.tight_layout()
    

###########################################################################################################
# get_distribution
###########################################################################################################
def get_distribution(
        distr: Union[pd.Series, np.ndarray], 
        unit_per_bin: float = None,
        interval: Tuple[float, float] = None,
        n_bins: int = None,
    ):
    
    """
    """
    data = np.asarray(distr)
    
    if n_bins is not None:
        bins=n_bins
    elif unit_per_bin is not None and interval is not None:
        bins = ceil(interval[1] - interval[0]) / unit_per_bin
    else:
        raise ValueError('n_bins or (unit_per_bin and interval) must be not None!')
    
    counts, bin_edges = np.histogram(data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    errors = np.sqrt(counts)
    
    mask = counts > 0
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    errors = errors[mask]
    
    return counts, bin_centers, bin_edges, errors


###########################################################################################################
# get_mass_distributions
###########################################################################################################
def get_mass_distributions(sig_mass_distr: Union[pd.Series, np.ndarray], 
                     bg_mass_distr: Union[pd.Series, np.ndarray], 
                     have_sig_events: int,
                     have_bg_events: int,
                     gev_per_bin: float = 0.000773,# (2.32497 - 2.24763) / 100
                     mass_interval: Tuple[float, float] = (2.24763, 2.32497)):
    """
    """
    
    # Calculate scaling weights to match real event expectations
    sig_weight = CFG.real_sig_events / have_sig_events
    bg_weight = CFG.real_bg_events / have_bg_events
    
    sig_mass_distr = np.asarray(sig_mass_distr)
    bg_mass_distr = np.asarray(bg_mass_distr)
    
    # Create weight arrays for signal and background
    weights_sig = np.full_like(sig_mass_distr, sig_weight, dtype=float)
    weights_bg = np.full_like(bg_mass_distr, bg_weight, dtype=float)
    
    mass_combined = np.concatenate([sig_mass_distr, bg_mass_distr])
    weights_combined = np.concatenate([weights_sig, weights_bg])
    labels_combined = np.concatenate([np.full(len(sig_mass_distr), 'Signal'), 
                                      np.full(len(bg_mass_distr), 'Background')])

    mass_df = pd.DataFrame({
        'mass': mass_combined,
        'weights': weights_combined,
        'type': labels_combined,
    })

    mass_mask = (mass_df['mass'] >= mass_interval[0]) & (mass_df['mass'] <= mass_interval[1])
    mass_df = mass_df[mass_mask]

    bins = ceil( (mass_interval[1] - mass_interval[0]) / gev_per_bin )
    binrange = (mass_df['mass'].min(), mass_df['mass'].max())

    bin_edges = np.linspace(binrange[0], binrange[1], bins + 1)
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

    total_heights = sig_heights + bg_heights

    # Errors
    sig_counts, _ = np.histogram(
        signal_data['mass'],
        bins=bin_edges
    )
    
    bg_counts, _ = np.histogram(
        background_data['mass'],
        bins=bin_edges
    )
    
    sig_heights_errors = np.sqrt(sig_counts) * sig_weight
    bg_heights_errors = np.sqrt(bg_counts) * bg_weight
    
    total_heights_errors = np.sqrt(sig_heights_errors**2 + bg_heights_errors**2)
    
    return total_heights, total_heights_errors, sig_heights, sig_heights_errors, bg_heights, bg_heights_errors, bin_edges, bin_centers



###########################################################################################################
# calculate_sb_ratio_calculate_sb_ratio_2gauss_const
###########################################################################################################
def calculate_sb_ratio_2gauss_const(
    fit_params,
    bin_edges,
    mass_interval=None,
    cov_matrix=None
):
    """
    Calculate S/B and its uncertainty for:
 
    f(x) = A1 * exp(-0.5*((x-mu1)/sigma1)**2)
         + A2 * exp(-0.5*((x-mu2)/sigma2)**2)
         + C
 
    Parameters
    ----------
    fit_params : array
        [A1, mu1, sigma1, A2, mu2, sigma2, C]
 
    bin_edges : array
        Histogram bin edges.
 
    mass_interval : tuple or None
        Integration interval.
 
    cov_matrix : 2D array
        Full covariance matrix from the fit.
    """
 
    version = 1.4
        
    if version != CFG.fit_version:
        raise 'Versions of fit and s/b calculations might be different!!'
 
    fit_params = np.asarray(fit_params, dtype=float)
 
    if mass_interval is None:
        x_min = bin_edges[0]
        x_max = bin_edges[-1]
    else:
        x_min, x_max = mass_interval
 
    def gaussian_integral(A, mu, sigma):
 
        z_min = (x_min - mu) / (np.sqrt(2) * sigma)
        z_max = (x_max - mu) / (np.sqrt(2) * sigma)
 
        return (
            A
            * sigma
            * np.sqrt(np.pi / 2)
            * (erf(z_max) - erf(z_min))
        )
 
    def calculate_R(params):
 
        A1, mu1, sigma1, A2, mu2, sigma2, C = params
 
        S1 = gaussian_integral(A1, mu1, sigma1)
        S2 = gaussian_integral(A2, mu2, sigma2)
 
        S = S1 + S2
 
        B = C * (x_max - x_min)
 
        R = S / B
 
        return S, B, R
 
    S, B, R = calculate_R(fit_params)
 
    result = {
        'signal': S,
        'background': B,
        'S_B_ratio': R,
        'integration_range': (x_min, x_max)
    }
 
    if cov_matrix is not None:
 
        cov_matrix = np.asarray(cov_matrix, dtype=float)
 
        gradient = np.zeros(len(fit_params))
 
        for i in range(len(fit_params)):
 
            params_plus = fit_params.copy()
            params_minus = fit_params.copy()
 
            step = 1e-5 * max(abs(fit_params[i]), 1.0)
 
            params_plus[i] += step
            params_minus[i] -= step
 
            _, _, R_plus = calculate_R(params_plus)
            _, _, R_minus = calculate_R(params_minus)
 
            gradient[i] = (
                R_plus - R_minus
            ) / (2 * step)
 
        variance_R = gradient @ cov_matrix @ gradient
        
        variance_R = max(variance_R, 0.0)
 
        result['S_B_uncertainty'] = np.sqrt(variance_R)
 
    return result


###########################################################################################################
# calculate_sb_ratio_calculate_sb_ratio_2gauss_pol1
###########################################################################################################
def calculate_sb_ratio_2gauss_pol1(
    fit_params,
    bin_edges,
    mass_interval=None,
    cov_matrix=None
):
    """
    Calculate S/B and its uncertainty for:
 
    f(x) = (A1 / (sqrt(2π) * sigma1)) * exp(-0.5*((x-mu1)/sigma1)^2)
         + (A2 / (sqrt(2π) * sigma2)) * exp(-0.5*((x-mu2)/sigma2)^2)
         + A3 + b * x
 
    Parameters
    ----------
    fit_params : array
        [A1, mu1, sigma1, A2, mu2, sigma2, A3, b]
 
    bin_edges : array
        Histogram bin edges.
 
    mass_interval : tuple or None
        Integration interval.
 
    cov_matrix : 2D array
        Full covariance matrix from the fit.
    """
 
    version = 1.6
    
    if version != CFG.fit_version:
        raise Exception('Versions of fit and s/b calculations might be different!!')
 
    fit_params = np.asarray(fit_params, dtype=float)
 
    if mass_interval is None:
        x_min = bin_edges[0]
        x_max = bin_edges[-1]
    else:
        x_min, x_max = mass_interval
 
 
    def gaussian_integral(A, mu, sigma):
        z_min = (x_min - mu) / (np.sqrt(2) * sigma)
        z_max = (x_max - mu) / (np.sqrt(2) * sigma)
        
        return A * (erf(z_max) - erf(z_min)) / 2


    def pol1_integral(A, b):
        func = lambda x: A * x + 0.5 * b * x**2
        return func(x_max) - func(x_min)
 
 
    def calculate_R(params):
        A1, mu1, sigma1, A2, mu2, sigma2, A3, b = params
 
        S1 = gaussian_integral(A1, mu1, sigma1)
        S2 = gaussian_integral(A2, mu2, sigma2)
 
        S = S1 + S2
        B = pol1_integral(A3, b)
        R = S / B
 
        return S, B, R
 
    S, B, R = calculate_R(fit_params)
 
    result = {
        'signal': S,
        'background': B,
        'S_B_ratio': R,
        'integration_range': (x_min, x_max)
    }
 
    if cov_matrix is not None:
        cov_matrix = np.asarray(cov_matrix, dtype=float)
 
        gradient_S = np.zeros(len(fit_params))
        gradient_B = np.zeros(len(fit_params))
        gradient_R = np.zeros(len(fit_params))
 
        for i in range(len(fit_params)):
            params_plus = fit_params.copy()
            params_minus = fit_params.copy()
 
            step = 1e-5 * max(abs(fit_params[i]), 1.0)
 
            params_plus[i] += step
            params_minus[i] -= step
 
            S_plus, B_plus, R_plus = calculate_R(params_plus)
            S_minus, B_minus, R_minus = calculate_R(params_minus)
 
            gradient_S[i] = (S_plus - S_minus) / (2 * step)
            gradient_B[i] = (B_plus - B_minus) / (2 * step)
            gradient_R[i] = (R_plus - R_minus) / (2 * step)
 
        variance_S = gradient_S @ cov_matrix @ gradient_S
        variance_B = gradient_B @ cov_matrix @ gradient_B
        variance_R = gradient_R @ cov_matrix @ gradient_R
 
        variance_S = max(variance_S, 0.0)
        variance_B = max(variance_B, 0.0)
        variance_R = max(variance_R, 0.0)
 
        result['signal_uncertainty'] = np.sqrt(variance_S)
        result['background_uncertainty'] = np.sqrt(variance_B)
        result['S_B_uncertainty'] = np.sqrt(variance_R)
 
    return result


###########################################################################################################
# calculate_tssa_errors_local_from_fit
###########################################################################################################
def calculate_tssa_errors_local(
    N_phi: int,
    N_phi_pi: int,
    phi_1: float,
    phi_2: float,
):
    
    """
    Calculates TSSA Uncertainties in local phi angle range.
    """
    
    def mean_abs_cos(phi_1, phi_2, n_points=10000):

        if phi_1 == phi_2:
            return abs(np.cos(phi_1))
    
        if phi_1 > phi_2:
            raise "Wrong boundaries"
    
        x = np.linspace(phi_1, phi_2, n_points)
        y = np.abs(np.cos(x))
    
        integral = np.trapz(y, x)
    
        return integral / (phi_2 - phi_1)
    
    
    A = 2 * N_phi * N_phi_pi / (N_phi + N_phi_pi)**2
    mean_abs_cos_val = mean_abs_cos(phi_1, phi_2)
    sigma_N_phi = np.sqrt(N_phi)
    sigma_N_phi_pi = np.sqrt(N_phi_pi)
    B = np.sqrt((sigma_N_phi / N_phi)**2 + (sigma_N_phi_pi / N_phi_pi)**2)
    
    sigma_tssa = A * B / mean_abs_cos_val
    
    return sigma_tssa


###########################################################################################################
# calculate_tssa_sig_errors_local
###########################################################################################################
def calculate_tssa_sig_errors_local(
    N_raw_phi: int,
    N_raw_phi_pi: int,
    N_bg_peak_phi: int,
    N_bg_peak_phi_pi: int,
    N_bg_side_phi: int,
    N_bg_side_phi_pi: int,
    phi_1: float,
    phi_2: float,
):
    
    r = (N_bg_peak_phi + N_bg_peak_phi_pi) / (N_raw_phi + N_raw_phi_pi)
    
    sigma_ann_raw = calculate_tssa_errors_local(
        N_phi=N_raw_phi,
        N_phi_pi=N_raw_phi_pi,
        phi_1=phi_1,
        phi_2=phi_2,
    )
    
    sigma_ann_bg_side = calculate_tssa_errors_local(
        N_phi=N_bg_side_phi,
        N_phi_pi=N_bg_side_phi_pi,
        phi_1=phi_1,
        phi_2=phi_2,
    )
    
    # print(f'{r=}')
    # print(f'{sigma_ann_raw=}')
    # print(f'{sigma_ann_bg_side=}')
    
    tssa_error = np.sqrt(sigma_ann_raw**2 + r**2 * sigma_ann_bg_side**2) / (1 - r)
    
    return tssa_error


###########################################################################################################
# calculate_tssa_sig_errors_local
###########################################################################################################
def epsilon_up(N_before, N_after, CL=0.95):
    """
    Односторонний верхний доверительный предел для эффективности отбора.
 
    Parameters
    ----------
    N_before : int
        Число событий до отбора.
    N_after : int
        Число событий после отбора.
    CL : float, optional
        Confidence level. По умолчанию 0.95.
 
    Returns
    -------
    float
        Upper limit for efficiency.
    """
 
    if N_before <= 0:
        raise ValueError("N_before must be > 0")
 
    if N_after < 0 or N_after > N_before:
        raise ValueError(
            "Must be 0 <= N_after <= N_before"
        )
 
    if N_after == 0:
        eps_up = 1.0 - (1.0 - CL)**(1.0 / N_before)
    else:
        eps_up = beta.ppf(
            CL,
            N_after + 1,
            N_before - N_after
        )
 
    return eps_up