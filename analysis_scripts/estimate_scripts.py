import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, Union
import sys


sys.path.append('../analysis_scripts')

from fit_scripts import fit_distr_flat_bg


###########################################################################################################
# signal_estimates
###########################################################################################################

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
    
    ax.grid(visible=True, alpha=0.3, linestyle='--')
    ax.set_xlabel('Threshold Set')
    ax.set_ylabel('Classification Rate')
    ax.set_title('Classification Performance vs. Probability Thresholds')
    
    ax.set_xticks(np.arange(len(thresholds)))
    ax.set_xticklabels([str(thresh) for thresh in thresholds], rotation=45)
    
    ax.legend()
    plt.tight_layout()