import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
from sklearn.metrics import auc
from typing import Tuple, List, Union
import os
import sys

# Path setup
try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.getcwd()
    
parent_dir = os.path.join(current_dir, '..', '..')
sys.path.insert(0, parent_dir)

from analysis_scripts.config import CFG
from analysis_scripts.estimate_scripts import get_mass_distributions

###########################################################################################################
# draw_distribution_pair
###########################################################################################################

def draw_distribution_pair(df: pd.DataFrame,
                       distr_name: str,
                       bins: int = 40) -> Tuple:
    """Draw combined boxplot and histogram for feature distribution analysis."""
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

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


###########################################################################################################
# draw_feature_distribution
###########################################################################################################
def draw_feature_distribution(
    df: pd.DataFrame,
    distr_name: str,
    axes: List, # [ax_bottom, ax_top]
    tag: str = None,
    hue: str = None,
    bins: int = 40,
    norma: bool = False,
    cut_point: float = None,
    select_direction: str = None
) -> Tuple:
    
    """Draw combined boxplot and histogram for feature distribution analysis."""

    df = df.copy()

    x_label = CFG.official_names[distr_name] if distr_name in CFG.official_names.keys() else distr_name

    if tag is not None:
        df = df[df.tag == tag]

    sns.boxplot(
        data=df,
        x=distr_name,
        y='tag',
        hue=hue,
        color=None if hue is not None else '#2E86AB',
        whis=1.5,
        ax=axes[1]
    )

    axes[1].set_xlabel(x_label, fontsize=8)
    axes[1].set_ylabel('Tag', fontsize=8)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].grid(True, alpha=0.8, linestyle='--')
    axes[1].tick_params(axis='both', labelsize=8)

    hist_ax = sns.histplot(
        data=df,
        x=distr_name,
        hue=hue,
        color=None if hue is not None else '#2E86AB',
        stat="density" if norma else 'count',
        common_norm=False,
        alpha=1,
        element='step',
        fill=False,
        linewidth=2,
        bins=bins,
        ax=axes[0]
    )

    # axes[0].set_title(f'{distr_name} Distribution', fontsize=14, fontweight='bold', pad=20)
    axes[0].set_xlabel(x_label, fontsize=8)
    y_label = 'Density' if norma else 'Counts'
    axes[0].set_ylabel(y_label, fontsize=8)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].grid(True, alpha=0.8, linestyle='--')
    axes[0].tick_params(axis='both', labelsize=8)

    # Threshold Line & arrow
    if cut_point is not None and select_direction is not None:

        y_min, y_max = axes[0].get_ylim()
        axes[0].axvline(cut_point, color='black', linestyle='dashed', 
                       linewidth=1.5, ymin=0, ymax=0.95)
    
        arrow_y_pos = y_max * 0.9    
        data = df[distr_name].dropna()
        counts, bin_edges = np.histogram(data, bins=bins)
        bin_width = bin_edges[1] - bin_edges[0]

        if select_direction == 'right':
            arrow_direction = cut_point + bin_width * 3
        else:
            arrow_direction = cut_point - bin_width * 3
            
        axes[0].annotate(
            '', 
            xy=(cut_point, arrow_y_pos), 
            xytext=(arrow_direction, arrow_y_pos),
            arrowprops=dict(arrowstyle='<-', color='black', lw=2)
        )
    
    return None


###########################################################################################################
# only_hist_plot
###########################################################################################################
def only_hist_plot(
    df: pd.DataFrame,
    distr_name: str,
    tag: str = None,
    hue: str = None,
    bins: int = 40,
    norma: bool = False,
    cut_point: float = None,
    select_direction: str = None,
    x_limits: List = None,
    fig = None,
    ax = None,
) -> Tuple:
    
    """Draw histogram for feature distribution analysis."""

    if (ax is None) or (fig is None):
        fig, ax = plt.subplots(figsize=(5, 5))

    df = df.copy()

    x_label = CFG.official_names[distr_name] if distr_name in CFG.official_names.keys() else distr_name

    if tag is not None:
        df = df[df.tag == tag]

    hist_ax = sns.histplot(
        data=df,
        x=distr_name,
        hue=hue,
        color=None if hue is not None else '#2E86AB',
        stat="density" if norma else 'count',
        common_norm=False,
        alpha=1,
        element='step',
        fill=False,
        linewidth=2,
        bins=bins,
        ax=ax
    )

    ax.set_title('Feature Distribution')
    ax.set_xlabel(x_label, fontsize=10)
    y_label = 'Density' if norma else 'Counts'
    ax.set_ylabel(y_label, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.8, linestyle='--')
    ax.tick_params(axis='both', labelsize=8)
    
    if x_limits is not None:
        ax.set_xlim(x_limits[0], x_limits[1])

    # Threshold Line & arrow
    if cut_point is not None and select_direction is not None:

        y_min, y_max = ax.get_ylim()
        ax.axvline(cut_point, color='black', linestyle='dashed', 
                       linewidth=1.5, ymin=0, ymax=0.95)
    
        arrow_y_pos = y_max * 0.9    
        data = df[distr_name].dropna()
        counts, bin_edges = np.histogram(data, bins=bins)
        bin_width = bin_edges[1] - bin_edges[0]

        if select_direction == 'right':
            arrow_direction = cut_point + bin_width * 3
        else:
            arrow_direction = cut_point - bin_width * 3
            
        ax.annotate(
            '', 
            xy=(cut_point, arrow_y_pos), 
            xytext=(arrow_direction, arrow_y_pos),
            arrowprops=dict(arrowstyle='<-', color='black', lw=2)
        )
    
    return fig, ax


###########################################################################################################
# draw_roc
###########################################################################################################

def draw_roc(tpr: List,
             fpr: List,
             select_direction: str,
             ax,
             best_cut_x: float = None,
             best_tpr: float = None,
             best_fpr: float = None,
             mask: np.ndarray = None,
             min_sig_sel: float = None,
    ) -> Tuple:
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
        None
    """

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

    roc_auc = auc(fpr, tpr)

    # Plot ROC curve
    sns.lineplot(
        data=roc_data,
        x='False Positive Rate',
        y='True Positive Rate',
        color='darkorange',
        linewidth=2,
        # label=f'ROC curve (AUC = {roc_auc:.3f})',
        ax=ax
    )
    
    # Plot dashed section
    if mask is not None and min_sig_sel is not None:
        
        # Dataframe for dashed section
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
            # label=f'Zone where Sig Eff < {min_sig_sel}',
            ax=ax
        )

    # Plot optimal threshold point
    sns.scatterplot(
        x=[best_fpr],
        y=[best_tpr],
        color='red',
        s=100,
        # label='Optimal threshold',
        zorder=5,
        ax=ax
    )

    # Add diagonal line
    ax.plot(
        [0, 1],
        [0, 1],
        color='navy',
        linestyle='--', 
        # label='Random classifier (AUC = 0.5)',
        alpha=0.8
    )

    # Define selection direction symbols
    sel_dir = {
        'right': '>',
        'left': '<'
    }

    if select_direction not in sel_dir:
        raise ValueError(f"select_direction must be 'right' or 'left', got '{select_direction}'")

    if (best_cut_x is not None) and (best_tpr is not None) and (best_fpr is not None):
        # Threshold annotation
        ax.text(best_fpr + 0.02, best_tpr - 0.05,
                f'Threshold: {sel_dir[select_direction]} {best_cut_x:.3f}\n',
                # f'TPR: {best_tpr:.3f}, FPR: {best_fpr:.3f}',
                fontsize=8, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        ax.text(
            0.7, 0.1,
            f'AUC: {roc_auc:.3f}', 
            fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
        )

    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=8)
    ax.set_ylabel('True Positive Rate', fontsize=8)
    # ax.set_title('Receiver Operating Characteristic (ROC) Curve', fontsize=14, pad=20)
    ax.set_title('ROC Curve', fontsize=8, pad=10)
    ax.grid(True, alpha=0.8, linestyle='--')
    # ax.legend(loc='lower right', framealpha=0.9)
    ax.tick_params(axis='both', which='major', labelsize=5)
    
    plt.tight_layout()


###########################################################################################################
# create_subplots_dynamic
###########################################################################################################

def create_subplots_dynamic(ncols, df, vars: list, bins: list = None, plot_type='hist', hue=False, norm=False):

    temp_df = df[['mc_pid'] + vars].copy()

    nrows = (temp_df.shape[1] + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5*ncols, 4*nrows))
    fig.suptitle('Initial track parameters distributions', fontsize=16, fontweight='bold')

    # Flatten axes for easy iteration
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


###########################################################################################################
# bar_plot
###########################################################################################################

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
    
    
###########################################################################################################
# hist_plot
###########################################################################################################
def hist_plot(
    data: pd.DataFrame,
    distr_name: str,
    ax,
    bins=40
) -> None:

    sns.histplot(
        data=data,
        x=distr_name,
        color='#2E86AB',
        alpha=1,
        element='step',
        fill=False,
        linewidth=2,
        bins=bins,
        ax=ax
    )

    ax.set_title(f'{distr_name} Distribution', fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel(f'{distr_name}', fontweight='bold')
    ax.set_ylabel('Counts', fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(True, alpha=0.8, linestyle='--')
    
    
###########################################################################################################
# total_spectrum
###########################################################################################################

def total_spectrum(sig_mass_distr: Union[pd.Series, np.ndarray], 
                     bg_mass_distr: Union[pd.Series, np.ndarray], 
                     have_sig_events: int,
                     have_bg_events: int,
                     mass_interval: Tuple[float, float] = (2.24763, 2.32497),
                     gev_per_bin: float = 0.000773,# (2.32497 - 2.24763) / 100
                     log_y: bool = True,
                     total: bool = True,
                     only_total: bool = False,
                     line_style: str = '--'):
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

    fig, ax = plt.subplots(figsize=(6, 6))

    colors = {'Signal': '#1f77b4', 'Background': 'black', 'total': '#2ca02c'}
    
    if not total or (total and not only_total):
        
        ax.step(
            bin_edges[:-1],
            sig_heights,
            where='post',
            color=colors['Signal'],
            linewidth=2,
            linestyle=line_style,
            label='Signal'
        )
        
        ax.step(
            bin_edges[:-1],
            bg_heights,
            where='post',
            color=colors['Background'],
            linewidth=2,
            linestyle=line_style,
            label='Background'
        )

    if total:
        ax.step(
            bin_edges[:-1],
            total_heights,
            where='post',
            color=colors['total'],
            linewidth=2.5,
            linestyle=line_style,
            label='Total (S+B)'
        )

    if log_y:
        ax.set_yscale('log')
        # ax.set_ylim(bottom=0.1)

    ax.set_title('Yield Mass Distribution', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Mass (GeV)', fontweight='bold', fontsize=12)
    y_label = 'Counts (log scale)' if log_y else 'Counts'
    ax.set_ylabel(y_label, fontweight='bold', fontsize=12)

    ax.legend(
        loc='best',
        title='Mass Spectrum',
        frameon=True,
        fancybox=True,
        shadow=True
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.8, linestyle='--')
    
    plt.tight_layout()
    return fig, ax
    

###########################################################################################################
# tssa_area_plot
###########################################################################################################

def tssa_area_plot(
    sig_mass_distr: Union[pd.Series, np.ndarray], 
    bg_mass_distr: Union[pd.Series, np.ndarray], 
    have_sig_events: int,
    have_bg_events: int,
    fit_params: Union[List, np.ndarray],
    sig_eff: float,
    bg_eff: float,
    borders_pos: Union[List, np.ndarray],
    mass_interval: Tuple[float, float] = (2.24763, 2.32497),
    gev_per_bin: float = 0.000773,# (2.32497 - 2.24763) / 100
    log_y: bool = True,
    total: bool = True,
    only_total: bool = False,
    line_style: str = '--'):
    
    """
    params: 8 parameters [A1, mu1, sigma1, A2, mu2, sigma2, A3, b]
    borders_pos: 6 params [left_left_bg, left_right_bg, left_sif, right_sig, right_left_bg, right_right_bg]
    """
    
    fig, ax = total_spectrum(
        sig_mass_distr=sig_mass_distr, 
        bg_mass_distr=bg_mass_distr, 
        have_sig_events=have_sig_events / sig_eff,
        have_bg_events=have_bg_events / bg_eff,
        mass_interval=mass_interval,
        gev_per_bin=gev_per_bin,
        log_y=log_y,
        total=total,
        only_total=only_total,
        line_style=line_style
    )
    
    total_heights, total_heights_errors, sig_heights, sig_heights_errors, bg_heights, bg_heights_errors, bin_edges, bin_centers = get_mass_distributions(
        sig_mass_distr=sig_mass_distr, 
        bg_mass_distr=bg_mass_distr, 
        have_sig_events=have_sig_events / sig_eff, 
        have_bg_events=have_bg_events / bg_eff,
        mass_interval=mass_interval,
        gev_per_bin=0.000773
    )
            
    x_bins = np.linspace(bin_edges[0], bin_edges[-1], 200)
            
    g1 = fit_params[0] / (np.sqrt(2*np.pi) * fit_params[2]) * np.exp( -0.5 * ((x_bins - fit_params[1]) / fit_params[2])**2 )
    g2 = fit_params[3] / (np.sqrt(2*np.pi) * fit_params[5]) * np.exp( -0.5 * ((x_bins - fit_params[4]) / fit_params[5])**2 )            
    pol1 = fit_params[6] + fit_params[7] * x_bins

    g1_label = 'Gaussian 1'
    g2_label = 'Gaussian 2'
    pol2_label = 'Pol1'
    total_label = "Fit"
    
    ax.plot(x_bins, g1, 'g:', label=g1_label, alpha=0.7, linewidth=2.5)
    ax.plot(x_bins, g2, 'b:', label=g2_label, alpha=0.7, linewidth=2.5)
    ax.plot(x_bins, pol1, 'm:', label=pol2_label, alpha=0.7, linewidth=2.5)
    ax.plot(x_bins, g1 + g2 + pol1, 'r:', label=total_label, linestyle='-', linewidth=3, alpha=0.7)
    
    ax.axvline(x=borders_pos[0], ymin=0., ymax=0.95, color='brown', linestyle='--', linewidth=1.5)
    ax.axvline(x=borders_pos[1], ymin=0., ymax=0.95, color='brown', linestyle='--', linewidth=1.5)
    
    ax.axvline(x=borders_pos[2], ymin=0., ymax=0.95, color='blue', linestyle='--', linewidth=1.5)
    ax.axvline(x=borders_pos[3], ymin=0., ymax=0.95, color='blue', linestyle='--', linewidth=1.5)
    
    ax.axvline(x=borders_pos[4], ymin=0., ymax=0.95, color='brown', linestyle='--', linewidth=1.5)
    ax.axvline(x=borders_pos[5], ymin=0., ymax=0.95, color='brown', linestyle='--', linewidth=1.5) 
    
    ax.legend()
    
    plt.tight_layout()
    return fig, ax
    
    
###########################################################################################################
# tssa_error_plot
###########################################################################################################

def tssa_error_plot(
    gaps: Union[np.ndarray, List],
    tssa_errors_list: Union[np.ndarray, List],
    x_label: str,
):
    
    fig, ax = plt.subplots(figsize=(12, 6))

    max_error = np.max(tssa_errors_list)

    for i, ((start, end), error) in enumerate(zip(gaps, tssa_errors_list)):
    
        width = end - start
        color = 'steelblue'
        edgecolor = 'black'
        linewidth = 1.5
        alpha = 0.7
        
        rect = plt.Rectangle(
            (start, -error),
            width,
            2 * error,
            fill=True,
            alpha=alpha,
            edgecolor=edgecolor,
            linewidth=linewidth,
            facecolor=color
        )
        ax.add_patch(rect)
        
        mid_x = start + width/2
        ax.text(mid_x, 0.1 * max_error, f'+-{error:.6f}', ha='center', va='center', fontsize=10, color='black')
        
        ax.text(mid_x, -max(tssa_errors_list)*0.1, f'[{start:.1f}, {end:.1f}]', ha='center', va='top', fontsize=10, color='black')
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.5, alpha=0.7)

    ax.set_xlim(gaps[0][0] - 0.02, gaps[-1][1] + 0.02)
    ylim_max = max(tssa_errors_list) * 1.4
    ax.set_ylim(-ylim_max, ylim_max)

    ax.set_xlabel(f'{x_label} Interval', fontsize=12)
    ax.set_ylabel('TSSA Error', fontsize=12)
    ax.set_title(f'TSSA Error by {x_label} Interval', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.8, linestyle='--')

    plt.tight_layout()
    
    return fig, ax


###########################################################################################################
# selection_eff_gap_plot
###########################################################################################################

def selection_eff_gap_plot(
    gaps: Union[np.ndarray, List],
    sig_eff_list: Union[np.ndarray, List],
    bg_eff_list: Union[np.ndarray, List],
    x_label: str,
):
    
    fig, ax = plt.subplots(figsize=(12, 6))

    max_eff = np.max(np.concatenate((sig_eff_list, bg_eff_list), axis=0))
    min_eff = np.min(np.concatenate((sig_eff_list, bg_eff_list), axis=0))

    for i, ((start, end), sig_eff, bg_eff) in enumerate(zip(gaps, sig_eff_list, bg_eff_list)):
    
        width = end - start
        
        rect_sig = plt.Rectangle(
            (start, 0),
            width,
            sig_eff,
            fill=True,
            alpha=0.7,
            edgecolor='black',
            linewidth=1.5,
            facecolor='steelblue',
            label='Signal Efficiency' if i == 0 else ""
        )
        ax.add_patch(rect_sig)
        
        rect_bg = plt.Rectangle(
            (start, 0),
            width,
            bg_eff,
            fill=True,
            alpha=0.7,
            edgecolor='black',
            linewidth=1.5,
            facecolor='coral',
            label='Background Efficiency' if i == 0 else ""
        )
        ax.add_patch(rect_bg)
        
        mid_x = start + width / 2

        ax.text(mid_x, sig_eff / 2, f'{sig_eff:.2e}', 
                ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        ax.text(mid_x, bg_eff / 2, f'{bg_eff:.2e}', 
                ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        
        ax.text(mid_x, 1e-3 * min_eff, f'[{start:.1f}, {end:.1f}]', 
                ha='center', va='top', fontsize=10, color='black')
    

    ax.set_xlim(gaps[0][0] - 0.02, gaps[-1][1] + 0.02)
    ylim_max = max_eff * 1e2
    ax.set_ylim(1e-8, ylim_max)

    ax.set_xlabel(f'{x_label} Interval', fontsize=12)
    ax.set_ylabel('Efficiency', fontsize=12)
    ax.set_title(f'Selection Efficiency by {x_label} Interval', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.8, linestyle='--')
    
    ax.set_yscale('log')
    
    ax.legend(loc='upper right', framealpha=0.9)

    plt.tight_layout()
    
    return fig, ax