import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from typing import Tuple, List


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

    axes[1].set_xlabel(f'{distr_name}', fontsize=8)
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
    axes[0].set_xlabel(f'{distr_name}', fontsize=8)
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