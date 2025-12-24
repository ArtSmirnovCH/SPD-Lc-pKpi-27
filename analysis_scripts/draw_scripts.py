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


###########################################################################################################
# draw_feature_distribution
###########################################################################################################

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
        bins=bins,
        ax=axes[1]
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
        # Get the current y-limits in data coordinates
        y_min, y_max = axes[1].get_ylim()
        
        # Set the vline to go from bottom to near the top (95% of the way up in axes coordinates)
        # ymax=0.95 means 95% from bottom to top in axes coordinates
        axes[1].axvline(cut_point, color='black', linestyle='dashed', 
                       linewidth=2, ymin=0, ymax=0.95)
        
        # Calculate position for arrow (in data coordinates)
        arrow_y_pos = y_max * 0.9  # 90% of the way up in data coordinates
        
        data = df[distr_name].dropna()
        counts, bin_edges = np.histogram(data, bins=bins)
        bin_width = bin_edges[1] - bin_edges[0]

        if select_direction == 'right':
            arrow_direction = cut_point + bin_width * 3
        else:
            arrow_direction = cut_point - bin_width * 3
            
        # Draw arrow in data coordinates
        axes[1].annotate('', 
                        xy=(cut_point, arrow_y_pos), 
                        xytext=(arrow_direction, arrow_y_pos),
                        arrowprops=dict(arrowstyle='<-', color='black', lw=2))

    plt.tight_layout()
    
    return fig, axes


###########################################################################################################
# draw_roc
###########################################################################################################

def draw_roc(tpr: List,
             fpr: List,
             select_direction: str,
             best_cut_x: float = None,
             best_tpr: float = None,
             best_fpr: float = None,
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

    if (best_cut_x is not None) and (best_tpr is not None) and (best_fpr is not None):
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


###########################################################################################################
# create_subplots_dynamic
###########################################################################################################

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