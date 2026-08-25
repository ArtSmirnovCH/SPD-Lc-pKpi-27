import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from typing import Tuple, Union
from matplotlib.ticker import FuncFormatter

import sys
import os

# Path setup
try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.getcwd()
    
parent_dir = os.path.join(current_dir, '..', '..')
sys.path.insert(0, parent_dir)

from analysis_scripts.config import CFG


###########################################################################################################
# fit_distr_triple_gauss
###########################################################################################################

def fit_distr_triple_gauss(distr: Union[pd.Series, np.ndarray],
              distr_name: str,
              title: str,
              x_label: str,
              optimize_method: str = 'L-BFGS-B',
              plot: bool = False,
              bins: int = 40,
              max_iter: int=1000) -> Tuple:

    version = 1.4
        
    if version != CFG.fit_version:
        raise 'Versions of fit and s/b calculations might be different!!'

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

    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])
    
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min

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
        sigma^2_eff = w1 * sigma1**2 + w2 * sigma2**2 + w3 * sigma3**2
        """
        A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = fit_params
        
        # Weights
        total_amplitude = A1 + A2 + A3
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
        w3 = A3 / total_amplitude
        
        effective_mean = w1 * mu1 + w2 * mu2 + w3 * mu3
        
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
            
            # Data and fit
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
            
            # Fitted curve
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
            
            # Add fit parameters
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
            
            # Pulls
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
            
            # Pull statistics
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



###########################################################################################################
# fit_distr_flat_bg
###########################################################################################################

def fit_distr_flat_bg(distr: Union[pd.Series, np.ndarray],
              title: str,
              x_label: str,
              optimize_method: str = 'L-BFGS-B',
              plot: bool = False,
              bins: int = 40,
              max_iter: int = 1000) -> Tuple:

    version = 1.4
        
    if version != CFG.fit_version:
        raise 'Versions of fit and s/b calculations might be different!!'

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
    
    def constant(x, *params):
        """
        Composite model of 2 Gaussian functions with constant background
        
        Args:
            x: array-like, input values
            params: 7 parameters [A1, mu1, sigma1, A2, mu2, sigma2, C]
        
        Returns:
            y: sum of 2 Gaussian functions and Constant
        """
        C = params

        return np.ones_like(x) * C
    
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])
    
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    initial_guess = [initial_amplitude]   # Constant
    
    bounds = [(0, None)]   # Constant
    
    
    #=======================================================================================
    # Loss Function and Optimization
    #=======================================================================================
    
    def weighted_loss(params, x, y, errors):
        """Objective function to minimize (sum of weighted squared residuals)"""
        predictions = constant(x, *params)
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

    n_params = len(fit_params)
    n_points = len(errors)
    n_dof = n_points - n_params  # degrees of freedom

    #=======================================================================================
    # Calculate Parameter Errors
    #=======================================================================================
    
    def calculate_parameter_errors(result, x, y, errors):
        """
        Calculate errors for fitted parameters using the covariance matrix
        
        The covariance matrix is approximated by the inverse of the Hessian matrix
        at the minimum, scaled by the reduced chi-squared
        """
        
        # Calculate Hessian matrix using finite differences
        def hessian_finite_diff(fun, x0, args=(), epsilon=1e-4):
            """Calculate Hessian matrix using central finite differences"""
            n = len(x0)
            hess = np.zeros((n, n))
            
            for i in range(n):
                for j in range(i, n):
                    # Create basis vectors
                    ei = np.zeros(n)
                    ej = np.zeros(n)
                    ei[i] = 1.0
                    ej[j] = 1.0
                    
                    # Central difference for Hessian
                    if i == j:
                        # Diagonal elements
                        f_plus = fun(x0 + epsilon * ei, *args)
                        f_minus = fun(x0 - epsilon * ei, *args)
                        f0 = fun(x0, *args)
                        hess[i, i] = (f_plus - 2*f0 + f_minus) / (epsilon**2)
                    else:
                        # Off-diagonal elements
                        f_plus_plus = fun(x0 + epsilon * ei + epsilon * ej, *args)
                        f_plus_minus = fun(x0 + epsilon * ei - epsilon * ej, *args)
                        f_minus_plus = fun(x0 - epsilon * ei + epsilon * ej, *args)
                        f_minus_minus = fun(x0 - epsilon * ei - epsilon * ej, *args)
                        hess[i, j] = (f_plus_plus - f_plus_minus - f_minus_plus + f_minus_minus) / (4 * epsilon**2)
                        hess[j, i] = hess[i, j]
            
            return hess
        
        # Calculate Hessian at minimum
        hess = hessian_finite_diff(
            lambda p: weighted_loss(p, x, y, errors), 
            result.x, 
            epsilon=1e-4
        )
        
        # Calculate covariance matrix (inverse of Hessian)
        try:
            cov_matrix = np.linalg.inv(hess)
            
            # Calculate chi-squared and reduced chi-squared
            predictions = constant(x, *result.x)
            residuals = y - predictions
            weights = 1.0 / (errors**2 + 1e-6)
            chi2 = np.sum(weights * residuals**2)
            reduced_chi2 = chi2 / n_dof if n_dof > 0 else chi2
            
            # Scale covariance matrix by reduced chi-squared
            cov_matrix *= reduced_chi2
            
            # Parameter errors are sqrt of diagonal elements
            param_errors = np.sqrt(np.diag(cov_matrix))
            
            return param_errors, cov_matrix, chi2, reduced_chi2
            
        except np.linalg.LinAlgError:
            # If matrix inversion fails, return NaNs
            print("Warning: Could not calculate covariance matrix (singular Hessian)")
            n_params = len(result.x)
            return np.full(n_params, np.nan), np.full((n_params, n_params), np.nan), np.nan, np.nan
    
    # Calculate parameter errors
    param_errors, cov_matrix, chi2, reduced_chi2 = calculate_parameter_errors(
        result, bin_centers, counts, errors
    )
    

    print("=" * 60)
    print("Optimization result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of iterations: {result.nit}")
    print(f"Final objective value: {result.fun:.6f}")
    print(f"Chi-squared: {chi2:.2f}")
    print(f"Reduced chi-squared: {reduced_chi2:.3f}")
    print("-" * 60)
    print("Fit parameters with errors:")
    print(f"Constant: C = {fit_params[0]:.4f} +- {param_errors[0]:.4f}")
    print("=" * 60)

    # Draw Data and Fit Curve
    if plot:
        
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [2, 1]})
        
        if result.success:
            
            # Calculate fit values and pulls
            fit_values = constant(bin_centers, *fit_params)
            residuals = counts - fit_values
            pulls = residuals / errors  # (data - fit) / error
            
            # Enable LaTeX
            plt.rcParams.update({
                "text.usetex": True,
                "font.family": "serif",
                "font.size": 12
            })
            
            # Data and fit
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
            y_fit = constant(x_fit, *fit_params)
            
            axes[0].plot(
                x_fit,
                y_fit,
                color='red',
                linestyle='-',
                linewidth=2, 
                label='Double Gaussian Fit'
            )
            
            C = fit_params[0]
            
            # Format parameter strings with errors
            if not np.isnan(param_errors[0]):
                const_label = fr'Const: C={fit_params[0]:.2f} $\pm$ {param_errors[0]:.2f}'
            else:
                const_label = fr'Const: C={fit_params[0]:.2f}'
            
            axes[0].plot(x_fit, constant(x_fit, *fit_params), 'm:', label=const_label, alpha=0.7)
            
            # Add fit parameters as text with errors
            if not np.isnan(param_errors[0]):
                fit_text = (fr'Fit Parameters:'
                          fr'\\Const: C={fit_params[0]:.1f} $\pm$ {param_errors[0]:.1f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            else:
                fit_text = (fr'Fit Parameters:'
                          fr'\\Const: C={fit_params[0]:.1f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
    
            
            axes[0].text(0.02, 0.98, fit_text, transform=axes[0].transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)
            
            axes[0].set_title(f'{title}', fontsize=14, fontweight='bold', pad=20)
            axes[0].set_xlabel(x_label, fontweight='bold')
            axes[0].set_ylabel('Counts', fontweight='bold')
            axes[0].legend(loc='best')
            axes[0].grid(True, alpha=0.3, linestyle='--')
            
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            
            # Pulls
            axes[1].axhline(y=0, color='red', linestyle='-', linewidth=1, alpha=0.7)
            axes[1].axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=-1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            axes[1].axhline(y=-2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            
            bin_width = bin_edges[1] - bin_edges[0]
            axes[1].bar(bin_centers, pulls, 
                        width=bin_width*0.8,  # 80% of bin width for spacing
                        color='skyblue', 
                        edgecolor='darkblue',
                        alpha=0.7, 
                        label='Pulls')
            
            # Calculate pull statistics
            mean_pull = np.mean(pulls)
            std_pull = np.std(pulls)
            
            pull_stats = (fr'Pull Statistics:'
                          fr'\\$\mu_{{pull}}$ = {mean_pull:.3f}'
                          fr'\\$\sigma_{{pull}}$ = {std_pull:.3f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            
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
        
        # Return additional information
        fit_results = {
            'params': fit_params,
            'errors': param_errors,
            'cov_matrix': cov_matrix,
            'chi2': chi2,
            'reduced_chi2': reduced_chi2,
            'result': result
        }
        
        return fit_results, fig, axes
    
    fit_results = {
        'params': fit_params,
        'errors': param_errors,
        'cov_matrix': cov_matrix,
        'chi2': chi2,
        'reduced_chi2': reduced_chi2,
        'result': result
    }
    
    return fit_results


###########################################################################################################
# fit_distr_double_gauss_and_pol2
###########################################################################################################

def fit_distr_double_gauss_and_pol2(
    bin_centers,
    bin_edges,
    counts,
    errors,
    title: str,
    x_label: str,
    optimize_method: str = 'L-BFGS-B',
    plot: bool = False,
    log_y: bool = False,
    max_iter: int = 1000,
    scaler: float = 1e4,
    scaler_round: int = 0,
    unit_name: str = '$\mu m$'
) -> Tuple:
    
    version = 1.5
        
    if version != CFG.fit_version:
        raise 'Versions of fit and s/b calculations might be different!!'
    
    mask = counts > 0
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    errors = errors[mask]

    #=======================================================================================
    # Fit Model - Double Gaussian and Pol2
    #=======================================================================================
    
    def double_gaussian_and_pol_2(x, params):
        """
        Composite model of 2 Gaussian functions with pol2 background
        
        Args:
            x: array-like, input values
            params: 7 parameters [A1, mu1, sigma1, A2, mu2, sigma2, A3, b, c]
        
        Returns:
            y: sum of 2 Gaussian functions and pol2
        """
        A1, mu1, sigma1, A2, mu2, sigma2, A3, b, c = params
        
        # g1 = A1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
        # g2 = A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
        # pol2 = A3 + b * x + c * x**2
        
        g1 = A1 / (np.sqrt(2*np.pi) * sigma1) * np.exp( -0.5 * ((x - mu1) / sigma1)**2 )
        g2 = A2 / (np.sqrt(2*np.pi) * sigma2) * np.exp( -0.5 * ((x - mu2) / sigma2)**2 )
        # z = x - np.mean(bin_centers)
        pol2 = A3 + b * x + c * x**2
        
        return g1 + g2 + pol2
    
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])
    
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    initial_guess = [
        initial_amplitude, data_min + 0.25 * data_range, 0.1 * data_range,  # Gaussian 1
        initial_amplitude, data_min + 0.75 * data_range, 0.1 * data_range,  # Gaussian 2  
        initial_amplitude, 1, -2,   # Pol2
    ]
    
    # Set bounds for parameters
    bounds = [
        (0, None), (2.284, 2.287), (1e-6, 0.5 * data_range),  # Gaussian 1
        (0, None), (2.284, 2.287), (1e-6, 0.5 * data_range),  # Gaussian 2
        (0, None), (None, None), (None, 0)   # Pol2
    ]
    
    def chi2_function(params, x, y, errors):
 
        predictions = double_gaussian_and_pol_2(x, params)
 
        residuals = y - predictions
 
        return np.sum((residuals / errors) ** 2)
        
        
    # =========================================================================
    # Fit
    # =========================================================================
 
    result = minimize(
        fun=chi2_function,
        x0=initial_guess,
        args=(bin_centers, counts, errors),
        method=optimize_method,
        bounds=bounds,
        options={"maxiter": max_iter}
    )
 
    if not result.success:
        print("Warning: fit did not fully converge:")
        print(result.message)
 
    fit_params = result.x
 
    n_points = len(counts)
    n_params = len(fit_params)
    n_dof = n_points - n_params
        
    # =========================================================================
    # Numerical Hessian
    # =========================================================================
 
    def numerical_hessian(fun, x0, epsilon=1e-4):
        """
        Numerical Hessian using central finite differences.
        """
 
        x0 = np.asarray(x0, dtype=float)
 
        n = len(x0)
        hess = np.zeros((n, n))
 
        steps = epsilon * np.maximum(np.abs(x0), 1.0)
 
        f0 = fun(x0)
 
        for i in range(n):
 
            hi = steps[i]
 
            ei = np.zeros(n)
            ei[i] = hi
 
            # Diagonal second derivative
            f_plus = fun(x0 + ei)
            f_minus = fun(x0 - ei)
 
            hess[i, i] = (
                f_plus - 2.0 * f0 + f_minus
            ) / hi**2
 
            for j in range(i + 1, n):
 
                hj = steps[j]
 
                ej = np.zeros(n)
                ej[j] = hj
 
                f_pp = fun(x0 + ei + ej)
                f_pm = fun(x0 + ei - ej)
                f_mp = fun(x0 - ei + ej)
                f_mm = fun(x0 - ei - ej)
 
                value = (
                    f_pp
                    - f_pm
                    - f_mp
                    + f_mm
                ) / (4.0 * hi * hj)
 
                hess[i, j] = value
                hess[j, i] = value
 
        return hess
        
    # =========================================================================
    # Covariance matrix and parameter errors
    # =========================================================================
 
    def calculate_parameter_errors(
        result,
        x,
        y,
        errors,
        n_dof
    ):
        """
        Calculate parameter covariance matrix.
 
        Since the Hessian is calculated for chi2:
 
            H = d²(chi2) / dtheta_i dtheta_j
 
        the covariance matrix in the quadratic approximation is:
 
            Cov = 2 * H^{-1}
 
        If the uncertainties are only known up to a common scale,
        the covariance is additionally scaled by reduced chi2.
        """
 
        objective = lambda p: chi2_function(
            p, x, y, errors
        )
 
        hess = numerical_hessian(
            objective,
            result.x,
            epsilon=1e-4
        )
 
        # Symmetrize to reduce numerical noise
        hess = 0.5 * (hess + hess.T)
 
        chi2 = objective(result.x)
 
        reduced_chi2 = (
            chi2 / n_dof
            if n_dof > 0
            else np.nan
        )
 
        try:
 
            # Covariance for Hessian of chi2
            cov_matrix = 2.0 * np.linalg.inv(hess)
 
            # Optional scale correction
            #
            # Keep this if errors are estimates and the overall
            # noise scale may be incorrect.
            #
            # If errors are known exact Poisson errors, one may
            # choose NOT to multiply by reduced_chi2.
            
            # Maybe remove
            if n_dof > 0:
                cov_matrix *= reduced_chi2
 
            # Symmetrize covariance
            cov_matrix = 0.5 * (
                cov_matrix + cov_matrix.T
            )
 
            variances = np.diag(cov_matrix)
 
            # Negative values can appear from numerical noise
            param_errors = np.sqrt(
                np.maximum(variances, 0.0)
            )
 
            return (
                param_errors,
                cov_matrix,
                chi2,
                reduced_chi2,
                hess
            )
 
        except np.linalg.LinAlgError:
 
            print(
                "Warning: Hessian is singular. "
                "Covariance matrix cannot be calculated."
            )
 
            n = len(result.x)
 
            return (
                np.full(n, np.nan),
                np.full((n, n), np.nan),
                chi2,
                reduced_chi2,
                hess
            )
 
    param_errors, cov_matrix, chi2, reduced_chi2, hess = calculate_parameter_errors(
        result=result,
        x=bin_centers,
        y=counts,
        errors=errors,
        n_dof=n_dof
    )
    
    # =========================================================================
    # Effective parameters of the Gaussian mixture
    # =========================================================================
 
    def calculate_effective_parameters(
        fit_params,
        cov_matrix
    ):
        """
        Calculate mean and standard deviation of the full
        two-Gaussian mixture.
 
        Parameters:
            [A1, mu1, sigma1, A2, mu2, sigma2, C]
 
        For the mixture weights:
 
            w1 = A1 / (A1 + A2)
            w2 = A2 / (A1 + A2)
 
        Mean:
 
            mu_eff = w1*mu1 + w2*mu2
 
        Full variance of the mixture:
 
            Var_eff =
                w1*sigma1²
                + w2*sigma2²
                + w1*w2*(mu1 - mu2)²
 
        Equivalently:
 
            Var_eff =
                w1*(sigma1² + (mu1 - mu_eff)²)
                + w2*(sigma2² + (mu2 - mu_eff)²)
 
        Errors are propagated using the full covariance matrix:
 
            Var(f) = grad(f)^T Cov grad(f)
        """
 
        A1, mu1, sigma1, A2, mu2, sigma2, A3, b, c = fit_params
 
        total_amplitude = A1 + A2
 
        if total_amplitude <= 0:
            return np.nan, np.nan, np.nan, np.nan
 
        # ---------------------------------------------------------------------
        # Weights
        # ---------------------------------------------------------------------
 
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
 
        # ---------------------------------------------------------------------
        # Effective mean
        # ---------------------------------------------------------------------
 
        effective_mean = (
            w1 * mu1
            + w2 * mu2
        )
 
        # ---------------------------------------------------------------------
        # Full effective variance
        # ---------------------------------------------------------------------
 
        delta_mu = mu1 - mu2
 
        effective_variance = (
            w1 * sigma1**2
            + w2 * sigma2**2
            + w1 * w2 * delta_mu**2
        )
 
        effective_variance = max(
            effective_variance,
            0.0
        )
 
        effective_sigma = np.sqrt(
            effective_variance
        )
 
        # ---------------------------------------------------------------------
        # Error propagation
        # ---------------------------------------------------------------------
 
        if (
            cov_matrix is None
            or np.any(~np.isfinite(cov_matrix))
        ):
            return (
                effective_sigma,
                effective_mean,
                np.nan,
                np.nan
            )
 
        # =====================================================================
        # Gradient of effective mean
        # =====================================================================
 
        grad_mean = np.zeros(len(fit_params))
 
        grad_mean[0] = (
            A2 * (mu1 - mu2)
            / total_amplitude**2
        )
 
        grad_mean[1] = w1
 
        grad_mean[2] = 0.0
 
        grad_mean[3] = (
            A1 * (mu2 - mu1)
            / total_amplitude**2
        )
 
        grad_mean[4] = w2
 
        grad_mean[5] = 0.0
 
        grad_mean[6] = 0.0
 
        # Full covariance propagation
        mean_variance = (
            grad_mean
            @ cov_matrix
            @ grad_mean
        )
 
        effective_mean_err = np.sqrt(
            max(mean_variance, 0.0)
        )
 
        # =====================================================================
        # Gradient of full effective variance
        #
        # V =
        #   w1*sigma1²
        #   + w2*sigma2²
        #   + w1*w2*(mu1-mu2)²
        #
        # =====================================================================
 
        grad_var = np.zeros(len(fit_params))
 
        dV_dw1 = (
            sigma1**2
            + w2 * delta_mu**2
        )
 
        dV_dw2 = (
            sigma2**2
            + w1 * delta_mu**2
        )
 
        grad_var[0] = (
            A2 / total_amplitude**2
        ) * (dV_dw1 - dV_dw2)
 
        grad_var[3] = (
            A1 / total_amplitude**2
        ) * (dV_dw2 - dV_dw1)
 
        # Derivatives with respect to mu1 and mu2
 
        grad_var[1] = (
            2.0
            * w1
            * w2
            * delta_mu
        )
 
        grad_var[4] = (
            -2.0
            * w1
            * w2
            * delta_mu
        )
 
        # Derivatives with respect to sigma1 and sigma2
 
        grad_var[2] = (
            2.0 * w1 * sigma1
        )
 
        grad_var[5] = (
            2.0 * w2 * sigma2
        )
 
        # Constant background does not enter effective moments
        grad_var[6] = 0.0
 
        # Error of variance
        variance_variance = (
            grad_var
            @ cov_matrix
            @ grad_var
        )
 
        variance_err = np.sqrt(
            max(variance_variance, 0.0)
        )
 
        if effective_sigma > 0:
 
            effective_sigma_err = (
                variance_err
                / (2.0 * effective_sigma)
            )
 
        else:
 
            effective_sigma_err = np.nan
 
        return (
            effective_sigma,
            effective_mean,
            effective_sigma_err,
            effective_mean_err
        )
        
    effective_sigma, effective_mean, effective_sigma_err, effective_mean_err = calculate_effective_parameters(
        fit_params,
        cov_matrix
    )

    print("=" * 60)
    print("Optimization result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of iterations: {result.nit}")
    print(f"Final objective value: {result.fun:.6f}")
    print(f"Chi-squared: {chi2:.2f}")
    print(f"Reduced chi-squared: {reduced_chi2:.3f}")
    print("-" * 60)
    print("Fit parameters with errors:")
    print(f"Gauss 1: A = {fit_params[0]:.4f} +- {param_errors[0]:.4f}")
    print(f"         mu = {fit_params[1]:.6f} +- {param_errors[1]:.6f}")
    print(f"         sigma = {fit_params[2]:.6f} +- {param_errors[2]:.6f}")
    print(f"Gauss 2: A = {fit_params[3]:.4f} +- {param_errors[3]:.4f}")
    print(f"         mu = {fit_params[4]:.6f} +- {param_errors[4]:.6f}")
    print(f"         sigma = {fit_params[5]:.6f} +- {param_errors[5]:.6f}")
    print(f"Pol2: A = {fit_params[6]:.4f} +- {param_errors[6]:.4f}")
    print(f"      b = {fit_params[7]:.6f} +- {param_errors[7]:.6f}")
    print(f"      c = {fit_params[8]:.6f} +- {param_errors[8]:.6f}")
    print("-" * 60)
    print("Derived quantities:")
    if effective_mean_err is not None:
        print(f"Effective Mean: {effective_mean:.6f} +- {effective_mean_err:.6f}")
    else:
        print(f"Effective Mean: {effective_mean:.6f}")
    
    if effective_sigma_err is not None:
        print(f"Effective Sigma: {effective_sigma:.6f} +- {effective_sigma_err:.6f}")
    else:
        print(f"Effective Sigma: {effective_sigma:.6f}")
    print("=" * 60)


    if plot:
        
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [2, 1]})
        
        if result.success:
            
            fit_values = double_gaussian_and_pol_2(bin_centers, fit_params)
            residuals = counts - fit_values
            pulls = residuals / errors

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
            
            x_fit = np.linspace(bin_edges[0], bin_edges[-1], 200)
            y_fit = double_gaussian_and_pol_2(x_fit, fit_params)
            
            axes[0].plot(
                x_fit,
                y_fit,
                color='red',
                linestyle='-',
                linewidth=2, 
                label='Double Gaussian Fit'
            )
            
            g1 = fit_params[0] / (np.sqrt(2*np.pi) * fit_params[2]) * np.exp( -0.5 * ((x_fit - fit_params[1]) / fit_params[2])**2 )
            g2 = fit_params[3] / (np.sqrt(2*np.pi) * fit_params[5]) * np.exp( -0.5 * ((x_fit - fit_params[4]) / fit_params[5])**2 )
            pol2 = fit_params[6] + fit_params[7] * x_fit + fit_params[8] * x_fit**2

            g1_label = 'Gaussian 1'
            g2_label = 'Gaussian 2'
            pol2_label = 'Pol2'
            
            axes[0].plot(x_fit, g1, 'g:', label=g1_label, alpha=0.7)
            axes[0].plot(x_fit, g2, 'b:', label=g2_label, alpha=0.7)
            axes[0].plot(x_fit, pol2, 'm:', label=pol2_label, alpha=0.7)
            
            def x_scaler(x, pos):
                return f'{x * scaler:.{scaler_round}f}'
            
            axes[0].xaxis.set_major_formatter(FuncFormatter(x_scaler))
            axes[1].xaxis.set_major_formatter(FuncFormatter(x_scaler))
            
            fit_text = ''

            if effective_sigma_err is not None and not np.isnan(effective_sigma_err):
                fit_text += f'Effective $\sigma$: ({effective_sigma * scaler:.4f} $\pm$ {effective_sigma_err * scaler:.4f}) $\mu$m'
            else:
                fit_text += f'Effective $\sigma$: {effective_sigma * scaler:.4f} $\mu$m'

            if effective_mean_err is not None and not np.isnan(effective_mean_err):
                fit_text += f'\nEffective $\mu$: ({effective_mean * scaler:.4f} $\pm$ {effective_mean_err * scaler:.4f}) $\mu$m'
            else:
                fit_text += f'\nEffective $\mu$: {effective_mean * scaler:.4f} $\mu$m'

            axes[0].text(0.02, 0.7, fit_text, transform=axes[0].transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=10)
            
            axes[0].set_title(f'{title}', fontsize=14, fontweight='bold', pad=20)
            axes[0].set_ylabel('Counts')
            axes[0].legend(loc='best')
            axes[0].grid(True, alpha=0.5, linestyle='--')
            
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            
            if log_y:
                axes[0].set_yscale('log')
            
            axes[1].axhline(y=0, color='red', linestyle='-', linewidth=1, alpha=0.7)
            axes[1].axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=-1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            axes[1].axhline(y=-2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            
            bin_width = bin_edges[1] - bin_edges[0]
            axes[1].bar(bin_centers, pulls, 
                        width=bin_width*0.8,
                        color='skyblue', 
                        edgecolor='darkblue',
                        alpha=0.7, 
                        label='Pulls')
            
            mean_pull = np.mean(pulls)
            std_pull = np.std(pulls)
            
            pull_stats = (
                f'Pull Statistics:'
                f'\n$\mu_{{pull}}$ = {mean_pull:.3f}'
                f'\n$\sigma_{{pull}}$ = {std_pull:.3f}'
            )
            
            axes[1].text(0.02, 0.98, pull_stats, transform=axes[1].transAxes,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8), fontsize=9)
            
            axes[1].set_ylabel(r'$(data - fit) / \sqrt{data}$')
            axes[1].legend(loc='upper right')
            
            axes[1].set_xlabel(f'{x_label}')
            axes[1].grid(True, alpha=0.5, linestyle='--')
            
            axes[1].spines['top'].set_visible(False)
            axes[1].spines['right'].set_visible(False)
        
        else:
            axes[0].text(0.5, 0.5, 'Fit failed!', transform=axes[0].transAxes, 
                    ha='center', va='center', fontsize=12, color='red')    
        
        plt.tight_layout()
        
        fit_results = {
            'params': fit_params,
            'errors': param_errors,
            'cov_matrix': cov_matrix,
            'chi2': chi2,
            'reduced_chi2': reduced_chi2,
            'effective_sigma': effective_sigma,
            'effective_sigma_err': effective_sigma_err,
            'effective_mean': effective_mean,
            'effective_mean_err': effective_mean_err,
            'result': result
        }
        
        return fit_results, fig, axes
    
    fit_results = {
        'params': fit_params,
        'errors': param_errors,
        'cov_matrix': cov_matrix,
        'chi2': chi2,
        'reduced_chi2': reduced_chi2,
        'effective_sigma': effective_sigma,
        'effective_sigma_err': effective_sigma_err,
        'effective_mean': effective_mean,
        'effective_mean_err': effective_mean_err,
        'result': result
    }
    
    return fit_results


###########################################################################################################
# fit_distr_double_gauss_and_const
###########################################################################################################

def fit_distr_double_gauss_and_const(
    bin_centers,
    bin_edges,
    counts,
    errors,
    title: str,
    x_label: str,
    optimize_method: str = 'L-BFGS-B',
    plot: bool = False,
    log_y: bool = False,
    max_iter: int = 1000,
    scaler: float = 1e4,
    scaler_round: int = 0,
    unit_name: str = '$\mu m$'
) -> Tuple:
    
    version = 1.5
        
    if version != CFG.fit_version:
        raise 'Versions of fit and s/b calculations might be different!!'
    
    mask = counts > 0
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    errors = errors[mask]

    #=======================================================================================
    # Fit Model - Double Gaussian and Pol2
    #=======================================================================================
    
    def double_gaussian_and_const(x, params):
        """
        Composite model of 2 Gaussian functions with const background
        
        Args:
            x: array-like, input values
            params: 7 parameters [A1, mu1, sigma1, A2, mu2, sigma2, C]
        
        Returns:
            y: sum of 2 Gaussian functions and const
        """
        A1, mu1, sigma1, A2, mu2, sigma2, C = params
        
        # g1 = A1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
        # g2 = A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
        # C
        
        g1 = A1 / (np.sqrt(2*np.pi) * sigma1) * np.exp( -0.5 * ((x - mu1) / sigma1)**2 )
        g2 = A2 / (np.sqrt(2*np.pi) * sigma2) * np.exp( -0.5 * ((x - mu2) / sigma2)**2 )
        
        return g1 + g2 + C
    
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])
    
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    initial_guess = [
        initial_amplitude, data_min + 0.25 * data_range, 0.1 * data_range,  # Gaussian 1
        initial_amplitude, data_min + 0.75 * data_range, 0.1 * data_range,  # Gaussian 2  
        initial_amplitude   # const
    ]
    
    # Set bounds for parameters
    bounds = [
        (0, None), (data_min, data_max), (1e-6, 0.5 * data_range),  # Gaussian 1
        (0, None), (data_min, data_max), (1e-6, 0.5 * data_range),  # Gaussian 2
        (0, None)   # Const
    ]
    
    def chi2_function(params, x, y, errors):
 
        predictions = double_gaussian_and_const(x, params)
 
        residuals = y - predictions
 
        return np.sum((residuals / errors) ** 2)
        
        
    # =========================================================================
    # Fit
    # =========================================================================
 
    result = minimize(
        fun=chi2_function,
        x0=initial_guess,
        args=(bin_centers, counts, errors),
        method=optimize_method,
        bounds=bounds,
        options={"maxiter": max_iter}
    )
 
    if not result.success:
        print("Warning: fit did not fully converge:")
        print(result.message)
 
    fit_params = result.x
 
    n_points = len(counts)
    n_params = len(fit_params)
    n_dof = n_points - n_params
        
    # =========================================================================
    # Numerical Hessian
    # =========================================================================
 
    def numerical_hessian(fun, x0, epsilon=1e-4):
        """
        Numerical Hessian using central finite differences.
        """
 
        x0 = np.asarray(x0, dtype=float)
 
        n = len(x0)
        hess = np.zeros((n, n))
 
        steps = epsilon * np.maximum(np.abs(x0), 1.0)
 
        f0 = fun(x0)
 
        for i in range(n):
 
            hi = steps[i]
 
            ei = np.zeros(n)
            ei[i] = hi
 
            # Diagonal second derivative
            f_plus = fun(x0 + ei)
            f_minus = fun(x0 - ei)
 
            hess[i, i] = (
                f_plus - 2.0 * f0 + f_minus
            ) / hi**2
 
            for j in range(i + 1, n):
 
                hj = steps[j]
 
                ej = np.zeros(n)
                ej[j] = hj
 
                f_pp = fun(x0 + ei + ej)
                f_pm = fun(x0 + ei - ej)
                f_mp = fun(x0 - ei + ej)
                f_mm = fun(x0 - ei - ej)
 
                value = (
                    f_pp
                    - f_pm
                    - f_mp
                    + f_mm
                ) / (4.0 * hi * hj)
 
                hess[i, j] = value
                hess[j, i] = value
 
        return hess
        
    # =========================================================================
    # Covariance matrix and parameter errors
    # =========================================================================
 
    def calculate_parameter_errors(
        result,
        x,
        y,
        errors,
        n_dof
    ):
        """
        Calculate parameter covariance matrix.
 
        Since the Hessian is calculated for chi2:
 
            H = d²(chi2) / dtheta_i dtheta_j
 
        the covariance matrix in the quadratic approximation is:
 
            Cov = 2 * H^{-1}
 
        If the uncertainties are only known up to a common scale,
        the covariance is additionally scaled by reduced chi2.
        """
 
        objective = lambda p: chi2_function(
            p, x, y, errors
        )
 
        hess = numerical_hessian(
            objective,
            result.x,
            epsilon=1e-4
        )
 
        # Symmetrize to reduce numerical noise
        hess = 0.5 * (hess + hess.T)
 
        chi2 = objective(result.x)
 
        reduced_chi2 = (
            chi2 / n_dof
            if n_dof > 0
            else np.nan
        )
 
        try:
 
            # Covariance for Hessian of chi2
            cov_matrix = 2.0 * np.linalg.inv(hess)
 
            # Optional scale correction
            #
            # Keep this if errors are estimates and the overall
            # noise scale may be incorrect.
            #
            # If errors are known exact Poisson errors, one may
            # choose NOT to multiply by reduced_chi2.
            
            # Maybe remove
            if n_dof > 0:
                cov_matrix *= reduced_chi2
 
            # Symmetrize covariance
            cov_matrix = 0.5 * (
                cov_matrix + cov_matrix.T
            )
 
            variances = np.diag(cov_matrix)
 
            # Negative values can appear from numerical noise
            param_errors = np.sqrt(
                np.maximum(variances, 0.0)
            )
 
            return (
                param_errors,
                cov_matrix,
                chi2,
                reduced_chi2,
                hess
            )
 
        except np.linalg.LinAlgError:
 
            print(
                "Warning: Hessian is singular. "
                "Covariance matrix cannot be calculated."
            )
 
            n = len(result.x)
 
            return (
                np.full(n, np.nan),
                np.full((n, n), np.nan),
                chi2,
                reduced_chi2,
                hess
            )
 
    param_errors, cov_matrix, chi2, reduced_chi2, hess = calculate_parameter_errors(
        result=result,
        x=bin_centers,
        y=counts,
        errors=errors,
        n_dof=n_dof
    )
    
    # =========================================================================
    # Effective parameters of the Gaussian mixture
    # =========================================================================
 
    def calculate_effective_parameters(
        fit_params,
        cov_matrix
    ):
        """
        Calculate mean and standard deviation of the full
        two-Gaussian mixture.
 
        Parameters:
            [A1, mu1, sigma1, A2, mu2, sigma2, C]
 
        For the mixture weights:
 
            w1 = A1 / (A1 + A2)
            w2 = A2 / (A1 + A2)
 
        Mean:
 
            mu_eff = w1*mu1 + w2*mu2
 
        Full variance of the mixture:
 
            Var_eff =
                w1*sigma1²
                + w2*sigma2²
                + w1*w2*(mu1 - mu2)²
 
        Equivalently:
 
            Var_eff =
                w1*(sigma1² + (mu1 - mu_eff)²)
                + w2*(sigma2² + (mu2 - mu_eff)²)
 
        Errors are propagated using the full covariance matrix:
 
            Var(f) = grad(f)^T Cov grad(f)
        """
 
        A1, mu1, sigma1, A2, mu2, sigma2, C = fit_params
 
        total_amplitude = A1 + A2
 
        if total_amplitude <= 0:
            return np.nan, np.nan, np.nan, np.nan
 
        # ---------------------------------------------------------------------
        # Weights
        # ---------------------------------------------------------------------
 
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
 
        # ---------------------------------------------------------------------
        # Effective mean
        # ---------------------------------------------------------------------
 
        effective_mean = (
            w1 * mu1
            + w2 * mu2
        )
 
        # ---------------------------------------------------------------------
        # Full effective variance
        # ---------------------------------------------------------------------
 
        delta_mu = mu1 - mu2
 
        effective_variance = (
            w1 * sigma1**2
            + w2 * sigma2**2
            + w1 * w2 * delta_mu**2
        )
 
        effective_variance = max(
            effective_variance,
            0.0
        )
 
        effective_sigma = np.sqrt(
            effective_variance
        )
 
        # ---------------------------------------------------------------------
        # Error propagation
        # ---------------------------------------------------------------------
 
        if (
            cov_matrix is None
            or np.any(~np.isfinite(cov_matrix))
        ):
            return (
                effective_sigma,
                effective_mean,
                np.nan,
                np.nan
            )
 
        # =====================================================================
        # Gradient of effective mean
        # =====================================================================
 
        grad_mean = np.zeros(len(fit_params))
 
        grad_mean[0] = (
            A2 * (mu1 - mu2)
            / total_amplitude**2
        )
 
        grad_mean[1] = w1
 
        grad_mean[2] = 0.0
 
        grad_mean[3] = (
            A1 * (mu2 - mu1)
            / total_amplitude**2
        )
 
        grad_mean[4] = w2
 
        grad_mean[5] = 0.0
 
        grad_mean[6] = 0.0
 
        # Full covariance propagation
        mean_variance = (
            grad_mean
            @ cov_matrix
            @ grad_mean
        )
 
        effective_mean_err = np.sqrt(
            max(mean_variance, 0.0)
        )
 
        # =====================================================================
        # Gradient of full effective variance
        #
        # V =
        #   w1*sigma1²
        #   + w2*sigma2²
        #   + w1*w2*(mu1-mu2)²
        #
        # =====================================================================
 
        grad_var = np.zeros(len(fit_params))
 
        # Derivatives with respect to A1 and A2
 
        # dw1/dA1 = A2 / (A1+A2)^2
        # dw2/dA1 = -A2 / (A1+A2)^2
 
        # dw1/dA2 = -A1 / (A1+A2)^2
        # dw2/dA2 = A1 / (A1+A2)^2
 
        dV_dw1 = (
            sigma1**2
            + w2 * delta_mu**2
        )
 
        dV_dw2 = (
            sigma2**2
            + w1 * delta_mu**2
        )
 
        grad_var[0] = (
            A2 / total_amplitude**2
        ) * (dV_dw1 - dV_dw2)
 
        grad_var[3] = (
            A1 / total_amplitude**2
        ) * (dV_dw2 - dV_dw1)
 
        # Derivatives with respect to mu1 and mu2
 
        grad_var[1] = (
            2.0
            * w1
            * w2
            * delta_mu
        )
 
        grad_var[4] = (
            -2.0
            * w1
            * w2
            * delta_mu
        )
 
        # Derivatives with respect to sigma1 and sigma2
 
        grad_var[2] = (
            2.0 * w1 * sigma1
        )
 
        grad_var[5] = (
            2.0 * w2 * sigma2
        )
 
        # Constant background does not enter effective moments
        grad_var[6] = 0.0
 
        # Error of variance
        variance_variance = (
            grad_var
            @ cov_matrix
            @ grad_var
        )
 
        variance_err = np.sqrt(
            max(variance_variance, 0.0)
        )
 
        # =====================================================================
        # sigma = sqrt(V)
        #
        # d sigma / d V = 1 / (2 sqrt(V))
        # =====================================================================
 
        if effective_sigma > 0:
 
            effective_sigma_err = (
                variance_err
                / (2.0 * effective_sigma)
            )
 
        else:
 
            effective_sigma_err = np.nan
 
        return (
            effective_sigma,
            effective_mean,
            effective_sigma_err,
            effective_mean_err
        )
        
    effective_sigma, effective_mean, effective_sigma_err, effective_mean_err = calculate_effective_parameters(
        fit_params,
        cov_matrix
    )

    print("=" * 60)
    print("Optimization result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of iterations: {result.nit}")
    print(f"Final objective value: {result.fun:.6f}")
    print(f"Chi-squared: {chi2:.2f}")
    print(f"Reduced chi-squared: {reduced_chi2:.3f}")
    print("-" * 60)
    print("Fit parameters with errors:")
    print(f"Gauss 1: A = {fit_params[0]:.4f} +- {param_errors[0]:.4f}")
    print(f"         mu = {fit_params[1]:.6f} +- {param_errors[1]:.6f}")
    print(f"         sigma = {fit_params[2]:.6f} +- {param_errors[2]:.6f}")
    print(f"Gauss 2: A = {fit_params[3]:.4f} +- {param_errors[3]:.4f}")
    print(f"         mu = {fit_params[4]:.6f} +- {param_errors[4]:.6f}")
    print(f"         sigma = {fit_params[5]:.6f} +- {param_errors[5]:.6f}")
    print(f"Const: C = {fit_params[6]:.4f} +- {param_errors[6]:.4f}")
    print("-" * 60)
    print("Derived quantities:")
    if effective_mean_err is not None:
        print(f"Effective Mean: {effective_mean:.6f} +- {effective_mean_err:.6f}")
    else:
        print(f"Effective Mean: {effective_mean:.6f}")
    
    if effective_sigma_err is not None:
        print(f"Effective Sigma: {effective_sigma:.6f} +- {effective_sigma_err:.6f}")
    else:
        print(f"Effective Sigma: {effective_sigma:.6f}")
    print("=" * 60)


    if plot:
        
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [2, 1]})
        
        if result.success:
            
            fit_values = double_gaussian_and_const(bin_centers, fit_params)
            residuals = counts - fit_values
            pulls = residuals / errors
            
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
            
            x_fit = np.linspace(bin_edges[0], bin_edges[-1], 200)
            y_fit = double_gaussian_and_const(x_fit, fit_params)
            
            axes[0].plot(
                x_fit,
                y_fit,
                color='red',
                linestyle='-',
                linewidth=2, 
                label='Double Gaussian Fit'
            )
            
            g1 = fit_params[0] / (np.sqrt(2*np.pi) * fit_params[2]) * np.exp( -0.5 * ((x_fit - fit_params[1]) / fit_params[2])**2 )
            g2 = fit_params[3] / (np.sqrt(2*np.pi) * fit_params[5]) * np.exp( -0.5 * ((x_fit - fit_params[4]) / fit_params[5])**2 )
            const = fit_params[6] * np.ones(shape=g1.size)

            g1_label = 'Gaussian 1'
            g2_label = 'Gaussian 2'
            const_label = 'Constant'
            
            axes[0].plot(x_fit, g1, 'g:', label=g1_label, alpha=0.7)
            axes[0].plot(x_fit, g2, 'b:', label=g2_label, alpha=0.7)
            axes[0].plot(x_fit, const, 'm:', label=const_label, alpha=0.7)
            
            def x_scaler(x, pos):
                return f'{x * scaler:.{scaler_round}f}'

            axes[0].xaxis.set_major_formatter(FuncFormatter(x_scaler))
            axes[1].xaxis.set_major_formatter(FuncFormatter(x_scaler))
            
            fit_text = ''

            if effective_sigma_err is not None and not np.isnan(effective_sigma_err):
                fit_text += f'Effective $\sigma$: ({effective_sigma * scaler:.4f} $\pm$ {effective_sigma_err * scaler:.4f}) $\mu$m'
            else:
                fit_text += f'Effective $\sigma$: {effective_sigma * scaler:.4f} $\mu$m'

            if effective_mean_err is not None and not np.isnan(effective_mean_err):
                fit_text += f'\nEffective $\mu$: ({effective_mean * scaler:.4f} $\pm$ {effective_mean_err * scaler:.4f}) $\mu$m'
            else:
                fit_text += f'\nEffective $\mu$: {effective_mean * scaler:.4f} $\mu$m'

            axes[0].text(0.02, 0.7, fit_text, transform=axes[0].transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=10)
            
            axes[0].set_title(f'{title}', fontsize=14, fontweight='bold', pad=20)
            axes[0].set_ylabel('Counts')
            axes[0].legend(loc='best')
            axes[0].grid(True, alpha=0.5, linestyle='--')
            
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            
            if log_y:
                axes[0].set_yscale('log')
            
            axes[1].axhline(y=0, color='red', linestyle='-', linewidth=1, alpha=0.7)
            axes[1].axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=-1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            axes[1].axhline(y=-2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            
            bin_width = bin_edges[1] - bin_edges[0]
            axes[1].bar(bin_centers, pulls, 
                        width=bin_width*0.8,
                        color='skyblue', 
                        edgecolor='darkblue',
                        alpha=0.7, 
                        label='Pulls')
            
            mean_pull = np.mean(pulls)
            std_pull = np.std(pulls)
            
            pull_stats = (
                f'Pull Statistics:'
                f'\n$\mu_{{pull}}$ = {mean_pull:.3f}'
                f'\n$\sigma_{{pull}}$ = {std_pull:.3f}'
                )
            
            axes[1].text(0.02, 0.98, pull_stats, transform=axes[1].transAxes,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8), fontsize=9)
            
            axes[1].set_ylabel(r'$(data - fit) / \sqrt{data}$')
            axes[1].legend(loc='upper right')
            
            axes[1].set_xlabel(f'{x_label}')
            axes[1].grid(True, alpha=0.5, linestyle='--')
            
            axes[1].spines['top'].set_visible(False)
            axes[1].spines['right'].set_visible(False)
        
        else:
            axes[0].text(0.5, 0.5, 'Fit failed!', transform=axes[0].transAxes, 
                    ha='center', va='center', fontsize=12, color='red')    
        
        plt.tight_layout()
        

        fit_results = {
            'params': fit_params,
            'errors': param_errors,
            'cov_matrix': cov_matrix,
            'chi2': chi2,
            'reduced_chi2': reduced_chi2,
            'effective_sigma': effective_sigma,
            'effective_sigma_err': effective_sigma_err,
            'effective_mean': effective_mean,
            'effective_mean_err': effective_mean_err,
            'result': result
        }
        
        return fit_results, fig, axes
    
    fit_results = {
        'params': fit_params,
        'errors': param_errors,
        'cov_matrix': cov_matrix,
        'chi2': chi2,
        'reduced_chi2': reduced_chi2,
        'effective_sigma': effective_sigma,
        'effective_sigma_err': effective_sigma_err,
        'effective_mean': effective_mean,
        'effective_mean_err': effective_mean_err,
        'result': result
    }
    
    return fit_results


###########################################################################################################
# fit_distr_double_gauss_and_pol2_and_exp
###########################################################################################################

def fit_distr_double_gauss_and_pol2_end_exp(
    bin_centers,
    bin_edges,
    counts,
    errors,
    title: str,
    x_label: str,
    optimize_method: str = 'L-BFGS-B',
    plot: bool = False,
    log_y: bool = False,
    max_iter: int = 1000,
    scaler: float = 1e4,
    unit_name: str = '$\mu m$') -> Tuple:
    
    version = 1.4
        
    if version != CFG.fit_version:
        raise 'Versions of fit and s/b calculations might be different!!'
    
    mask = counts > 0
    bin_centers = bin_centers[mask]
    counts = counts[mask]
    errors = errors[mask]

    #=======================================================================================
    # Fit Model - Double Gaussian with Constant
    #=======================================================================================
    
    def double_gaussian_and_pol_2_and_exp(x, *params):
        """
        Composite model of 2 Gaussian functions with pol2 background
        
        Args:
            x: array-like, input values
            params: 7 parameters [A1, mu1, sigma1, A2, mu2, sigma2, A3, b, c, A4, e]
        
        Returns:
            y: sum of 2 Gaussian functions and pol2 and exp
        """
        A1, mu1, sigma1, A2, mu2, sigma2, A3, b, c, A4, e = params
        
        g1 = A1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
        g2 = A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
        pol2 = A3 + b * x + c * x**2
        exp = A4 * np.exp(-x * e)
        
        return g1 + g2 + pol2 + exp
    
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])
    
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    initial_guess = [
        initial_amplitude, data_min + 0.25 * data_range, 0.1 * data_range,  # Gaussian 1
        initial_amplitude, data_min + 0.75 * data_range, 0.1 * data_range,  # Gaussian 2  
        initial_amplitude, 1, -2,   # Pol2
        initial_amplitude, 1,   # Exp
    ]
    
    # Set bounds for parameters
    bounds = [
        # (0, None), (2.284, 2.287), (1e-6, 0.5 * data_range),  # Gaussian 1
        # (0, None), (2.284, 2.287), (1e-6, 0.5 * data_range),  # Gaussian 2
        (0, None), (2.284, 2.287), (1e-6, 0.02),  # Gaussian 1
        (0, None), (2.284, 2.287), (1e-6, 0.02),  # Gaussian 2
        (0, None), (None, None), (None, None),   # Pol2
        (0, None), (None, None) # Exp
    ]
    
    
    def calculate_effective_sigma_double_gaussian(fit_params, fit_errors=None):
        """
        Calculate effective sigma for double Gaussian mixture
        
        For a mixture distribution: f(x) = w1*g1(x) + w2*g2(x)
        where w_i = A_i / (A1 + A2) are the weights
        
        The variance of the mixture is:
        sigma_eff = sqrt(w1 * sigma1**2 + w2*sigma2**2)
        
        Optionally propagates errors if fit_errors is provided
        """
        A1, mu1, sigma1, A2, mu2, sigma2, A3, b, c, A4, e = fit_params
        
        # Weights
        total_amplitude = A1 + A2
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
        
        effective_mean = w1 * mu1 + w2 * mu2
        effective_variance = w1 * sigma1**2 + w2 * sigma2**2
        effective_sigma = np.sqrt(effective_variance)
        
        # Uncertainties
        effective_mean_err = None
        effective_sigma_err = None
        
        if fit_errors is not None:
            A1_err, mu1_err, sigma1_err, A2_err, mu2_err, sigma2_err, A3_err, b_err, c_err, A4_err, e_err = fit_errors
            
            # Partial derivatives for effective mean
            dmean_dA1 = (mu1 - mu2) * A2 / (total_amplitude**2)
            dmean_dA2 = (mu2 - mu1) * A1 / (total_amplitude**2)
            dmean_dmu1 = w1
            dmean_dmu2 = w2
            
            effective_mean_err = np.sqrt(
                (dmean_dA1 * A1_err)**2 + 
                (dmean_dA2 * A2_err)**2 + 
                (dmean_dmu1 * mu1_err)**2 + 
                (dmean_dmu2 * mu2_err)**2
            )
            
            # Partial derivatives for effective sigma^2
            # f_sigma2 = w1*sigma1^2 + w2*sigma2^2
            df_sigma2_dA1 = (sigma1**2 - sigma2**2) * A2 / (total_amplitude**2)
            df_sigma2_dA2 = (sigma2**2 - sigma1**2) * A1 / (total_amplitude**2)
            df_sigma2_dsigma1 = 2 * w1 * sigma1
            df_sigma2_dsigma2 = 2 * w2 * sigma2
            
            sigma2_err_sq = (
                (df_sigma2_dA1 * A1_err)**2 + 
                (df_sigma2_dA2 * A2_err)**2 + 
                (df_sigma2_dsigma1 * sigma1_err)**2 + 
                (df_sigma2_dsigma2 * sigma2_err)**2
            )
            
            # Error for effective sigma: d(sqrt(f))/dx = 0.5 * (1/sqrt(f)) * df/dx
            if effective_variance > 0:
                effective_sigma_err = 0.5 * np.sqrt(sigma2_err_sq) / np.sqrt(effective_variance)
        
        return effective_sigma, effective_mean, effective_sigma_err, effective_mean_err
    
    #=======================================================================================
    # Loss Function and Optimization
    #=======================================================================================
    
    def weighted_loss(params, x, y, errors):
        """Objective function to minimize (sum of weighted squared residuals)"""
        predictions = double_gaussian_and_pol_2_and_exp(x, *params)
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

    n_params = len(fit_params)
    n_points = len(errors)
    n_dof = n_points - n_params

    #=======================================================================================
    # Calculate Parameter Errors
    #=======================================================================================
    
    def calculate_parameter_errors(result, x, y, errors):
        """
        Calculate errors for fitted parameters using the covariance matrix
        
        The covariance matrix is approximated by the inverse of the Hessian matrix
        at the minimum, scaled by the reduced chi-squared
        """
        
        # Calculate Hessian matrix using finite differences
        def hessian_finite_diff(fun, x0, args=(), epsilon=1e-4):
            """Calculate Hessian matrix using central finite differences"""
            n = len(x0)
            hess = np.zeros((n, n))
            
            for i in range(n):
                for j in range(i, n):
                    # Create basis vectors
                    ei = np.zeros(n)
                    ej = np.zeros(n)
                    ei[i] = 1.0
                    ej[j] = 1.0
                    
                    # Central difference for Hessian
                    if i == j:
                        # Diagonal elements
                        f_plus = fun(x0 + epsilon * ei, *args)
                        f_minus = fun(x0 - epsilon * ei, *args)
                        f0 = fun(x0, *args)
                        hess[i, i] = (f_plus - 2*f0 + f_minus) / (epsilon**2)
                    else:
                        # Off-diagonal elements
                        f_plus_plus = fun(x0 + epsilon * ei + epsilon * ej, *args)
                        f_plus_minus = fun(x0 + epsilon * ei - epsilon * ej, *args)
                        f_minus_plus = fun(x0 - epsilon * ei + epsilon * ej, *args)
                        f_minus_minus = fun(x0 - epsilon * ei - epsilon * ej, *args)
                        hess[i, j] = (f_plus_plus - f_plus_minus - f_minus_plus + f_minus_minus) / (4 * epsilon**2)
                        hess[j, i] = hess[i, j]
            
            return hess
        
        # Calculate Hessian at minimum
        hess = hessian_finite_diff(
            lambda p: weighted_loss(p, x, y, errors), 
            result.x, 
            epsilon=1e-4
        )
        
        # Calculate covariance matrix (inverse of Hessian)
        try:
            cov_matrix = np.linalg.inv(hess)
            
            # Calculate chi-squared and reduced chi-squared
            predictions = double_gaussian_and_pol_2_and_exp(x, *result.x)
            residuals = y - predictions
            weights = 1.0 / (errors**2 + 1e-6)
            chi2 = np.sum(weights * residuals**2)
            reduced_chi2 = chi2 / n_dof if n_dof > 0 else chi2
            
            # Scale covariance matrix by reduced chi-squared
            cov_matrix *= reduced_chi2
            
            # Parameter errors are sqrt of diagonal elements
            param_errors = np.sqrt(np.diag(cov_matrix))
            
            return param_errors, cov_matrix, chi2, reduced_chi2
            
        except np.linalg.LinAlgError:
            # If matrix inversion fails, return NaNs
            print("Warning: Could not calculate covariance matrix (singular Hessian)")
            n_params = len(result.x)
            return np.full(n_params, np.nan), np.full((n_params, n_params), np.nan), np.nan, np.nan
    
    # Calculate parameter errors
    param_errors, cov_matrix, chi2, reduced_chi2 = calculate_parameter_errors(
        result, bin_centers, counts, errors
    )
    
    # Calculate effective sigma with error propagation
    effective_sigma, effective_mean, effective_sigma_err, effective_mean_err = calculate_effective_sigma_double_gaussian(
        fit_params, param_errors
    )

    print("=" * 60)
    print("Optimization result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Number of iterations: {result.nit}")
    print(f"Final objective value: {result.fun:.6f}")
    print(f"Chi-squared: {chi2:.2f}")
    print(f"Reduced chi-squared: {reduced_chi2:.3f}")
    print("-" * 60)
    print("Fit parameters with errors:")
    print(f"Gauss 1: A = {fit_params[0]:.4f} +- {param_errors[0]:.4f}")
    print(f"         mu = {fit_params[1]:.6f} +- {param_errors[1]:.6f}")
    print(f"         sigma = {fit_params[2]:.6f} +- {param_errors[2]:.6f}")
    print(f"Gauss 2: A = {fit_params[3]:.4f} +- {param_errors[3]:.4f}")
    print(f"         mu = {fit_params[4]:.6f} +- {param_errors[4]:.6f}")
    print(f"         sigma = {fit_params[5]:.6f} +- {param_errors[5]:.6f}")
    print(f"Pol2: A = {fit_params[6]:.4f} +- {param_errors[6]:.4f}")
    print(f"      b = {fit_params[7]:.6f} +- {param_errors[7]:.6f}")
    print(f"      c = {fit_params[8]:.6f} +- {param_errors[8]:.6f}")
    print(f"Exp: A = {fit_params[9]:.4f} +- {param_errors[9]:.4f}")
    print(f"      b = {fit_params[10]:.6f} +- {param_errors[10]:.6f}")
    print("-" * 60)
    print("Derived quantities:")
    if effective_mean_err is not None:
        print(f"Effective Mean: {effective_mean:.6f} +- {effective_mean_err:.6f}")
    else:
        print(f"Effective Mean: {effective_mean:.6f}")
    
    if effective_sigma_err is not None:
        print(f"Effective Sigma: {effective_sigma:.6f} +- {effective_sigma_err:.6f}")
    else:
        print(f"Effective Sigma: {effective_sigma:.6f}")
    print("=" * 60)


    if plot:
        
        fig, axes = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [2, 1]})
        
        if result.success:
            
            # Calculate fit values and pulls
            fit_values = double_gaussian_and_pol_2_and_exp(bin_centers, *fit_params)
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
            y_fit = double_gaussian_and_pol_2_and_exp(x_fit, *fit_params)
            
            axes[0].plot(
                x_fit,
                y_fit,
                color='red',
                linestyle='-',
                linewidth=2, 
                label='Double Gaussian Fit'
            )
            
            # Plot individual Gaussians
            g1 = fit_params[0] * np.exp(-0.5 * ((x_fit - fit_params[1]) / fit_params[2]) ** 2)
            g2 = fit_params[3] * np.exp(-0.5 * ((x_fit - fit_params[4]) / fit_params[5]) ** 2)
            pol2 = fit_params[6] + fit_params[7] * x_fit + fit_params[8] * x_fit**2
            exp = fit_params[9] * np.exp(-x_fit * fit_params[10])

            g1_label = 'Gaussian 1'
            g2_label = 'Gaussian 2'
            pol2_label = 'Pol2'
            exp_label = 'Exp'
            
            axes[0].plot(x_fit, g1, 'g:', label=g1_label, alpha=0.7)
            axes[0].plot(x_fit, g2, 'b:', label=g2_label, alpha=0.7)
            axes[0].plot(x_fit, pol2, 'm:', label=pol2_label, alpha=0.7)
            axes[0].plot(x_fit, exp, 'y:', label=exp_label, alpha=0.7)
                        
            # Add fit parameters as text with errors
            if not np.isnan(param_errors[0]):
                fit_text = (fr'Fit Parameters:'
                          fr'\\G1: A={fit_params[0]:.1f} $\pm$ {param_errors[0]:.1f},\\ $\mu$={fit_params[1] * scaler:.3f} $\pm$ {param_errors[1] * scaler:.3f} {unit_name},\\ $\sigma$={fit_params[2] * scaler:.3f} $\pm$ {param_errors[2] * scaler:.3f} {unit_name}'
                          fr'\\G2: A={fit_params[3]:.1f} $\pm$ {param_errors[3]:.1f},\\ $\mu$={fit_params[4] * scaler:.3f} $\pm$ {param_errors[4] * scaler:.3f} {unit_name},\\ $\sigma$={fit_params[5] * scaler:.3f} $\pm$ {param_errors[5] * scaler:.3f} {unit_name}'
                          fr'\\Pol2: A={fit_params[6]:.1f} $\pm$ {param_errors[6]:.1f},\\ $b$={fit_params[7] * scaler:.3f} $\pm$ {param_errors[7] * scaler:.3f} {unit_name},\\ c={fit_params[8] * scaler:.3f} $\pm$ {param_errors[8] * scaler:.3f} {unit_name}'
                          fr'\\Exp: A={fit_params[9]:.1f} $\pm$ {param_errors[9]:.1f},\\ $e$={fit_params[10] * scaler:.3f} $\pm$ {param_errors[10] * scaler:.3f} {unit_name}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            else:
                fit_text = (fr'Fit Parameters:'
                          fr'\\G1: A={fit_params[0]:.1f},\\ $\mu$={fit_params[1] * scaler:.3f} {unit_name},\\ $\sigma$={fit_params[2] * scaler:.3f} {unit_name}'
                          fr'\\G2: A={fit_params[3]:.1f},\\ $\mu$={fit_params[4] * scaler:.3f} {unit_name},\\ $\sigma$={fit_params[5] * scaler:.3f} {unit_name}'
                          fr'\\Pol2: A={fit_params[6]:.1f},\\ $b$={fit_params[7] * scaler:.3f} {unit_name},\\ c={fit_params[8] * scaler:.3f} {unit_name}'
                          fr'\\Exp: A={fit_params[9]:.1f},\\ $e$={fit_params[10] * scaler:.3f} {unit_name}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            
            if effective_sigma_err is not None and not np.isnan(effective_sigma_err):
                fit_text += fr'\\Effective $\sigma$: {effective_sigma * scaler:.4f} $\pm$ {effective_sigma_err * scaler:.4f} {unit_name}'
            else:
                fit_text += fr'\\Effective $\sigma$: {effective_sigma * scaler:.4f} {unit_name}'
            
            if effective_sigma_err is not None and not np.isnan(effective_sigma_err):
                fit_text += fr'\\Effective mean: {effective_mean * scaler:.4f} $\pm$ {effective_mean_err * scaler:.4f} {unit_name}'
            else:
                fit_text += fr'\\Effective mean: {effective_mean * scaler:.4f} {unit_name}'

            axes[0].text(0.02, 0.7, fit_text, transform=axes[0].transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=12)
            
            axes[0].set_title(f'{title}', fontsize=14, fontweight='bold', pad=20)
            axes[0].set_ylabel('Counts', fontweight='bold')
            axes[0].legend(loc='best')
            axes[0].grid(True, alpha=0.5, linestyle='--')
            
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            
            if log_y:
                axes[0].set_yscale('log')
            
            # Lower plot: Pulls
            axes[1].axhline(y=0, color='red', linestyle='-', linewidth=1, alpha=0.7)
            axes[1].axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=-1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            axes[1].axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            axes[1].axhline(y=-2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
            
            bin_width = bin_edges[1] - bin_edges[0]
            axes[1].bar(bin_centers, pulls, 
                        width=bin_width*0.8,  # 80% of bin width for spacing
                        color='skyblue', 
                        edgecolor='darkblue',
                        alpha=0.7, 
                        label='Pulls')
            
            # Calculate pull statistics
            mean_pull = np.mean(pulls)
            std_pull = np.std(pulls)
            
            pull_stats = (fr'Pull Statistics:'
                          fr'\\$\mu_{{pull}}$ = {mean_pull:.3f}'
                          fr'\\$\sigma_{{pull}}$ = {std_pull:.3f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            
            axes[1].text(0.02, 0.98, pull_stats, transform=axes[1].transAxes,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8), fontsize=9)
            
            axes[1].set_ylabel('(data - fit) / sqrt(data) at bin', fontweight='bold')
            axes[1].legend(loc='upper right')
            
            axes[1].set_xlabel(f'{x_label}', fontweight='bold')
            axes[1].grid(True, alpha=0.5, linestyle='--')
            
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
        
        # Return additional information
        fit_results = {
            'params': fit_params,
            'errors': param_errors,
            'cov_matrix': cov_matrix,
            'chi2': chi2,
            'reduced_chi2': reduced_chi2,
            'effective_sigma': effective_sigma,
            'effective_sigma_err': effective_sigma_err,
            'effective_mean': effective_mean,
            'effective_mean_err': effective_mean_err,
            'result': result
        }
        
        return fit_results, fig, axes
    
    fit_results = {
        'params': fit_params,
        'errors': param_errors,
        'cov_matrix': cov_matrix,
        'chi2': chi2,
        'reduced_chi2': reduced_chi2,
        'effective_sigma': effective_sigma,
        'effective_sigma_err': effective_sigma_err,
        'effective_mean': effective_mean,
        'effective_mean_err': effective_mean_err,
        'result': result
    }
    
    return fit_results