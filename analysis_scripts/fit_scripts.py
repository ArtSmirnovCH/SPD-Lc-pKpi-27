import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from typing import Tuple, Union


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


###########################################################################################################
# fit_distr_double_gauss
###########################################################################################################

def fit_distr_double_gauss(distr: Union[pd.Series, np.ndarray],
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
    
    
    def calculate_effective_sigma_double_gaussian(fit_params, fit_errors=None):
        """
        Calculate effective sigma for double Gaussian mixture
        
        For a mixture distribution: f(x) = w1*g1(x) + w2*g2(x)
        where w_i = A_i / (A1 + A2) are the weights
        
        The variance of the mixture is:
        sigma_eff = sqrt(w1 * sigma1**2 + w2*sigma2**2)
        
        Optionally propagates errors if fit_errors is provided
        """
        A1, mu1, sigma1, A2, mu2, sigma2, C = fit_params
        
        # Calculate weights
        total_amplitude = A1 + A2
        w1 = A1 / total_amplitude
        w2 = A2 / total_amplitude
        
        # Calculate effective mean
        effective_mean = w1 * mu1 + w2 * mu2
        
        # Calculate effective variance
        effective_variance = w1 * sigma1**2 + w2 * sigma2**2
        effective_sigma = np.sqrt(effective_variance)
        
        # Error propagation if errors are provided
        effective_mean_err = None
        effective_sigma_err = None
        
        if fit_errors is not None:
            A1_err, mu1_err, sigma1_err, A2_err, mu2_err, sigma2_err, C_err = fit_errors
            
            # Partial derivatives for effective mean
            dmean_dA1 = (mu1 - mu2) * A2 / (total_amplitude**2)
            dmean_dA2 = (mu2 - mu1) * A1 / (total_amplitude**2)
            dmean_dmu1 = w1
            dmean_dmu2 = w2
            
            # Error propagation for effective mean
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
            
            # Error propagation for effective variance
            sigma2_err_sq = (
                (df_sigma2_dA1 * A1_err)**2 + 
                (df_sigma2_dA2 * A2_err)**2 + 
                (df_sigma2_dsigma1 * sigma1_err)**2 + 
                (df_sigma2_dsigma2 * sigma2_err)**2
            )
            
            # Error propagation for effective sigma: d(sqrt(f))/dx = 0.5 * (1/sqrt(f)) * df/dx
            if effective_variance > 0:
                effective_sigma_err = 0.5 * np.sqrt(sigma2_err_sq) / np.sqrt(effective_variance)
        
        return effective_sigma, effective_mean, effective_sigma_err, effective_mean_err
    
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
            predictions = double_gaussian(x, *result.x)
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
    print(f"Constant: C = {fit_params[6]:.4f} +- {param_errors[6]:.4f}")
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
                label='Double Gaussian Fit'
            )
            
            # Plot individual Gaussians
            g1 = fit_params[0] * np.exp(-0.5 * ((x_fit - fit_params[1]) / fit_params[2]) ** 2)
            g2 = fit_params[3] * np.exp(-0.5 * ((x_fit - fit_params[4]) / fit_params[5]) ** 2)
            C = fit_params[6]
            
            # Format parameter strings with errors
            if not np.isnan(param_errors[1]):
                g1_label = fr'Gaussian 1: $\mu$={fit_params[1] * 1e4:.2f} $\pm$ {param_errors[1] * 1e4:.2f} $\mu m$'
                g2_label = fr'Gaussian 2: $\mu$={fit_params[4] * 1e4:.2f} $\pm$ {param_errors[4] * 1e4:.2f} $\mu m$'
                const_label = fr'Const: C={fit_params[6]:.2f} $\pm$ {param_errors[6]:.2f}'
            else:
                g1_label = fr'Gaussian 1: $\mu$={fit_params[1] * 1e4:.2f} $\mu m$'
                g2_label = fr'Gaussian 2: $\mu$={fit_params[4] * 1e4:.2f} $\mu m$'
                const_label = fr'Const: C={fit_params[6]:.2f}'
            
            axes[0].plot(x_fit, g1, 'g:', label=g1_label, alpha=0.7)
            axes[0].plot(x_fit, g2, 'b:', label=g2_label, alpha=0.7)
            axes[0].plot(x_fit, C * np.ones(shape=g1.size), 'm:', label=const_label, alpha=0.7)
            
            # Add fit parameters as text with errors
            if not np.isnan(param_errors[0]):
                fit_text = (fr'Fit Parameters:'
                          fr'\\G1: A={fit_params[0]:.1f} $\pm$ {param_errors[0]:.1f},\\ $\mu$={fit_params[1] * 1e4:.3f} $\pm$ {param_errors[1] * 1e4:.3f} $\mu m$,\\ $\sigma$={fit_params[2] * 1e4:.3f} $\pm$ {param_errors[2] * 1e4:.3f} $\mu m$'
                          fr'\\G2: A={fit_params[3]:.1f} $\pm$ {param_errors[3]:.1f},\\ $\mu$={fit_params[4] * 1e4:.3f} $\pm$ {param_errors[4] * 1e4:.3f} $\mu m$,\\ $\sigma$={fit_params[5] * 1e4:.3f} $\pm$ {param_errors[5] * 1e4:.3f} $\mu m$'
                          fr'\\Const: C={fit_params[6]:.1f} $\pm$ {param_errors[6]:.1f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            else:
                fit_text = (fr'Fit Parameters:'
                          fr'\\G1: A={fit_params[0]:.1f},\\ $\mu$={fit_params[1] * 1e4:.3f} $\mu m$,\\ $\sigma$={fit_params[2] * 1e4:.3f} $\mu m$'
                          fr'\\G2: A={fit_params[3]:.1f},\\ $\mu$={fit_params[4] * 1e4:.3f} $\mu m$,\\ $\sigma$={fit_params[5] * 1e4:.3f} $\mu m$'
                          fr'\\Const: C={fit_params[6]:.1f}'
                          fr'\\$\chi^2$/ndf = {chi2:.1f}/{n_dof}')
            
            if effective_sigma_err is not None and not np.isnan(effective_sigma_err):
                fit_text += fr'\\Effective $\sigma$: {effective_sigma * 1e4:.4f} $\pm$ {effective_sigma_err * 1e4:.4f} $\mu m$'
            else:
                fit_text += fr'\\Effective $\sigma$: {effective_sigma * 1e4:.4f} $\mu m$'
            
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
# fit_distr_flat_bg
###########################################################################################################

def fit_distr_flat_bg(distr: Union[pd.Series, np.ndarray],
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
    
    # Better Initial Guess
    initial_amplitude = np.max(counts) * (bin_edges[1] - bin_edges[0])  # Scale by bin width
    
    # Use percentiles to ensure means are within bin range
    data_min = np.min(bin_centers)
    data_max = np.max(bin_centers)
    data_range = data_max - data_min
    
    # Spread initial means across the data range
    initial_guess = [initial_amplitude]   # Constant
    
    # Set bounds for parameters
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