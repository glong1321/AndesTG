"""
River Profile Inversion Analysis
Python implementation of Goren et al. (2014) linear inversion methods
for inferring uplift rate history from river profiles.

Author: Adapted for TopoAnalysis framework
Date: November 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.optimize import curve_fit
import pandas as pd
from typing import Tuple, Dict, Optional, List


# ============================================================================
# PART 1: CORE CHI COMPUTATION AND NETWORK EXTRACTION
# ============================================================================

def compute_receiver_array(fd_grid):
    """
    Convert D8 flow direction grid to receiver index array.
    
    Parameters
    ----------
    fd_grid : np.ndarray
        Flow direction grid using D8 encoding (1,2,4,8,16,32,64,128)
        
    Returns
    -------
    rec_array : np.ndarray
        Flattened receiver indices for each pixel
    """
    nrows, ncols = fd_grid.shape
    
    # D8 direction to offset mapping
    d8_offsets = {
        1: (0, 1),    # East
        2: (1, 1),    # SE
        4: (1, 0),    # South
        8: (1, -1),   # SW
        16: (0, -1),  # West
        32: (-1, -1), # NW
        64: (-1, 0),  # North
        128: (-1, 1)  # NE
    }
    
    rec_array = np.arange(nrows * ncols)  # Initialize as self-receiver
    
    for idx in range(nrows * ncols):
        row, col = idx // ncols, idx % ncols
        flow_dir = fd_grid[row, col]
        
        if flow_dir in d8_offsets:
            dr, dc = d8_offsets[flow_dir]
            rec_row, rec_col = row + dr, col + dc
            
            # Check bounds
            if 0 <= rec_row < nrows and 0 <= rec_col < ncols:
                rec_array[idx] = rec_row * ncols + rec_col
    
    return rec_array


def extract_tributary(area, fd, elevation, outlet_xy, min_area=1e6, m=0.45, A0=1e7):
    """
    Extract tributary network given outlet coordinates using TopoAnalysis Chi class.
    
    Parameters
    ----------
    area : TopoAnalysis.Area
        Drainage area object
    fd : TopoAnalysis.FlowDirectionD8
        Flow direction object
    elevation : TopoAnalysis.Elevation
        Elevation object
    outlet_xy : tuple
        (x, y) UTM coordinates of outlet
    min_area : float
        Minimum drainage area threshold (m²)
    m : float
        Area exponent for chi calculation
    A0 : float
        Reference drainage area (m²)
        
    Returns
    -------
    dict with keys:
        - indices: flat indices of tributary pixels
        - x, y: coordinates
        - z: elevations
        - area: drainage areas
        - chi: chi values
        - rec_array: receiver indices for tributary only
    """
    # Import TopoAnalysis dem module
    try:
        from TopoAnalysis import dem as d
    except ImportError:
        raise ImportError("TopoAnalysis package not found. Make sure it's installed.")
    
    # Create Chi object for this outlet
    chi_obj = d.Chi(flow_direction=fd, area=area, theta=m, Ao=A0, outlets=[outlet_xy])
    
    # Extract data where chi > 0 (the tributary network)
    chi_mask = chi_obj._griddata > 0
    chi_data = chi_obj._griddata[chi_mask]
    elev_data = elevation._griddata[chi_mask]
    area_data = area._griddata[chi_mask]
    
    # Apply area threshold
    area_mask = area_data >= min_area
    
    chi_filtered = chi_data[area_mask]
    elev_filtered = elev_data[area_mask]
    area_filtered = area_data[area_mask]
    
    # Get indices and coordinates
    nrows, ncols = area._griddata.shape
    all_indices = np.where(chi_mask)
    rows_all = all_indices[0]
    cols_all = all_indices[1]
    
    # Filter by area threshold
    rows = rows_all[area_mask]
    cols = cols_all[area_mask]
    tributary_indices = rows * ncols + cols
    
    # Get coordinates
    georef = area._georef_info
    x_coords = georef.left + cols * georef.dx + georef.dx / 2
    y_coords = georef.top - rows * georef.dy - georef.dy / 2
    
    # Handle NaN values
    valid = ~np.isnan(elev_filtered)
    tributary_indices = tributary_indices[valid]
    x_coords = x_coords[valid]
    y_coords = y_coords[valid]
    chi_filtered = chi_filtered[valid]
    elev_filtered = elev_filtered[valid]
    area_filtered = area_filtered[valid]
    
    if len(tributary_indices) == 0:
        raise ValueError(f"No valid pixels found. Check outlet coordinates {outlet_xy} and min_area threshold.")
    
    # Compute receiver array for tributary
    full_rec_array = compute_receiver_array(fd._griddata)
    
    # Map global indices to local tributary indices
    index_map = {global_idx: local_idx for local_idx, global_idx in enumerate(tributary_indices)}
    local_rec_array = np.arange(len(tributary_indices))
    
    for i, global_idx in enumerate(tributary_indices):
        global_rec = full_rec_array[global_idx]
        if global_rec in index_map:
            local_rec_array[i] = index_map[global_rec]
        else:
            local_rec_array[i] = i  # Self-receiver if receiver not in tributary
    
    return {
        'indices': tributary_indices,
        'x': x_coords,
        'y': y_coords,
        'z': elev_filtered,
        'area': area_filtered,
        'chi': chi_filtered,
        'rec_array': local_rec_array
    }


def compute_chi(x, y, rec_array, area_array, m=0.45, A0=1e7):
    """
    Calculate chi parameter for river pixels.
    
    Parameters
    ----------
    x, y : np.ndarray
        Coordinates of each pixel (m)
    rec_array : np.ndarray
        Receiver indices for each pixel
    area_array : np.ndarray
        Upstream drainage area for each pixel (m²)
    m : float
        Area exponent (assuming n=1)
    A0 : float
        Reference drainage area (m²), default 10 km²
        
    Returns
    -------
    chi : np.ndarray
        Chi values for each pixel (m)
    """
    n = len(x)
    
    # Sort by drainage area (largest first) to process downstream to upstream
    data_mat = np.column_stack([np.arange(n), x, y, area_array])
    sorted_mat = data_mat[np.argsort(-data_mat[:, 3])]
    
    chi = np.zeros(n)
    
    for i in range(n):
        my_id = int(sorted_mat[i, 0])
        my_rec = rec_array[my_id]
        
        if my_id == my_rec:  # Outlet
            chi[i] = 0
        else:
            # Find receiver in sorted array
            rec_idx = np.where(sorted_mat[:, 0] == my_rec)[0][0]
            
            # Calculate distance
            dist = np.sqrt((sorted_mat[i, 1] - sorted_mat[rec_idx, 1])**2 +
                          (sorted_mat[i, 2] - sorted_mat[rec_idx, 2])**2)
            
            # Chi integral
            chi[i] = chi[rec_idx] + dist * (A0 / sorted_mat[i, 3])**m
    
    # Resort back to original order
    sorted_mat = np.column_stack([sorted_mat, chi])
    desorted_mat = sorted_mat[np.argsort(sorted_mat[:, 0])]
    
    return desorted_mat[:, -1]


# ============================================================================
# PART 2: OPTIMIZATION AND DIAGNOSTICS
# ============================================================================

def find_m_collapse_chi(x, y, z, rec_array, area_array, m_range=None, n_bins=25):
    """
    Find optimal m that minimizes scatter in chi-z domain.
    
    Parameters
    ----------
    x, y, z : np.ndarray
        Pixel coordinates and elevations
    rec_array : np.ndarray
        Receiver indices
    area_array : np.ndarray
        Drainage areas
    m_range : tuple, optional
        (min_m, max_m, step), default (0.1, 0.95, 0.05)
    n_bins : int
        Number of bins for scatter calculation
        
    Returns
    -------
    m_opt : float
        Optimal m value
    chi_opt : np.ndarray
        Chi values using optimal m
    scatter_metric : np.ndarray
        Scatter metric for each m value tested
    """
    if m_range is None:
        m_range = (0.1, 0.95, 0.05)
    
    m_try = np.arange(m_range[0], m_range[1], m_range[2])
    scatter_metric = np.zeros(len(m_try))
    A0 = 1e7  # 10 km²
    
    for i, m in enumerate(m_try):
        chi = compute_chi(x, y, rec_array, area_array, m=m, A0=A0)
        
        # Remove zeros
        valid = z != 0
        valid_chi = chi[valid]
        valid_z = z[valid]
        
        # Calculate scatter in bins
        val_in_bins = len(valid_chi) // n_bins
        chi_scatter = np.zeros(n_bins)
        
        # Sort by chi
        sorted_idx = np.argsort(valid_chi)
        sorted_chi = valid_chi[sorted_idx]
        sorted_z = valid_z[sorted_idx]
        
        for k in range(n_bins - 1):
            start_idx = val_in_bins * k
            end_idx = val_in_bins * (k + 1)
            chi_scatter[k] = np.std(sorted_z[start_idx:end_idx])
        
        chi_scatter[-1] = np.std(sorted_z[val_in_bins * (n_bins - 1):])
        scatter_metric[i] = np.mean(chi_scatter)
    
    # Find minimum
    min_idx = np.argmin(scatter_metric)
    m_opt = m_try[min_idx]
    chi_opt = compute_chi(x, y, rec_array, area_array, m=m_opt, A0=A0)
    
    return m_opt, chi_opt, scatter_metric, m_try


def find_m_slope_area(slope_array, area_array):
    """
    Calculate concavity index from slope-area relationship.
    
    Parameters
    ----------
    slope_array : np.ndarray
        Slope values (dimensionless)
    area_array : np.ndarray
        Drainage area values (m²)
        
    Returns
    -------
    m : float
        Best-fit concavity index
    lower_bound, upper_bound : float
        95% confidence interval bounds
    R2 : float
        R-squared value
    """
    # Use only positive slopes
    valid = slope_array > 0
    log_area = np.log(area_array[valid])
    log_slope = np.log(slope_array[valid])
    
    # Linear regression
    coeffs = np.polyfit(log_area, log_slope, 1)
    p = np.poly1d(coeffs)
    
    # R-squared
    yhat = p(log_area)
    ybar = np.mean(log_slope)
    ssreg = np.sum((yhat - ybar)**2)
    sstot = np.sum((log_slope - ybar)**2)
    R2 = ssreg / sstot if sstot > 0 else 0
    
    m = -coeffs[0]
    
    # Simple confidence interval estimate (approximate)
    residuals = log_slope - yhat
    std_err = np.std(residuals) / np.sqrt(len(log_area))
    margin = 1.96 * std_err  # 95% CI
    
    return m, m - margin, m + margin, R2


# ============================================================================
# PART 3: BLOCK UPLIFT INVERSION
# ============================================================================

def invert_block_uplift(chi, z, Gamma, q, plot=False):
    """
    Block uplift linear inversion following Goren et al. (2014).
    
    Parameters
    ----------
    chi : np.ndarray
        Chi values (m)
    z : np.ndarray
        Elevation values (m)
    Gamma : float
        Damping coefficient
    q : int
        Number of time intervals
    plot : bool
        Whether to plot results
        
    Returns
    -------
    Ustar : np.ndarray
        Non-dimensional uplift rate history (length q)
    tstar : np.ndarray
        Scaled time boundaries (length q+1)
    misfit : float
        RMS misfit between data and model
    """
    # Remove zeros
    valid = z != 0
    nz_chi = chi[valid]
    nz_z = z[valid]
    N = len(nz_chi)
    
    # Sort by chi
    sort_idx = np.argsort(nz_chi)
    sorted_chi = nz_chi[sort_idx]
    sorted_z = nz_z[sort_idx]
    
    # Construct time intervals (equal number of pixels per interval)
    val_per_dt = N // q
    scaled_t_vec = np.zeros(q + 1)
    scaled_t_vec[0] = 0
    
    for i in range(1, q):
        scaled_t_vec[i] = sorted_chi[(i - 1) * val_per_dt]
    scaled_t_vec[q] = sorted_chi[-1]
    
    scaled_dt_vec = np.diff(scaled_t_vec)
    
    # Construct forward model matrix A*
    Astar = np.zeros((N, q))
    
    for i in range(N):
        filled_full = np.where(scaled_t_vec >= sorted_chi[i])[0][0] - 1
        
        if filled_full > 0:
            Astar[i, :filled_full] = scaled_dt_vec[:filled_full]
        
        if filled_full < q:
            Astar[i, filled_full] = sorted_chi[i] - np.sum(scaled_dt_vec[:filled_full])
    
    # Prior model
    U_pri_star = np.ones(q) * np.mean(sorted_z / sorted_chi)
    
    # Least squares with Tikhonov regularization
    ATA = Astar.T @ Astar
    damping = Gamma**2 * np.eye(q)
    denom = ATA + damping
    nom = Astar.T @ (sorted_z - Astar @ U_pri_star)
    
    Ustar = U_pri_star + np.linalg.solve(denom, nom)
    tstar = scaled_t_vec
    
    # Misfit
    residuals = sorted_z - Astar @ Ustar
    misfit = np.sqrt(np.sum(residuals**2) / (N - q))
    
    if plot:
        plot_uplift_history(tstar, Ustar, xlabel='t* [m]', ylabel='U*')
    
    return Ustar, tstar, misfit


def bootstrap_invert_block_uplift(chi, z, Gamma, q, percent_sample=0.8, 
                                   num_iterations=100, K=1, plot=False):
    """
    Bootstrap block uplift inversion for uncertainty estimation.
    
    Parameters
    ----------
    chi, z : np.ndarray
        Chi and elevation values
    Gamma : float
        Damping coefficient
    q : int
        Number of time intervals
    percent_sample : float
        Fraction of data to sample (0-1)
    num_iterations : int
        Number of bootstrap iterations
    K : float
        Erodibility coefficient for dimensional plotting (default 1 for non-dimensional)
    plot : bool
        Whether to plot results
        
    Returns
    -------
    Ustar_mat : np.ndarray
        Matrix of uplift histories (num_iterations x q)
    tstar_best : np.ndarray
        Time boundaries from full dataset inversion
    Ustar_best : np.ndarray
        Best-fit uplift history from full dataset
    """
    # Remove zeros
    valid = z != 0
    nz_chi = chi[valid]
    nz_z = z[valid]
    
    data_length = len(nz_chi)
    sample_length = int(data_length * percent_sample)
    
    Ustar_mat = np.zeros((num_iterations, q))
    tstar_mat = np.zeros((num_iterations, q + 1))
    
    # Bootstrap iterations
    for i in range(num_iterations):
        sample_idx = np.sort(np.random.choice(data_length, sample_length, replace=False))
        Ustar, tstar, _ = invert_block_uplift(nz_chi[sample_idx], nz_z[sample_idx], 
                                               Gamma, q, plot=False)
        Ustar_mat[i, :] = Ustar
        tstar_mat[i, :] = tstar
    
    # Full dataset inversion
    Ustar_best, tstar_best, _ = invert_block_uplift(nz_chi, nz_z, Gamma, q, plot=False)
    
    if plot:
        plot_bootstrap_results(tstar_mat, Ustar_mat, tstar_best, Ustar_best, K)
    
    return Ustar_mat, tstar_best, Ustar_best


def invert_with_different_gamma(chi, z, q, gamma_range=None):
    """
    Generate L-curve for selecting optimal damping coefficient.
    
    Parameters
    ----------
    chi, z : np.ndarray
        Chi and elevation values
    q : int
        Number of time intervals
    gamma_range : tuple, optional
        (min, max, num_points), default (0.1, 100, 100)
        
    Returns
    -------
    Gamma_vec : np.ndarray
        Gamma values tested
    Misfit_vec : np.ndarray
        Corresponding misfit values
    """
    if gamma_range is None:
        Gamma_vec = np.logspace(-1, 2, 100)
    else:
        Gamma_vec = np.logspace(np.log10(gamma_range[0]), 
                                np.log10(gamma_range[1]), 
                                gamma_range[2])
    
    Misfit_vec = np.zeros(len(Gamma_vec))
    
    for i, Gamma in enumerate(Gamma_vec):
        _, _, Misfit_vec[i] = invert_block_uplift(chi, z, Gamma, q, plot=False)
    
    return Gamma_vec, Misfit_vec


# ============================================================================
# PART 4: SPACE-TIME INVERSION (PARABOLIC UPLIFT)
# ============================================================================

def invert_parabola(chi, z, x, rec_array, Gamma, q, K=1, plot=False):
    """
    Space-time inversion with parabolic spatial uplift pattern.
    
    Parameters
    ----------
    chi, z, x : np.ndarray
        Chi, elevation, and x-coordinate values
    rec_array : np.ndarray
        Receiver indices
    Gamma : float
        Damping coefficient
    q : int
        Number of time intervals
    K : float
        Erodibility coefficient (default 1 for non-dimensional)
    plot : bool
        Whether to plot results
        
    Returns
    -------
    Up : np.ndarray
        Parabolic uplift parameters (3q values: a, b, c for each time interval)
    tstar : np.ndarray
        Scaled time boundaries
    misfit : float
        RMS misfit
    """
    # Remove zeros
    valid = z != 0
    nz_chi = chi[valid]
    nz_z = z[valid]
    nz_x = x[valid]
    nz_id = np.arange(len(chi))[valid]
    
    N = len(nz_chi)
    
    # Sort by chi
    sort_data = np.column_stack([nz_id, nz_chi, nz_z])
    sort_idx = np.argsort(sort_data[:, 1])
    sorted_data = sort_data[sort_idx]
    sorted_chi = sorted_data[:, 1]
    
    # Time intervals
    val_per_dt = N // q
    scaled_t_vec = np.zeros(q + 1)
    scaled_t_vec[0] = 0
    
    for i in range(1, q):
        scaled_t_vec[i] = sorted_chi[(i - 1) * val_per_dt]
    scaled_t_vec[q] = sorted_chi[-1]
    
    scaled_dt_vec = np.diff(scaled_t_vec)
    
    # Build space-time matrix A
    A = np.zeros((N, N * q))
    
    for i in range(N):
        global_time = nz_chi[i]
        time_int = 0
        local_time = scaled_dt_vec[time_int]
        curr_id = nz_id[i]
        curr_ind = i
        id_rec = rec_array[curr_id]
        next_ind = np.where(nz_id == id_rec)[0]
        next_ind = next_ind[0] if len(next_ind) > 0 else None
        
        while global_time > 1e-12 and time_int < q:
            col = time_int * N + curr_ind
            
            if next_ind is None:
                chi_next = 0
            else:
                chi_next = nz_chi[next_ind]
            
            if global_time - chi_next >= local_time + 1e-10:
                A[i, col] = local_time
                global_time -= local_time
                time_int += 1
                if time_int < len(scaled_dt_vec):
                    local_time = scaled_dt_vec[time_int]
                else:
                    break
            else:
                A[i, col] = global_time - chi_next
                local_time -= (global_time - chi_next)
                global_time = chi_next
                curr_ind = next_ind
                curr_id = id_rec
                id_rec = rec_array[curr_id] if curr_id < len(rec_array) else curr_id
                next_ind = np.where(nz_id == id_rec)[0]
                next_ind = next_ind[0] if len(next_ind) > 0 else None
    
    # Build parabola matrix Bp
    Bp = np.zeros((q * N, 3 * q))
    x_km = nz_x / 1e3  # Convert to km
    
    for i in range(q):
        for j in range(N):
            Bp[i * N + j, i * 3] = x_km[j]**2
            Bp[i * N + j, i * 3 + 1] = x_km[j]
            Bp[i * N + j, i * 3 + 2] = 1
    
    Ap = A @ Bp
    
    # Prior model
    U_pri_star = np.ones(3 * q) * np.mean(nz_z / nz_chi)
    
    # Least squares
    ApTAp = Ap.T @ Ap
    damping = Gamma**2 * np.eye(3 * q)
    denom = ApTAp + damping
    nom = Ap.T @ (nz_z - Ap @ U_pri_star)
    
    Up = U_pri_star + np.linalg.solve(denom, nom)
    tstar = scaled_t_vec
    
    # Misfit
    residuals = nz_z - Ap @ Up
    misfit = np.sqrt(np.sum(residuals**2) / (N - 3 * q))
    
    if plot:
        plot_spacetime_uplift(Up, tstar, x, K)
    
    return Up, tstar, misfit


# ============================================================================
# PART 5: CALIBRATION
# ============================================================================

def calibrate_k_total_uplift(H, t_H, Ustar, tstar, A0=1e7, m=0.45, plot=False):
    """
    Calibrate erodibility coefficient K from total uplift data.
    
    Parameters
    ----------
    H : float
        Total uplift (m) from present to time t_H
    t_H : float
        Age of dated uplifted feature (years)
    Ustar : np.ndarray
        Non-dimensional uplift rate history
    tstar : np.ndarray
        Scaled time boundaries
    A0 : float
        Reference drainage area (m²)
    m : float
        Area exponent
    plot : bool
        Whether to plot dimensional uplift history
        
    Returns
    -------
    K : float
        Calibrated erodibility coefficient (m^(1-2m)/yr)
    U : np.ndarray
        Dimensional uplift rate history (m/yr)
    t : np.ndarray
        Dimensional time boundaries (years)
    """
    del_t_star = np.diff(tstar)
    testH_vec = np.cumsum(Ustar * del_t_star)
    
    # Find time interval containing H
    time_ind = np.where(testH_vec > H)[0]
    if len(time_ind) == 0:
        time_ind = len(testH_vec) - 1
    else:
        time_ind = time_ind[0] - 1
    
    rem_H = H - testH_vec[time_ind]
    del_scaled_time = rem_H / Ustar[time_ind + 1]
    t_H_star = tstar[time_ind + 1] + del_scaled_time
    
    K = t_H_star / (t_H * A0**m)
    U = Ustar * K
    t = tstar / K
    
    if plot:
        plot_uplift_history(t, U, xlabel='t [yr]', ylabel='U [m/yr]')
    
    return K, U, t


# ============================================================================
# PART 6: PLOTTING UTILITIES
# ============================================================================

def plot_chi_z(chi, z, rec_array, ax=None):
    """Plot chi-z relationship."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot river segments
    for i in range(len(chi)):
        j = rec_array[i]
        if j != i:
            ax.plot([chi[j], chi[i]], [z[j], z[i]], 'b-', alpha=0.5, linewidth=0.5)
    
    ax.set_xlabel('χ [m]', fontsize=12)
    ax.set_ylabel('Elevation [m]', fontsize=12)
    ax.set_title('χ-z Profile')
    ax.grid(True, alpha=0.3)
    
    return ax


def plot_slope_area(slope, area, m=None, ax=None):
    """Plot slope-area relationship."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    valid = slope > 0
    ax.loglog(area[valid], slope[valid], 'o', alpha=0.5, markersize=3)
    
    if m is not None:
        # Add best-fit line
        log_area = np.log10(area[valid])
        log_area_line = np.linspace(log_area.min(), log_area.max(), 100)
        area_line = 10**log_area_line
        slope_line = 10**(log_area_line.max()) * (area_line / area_line.max())**(-m)
        ax.loglog(area_line, slope_line, 'r-', linewidth=2, 
                 label=f'm = {m:.3f}')
        ax.legend()
    
    ax.set_xlabel('Drainage Area [m²]', fontsize=12)
    ax.set_ylabel('Slope', fontsize=12)
    ax.set_title('Slope-Area Relationship')
    ax.grid(True, alpha=0.3)
    
    return ax


def plot_uplift_history(tstar, Ustar, xlabel='t*', ylabel='U*', ax=None):
    """Plot uplift rate history as staircase."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    q = len(Ustar)
    t_plot = []
    U_plot = []
    
    for i in range(q):
        t_plot.extend([tstar[i], tstar[i+1]])
        U_plot.extend([Ustar[i], Ustar[i]])
    
    ax.plot(t_plot, U_plot, 'b-', linewidth=2)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title('Uplift Rate History')
    ax.grid(True, alpha=0.3)
    
    return ax


def plot_bootstrap_results(tstar_mat, Ustar_mat, tstar_best, Ustar_best, K=1):
    """Plot bootstrap inversion results with uncertainty."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    num_iterations, q = Ustar_mat.shape
    
    # Plot all bootstrap realizations in light gray
    for i in range(num_iterations):
        t_plot = []
        U_plot = []
        for j in range(q):
            t_plot.extend([tstar_mat[i, j], tstar_mat[i, j+1]])
            U_plot.extend([Ustar_mat[i, j], Ustar_mat[i, j]])
        
        ax.plot(np.array(t_plot) / K / 1e6, 
               np.array(U_plot) * K / 1e-3,
               color=[0.9, 0.9, 0.9], linewidth=1)
    
    # Plot best-fit solution
    t_best = []
    U_best = []
    for j in range(q):
        t_best.extend([tstar_best[j], tstar_best[j+1]])
        U_best.extend([Ustar_best[j], Ustar_best[j]])
    
    ax.plot(np.array(t_best) / K / 1e6,
           np.array(U_best) * K / 1e-3,
           'k-', linewidth=2, label='Best fit')
    
    # Plot mean and std
    mean_U = np.mean(Ustar_mat, axis=0)
    std_U = np.std(Ustar_mat, axis=0)
    
    t_mean = []
    U_mean = []
    U_up = []
    U_down = []
    for j in range(q):
        t_mean.extend([tstar_best[j], tstar_best[j+1]])
        U_mean.extend([mean_U[j], mean_U[j]])
        U_up.extend([mean_U[j] + std_U[j], mean_U[j] + std_U[j]])
        U_down.extend([mean_U[j] - std_U[j], mean_U[j] - std_U[j]])
    
    ax.plot(np.array(t_mean) / K / 1e6,
           np.array(U_mean) * K / 1e-3,
           'm-', linewidth=1.5, label='Mean')
    ax.plot(np.array(t_mean) / K / 1e6,
           np.array(U_up) * K / 1e-3,
           'm:', linewidth=1, label='Mean ± 1σ')
    ax.plot(np.array(t_mean) / K / 1e6,
           np.array(U_down) * K / 1e-3,
           'm:', linewidth=1)
    
    ax.set_xlabel('t [Ma]', fontsize=14)
    ax.set_ylabel('U [mm/yr]', fontsize=14)
    ax.set_title('Bootstrap Block Uplift Inversion')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    return ax


def plot_spacetime_uplift(Up, tstar, x, K=1):
    """Plot space-time uplift pattern."""
    q = len(tstar) - 1
    x_domain = np.linspace(x.min(), x.max(), 100) / 1e3  # Convert to km
    
    U_space_time = np.zeros((q, len(x_domain)))
    
    for i in range(q):
        a = Up[i * 3]
        b = Up[i * 3 + 1]
        c = Up[i * 3 + 2]
        U_space_time[i, :] = a * x_domain**2 + b * x_domain + c
    
    # Create staircase in time
    t_plot = []
    U_plot = []
    for i in range(q):
        t_plot.extend([tstar[i], tstar[i+1]])
        U_plot.append(U_space_time[i, :])
        U_plot.append(U_space_time[i, :])
    
    U_plot = np.array(U_plot)
    t_plot = np.array(t_plot)
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    X, T = np.meshgrid(x_domain, t_plot)
    
    if K == 1:
        surf = ax.plot_surface(X, T, U_plot, cmap='viridis', 
                               edgecolor='none', alpha=0.9)
        ax.set_ylabel('Scaled t [m]')
        ax.set_zlabel('Non-dimensional Uplift Rate')
    else:
        surf = ax.plot_surface(X, T / K / 1e6, U_plot * K / 1e-3, 
                               cmap='viridis', edgecolor='none', alpha=0.9)
        ax.set_ylabel('Time [Ma]')
        ax.set_zlabel('Uplift Rate [mm/yr]')
    
    ax.set_xlabel('x [km]')
    ax.set_title('Space-Time Uplift Pattern')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    return ax


def plot_l_curve(Gamma_vec, Misfit_vec, ax=None):
    """Plot L-curve for damping parameter selection."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(1 / Gamma_vec, Misfit_vec, 'b-', linewidth=2)
    ax.set_xlabel('1/Γ', fontsize=12)
    ax.set_ylabel('Misfit [m]', fontsize=12)
    ax.set_title('L-Curve for Damping Parameter Selection')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([-1, np.max(1 / Gamma_vec) * 1.1])
    ax.set_ylim([0, np.max(Misfit_vec) * 1.1])
    
    return ax


def plot_m_optimization(m_try, scatter_metric, chi, z, rec_array, m_opt):
    """Plot m optimization results."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Scatter metric vs m
    ax1.plot(m_try, scatter_metric, 'ob', markersize=8)
    ax1.axvline(m_opt, color='r', linestyle='--', linewidth=2, 
                label=f'Optimal m = {m_opt:.3f}')
    ax1.set_xlabel('m', fontsize=12)
    ax1.set_ylabel('Mean of z scatter in bins', fontsize=12)
    ax1.set_title('Concavity Index Optimization')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Chi-z with optimal m
    for i in range(len(chi)):
        j = rec_array[i]
        if j != i:
            ax2.plot([chi[j], chi[i]], [z[j], z[i]], 'b-', 
                    alpha=0.5, linewidth=0.5)
    
    ax2.set_xlabel('χ [m]', fontsize=12)
    ax2.set_ylabel('z [m]', fontsize=12)
    ax2.set_title(f'χ-z relation using optimal m = {m_opt:.3f}')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


# ============================================================================
# PART 7: HIGH-LEVEL WRAPPER FUNCTION
# ============================================================================

def run_river_inversion(area, fd, elevation, outlet_coords, 
                        m=0.45, A0=1e7, min_area=1e6,
                        Gamma=10, q=10, 
                        optimize_m=False,
                        block_uplift=True,
                        space_time=False,
                        bootstrap=False,
                        bootstrap_iterations=100,
                        bootstrap_sample_pct=0.8,
                        calibration_data=None,
                        plot_results=True,
                        save_results=None):
    """
    Complete river profile inversion workflow.
    
    Parameters
    ----------
    area : TopoAnalysis.Area
        Drainage area object
    fd : TopoAnalysis.FlowDirectionD8
        Flow direction object
    elevation : TopoAnalysis.Elevation
        Elevation object
    outlet_coords : tuple
        (x, y) UTM coordinates of outlet
    m : float
        Area exponent (default 0.45)
    A0 : float
        Reference drainage area in m² (default 1e7 = 10 km²)
    min_area : float
        Minimum drainage area threshold (m²)
    Gamma : float
        Damping coefficient for inversion
    q : int
        Number of time intervals
    optimize_m : bool
        Whether to optimize m value
    block_uplift : bool
        Whether to run block uplift inversion
    space_time : bool
        Whether to run space-time inversion
    bootstrap : bool
        Whether to run bootstrap analysis
    bootstrap_iterations : int
        Number of bootstrap iterations
    bootstrap_sample_pct : float
        Fraction of data to sample in bootstrap (0-1)
    calibration_data : dict, optional
        Dict with keys 'H' (total uplift, m), 't_H' (age, yr)
    plot_results : bool
        Whether to generate plots
    save_results : str, optional
        Path to save results (as .npz file)
        
    Returns
    -------
    results : dict
        Dictionary containing all results and intermediate data
    """
    print("="*60)
    print("RIVER PROFILE INVERSION ANALYSIS")
    print("="*60)
    
    results = {}
    
    # Step 1: Extract tributary network
    print("\n[1/7] Extracting tributary network...")
    tributary = extract_tributary(area, fd, elevation, outlet_coords, min_area, m=m, A0=A0)    results['tributary'] = tributary
    print(f"  - Found {len(tributary['x'])} channel pixels")
    print(f"  - Elevation range: {tributary['z'].min():.1f} - {tributary['z'].max():.1f} m")
    print(f"  - Area range: {tributary['area'].min()/1e6:.2f} - {tributary['area'].max()/1e6:.2f} km²")
    
    # Step 2: Optimize m if requested (skip if chi already computed with specified m)
    if optimize_m:
        print("\n[2/7] Optimizing concavity index m...")
        print("  - Note: Chi already computed with m={:.3f} during extraction".format(m))
        print("  - Re-optimizing m will require recomputing chi...")
        m_opt, chi_opt, scatter, m_try = find_m_collapse_chi(
            tributary['x'], tributary['y'], tributary['z'],
            tributary['rec_array'], tributary['area']
        )
        print(f"  - Optimal m = {m_opt:.3f}")
        chi = chi_opt
        m = m_opt
        results['m_optimization'] = {
            'm_opt': m_opt,
            'm_try': m_try,
            'scatter_metric': scatter
        }
    else:
        print(f"\n[2/7] Using specified m = {m:.3f}")
        # Use chi values already computed during extraction
        chi = tributary['chi']
    
    # Step 3: Store chi values
    print("\n[3/7] Chi parameter ready...")
    results['chi'] = chi
    results['m'] = m
    results['A0'] = A0
    print(f"  - χ range: {chi.min():.1f} - {chi.max():.1f} m")
    
    # Step 4: Block uplift inversion
    if block_uplift:
        print("\n[4/7] Running block uplift inversion...")
        Ustar, tstar, misfit = invert_block_uplift(
            chi, tributary['z'], Gamma, q, plot=False
        )
        results['block_uplift'] = {
            'Ustar': Ustar,
            'tstar': tstar,
            'misfit': misfit,
            'Gamma': Gamma,
            'q': q
        }
        print(f"  - RMS misfit: {misfit:.2f} m")
        print(f"  - Time range: {tstar[0]:.1f} - {tstar[-1]:.1f} m")
        
        # Bootstrap if requested
        if bootstrap:
            print(f"\n  Running bootstrap analysis ({bootstrap_iterations} iterations)...")
            Ustar_mat, tstar_best, Ustar_best = bootstrap_invert_block_uplift(
                chi, tributary['z'], Gamma, q,
                percent_sample=bootstrap_sample_pct,
                num_iterations=bootstrap_iterations,
                K=1, plot=False
            )
            results['block_uplift']['bootstrap'] = {
                'Ustar_mat': Ustar_mat,
                'mean_Ustar': np.mean(Ustar_mat, axis=0),
                'std_Ustar': np.std(Ustar_mat, axis=0)
            }
            print(f"  - Bootstrap complete")
    else:
        print("\n[4/7] Skipping block uplift inversion")
    
    # Step 5: Space-time inversion
    if space_time:
        print("\n[5/7] Running space-time inversion...")
        Up, tstar_st, misfit_st = invert_parabola(
            chi, tributary['z'], tributary['x'],
            tributary['rec_array'], Gamma, q, K=1, plot=False
        )
        results['space_time'] = {
            'Up': Up,
            'tstar': tstar_st,
            'misfit': misfit_st
        }
        print(f"  - RMS misfit: {misfit_st:.2f} m")
    else:
        print("\n[5/7] Skipping space-time inversion")
    
    # Step 6: Calibration
    if calibration_data is not None and block_uplift:
        print("\n[6/7] Calibrating with erosion/uplift data...")
        K, U, t = calibrate_k_total_uplift(
            calibration_data['H'],
            calibration_data['t_H'],
            results['block_uplift']['Ustar'],
            results['block_uplift']['tstar'],
            A0=A0, m=m, plot=False
        )
        results['calibration'] = {
            'K': K,
            'U': U,
            't': t,
            'H': calibration_data['H'],
            't_H': calibration_data['t_H']
        }
        print(f"  - Calibrated K = {K:.2e} m^(1-2m)/yr")
        print(f"  - Time range: {t[0]/1e6:.2f} - {t[-1]/1e6:.2f} Ma")
        print(f"  - Uplift rate range: {U.min()*1e3:.3f} - {U.max()*1e3:.3f} mm/yr")
    else:
        print("\n[6/7] No calibration data provided")
    
    # Step 7: Plotting
    if plot_results:
        print("\n[7/7] Generating plots...")
        
        # Plot 1: Chi-z profile
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        plot_chi_z(chi, tributary['z'], tributary['rec_array'], ax=ax1)
        plt.tight_layout()
        
        if optimize_m:
            # Plot 2: m optimization
            fig2 = plot_m_optimization(
                results['m_optimization']['m_try'],
                results['m_optimization']['scatter_metric'],
                chi, tributary['z'], tributary['rec_array'],
                results['m_optimization']['m_opt']
            )
        
        if block_uplift:
            # Plot 3: Block uplift results
            fig3, ax3 = plt.subplots(figsize=(10, 6))
            if bootstrap and calibration_data is not None:
                K_plot = results['calibration']['K']
                plot_bootstrap_results(
                    np.tile(results['block_uplift']['tstar'], 
                           (bootstrap_iterations, 1)),
                    results['block_uplift']['bootstrap']['Ustar_mat'],
                    results['block_uplift']['tstar'],
                    results['block_uplift']['Ustar'],
                    K=K_plot
                )
            elif calibration_data is not None:
                plot_uplift_history(
                    results['calibration']['t'],
                    results['calibration']['U'],
                    xlabel='Time [yr]',
                    ylabel='Uplift Rate [m/yr]',
                    ax=ax3
                )
            else:
                plot_uplift_history(
                    results['block_uplift']['tstar'],
                    results['block_uplift']['Ustar'],
                    xlabel='t* [m]',
                    ylabel='U*',
                    ax=ax3
                )
            plt.tight_layout()
        
        if space_time:
            # Plot 4: Space-time uplift
            K_plot = results['calibration']['K'] if calibration_data else 1
            fig4 = plt.figure(figsize=(12, 8))
            plot_spacetime_uplift(
                results['space_time']['Up'],
                results['space_time']['tstar'],
                tributary['x'],
                K=K_plot
            )
        
        plt.show()
        print("  - Plots generated")
    else:
        print("\n[7/7] Skipping plots")
    
    # Save results
    if save_results is not None:
        print(f"\nSaving results to {save_results}...")
        np.savez_compressed(save_results, **{k: v for k, v in results.items() 
                                              if not isinstance(v, dict)})
        print("  - Results saved")
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    
    return results


# ============================================================================
# PART 8: UTILITY FUNCTIONS
# ============================================================================

def compute_steady_state_uplift(chi, z):
    """
    Compute steady-state uplift pattern from dz/dχ.
    
    Parameters
    ----------
    chi, z : np.ndarray
        Chi and elevation values
        
    Returns
    -------
    dz_dchi : np.ndarray
        Derivative dz/dχ (proportional to U*)
    """
    # Sort by chi
    sort_idx = np.argsort(chi)
    chi_sorted = chi[sort_idx]
    z_sorted = z[sort_idx]
    
    # Compute derivative
    dz_dchi = np.gradient(z_sorted, chi_sorted)
    
    # Unsort
    unsort_idx = np.argsort(sort_idx)
    dz_dchi_unsorted = dz_dchi[unsort_idx]
    
    return dz_dchi_unsorted


def detect_knickpoints(chi, z, rec_array, threshold=0.5):
    """
    Detect knickpoints as breaks in slope (dz/dχ).
    
    Parameters
    ----------
    chi, z : np.ndarray
        Chi and elevation values
    rec_array : np.ndarray
        Receiver indices
    threshold : float
        Minimum relative change in slope to detect knickpoint
        
    Returns
    -------
    knickpoint_indices : np.ndarray
        Indices of detected knickpoints
    slope_changes : np.ndarray
        Relative slope change at each knickpoint
    """
    dz_dchi = compute_steady_state_uplift(chi, z)
    
    # Compute slope changes
    slope_changes = []
    knickpoint_idx = []
    
    for i in range(len(chi)):
        rec = rec_array[i]
        if rec != i:  # Not an outlet
            slope_here = dz_dchi[i]
            slope_downstream = dz_dchi[rec]
            
            if slope_downstream > 0:
                rel_change = abs(slope_here - slope_downstream) / slope_downstream
                
                if rel_change > threshold:
                    knickpoint_idx.append(i)
                    slope_changes.append(rel_change)
    
    return np.array(knickpoint_idx), np.array(slope_changes)


def estimate_knickpoint_timing(chi, tstar, knickpoint_indices):
    """
    Estimate formation timing of knickpoints using chi-time mapping.
    
    Parameters
    ----------
    chi : np.ndarray
        Chi values
    tstar : np.ndarray
        Time boundaries from inversion
    knickpoint_indices : np.ndarray
        Indices of knickpoints
        
    Returns
    -------
    knickpoint_times : np.ndarray
        Estimated formation times for each knickpoint
    """
    chi_knicks = chi[knickpoint_indices]
    
    # Map chi to time (linear interpolation)
    knickpoint_times = np.interp(chi_knicks, 
                                  [tstar[0], tstar[-1]], 
                                  [0, tstar[-1]])
    
    return knickpoint_times


def export_results_to_csv(results, output_path):
    """
    Export key results to CSV file.
    
    Parameters
    ----------
    results : dict
        Results dictionary from run_river_inversion
    output_path : str
        Path to save CSV file
    """
    if 'block_uplift' in results:
        df = pd.DataFrame({
            'time_interval_start': results['block_uplift']['tstar'][:-1],
            'time_interval_end': results['block_uplift']['tstar'][1:],
            'uplift_rate': results['block_uplift']['Ustar']
        })
        
        if 'calibration' in results:
            df['time_start_yr'] = results['calibration']['t'][:-1]
            df['time_end_yr'] = results['calibration']['t'][1:]
            df['uplift_rate_m_per_yr'] = results['calibration']['U']
        
        df.to_csv(output_path, index=False)
        print(f"Results exported to {output_path}")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

def example_workflow():
    """
    Example workflow showing how to use the river inversion code.
    """
    
    # Import your TopoAnalysis package
    from TopoAnalysis import dem as d
    
    # Load data
    area = d.Area.load('/path/to/rapel_area_utm30m')
    fd = d.FlowDirectionD8.load('/path/to/rapel_fd_utm30m')
    elevation = d.Elevation.load('/path/to/rapel_SRTMGL130m_dem_utm.tif')
    
    # Define outlet coordinates (example)
    outlet_coords = (300000, 6200000)  # Replace with your coordinates
    
    # Run basic inversion (no calibration)
    results = run_river_inversion(
        area=area,
        fd=fd,
        elevation=elevation,
        outlet_coords=outlet_coords,
        m=0.45,
        A0=1e7,
        min_area=1e6,
        Gamma=10,
        q=10,
        optimize_m=False,
        block_uplift=True,
        space_time=False,
        bootstrap=True,
        bootstrap_iterations=100,
        plot_results=True
    )
    
    # Run with calibration data
    calibration = {
        'H': 1000,  # 1000 m total uplift
        't_H': 5e6  # 5 Ma age
    }
    
    results_calibrated = run_river_inversion(
        area=area,
        fd=fd,
        elevation=elevation,
        outlet_coords=outlet_coords,
        m=0.45,
        Gamma=10,
        q=10,
        block_uplift=True,
        bootstrap=True,
        calibration_data=calibration,
        plot_results=True,
        save_results='river_inversion_results.npz'
    )
    
    # Export to CSV
    export_results_to_csv(results_calibrated, 'uplift_history.csv')
    
    return results_calibrated


if __name__ == "__main__":
    print(__doc__)
    print("\nThis module provides tools for river profile inversion analysis.")
    print("See example_workflow() for usage examples.")
