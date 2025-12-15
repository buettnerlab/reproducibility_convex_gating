import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, interpolate
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
import os


def _generate_biex_lut(channel_range=4096, pos=4.418540, neg=0.0, width_basis=-10, max_value=262144.000029):
    """Creates a FlowJo compatible biex lookup table.

    Code adopted from FlowKit Python package.
    Creates a FlowJo compatible biex lookup table.

    Implementation ported from the R library cytolib, which claims to be
    directly ported from the legacy Java code from TreeStar.

    Parameters
    ----------
    channel_range
        Maximum positive value of the output range.
    pos
        Number of decades.
    neg
        Number of extra negative decades.
    width_basis
        Controls the amount of input range compressed in the zero / linear region. A higher width basis value will include more input values in the zero / linear region.
    max_value
        Maximum input value to scale.

    Returns
    -------
    2-column NumPy array of the LUT (column order: input, output).
    """
    ln10 = np.log(10.0)
    decades = pos
    low_scale = width_basis
    width = np.log10(-low_scale)

    decades = decades - (width / 2)

    extra = neg

    if extra < 0:
        extra = 0

    extra = extra + (width / 2)

    zero_point = int((extra * channel_range) / (extra + decades))
    zero_point = int(np.min([zero_point, channel_range / 2]))

    if zero_point > 0:
        decades = extra * channel_range / zero_point

    width = width / (2 * decades)

    maximum = max_value
    positive_range = ln10 * decades
    minimum = maximum / np.exp(positive_range)

    negative_range = _log_root(positive_range, width)

    max_channel_value = channel_range + 1
    n_points = max_channel_value

    step = (max_channel_value - 1) / (n_points - 1)

    values = np.arange(n_points)
    positive = np.exp(values / float(n_points) * positive_range)
    negative = np.exp(values / float(n_points) * -negative_range)

    # apply step to values
    values = values * step

    s = np.exp((positive_range + negative_range) * (width + extra / decades))

    negative *= s
    s = positive[zero_point] - negative[zero_point]

    positive[zero_point:n_points] = positive[zero_point:n_points] - negative[zero_point:n_points]
    positive[zero_point:n_points] = minimum * (positive[zero_point:n_points] - s)

    neg_range = np.arange(zero_point)
    m = 2 * zero_point - neg_range

    positive[neg_range] = -positive[m]

    return positive, values

def _log_root(b: float, w: float) -> float:
    """Helper function.

    Parameters
    ----------
    b
        Upper bound
    w
        Step parameter

    Returns
    -------
    Solution to interpolation
    """
    # Code adopted from FlowKit Python package
    x_lo = 0.0
    x_hi = b
    d = (x_lo + x_hi) / 2
    dx = abs(int(x_lo - x_hi))  # type: float
    dx_last = dx
    fb = -2 * np.log(b) + w * b
    f = 2.0 * np.log(d) + w * b + fb
    df = 2 / d + w

    if w == 0:
        return b

    for _ in range(100):
        if (((d - x_hi) * df - f) - ((d - x_lo) * df - f)) > 0 or abs(2 * f) > abs(dx_last * df):
            dx = (x_hi - x_lo) / 2
            d = x_lo + dx
            if d == x_lo:
                return d
        else:
            dx = f / df
            t = d
            d -= dx
            # if dx is smaller than some precision threshold
            if d == t:
                return d

        # if dx is smaller than some precision threshold
        if abs(dx) < 1.0e-12:
            return d

        dx_last = dx
        f = 2 * np.log(d) + w * d + fb
        df = 2 / d + w
        if f < 0:
            x_lo = d
        else:
            x_hi = d

    return d



def generate_lut_func(channel_range = None, 
                      neg = None, 
                      width_basis = None, 
                      pos = None, 
                      max_value = None):
    """
    Generate biexponential transformation functions.
    
    Linear interpolation used for function values and constant outside the original range.
    
    Parameters
    ----------
    channel_range
        Maximum positive value of the output range.
    pos
        Number of decades.
    neg
        Number of extra negative decades.
    width_basis
        Controls the amount of input range compressed in the zero / linear region. A higher width basis value will include more input values in the zero / linear region.
    max_value
        Maximum input value to scale.

    Returns
    -------
    tuple
        (lut_func, inverse_lut_func) - biexponential transformation function and its inverse
    
    """
    # Define default values within the function
    if channel_range is None:
        channel_range = 4096  # default value
    if neg is None:
        neg = 0  # default value
    
    if width_basis is None:
        width_basis = -10  # default value
    
    if pos is None:
        pos = 4.418540  # default value
    
    if max_value is None:
        max_value = 262144  # default value
    
    x, y = _generate_biex_lut(
                         channel_range= channel_range,
                         neg=neg, 
                         width_basis=width_basis, 
                         pos=pos, 
                         max_value=max_value)
    
    lut_func = interpolate.interp1d(
            x, y, kind='linear', bounds_error=False, fill_value=(np.min(y), np.max(y))
        )
    inverse_lut_func = interpolate.interp1d(
            y, x, kind='linear', bounds_error=False, fill_value=(np.min(x), np.max(x))
        )
    
    return lut_func, inverse_lut_func

def generate_hybrid_ticks_and_labels(lims, transform_func=None, inverse_transform_func=None, 
                                    linear_threshold=100):
    """
    Generate hybrid ticks: minimal linear ticks (0, ±threshold), log for large values
    
    Parameters:
    -----------
    lims : array_like
        [min, max] limits in transformed space
    transform_func : callable, optional
        Function to transform from original to display space (e.g., lut_func, arcsinh)
        Default: lut_func, which refers to the biexp transformation used in pytometry and flowkit
    inverse_transform_func : callable, optional
        Function to transform from display to original space (e.g., inverse_lut_func, sinh)
        Default: inverse_lut_func, which refers to the inverse biexp transformation used in pytometry and flowkit
    linear_threshold : float, optional
        Threshold in original space below which to use minimal linear ticks
        Default: 100
        
    Returns:
    --------
    tuple
        (all_ticks, all_labels) - lists of tick positions and corresponding labels
        
    Notes:
    ------
    Example usage with different transforms:
    1. Biexponential (current):
       generate_hybrid_ticks_and_labels(lims, lut_func, inverse_lut_func)

    2. Arcsinh transformation:
       generate_hybrid_ticks_and_labels(lims, np.arcsinh, np.sinh)
    
    3. Log transformation (positive values only):
       generate_hybrid_ticks_and_labels(lims, np.log10, lambda x: 10**x)

    4. Custom transformation:
       def my_transform(x): return x**2
       def my_inverse(y): return np.sqrt(np.abs(y)) * np.sign(y)
       generate_hybrid_ticks_and_labels(lims, my_transform, my_inverse)
    """
    
    # Use default functions if none provided
    if transform_func is None:
        transform_func = lut_func
    if inverse_transform_func is None:
        inverse_transform_func = inverse_lut_func
    
    orig_lims = inverse_transform_func(lims)
    
    all_ticks = []
    all_labels = []
    
    # 1. Logarithmic ticks for large values
    power_lims = np.round(np.log10(np.abs(inverse_transform_func(lims))))
    power_lims = power_lims[~np.isnan(power_lims)]
    
    if len(power_lims) > 0:
        # Only add log ticks for values >= linear_threshold
        min_log_power = max(2, int(np.log10(linear_threshold)))  # Start from 10^2 = 100
        max_power = int(max(power_lims))
        
        for power in range(min_log_power, max_power + 1):
            # Positive powers
            val_pos = 10**power
            if orig_lims[0] <= val_pos <= orig_lims[1]:
                tick_pos = transform_func(val_pos)
                if lims[0] <= tick_pos <= lims[1]:
                    all_ticks.append(tick_pos)
                    all_labels.append(f"$10^{power}$")
            
            # Negative powers
            val_neg = -10**power
            if orig_lims[0] <= val_neg <= orig_lims[1]:
                tick_neg = transform_func(val_neg)
                if lims[0] <= tick_neg <= lims[1]:
                    all_ticks.append(tick_neg)
                    all_labels.append(f"$-10^{power}$")
    
    # 2. Minimal linear ticks: only 0 and ±threshold for reference
    reference_values = [0, linear_threshold, -linear_threshold]
    
    for val in reference_values:
        if orig_lims[0] <= val <= orig_lims[1]:
            tick_pos = transform_func(val)
            if lims[0] <= tick_pos <= lims[1]:
                all_ticks.append(tick_pos)
                if val == 0:
                    all_labels.append("0")
                elif val == linear_threshold:
                    all_labels.append("")
                elif val == -linear_threshold:
                    all_labels.append("")
    
    # Sort ticks
    if all_ticks:
        sorted_indices = np.argsort(all_ticks)
        all_ticks = [all_ticks[i] for i in sorted_indices]
        all_labels = [all_labels[i] for i in sorted_indices]
    
    return all_ticks, all_labels

def contour_lines_separate_outliers(
    meta_info, run_ID, base_df, population, population_path,
    num_hierarchies=1,
    pop_min=100,
    mode=None,  
    transform_func=None, 
    inverse_transform_func=None, 
    linear_threshold=100,
    outlier_percentile=1, 
    grid_size=50,
    target_color='#FF4F00', 
    non_target_color='#2d5380', 
    add_grid=False,
    save=None,
    save_png=None,
):
    """
    Create contour plots with outlier detection for flow cytometry gating analysis.
    
    This function generates publication-ready contour plots showing density distributions
    of target and non-target cell populations, with automatic outlier detection and 
    intelligent handling of small populations. Uses ultra-fast histogram-based density
    estimation for optimal performance with large datasets.
    
    Parameters:
    -----------
    meta_info : dict
        Metadata dictionary containing gating information with structure:
        meta_info['gating_summary'][run_ID][1][hierarchy]['marker_combo']
        
    run_ID : str
        Identifier for the current experimental run
        
    base_df : pandas.DataFrame
        Main dataframe containing cell measurements with columns:
        - marker columns (e.g., 'CD3', 'CD4', etc.)
        - 'true_label' : binary column (1=target, 0=non-target)
        - 'gate_hull_{i}' : binary gating results for hierarchy i
        
    population : str
        Name of the cell population being analyzed (used in plot title)
        
    population_path : str
        Path to directory containing gate edge CSV files with format:
        'cluster_{population}_gate_edges_hierarchy_{hierarchy}.csv'
        
    num_hierarchies : int, optional
        Number of gating hierarchies to plot (creates subplots). Default: 1
        
    pop_min : int, optional
        Minimum population size for contour plotting. Populations smaller than
        this threshold are displayed as scatter plots only. Default: 100
        
    mode : str or None, optional
        Axis scaling mode for transformed data:
        - None: Use original axis scaling
        - 'both': Add secondary axes with original scale
        - 'replace': Replace axes with original scale
        Default: None
        
    transform_func : callable or dictionary, optional
        Function to transform data from original to display space.
        If the transform function is different per marker, pass a dictionary 
        of functions with the marker name as keys.
        (e.g., arcsinh, biexp transformation). Default: None
        
    inverse_transform_func : callable or dictionary, optional
        Inverse transformation function for axis labeling. 
        If the transform function is different per marker, pass a dictionary 
        of functions with the marker name as keys.
        Default: None
        
    linear_threshold : float, optional
        Threshold for hybrid tick generation in original units. Default: 100
        
    outlier_percentile : float, optional
        Percentile threshold for outlier detection (lower = more outliers).
        Default: 1 (bottom 1% of density)
        
    grid_size : int, optional
        Resolution of density estimation grid. Higher values give smoother
        contours but slower computation. Default: 50
        
    target_color : str, optional
        Color for target population (hex or named color). Default: '#FF4F00'
        
    non_target_color : str, optional
        Color for non-target population (hex or named color). Default: '#2d5380'
        
    add_grid : bool, optional
        Whether to add grid lines to plots. Default: False
        
    save : str, optional
        File path to save plot (format determined by extension). Default: None
        
    save_png : str, optional
        File path to save PNG version of plot. Default: None
    
    Returns:
    --------
    None
        Displays plot and optionally saves to file
    
    Notes:
    ------
    **Contour Generation:**
    - Uses 2D histogram with Gaussian smoothing for fast density estimation
    - Adaptive contour levels: more levels for larger populations
    - Population-specific coloring: targets and non-targets get distinct colors
    
    **Outlier Detection:**
    - Identifies cells in lowest 5% of population density
    - Only applied to populations ≥ pop_min events
    - Uses histogram bin assignment for computational efficiency
    
    **Mixed Population Handling:**
    - Populations ≥ pop_min: Shows contours + outlier detection
    - Populations < pop_min: Shows scatter points only
    - Each population handled independently
    
    **Performance:**
    - Optimized for datasets with 10K-1M+ cells
    - O(n) complexity for outlier detection
    - Minimal memory footprint through histogram approach
    
    Examples:
    ---------
    Basic usage:
    >>> contour_lines_separate_outliers(
    ...     meta_info, run_ID, df, 'T_cells', './gates/',
    ...     num_hierarchies=2
    ... )
    
    With custom styling:
    >>> contour_lines_separate_outliers(
    ...     meta_info, run_ID, df, 'B_cells', './gates/',
    ...     target_color='red', non_target_color='blue',
    ...     grid_size=100, add_grid=True
    ... )
    
    Save publication figure:
    >>> contour_lines_separate_outliers(
    ...     meta_info, run_ID, df, 'NK_cells', './gates/',
    ...     save='./figures/NK_analysis.pdf',
    ...     save_png='./figures/NK_analysis.png'
    ... )
    """
    
    # Define a fixed size for each subplot
    fig_width = 5*num_hierarchies  # Adjust as needed
    fig_height = 5  # Adjust as needed
    fig, axes = plt.subplots(1, num_hierarchies, figsize=(fig_width, fig_height))  # Always create num_hierarchies subplots
    
    for idx, hierarchy in enumerate(map(str, range(1, num_hierarchies+1))):  # Loop through up to 5 hierarchies
        if idx >= num_hierarchies:  # Skip plotting if hierarchy number exceeds num_hierarchies
            axes[idx].axis('off')  # Turn off unused subplots
            continue
        
        # Extract marker information for current hierarchy
        marker = meta_info['gating_summary'][run_ID][1][hierarchy]['marker_combo']
        marker1 = marker[0]
        marker2 = marker[1]
        
        # Filter data for current gating hierarchy
        points_df = base_df[base_df['gate_hull_' + str(int(hierarchy)-1)] == 1][[marker1, marker2, 'true_label']]
        points_targets = points_df[points_df['true_label'] == 1]
        points_non_targets = points_df[points_df['true_label'] == 0]
        nr_targets = len(points_targets)
        nr_non_targets = len(points_non_targets)
        
        print('nr_targets :' +str(nr_targets))
        print('nr_non_targets :' +str(nr_non_targets))
        
        # Load gate edges from ConvexGating
        edges = pd.read_csv(os.path.join(population_path, 
                                     f'cluster_{population}_gate_edges_hierarchy_{hierarchy}.csv'), index_col=0)
        
        # Determine axes limits to display both cells and gate edges
        min_x = min(np.min(points_df[marker1]), np.min(edges['x_coordinate']))
        max_x = max(np.max(points_df[marker1]), np.max(edges['x_coordinate']))
        
        min_y = min(np.min(points_df[marker2]), np.min(edges['y_coordinate']))
        max_y = max(np.max(points_df[marker2]), np.max(edges['y_coordinate']))
        
        range_y = max_y - min_y
        range_x = max_x - min_x
        shift_y = 0.05 * range_y 
        shift_x = 0.05 * range_x
        
        x_plot_range = [min_x - shift_x, max_x + shift_x]
        y_plot_range = [min_y - shift_y, max_y + shift_y]
        
        # Create 2D histograms for each population separately
        all_x = points_df[marker1].values
        all_y = points_df[marker2].values
    
        # Define common bin edges for both populations
        x_min, x_max = all_x.min(), all_x.max()
        y_min, y_max = all_y.min(), all_y.max()
        x_edges = np.linspace(x_min, x_max, grid_size + 1)
        y_edges = np.linspace(y_min, y_max, grid_size + 1)
    
        # Create histograms for density estimation
        H_non_targets, _, _ = np.histogram2d(points_non_targets[marker1], points_non_targets[marker2], 
                                        bins=[x_edges, y_edges])
        H_targets, _, _ = np.histogram2d(points_targets[marker1], points_targets[marker2], 
                                    bins=[x_edges, y_edges])
    
        H_non_targets = H_non_targets.T
        H_targets = H_targets.T
    
        # Smooth histograms for better contour appearance
        H_non_targets_smooth = gaussian_filter(H_non_targets, sigma=0.8)
        H_targets_smooth = gaussian_filter(H_targets, sigma=0.8)
    
        # Create coordinate arrays for contour plotting
        X_centers = (x_edges[:-1] + x_edges[1:]) / 2
        Y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        X, Y = np.meshgrid(X_centers, Y_centers)
    
        # Select appropriate subplot
        ax = axes[idx] if num_hierarchies > 1 else axes
        
        # Adaptive contour levels based on population size
        if nr_non_targets > 0 :
            n_levels_non_targets = int(np.round(nr_non_targets/np.min([3000, nr_non_targets/4])))
        else:
            n_levels_non_targets = 1
        n_levels_targets = int(np.round(nr_targets/np.min([400, nr_targets/4])))
        
        levels_non_targets = np.linspace(H_non_targets_smooth.min(), H_non_targets_smooth.max(), 
                                         n_levels_non_targets)[1:]  # More detail for larger population
        levels_targets = np.linspace(H_targets_smooth.min(), H_targets_smooth.max(), 
                                     n_levels_targets)[1:]  # Fewer levels for smaller population
        
        # Handle non-target population
        if nr_non_targets > pop_min:
            # Draw contour lines for populations above threshold
            ax.contour(X, Y, H_non_targets_smooth, levels=levels_non_targets, 
                  colors=non_target_color, alpha=0.7, linewidths=1.0, linestyles='solid')
            
            # Fast outlier detection using histogram bins
            x_bin_indices_nt = np.digitize(points_non_targets[marker1], x_edges) - 1
            y_bin_indices_nt = np.digitize(points_non_targets[marker2], y_edges) - 1
            x_bin_indices_nt = np.clip(x_bin_indices_nt, 0, grid_size - 1)
            y_bin_indices_nt = np.clip(y_bin_indices_nt, 0, grid_size - 1)
    
            # Outlier detection: cells in lowest 5% of density
            point_densities_nt = H_non_targets_smooth[y_bin_indices_nt, x_bin_indices_nt]
            density_threshold_nt = np.percentile(point_densities_nt[point_densities_nt > 0], 5)
            non_target_outlier_mask = point_densities_nt <= density_threshold_nt
            
            # Plot outliers
            if np.any(non_target_outlier_mask):
                outlier_x_nt = points_non_targets[marker1][non_target_outlier_mask]
                outlier_y_nt = points_non_targets[marker2][non_target_outlier_mask]
                ax.scatter(outlier_x_nt, outlier_y_nt, c=non_target_color, s=2, alpha=0.1, 
                      label=f'Non-target outliers ({np.sum(non_target_outlier_mask)})', zorder=5)
        else:
            # Scatter plot for small populations
            ax.scatter(points_non_targets[marker1], points_non_targets[marker2], 
                  c=non_target_color, s=2, alpha=0.3, 
                  label=f'Non-targets ({nr_non_targets})', zorder=3)
            
        # Handle target population    
        if nr_targets > pop_min:
            # Draw contour lines for populations above threshold
            ax.contour(X, Y, H_targets_smooth, levels=levels_targets, 
                  colors=target_color, alpha=0.7, linewidths=1.0, linestyles='solid')
            
            # Fast outlier detection using histogram bins
            x_bin_indices_t = np.digitize(points_targets[marker1], x_edges) - 1
            y_bin_indices_t = np.digitize(points_targets[marker2], y_edges) - 1
            x_bin_indices_t = np.clip(x_bin_indices_t, 0, grid_size - 1)
            y_bin_indices_t = np.clip(y_bin_indices_t, 0, grid_size - 1)
    
            # Outlier detection: cells in lowest 5% of density
            point_densities_t = H_targets_smooth[y_bin_indices_t, x_bin_indices_t]
            density_threshold_t = np.percentile(point_densities_t[point_densities_t > 0], 5)
            target_outlier_mask = point_densities_t <= density_threshold_t
            
            # Plot outliers
            if np.any(target_outlier_mask):
                outlier_x_t = points_targets[marker1][target_outlier_mask]
                outlier_y_t = points_targets[marker2][target_outlier_mask]
                ax.scatter(outlier_x_t, outlier_y_t, c=target_color, s=2, alpha=0.1, 
                          label=f'Target outliers ({np.sum(target_outlier_mask)})', zorder=5)
        else: 
            # Scatter plot for small populations
            ax.scatter(points_targets[marker1], points_targets[marker2], 
                  c=target_color, s=2, alpha=0.3, 
                  label=f'Targets ({nr_targets})', zorder=3) 

        # Configure subplot appearance
        ax.set_xlabel(marker1)
        ax.set_ylabel(marker2)
        ax.set_title(f'Hierarchy {hierarchy}')
        
        if add_grid:
            ax.grid(True, alpha=0.3)
        
        # Set plot limits and add gate boundaries
        ax.set_xlim(x_plot_range)
        ax.set_ylim(y_plot_range)
        
        # Plot gate edges (close the polygon)
        edges = edges.append(edges.iloc[0], ignore_index=True)
        ax.plot(edges['x_coordinate'], edges['y_coordinate'], 
                linestyle='-', c='tab:red', alpha=1, linewidth=3)
        
        # Handle axis transformations for special scales (e.g., arcsinh)
        if mode is not None:
            # Apply custom axis scaling if transformation functions provided
            lims_y = np.array(ax.get_ylim())
            lims_x = np.array(ax.get_xlim())
            if (callable(transform_func) & callable(inverse_transform_func)):
                iticks_y, iticklabels_y = generate_hybrid_ticks_and_labels(
                    lims_y, linear_threshold=linear_threshold,
                    transform_func=transform_func,
                    inverse_transform_func=inverse_transform_func)

            
                iticks_x, iticklabels_x = generate_hybrid_ticks_and_labels(
                    lims_x, linear_threshold=linear_threshold,
                    transform_func=transform_func,
                    inverse_transform_func=inverse_transform_func)
            elif (isinstance(transform_func, dict) & isinstance(inverse_transform_func, dict)):
                iticks_y, iticklabels_y = generate_hybrid_ticks_and_labels(
                    lims_y, linear_threshold=linear_threshold,
                    transform_func=transform_func[marker2],
                    inverse_transform_func=inverse_transform_func[marker2])

            
                iticks_x, iticklabels_x = generate_hybrid_ticks_and_labels(
                    lims_x, linear_threshold=linear_threshold,
                    transform_func=transform_func[marker1],
                    inverse_transform_func=inverse_transform_func[marker1])
            
            if mode == 'both':
                # Add secondary axes with original scale
                ax1 = ax.twinx()
                ax1.set_ylim(lims_y[0], lims_y[1])
                ax1.set_yticks(iticks_y)
                ax1.set_yticklabels(iticklabels_y)
                
                ax1 = ax.twiny()
                ax1.set_xlim(lims_x[0],lims_x[1])
                ax1.set_xticks(iticks_x)
                ax1.set_xticklabels(iticklabels_x)
                
            elif mode == 'replace':
                # Replace current axis scaling
                ax.set_yticks(iticks_y)
                ax.set_yticklabels(iticklabels_y)
                ax.set_xticks(iticks_x)
                ax.set_xticklabels(iticklabels_x)
        
    # Finalize plot
    plt.suptitle(population, fontsize=12)
    plt.tight_layout()
    
    # Save plots if requested
    if save is not None:
        plt.savefig(save, bbox_inches='tight')
    if save_png is not None:
        plt.savefig(save_png, bbox_inches='tight')
    
    plt.show()
    
    return

def prepare_contour_plot(meta_info,add_hull_path,run_ID = 0):
    base_df = meta_info['general_summary'][run_ID]
    population = meta_info['clusterkeys'][run_ID]
    population_path = os.path.join(add_hull_path,population)
    if not os.path.exists(population_path):
        os.mkdir(population_path)
    base_df = add_tight_analysis(meta_info=meta_info,
                                run_ID=run_ID,
                                save_loc=population_path)
    return population,population_path,base_df