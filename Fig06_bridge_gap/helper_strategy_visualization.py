import os 
import numpy as np
import convexgating as cg
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kde

def plot_hierarchies_v3(meta_info, run_ID, base_df, population, population_path,plot_config_targets,plot_config_non_targets, num_hierarchies=1,min_contour = 200,save = None,save_png=None):
    # Define a fixed size for each subplot
    fig_width = 15  # Adjust as needed
    fig_height = 3  # Adjust as needed
    fig, axes = plt.subplots(1, 5, figsize=(fig_width, fig_height))  # Always create 5 subplots
    
    for idx, hierarchy in enumerate(map(str, range(1, 6))):  # Loop through up to 5 hierarchies
        if idx >= num_hierarchies:  # Skip plotting if hierarchy number exceeds num_hierarchies
            axes[idx].axis('off')  # Turn off unused subplots
            continue
        
        marker = meta_info['gating_summary'][run_ID][1][hierarchy]['marker_combo']
        marker1 = marker[0]
        marker2 = marker[1]
        
        points_df = base_df[base_df['gate_hull_' + str(int(hierarchy)-1)] == 1][[marker1, marker2, 'true_label']]
        points_targets = points_df[points_df['true_label'] == 1]
        points_non_targets = points_df[points_df['true_label'] == 0]
        nr_targets = len(points_targets)
        nr_non_targets = len(points_non_targets)
        print('nr_targets :' +str(nr_targets))
        print('nr_non_targets :' +str(nr_non_targets))
        
        
        x = np.linspace(points_df[marker1].min(), points_df[marker1].max(), 150)
        y = np.linspace(points_df[marker2].min(), points_df[marker2].max(), 150)
        X, Y = np.meshgrid(x, y)
        
        kde_0 = kde.gaussian_kde([points_non_targets[marker1], points_non_targets[marker2]])
        Z_0 = kde_0([X.flatten(), Y.flatten()]).reshape(X.shape)
        
        kde_1 = kde.gaussian_kde([points_targets[marker1], points_targets[marker2]])
        Z_1 = kde_1([X.flatten(), Y.flatten()]).reshape(X.shape)
        
        edges = pd.read_csv(os.path.join(population_path, f'cluster_{population}_gate_edges_hierarchy_{hierarchy}.csv'), index_col=0)
        
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
        
        ax = axes[idx]  # Select the appropriate subplot
        if plot_config_non_targets['hierarchy_' + str(hierarchy)] == 'contour':
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
        if plot_config_targets['hierarchy_' + str(hierarchy)] == 'contour':
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')

        if plot_config_non_targets['hierarchy_' + str(hierarchy)] == 'scatter':
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.4, markersize=1.5)
        if plot_config_targets['hierarchy_' + str(hierarchy)] == 'scatter':
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.4, markersize=1.5)

        if plot_config_non_targets['hierarchy_' + str(hierarchy)] == 'both':
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.05, markersize=1.5)
            
        if plot_config_targets['hierarchy_' + str(hierarchy)] == 'both':
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.05, markersize=1.5)
           

        '''
        # Contour plot for label 0
        if (nr_non_targets >= min_contour) & (nr_non_targets < 2000):
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.15, markersize=1.5)
        if nr_non_targets >= 2000:
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.05, markersize=1.5)
        if (nr_non_targets < min_contour):
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.4, markersize=1.5)
        
        # Contour plot for label 1
        if (nr_targets >= min_contour) & (nr_targets < 2000):
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.15, markersize=1.5)
        if (nr_targets >= 2000):
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.05, markersize=1.5)
        if (nr_targets < min_contour):
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.4, markersize=1.5)
        '''
        
        # Adding labels and title
        ax.set_xlabel(marker1)
        ax.set_ylabel(marker2)
        ax.set_title(f'Hierarchy {hierarchy}')
        
        #ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.2, markersize=1.5)
        #ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.2, markersize=1.5)
        ax.set_xlim(x_plot_range)
        ax.set_ylim(y_plot_range)
        
        edges = edges.append(edges.iloc[0], ignore_index=True)
        ax.plot(edges['x_coordinate'], edges['y_coordinate'], linestyle='-', c='tab:red', alpha=1, linewidth=3)
    
    plt.suptitle(population, fontsize=12)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    if save_png is not None:
        plt.savefig(save_png,bbox_inches='tight')
    plt.show()

def plot_hierarchies_v4(meta_info, run_ID, base_df, population, population_path,plot_config_targets,plot_config_non_targets, num_hierarchies=1,min_contour = 200,save = None,save_png=None,hierarchies_exclude_kde = None):
    # Define a fixed size for each subplot
    fig_width = 15  # Adjust as needed
    fig_height = 3  # Adjust as needed
    fig, axes = plt.subplots(1, 5, figsize=(fig_width, fig_height))  # Always create 5 subplots
    
    for idx, hierarchy in enumerate(map(str, range(1, 6))):  # Loop through up to 5 hierarchies
        if idx >= num_hierarchies:  # Skip plotting if hierarchy number exceeds num_hierarchies
            axes[idx].axis('off')  # Turn off unused subplots
            continue
        
        marker = meta_info['gating_summary'][run_ID][1][hierarchy]['marker_combo']
        marker1 = marker[0]
        marker2 = marker[1]
        
        points_df = base_df[base_df['gate_hull_' + str(int(hierarchy)-1)] == 1][[marker1, marker2, 'true_label']]
        points_targets = points_df[points_df['true_label'] == 1]
        points_non_targets = points_df[points_df['true_label'] == 0]
        nr_targets = len(points_targets)
        nr_non_targets = len(points_non_targets)
        print('nr_targets :' +str(nr_targets))
        print('nr_non_targets :' +str(nr_non_targets))
        
        
        x = np.linspace(points_df[marker1].min(), points_df[marker1].max(), 150)
        y = np.linspace(points_df[marker2].min(), points_df[marker2].max(), 150)
        X, Y = np.meshgrid(x, y)

        if hierarchy not in hierarchies_exclude_kde:
            kde_0 = kde.gaussian_kde([points_non_targets[marker1], points_non_targets[marker2]])
        Z_0 = kde_0([X.flatten(), Y.flatten()]).reshape(X.shape)

        if hierarchy not in hierarchies_exclude_kde:
            kde_1 = kde.gaussian_kde([points_targets[marker1], points_targets[marker2]])
        Z_1 = kde_1([X.flatten(), Y.flatten()]).reshape(X.shape)
        
        edges = pd.read_csv(os.path.join(population_path, f'cluster_{population}_gate_edges_hierarchy_{hierarchy}.csv'), index_col=0)
        
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
        
        ax = axes[idx]  # Select the appropriate subplot
        if plot_config_non_targets['hierarchy_' + str(hierarchy)] == 'contour':
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
        if plot_config_targets['hierarchy_' + str(hierarchy)] == 'contour':
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')

        if plot_config_non_targets['hierarchy_' + str(hierarchy)] == 'scatter':
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.4, markersize=1.5)
        if plot_config_targets['hierarchy_' + str(hierarchy)] == 'scatter':
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.4, markersize=1.5)

        if plot_config_non_targets['hierarchy_' + str(hierarchy)] == 'both':
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.05, markersize=1.5)
            
        if plot_config_targets['hierarchy_' + str(hierarchy)] == 'both':
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.05, markersize=1.5)
           

        '''
        # Contour plot for label 0
        if (nr_non_targets >= min_contour) & (nr_non_targets < 2000):
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.15, markersize=1.5)
        if nr_non_targets >= 2000:
            ax.contour(X, Y, Z_0, colors='#2d5380', linestyles='solid', label='Label 0')
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.05, markersize=1.5)
        if (nr_non_targets < min_contour):
            ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.4, markersize=1.5)
        
        # Contour plot for label 1
        if (nr_targets >= min_contour) & (nr_targets < 2000):
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.15, markersize=1.5)
        if (nr_targets >= 2000):
            ax.contour(X, Y, Z_1, colors='#FF4F00', linestyles='solid', label='Label 1')
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.05, markersize=1.5)
        if (nr_targets < min_contour):
            ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.4, markersize=1.5)
        '''
        
        # Adding labels and title
        ax.set_xlabel(marker1)
        ax.set_ylabel(marker2)
        ax.set_title(f'Hierarchy {hierarchy}')
        
        #ax.plot(points_non_targets[marker1], points_non_targets[marker2], 'o', color='#2d5380', alpha=0.2, markersize=1.5)
        #ax.plot(points_targets[marker1], points_targets[marker2], 'o', color='#FF4F00', alpha=0.2, markersize=1.5)
        ax.set_xlim(x_plot_range)
        ax.set_ylim(y_plot_range)
        
        edges = edges.append(edges.iloc[0], ignore_index=True)
        ax.plot(edges['x_coordinate'], edges['y_coordinate'], linestyle='-', c='tab:red', alpha=1, linewidth=3)
    
    plt.suptitle(population, fontsize=12)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    if save_png is not None:
        plt.savefig(save_png,bbox_inches='tight')
    plt.show()

def prepare_contour_plot(meta_info,add_hull_path,run_ID =0):
    base_df = meta_info['general_summary'][run_ID]
    population = meta_info['clusterkeys'][run_ID]
    population_path = os.path.join(add_hull_path,population)
    if not os.path.exists(population_path):
        os.mkdir(population_path)
    base_df = cg.tools.add_tight_analysis(meta_info=meta_info,
                                run_ID=run_ID,
                                save_loc=population_path)
    return population,population_path,base_df