import os
from pathlib import Path
import scanpy as sc
import anndata as ann
import pandas as pd
import numpy as np
import json

from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
import time

import convexgating as cg
import convex_hull_add_on as ch
from apply_strategy_to_new_anndata import apply_gating_strategy


def run_convex_gating(data: ann.AnnData = 'adata.h5ad', 
                      key: 'str' = 'louvain',
                      group: 'str' = '0',
                      batch: 'str' = None,
                      channels: 'list' = [],
                      subsample: 'int' = 0,
                      add_noise: 'bool' = True,
                      save_path: 'str' = None,
                      verbose: 'int' = 0
                     ):
    
    #read data 
    adata = sc.read(data)
    
    #check inputs
    if verbose>0:
        print(f'Key {key}')
        print(f'Group {group}')
    
    
    if key not in adata.obs_keys():
        raise KeyError(f"Did not find '{key}' in `.obs_keys()`.")
    
    if group not in adata.obs[key].unique():
        raise KeyError(f"Did not find '{group}' in `.obs['{key}']`.")
    
    if channels == []:
        channels = list(adata.var.index)
    
    for channel in channels:
        if channel not in list(adata.var.index):
            raise KeyError(f"Did not find '{channel}' in `.var.index`.") 
    
    #transform group into list
    groups = []
    groups.append(group)
    
    if add_noise:
        adata.X = adata.X + (np.random.rand(adata.X.shape[0], adata.X.shape[1]) - 0.5) / 10000
    
    if batch is not None:
        if batch not in adata.obs_keys():
            raise KeyError(f"Did not find '{batch}' in `.obs_keys()`.")
        
        for batch_key in adata.obs[batch].unique():
            #subset to batch key
            adata_batch = adata[adata.obs[batch]==batch_key].copy()
            #skip iteration if cell type not present
            if group not in adata_batch.obs[key].unique():
                continue
            #subsample anndata object
            if subsample>0:
                sc.pp.subsample(adata_batch, n_obs=subsample)
    
            #pre-process data
            cell_data = cg.preprocess_adata_gating(adata_batch, cluster_string = key)
            #print(cell_data.head())
    
            if save_path is None:
                save_path = os.getcwd()
            #create output directory
            save_to_path='{}/cluster_{}_batch_{}/'.format(save_path, group, batch_key)
            Path(save_to_path).mkdir(parents=True, exist_ok=True)
    
            keys, gating_core, gating_overview = cg.FIND_GATING_STRATEGY(cell_data = cell_data,
                                                            channels = channels,
                                                            cluster_numbers = groups,
                                                            cluster_string = key,
                                                            save_path = save_to_path
                                                            )
            meta_info = {}
            meta_info["clusterkeys"] = keys
            meta_info["gating_summary"] = gating_core
            meta_info["general_summary"] = gating_overview
            np.save(os.path.join(save_to_path, "meta_info.npy"), meta_info)
            
            cg.tools.convex_hull_add_on(meta_info_path = os.path.join(save_to_path, 'meta_info.npy'),
                              target_location=save_to_path)
            
            for keyID in keys:
                save_to_path2='{}/cluster_{}_batch_{}/'.format(save_path, keys[keyID], batch_key)
                Path(save_to_path2).mkdir(parents=True, exist_ok=True)
                gating_overview[keyID].to_csv('{}gate_overview.csv'.format(save_to_path2))
    else:
        #subsample anndata object
        if subsample>0:
            sc.pp.subsample(adata, n_obs=subsample)
    
        #pre-process data
        cell_data = cg.preprocess_adata_gating(adata, cluster_string = key)
        #print(cell_data.head())
        #add batch info
        sample_ids = adata.obs['sample'].reset_index().drop(columns='index')
    
        if save_path is None:
            save_path = os.getcwd()
        #create output directory
        save_to_path='{}/cluster_{}/'.format(save_path, group)
        Path(save_to_path).mkdir(parents=True, exist_ok=True)
        
        keys, gating_core, gating_overview = cg.FIND_GATING_STRATEGY(cell_data = cell_data,
                                                            channels = channels,
                                                            cluster_numbers = groups,
                                                            cluster_string = key,
                                                            save_path = save_to_path
                                                            )
        #add batch info to metadata output
        gating_overview['sample'] = sample_ids
        meta_info = {}
        meta_info["clusterkeys"] = keys
        meta_info["gating_summary"] = gating_core
        meta_info["general_summary"] = gating_overview
        np.save(os.path.join(save_to_path, "meta_info.npy"), meta_info)
        
        cg.tools.convex_hull_add_on(meta_info_path = os.path.join(save_to_path, 'meta_info.npy'),
                              target_location=save_to_path)
        
        for keyID in keys:
            save_to_path2='{}/cluster_{}/'.format(save_path, keys[keyID])
            Path(save_to_path2).mkdir(parents=True, exist_ok=True)
            gating_overview[keyID].to_csv('{}gate_overview.csv'.format(save_to_path2))
    
    return

def run_svm(data: ann.AnnData = 'adata.h5ad', 
            key: 'str' = 'louvain',
            batch: 'str' = None,
            subsample: 'int' = 0,
            svm_type: 'str' = 'linear', #choose from linear and rbf
            save_path: 'str' = None,
            verbose: 'int' = 0):
    def _run_svm(adata, key, batch_key, targets, subsample, svm_type):
        #initialize result table
        cols = ['dataset','sample','celltype', 'subsample','f1','recall','precision']
        svm_df = pd.DataFrame(columns=cols)
        
        #subsample anndata object
        if subsample>0:
            sc.pp.subsample(adata_batch, n_obs=subsample)
        
        #prepare data
        cell_data = cg.preprocess_adata_gating(adata_batch, cluster_string = key)
        cell_data = cell_data.drop("cell_ID",axis=1)
        #iterate over all cell types
        for target_cluster in targets:
            y_true = (cell_data[key] == target_cluster).values*1
            if sum(y_true) != 0:
                #create data matrix for SVM
                X_true = cell_data.drop(key,axis=1).values
                #initialize SVM and fit model
                if svm_type=='linear':
                    svm = LinearSVC()
                elif svm_type=='rbf':
                    svm = SVC(kernel = 'rbf')
                svm.fit(X_true,y_true)
                #compute predictions
                y_pred = svm.predict(X_true)
                #compute resulting metrics
                f1 = f1_score(y_true, y_pred)
                recall = recall_score(y_true, y_pred)
                precision = precision_score(y_true, y_pred)
                #generate output tables
                ser = pd.Series((prefix, batch_key, target_cluster, subsample ,f1,recall,precision), 
                                    index = cols)
                svm_df = svm_df.append(ser,ignore_index=True)
            else:
                ser = pd.Series((prefix, batch_key, target_cluster, subsample,'na','na','na'),
                                  index = cols)
                svm_df = svm_df.append(ser,ignore_index=True)
        
        return svm_df
        
    #read data 
    adata = sc.read(data)
    prefix = Path(data).stem
    
   
    #check inputs
    if verbose>0:
        print(f'Key {key}')
        print(f'Batch {batch}')
        print(f'SVM type {svm_type}')
    
    svm_type = svm_type.lower()
    
    if svm_type not in ['linear', 'rbf']:
        raise KeyError(f"The argument '{svm_type}' needs to be either 'linear' or 'rbf'.")
        
    if key not in adata.obs_keys():
        raise KeyError(f"Did not find '{key}' in `.obs_keys()`.")
    
    if save_path is None:
            save_path = os.getcwd()
    
    #get groups 
    targets = adata.obs[key].unique()
      
    if batch is not None:
        if batch not in adata.obs_keys():
            raise KeyError(f"Did not find '{batch}' in `.obs_keys()`.")
        
        #initialize result table
        cols = ['dataset','sample','celltype', 'subsample','f1','recall','precision']
        svm_df = pd.DataFrame(columns=cols)
        
        for batch_key in adata.obs[batch].unique():
            #subset to batch key
            adata_batch = adata[adata.obs[batch]==batch_key].copy()
            svm_batch = _run_svm(adata = adata_batch, 
                                     key = key,
                                     batch_key = batch_key, 
                                     targets = targets, 
                                     subsample = subsample,
                                 svm_type = svm_type
                                )
            svm_df = svm_df.append(svm_batch, ignore_index=True)
            
    else:
        svm_df = _run_svm(adata = adata, 
                                     key = key,
                                     batch_key = batch, 
                                     targets = targets, 
                                     subsample = subsample,
                             svm_type = svm_type)
    
    #save to file
    ## check if save path exists and create if necessary
    Path(save_path).mkdir(parents=True, exist_ok=True)     
    svm_df.to_csv('{}/{}_{}_{}_svm_overview.csv'.format(save_path, prefix, key, svm_type))
    
    return





def load_cell_names_from_r(cell_names_file):
    """
    Load train/test cell names saved from R analysis.
    
    Args:
        cell_names_file: Path to .json file containing cell names
        
    Returns:
        Dictionary containing train/test cell names for each sample/cell_type combination
    """
    cell_names_path = Path(cell_names_file)
    
    if cell_names_path.suffix == '.json':
        # Load from JSON file (recommended format)
        with open(cell_names_file, 'r') as f:
            cell_names = json.load(f)
        return cell_names

    else:
        raise ValueError("Cell names file must be either .json or .rds format. JSON is recommended.")
    
    return

def load_h5ad_data(file_path):
    """Load h5ad data using scanpy."""
    adata = sc.read_h5ad(file_path)
    return adata

def run_convexgating_with_hypergate_cell_names(adata, cell_names, sample_key, cell_type_key, 
                                              sample_name, cell_type_name, 
                                              save_path, add_noise=True, verbose=True):
    """
    Run ConvexGating using the same train/test split as HyperGate.
    
    Args:
        adata: AnnData object containing the single-cell data
        cell_names: Dictionary of train/test cell names from R HyperGate analysis
        sample_key: Key for sample information in adata.obs
        cell_type_key: Key for cell type information in adata.obs
        sample_name: Name of the sample to analyze
        cell_type_name: Name of the cell type to analyze
        save_path: Base path to save ConvexGating results
        add_noise: Whether to add noise to training data
        verbose: Whether to print progress information
        
    Returns:
        Dictionary containing ConvexGating results and performance metrics
    """
    
    # Get the specific cell names for this sample/cell_type combination
    split_key = f"{sample_name}_{cell_type_name}"
    if split_key not in cell_names:
        raise ValueError(f"No cell names found for {split_key}")
    
    cell_names_group = cell_names[split_key]
    
    # Extract train, test, and full sample data using cell names
    train_cell_names = set(cell_names_group['all_train_cell_names'])
    test_cell_names = set(cell_names_group['all_test_cell_names'])
    sample_cell_names = set(cell_names_group['sample_cell_names'])
    
    # Filter adata using cell name intersections
    adata_train = adata[list(adata.obs_names & train_cell_names)].copy()
    adata_test = adata[list(adata.obs_names & test_cell_names)].copy()
    adata_full = adata[adata.obs[sample_key] == sample_name].copy()
    
    if verbose:
        print(f"Processing {sample_name} - {cell_type_name}")
        print(f"Training cells: {len(adata_train)}")
        print(f"Test cells: {len(adata_test)}")
        print(f"Full sample cells: {len(adata_full)}")
        print(f"Train ratio: {cell_names_group['train_ratio']}")
        print(f"Ratio non-target/target: {cell_names_group['ratio_nontarget_target']}")
    
    # Transform group into list
    groups = []
    groups.append(cell_type_name)
    
    # Add noise to training data if requested
    if add_noise:
        adata_train.X = adata_train.X + (np.random.rand(adata_train.X.shape[0], adata_train.X.shape[1]) - 0.5) / 10000
        if verbose:
            print("Added noise to training data")
    
    # Pre-process data using ConvexGating
    if verbose:
        print("Preprocessing data with ConvexGating...")
    
    cell_data = cg.preprocess_adata_gating(adata_train, cluster_string=cell_type_key)
    
    # Create output directory
    save_to_path = f'{save_path}/cluster_{cell_type_name}_batch_{sample_name}/'
    Path(save_to_path).mkdir(parents=True, exist_ok=True)
    
    # Train ConvexGating model
    if verbose:
        print("Training ConvexGating model...")
    
    start_time = time.time()
    
    keys, gating_core, gating_overview = cg.FIND_GATING_STRATEGY(
        cell_data=cell_data,
        channels=list(adata.var.index),
        cluster_numbers=groups,
        cluster_string=cell_type_key,
        save_path=save_to_path,
        show_HEAT=False,
        save_HEAT=False,
        show_SCATTER=False,
        save_SCATTER=False,
        show_metrics_df=False,
        save_metrics_df=False,
        save_metrics_plot=False,
        show_metrics_plot=False,
    )
    
    # Create and save meta info
    meta_info = {
        "clusterkeys": keys,
        "gating_summary": gating_core,
        "general_summary": gating_overview
    }
    np.save(os.path.join(save_to_path, "meta_info.npy"), meta_info)
    
    # Adjust convex hull
    cg.tools.convex_hull_add_on(
        meta_info_path=os.path.join(save_to_path, 'meta_info.npy'),
        target_location=save_to_path
    )
    
    train_time = time.time() - start_time
    
    # Save convex hull results
    for keyID in keys:
        save_to_path2 = f'{save_path}/cluster_{keys[keyID]}_batch_{sample_name}/'
        Path(save_to_path2).mkdir(parents=True, exist_ok=True)
        gating_overview[keyID].to_csv(f'{save_to_path2}gate_overview.csv')
    
    if verbose:
        print(f"Training completed in {train_time:.3f} seconds")
    
    # Apply gates to test data
    if verbose:
        print("Applying gates to test data...")
    
    start_time = time.time()
    
    try:
        f1_test, recall_test, precision_test = apply_gating_strategy(
            adata=adata_test,
            base_path=f'{save_path}/cluster_{cell_type_name}_batch_{sample_name}/',
            cluster=cell_type_name,
            cluster_string=cell_type_key,
            output_scores=True
        )
    except NameError:
        print("Warning: apply_gating_strategy function not found. Please import it.")
        f1_test = recall_test = precision_test = 0.0
    
    test_time = time.time() - start_time
    
    # Apply gates to full sample
    if verbose:
        print("Applying gates to full sample...")
    
    start_time = time.time()
    
    try:
        f1_full, recall_full, precision_full = apply_gating_strategy(
            adata=adata_full,
            base_path=f'{save_path}/cluster_{cell_type_name}_batch_{sample_name}/',
            cluster=cell_type_name,
            cluster_string=cell_type_key,
            output_scores=True
        )
    except NameError:
        print("Warning: apply_gating_strategy function not found. Please import it.")
        f1_full = recall_full = precision_full = 0.0
    
    full_time = time.time() - start_time
    
    if verbose:
        print(f"Test evaluation completed in {test_time:.3f} seconds")
        print(f"Full evaluation completed in {full_time:.3f} seconds")
        print(f"Test set - Precision: {precision_test:.3f}, Recall: {recall_test:.3f}, F1: {f1_test:.3f}")
        print(f"Full sample - Precision: {precision_full:.3f}, Recall: {recall_full:.3f}, F1: {f1_full:.3f}")
    
    # Prepare results dictionary
    result_dict_cluster_sample = {}
    result_dict_cluster_sample['test'] = {
        'cluster': cell_type_name,
        'f1': f1_test,
        'recall': recall_test,
        'precision': precision_test
    }
    
    result_dict_cluster_sample['full'] = {
        'cluster': cell_type_name,
        'f1': f1_full,
        'recall': recall_full,
        'precision': precision_full
    }
    
    # Save timing
    timing = {}
    timing['train'] = {
        'cluster': cell_type_name,
        'time': train_time
    }
    timing['test'] = {
        'cluster': cell_type_name,
        'time': test_time
    }
    timing['full'] = {
        'cluster': cell_type_name,
        'time': full_time
    }
    
    df_timing = pd.DataFrame(timing).T.reset_index().rename(columns={'index': 'set'})
    df_timing.to_csv(f'{save_to_path}timing.csv', index=False)
    
    # Save test performance to file
    test_perf = pd.DataFrame(result_dict_cluster_sample).T.reset_index().rename(
        columns={'index': 'set'}).melt(['set', 'cluster'])
    test_perf.to_csv(f'{save_to_path}performance_test.csv', index=False)
    
    return {
        'results': result_dict_cluster_sample,
        'timing': timing,
        'keys': keys,
        'train_ratio': cell_names_group['train_ratio'],
        'ratio_nontarget_target': cell_names_group['ratio_nontarget_target'],
        'cell_counts': {
            'train': adata_train.n_obs,
            'test': adata_test.n_obs,
            'full': adata_full.n_obs
        }
    }

def run_svm_with_hypergate_cell_names(adata, cell_names, sample_key, cell_type_key, 
                                      sample_name, cell_type_name, 
                                     svm_type: 'str' = 'linear', #choose from linear and rbf
                                     verbose = True
                                     ):
    """
    Run ConvexGating using the same train/test split as HyperGate.
    
    Args:
        adata: AnnData object containing the single-cell data
        cell_names: Dictionary of train/test cell names from R HyperGate analysis
        sample_key: Key for sample information in adata.obs
        cell_type_key: Key for cell type information in adata.obs
        sample_name: Name of the sample to analyze
        cell_type_name: Name of the cell type to analyze
        svm_type: Model type of the svm, choose from linear and rbf
        verbose: Whether to print progress information
        
    Returns:
        Dictionary containing ConvexGating results and performance metrics
    """
    
    # Get the specific cell names for this sample/cell_type combination
    split_key = f"{sample_name}_{cell_type_name}"
    if split_key not in cell_names:
        raise ValueError(f"No cell names found for {split_key}")
    
    cell_names_group = cell_names[split_key]
    
    # Extract train, test, and full sample data using cell names
    train_cell_names = set(cell_names_group['all_train_cell_names'])
    test_cell_names = set(cell_names_group['all_test_cell_names'])
    sample_cell_names = set(cell_names_group['sample_cell_names'])
    
    # Filter adata using cell name intersections
    adata_train = adata[list(adata.obs_names & train_cell_names)].copy()
    adata_test = adata[list(adata.obs_names & test_cell_names)].copy()
    adata_full = adata[adata.obs[sample_key] == sample_name].copy()
    
    train_ratio = cell_names_group['train_ratio']
    ratio_nontarget_target = cell_names_group['ratio_nontarget_target']
    actual_training_size = cell_names_group['actual_training_size']
    requested_training_size = cell_names_group['requested_training_size']
    
    if verbose:
        print(f"Processing {sample_name} - {cell_type_name}")
        print(f"Train ratio: {train_ratio}")
        print(f"Ratio non-target/target: {ratio_nontarget_target}")
        print(f"Requested training size: {requested_training_size}")
        print(f"Actual training size: {actual_training_size}")
        print(f"Training target cells: {len(train_target_cell_names)}")
        #print(f"Training non-target cells: {len(train_nontarget_cell_names)}")
        print(f"Test target cells: {len(test_target_cell_names)}")
       # print(f"Test non-target cells: {len(test_nontarget_cell_names)}")
        print(f"Total sample cells: {len(sample_cell_names)}")
    
    # Initialize result table
    cols = ['dataset','sample','celltype', 'time', 'set','f1','recall','precision']
    svm_df = pd.DataFrame(columns=cols)
    
    # Prepare data
    cell_data = cg.preprocess_adata_gating(adata_train, cluster_string = cell_type_key)
    cell_data = cell_data.drop("cell_ID",axis=1)
    X_train = cell_data.drop(cell_type_key,axis=1).values
    y_train = (cell_data[cell_type_key] == cell_type_name).values*1
    X_test = adata_test.X
    y_test = (adata_test.obs[cell_type_key] == cell_type_name).values*1
    X_full = adata_full.X
    y_full = (adata_full.obs[cell_type_key] == cell_type_name).values*1
    
    # Train svm per cell type

    if sum(y_train) != 0:
        # Create data matrix for SVM

        # Initialize SVM and fit model
        if verbose:
            print("Fitting SVM model...")
        start_time = time.time()
        if svm_type=='linear':
            svm = LinearSVC()
        elif svm_type=='rbf':
            svm = SVC(kernel = 'rbf')
        svm.fit(X_train, y_train)
        # Compute prediction on training data
        y_pred = svm.predict(X_train)
        # Get train time
        train_time = time.time() - start_time
        # Compute performance metrics
        f1_train = f1_score(y_train, y_pred)
        recall_train = recall_score(y_train, y_pred)
        precision_train = precision_score(y_train, y_pred)
        # Generate output tables
        ser_train = pd.Series((prefix, sample_name, target_cluster, train_time, 'train',
                         f1_train, recall_train, precision_train), 
                                    index = cols)
        svm_df = svm_df.append(ser_train,ignore_index=True)
         
        # Compute prediction on test data
        start_time = time.time()
        y_pred = svm.predict(X_test)
        # Get train time
        test_time = time.time() - start_time
        
        # Compute performance metrics on test data
        f1_test = f1_score(y_test, y_pred)
        recall_test = recall_score(y_test, y_pred)
        precision_test = precision_score(y_test, y_pred)
        # Generate output tables
        ser_test = pd.Series((prefix, sample_name, target_cluster, test_time, 'test',
                         f1_test, recall_test, precision_test), 
                                    index = cols)
        svm_df = svm_df.append(ser_test,ignore_index=True)
        
        # Compute prediction on full data
        start_time = time.time()
        y_pred = svm.predict(X_full)
        # Get train time
        full_time = time.time() - start_time
        
        # Compute performance metrics on test data
        f1_full = f1_score(y_full, y_pred)
        recall_full = recall_score(y_full, y_pred)
        precision_full = precision_score(y_full, y_pred)
        # Generate output tables
        ser_full = pd.Series((prefix, sample_name, target_cluster, full_time, 'full',
                         f1_full, recall_full, precision_full), 
                                    index = cols)
        svm_df = svm_df.append(ser_full, ignore_index=True)
        
    else:
        ser_train = pd.Series((prefix, sample_name, target_cluster, 0, 'train',0,0,0),
                                  index = cols)
        svm_df = svm_df.append(ser_train, ignore_index=True)
    
    if verbose:
        print(f"Training time: {train_time:.3f} seconds")
        print(f"Application time: {test_time:.3f} seconds")
        print(f"Total time: {full_time:.3f} seconds")
        print(f"Full sample - Precision: {precision_full:.3f}, Recall: {recall_full:.3f}, F1: {f1_full:.3f}")
        print(f"Training set - Precision: {precision_train:.3f}, Recall: {recall_train:.3f}, F1: {f1_train:.3f}")
        print(f"Test set - Precision: {precision_test:.3f}, Recall: {recall_test:.3f}, F1: {f1_test:.3f}")
    
    # save to file
    prefix = Path(adata).stem
    ## check if save path exists and create if necessary
    Path(save_path).mkdir(parents=True, exist_ok=True)     
    svm_df.to_csv(f'{save_path}/{prefix}_{sample_name}_{cell_type_key}_{cell_type_name}_{svm_type}_svm_performance.csv')
    
    return