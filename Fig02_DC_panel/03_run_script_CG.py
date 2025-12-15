#!/usr/bin/env python
# coding: utf-8

import warnings

import os
from pathlib import Path
import scanpy as sc
import anndata as ann
import pandas as pd
import numpy as np
import json

import time

import convexgating as cg

warnings.filterwarnings('ignore')

if __name__ == '__main__':
    
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Compute convex gating')

    parser.add_argument('-u', '--data', required=True, help='Path to data file in h5ad format')
    parser.add_argument('-o', '--output', required=True, help='Output path')
    parser.add_argument('-s', '--subsample', required=False, default=None, help='Number of cells to subsample to')

    parser.add_argument('-k', '--key', required=True, help='Key of annotated labels, e.g. "cell_type"')
    parser.add_argument('-b', '--batch', required=False, default=None, help='Key of batch labels when data contains several patients, e.g. "patient_ID"')
    parser.add_argument('-g', '--group', required=True, help='Group of annotated cells, e.g. "B cell"')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    
    verbose = args.verbose
    key = args.key
    group = args.group
    batch = args.batch
    subsample = int(args.subsample) if not None else 0
    save_path = args.output
    verbose = args.verbose
    
    
    run_convex_gating(data = args.data, 
                      key = key,
                      group = group,
                      batch = batch,
                      subsample = subsample,
                      save_path = save_path,
                      verbose = verbose
                     )
    
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