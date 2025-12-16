import json
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
import time

import sys
import os
import re
from pathlib import Path

import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd

import convexgating as cg
import convex_hull_add_on as ch
from apply_strategy_to_new_anndata import apply_gating_strategy
from utils import run_convexgating_with_hypergate_cell_names, load_cell_names_from_r, load_h5ad_data

# File paths (adjust these to your actual file paths)
train_dir = '../data/'
cluster_key = "cell_type_lvl"
project_dir = f"{train_dir}train_test_split/hypergate/Oetjen/{cluster_key}/"
h5ad_file = train_dir + 'Oetjen_2018/anndata/cytof_data_tmp.h5ad'
cell_names_file = os.path.join(project_dir, 
                                f'cytof_data_tmp_{cluster_key}_train_test_cell_names.json')  

sample_key = "sample"
save_path = f"{train_dir}train_test_split_1to1/convex_gating/Oetjen/{cluster_key}/"

add_noise = True

# Load data and cell names
adata = load_h5ad_data(h5ad_file)
cell_names = load_cell_names_from_r(cell_names_file)
    
# Print structure of loaded cell names for debugging
print("Loaded cell names structure:")
for key in list(cell_names.keys())[:2]:  # Show first 2 entries
    print(f"\nKey: {key}")
    entry = cell_names[key]
    print(f"Keys in entry: {list(entry.keys())}")
    print(f"Number of train target cells: {len(entry['train_target_cell_names'])}")
    print(f"Number of train non-target cells: {len(entry['train_nontarget_cell_names'])}")
    break  # Just show one example
    
# Example: Run ConvexGating for a specific sample and cell type 
sample_cell_type_combo = list(cell_names.keys())
for sample_name, cell_type_name in [item.split('_') for item in sample_cell_type_combo]:
    # Run ConvexGating with the same cell names as HyperGate 
    cg_results = run_convexgating_with_hypergate_cell_names(
        adata=adata,
        cell_names=cell_names,
        sample_key=sample_key,
        cell_type_key=cluster_key,
        sample_name=sample_name,
        cell_type_name=cell_type_name,
        save_path=save_path,
        add_noise=add_noise,
        verbose=True
    )
