import convexgating as cg
import numpy as np
import pandas as pd
import os
import scanpy as sc

base_data = '/work/users/mh823zote/projects/gating/data/cytof/all/oetjen_all.h5ad'
level = 'cell_type_lvl5'
target_folder = '/work/users/mh823zote/projects/gating/data/oetjen/Figure_5/'
oetjen_all = sc.read_h5ad(base_data)
samples = ['B']
#samples = ['J','B']
for sample in samples:
    oetjen = oetjen_all[oetjen_all.obs['sample'] == 'B']

    oetjen_TRM = oetjen[oetjen.obs['cell_type_lvl5'] == 'CD8+ TRM T cell']

    oetjen_no_TRM = oetjen[oetjen.obs['cell_type_lvl5'] != 'CD8+ TRM T cell']

    sc.pp.subsample(oetjen_TRM,n_obs=100)
    sc.pp.subsample(oetjen_no_TRM,n_obs = 49900)
    oetjen =oetjen_TRM.concatenate(oetjen_no_TRM)
    cts = ['CD8+ TRM T cell']
    if 'not annotated' in cts:
        cts = [item for item in cts if item != 'not annotated']
    if not os.path.exists(target_folder + 'P01_' + level + 'TRM_' + sample):
        os.mkdir(target_folder + 'P01_' + level + 'TRM_' + sample)
    subpath = os.path.join(target_folder + 'P01_' + level + 'TRM_' + sample, sample)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    cg.tools.gating_strategy(oetjen,cts,cluster_string = level, save_path = subpath, add_noise=True)
