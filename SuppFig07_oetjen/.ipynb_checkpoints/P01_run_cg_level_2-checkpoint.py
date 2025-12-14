import convexgating as cg
import numpy as np
import pandas as pd
import os
import scanpy as sc

base_data = '/work/users/mh823zote/projects/gating/data/cytof/all/oetjen_all.h5ad'
level = 'cell_type_lvl2'
target_folder = '/work/users/mh823zote/projects/gating/data/oetjen/Figure_5/'
oetjen_all = sc.read_h5ad(base_data)
samples = list(pd.unique(oetjen_all.obs['sample']))
for sample in samples:
    oetjen = oetjen_all[oetjen_all.obs['sample'] == sample]
    sc.pp.subsample(oetjen,n_obs=50000)
    cts = pd.unique(oetjen.obs[level])
    if 'not annotated' in cts:
        cts = [item for item in cts if item != 'not annotated']
    if not os.path.exists(target_folder + 'P01_' + level):
        os.mkdir(target_folder + 'P01_' + level)
    subpath = os.path.join(target_folder + 'P01_' + level, sample)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    cg.tools.gating_strategy(oetjen,cts,cluster_string = level, save_path = subpath, add_noise=True)