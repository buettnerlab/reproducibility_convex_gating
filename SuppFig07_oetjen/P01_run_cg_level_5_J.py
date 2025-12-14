import convexgating as cg
import numpy as np
import pandas as pd
import os
import scanpy as sc

base_data = '/work/users/mh823zote/projects/gating/data/cytof/all/oetjen_all.h5ad'
level = 'cell_type_lvl5'
target_folder = '/work/users/mh823zote/projects/gating/data/oetjen/Figure_5/'
oetjen_all = sc.read_h5ad(base_data)
samples = ['J']
#samples = ['J','B']
for sample in samples:
    oetjen = oetjen_all[oetjen_all.obs['sample'] == sample]
    sc.pp.subsample(oetjen,n_obs=50000)
    cts = ['NK cell', 'Naive CD4+ T cell', 'CD8+ TE T cell', 'CD4+ CM T cell','Naive CD8+ T cell','CD8+ EM T cell','Double positive T cell', 
'CD4+ EM T cell','CD8+ CM T cell', 'CD4+ TRM T cell', 'Double negative T cell', 'CD8+ TRM T cell']
    if 'not annotated' in cts:
        cts = [item for item in cts if item != 'not annotated']
    if not os.path.exists(target_folder + 'P01_' + level):
        os.mkdir(target_folder + 'P01_' + level)
    subpath = os.path.join(target_folder + 'P01_' + level, sample)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    cg.tools.gating_strategy(oetjen,cts,cluster_string = level, save_path = subpath, add_noise=True)
