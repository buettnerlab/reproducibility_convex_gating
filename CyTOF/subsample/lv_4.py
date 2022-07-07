import scanpy as sc
import convex_gating as cg
from convex_gating import tools as ct
from convex_gating import helper as ch
import numpy as np
import pandas as pd
import os
path = os.getcwd()
sample = os.path.basename(path)
cluster_string = 'cell_type_lvl4'
#-----------------------------------------------------------------
targets = ['CCR7+ CD4+ T cell', 
    'CCR7+ CD8+ T cell', 
    'CCR7- CD4+ T cell',
    'CCR7- CD8+ T cell',
    'Double negative T cell', 
    'Double positive T cell', 
    'NK cell']
save_path = os.getcwd() + '/level_4'
tot_cells = 20000
#-----------------------------------------------------------------
for targ in targets:
    oetjen = sc.read_h5ad('oetjen_' + sample + '.h5ad')
    oetjen.X = oetjen.X + (np.random.rand(oetjen.X.shape[0],oetjen.X.shape[1])-0.5)/10000
    oetjen_targets =  oetjen[oetjen.obs[cluster_string] == targ]
    oetjen_non_targets =  oetjen[oetjen.obs[cluster_string] != targ]
    if len(oetjen) > tot_cells:
        ratio = len(oetjen_targets)/len(oetjen_non_targets)
        nr_non_targ_sub = int(np.floor(tot_cells/(1+ratio)))
        nr_targ_sub = tot_cells - nr_non_targ_sub
        sc.pp.subsample(oetjen_non_targets,n_obs = nr_non_targ_sub)
        sc.pp.subsample(oetjen_targets,n_obs = nr_targ_sub)
        oetjen = oetjen_targets.concatenate(oetjen_non_targets) 
    ct.gating_strategy(oetjen,targets,cluster_string,save_path=save_path)







