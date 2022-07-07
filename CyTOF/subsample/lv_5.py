import scanpy as sc
import convex_gating as cg
from convex_gating import tools as ct
from convex_gating import helper as ch
import numpy as np
import pandas as pd
import os
path = os.getcwd()
sample = os.path.basename(path)
cluster_string = 'cell_type_lvl5'
#-----------------------------------------------------------------
targets = ['CD8+ EM T cell',
'CD4+ TEMRA',
'Naive CD4+ T cell',
'CD8+ TE T cell',
'CD4+ CM T cell',
'CD4+ EM T cell',
'Naive CD8+ T cell',
'NK cell',
'Double negative T cell',
'CD4+ TRM T cell',
'CD8+ CM T cell',
'Double positive T cell',
'CD8+ TRM T cell']
save_path = os.getcwd() + '/level_5'
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







