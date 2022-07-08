import scanpy as sc
import convex_gating as cg
from convex_gating import tools as ct
from convex_gating import helper as ch
import pandas as pd

#----------------------------------------------------------------------------------------------------------------#
file = sc.read_h5ad('Crtitical/crtitical.h5ad')
save_path = '/work/users/mh823zote/projects/gating/data/Stephenson_2/no_subsampling/Crtitical'
#----------------------------------------------------------------------------------------------------------------#
targets = list(pd.read_csv('target_clusters_sccoda.csv')['credible_effect_clusters'])
cluster_string = 'initial_clustering'
ct.gating_strategy(file,targets,cluster_string,save_path = save_path)