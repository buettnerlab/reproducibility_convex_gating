import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import mudata as md
import muon as mu
#import scvi

sc.logging.print_header()

file_path = '/work/users/mh823zote/projects/gating/data/CITEseq_Covid/' 

data_dir = file_path + 'data/' 

data_raw_dir = file_path + 'input/' 

sc.settings.figdir = file_path + 'figures/'

import datetime

today = datetime.date.today().strftime('%y%m%d') #creates a YYMMDD string of today's date

adata = md.read(f'{data_dir}haniffa21_totalVI_Covid_Healthy.h5mu/rna')

pdata = md.read(f'{data_dir}haniffa21_totalVI_Covid_Healthy.h5mu/protein')

adata_subset = md.read_h5ad(f'{data_dir}haniffa21_totalVI_Covid_Healthy.h5mu', mod='rna_subset')

mdata = md.MuData({"rna": adata, "protein": pdata, "rna_subset": adata_subset})

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    color=["rna_subset:leiden_totalVI", "rna_subset:full_clustering", 'rna_subset:initial_clustering'],
    frameon=False,
    ncols=1,
)

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    color=["rna_subset:leiden_totalVI", "rna_subset:full_clustering"],
    frameon=False,
    ncols=1,
    legend_loc='on data'
)

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    color=['AB_CD3', 'AB_CD4', 'AB_CD8', 'AB_CD16', 'AB_CD14'],
    frameon=False,
    ncols=3,
    vmax=300,
    wspace=0.1,
    cmap='Reds',
    size=5,
    layer="denoised_protein",
)



pdata.obs['leiden_totalVI']=mdata.mod['rna_subset'].obs['leiden_totalVI'].copy()

sc.pl.matrixplot(pdata, groupby='full_clustering', dendrogram=True,
                 var_names=pdata.var_names, standard_scale='var', log=True,
                 cmap='Reds', layer='denoised_protein')

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    color=pdata.var_names[:39],
    frameon=False,
    ncols=3,
    vmax="p99",
    wspace=0.1,
    cmap='Reds',
    size=5,
    layer="denoised_protein",
)

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    color=pdata.var_names[39:][:39],
    frameon=False,
    ncols=3,
    vmax="p99",
    wspace=0.1,
    cmap='Reds',
    size=5,
    layer="denoised_protein",
)

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    color=pdata.var_names[78:][:39],
    frameon=False,
    ncols=3,
    vmax="p99",
    wspace=0.1,
    cmap='Reds',
    size=5,
    layer="denoised_protein",
)

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    layer="protein_foreground_prob",
    color=['AB_CD3', 'AB_CD16', 'AB_CD4', 'AB_CD8'],
    frameon=False,
    ncols=2,
    vmax=0.9,
    size=3,
    wspace=0.1,
    color_map="cividis",
)

ax = sc.pl.scatter(pdata, x='AB_CD3', y='AB_CD16', color='full_clustering', show=False, size=5, layers='denoised_protein')
ax.set_xlim([-10, 500])
ax.set_ylim([-10, 500])
plt.show()

cd4_comp = []
cd8_comp = []

for group in adata.obs['initial_clustering'].cat.categories:
    
    if group.startswith('CD4'):
        cd4_comp.append(group)
    elif group.startswith('CD8'):
        cd8_comp.append(group)
        
pdata_t_cell = pdata[pdata.obs['initial_clustering'].isin(cd4_comp + cd8_comp)].copy()
prob_cd16 = pdata_t_cell.layers['protein_foreground_prob'][:, pdata.var_names =='AB_CD16'].squeeze()
pdata_t_cell.obs['prob_cd16'] = prob_cd16

protein = mdata.mod['protein']

prob_cd16 = pdata.layers['protein_foreground_prob'][:, pdata.var_names =='AB_CD16'].squeeze()

bool_mask = np.logical_and(prob_cd16>0.5,  pdata.obs['initial_clustering'].isin(cd4_comp))
bool_mask2 = np.logical_and(prob_cd16>0.5,  pdata.obs['initial_clustering'].isin(cd8_comp))

protein.obs['CD16+ T cells'] = 'other'
protein.obs['CD16+ T cells'][bool_mask] = 'CD16+ CD4+ T cell'
protein.obs['CD16+ T cells'][bool_mask2] = 'CD16+ CD8+ T cell'
protein.obs['CD16+ T cells'] = protein.obs['CD16+ T cells'].astype('category')

protein.uns['CD16+ T cells_colors'] = ['#fd349a', '#92fab1', '#bbbbbb']
mdata.update()

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    #layer="protein",
    color='protein:CD16+ T cells',
    groups=['CD16+ CD8+ T cell', 'CD16+ CD4+ T cell'],
    frameon=False,
    ncols=3,
    size=20,
    wspace=0.1)

mu.pl.embedding(
    mdata,
    basis="rna_subset:X_umap",
    #layer="protein",
    color='rna_subset:initial_clustering',
    groups=['CD8', 'CD4'],
    frameon=False,
    ncols=3,
    size=20,
    wspace=0.1)
mdata.update()
mdata.mod['protein'] = protein
mdata.write(f'{data_dir}haniffa21_totalVI_Covid_Healthy_annotated.h5mu')
