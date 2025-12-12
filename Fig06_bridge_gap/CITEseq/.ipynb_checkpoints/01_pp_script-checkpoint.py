import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rcParams
import mudata as md
#import muon as mu
import scvi

sc.logging.print_header()

file_path = '/work/users/mh823zote/projects/gating/data/CITEseq_Covid/' 
data_dir = file_path + 'data/' 
data_raw_dir = file_path + 'input/' 
sc.settings.figdir = file_path + 'figures/'

import datetime

today = datetime.date.today().strftime('%y%m%d') #creates a YYMMDD string of today's date

adata = sc.read(data_raw_dir + 'covid_portal_210320_with_raw.h5ad')

cd4_comp = []
cd8_comp = []

for group in adata.obs['full_clustering'].cat.categories:
    
    if group.startswith('CD4.'):
        cd4_comp.append(group)
    elif group.startswith('CD8.'):
        cd8_comp.append(group)
        
adata_t_cell = adata[adata.obs['full_clustering'].isin(cd4_comp + cd8_comp)].copy()

cd16pos_cells = adata[(adata[:, adata.var_names=='AB_CD16'].layers['raw']>200).todense()].copy()

cd16pos_tcells = adata_t_cell[(adata_t_cell[:, adata_t_cell.var_names=='AB_CD16'].layers['raw']>200).todense()].copy()

#TotalVI integration
adata = adata[adata.obs['Status'].isin(['Covid', 'Healthy'])]
protein_adata = adata[:,adata.var['feature_types']=='Antibody Capture'].copy()
rna_adata = adata[:,adata.var['feature_types']=='Gene Expression'].copy()

protein_adata.X = protein_adata.layers['raw'].A.copy()
rna_adata.X = rna_adata.layers['raw'].copy()
sc.pp.normalize_total(rna_adata, target_sum=1e4)
sc.pp.log1p(rna_adata)
rna_adata.obs_names_make_unique()

mdata = md.MuData({"rna": rna_adata, "protein": protein_adata})

sc.pp.highly_variable_genes(
    mdata.mod["rna"],
    n_top_genes=3000,
    flavor="cell_ranger",
    batch_key="sample_id",
)

# Place subsetted counts in a new modality
mdata.mod["rna_subset"] = mdata.mod["rna"][
    :, mdata.mod["rna"].var["highly_variable"]
].copy()

mdata.update()

scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="raw",
    protein_layer=None,
    batch_key="sample_id",
    modalities={
        "rna_layer": "rna_subset",
        "protein_layer": "protein",
        "batch_key": "rna_subset",
    },
)

vae = scvi.model.TOTALVI(mdata,override_missing_proteins=True)

vae.train()

fig, ax = plt.subplots(1, 1)
vae.history["elbo_train"].plot(ax=ax, label="train")
vae.history["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs", 
       ylim=(1000, 3500)
      )
ax.legend()

vae.save(f'{data_dir}totalVI_model_Covid_Healthy')

vae = scvi.model.TOTALVI.load(f'{data_dir}totalVI_model_Covid_Healthy/', adata=mdata)

rna = mdata.mod["rna_subset"]
protein = mdata.mod["protein"]

# arbitrarily store latent in rna modality
rna.obsm["X_totalVI"] = vae.get_latent_representation()

chunk_size = 10000
n_full_chunks = rna.n_obs // chunk_size
rem_chunk_size = rna.n_obs % chunk_size

if n_full_chunks > 0:
    index_ranges = np.arange(n_full_chunks + 1) * chunk_size
    if rem_chunk_size > 0:
        full_ranges = np.concatenate([index_ranges, np.array([rna.n_obs])]) 
    else:
        full_ranges = index_ranges    
else:
    full_ranges = np.arange(0, 2) * rem_chunk_size
    
#run in several batches to lower cpu and memory load
rna_list = []
protein_list = []
for idx in range(1, len(full_ranges)):
    print(idx)
    rna_denoised, protein_denoised = vae.get_normalized_expression(adata=mdata, 
                            indices=range(full_ranges[idx-1], full_ranges[idx]), #use only subset if kernel dies
        n_samples=25, return_mean=True, transform_batch=None, 
    )
    rna_list.append(rna_denoised)
    protein_list.append(protein_denoised)
    
rna_denoised_all = pd.concat(rna_list)

protein_denoised_all = pd.concat(protein_list)

(
    rna.layers["denoised_rna"],
    protein.layers["denoised_protein"],
) = (rna_denoised_all, protein_denoised_all)

protein.layers["protein_foreground_prob"] = vae.get_protein_foreground_probability(
    n_samples=25, return_mean=True, transform_batch=None
)
#parsed_protein_names = [p.split("_")[1] for p in protein.var_names]
#protein.var["clean_names"] = parsed_protein_names
mdata.update()

sc.pp.neighbors(rna, use_rep="X_totalVI")
sc.tl.umap(rna)
sc.tl.leiden(rna, key_added="leiden_totalVI")

mdata.update()

mdata.write(f'{data_dir}haniffa21_totalVI_Covid_Healthy.h5mu')