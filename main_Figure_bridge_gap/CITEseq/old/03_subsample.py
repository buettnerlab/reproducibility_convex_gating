import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import mudata as md
import muon as mu
import datetime
import anndata as ad

sc.logging.print_header()
file_path = '/work/users/mh823zote/projects/gating/data/CITEseq_Covid/' 
data_dir = file_path + 'data/' 
data_raw_dir = file_path + 'input/' 
sc.settings.figdir = file_path + 'figures/'
today = datetime.date.today().strftime('%y%m%d') #creates a YYMMDD string of today's date
pdata = md.read(f'{data_dir}haniffa21_totalVI_Covid_Healthy_annotated.h5mu/protein')

#CD16 CD4 - denoised
pdata_CD16_CD4_target = pdata[pdata.obs['CD16+ T cells'] == 'CD16+ CD4+ T cell']
pdata_CD16_CD4_non_target = pdata[pdata.obs['CD16+ T cells'] != 'CD16+ CD4+ T cell']
sc.pp.subsample(pdata_CD16_CD4_non_target ,n_obs = 15*len(pdata_CD16_CD4_target))
pdata_CD16_CD4_sub_denoised = pdata_CD16_CD4_non_target.concatenate(pdata_CD16_CD4_target)

#CD16 CD4 - raw
pdata_CD16_CD4_sub_raw = ad.AnnData(pdata_CD16_CD4_sub_denoised.layers['raw'], 
                                    obs=pdata_CD16_CD4_sub_denoised.obs,
                                    var=pdata_CD16_CD4_sub_denoised.var)


#CD16 CD8 - denoised
pdata_CD16_CD8_target = pdata[pdata.obs['CD16+ T cells'] == 'CD16+ CD8+ T cell']
pdata_CD16_CD8_non_target = pdata[pdata.obs['CD16+ T cells'] != 'CD16+ CD8+ T cell']
sc.pp.subsample(pdata_CD16_CD8_non_target ,n_obs = 15*len(pdata_CD16_CD8_target))
pdata_CD16_CD8_sub_denoised = pdata_CD16_CD8_non_target.concatenate(pdata_CD16_CD8_target)

#CD16 CD8 - raw
pdata_CD16_CD8_sub_raw = ad.AnnData(pdata_CD16_CD8_sub_denoised.layers['raw'], 
                                    obs=pdata_CD16_CD8_sub_denoised.obs,
                                    var=pdata_CD16_CD8_sub_denoised.var)


pdata_CD16_CD4_sub_denoised.write_h5ad(data_dir + 'pdata_CD16_CD4_sub_denoised.h5ad')
pdata_CD16_CD4_sub_raw.write_h5ad(data_dir + 'pdata_CD16_CD4_sub_raw.h5ad')
pdata_CD16_CD8_sub_denoised.write_h5ad(data_dir + 'pdata_CD16_CD8_sub_denoised.h5ad')
pdata_CD16_CD8_sub_raw.write_h5ad(data_dir + 'pdata_CD16_CD8_sub_raw.h5ad')