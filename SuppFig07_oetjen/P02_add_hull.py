import convex_hull_add_on as ch
import numpy as np
import convexgating as cg
import scanpy as sc
import anndata as ann
import pandas as pd
import os

#level2
os.mkdir('P02_cell_type_lvl2_hull')
samples = ['A','B','C','H','J','O','T','U']
for sample in samples:
    meta_info_path = os.path.join(os.getcwd(),'P01_cell_type_lvl2',sample,'meta_info.npy')
    #os.mkdir(os.path.join(os.getcwd(),'P02_cell_type_lvl2_hull',sample))
    target_location = os.path.join(os.getcwd(),'P02_cell_type_lvl2_hull',sample)
    ch.convex_hull_add_on(meta_info_path,target_location)

#level3
os.mkdir('P02_cell_type_lvl3_hull')
samples = ['A','B','C','H','J','O','T','U']
for sample in samples:
    meta_info_path = os.path.join(os.getcwd(),'P01_cell_type_lvl3',sample,'meta_info.npy')
    os.mkdir(os.path.join(os.getcwd(),'P02_cell_type_lvl3_hull',sample))
    target_location = os.path.join(os.getcwd(),'P02_cell_type_lvl3_hull',sample)
    ch.convex_hull_add_on(meta_info_path,target_location)

#level4
os.mkdir('P02_cell_type_lvl4_hull')
samples = ['A','B','C','H','J','O','T','U']
for sample in samples:
    meta_info_path = os.path.join(os.getcwd(),'P01_cell_type_lvl4',sample,'meta_info.npy')
    os.mkdir(os.path.join(os.getcwd(),'P02_cell_type_lvl4_hull',sample))
    target_location = os.path.join(os.getcwd(),'P02_cell_type_lvl4_hull',sample)
    ch.convex_hull_add_on(meta_info_path,target_location)

#level5
os.mkdir('P02_cell_type_lvl5_hull')
samples = ['A','B','C','H','J','O','T','U']
for sample in samples:
    meta_info_path = os.path.join(os.getcwd(),'P01_cell_type_lvl5',sample,'meta_info.npy')
    os.mkdir(os.path.join(os.getcwd(),'P02_cell_type_lvl5_hull',sample))
    target_location = os.path.join(os.getcwd(),'P02_cell_type_lvl5_hull',sample)
    ch.convex_hull_add_on(meta_info_path,target_location)