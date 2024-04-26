import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pkg_resources
import os
import scanpy as sc
import re
import warnings
from contextlib import suppress
import anndata
import datetime

def print_main_versions():
    sc.logging.print_header()

def save_package_versions(base_package_version_path,pre,do_print = True):
    #base_package_version_path: path to folder where to save all package versions
    #pre: unique script/notebook identifier, e.g. 'M10'
    with open(os.path.join(base_package_version_path, pre + '_package_versions.txt'), "w") as file:
        for package in pkg_resources.working_set:
            file.write(f"{package.key}=={package.version}\n")
            if do_print:
                print(f"{package.key}=={package.version}")

def return_base_package_version_path():
    path = os.path.join(os.getcwd(),'package_versions')
    if not os.path.exists(path):
                os.mkdir(path)
    return path

def return_input_data_path():
    path = os.path.join(os.getcwd(),'input_data')
    if not os.path.exists(path):
                os.mkdir(path)
    return path

def make_path(path):
    if not os.path.exists(path):
                os.mkdir(path)
    else:
        print('path exists')
    

def return_output_data_path():
    path = os.path.join(os.getcwd(),'output_data')
    if not os.path.exists(path):
                os.mkdir(path)
    return path


def get_time():
    return datetime.datetime.now()

def get_time_delta(start,end,save = False,pre=None):
    time_delta = end - start
    path = os.path.join(os.getcwd(),'time_consumption')
    if save:
        if not os.path.exists(path):
                os.mkdir(path)
        with open(os.path.join(path,pre + '_time_consumption.txt'), "w") as file:
            file.write(str(time_delta))
    return time_delta

def get_subsample_fraction(adata,cluster_string,target_cluster,subsample_factor):
    #adata: AnnData object
    #cluster_string: string in adata.obs
    #target_cluster: category in adata.obs['cluster_string']
    #subsample_factor: int, to obtain ratio of subsample_factor:1 for non_targets:targets 
    adata_target = adata[adata.obs[cluster_string] == target_cluster]
    adata_rest = adata[adata.obs[cluster_string] != target_cluster]
    sc.pp.subsample(adata_rest,n_obs = subsample_factor*len(adata_target))
    adata_sub = adata_rest.concatenate(adata_target)
    return adata_sub


def get_subsample_fraction_final_obs(adata,cluster_string,target_cluster,fraction,final_obs):
    adata_target = adata[adata.obs[cluster_string] == target_cluster]
    adata_rest = adata[adata.obs[cluster_string] != target_cluster]
    n_obs_targets = int(fraction*final_obs)
    n_obs_rest = int((1-fraction)*final_obs)
    sc.pp.subsample(adata_target,n_obs = n_obs_targets)
    sc.pp.subsample(adata_rest,n_obs = n_obs_rest)
    adata_sub = adata_rest.concatenate(adata_target)
    return adata_sub

def get_fraction(adata,cluster_string,cluster_name):
    total_all = len(adata)
    total_cluster = adata.obs[cluster_string].value_counts()[cluster_name]
    return total_cluster/total_all
   
    
    
