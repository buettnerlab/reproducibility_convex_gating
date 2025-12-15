import numpy as np
import os
import pandas as pd
from shapely.geometry import Polygon, Point
import anndata as ann
from sklearn.metrics import precision_score, recall_score, f1_score

def adata_to_df_gating(adata, cluster_string):
    """
    transform adata such that in right form for gating algorithm

    Parameters
    ----------
    adata : anndata object
        annotated data frame as data input.
    cluster_string : str
        needs to be a column name of the adata.obs metadata.

    Returns
    -------
    cell_data : pandas DataFrame
        prepared DataFrame as input for the gating algorithm.

    """

    channels = list(adata.var.index)
    if type(adata.X) != np.ndarray:
        cell_data = pd.DataFrame(adata.X.todense(), columns=channels)
        cell_data[cluster_string] = adata.obs[cluster_string].values
    else:
        cell_data = pd.DataFrame(adata.X, columns=channels)
        cell_data[cluster_string] = adata.obs[cluster_string].values
    return cell_data

def apply_gating_strategy(adata,base_path,cluster,cluster_string,output_scores =True):
    #adata: anndata object 
    #base_path: output path of convexgating (directory with 'marker_summary.csv', 'performance_summary.csv' and  the subdirectory with the gate edge coordinates per hierarchy)
    #cluster: cluster number
    #cluster_string: column in adata.obs where cluster annotation
    cell_data = adata_to_df_gating(adata, cluster_string=cluster_string)
    cluster = str(cluster)
    marker_summary_base = pd.read_csv(os.path.join(os.path.join(base_path,'marker_summary.csv')))
    performance_summary_base = pd.read_csv(os.path.join(os.path.join(base_path,'performance_summary.csv')))
    n_hierarchies = performance_summary_base[performance_summary_base['cluster'].astype(str) == cluster]['hierarchy'].values[0]
    for hierarchy in range(1,n_hierarchies+1):
        base_markers = marker_summary_base[
                    (marker_summary_base['hierarchy'] == hierarchy) &
                    (marker_summary_base['cluster'].astype(str) == cluster)
                ]['marker'].values
        coordinates_base = pd.read_csv(os.path.join(base_path,f'cluster_{cluster}', f'cluster_{cluster}_gate_edges_hierarchy_{hierarchy}.csv'),
                                       index_col=0)
        marker1_base, marker2_base = base_markers
        coordinates_base.columns = [marker1_base, marker2_base]
        poly_base = Polygon(coordinates_base.values)
        is_in_gate = []
        for x, y in zip(cell_data[base_markers[0]], cell_data[base_markers[1]]):
            point = Point(x, y)
            is_in_gate.append(poly_base.covers(point))
        cell_data['hierarchy_' +str(hierarchy) +'_in_gate'] = is_in_gate
    hierarchy_cols = [f'hierarchy_{h}_in_gate' for h in range(1, n_hierarchies + 1)]
    cell_data['gated'] = cell_data[hierarchy_cols].all(axis=1)
    cell_data['target'] = cell_data[cluster_string] ==cluster
    f1 = f1_score(cell_data['target'], cell_data['gated'])
    recall = recall_score(cell_data['target'], cell_data['gated'])
    precision = precision_score(cell_data['target'], cell_data['gated'])
    if output_scores:
        return f1,recall,precision
    else:
        return cell_data,f1,recall,precision

#example usage
#f1,recall,precision = apply_gating_strategy(adata=adata_to_test,base_path = save_path,cluster = '1',cluster_string='KMeans')