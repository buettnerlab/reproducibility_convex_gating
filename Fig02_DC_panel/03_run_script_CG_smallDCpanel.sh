#!/bin/bash

#source $HOME/.bashrc

# declare an array called array and define 3 values
filename="/groups/NovaSeq-01/bioinformatics/users/buettnerm/convex_gating/data/Hofer_IMA/mono_merge_annotated.h5ad"
clusterid="cell_types_lvl2"
outpath="/groups/NovaSeq-01/bioinformatics/users/buettnerm/convex_gating/data/performance_test/convex_gating/Hofer_data/${clusterid}/"
sampleid="sample"
array=("B cell" "classical monocyte" "non-classical monocyte" "intermediate monocyte" 
       "DC Q1(CD14-, CD5+)" "DC Q3(CD14+, CD5-)" "DC Q4(CD14-, CD5-)" )
#("B cell" "monocyte" "dendritic cell")
#("B cell" "classical monocyte" "non-classical monocyte" "intermediate monocyte" 
#      "myeloid-derived suppressor cell" "DC Q1(CD14-, CD5+)" "DC Q3(CD14+, CD5-)" "DC Q4(CD14-, CD5-)" )#("DC Q1(CD14-, CD5+)")
#
#("B cell" "monocyte" "dendritic cell") 
#("B cell" "T cell" "NK cell" "CD4+ T cell" "CD8+ T cell" "NK T cell" 
# "Classical monocyte" "Intermediate monocyte" "Non-classical monocyte" "cDC" "pDC" "pre-DC" ) 
#("B cell" "Dendritic cell" "Monocyte" "NK cell" "T cell" ) 
#("B cell" "monocyte" "dendritic cell") 
#("B cell" "classical monocyte" "non-classical monocyte" "intermediate monocyte" "dendritic cell Q1"        "dendritic cell Q3" "dendritic cell Q4" ) #"CD4+ T cell" "CD8+ T cell" "NK T cell" "Classical monocyte" "Non-classical monocyte" 

for i in "${array[@]}"
do
  /opt/python/bin/python run_script_cg.py -u $filename -o $outpath -k "$clusterid" -g "$i" -b $sampleid -s 50000 
done