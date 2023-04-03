# CITE-seq COVID19

In this example, we demonstrate the applicability of Convex Gating across modalities. In the notebooks, we first denoise the CITE-seq data using totalVI. Th
Then, we use the denoised protein data to define a CD16+ T cell type. Here we rely on the previously published cell type definition for CD4 and CD8 T cells. 
We denote CD16+ T cells as T cells with a foreground probability of CD16 higher than 0.5 and being clustered with the CD4 and CD8 T cell compartment, respectively.

## Define a target population

We aim to derive a gating strategy for CD16+ CD4+ and CD8+ T cells, respectively, in COVID-19 patients. 
We therefore subset the data to the COVID-19 samples. The target populations contain very few cells (less than 1,000 each). 
To obtain a more favorable ratio of target to non-target cells, we subsample the non-target population.

## Data normalization

In order to run Convex Gating, we have to normalize the protein expression data. In general, CITE-seq data is count-based with values ranging from 0 - ca. 2.000.
We will compare the impact of denoising on the performance of Convex Gating. 
Thus, we run Convex Gating on both denoised and raw protein data. In both settings, we scale each feature by log+1 and min-max. 

## Tasks

1. Derive a gating strategy for CD16+ T cells from CITE-seq data (both denoised and raw).
2. Compare gating strategy to gating strategy on cyTOF data: Did we use identical markers? Which markers are different? 
3. Subset feature spaces and apply gating strategies across data types: What is the performance of Convex Gating if we subset to the cyTOF panel? Vice versa: What is the performance of Convex Gating if we subset the cyTOF panel to CITEseq data?
4. Can we directly apply the gating strategy derived from CITE-seq and cyTOF to the corresponding other modality?    
