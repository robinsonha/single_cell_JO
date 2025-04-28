## Update of scanpy_data_prep_QC_10x_h5

## Preparing conda environment 
## This is done on commandline ie Ubuntu and assumes that you have conda ready 

## Create conda environment and install necessary pre-requisites 
conda create -n scRNA 

## Activate and install pre-requisites 
conda activate scRNA 

conda install -c conda-forge scanpy python-igraph leidenalg
conda install pandas numpy

sudo apt-get update
sudo apt-get install g++

pip install scrublet 
pip install scvi-tools
pip install pybiomart

## More advanced tools 
# DEG analysis 
pip install pydeseq2

# GO analysis 
pip install goatools

python 

################################### Now run in Python 

## Import libraries  
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import seaborn as sns
import anndata as ad
import pickle
from sklearn.preprocessing import StandardScaler
import scanpy as sc
import scrublet as scr
import scvi
import os, re

## Change location to where data has been stored 
# os.chdir('/Path_to_sc_experiment/')

adata = ad.read_h5ad("sc_experiment.h5ad")

## The data uses ensembl IDs as gene features -> Switch to symbol instead 
annot = sc.queries.biomart_annotations(
    "hsapiens",
    ["ensembl_gene_id", "hgnc_symbol"],
).set_index("ensembl_gene_id")

## Remove anything that is NAN 
annot = annot.dropna()

## Retain only ensembl IDs that are consensus in both 
## Some of the annot has duplicate enmsebl IDs 
## Remove them 

annot_ids = pd.DataFrame(annot.index.to_list()).value_counts()
annot_pass = annot_ids[annot_ids == 1]

gene_table = pd.concat([pd.DataFrame(annot_pass.index.to_list()), pd.DataFrame(adata.var_names.to_list())], axis = 0, ignore_index=True).value_counts()
gene_table = gene_table.reset_index()
gene_table = gene_table.rename(columns = {0:'ENSEMBL'})
pass_genes = gene_table.loc[gene_table['count'] == 2,'ENSEMBL']

annot = annot.loc[pass_genes,]
adata_filtered = adata[:, adata.var_names.isin(pass_genes)]
adata_filtered.var[annot.columns] = annot.loc[adata_filtered.var_names.tolist()]
adata_filtered.var['Original_ENSEMBL'] = adata_filtered.var_names
adata_filtered.var_names = adata_filtered.var['hgnc_symbol']


## Now process 

## This MT threshold is non-standard! 
## Left at 100% to retrieve all cells in the original h5ad file
mt_threshold = 100

minimum_cell_per_gene = 3
min_count_threshold = 2

adata_filtered.var["mt"] = adata_filtered.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_filtered, qc_vars= ["mt"], inplace=True)

## MT and gene count filtering
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt <= mt_threshold,:]
sc.pp.filter_genes(adata_filtered, min_counts=minimum_cell_per_gene)
sc.pp.filter_cells(adata_filtered, min_genes=min_count_threshold)

## Remove guides 
adata_filtered = adata_filtered[adata_filtered.obs['guide_identity'].apply(lambda x: re.search('_pos', x)).isna()] 

## Normalise and log the data
sc.pp.normalize_total(adata_filtered, target_sum= 10000)
sc.pp.log1p(adata_filtered)


##### Doublet detection and removal 
scrub = scr.Scrublet(adata_filtered.X, expected_doublet_rate=0.06, random_state= 0)

## The min_counts and min_cells is used for gene filtering
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                           min_cells=3,
                                                           min_gene_variability_pctl=85,
                                                           n_prin_comps=30)
adata_filtered.obs["scr_score"] = doublet_scores
adata_filtered.obs["scr_prediction"] = predicted_doublets

## Remove doublets
adata_filtered = adata_filtered[adata_filtered.obs.scr_prediction == False, :]

## The data is now clean and can go through various steps of interest i.e. dimension reduction 
sc.pp.highly_variable_genes(adata_filtered)
sc.pp.pca(adata_filtered)
sc.pp.neighbors(adata_filtered)
sc.tl.leiden(adata_filtered)
sc.tl.umap(adata_filtered)

## Visualise MT across different clusters 

## Replicates analysis6 in R script (looking at MT% and genes detected%)
## Can see a clear cluster where cells have high MT% and low genes detected 
sc.pl.violin(
    adata_filtered,
    ["n_genes_by_counts", "mito_frac"],
    groupby="leiden",
    stripplot=False,  # remove the internal dots
    inner="box",  # adds a boxplot inside violins
)

## Save h5ad output 
adata_filtered.write_h5ad("sc_experiment_mt_check.h5ad")

