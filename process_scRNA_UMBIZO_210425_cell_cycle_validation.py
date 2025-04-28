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
import goatools as goat




import os, re

## Change location to where data has been stored 
os.chdir('/Path_to_sc_experiment/')

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
mt_threshold = 20
minimum_cell_per_gene = 3
min_count_threshold = 2

adata_filtered.var["mt"] = adata_filtered.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_filtered, qc_vars= ["mt"], inplace=True)

## MT and gene count filtering
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt <= mt_threshold,:]
sc.pp.filter_genes(adata_filtered, min_counts=minimum_cell_per_gene)
sc.pp.filter_cells(adata_filtered, min_genes=min_count_threshold)

## Normalise and log the data
sc.pp.normalize_total(adata_filtered, target_sum= 10000)
sc.pp.log1p(adata_filtered)

## Retain cells that have a score of atlas 0.5 in one of the cell cycles 
adata_filtered = adata_filtered[adata_filtered.obs[['G1-S', 'G2-M', 'M', 'M-G1', 'S']].max(axis = 1) > 0.5]
adata_filtered.obs['Score_Delta'] = adata_filtered.obs[['G1-S', 'G2-M', 'M', 'M-G1', 'S']].apply(lambda x:sorted(x,reverse = True), axis=1).apply(lambda y:y[0] - y[1])
adata_guide = adata_filtered[adata_filtered.obs['guide_identity'].apply(lambda x: re.search('_pos', x)).isna()] 

## Included a 'delta' score as another filter -> This is the difference between the current best-scoring cell cycle phase versus the next best -> 
## The greater the delta, the more confident in the phase being unambigious 
## Left as 0 for now (meaning no delta) but can be altered  
adata_delta = adata_guide[adata_guide.obs['Score_Delta'] > 0]

##### Doublet detection and removal 
scrub = scr.Scrublet(adata_delta.X, expected_doublet_rate=0.06, random_state= 0)

## The min_counts and min_cells is used for gene filtering
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                           min_cells=3,
                                                           min_gene_variability_pctl=85,
                                                           n_prin_comps=30)
adata_delta.obs["scr_score"] = doublet_scores
adata_delta.obs["scr_prediction"] = predicted_doublets

## Remove doublets
adata_delta = adata_delta[adata_delta.obs.scr_prediction == False, :]

## The data is now clean and can go through various steps of interest i.e. dimension reduction 
sc.pp.highly_variable_genes(adata_delta)
sc.pp.pca(adata_delta)
sc.pp.neighbors(adata_delta)
sc.tl.leiden(adata_delta)
sc.tl.umap(adata_delta)

## UMAP looks fairly nice overall -> There is separation of cells according to phases 
sc.pl.umap(adata_delta, color = 'cell_cycle_phase')

## Identification of cell markers that differentiate between phases  
sc.tl.rank_genes_groups(adata_delta, groupby="cell_cycle_phase", method = "wilcoxon")

def rank_gene_extract(in_adata, group_info, logfoldchanges = 0.58, pvals_adj = 0.001):
    out_df = sc.get.rank_genes_groups_df(in_adata, group = group_info)
    out_df = out_df[out_df['logfoldchanges'].abs() >= logfoldchanges]
    out_df = out_df[out_df['pvals_adj'] <= pvals_adj]
    return(out_df)

g1_s_df = rank_gene_extract(adata_delta, group_info = "G1-S", logfoldchanges = 0.58, pvals_adj= 0.001)
s_df = rank_gene_extract(adata_delta, group_info = "S", logfoldchanges = 0.58, pvals_adj= 0.001)
g2_m_df = rank_gene_extract(adata_delta, group_info = "G2-M", logfoldchanges = 0.58, pvals_adj= 0.001)
m_df = rank_gene_extract(adata_delta, group_info = "M", logfoldchanges = 0.58, pvals_adj= 0.001)
m_g1_df = rank_gene_extract(adata_delta, group_info = "M-G1", logfoldchanges = 0.58, pvals_adj= 0.001)

## Writing out results for GO and pathway analysis  
g1_s_df.to_csv("g1_s.csv", index = False)
s_df.to_csv("s.csv", index = False)
g2_m_df.to_csv("g2_m.csv", index = False)
m_df.to_csv("m_df.csv", index = False)
m_g1_df.to_csv("m_g1.csv", index = False)

g1_s_df['names'].to_csv("g1_s_gene.csv", index = False)
s_df['names'].to_csv("s_gene.csv", index = False)
g2_m_df['names'].to_csv("g2_m_gene.csv", index = False)
m_df['names'].to_csv("m_df_gene.csv", index = False)
m_g1_df['names'].to_csv("m_g1_gene.csv", index = False)

## Can export these genes to an enrichment server e.g. DAVID 
# https://davidbioinformatics.nih.gov/home.jsp


## Save the h5ad file 
adata_delta = sc.read_h5ad("adata_cell_cycle_validation.h5ad")

