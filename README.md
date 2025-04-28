# single_cell_JO
Repository for single cell analysis relating to a specific project.
main_analysis.R details the R process followed.

# For python users: 
process_sRNA.py has been included to demonstrate basic processing of scRNA data via scanpy for python users 
It includes a short set of commands initially to create and install pre-requisites in a conda environment 
Then moves onto pre-processing, removal of doublets, dimension reduction and clustering of scRNA data. 
process_scRNA_MT.py details the mitochondrial population checks but leaves the usual mitochondrial filter off 
so that all cells (originally provided and containing an off-shoot cluster) can be used.
process_scRNA_cell_cycle_validation.py details the cell cycle validation.
