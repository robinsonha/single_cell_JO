## R installation -> Assuming 4.4.2 

## Rtools 
## 
# https://cran.r-project.org/bin/windows/Rtools/rtools44/files/rtools44-6459-6401.exe
# 
# Rtools needed for Seurat and SeuratObject 

local.dir <- "/Path/to/LocationToInstallPackages"
dir.create(local.dir)

bioc.version <- "3.20"
.libPaths(new = local.dir)
install.packages("BiocManager", lib = local.dir, repos = "https://www.stats.bris.ac.uk/R/")

install.packages("devtools",lib=local.dir, repos = "https://www.stats.bris.ac.uk/R/")
BiocManager::install("org.Hs.eg.db", lib=local.dir, version = bioc.version)

install.packages("tidyverse",lib=local.dir, repos = "https://www.stats.bris.ac.uk/R/")
install.packages("dplyr",lib=local.dir, repos = "https://www.stats.bris.ac.uk/R/")

install.packages("data.table",lib=local.dir, repos = "https://www.stats.bris.ac.uk/R/")
install.packages("Seurat",lib=local.dir, repos = "https://www.stats.bris.ac.uk/R/")
install.packages("SeuratObject",lib=local.dir, repos = "https://www.stats.bris.ac.uk/R/")

BiocManager::install("org.Hs.eg.db", lib=local.dir, version = bioc.version)
BiocManager::install("biomaRt", lib=local.dir, version = bioc.version)
BiocManager::install('DESeq2', lib=local.dir, version = bioc.version)
BiocManager::install("clusterProfiler", lib=local.dir, version = bioc.version)



