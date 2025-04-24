
## Bespoke functions 

# Customised version of plotPCA function from DESeq2 package 
plotPCA_DESeq2.local <- function(object, intgroup="condition", inputGenes = NA, ntop = 500, returnData=FALSE)
{
  # calculate the variance for each gene
  if(is.na(inputGenes) == TRUE){
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  } 
  
  ## Added this conditional here 
  if(is.na(inputGenes[1]) == FALSE){
  select <- inputGenes
  }
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }

  gplot.out <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", label = "name")) +  geom_point(size = 4) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        coord_fixed()
                
return(gplot.out)               
}


## Takes in human or mouse identifiers and generates a data frame with ensembl, entrezID and symbol 
## Input -> 
# in.gene: character vector of identifier (can be ensembl, entrez or symbol)
# in.geneSource: character stating what format the in.gene input is in -> Default is ensembl 
# in.org -> R organism database -> Works for human (org.Hs.eg.db) or mouse (org.Mm.eg.db) -> org.Hs.eg.db by default 
# output -> 3 column table with ensembl, entrez and symbol 

conversion.function.clusterProfiler <- function(in.gene, in.geneSource = "ENSEMBL", in.org = "org.Hs.eg.db"){
gene.sources <- c("ENSEMBL", "ENTREZID", "SYMBOL")
input.source <- in.geneSource 
out.source <- gene.sources[which(input.source != gene.sources)]
id.table <- clusterProfiler::bitr(in.gene, fromType = input.source, toType = out.source, OrgDb = in.org, drop = FALSE)
id.table <- id.table[,gene.sources]
return(id.table)
}

## Over-representation analysis for GO and KEGG 
# Input is the 3-column table from conversion.function.clusterProfiler()
# in.universe -> The background from which to perform GO analysis against 
# The 'universe' by default is all available genes in a given genome from the input organism database 
# If user wants to specify; they can put in a 3-column table from conversion.function.clusterProfiler() of background 
# Conventionally, a universe is either the total number of genes in the organism or the total number of genes detected in 
# the experiment (not just DEGs). 
# 
# in.org -> R organism database -> Default to 'org.Hs.eg.db' for human 
# in.species -> Required for the KEGG portion of code -> Either "HUMAN" or "MOUSE" 
# code.name -> A pre-fix for the output files you wish to add 

overrepresentation.analysis <- function(in.set, in.universe = NA, in.org = "org.Hs.eg.db", in.species = "HUMAN", code.name){
## Same as set.enrichment.analysis v1 but universe is set to default 
## GO 
gene <- unique(na.omit(in.set[,"ENTREZID"]))
gene.sym <- unique(na.omit(in.set[,"SYMBOL"]))
gene.ensembl <- unique(na.omit(in.set[,"ENSEMBL"]))

if(is.na(in.universe) == FALSE){
universe.genes <- unique(na.omit(in.universe[,"ENSEMBL"]))
ego <- clusterProfiler::enrichGO(gene          = gene.ensembl,
                OrgDb         = in.org,
				universe = universe.genes,
				keyType = "ENSEMBL",
                pAdjustMethod = "BH",
				ont = "ALL",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)

}
else{
## Setting the ID to ENSEMBL 

ego <- clusterProfiler::enrichGO(gene          = gene.ensembl,
                OrgDb         = in.org,
				keyType = "ENSEMBL",
                pAdjustMethod = "BH",
				ont = "ALL",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
}

ego.df <- data.frame(ego)

## KEGG needs entrez format so these are converted 
## KEGG analysis 
if(in.species == "HUMAN"){

kk <- clusterProfiler::enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
}
if(in.species == "MOUSE"){
kk <- clusterProfiler::enrichKEGG(gene         = gene,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
}				 
				 
kk.df <- data.frame(kk)
				 
write.table(ego, paste(code.name, "_GO_Enrichment.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(kk, paste(code.name, "_KEGG_Enrichment.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(in.set, paste(code.name, "_GeneTable_Enrichment.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
writeLines(gene.sym, paste(code.name, "_GeneSymbol.txt", sep = ""))

## Including plots 
ego.bar <- barplot(ego, showCategory = 15) + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10))
ego.dot <- dotplot(ego, showCategory = 15) + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10))
kk.bar <- barplot(kk, showCategory = 15) + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10))
kk.dot <- dotplot(kk, showCategory = 15) + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10))

pdf(paste(code.name, "_plots.pdf", sep = ""), width = 15)
plot(ego.bar, cex.axis = 1.5)
plot(ego.dot)
plot(kk.bar)
plot(kk.dot)
dev.off()
}


## Function for reference mapping
#
# Seurat has a handy function which allows the user to more-easily annotate a dataset by taking an already well-defined dataset 
# and use it as a reference to which one can overlay a new dataset on top of -> 
# The principle here is that the cells in the new/query dataset that are transcriptomically more similar to the reference will 
# closely align/map onto the same reduced space. 
# This assumption allows for transferring of the reference annotations onto the new/query data 

# Input 
# reference.obj -> a seurat object representing the reference 
# query.obj -> a seurat object representing the query/new dataset 
#
# It needs to be noted that the reference and query objects should be processed and normalised in the same way (in this case 'SCT' normalisation) 
# SCT (single-cell transformation) is the recommended approach from seurat 
#
# reference.name -> The name of mapping annotations in the reference that you want to transfer over to the 
# query -> This is a column name in the reference meta data and defaults to "seurat_clusters"
# Different references/atlases have been produced with different pipelines/labs/people and naming conventions are 
# therefore inconsistent
# As such, the seurat default column name 'seurat_clusters' is used as default to make sure that this function runs. 
#
# sample.name -> Pre-fix for the output files (figures and annotated query)
reference.mapping.seurat <- function(reference.obj, query.obj, reference.name = "seurat_clusters", sample.name)
{

## Auto check the dimension of PCA reduction slot 

#ref.dim <- ncol(reference.obj@reductions$pca)
#query.dim <- ncol(query.obj@reductions$pca)

min.dim <- 30

# min.dim <- min(c(ref.dim, query.dim))

if(min.dim >= 30){
dim.param <- 30 
}
else{
dim.param <- min.dim
}

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference.obj,
  query = query.obj,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:dim.param
)

# set the cell identities to the cell type predictions
reference.obj <- RunUMAP(reference.obj, dims = 1:dim.param, reduction = "pca", return.model = TRUE)

query.obj <- MapQuery(anchorset = transfer_anchors, reference = reference.obj, query = query.obj,
    refdata = list(seurat_clusters = reference.name), reference.reduction = "pca", reduction.model = "umap", transferdata = list(prediction.assay = TRUE))

## Prediction information is now in: 
# query.obj@assays$prediction.score.seurat_clusters@data

prediction.scores <- t(query.obj@assays$prediction.score.seurat_clusters@data)
colnames(prediction.scores) <- paste0("predicted.", colnames(prediction.scores), ".scores")
query.obj@meta.data <- cbind(query.obj@meta.data, prediction.scores)

## Plot next to each other 
p1 <- DimPlot(reference.obj, reduction = "umap", group.by = reference.name, label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query.obj, reduction = "ref.umap", group.by = "predicted.seurat_clusters", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p3 <- p1 + p2
ggsave(paste0(sample.name, "_Reference_Mapped.pdf"), p3, width = 15, height = 10)

## Added to save object 
base::saveRDS(query.obj, paste(sample.name, "_seuratObject.rds", sep = "", collapse = ""))

return(query.obj)
}


package.dir <- "/Path/to/packages"
.libPaths(new = package.dir)

work.dir <- "/Path/to/workingDirectory"
setwd(work.dir)

library(DESeq2)
library(Seurat)
library(data.table)
library(clusterProfiler)

## This option is needed to increase RAM allowed to be used 
## In this case -> 30 GB 
options(future.globals.maxSize= 30.0 * 1e9)

dir.create("obj1")
setwd("obj1")


## Read in and generate seurat object from counts and meta data 
# analysis1 <- function(){
# dir.create("analysis1")
# setwd("analysis1")

# count.mat <- data.frame(fread("../../matrix.csv"))
# meta.info <- data.frame(fread("../../count_mat/obs.csv"))

# gct.mat <- count.mat[,2:ncol(count.mat)]
# rownames(gct.mat) <- count.mat[,1]

# ## Change back to gene names 
# gct.mat <- ConvertEnsembleToSymbol.local(gct.mat, species = "human")


# s.obj <- CreateSeuratObject(counts = gct.mat, project = "UM_04", meta.data = meta.info)
# s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^MT-")

# in.cells <- s.obj 

# ## non-targeting
# ## Targeting 

# in.cells <- SCTransform(in.cells, verbose = FALSE) %>% RunPCA() %>% 
# RunUMAP(dims = 1:30, reduction.name = 'umap', reduction.key = 'rnaUMAP_') %>%
# RunTSNE(dims = 1:30, method = "FIt-SNE") %>% FindNeighbors(reduction = "pca", dims = 1:30) %>% 
# FindClusters(resolution = 1)

# saveRDS(in.cells, "seuratObject.rds")

# }

## Generate seurat object for targeted cells 
analysis1b <- function(){
dir.create("analysis1b")
setwd("analysis1b")

count.mat <- data.frame(fread("../../matrix.csv"))
meta.info <- data.frame(fread("../../count_mat/obs.csv"))

gct.mat <- count.mat[,2:ncol(count.mat)]
rownames(gct.mat) <- count.mat[,1]

## Change back to gene names 
gct.mat <- ConvertEnsembleToSymbol.local(gct.mat, species = "human")


s.obj <- CreateSeuratObject(counts = gct.mat, project = "UM_04", meta.data = meta.info)
s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^MT-")

in.cells <- s.obj 

in.cells <- subset(in.cells, gene != "non-targeting")

## non-targeting
## Targeting 

in.cells <- SCTransform(in.cells, verbose = FALSE) %>% RunPCA() %>% 
RunUMAP(dims = 1:30, reduction.name = 'umap', reduction.key = 'rnaUMAP_') %>%
RunTSNE(dims = 1:30, method = "FIt-SNE") %>% FindNeighbors(reduction = "pca", dims = 1:30) %>% 
FindClusters(resolution = 1)

saveRDS(in.cells, "seuratObject.rds")

## Cluster 17 appears to be MT related 
## Cluster 21 -> Enriched for one guide type? 
## Cluster 19, 10 and 15 -> Enriched for another guide type? 


}

## Generate seurat object for non-targeted/control cells 
analysis1c <- function(){
dir.create("analysis1c")
setwd("analysis1c")

count.mat <- data.frame(fread("../../matrix.csv"))
meta.info <- data.frame(fread("../../count_mat/obs.csv"))

gct.mat <- count.mat[,2:ncol(count.mat)]
rownames(gct.mat) <- count.mat[,1]

## Change back to gene names 
gct.mat <- ConvertEnsembleToSymbol.local(gct.mat, species = "human")


s.obj <- CreateSeuratObject(counts = gct.mat, project = "UM_04", meta.data = meta.info)
s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^MT-")

in.cells <- s.obj 

in.cells <- subset(in.cells, gene == "non-targeting")

## non-targeting
## Targeting 

in.cells <- SCTransform(in.cells, verbose = FALSE) %>% RunPCA() %>% 
RunUMAP(dims = 1:30, reduction.name = 'umap', reduction.key = 'rnaUMAP_') %>%
RunTSNE(dims = 1:30, method = "FIt-SNE") %>% FindNeighbors(reduction = "pca", dims = 1:30) %>% 
FindClusters(resolution = 1)

saveRDS(in.cells, "seuratObject.rds")

}


analysis2 <- function(){
dir.create("analysis2")
setwd("analysis2")

meta.info <- data.frame(fread("../../count_mat/obs.csv"))

meta.info.control <- meta.info[which(meta.info[,"gene"] == "non-targeting"),]
meta.info.nc <- meta.info[which(meta.info[,"gene"] != "non-targeting"),]

meta.info.nc[,"guide_pos"] <- "N"
meta.info.nc[grep("_posA", fixed = TRUE, meta.info.nc[,"guide_identity"]),"guide_pos"] <- "A"
meta.info.nc[grep("_posB", fixed = TRUE, meta.info.nc[,"guide_identity"]),"guide_pos"] <- "B"

## Do genes have multiple guides? 

meta.info.nc.N <- data.frame(setorder(data.table(meta.info.nc[which(meta.info.nc[,"guide_pos"] == "N"),]), gene, guide_identity))
gene.table <- unique(meta.info.nc.N[,c("gene_id", "gene_transcript", "guide_identity")])
## one guide per gene 

}

## Pseudo-bulk assessment of data
## Targeting only and using each individual guide_identity set of cells as way to group cells (assuming more than 20 cells available) 
## Split according to cell phase and then draw PCA 
## 
## Pseudo-bulking is a process by which we convert single-cell RNA data into a conventional 'bulk' dataset -> This alleviates the data-sparsity 
## in single cell (at the cost of single-cell granularity) and allows the option to use well-established bulk-RNAseq algorithms (like DESeq2)

analysis3b <- function(){
dir.create("analysis3b")
setwd("analysis3b")

count.mat <- data.frame(fread("../../matrix.csv"))
meta.info <- data.frame(fread("../../count_mat/obs.csv"))

gct.mat <- count.mat[,2:ncol(count.mat)]
rownames(gct.mat) <- count.mat[,1]

## Change back to gene names 
gct.mat <- ConvertEnsembleToSymbol.local(gct.mat, species = "human")

s.obj <- CreateSeuratObject(counts = gct.mat, project = "UM_04", meta.data = meta.info)
s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^MT-")

in.cells <- s.obj 
in.cells <- subset(in.cells, gene != "non-targeting")

## Remove the posA and posB as they are positive controls
cell.phase <- unique(in.cells[[]][, "cell_cycle_phase"])

## Use the cell phase score as well -> Retain only those that have a high/confident score 
# Normal distribution between -1 and 1 -> Choosing 0.5 for now 
for(i in 1:length(cell.phase)){

tmp.i <- subset(in.cells, cell_cycle_phase == cell.phase[i])
s.meta <- tmp.i[[]]
s.meta[,"guide_pos"] <- "N"
s.meta[grep("_posA", fixed = TRUE, s.meta[,"guide_identity"]),"guide_pos"] <- "A"
s.meta[grep("_posB", fixed = TRUE, s.meta[,"guide_identity"]),"guide_pos"] <- "B"
int.meta <- s.meta[which(s.meta[,"guide_pos"] == "N"),]

## Subset based on cell phase score (0.5)
int.meta <- int.meta[which(int.meta[, gsub("-", ".", fixed = TRUE, cell.phase[i])] >= 0.5),]

id.table <- table(int.meta[,"guide_identity"])

# Min of 20 cells for pseudo-bulking 
int.guides <- names(id.table)[which(id.table >= 20)]
sub.int.meta <- int.meta[which(is.na(match(int.meta[,"guide_identity"], int.guides)) == FALSE),]
s.subset <- subset(in.cells, cells = rownames(sub.int.meta))

s.list <- SplitObject(s.subset, split.by = "guide_identity")

agg.df <- do.call(cbind, lapply(s.list, function(x){
tmp.x <- AggregateExpression(x)[["RNA"]]
return(tmp.x)
}))
colnames(agg.df) <- names(s.list)
colnames(agg.df) <- apply(as.matrix(colnames(agg.df)),1,function(x){strsplit(x,split = "_",fixed = TRUE)[[1]][[1]]})

saveRDS(agg.df, paste(cell.phase[i], "_agg.rds", sep = "", collapse = ""))

rm(tmp.i, s.meta, int.meta, id.table, int.guides, sub.int.meta, s.subset, s.list, agg.df)
}
rm(i) 

# "G1-S_agg.rds" "G2-M_agg.rds" "M-G1_agg.rds" "M_agg.rds"    "S_agg.rds" 
g1.s.agg <- readRDS("G1-S_agg.rds")
g2.m.agg <- readRDS("G2-M_agg.rds")
m.g1.agg <- readRDS("M-G1_agg.rds")
m.agg <- readRDS("M_agg.rds")
s.agg <- readRDS("S_agg.rds")

colnames(g1.s.agg) <- paste(colnames(g1.s.agg), "__G1-S", sep = "")
colnames(g2.m.agg) <- paste(colnames(g2.m.agg), "__G2-M", sep = "")
colnames(m.g1.agg) <- paste(colnames(m.g1.agg), "__M-G1", sep = "")
colnames(m.agg) <- paste(colnames(m.agg), "__M", sep = "")
colnames(s.agg) <- paste(colnames(s.agg), "__S", sep = "")

agg.df <- cbind(
g1.s.agg, 
g2.m.agg,
m.g1.agg, 
m.agg, 
s.agg
)

gct.mat <- agg.df 
gct.design <- data.frame(
SampleID = colnames(gct.mat), 
SampleType = apply(as.matrix(colnames(gct.mat)),1,function(x){strsplit(x,split = "__",fixed = TRUE)[[1]][[1]]}),
Phase = apply(as.matrix(colnames(gct.mat)),1,function(x){strsplit(x,split = "__",fixed = TRUE)[[1]][[2]]})
)

dds.obj <- DESeq(DESeqDataSetFromMatrix(
countData = gct.mat,
colData = gct.design,
design = ~1
))

## Using variance stabilisation transformation vst() for PCA 
g1 <- plotPCA_DESeq2.local(vst(dds.obj), intgroup = "SampleType")
g2 <- plotPCA_DESeq2.local(vst(dds.obj), intgroup = "Phase")

}


analysis4 <- function(){
dir.create("analysis4")
setwd("analysis4")

control.obj <- readRDS("../analysis1c/seuratObject.rds")

## CLean up the control object 
control.meta <- control.obj@meta.data 
control.cell.phase <- gsub("-", ".", unique(control.meta[, "cell_cycle_phase"]))
control.meta[,"top_cycle_score"] <- apply(control.meta[,control.cell.phase],1,max)
control.meta[,"guide_pos"] <- "N"
control.meta[grep("_posA", fixed = TRUE, control.meta[,"guide_identity"]),"guide_pos"] <- "A"
control.meta[grep("_posB", fixed = TRUE, control.meta[,"guide_identity"]),"guide_pos"] <- "B"

control.obj[["top_cycle_score"]] <- control.meta[,"top_cycle_score"]
control.obj[["Guide"]] <- control.meta[,"guide_pos"]

control.s.obj <- subset(control.obj, top_cycle_score >= 0.5)
control.s.obj <- subset(control.s.obj, Guide == "N")

control.s.obj <- CreateSeuratObject(counts = GetAssayData(control.s.obj, assay = "RNA", slot = "counts"), meta.data = control.s.obj[[]])

c.s.obj <- SCTransform(control.s.obj, verbose = FALSE) %>% RunPCA() %>% 
RunUMAP(dims = 1:30, reduction.name = 'umap', reduction.key = 'rnaUMAP_') %>%
RunTSNE(dims = 1:30, method = "FIt-SNE") %>% FindNeighbors(reduction = "pca", dims = 1:30) %>% 
FindClusters(resolution = 1)


target.obj <- readRDS("../analysis1b/seuratObject.rds")

## Clean up the target list (retain only those cells where the cell phase score > 0.5), split and then 
## perform projection mapping 

target.meta <- target.obj@meta.data 
cell.phase <- gsub("-", ".", unique(target.meta[, "cell_cycle_phase"]))
target.meta[,"top_cycle_score"] <- apply(target.meta[,cell.phase],1,max)
target.meta[,"guide_pos"] <- "N"
target.meta[grep("_posA", fixed = TRUE, target.meta[,"guide_identity"]),"guide_pos"] <- "A"
target.meta[grep("_posB", fixed = TRUE, target.meta[,"guide_identity"]),"guide_pos"] <- "B"

target.obj[["top_cycle_score"]] <- target.meta[,"top_cycle_score"]
target.obj[["Guide"]] <- target.meta[,"guide_pos"]

target.s.obj <- subset(target.obj, top_cycle_score >= 0.5)
target.s.obj <- subset(target.s.obj, Guide == "N")

target.mapped <- reference.mapping.seurat(c.s.obj, target.s.obj, reference.name = "cell_cycle_phase", sample.name = "vs_target")

base::saveRDS(c.s.obj, "control_reference.rds")
base::saveRDS(target.mapped, "target_mapped.rds")

}

## GO and KEGG analysis 
## Performing over-representation analysis on the data
## This is done in context of the high mitochondrial population (cluster 17) which is high MT% vs other cells 
## and cell-cycle phase 
## The MT is to check whether certain pathways/GO terms are enriched for this cluster
## The cell-cycle analysis is performed for validation purposes  

analysis5 <- function(){
dir.create("analysis5")
setwd("analysis5")

control.obj <- readRDS("../analysis4/control_reference.rds")
target.obj <- readRDS("../analysis4/target_mapped.rds")

## The mito population and library size are invertly proportional to each other 
# DimPlot(target.obj, reduction = "ref.umap", label = TRUE)
# cluster 17 is the MT cluster of interest 

target.meta <- target.obj[[]]
target.meta[,"MT_Status"] <- NA 
target.meta[which(target.meta[,"seurat_clusters"] == "17"), "MT_Status"] <- "MT_Interest"
target.meta[which(target.meta[,"seurat_clusters"] != "17"), "MT_Status"] <- "nMT_Interest"

target.obj[["MT_Status"]] <- target.meta[,"MT_Status"]

DimPlot(target.obj, group.by = "MT_Status", label = TRUE)

Idents(target.obj) <- "MT_Status"
levels(Idents(target.obj)) <- c("MT_Interest", "nMT_Interest")
marker.set <- FindAllMarkers(target.obj, only.pos = TRUE)
cc.marker.set <- FindAllMarkers(target.obj, only.pos = TRUE, group.by = "cell_cycle_phase")

marker.set <- marker.set[which(marker.set[,"p_val_adj"] <= 0.001),]
cc.marker.set <- cc.marker.set[which(cc.marker.set[,"p_val_adj"] <= 0.001),]

library(ggplot2)
mt.converted <- conversion.function.clusterProfiler(marker.set[which(marker.set[,"avg_log2FC"] > 0.5),"gene"], "SYMBOL", "org.Hs.eg.db")
overrepresentation.analysis(mt.converted, in.universe = NA, "org.Hs.eg.db", "HUMAN", "MT_Comparison")

cc.converted <- conversion.function.clusterProfiler(cc.marker.set[which(cc.marker.set[,"avg_log2FC"] > 0.5),"gene"], "SYMBOL", "org.Hs.eg.db")
overrepresentation.analysis(cc.converted, in.universe = NA, "org.Hs.eg.db", "HUMAN", "CC_Comparison")

}


## Sanity-checking 
## To ensure that the high MT% cluster matches the bimodal cluster (high MT% from clients analysis); 
## A sanity-check is done where the annotation of this cluster of cells was done blind (as is seen above in analysis1b and analysis1c) 
## and the same UMAP structure is re-coloured according to the client's bimodal annotated cells 
## To ensure that they are the same cells/cluster 

## In addition, further QC checks were done on this high MT% cluster
## Specifically, looking at number of genes identified in the high MT% cluster vs other clusters
## Test here is that dying cells will have high MT% content but also low genes detected -> principle here being that a dying cell will 
## leak out their transcript -> leading to loss of detected genes 

analysis6 <- function(){
dir.create("analysis6")
setwd("analysis6")

factor.vector <- data.frame(fread("../../factor_output.csv"))
mt.vector <- data.frame(fread("../../mt_output.csv"))

fact.vec <- factor.vector[2:nrow(factor.vector),2]
mt.vec <- mt.vector[,2]

names(fact.vec) <- factor.vector[2:nrow(factor.vector),1]
names(mt.vec) <- mt.vector[,1]

bimodal.df <- cbind(fact.vec, mt.vec[names(fact.vec)])
colnames(bimodal.df) <- c("Factor", "Mitochondria")
rownames(bimodal.df) <- gsub("-", ".", fixed = TRUE,rownames(bimodal.df))

int.set <- bimodal.df[which(bimodal.df[,"Mitochondria"] >= 0.2),]

target.obj <- readRDS("../analysis4/target_mapped.rds")


## Annotation based on UMBIZO clustering (cluster 17)
target.meta <- target.obj[[]]
target.meta[,"MT_Status"] <- NA 
target.meta[which(target.meta[,"seurat_clusters"] == "17"), "MT_Status"] <- "MT_Interest"
target.meta[which(target.meta[,"seurat_clusters"] != "17"), "MT_Status"] <- "nMT_Interest"

## Annotation based on bimodal annotation 
target.meta[,"MT_Factor"] <- "nInt" 
int.set.subset <- int.set[which(is.na(match(rownames(int.set), rownames(target.meta))) == FALSE),]
target.meta[rownames(int.set.subset),"MT_Factor"] <- "Interest"

target.obj[["MT_Status"]] <- target.meta[,"MT_Status"]
target.obj[["MT_Factor"]] <- target.meta[,"MT_Factor"]

rm(target.meta)
target.meta <- target.obj[[]]

library(ggplot2)
d1 <- DimPlot(target.obj, group.by = "MT_Status", label = TRUE) + ggtitle("UMBIZO Clustering")
d2 <- DimPlot(target.obj, group.by = "MT_Factor", label = TRUE) + ggtitle("Bimodal Annotation")

d3 <- d1 + d2

## Appears to be the same cluster 


## Additional QC testing 
pdf("Preliminary_QC_Plots.pdf", width = 14)
boxplot(nFeature_RNA ~ MT_Factor, target.meta, outline = FALSE, xlab = "Cluster", ylab = "Number of genes", main = "Gene number discrepancies for bimodal annotated cells")
boxplot(nFeature_RNA ~ seurat_clusters, target.meta, outline = FALSE, xlab = "Cluster", ylab = "Number of genes", main = "Gene number discrepancies for UMBIZO clusters")
boxplot(nCount_RNA ~ MT_Factor, target.meta, outline = FALSE, xlab = "Cluster", ylab = "Number of genes", main = "Count distributions for bimodal annotated cells")
boxplot(nCount_RNA ~ seurat_clusters, target.meta, outline = FALSE, xlab = "Cluster", ylab = "Number of genes", main = "Count distributions for UMBIZO clusters")
plot(d3)
dev.off()

## High MT% cluster has fewer genes compared to other clusters 
## Suggestive that this cluster is mainly composed of dying cells

}

