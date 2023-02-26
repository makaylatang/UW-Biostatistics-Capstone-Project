library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(InSituType)
library(clue)
library(pheatmap)
library(SingleCellExperiment)
library(scater)
library(mclust)
library(scran)
library(igraph)
library(SC3)
library(Seurat)
library(CHETAH)

# load data
setwd("~/Desktop/Nanostring")
load("lupus nephritis SP11_1139 data.RData")
cols = readRDS("kidney cols.RDS")
refprofiles = readRDS("HCA kidney profiles plus nstg immune profiles.RDS")
in.glom = readRDS("in.glom.RDS")

# create HCA reference
HCA_reference <- readRDS("local.rds")
gene_key <- data.frame(ENSG = rownames(HCA_reference@assays[["RNA"]]@meta.features),
                       gene = as.character(HCA_reference@assays[["RNA"]]@meta.features$feature_name))
# correct gene features for CCL3L3 and RGS5
gene_key$gene[gene_key$gene == "CCL3L1"] <- "CCL3L3"
gene_key$gene[substr(gene_key$gene, 1, 4) == "RGS5"][1] <- "RGS5"
rownames(HCA_reference@assays[["RNA"]]@counts) <- gene_key$gene
raw_ref <- HCA_reference@assays[["RNA"]]@counts[rownames(HCA_reference@assays[["RNA"]]@counts) %in% colnames(counts), ]
counts_ref <- t(as.matrix(raw_ref))
rm(raw_ref)
annot_ref <- as.character(HCA_reference@meta.data[["cell_type"]])
meanprofiles <- InSituType:::Estep(counts = counts_ref,
                                   clust = annot_ref,
                                   neg = rep(0, 40268))
## Half genes
set.seed(101)
subsample_half = sample.int(ncol(counts), ncol(counts)/2)  #480 genes
raw_ref_half <- HCA_reference@assays[["RNA"]]@counts[rownames(HCA_reference@assays[["RNA"]]@counts) %in% colnames(counts[,subsample_half]), ]
counts_ref_half <- t(as.matrix(raw_ref_half))
rm(raw_ref_half)

# 96 HVGs (Totoal genes: 960)
## Preprocess
sce_hvg = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                               colData = annot)
sce_hvg = logNormCounts(sce_hvg) #normalize
top = getTopHVGs(sce_hvg, n = 96)
sce_hvg = sce_hvg[top, ]
sce_hvg = scater::runPCA(sce_hvg, subset_row = top)
#Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#You're computing too large a percentage of total singular values, use a standard svd instead.

counts_hvg = counts[,top]
hvg_filter = rowSums(counts_hvg) >= 5
annot_hvg = annot[hvg_filter,]
counts_hvg = counts_hvg[hvg_filter,]
#sce_hvg = sce_hvg[,hvg_filter]
pcs_hvg = reducedDim(sce_hvg)

raw_ref_hvg <- HCA_reference@assays[["RNA"]]@counts[rownames(HCA_reference@assays[["RNA"]]@counts) %in% colnames(counts[,top]), ]
counts_ref_hvg <- t(as.matrix(raw_ref_hvg))
rm(raw_ref_hvg)

# without cohorts
########################### Single R ###########################################
sce = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                           colData = annot)
sce = logNormCounts(sce) #normalize
sce = scater::runPCA(sce)

library(SingleR)
singleR = SingleR(test = sce,
                  ref = meanprofiles,
                  labels = colnames(meanprofiles))

system.time(SingleR(test = sce,
                    ref = meanprofiles,
                    labels = colnames(meanprofiles)))

# user  system elapsed 
#80.002   5.437  88.908 

in.glom <- tibble::rownames_to_column(as.data.frame(in.glom), "cell_ID")
singleR.nocohorts <- data.frame(cell_ID = rownames(singleR),
                                singleR.labels = singleR$pruned.labels)
annot <- merge(merge(annot, in.glom, by = "cell_ID"), singleR.nocohorts, by = "cell_ID")

# half genes
sce_half = SingleCellExperiment(assays = list(counts = t(counts[,subsample_half])), #make SCE
                                colData = annot)
sce_half = logNormCounts(sce_half) #normalize
sce_half = scater::runPCA(sce_half)
pcs_half = reducedDim(sce_half)
singleR_half = SingleR(test = sce_half,
                       ref = meanprofiles,
                       labels = colnames(meanprofiles))

system.time(SingleR(test = sce_half,
                    ref = meanprofiles,
                    labels = colnames(meanprofiles)))
#user  system elapsed 
#49.765   2.467  53.266 

singleR.nocohorts_half <- data.frame(cell_ID = rownames(singleR_half),
                                     SingleR_half.labels = singleR_half$pruned.labels)
annot <- merge(annot, singleR.nocohorts_half, by = "cell_ID")

# 96 HVGs
singleR_hvg = SingleR(test = sce_hvg,
                      ref = meanprofiles,
                      labels = colnames(meanprofiles))
system.time(SingleR(test = sce_hvg,
                    ref = meanprofiles,
                    labels = colnames(meanprofiles)))
#user  system elapsed 
#4.282   0.599   5.227 

singleR.nocohorts_hvg <- data.frame(cell_ID = rownames(singleR_hvg),
                                    SingleR_hvg.labels = singleR_hvg$pruned.labels)
annot <- merge(annot, singleR.nocohorts_hvg, by = "cell_ID")

########################### InsituType #########################################
library(tictoc)
tic()
sup.nocohorts = insitutypeML(x = counts,
                             neg = annot$negmean,
                             reference_profiles = meanprofiles,
                             cohort = NULL)$clust
elapsed_time = toc(quiet = TRUE)
# "99.579 sec elapsed"
sup.nocohorts <- tibble::rownames_to_column(as.data.frame(sup.nocohorts), "cell_ID")
sup.nocohorts <- sup.nocohorts %>% mutate( Insitutype.labels = sup.nocohorts) %>% select(-sup.nocohorts)
annot <- merge(merge(annot, in.glom, by = "cell_ID"), sup.nocohorts, by = "cell_ID")

# half genes
tic()
sup.nocohorts_half = insitutypeML(x = counts[,subsample_half],
                                  neg = annot$negmean,
                                  reference_profiles = meanprofiles,
                                  cohort = NULL)$clust
elapsed_time = toc(quiet = TRUE)
# "36.796 sec elapsed"
sup.nocohorts_half <- tibble::rownames_to_column(as.data.frame(sup.nocohorts_half), "cell_ID")
sup.nocohorts_half <- sup.nocohorts_half %>% mutate(Insitutype_half.labels = sup.nocohorts_half) %>% select(-sup.nocohorts_half)
annot <- merge(annot, sup.nocohorts_half, by = "cell_ID")

# 96 HVGs
tic()
sup.nocohorts_hvg = insitutypeML(x = counts[,top],
                                 neg = annot$negmean,
                                 reference_profiles = meanprofiles,
                                 cohort = NULL)$clust
elapsed_time = toc(quiet = TRUE)
# "19.027 sec elapsed"
sup.nocohorts_hvg <- tibble::rownames_to_column(as.data.frame(sup.nocohorts_hvg), "cell_ID")
sup.nocohorts_hvg <- sup.nocohorts_hvg %>% mutate(Insitutype_hvg.labels = sup.nocohorts_hvg) %>% select(-sup.nocohorts_hvg)
annot <- merge(annot, sup.nocohorts_hvg, by = "cell_ID")

########################### Seurat #############################################
library(Seurat)
rownames(annot) = rownames(counts)
seurat.test = CreateSeuratObject(t(counts), meta.data = annot)
# Warning: Some cells in meta.data not present in provided counts matrix
annot_ref <- as.data.frame(annot_ref)
rownames(annot_ref) = rownames(counts_ref)
seurat.ref = CreateSeuratObject(t(counts_ref), meta.data = annot_ref)
seurat.test = NormalizeData(seurat.test, verbose = FALSE)
seurat.ref = NormalizeData(seurat.ref, verbose = FALSE)
gene_intersect <- intersect(colnames(counts_ref), colnames(counts))
tic()
anchors <- FindTransferAnchors(reference = seurat.ref, query = seurat.test,
                               features = gene_intersect)
elapsed_time = toc(quiet = TRUE)
#"281.799 sec elapsed"
tic()
seurat.predictions <- TransferData(anchorset = anchors, refdata = seurat.ref@meta.data[["annot_ref"]])
elapsed_time = toc(quiet = TRUE)
# "24.397 sec elapsed"

seurat.nocohorts <- data.frame(cell_ID = rownames(seurat.predictions),
                               Seurat.labels = seurat.predictions$predicted.id)
annot <- merge(annot, seurat.nocohorts, by = "cell_ID")

# half genes
seurat.test_half = CreateSeuratObject(t(counts[,subsample_half]), meta.data = annot)
# Warning: Some cells in meta.data not present in provided counts matrix
seurat.ref_half = CreateSeuratObject(t(counts_ref_half), meta.data = annot_ref)
seurat.test_half = NormalizeData(seurat.test_half, verbose = FALSE)
seurat.ref_half = NormalizeData(seurat.ref_half, verbose = FALSE)
gene_intersect_half <- intersect(colnames(counts_ref_half), colnames(counts[,subsample_half]))
tic()
anchors_half <- FindTransferAnchors(reference = seurat.ref_half, query = seurat.test_half,
                                    features = gene_intersect_half)
elapsed_time = toc(quiet = TRUE)
# "360.041 sec elapsed"
tic()
seurat.predictions_half <- TransferData(anchorset = anchors_half, refdata = seurat.ref@meta.data[["annot_ref"]])
elapsed_time = toc(quiet = TRUE)
# "22.933 sec elapsed"
seurat.nocohorts_half <- data.frame(cell_ID = rownames(seurat.predictions_half),
                                    Seurat_half.labels = seurat.predictions_half$predicted.id)
annot <- merge(annot, seurat.nocohorts_half, by = "cell_ID")

# 96 HVGs
seurat.test_hvg = CreateSeuratObject(t(counts_hvg), meta.data = annot)
seurat.ref_hvg = CreateSeuratObject(t(counts_ref_hvg), meta.data = annot_ref)
seurat.test_hvg = NormalizeData(seurat.test_hvg, verbose = FALSE)
seurat.ref_hvg = NormalizeData(seurat.ref_hvg, verbose = FALSE)
gene_intersect_hvg <- intersect(colnames(counts_ref_hvg), colnames(counts_hvg))
tic()
anchors_hvg <- FindTransferAnchors(reference = seurat.ref_hvg, query = seurat.test_hvg,
                                   features = gene_intersect_hvg)
elapsed_time = toc(quiet = TRUE)
# "578.547 sec elapsed"
seurat.predictions_hvg <- TransferData(anchorset = anchors_hvg, refdata = seurat.ref@meta.data[["annot_ref"]])
seurat.nocohorts_hvg <- data.frame(cell_ID = rownames(seurat.predictions_hvg),
                                   Seurat_hvg.labels = seurat.predictions_hvg$predicted.id)
annot <- merge(annot, seurat.nocohorts_hvg, by = "cell_ID")


########################### CHETAH #############################################
## For the reference we define a "counts" assay and "celltypes" metadata
sce_ref = SingleCellExperiment(assays = list(counts = t(counts_ref)),
                               colData = annot_ref)
sce_ref = sce_ref[,colSums(counts(sce_ref))>0]
sce_ref = logNormCounts(sce_ref)
## For the input we define a "counts" assay and "TSNE" reduced dimensions
tic()
sce = runTSNE(sce)
elapsed_time = toc(quiet = TRUE)
# "492.171 sec elapsed"
tic()
sce <- CHETAH::CHETAHclassifier(input = sce,
                                ref_cells = sce_ref,
                                ref_ct = "annot_ref")
elapsed_time = toc(quiet = TRUE)
# "2817.674 sec elapsed"
CHETAH.nocohorts <- data.frame(cell_ID = names(sce$celltype_CHETAH),
                               CHETAH.labels = sce$celltype_CHETAH)
annot <- merge(annot, CHETAH.nocohorts, by = "cell_ID")

# half genes
sce_ref_half = SingleCellExperiment(assays = list(counts = t(counts_ref_half)),
                                    colData = annot_ref)
sce_ref_half = sce_ref_half[,colSums(counts(sce_ref_half))>0]
sce_ref_half = logNormCounts(sce_ref_half)
tic()
sce_half = runTSNE(sce_half)
elapsed_time = toc(quiet = TRUE)
# "604.379 sec elapsed"
tic()
sce_half <- CHETAH::CHETAHclassifier(input = sce_half,
                                     ref_cells = sce_ref_half,
                                     ref_ct = "annot_ref")
elapsed_time = toc(quiet = TRUE)
# "3780.206 sec elapsed"
CHETAH.nocohorts_half <- data.frame(cell_ID = names(sce_half$celltype_CHETAH),
                                    CHETAH_half.labels = sce_half$celltype_CHETAH)
annot <- merge(annot, CHETAH.nocohorts_half, by = "cell_ID")

# 96 HVGs 
## Note: CHETAH is not compatible with the top 10% of the most highly variable genes. 
sce_ref_hvg = SingleCellExperiment(assays = list(counts = t(counts_ref_hvg)),
                                   colData = annot_ref)
sce_ref_hvg = sce_ref_hvg[,colSums(counts(sce_ref_hvg))>0]
sce_ref_hvg = logNormCounts(sce_ref_hvg)
tic()
sce_hvg = runTSNE(sce_hvg)
elapsed_time = toc(quiet = TRUE)
# "419.214 sec elapsed"
tic()
sce_hvg <- CHETAH::CHETAHclassifier(input = sce_hvg,
                                    ref_cells = sce_ref,
                                    ref_ct = "annot_ref")
# Error in ref_profiles[genes, type, drop = FALSE] : subscript out of bounds
