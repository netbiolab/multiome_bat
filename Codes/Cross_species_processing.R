### library loading ###
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(future)
library(harmony)
library(dplyr)
library(patchwork)
library(reshape2)
set.seed(1234)

### Caution
## Run after processing each species' samples with inter_species_processing.R

### Sample loading 
obj_ls_species1 <- readRDS("~/processed_Seurat_object_list_species1.rds")
obj_ls_species2 <- readRDS("~/processed_Seurat_object_list_species2.rds")

### Merging
Species1 <- merge(x=obj_ls_species1[[1]],y=obj_ls_species1[c(2:length(obj_ls_species1))], add.cell.ids = c("sample_name1","sample_name2","sample_name3","sample_name4"), merge.data=TRUE)

Species2 <- merge(x=obj_ls_species2[[1]],y=obj_ls_species2[c(2:length(obj_ls_species2))], add.cell.ids = c("sample_name5","sample_name6","sample_name7","sample_name8"), merge.data=TRUE)

## Find variable features for genescore by species
DefaultAssay(Species1) <- "genescore"
Species1 <- NormalizeData(Species1, normalization.method = 'LogNormalize', scale.factor = median(fruit$nCount_genescore))
Species1 <- FindVariableFeatures(Species1, selection.method = "vst", nfeatures = 2000)
Species1

DefaultAssay(Species2) <- "genescore"
Species2 <- NormalizeData(Species2, normalization.method = 'LogNormalize', scale.factor = median(Species2$nCount_genescore))
Species2 <- FindVariableFeatures(Species2, selection.method = "vst", nfeatures = 2000)
Species2

## Merge datasets
combined <- merge(x = Species1, y = Species2, merge.data = TRUE)
obj_ls_combined <- append(obj_ls_species1, obj_ls_species2)

## RNA integration via SCTransform
DefaultAssay(combined) <- "SCT"
features <- SelectIntegrationFeatures(object.list = obj_ls_combined, nfeatures = 3000) # Get variable genes for SCT
combined@assays$SCT@var.features <- features
combined <- RunPCA(combined, assay = 'SCT', features = features)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30, assay = 'SCT')

combined <- RunHarmony(combined, assay.use='SCT', max.iter.harmony=20, group.by.vars=c("orig.ident","species"), dims.use=1:30, reduction = "pca", reduction.save = "harmony", plot_convergence = T)

## ATAC integration via genescore
DefaultAssay(combined) <- "genescore"
species_ls <- list(Species1, Species2)
genescore_features <- SelectIntegrationFeatures(object.list = species_ls, nfeatures = 2000) # Get variable genes for SCT
combined@assays$genescore@var.features <- genescore_features
combined <- ScaleData(combined, verbose = F)
combined <- RunPCA(combined, verbose = F)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30, verbose = F)
combined <- RunHarmony(combined, assay.use = 'genescore', max.iter.harmony=20, group.by.vars = c("orig.ident","species"), dims.use = 1:30, reduction = 'pca', reduction.save = "harmony_genescore", plot_convergence = T)

## WNN integration
combined <- FindMultiModalNeighbors(combined, 
                               reduction.list = list("harmony", "harmony_genescore"), 
                               dims.list = list(1:30, 1:30),
                               modality.weight.name = c("SCT_har","genescore_har"),
                               weighted.nn.name = "harmony.nn",
                               verbose = TRUE)
combined <- RunUMAP(combined, nn.name = "harmony.nn", reduction.name = "harmonywnn.umap", reduction.key = "harmonyWNNUMAP_", assay = "RNA", verbose = TRUE)

## clustering
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, resolution = .6, verbose = FALSE)

saveRDS(combined, file="combined_Seurat_object.rds")
