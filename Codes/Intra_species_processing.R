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

### genome loading ###
## Ajam (fruit bat) genome

ajam.gtf <- rtracklayer::import('~/Ajam_genes_modified_final.gtf')
ajam.gene.coords <- ajam.gtf[ajam.gtf$type == 'exon']
ajam.gene.coords$gene_name <- ajam.gene.coords$gene_id

## Efus (insect bat) genome
efus.gtf <- rtracklayer::import('~/Efus_genes_modified_final.gtf')
efus.gene.coords <- efus.gtf[efus.gtf$type == 'exon']
efus.gene.coords$gene_name <- efus.gene.coords$gene_id

### Create Seurat/Signac Object
## Create ATAC object
makeObject <- function(h5.file, fragpath, sample) {
  inputdata.10x <- Read10X_h5(h5.file)
  atac_counts <- inputdata.10x$Peaks
  # Create Seurat object
  chrom_assay <- CreateChromatinAssay(counts = atac_counts, 
                                      sep = c(":", "-"), 
                                      fragments = fragpath, 
                                      annotation = ajam.gene.coords)
  x <- CreateSeuratObject(chrom_assay, assay = "peaks", project = sample)
  # delete to make space
  return(x)
}

ATAC_obj1 <- makeObject("~/sample1/raw_feature_bc_matrix.h5","~/sample1/atac_fragments.tsv.gz","sample_name1") # repeat for all samples
ATAC_obj2 <- makeObject("~/sample2/raw_feature_bc_matrix.h5","~/sample2/atac_fragments.tsv.gz","sample_name2") # repeat for all samples
ATAC_obj3 <- makeObject("~/sample3/raw_feature_bc_matrix.h5","~/sample3/atac_fragments.tsv.gz","sample_name3") # repeat for all samples
ATAC_obj4 <- makeObject("~/sample4/raw_feature_bc_matrix.h5","~/sample4/atac_fragments.tsv.gz","sample_name4") # repeat for all samples

obj_ls <- list(ATAC_obj1, ATAC_obj2, ATAC_obj3, ATAC_obj4) # list of ATAC objects

## Convert ATAC objects to genomic ranges
gr_obj_ls <- lapply(obj_ls, function(x) {
    return(granges(x))
})

## Create a unified set of peaks to quantify in each dataset (function will only return scaffolds that are common between the replicates, so warning message lets you know which scaffolds were not shared and thus thrown out)
combined_peaks <- reduce(x=unlist(gr_obj_ls))

## quantify peaks in each dataset
fragment_count_ls <- lapply(obj_ls, function(x) {
    fragment_count_obj <- FeatureMatrix(
        fragments = Fragments(x),
        features = combined_peaks,
        cells = Cells(x))
})

## Create new ATAC objects using the common feature peaks
new_obj_ls <- list()
for (i in seq_along(fragment_count_ls)) {
    assay <- CreateChromatinAssay(
        counts = fragment_count_ls[[i]],
        min.features = 100,
        fragments = Fragments(obj_ls[[i]]),
        annotation = ajam.gene.coords # choose the correct species' genome
    )
    object <- CreateSeuratObject(counts = assay, assay = "peaks", project = paste0("sample_name",i))
    new_obj_ls[[i]] <- object
}

## Create RNA object
h5_dir_ls <- c( # use the same sets of samples as ATAC objects with the same order
    "~/sample1/raw_feature_bc_matrix.h5",
    "~/sample2/raw_feature_bc_matrix.h5",
    "~/sample3/raw_feature_bc_matrix.h5",
    "~/sample4/raw_feature_bc_matrix.h5"
)

for (i in seq_along(h5_dir_ls)) {
    ATAC_obj <- new_obj_ls[[i]]
    input_10x <- Read10X_h5(file=h5_dir_ls[i])
    rna_counts <- input_10x$"Gene Expression"
    rna_counts <- rna_counts[, colnames(ATAC_obj)]
    rna_assay <- CreateAssayObject(counts = rna_counts, 
                          min.cells = 3, 
                          min.features = 10)
    # use only cells that co-exist for both assay
    if (ncol(ATAC_obj) != ncol(rna_assay)) {
        ATAC_obj <- ATAC_obj[,colnames(rna_assay@counts)]
    }
    # add RNA assay to existing object
    ATAC_obj[["RNA"]] <- rna_assay
    new_obj_ls[[i]] <- ATAC_obj
}

### Quality Control
mitochondrial <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB") # Only for fruit bat genome

run_qc <- function(x) {
  DefaultAssay(x) <- "RNA"
  # Choose the right species
  x[["percent.mt"]] <- PercentageFeatureSet(x, features = mitochondrial) # for fruit bat
#   x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-") # for insect bat
  DefaultAssay(x) <- "peaks"
  x <- NucleosomeSignal(x)
  x <- TSSEnrichment(x)
  # adjust QC parameters
  x_sub <- subset(x = x, subset = percent.mt < 25 &
                  nCount_peaks < 100000 &
                  nCount_peaks > 500 &
                  nCount_RNA < 25000 &
                  nCount_RNA > 200 &
                  nucleosome_signal < 2 &
                  TSS.enrichment > 1)
  x <- x_sub
  return(x)
}

new_obj_ls <- lapply(new_object_ls, function(Seurat) {
    Seurat_sub <- run_qc(Seurat)
    return(Seurat_sub)
})

### Add Gene score
add_genescore <- function(x) {
    DefaultAssay(x) <- "peaks"
    x.geneactivities <- GeneActivity(x)
    x[['genescore']] <- CreateAssayObject(counts = x.geneactivities)
    DefaultAssay(x) <- "genescore"
    x[["peaks"]] <- NULL
    return(x)
}

new_obj_ls <- lapply(new_object_ls, function(Seurat) {
    Seurat <- add_genescore(Seurat)
    return(Seurat)
})

### Run SCTransform for each sample
Run_SCTransform <- function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
  return(x)
}

new_obj_ls <- lapply(new_object_ls, function(Seurat) {
    Seurat <- Run_SCTransform(Seurat)
    return(Seurat)
})

### Save the Seurat object list
saveRDS(new_obj_ls, file="processed_Seurat_object_list_species1.rds")
