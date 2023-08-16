library(Seurat)
library(ggplot2)
library(patchwork)
library(Pando)
library(rtracklayer)
library(doParallel)
registerDoParallel(10) ## adjust core count
library(tidyr)
# custom BSgenome
library(BSgenome.Ajamaicensis.NCBI.WHUAjamv2)
library(BSgenome.Efuscus.NCBI.EptFus1.0)
data(motifs)
data('motif2tf')

## genome
fruit_gtf <- rtracklayer::import("~/Ajam_genes_modified_final.gtf")
fruit_genome <- fruit_gtf[fruit_gtf$type == 'exon']
fruit_genome$gene_name <- fruit_genome$gene_id

insect_gtf <- rtracklayer::import("~/Efus_genes_modified_final.gtf")
insect_genome <- insect_gtf[insect_gtf$type == 'exon']
insect_genome$gene_name <- insect_genome$gene_id

## variable gene
var <- readRDS("~/variable_genes.rds")

## Seurat obj
Seurat <- readRDS("~/Seurat_object.rds")

Seurat_grn <- initiate_grn(
    Seurat,
    rna_assay = 'RNA',
    peak_assay = 'peaks',
    regions = fruit_genome, # or insect_genome
    exclude_exon=FALSE
)

Seurat_grn <- find_motifs(
    Seurat_grn,
    pfm = motifs,
    motif_tfs = motif2tf,
    genome = BSgenome.Ajamaicensis.NCBI.WHUAjamv2 # or BSgenome.Efuscus.NCBI.EptFus1.0
)

Seurat_grn <- infer_grn(
    Seurat_grn,
    peak_to_gene_method = 'Signac',
    # genes = VariableFeatures(Seurat_grn, assay='RNA'),
    genes = var, 
    parallel = TRUE
)

Seurat_grn <- find_modules(
    Seurat_grn, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

Seurat_grn <- get_network_graph(Seurat_grn)

saveRDS(Seurat_grn,file="./Seurat_grn.rds")
