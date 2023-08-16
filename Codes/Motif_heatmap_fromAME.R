library(Seurat)
library(ggplot2)
library(Signac)
library(ggpubr)
library(pheatmap)
library(tidyr)
library(ktools)
library(RColorBrewer)
library(ggmisc)
library(dplyr)
library(gdata)

convert_to_heatmap_rank <- function(mat,name) {
    mat <- mutate(mat, percent_rank = rank(-rank) / length(rank))
    df <- data.frame(motif=mat$motif_alt_ID,P_rank=mat$percent_rank)
    df <- df[!duplicated(df$motif),]
    colnames(df)[2] <- name
    return(df)
}

convert_to_heatmap_pval <- function(mat,name) {
    pvalue <- mat$adj_p.value
    if (min(pvalue) == 0) {
        pvalue[pvalue == 0] <- sort(pvalue)[sort(pvalue) != 0][1]
    }
    mat$adj_p.value_update <- pvalue
    # mat <- mutate(mat, percent_rank = rank(-rank) / length(rank))
    df <- data.frame(motif=mat$motif_alt_ID,p_val=-log(mat$adj_p.value_update,base=10))
    df <- df[!duplicated(df$motif),]
    colnames(df)[2] <- name
    return(df)
}

draw_motif_heatmap <- function(mat_ls,title) {
    mat_fin_rank <- mat_fin_pval <- as.data.frame(matrix(nrow=0,ncol=1))
    colnames(mat_fin_rank) <- colnames(mat_fin_pval) <- c("motif")
    for (i in seq_along(mat_ls)) {
        mat_tmp_rank <- convert_to_heatmap_rank(mat_ls[[i]], names(mat_ls)[i])
        mat_tmp_pval <- convert_to_heatmap_pval(mat_ls[[i]], names(mat_ls)[i])
        mat_fin_rank <- merge(mat_fin_rank, mat_tmp_rank,by="motif", all=TRUE)
        mat_fin_pval <- merge(mat_fin_pval, mat_tmp_pval,by="motif", all=TRUE)
    }
    rownames(mat_fin_rank) <- mat_fin_rank$motif
    mat_fin_rank$motif <- NULL
    rownames(mat_fin_pval) <- mat_fin_pval$motif
    mat_fin_pval$motif <- NULL
    ## TOP 10% and more than 1
    mat_fin_rank_fil <- mat_fin_rank[apply(mat_fin_rank,1,function(x) {
        return(sum(x > 0.9, na.rm=TRUE))
    }) > 1,]
    mat_fin_rank_fil[is.na(mat_fin_rank_fil)] <- 0
    new_motif_top10 <- rownames(mat_fin_rank_fil)

    mat_fin_pval_fil <- mat_fin_pval[rownames(mat_fin_rank_fil),]
    mat_fin_pval_fil[is.na(mat_fin_pval_fil)] <- 0
    colors <- list(species = c(Insectivore = "#6747A2",Frugivore = "#EDC458"))
    col_color <- data.frame(species = unlist(lapply(strsplit(colnames(mat_fin_rank_fil),"_"), function(x) {return(tail(x,1))})))
    rownames(col_color) <- colnames(mat_fin_rank_fil)
    col_color$species <- factor(col_color$species, levels=c("Insectivore","Frugivore"))
    # labels_col <- unlist(lapply(strsplit(colnames(mat_fin_rank_fil),"_"), function(x) {return(paste(head(x,-1), collapse=" "))}))

    a <- pheatmap(mat_fin_rank_fil, color = colorRampPalette(c("white",brewer.pal(n = 6, name ="Greens")))(100), cluster_rows = TRUE, cluster_cols = TRUE, annotation_colors = colors, annotation_col = col_color, angle_col = 45,silent=TRUE, cellwidth=15, cellheight=12.5,border_color='#c5c5c5', main=paste0(title," ranked 10% & non-unique(percentile rank)"))
    ggsave(a, file=paste0(title,"_motif_rank_heatmap_ranked10_fin.pdf"), width=8, height=(3.5 + nrow(mat_fin_rank_fil) * 0.25), dpi=600, useDingbats=FALSE)
    return(a)
}

## load data

insect_ame_CT1 <- read.table("~/insect_ame_celltype1.tsv", sep='\t', header=TRUE)
insect_ame_CT2 <- read.table("~/insect_ame_celltype2.tsv", sep='\t', header=TRUE)
fruit_ame_CT1 <- read.table("~/fruit_ame_celltype1.tsv", sep='\t', header=TRUE)
fruit_ame_CT2 <- read.table("~/fruit_ame_celltype2.tsv", sep='\t', header=TRUE)

mat_ls <- list(insect_ame_CT1,fruit_ame_CT1,insect_ame_CT2,fruit_ame_CT2)
names(mat_ls) <- c('insect_ame_CT1','fruit_ame_CT1','insect_ame_CT2','fruit_ame_CT2')

a <- draw_motif_heatmap(mat_ls,"celltype_comb1")