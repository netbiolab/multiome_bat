library(Seurat)
library(ggplot2)
library(patchwork)
library(ktools)
library(Pando)
library(rtracklayer)
library(tidyr)
library(igraph)
library(enrichR)
library(stringr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

## Functions
get_cent <- function(Seurat) {
    net_file <- NetworkGraph(Seurat)
    df <- igraph::as_data_frame(net_file)
    df_net <- df[,c("from_node","to_node","corr","estimate")]
    colnames(df_net) <- c("From","To","weight","estimate")
    net <- graph.data.frame(df_net)
    E(net)$regulation <- "positive"
    E(net)$regulation[which(E(net)$weight < 0)] <- "negative"
    E(net)$weight <- abs(E(net)$weight)
    bet <- betweenness(net)
    bet <- bet[order(bet,decreasing=TRUE)]
    deg <- degree(net)
    deg <- deg[order(deg,decreasing=TRUE)]
    str <- strength(net)
    str <- str[order(str,decreasing=TRUE)]
    bet <- as.data.frame(bet)
    deg <- as.data.frame(deg)
    str <- as.data.frame(str)
    bet <- mutate(bet, bet_pr = rank(bet,ties.method="min")/length(bet))
    deg <- mutate(deg, deg_pr = rank(deg,ties.method="min")/length(deg))
    str <- mutate(str, str_pr = rank(str,ties.method="min")/length(str))
    bet <- bet[rownames(bet),]
    deg <- deg[rownames(bet),]
    str <- str[rownames(bet),]
    cent_df <- do.call('cbind', list(bet,deg,str))
    result <- list(net,cent_df)
    return(result)
}

combine_cent_result <- function(fruit_cent,insect_cent,celltype) {
    # total_genes <- unique(c(rownames(fruit_cent), rownames(insect_cent)))    
    colnames(fruit_cent) <- paste0("fruit_",colnames(fruit_cent))
    colnames(insect_cent) <- paste0("insect_",colnames(insect_cent))
    fruit_cent$gene <- rownames(fruit_cent)
    insect_cent$gene <- rownames(insect_cent)
    comb_cent <- merge(fruit_cent,insect_cent,by="gene",all=TRUE)
    comb_cent$exist <- "both"
    comb_cent$exist[is.na(comb_cent$fruit_bet)] <- "insect"
    comb_cent$exist[is.na(comb_cent$insect_bet)] <- "fruit"
    comb_cent[is.na(comb_cent)] <- 0
    comb_cent$bet_diff <- comb_cent$fruit_bet_pr - comb_cent$insect_bet_pr
    comb_cent$deg_diff <- comb_cent$fruit_deg_pr - comb_cent$insect_deg_pr
    comb_cent$str_diff <- comb_cent$fruit_str_pr - comb_cent$insect_str_pr
    comb_cent <- comb_cent[order(abs(comb_cent$bet_diff), decreasing=TRUE),]
    comb_cent$gene_convert <- comb_cent$gene
    comb_cent$gene_convert[which(comb_cent$gene_convert %in% as.character(og_convert$OG))] <- og_convert$gene[match(comb_cent$gene_convert[which(comb_cent$gene_convert %in% as.character(og_convert$OG))], og_convert$OG)]
    comb_cent <- comb_cent[,c(1,18,14,2,8,3,9,15,4,10,5,11,16,6,12,7,13,17)]
    a <- ggplot(comb_cent, aes(x=fruit_bet_pr, y=insect_bet_pr, color=bet_diff)) +
    geom_point() +
    geom_text_repel(aes(label=gene_convert),size=2) +
    scale_color_gradient2(low="#6747A2",mid='#AA867D',high="#EDC458") +
    theme_bw() +
    labs(title=paste0(celltype," Betweenness Difference"))
    b <- ggplot(comb_cent, aes(x=fruit_deg_pr, y=insect_deg_pr, color=deg_diff)) +
    geom_point() +
    geom_text_repel(aes(label=gene_convert),size=2) +
    scale_color_gradient2(low="#6747A2",mid='#AA867D',high="#EDC458") +
    theme_bw() +
    labs(title=paste0(celltype," Degree Difference"))
    c <- ggplot(comb_cent, aes(x=fruit_str_pr, y=insect_str_pr, color=str_diff)) +
    geom_point() +
    geom_text_repel(aes(label=gene_convert),size=2) +
    scale_color_gradient2(low="#6747A2",mid='#AA867D',high="#EDC458") +
    theme_bw() +
    labs(title=paste0(celltype," Strength Difference"))
    plot <- a|b|c
    result <- list(comb_cent,plot)
}

run_pathway <- function(fruit_gene,insect_gene,celltype) {
    fruit_res <- enrichr(fruit_gene,sets)
    insect_res <- enrichr(insect_gene,sets)
    comb_df_all <- as.data.frame(matrix(nrow=0,ncol=0))
    for (i in c(1:length(sets))) {
        terms <- unique(c(
            head(fruit_res[[i]]$Term,term_count),
            head(insect_res[[i]]$Term,term_count)
        ))
        fruit_df <- data.frame(terms = fruit_res[[i]]$Term[which(fruit_res[[i]]$Term %in% terms)],fruit = fruit_res[[i]]$Adjusted.P.value[which(fruit_res[[i]]$Term %in% terms)])
        insect_df <- data.frame(terms = insect_res[[i]]$Term[which(insect_res[[i]]$Term %in% terms)],insect = insect_res[[i]]$Adjusted.P.value[which(insect_res[[i]]$Term %in% terms)])
        comb_df <- merge(fruit_df, insect_df, by='terms',all=TRUE)
        rownames(comb_df) <- comb_df$terms
        comb_df$terms <- NULL
        comb_df[is.na(comb_df)] <- 1
        comb_df <- -log(comb_df,base=10)
        comb_df$pathway <- sets[i]
        comb_df_all <- rbind(comb_df_all, comb_df)
    }
    row_color <- data.frame(pathway=comb_df_all$pathway)
    rownames(row_color) <- rownames(comb_df_all)
    comb_df_all$pathway <- NULL
    comb_df_all <- comb_df_all[comb_df_all$fruit > -log(0.05,base=10) | comb_df_all$insect > -log(0.05,base=10),]
    inc <- (max(comb_df_all)+log(0.05,base=10)) / 5
    a <- pheatmap(comb_df_all[apply(comb_df_all, 1, max, na.rm=TRUE) > -log(0.05,base=10),], cluster_cols=TRUE, cluster_rows=TRUE, treeheight_row=0, treeheight_col=0, cellwidth=10, cellheight=10, color=c('white',brewer.pal(n = 6, name ="Greens")[2:6]),main=paste0(celltype,"\nFruit:",length(fruit_gene)," insect:",length(insect_gene)),breaks= c(0,-log(0.05,base=10),seq(from=inc+(-log(0.05,base=10)), to=max(comb_df_all), by=inc)), silent=TRUE,annotation_row=row_color)
}


## define title
title <- "celltype"

## load pando grn results for each species
Seurat_f <- readRDS("~/Seurat_grn_fruit.rds")
Seurat_i <- readRDS("~/Seurat_grn_insect.rds")

## load OG to gene name file
og_convert <- read.csv("~/OG_convert.csv",header=TRUE)


## get centrality
print("Getting Centrality and compare")
result_fruit <- get_cent(Seurat_f)
result_insect <- get_cent(Seurat_i)

fruit_net <- result_fruit[[1]]
insect_net <- result_insect[[1]]

name <- V(fruit_net)$name
name[which(name %in% as.character(og_convert$OG))] <- og_convert$gene[match(name[which(name %in% as.character(og_convert$OG))], og_convert$OG)]
V(fruit_net)$name <- name
name <- V(insect_net)$name
name[which(name %in% as.character(og_convert$OG))] <- og_convert$gene[match(name[which(name %in% as.character(og_convert$OG))], og_convert$OG)]
V(insect_net)$name <- name

E(fruit_net)$color <- "#B71103"
E(fruit_net)$color[which(E(fruit_net)$regulation == "positive")]  <- "#378805"
E(insect_net)$color <- "#B71103"
E(insect_net)$color[which(E(insect_net)$regulation == "positive")]  <- "#378805"

pdf(file=paste0(title,"_Network.pdf"), width=28 ,height=14)
par(mfrow=c(1,2))
plot(fruit_net,vertex.size=5,vertex.label.cex=0.8,vertex.label.color = "black", layout=layout.fruchterman.reingold, vertex.color='grey')
title(paste0(title," Fruit Network\nNode:",length(V(fruit_net)$name)))
plot(insect_net,vertex.size=5,vertex.label.cex=0.8,vertex.label.color = "black", layout=layout.fruchterman.reingold, vertex.color='grey')
title(paste0(title," Insect Network\nNode:",length(V(insect_net)$name)))
dev.off()

fruit_cent <- result_fruit[[2]]
insect_cent <- result_insect[[2]]

comb_result <- combine_cent_result(fruit_cent, insect_cent,title)

pdf(file=paste0(title,"_centraility_comparison.pdf"), width=21, height=7)
print(comb_result[[2]])
dev.off()

## centrality pathway
print("Comparing pathways genes with high strength centrality")

cent <- comb_result[[1]]
fruit_gene <- cent$gene_convert[cent$fruit_str_pr >= max(0.5,sort(cent$fruit_str_pr,decreasing=TRUE)[100]) & cent$insect_str_pr < 0.75]
insect_gene <- cent$gene_convert[cent$fruit_str_pr < 0.75 & cent$insect_str_pr >= max(0.5,sort(cent$insect_str_pr,decreasing=TRUE)[100])]

sets <- c(
    "Reactome_2022",
    "BioPlanet_2019",
    "KEGG_2021_Human",
    "WikiPathway_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "MSigDB_Hallmark_2020",
    "HuBMAP_ASCTplusB_augmented_2022",
    "PanglaoDB_Augmented_2021"
)


if (length(fruit_gene) > 10 & length(insect_gene) > 10) {
    print("Running pathway analysis")
    term_count <- 5
    plot <- run_pathway(fruit_gene,insect_gene,title)
    pdf(file=paste0(title,"_pheatmap_pathway.pdf"),width=13, height=(1.5+length(plot$tree_row$height)*0.175))
    print(plot)
    dev.off()
} else {
    print("Not enough genes")
}