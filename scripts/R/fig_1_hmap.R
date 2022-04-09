library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)


config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"


weightmat_phmap_plot_i <- function(args) {
topgenes_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"top_5_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)


tcells_genes_cd48 = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"top_5_genes_topic_CDMarkers.tsv",sep="")
df_tcell_genes_cd48 = read.table(tcells_genes_cd48, header=TRUE,sep="\t")


# df_tg <- df_tg[order(df_tg$Genes),]
df_tg$Topic <- factor(df_tg$Topic,levels = c("k0", "k1", "k2", "k3", "k4","k5",
      "k6","k7","k8","k9","k10","k11"))
df_tg$Proportion = sqrt(df_tg$Proportion)
df_topic_i = as.matrix(cast( df_tg[df_tg$GeneType=="immune",] , Topic~Gene) )
df_topic_ni = as.matrix(cast( df_tg[df_tg$GeneType=="non-immune",] , Topic~Gene) )

i_clust = data.frame(df_tcell_genes_cd48[df_tcell_genes_cd48$Gene %in% colnames(df_topic_i),])
i_clust$CD4 = factor(i_clust$CD4)
i_clust$CD8 = factor(i_clust$CD8)
rownames(i_clust)= i_clust$Gene
i_clust[order(i_clust$Gene),]
i_clust = i_clust[,c('CD4','CD8')]

newCols <- colorRampPalette(brewer.pal(length(unique(i_clust$CD8)),"Paired"))
cd4c <- newCols(length(unique(i_clust$CD4)))
names(cd4c) <- unique(i_clust$CD4)
cd8c <- newCols(length(unique(i_clust$CD8)))
names(cd8c) <- unique(i_clust$CD8)

my_colour = list(
    CD4 = cd4c,
    CD8 = cd8c
    )
mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df_topic_i,legend=FALSE,color = mat_colors,annotation_col = i_clust ,annotation_colors = my_colour,
cluster_rows=FALSE,cluster_cols=FALSE,fontsize_row=6,fontsize_col=6)

ggsave("topic_imm.pdf",p1)

# plotlist = list()
# plotlist[[1]] = p1 
# plotlist[[2]] = p2 

# stplt <- grid.arrange(grobs=plotlist,ncol=1,
# heights = c(1/2, 1/2))
# ggsave("pp_weightmat_phmap_plot.pdf",stplt)
}

weightmat_phmap_plot_i(args)
