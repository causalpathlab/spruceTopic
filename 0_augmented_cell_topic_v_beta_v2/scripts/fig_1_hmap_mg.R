library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
setwd(box::file())
source("Util.R")

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)


weightmat_phmap_plot_i <- function(args) {
topgenes_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"marker_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)
df_tg = sqrt(df_tg)

rownames(df_tg) = seq(0,14,by=1)

# df_tg$Proportion = log(df_tg$Proportion)
# df_topic_i = as.matrix(cast( df_tg[df_tg$GeneType=="immune",] , Topic~Gene) )
# df_topic_ni = as.matrix(cast( df_tg[df_tg$GeneType=="non-immune",] , Topic~Gene) )

# xlab = paste("Marker", "genes")
mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df_tg,color = mat_colors,legend =FALSE,fontsize_row=12,fontsize_col=12)

f = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"marker_gene_weight_hmap.pdf",sep="")

ggsave(f,p1)

}

weightmat_phmap_plot_i(args)
