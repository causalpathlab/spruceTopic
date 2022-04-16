library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
source("Util.R")


config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"


weightmat_phmap_plot_i <- function(args) {
topgenes_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"top_1_lrpair_topic_cd4cd8_top_genes.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)


rownames(df_tg) = df_tg$X
df_tg$X = NULL

df_tg = sqrt(df_tg)
colnames(df_tg) <- gsub("X","",colnames(df_tg))

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df_tg,color = mat_colors,legend =FALSE,fontsize_row=4,fontsize_col=4)

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"topic_lr_topic_cd4cd8_top_genes.pdf",sep="")
ggsave(f,p1)

}

weightmat_phmap_plot_i(args)
