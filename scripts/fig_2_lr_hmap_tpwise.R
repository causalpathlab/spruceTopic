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
args_home ="/home/BCCRC.CA/ssubedi/projects/spruce_topic/"
config = paste(args_home,"/config/",args[1],".yaml",sep="") 
args = read_yaml(config)


weightmat_phmap_plot_i <- function(args) {
topgenes_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"top_5_lrpair_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)


rownames(df_tg) = df_tg$X
df_tg$X = NULL

df_tg = sqrt(df_tg)
colnames(df_tg) <- gsub("X","",colnames(df_tg))

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df_tg,color = mat_colors,legend =FALSE,fontsize_row=2,fontsize_col=4)

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"top5_lr_pair_topic_hmap.pdf",sep="")
ggsave(f,p1)

}

weightmat_phmap_plot_i(args)
