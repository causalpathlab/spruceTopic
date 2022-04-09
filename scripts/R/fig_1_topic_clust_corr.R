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


loss_plot <- function(args) {
loss_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"top_5_topic_cluster_corrCD4.tsv",sep="")
df = read.table(loss_file, sep = "\t", header=TRUE)

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p <- pheatmap(df,color=mat_colors,cluster_rows=FALSE,cluster_cols=FALSE,fontsize_row=10,fontsize_col=8)

ggsave("top_5_topic_cluster_corrCD4.pdf",p)
}

loss_plot(args)