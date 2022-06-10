library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
setwd(box::file())
source("Util.R")

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/spruce_topic/"
config = paste(args_home,"/config/",args[1],".yaml",sep="") 
args = read_yaml(config)

weightmat_lr_plot <- function(args) {

topgenes_file = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"_ietm_top_5_genes_topic.tsv.gz",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)
print(head(df_tg))
df_tg = df_tg[df_tg$GeneType=="ligands",]
top_genes = unique(df_tg$Gene)
print(top_genes)


beta1 = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"_ietm_beta2.tsv.gz",sep="")
beta1_cols = read.table(paste(args_home,args$input,'ligands.csv.gz',sep=''),header=TRUE)

df_beta1 = read.table(beta1, sep = "\t", header=TRUE)

colnames(df_beta1) = beta1_cols$X0

df_beta1 = select(df_beta1,all_of(top_genes))
print(colnames(df_beta1))

row_order = row.order(df_beta1)

df_beta_t = df_beta1
df_beta_t$topic = rownames(df_beta1)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
col_order = col.order(df_beta_t,row_order)

df_beta1 = df_beta1[,col_order]
df_beta1 = df_beta1[row_order,]

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)



p1 <- pheatmap(df_beta1,color = mat_colors,legend =FALSE,fontsize_row=2,fontsize_col=4,cluster_rows=FALSE,cluster_cols=FALSE)

f = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"hmap_clust_ligands_tg.pdf",sep="")
ggsave(f,p1)

}
weightmat_lr_plot(args)



