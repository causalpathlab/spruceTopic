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
args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/5_augmented_lr_multinomial_dirichlet_beta_spike/"
config = paste(args_home,"/config/",args[1],".yaml",sep="") 
args = read_yaml(config)

weightmat_lr_plot <- function(beta1,beta_cols,f) {


df_beta = read.table(beta, sep = "\t", header=TRUE)

colnames(df_beta) = beta_cols$X0

# topgenes_file = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"_ietm_top_5_genes_topic.tsv.gz",sep="")
# df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)
# print(head(df_tg))
# df_tg = df_tg[df_tg$GeneType=="receptors",]
# top_genes = unique(df_tg$Gene)
# print(top_genes)
# df_beta1 = select(df_beta1,all_of(top_genes))
# print(colnames(df_beta1))

row_order = row.order(df_beta)

df_beta_t = df_beta
df_beta_t$topic = rownames(df_beta)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
col_order = col.order(df_beta_t,row_order)

df_beta = df_beta[,col_order]
df_beta = df_beta[row_order,]

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "RdYlGn"))(100)

p1 <- pheatmap(df_beta,color = mat_colors,fontsize_row=2,fontsize_col=4,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE)

ggsave(f,p1)

}

beta = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"_ietm_beta1.tsv.gz",sep="")
beta_cols = read.table(paste(args$home,args$data,args$sample_id,'receptors.csv.gz',sep=''),header=TRUE)
f = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"hmap_clust_receptors_all.pdf",sep="")

weightmat_lr_plot(beta,beta_cols,f)

beta = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"_ietm_beta2.tsv.gz",sep="")
beta_cols = read.table(paste(args$home,args$data,args$sample_id,'ligands.csv.gz',sep=''),header=TRUE)
f = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"hmap_clust_ligands_all.pdf",sep="")
weightmat_lr_plot(beta,beta_cols,f)



