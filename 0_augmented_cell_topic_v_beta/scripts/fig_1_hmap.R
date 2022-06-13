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
args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/0_augmented_cell_topic_v_beta/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)

weightmat_phmap_plot_i <- function(args) {

# topgenes_file = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cell_topic_top_5_genes_topic.tsv.gz",sep="")
# df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)
# top_genes = unique(df_tg$Gene)



beta = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cell_topic_beta.tsv.gz",sep="")
beta_cols = read.table(paste(args$home,args$data,args$sample_id,'genes.csv.gz',sep=''),header=TRUE)

df_beta = read.table(beta, sep = "\t", header=TRUE)

colnames(df_beta) = beta_cols$X0

# df_beta = select(df_beta,all_of(top_genes))
# print(colnames(df_beta))

row_order = row.order(df_beta)

df_beta_t = df_beta
df_beta_t$topic = rownames(df_beta)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
col_order = col.order(df_beta_t,row_order)

df_beta = df_beta[,col_order]
df_beta = df_beta[row_order,]

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "Spectral"))(100)
p1 <- pheatmap(df_beta,color = mat_colors,fontsize_row=2,fontsize_col=4,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)


f = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cell_topic_gene_weight_hmap_all_genes.pdf",sep="")
ggsave(f,p1)

}
weightmat_phmap_plot_i(args)
