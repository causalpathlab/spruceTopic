library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
# library(pheatmap)
library(dplyr)
source("Util.R")
options(repr.plot.width = 15, repr.plot.height = 7, repr.plot.res = 300)


weightmat_plot <- function(df_beta,p_top_genes,df_tg,f) {

if (p_top_genes){
top_genes = unique(df_tg$Gene)
df_beta = select(df_beta,all_of(top_genes))
}

row_order = row.order(df_beta)

df_beta_t = df_beta
df_beta_t$topic = rownames(df_beta)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
col_order = col.order(df_beta_t,row_order)

df_beta = df_beta[,col_order]
df_beta = df_beta[row_order,]


df_beta[df_beta < -20] = -20
df_beta[df_beta > 20] = 20


ignore = c("26","14","45","43","16","4","42","48","17","19","50","30","37","11","9","6")

df_beta = df_beta[!(rownames(df_beta) %in% ignore ),]
if(p_top_genes){
p1 <- pheatmap(df_beta,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),fontsize_row=8,fontsize_col=3,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)
}
else{
p1 <- pheatmap(df_beta,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),fontsize_row=8,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)
}
ggsave(f,p1)
}

# args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
# config = paste(args_home,"config.yaml",sep="") 
# args = read_yaml(config)
# model_id = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,sep='')

# topgenes_file = paste(model_id,"_cell_topic_top_5_genes_topic.tsv.gz",sep="")
# df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

# beta = paste(model_id,"_cell_topic_beta.tsv.gz",sep="")
# beta_cols = read.table(paste(args_home,args$data,args$sample_id,'genes.csv.gz',sep=''),header=TRUE)
# df_beta = read.table(beta, sep = "\t", header=TRUE)
# colnames(df_beta) = beta_cols$X0


# top_genes=TRUE
# weightmat_lr_plot(df_beta,top_genes,topgenes_file)
