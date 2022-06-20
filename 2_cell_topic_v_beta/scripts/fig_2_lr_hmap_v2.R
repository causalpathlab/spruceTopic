library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(pheatmap)
library(dplyr)
setwd(box::file())
source("Util.R")



weightmat_lr_plot <- function(df_beta,tag,p_top_genes,topgenes_file,f) {

if (p_top_genes){
df_tg = read.table(topgenes_file, sep = ",", header=TRUE)
df_tg = df_tg[df_tg$GeneType==tag,]
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

if(p_top_genes){
p1 <- pheatmap(df_beta,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),fontsize_row=8,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=T)
}
else{
p1 <- pheatmap(df_beta,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),fontsize_row=8,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)
}

ggsave(f,p1)
}




