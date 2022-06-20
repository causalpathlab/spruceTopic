
# library(msigdbr)


# immunesigdb_human_db <- msigdbr::msigdbr(species = "human",
# 										category = "C7",
# 										subcategory = "IMMUNESIGDB")
# immunesigdb_human_db = as.data.frame(immunesigdb_human_db)
# write.csv(immunesigdb_human_db, file=gzfile("immunesigdb_human_db.csv.gz"))

# hallmark_human_db <- msigdbr::msigdbr(species = "human",
# 										category = "H")
# hallmark_human_db = as.data.frame(hallmark_human_db)
# write.csv(hallmark_human_db, file=gzfile("hallmark_human_db.csv.gz"))

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
args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)

weightmat_phmap_plot_i <- function(args) {



beta = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cmark_hypergeom_test.tsv.gz",sep="")

df_beta = read.table(beta, sep = "\t", header=TRUE)

df_beta = cast(df_beta,pathway~topic,value.var='pval')
df_beta = as.data.frame(df_beta)
rownames(df_beta) = df_beta$pathway
df_beta = df_beta[,2:dim(df_beta)[2]]
print(head(df_beta))
df_beta = t(df_beta)


row_order = row.order(df_beta)

df_beta_t = as.data.frame(df_beta)
df_beta_t$topic = rownames(df_beta)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
col_order = col.order(df_beta_t,row_order)

df_beta = df_beta[,col_order]
df_beta = df_beta[row_order,]

# df_beta[df_beta < 0.01] = 0
# df_beta[df_beta > 1e-5] = 1

p1 <- pheatmap(t(df_beta),colorRampPalette(c("navy", "grey","white"))(100),fontsize_row=8,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=T)


f = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cell_topic_cmark.pdf",sep="")
ggsave(f,p1)

}
weightmat_phmap_plot_i(args)
