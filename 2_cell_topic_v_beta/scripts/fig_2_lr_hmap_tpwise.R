library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(pheatmap)
# library(dendextend)
library(RColorBrewer)
source("Util.R")



weightmat_plot_lrpair <- function(df_lrpair,f) {

df_lrpair = as.data.frame(df_lrpair)
rownames(df_lrpair) = df_lrpair$X
df_lrpair$X = NULL

colnames(df_lrpair) = gsub("X","",as.character(colnames(df_lrpair)))

df_lrpair = t(df_lrpair)


row_order = row.order(df_lrpair)

df_lrpair_t = df_lrpair
# df_lrpair_t$topic = rownames(df_lrpair)
df_lrpair_t = melt(df_lrpair_t)
colnames(df_lrpair_t)=c('row','col','weight')
col_order = col.order(df_lrpair_t,row_order)

df_lrpair = df_lrpair[,col_order]
df_lrpair = df_lrpair[row_order,]

df_lrpair[df_lrpair < -20] = -20
df_lrpair[df_lrpair > 20] = 20

# df_lrpair = sqrt(df_lrpair)

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(t(df_lrpair),color =colorRampPalette(c("navy", "white", "firebrick3"))(100),legend =TRUE,fontsize_row=2,fontsize_col=5,cluster_rows=FALSE,cluster_cols=FALSE)

ggsave(f,p1)

}

# weightmat_phmap_plot_i(args)
