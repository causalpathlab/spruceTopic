library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
source("Util.R")



weightmat_plot_lrpair <- function(df_tg,f) {

df_tg = as.data.frame(df_tg)
rownames(df_tg) = df_tg$X
df_tg$X = NULL

row_order = row.order(df_tg)

df_tg_t = df_tg
df_tg_t$topic = rownames(df_tg)
df_tg_t = melt(df_tg_t)
colnames(df_tg_t)=c('row','col','weight')
col_order = col.order(df_tg_t,row_order)

df_tg = df_tg[,col_order]
df_tg = df_tg[row_order,]

df_tg[df_tg < -20] = -20
df_tg[df_tg > 20] = 20

colnames(df_tg) <- gsub("X","",colnames(df_tg))

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df_tg,color =colorRampPalette(c("navy", "white", "firebrick3"))(100),legend =TRUE,fontsize_row=4,fontsize_col=12)

ggsave(f,p1)

}

# weightmat_phmap_plot_i(args)
