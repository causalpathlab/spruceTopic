library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(pheatmap)
library(viridis)
library(dplyr)
source("Util.R")
options(repr.plot.width = 15, repr.plot.height = 7, repr.plot.res = 300)


beta_loadings_plot <- function(df_beta,p_top_genes,df_tg,f) {

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
p1 <- pheatmap(df_beta,color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),fontsize_row=8,fontsize_col=3,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)
}
else{
p1 <- pheatmap(df_beta,color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),fontsize_row=8,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)
}
ggsave(f,p1)
}



beta_loadings_plot_with_toplabels <- function(df_beta,df_top3,f) {

df_beta = as.data.frame(df_beta)
row_order = row.order(df_beta)
df_beta = df_beta[row_order,]

df_beta_t = df_beta
df_beta_t$topic = rownames(df_beta)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
df_beta_t = as.data.frame(df_beta_t)
col_order = col.order(df_beta_t,rownames(df_beta)) 
df_beta = df_beta[,col_order]

df_beta[ df_beta < -10] = -10
df_beta[ df_beta > 10] = 10

library(ComplexHeatmap)
library(viridis)



pdf(f,width=12,height=10)

ha = rowAnnotation(Gene = anno_text(df_top3$gene,
just = "center", 
location = unit(0.5, "npc")), 
annotation_name_rot = 0)

Heatmap(as.matrix(df_beta),
cluster_rows=FALSE,cluster_columns=FALSE,
name="Weight",
show_column_names=FALSE,
left_annotation = ha,
col = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis")
)

dev.off()
}

beta_loadings_plot_interaction_topic <- function(df_beta,tag,p_top_genes,topgenes_file,f) {

if (p_top_genes){
df_tg = read.table(topgenes_file, sep = ",", header=TRUE)
df_tg = df_tg[df_tg$GeneType==tag,]
top_genes = unique(df_tg$Gene)
print(top_genes)
df_beta = select(df_beta,all_of(top_genes))
}

# row_order = row.order(df_beta)
# df_beta = df_beta[row_order,]

df_beta_t = df_beta
df_beta_t$topic = rownames(df_beta)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
col_order = col.order(df_beta_t,rownames(df_beta))

df_beta = df_beta[,col_order]


df_beta[df_beta < -10] = -10
df_beta[df_beta > 10] = 10

if(p_top_genes){
p1 <- pheatmap(t(df_beta),color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),fontsize_row=6,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=T)
ggsave(f,p1,width = 2, height = 6)

}
else{
# p1 <- pheatmap(df_beta,color = colorRampPalette(c("blue", "white", "red"))(20),fontsize_row=6,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)
p1 <- pheatmap(df_beta,color =viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),fontsize_row=6,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=F)

ggsave(f,p1,width = 6, height = 2)

}


}



######### archive


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

df_lrpair[df_lrpair < -10] = -10
df_lrpair[df_lrpair > 10] = 10

# df_lrpair = sqrt(df_lrpair)

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(t(df_lrpair),color =colorRampPalette(c("navy", "white", "firebrick3"))(100),legend =TRUE,fontsize_row=2,fontsize_col=5,cluster_rows=FALSE,cluster_cols=FALSE)

ggsave(f,p1)

}


bias_hmap_plot <- function() {
args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/spruce_topic/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)


f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"_ietm_beta1_bias_v2.tsv",sep="")
df = read.table(f, sep = "\t", header=TRUE)
p1 <- ggplot(df[df$group=='l25',], aes(x=group, y=gene,fill = val)) +
    geom_tile() +
    scale_fill_distiller("",palette = "PuRd", direction = 1)+
    theme(
    strip.background = element_rect(fill='transparent'),
    panel.background = element_rect(fill='transparent', color=NA),
    plot.background = element_rect(fill='transparent', color=NA))
p2 <- ggplot(df[df$group=='h25',], aes(x=group, y=gene,fill = val)) + 
    geom_tile() +
    scale_fill_distiller("",palette = "PuRd", direction = 1) +
    theme(
    strip.background = element_rect(fill='transparent'),
    panel.background = element_rect(fill='transparent', color=NA),
    plot.background = element_rect(fill='transparent', color=NA))
plotlist = list()
plotlist[[1]] = p2 
plotlist[[2]] = p1 
stplt <- grid.arrange(grobs=plotlist,ncol=2,heights = c(1/2, 1/2))
f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"beta1_bias.pdf",sep="")
ggsave(f,stplt)

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"_ietm_beta2_bias_v2.tsv",sep="")
df = read.table(f, sep = "\t", header=TRUE)
p1 <- ggplot(df[df$group=='l25',], aes(x=group, y=gene,fill = val)) +
    geom_tile() +
    scale_fill_distiller("",palette = "PuRd", direction = 1)+
    theme(
    strip.background = element_rect(fill='transparent'),
    panel.background = element_rect(fill='transparent', color=NA),
    plot.background = element_rect(fill='transparent', color=NA))
p2 <- ggplot(df[df$group=='h25',], aes(x=group, y=gene,fill = val)) + 
    geom_tile() +
    scale_fill_distiller("",palette = "PuRd", direction = 1) +
    theme(
    strip.background = element_rect(fill='transparent'),
    panel.background = element_rect(fill='transparent', color=NA),
    plot.background = element_rect(fill='transparent', color=NA))
plotlist = list()
plotlist[[1]] = p2 
plotlist[[2]] = p1 
stplt <- grid.arrange(grobs=plotlist,ncol=2,heights = c(1/2, 1/2))
f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"beta2_bias.pdf",sep="")
ggsave(f,stplt)

}