#########################################################
################################ set up 
args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
# args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
setwd(paste(args_home,'scripts/',sep=''))
library(yaml)
config = paste(args_home,"config.yaml",sep="") 
args = read_yaml(config)
ct_model_id = paste(args_home,args$cell_topic$out,args$cell_topic$model_id,sep='')
it_model_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$model_id,sep='')
ct_id = paste(args_home,args$cell_topic$out,args$cell_topic$id,sep='')
it_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$id,sep='')


chordDiagram <- function(){

file = paste(it_id,"8_chorddata.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)
df = df[order(df$topic),]
df$topic = as.factor(df$topic)

library(reshape)

dfm = melt(df,id=c('topic'))
dfm = dfm[,c('value','variable')]
colnames(dfm) = c('gene','genetype')

dfl = df[,c(1,2)]
library("circlize")
f = paste(it_id,"8_chorddata_ligand.pdf",sep="")
pdf(f)
# df[,c(1,3,4)]
genes<- as.character(dfl[[2]])
othercol = structure(rep("grey", length(genes)), names = genes)
grid_col = c("2" = "#b3446c", "4" = "#dcd300", "7" = "#882d17",
"10"="#8db600","18"="#654522","22"="#e25822","24"="#2b3d26", othercol )




# chordDiagram(df,,grid.col=grid_col,col = grid_col[as.character(df[[1]])],annotationTrack=NULL)


# Plot chord diagram
chordDiagram(dfl,
             grid.col=grid_col,
             col = grid_col[as.character(dfl[[1]])],
             annotationTrack = c("grid"), 
             preAllocateTracks = 1, 

             )


highlight.sector(dfl$ligand,
                 track.index = 1, col = "white",
                 text = "Ligand", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
highlight.sector(dfl$topic,
                 track.index = 1, col = "white",
                 text = "Interaction topics", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.4)
#   circos.axis(h = "top", labels.cex = 0, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  circos.axis(h = "top",labels=NULL)
}, bg.border = NA)
dev.off()
circos.clear()


dfr = df[,c(1,3)]
library("circlize")
f = paste(it_id,"8_chorddata_receptor.pdf",sep="")
pdf(f)
# df[,c(1,3,4)]
genes<- as.character(dfr[[2]])
othercol = structure(rep("grey", length(genes)), names = genes)
grid_col = c("2" = "#b3446c", "4" = "#dcd300", "7" = "#882d17",
"10"="#8db600","18"="#654522","22"="#e25822","24"="#2b3d26", othercol )




# chordDiagram(df,,grid.col=grid_col,col = grid_col[as.character(df[[1]])],annotationTrack=NULL)


# Plot chord diagram
chordDiagram(dfr,
             grid.col=grid_col,
             col = grid_col[as.character(dfr[[1]])],
             annotationTrack = c("grid"), 
             preAllocateTracks = 1, 

             )
highlight.sector(df$receptor,
                 track.index = 1, col = "white",
                 text = "Receptor", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
highlight.sector(dfr$topic,
                 track.index = 1, col = "white",
                 text = "Interaction topics", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.4)
#   circos.axis(h = "top", labels.cex = 0, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  circos.axis(labels=NULL)
}, bg.border = NA)
dev.off()
circos.clear()

}

########
source('fig_1_stplot.R')
h_sample_file = paste(ct_id,"2_i_kmeans_celltype_sample.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(ct_id,'2_i_kmeans_celltype_sample_stplot.pdf',sep=''))


source('fig_1_stplot.R')
h_sample_file = paste(ct_model_id,"_ct_3_h_kmeans_cluster_sample.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(ct_model_id,'_ct_3_struct_plot_kmeans_cluster.pdf',sep=''))
########


source('fig_1_hmap.R')
topgenes_file = paste(ct_id,"_ct_4_beta_weight_top_25_genes.csv.gz",sep="")
df_tg = read.table(topgenes_file, sep = ",", header=TRUE)

beta = paste(ct_id,"_ct_beta_mean.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'genes.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0

# f = paste(ct_model_id,"_ct_beta_hmap_all.png",sep="")
# top_genes=FALSE
# weightmat_plot(df_beta,top_genes,df_tg,f)

f = paste(ct_id,"4_beta_hmap_tp_25.pdf",sep="")
top_genes=TRUE
weightmat_plot(df_beta,top_genes,df_tg,f)



source('fig_1_hmap.R')
topgenes_file = paste(ct_id,"4_beta_weight_top_25_genes_selected.csv.gz",sep="")
df_beta = read.table(topgenes_file, sep = ",", header=TRUE,row.names=1)

top3_file = paste(ct_id,"4_beta_weight_top_25_genes_top3.csv.gz",sep="")
df_top3 = read.table(top3_file, sep = ",", header=TRUE,row.names=1)


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

f = paste(ct_id,"4_beta_hmap_tp_selected.pdf",sep="")

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


#########################################################
##############    interaction topic
#########################################################


#########################################################
source('fig_2_lr_loss.R')

plot_loss(paste(it_model_id,'_it_model_lossm.txt.gz',sep=''),paste(it_model_id,'_it_model_loss.png',sep=''))



#########################################################
source('fig_2_lr_hmap_v2.R')
topgenes_file = paste(it_id,"2_beta_weight_top_10_genes.csv.gz",sep="")
beta = paste(it_model_id,"_it_beta_lm.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'receptors.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0
df_beta = as.data.frame(df_beta)
rownames(df_beta) = 0:24
df_beta = df_beta[c("2","4","7", "10","18","22","24"),]
tag='ligands'
f = paste(it_id,"2_beta_weight_top_10_receptors.pdf",sep="")
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

topgenes_file = paste(it_id,"2_beta_weight_top_10_genes.csv.gz",sep="")
beta = paste(it_model_id,"_it_beta_rm.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'ligands.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0
df_beta = as.data.frame(df_beta)
rownames(df_beta) = 0:24
df_beta = df_beta[c("2","4","7", "10","18","22","24"),]
tag='receptors'
f = paste(it_id,"2_beta_weight_top_10_ligands.pdf",sep="")
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)


# f = paste(it_model_id,"_it_beta_r_hmap_all_v2.png",sep="")
# top_genes=FALSE
# weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

#########################################################

beta = paste(it_model_id,"_it_beta_r.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'ligands.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0

rownames(df_beta) = 0:24

df_beta = df_beta[c("2","4","7", "10","18","22","24"),]

f = paste(it_model_id,"_it_beta_l_hmap_tp.pdf",sep="")
tag='receptors'
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

# f = paste(it_model_id,"_it_beta_l_hmap_all.png",sep="")
# top_genes=FALSE
# weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

#########################################################

# source('fig_2_lr_hmap_tpwise.R')
# lr_pair_file = paste(it_model_id,'_it_beta_top_25_lrpair.csv.gz',sep="")
# df_lrpair = read.table(lr_pair_file,sep=',', header=TRUE)

# # df_lrpair = df_lrpair[,c("X","X1","X2","X4","X7", "X10","X18","X22","X24")]

# f = paste(it_model_id,"_it_beta_lr_pair_hmap.png",sep="")
# weightmat_plot_lrpair(df_lrpair,f)


#########################################################

# source('fig_3_summary.R')
# summary_file = paste(it_model_id,'_it_model_summary.csv.gz',sep="")
# df_summary = read.table(summary_file,sep=',', header=TRUE)
# cancer_cells = c('Cancer_Basal_SC','Cancer_Cycling','Cancer_Her2_SC','Cancer_LumA_SC','Cancer_LumB_SC')
# df_cancer = df_summary[df_summary$celltype %in% cancer_cells, ]
# f = paste(it_model_id,"_it_ccview_summary_cancer.png",sep="")
# col=1
# summary_plot_all(df_cancer,f,col)

source('fig_3_summary.R')
summary_file = paste(it_id,'4_cell_nbr_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_it_summary.pdf',sep="")
col=9
tag='celltype'
summary_plot_cancer_v1(df_summary,f,tag,col)    

source('fig_3_summary.R')
summary_file = paste(it_id,'4_cell_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_ct_it_summary.pdf',sep="")
col=3
tag='celltype'
df_summary <- df_summary[order(df_summary$celltype),]
summary_plot_cancer(df_summary,f,tag,col)    

source('fig_3_summary.R')
summary_file = paste(it_id,'4_cell_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_ct_it_summary_boxplt.pdf',sep="")
summary_boxplot(df_summary,f)

source('fig_3_summary.R')
summary_file = paste(it_id,'5_cancer_stype_nbr_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'5_cancer_stype_nbr_it_summary.pdf',sep="")
col=4
tag='celltype'
summary_plot_cancer_v1(df_summary,f,tag,col)    

summary_file = paste(it_id,'5_cancer_stype_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'5_cancer_stype_nbr_ct_it_summary.pdf',sep="")
col=1
tag='celltype'
summary_plot_cancer(df_summary,f,tag,col)    

source('fig_3_summary.R')
summary_file = paste(it_id,'5_cancer_stype_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'5_cancer_stype_nbr_ct_it_summary_boxplot.pdf',sep="")
summary_boxplot(df_summary,f)

#### normal vs cancer dataset 

source('fig_3_summary.R')
summary_file = paste(it_model_id,'_it_3_cell_nbr_it_norm_cancer_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_model_id,'_it_3_cell_nbr_it_norm_cancer_summary.pdf',sep="")
col=1
tag='interact_topic_x'
df_summary$cells = 'Cancer neighbours'
summary_plot_cancer_v2(df_summary,f,tag,col)    

source('fig_3_summary.R')
summary_file = paste(it_model_id,'_it_3_cell_nbr_ct_norm_cancer_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_model_id,'_it_3_cell_nbr_ct_norm_cancer_summary.pdf',sep="")
col=1
tag='cell_topic'
df_summary$cells = 'Cells'
summary_plot_cancer_v2(df_summary,f,tag,col)    


#########################################################

source('fig_4_ccview.R')
df = read.table(paste(it_id,'4_it_celltypedist.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_id,"4_it_celltypedist.pdf",sep="")
ccv_struct_plot_v2(df,f)

source('fig_4_ccview.R')
df = read.table(paste(it_model_id,'_it_4_cancercells_it_cluster_celltypedist_normal.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_model_id,"_it_4_cancercells_it_cluster_celltypedist_normal.pdf",sep="")
ccv_struct_plot(df,f)

df = read.table(paste(it_model_id,'_it_4_cancercells_it_cluster_celltypedist_notnormal.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_model_id,"_it_4_cancercells_it_cluster_celltypedist_notnormal.pdf",sep="")
ccv_struct_plot(df,f)

source('fig_4_ccview.R')
df = read.table(paste(it_model_id,'_it_4_cancercells_it_cluster_subtypedist.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_model_id,"_it_4_cancercells_it_cluster_subtypedist.pdf",sep="")
ccv_struct_plot_subtye(df,f)




source('fig_1_stplot.R')
# h_sample_file = paste(ct_model_id,"_ct_h_sample_celltype.csv.gz",sep="")
# df_h = read.table(h_sample_file, sep = ",", header=TRUE)
# struct_plot(df_h,paste(ct_model_id,'_ct_struct_plot_celltype.png',sep=''))

# h_sample_file = paste(it_model_id,"_it_h_sample_kmeans.csv.gz",sep="")
# df_h = read.table(h_sample_file, sep = ",", header=TRUE)
# struct_plot(df_h,paste(it_model_id,'_it_struct_plot_kmeans.png',sep=''))

h_sample_file = paste(it_id,"6_topic_sample_kmeans_stplot.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot_interaction_topic(df_h,paste(it_id,'_it_6_topic_sample_kmeans_stplot.png',sep=''))


