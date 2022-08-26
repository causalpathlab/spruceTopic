################################ set up ##################
args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
# args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
setwd(paste(args_home,'scripts/',sep=''))
library(yaml)
options(rlib_downstream_check = FALSE)
config = paste(args_home,"config.yaml",sep="") 
args = read_yaml(config)
ct_model_id = paste(args_home,args$cell_topic$out,args$cell_topic$model_id,sep='')
it_model_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$model_id,sep='')
ct_id = paste(args_home,args$cell_topic$out,args$cell_topic$id,sep='')
it_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$id,sep='')



##### gsea 
source('fig_1_gsea.R')
file = paste(it_id,"12_gsea.csv",sep="")
df = read.csv(file, sep = ",", header=TRUE,check.names=FALSE)
run_gsea(df)



subtype_heatmap <- function(){
file = paste(it_id,"5_subtype_heatmap_r.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

dfm = melt(df,id=c('Cancer', 'subtype', 'cell_topic' ))

dfm$value = ifelse(dfm$value > 3.0,3.0,dfm$value)

p <-
ggplot(dfm, aes(x=Cancer, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred")+
  facet_grid(~ cell_topic + subtype, scales = "free", switch = "x", space = "free")+
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      text = element_text(size=2),
        panel.margin=unit(.01, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))
    #   axis.text.y=element_blank(),
    #   axis.ticks.y=element_blank())

ggsave(paste(it_id,"5_subtype_heatmap_r.pdf",sep=""),p,width = 8, height = 6,dpi=600)

file = paste(it_id,"5_subtype_heatmap_l.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

dfm = melt(df,id=c('Cancer', 'subtype', 'cell_topic' ))

dfm$value = ifelse(dfm$value > 3.0,3.0,dfm$value)

p <-
ggplot(dfm, aes(x=Cancer, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred")+
  facet_grid(~ cell_topic + subtype, scales = "free", switch = "x", space = "free")+
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      text = element_text(size=2),
        panel.margin=unit(.01, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))
    #   axis.text.y=element_blank(),
    #   axis.ticks.y=element_blank())

ggsave(paste(it_id,"5_subtype_heatmap_l.pdf",sep=""),p,width = 8, height = 6,dpi=600)

}
#### bipartite 
#https://rdrr.io/cran/bipartite/man/plotweb.html
library(bipartite)
library(ggplot2)
library(gridExtra)
library(reshape)

file = paste(it_id,"8_lrnetwork.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

for ( i in c(2,4,7,10,18,22,24)){
print(i)
dfm = df[df$topic==i,]
if (dim(dfm)[1]>100){dfm = dfm[c(1:100),]}
print(max(dfm$score))
print(min(dfm$score))
matrix = table(dfm[,c(1,2)])
edge_colors = colorRampPalette(c("grey80","tan2"))(length(dfm$score))
f=paste(it_id,"8_lrnetwork_",i,".pdf",sep="")
pdf(f)
p <- plotweb(matrix, method="normal", 
labsize=1,text.rot = 90,
col.high='salmon',col.low='cyan',
col.interaction = edge_colors,
bor.col.high='black',bor.col.low='black',bor.col.interaction='NA',
y.width.low=0.05,y.width.high=0.05,
arrow='no',ybig=0.1)
dev.off()
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


source('fig_2_lr_hmap_v2.R')
beta = paste(it_model_id,"_it_beta_lm.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'receptors.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0
df_beta = as.data.frame(df_beta)
rownames(df_beta) = 0:24
df_beta = df_beta[c("2","4","7", "10","18","22","24"),]
tag='ligands'
f = paste(it_id,"2_beta_weight_all_receptors.pdf",sep="")
top_genes=FALSE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

beta = paste(it_model_id,"_it_beta_rm.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'ligands.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0
df_beta = as.data.frame(df_beta)
rownames(df_beta) = 0:24
df_beta = df_beta[c("2","4","7", "10","18","22","24"),]
tag='receptors'
f = paste(it_id,"2_beta_weight_all_ligands.pdf",sep="")
top_genes=FALSE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

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





#### fgsea
library('fgsea')
data(examplePathways)
df_file = paste(it_id,"_gse.csv.gz",sep="")
df = read.table(df_file, sep = ",", header=TRUE)

er = df[1,]