################################ set up ##################
args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
# args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
setwd(paste(args_home,'r_scripts/',sep=''))
library(yaml)
options(rlib_downstream_check = FALSE)
config = paste(args_home,"config.yaml",sep="") 
args = read_yaml(config)
ct_model_id = paste(args_home,args$cell_topic$out,args$cell_topic$model_id,sep='')
it_model_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$model_id,sep='')
ct_id = paste(args_home,args$cell_topic$out,args$cell_topic$id,sep='')
it_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$id,sep='')

###############################################################################
# cell topic analysis 
###############################################################################

##loss plot
source('fig_loss_plot.R')
plot_cell_topic_loss(paste(ct_model_id,'_ct_loss2.txt.gz',sep=''),paste(ct_model_id,'_ct_model_loss.png',sep=''))



## cell topic structure plot
source('fig_struct_plot.R')
h_sample_file = paste(ct_id,"2_i_kmeans_celltype_sample.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(ct_id,'2_i_kmeans_celltype_sample_stplot.pdf',sep=''))


##cell topic heatmap plot with top genes label
source('fig_beta_hmap.R')
topgenes_file = paste(ct_id,"4_beta_weight_top_25_genes_selected.csv.gz",sep="")
df_beta = read.table(topgenes_file, sep = ",", header=TRUE,row.names=1)

top3_file = paste(ct_id,"4_beta_weight_top_25_genes_top3.csv.gz",sep="")
df_top3 = read.table(top3_file, sep = ",", header=TRUE,row.names=1)

f = paste(ct_id,"4_beta_hmap_tp_selected.pdf",sep="")
beta_loadings_plot_with_toplabels(df_beta,df_top3,f)


###############################################################################
# interaction topic analysis 
###############################################################################


source('fig_loss_plot.R')
plot_interaction_topic_loss(paste(it_model_id,'_it_model_lossm.txt.gz',sep=''),paste(it_model_id,'_it_model_loss.png',sep=''))



#########################################################
source('fig_beta_hmap.R')
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
beta_loadings_plot_interaction_topic(df_beta,tag,top_genes,topgenes_file,f)

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
beta_loadings_plot_interaction_topic(df_beta,tag,top_genes,topgenes_file,f)


source('fig_beta_hmap.R')
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
beta_loadings_plot_interaction_topic(df_beta,tag,top_genes,topgenes_file,f)

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
beta_loadings_plot_interaction_topic(df_beta,tag,top_genes,topgenes_file,f)

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
beta_loadings_plot_interaction_topic(df_beta,tag,top_genes,topgenes_file,f)

#########################################################

source('fig_topic_summary.R')
summary_file = paste(it_id,'4_cell_nbr_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_it_summary.pdf',sep="")
col=9
tag='celltype'
summary_plot_cancer_v1(df_summary,f,tag,col)    

source('fig_topic_summary.R')
summary_file = paste(it_id,'4_cell_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_ct_it_summary.pdf',sep="")
col=3
tag='celltype'
df_summary <- df_summary[order(df_summary$celltype),]
summary_plot_cancer(df_summary,f,tag,col)    

source('fig_topic_summary.R')
summary_file = paste(it_id,'4_cell_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_ct_it_summary_boxplt.pdf',sep="")
summary_boxplot(df_summary,f)

source('fig_topic_summary.R')
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

source('fig_topic_summary.R')
summary_file = paste(it_id,'5_cancer_stype_nbr_ct_it_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'5_cancer_stype_nbr_ct_it_summary_boxplot.pdf',sep="")
summary_boxplot(df_summary,f)

#### normal vs cancer dataset 

source('fig_topic_summary.R')
summary_file = paste(it_id,'4_cell_nbr_ct_norm_cancer_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_ct_norm_cancer_summary.pdf',sep="")
col=1
tag='cell_topic'
df_summary$cells = 'cells'
summary_plot_cancer_v2(df_summary,f,tag,col)    

source('fig_topic_summary.R')
summary_file = paste(it_id,'4_cell_nbr_it_norm_cancer_summary.csv.gz',sep="")
df_summary = read.table(summary_file,sep=',', header=TRUE)
f = paste(it_id,'4_cell_nbr_it_norm_cancer_summary.pdf',sep="")
col=1
tag='interact_topic_x'
df_summary$cells = 'Cancer neighbours'
summary_plot_cancer_v2(df_summary,f,tag,col)    



#########################################################

source('fig_celltype_summary.R')
df = read.table(paste(it_id,'4_it_celltypedist.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_id,"4_it_celltypedist.pdf",sep="")
ccv_struct_plot_v2(df,f)


source('fig_celltype_summary.R')
df = read.table(paste(it_id,'4_it_celltypedist.csv.gz',sep=""),sep=',', header=TRUE)
df$interact_topic = as.factor(df$interact_topic)

library(ggplot2)
library(Polychrome)

col_vector <- as.vector(kelly.colors(22))[16:22]

p <-
ggplot(df, aes(x=cluster_celltype, y=celltype,size=ncount, fill=as.factor(interact_topic))) +
  geom_point() +
  scale_fill_manual("Cell type ",values=col_vector)+
  facet_wrap(~cluster_celltype+interact_topic)

f = paste(it_id,"4_it_celltypedist_hmap.pdf",sep="")
ggsave(f,p,width = 30, height = 20,limitsize=F)

source('fig_celltype_summary.R')
df = read.table(paste(it_model_id,'_it_4_cancercells_it_cluster_celltypedist_normal.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_model_id,"_it_4_cancercells_it_cluster_celltypedist_normal.pdf",sep="")
ccv_struct_plot(df,f)

df = read.table(paste(it_model_id,'_it_4_cancercells_it_cluster_celltypedist_notnormal.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_model_id,"_it_4_cancercells_it_cluster_celltypedist_notnormal.pdf",sep="")
ccv_struct_plot(df,f)



### interaction topic struct plot 
source('fig_struct_plot.R')
h_sample_file = paste(it_id,"6_topic_sample_kmeans_stplot.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot_interaction_topic(df_h,paste(it_id,'_it_6_topic_sample_kmeans_stplot.png',sep=''))



##### gsea 
source('fig_gsea.R')
file = paste(it_id,"12_gsea.csv",sep="")
df = read.csv(file, sep = ",", header=TRUE,check.names=FALSE)
run_gsea(df)


####others
# subtype heatmap
# chord diagram
# survial
# bipartite
