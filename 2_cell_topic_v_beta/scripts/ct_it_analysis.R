#########################################################
################################ set up 
args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
setwd(paste(args_home,'scripts/',sep=''))
library(yaml)
config = paste(args_home,"config.yaml",sep="") 
args = read_yaml(config)
ct_model_id = paste(args_home,args$cell_topic$out,args$cell_topic$model_id,sep='')
it_model_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$model_id,sep='')



source('fig_1_stplot.R')
h_sample_file = paste(ct_model_id,"_ct_h_sample_celltype.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(ct_model_id,'_ct_struct_plot_celltype.png',sep=''))

h_sample_file = paste(ct_model_id,"_ct_h_sample_kmeans_celltype.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(ct_model_id,'_ct_struct_plot_kmeans.png',sep=''))

h_sample_file = paste(ct_model_id,"_ct_h_sample_kmeans_cluster.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(ct_model_id,'_ct_struct_plot_kmeans_cluster.png',sep=''))



source('../scripts/fig_1_hmap.R')
topgenes_file = paste(ct_model_id,"_ct_beta_weight_top_5_genes.csv.gz",sep="")
df_tg = read.table(topgenes_file, sep = ",", header=TRUE)

beta = paste(ct_model_id,"_ct_beta_mean.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'genes.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0

# f = paste(ct_model_id,"_ct_beta_hmap_all.png",sep="")
# top_genes=FALSE
# weightmat_plot(df_beta,top_genes,df_tg,f)

f = paste(ct_model_id,"_ct_beta_hmap_tp.png",sep="")
top_genes=TRUE
weightmat_plot(df_beta,top_genes,df_tg,f)



#########################################################
##############    interaction topic
#########################################################


#########################################################
source('fig_2_lr_loss.R')

plot_loss(paste(it_model_id,'_it_model_lossm.txt.gz',sep=''),paste(it_model_id,'_it_model_loss.png',sep=''))



#########################################################
source('fig_2_lr_hmap_v2.R')
topgenes_file = paste(it_model_id,"_it_beta_weight_top_10_genes.csv.gz",sep="")

beta = paste(it_model_id,"_it_beta_l.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'receptors.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0

f = paste(it_model_id,"_it_beta_r_hmap_tp.png",sep="")
tag='receptors'
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

f = paste(it_model_id,"_it_beta_r_hmap_all.png",sep="")
top_genes=FALSE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

#########################################################

beta = paste(it_model_id,"_it_beta_r.csv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'ligands.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = ",", header=TRUE)
colnames(df_beta) = beta_cols$X0

f = paste(it_model_id,"_it_beta_l_hmap_tp.png",sep="")
tag='ligands'
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

f = paste(it_model_id,"_it_beta_l_hmap_all.png",sep="")
top_genes=FALSE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

#########################################################

source('fig_2_lr_hmap_tpwise.R')
lr_pair_file = paste(it_model_id,'_it_beta_top_25_lrpair.csv.gz',sep="")
df_lrpair = read.table(lr_pair_file,sep=',', header=TRUE)
f = paste(it_model_id,"_it_beta_lr_pair_hmap.png",sep="")
weightmat_plot_lrpair(df_lrpair,f)


#########################################################

source('fig_3_summary.R')
summary_file = paste(it_model_id,'_it_model_summary.csv.gz',sep="")

source('fig_3_summary.R')
df_summary = read.table(summary_file,sep=',', header=TRUE)
df_cancer = df_summary[df_summary$celltype %in% c('Cancer Epithelial'), ]
f = paste(it_model_id,"_it_ccview_summary_cancer.png",sep="")
col=1
summary_plot_all(df_cancer,f,col)

# source('fig_3_summary.R')
# df_summary = read.table(summary_file,sep=',', header=TRUE)
# df_others = df_summary[df_summary$celltype %in% c('PVL','Endothelial','T-cells','B-cells','Myeloid','CAFs','Plasmablasts','pan-CD8'), ]
# f = paste(it_model_id,"_it_model_summary_others.png",sep="")
# col=4
# summary_plot_all(df_others,f,col)

# source('fig_3_summary.R')
# df_summary = read.table(summary_file,sep=',', header=TRUE)
# df_normal = df_summary[df_summary$celltype %in% c('n_B-cells','n_Endothelial',
# 'n_Epithelial','n_Fibroblasts','n_Others','n_T-cells','Normal Epithelial'), ]
# f = paste(it_model_id,"_it_model_summary_normal.png",sep="")
# col=4
# summary_plot_all(df_normal,f,col)



#########################################################

source('fig_4_ccview.R')
df = read.table(paste(it_model_id,'_it_cancercells_interactions.csv.gz',sep=""),sep=',', header=TRUE)
f = paste(it_model_id,"_it_ccview_stplot.png",sep="")
ccv_struct_plot(df,f)

# source('fig_4_ccview.R')
# df = read.table(paste(it_model_id,'_it_cancer_normal_cl.csv.gz',sep=""),sep=',', header=TRUE)
# f = paste(it_model_id,"_it_ccview_cl.png",sep="")
# cancer_nbr_lr_plot(df,f)