

#########################################################

#########################################################


################# set up 
args_home ="/home/BCCRC.CA/ssubedi/projects/experiments/spruce_topic/2_cell_topic_v_beta/"
setwd(paste(args_home,'scripts/',sep=''))
library(yaml)
config = paste(args_home,"config.yaml",sep="") 
args = read_yaml(config)
ct_model_id = paste(args_home,args$cell_topic$out,args$cell_topic$model_id,sep='')
it_model_id = paste(args_home,args$interaction_topic$out,args$interaction_topic$model_id,sep='')



source('fig_1_stplot.R')
h_sample_file = paste(model_id,"_ct_h_sample_celltype.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(model_id,'_ct_struct_plot_celltype.png',sep=''))

h_sample_file = paste(model_id,"_ct_h_sample_kmeans_celltype.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(model_id,'_ct_struct_plot_kmeans.png',sep=''))

h_sample_file = paste(model_id,"_ct_h_sample_kmeans_cluster.csv.gz",sep="")
df_h = read.table(h_sample_file, sep = ",", header=TRUE)
struct_plot(df_h,paste(model_id,'_ct_struct_plot_kmeans_cluster.png',sep=''))



source('../scripts/fig_1_hmap.R')
topgenes_file = paste(model_id,"_cell_topic_top_5_genes_topic.tsv.gz",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

beta = paste(model_id,"_cell_topic_beta.tsv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'genes.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = "\t", header=TRUE)
colnames(df_beta) = beta_cols$X0


top_genes=FALSE
weightmat_plot(df_beta,top_genes,df_tg)


top_genes=TRUE
weightmat_plot(df_beta,top_genes,df_tg)


#########################################################
source('fig_2_lr_loss.R')

plot_loss(paste(model_id,'_it_model_lossm.txt.gz',sep=''),paste(model_id,'_it_model_loss.png',sep=''))



#########################################################
source('fig_2_lr_hmap_v2.R')
topgenes_file = paste(model_id,"_it_beta_weight_top_5_genes.csv.gz",sep="")

beta = paste(model_id,"_it_beta_l.tsv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'receptors.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = "\t", header=TRUE)
colnames(df_beta) = beta_cols$X0

f = paste(model_id,"_it_beta_r_hmap_tp.png",sep="")
tag='receptors'
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

f = paste(model_id,"_it_beta_r_hmap_all.png",sep="")
top_genes=FALSE

#########################################################

beta = paste(model_id,"_it_beta_r.tsv.gz",sep="")
beta_cols = read.table(paste(args_home,args$data,args$sample_id,'ligands.csv.gz',sep=''),header=TRUE)
df_beta = read.table(beta, sep = "\t", header=TRUE)
colnames(df_beta) = beta_cols$X0

f = paste(model_id,"_it_beta_l_hmap_tp.png",sep="")
tag='ligands'
top_genes=TRUE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

f = paste(model_id,"_it_beta_l_hmap_all.png",sep="")
top_genes=FALSE
weightmat_lr_plot(df_beta,tag,top_genes,topgenes_file,f)

#########################################################

source('fig_2_lr_hmap_tpwise.R')
lr_pair_file = paste(it_model_id,'_it_beta_top_25_lrpair.csv.gz',sep="")
df_lrpair = read.table(lr_pair_file,sep=',', header=TRUE)
f = paste(it_model_id,"_it_beta_lr_pair_hmap_25.png",sep="")
weightmat_plot_lrpair(df_lrpair,f)