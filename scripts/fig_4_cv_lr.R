library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
setwd(box::file())
source("Util.R")

ctype = commandArgs(trailingOnly=TRUE)[1]
args_home ="/home/BCCRC.CA/ssubedi/projects/spruce_topic/"
config = paste(args_home,"/config/",ctype,".yaml",sep="") 
args = read_yaml(config)

summary_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"_cv_ligands.tsv.gz",sep="")


df = read.table(summary_file, sep = "\t", header=TRUE)

df = df[df$state==19,]
df = df[,3:dim(df)[2]]
df = sqrt(df)

mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df,color = mat_colors,legend =FALSE,fontsize_row=2,fontsize_col=4)
f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"_cv_ligands_s19.pdf",sep="")
ggsave(f,p1)


summary_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"_cv_receptors.tsv.gz",sep="")


df = read.table(summary_file, sep = "\t", header=TRUE)

df = df[df$state==19,]
df = df[,3:dim(df)[2]]
df = sqrt(df)


mat_colors <- colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)

p1 <- pheatmap(df,color = mat_colors,legend =FALSE,fontsize_row=2,fontsize_col=4)
f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"_cv_receptors_s19.pdf",sep="")
ggsave(f,p1)
