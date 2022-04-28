library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
library(data.table)
setwd(box::file())
source("Util.R")

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)


f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"etm_beta1_bias_data_v2.tsv",sep="")
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

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"etm_beta2_bias_data_v2.tsv",sep="")
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
