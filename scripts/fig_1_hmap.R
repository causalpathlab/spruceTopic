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

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/spruce_topic/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)

weightmat_phmap_plot_i <- function(args) {
topgenes_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"_netm_top_5_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

# df_tg <- df_tg[order(df_tg$Genes),]
# df_tg$Topic <- factor(df_tg$Topic)

df_tg$Topic <- gsub("k","",df_tg$Topic)
# df_tg$Topic <- factor(df_tg$Topic, levels = c("0", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14"))

# df_tg$Proportion = log(df_tg$Proportion)
# df_topic_i = as.matrix(cast( df_tg[df_tg$GeneType=="immune",] , Topic~Gene) )
# df_topic_ni = as.matrix(cast( df_tg[df_tg$GeneType=="non-immune",] , Topic~Gene) )

xlab = paste("Top","genes")
p1 <- .gg.plot(df_tg[df_tg$GeneType=="top_genes",], aes(x = Gene, y = Topic, fill = Proportion)) +
        theme(axis.text.x = element_text(angle=70, vjust=1, hjust=1, size=6)) +
        theme(legend.position = "none") +
        xlab(xlab) + ylab("Topics") +
        geom_tile(colour = "black", size = .1) +
        scale_fill_distiller("",palette = "PuRd", direction = 1,trans="sqrt")


f = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"_netm_topic_gene_weight_hmap.pdf",sep="")

plotlist = list()
plotlist[[1]] = p1 

stplt <- grid.arrange(grobs=plotlist,ncol=1,
heights = c(1))

ggsave(f,stplt)

}

weightmat_phmap_plot_i(args)
