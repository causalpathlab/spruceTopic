library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"

weightmat_lr_plot <- function(args) {
topgenes_file = paste(args$home,args$output,args$lr_model$out,"128_128_1000_128_ep300_lrnet_top_5_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

df_tg$Topic <- factor(df_tg$Topic)

df_tg$Proportion = log(df_tg$Proportion)
p1 <- ggplot(df_tg[df_tg$GeneType=="receptors",], aes(Gene,Topic)) + 
     geom_tile(aes(fill = Proportion), colour = "white") + 
#      geom_text(aes(label=Gene),size=2) +
     scale_fill_gradient(low = "white", high = "steelblue") + 
     xlab("Receptors") +
     theme(panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
      axis.text.x = element_text(size=4,face="bold",angle = 45,hjust=1),
      axis.text.y = element_text(size=4))

p2 <- ggplot(df_tg[df_tg$GeneType=="ligands",], aes( Gene, Topic)) + 
     geom_tile(aes(fill = Proportion), colour = "white") + 
#      geom_text(aes(label=Gene),size=2) +
     scale_fill_gradient(low = "white", high = "steelblue") +
     xlab("Ligands") +
      theme(panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
      axis.text.x = element_text(size=4,face="bold",angle = 45,hjust=1),
      axis.text.y = element_text(size=4))

plotlist = list()
plotlist[[1]] = p1 
plotlist[[2]] = p2 

stplt <- grid.arrange(grobs=plotlist,ncol=1,
heights = c(1/2, 1/2))
ggsave("pp_weightmat_plot_lr.pdf",stplt)
}

weightmat_lr_phmap_plot <- function(args) {

topgenes_file = paste(args_home,args$output,args$lr_model$out,"128_128_1000_128_ep300_lrnet_top_5_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

df_tg$Topic <- factor(df_tg$Topic)

df_tg$Proportion = log(df_tg$Proportion)

df_topic_l = as.matrix(cast( df_tg[df_tg$GeneType=="ligands",] , Topic~Gene) )
df_topic_r = as.matrix(cast( df_tg[df_tg$GeneType=="receptors",] , Topic~Gene) )

p1 <- pheatmap(df_topic_l,fontsize_row=6,fontsize_col=6)
p2 <- pheatmap(df_topic_r,fontsize_row=6,fontsize_col=6)

ggsave("topic_l.pdf",p1)
ggsave("topic_r.pdf",p2)

}


# weightmat_lr_plot(args)
# struct_plot(args)
weightmat_phmap_plot(args)
