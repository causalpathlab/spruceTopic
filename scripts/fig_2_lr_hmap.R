library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
source("/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/scripts/Util.R")

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/tcell.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"

weightmat_lr_plot <- function(args) {
topgenes_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"top_5_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

df_tg$Topic <- gsub("k","",df_tg$Topic)

df_tg$Topic <- factor(df_tg$Topic,levels = factor(seq(0,49)))

p1 <- .gg.plot(df_tg[df_tg$GeneType=="receptors",], aes(x = Gene, y = Topic, fill = Proportion)) +
        theme(axis.text.x = element_text(angle=70, vjust=1, hjust=1, size=4)) +
        theme(axis.text.y = element_text( vjust=1, hjust=1, size=5)) +
        theme(legend.position = "none") +
        xlab("Top genes - receptors") + ylab("Topic") +
        geom_tile(colour = "black", size = .1) +
        scale_fill_distiller("",palette = "PuRd", direction = 1,trans="sqrt")

p2 <- .gg.plot(df_tg[df_tg$GeneType=="ligands",], aes(x = Gene, y = Topic, fill = Proportion)) +
        theme(axis.text.x = element_text(angle=70, vjust=1, hjust=1, size=6)) +
        theme(axis.text.y = element_text( vjust=1, hjust=1, size=5)) +
        theme(legend.position = "none") +
        xlab("Top genes - ligands ") + ylab("Topic") +
        geom_tile(colour = "black", size = .1) +
        scale_fill_distiller("",palette = "PuRd", direction = 1,trans="sqrt")

plotlist = list()
plotlist[[1]] = p2 
plotlist[[2]] = p1 

stplt <- grid.arrange(grobs=plotlist,ncol=1,
heights = c(1/2, 1/2))
f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"pp_weightmat_plot_lr.pdf",sep="")
ggsave(f,stplt)
}
weightmat_lr_plot(args)



