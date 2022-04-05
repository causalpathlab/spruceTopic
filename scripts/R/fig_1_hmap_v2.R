library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
source("Util.R")


config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"


weightmat_phmap_plot_i <- function(args) {
topgenes_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"top_5_genes_topic.tsv",sep="")
df_tg = read.table(topgenes_file, sep = "\t", header=TRUE)

# df_tg <- df_tg[order(df_tg$Genes),]
df_tg$Topic <- factor(df_tg$Topic,levels = c("k0", "k1", "k2", "k3", "k4","k5",
      "k6","k7","k8","k9","k10","k11"))
# df_tg$Proportion = log(df_tg$Proportion)
# df_topic_i = as.matrix(cast( df_tg[df_tg$GeneType=="immune",] , Topic~Gene) )
# df_topic_ni = as.matrix(cast( df_tg[df_tg$GeneType=="non-immune",] , Topic~Gene) )

p1 <- .gg.plot(df_tg[df_tg$GeneType=="immune",], aes(x = Gene, y = Topic, fill = Proportion)) +
        theme(axis.text.x = element_text(angle=70, vjust=1, hjust=1, size=6)) +
        theme(legend.position = "none") +
        xlab("Top genes") + ylab("Topics") +
        geom_tile(colour = "black", size = .1) +
        scale_fill_distiller("",palette = "PuRd", direction = 1,trans="sqrt")

ggsave("topic_imm_v2.pdf",p1)

}

weightmat_phmap_plot_i(args)
