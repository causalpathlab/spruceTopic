library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
source("/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/scripts/Util.R")

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/pbmc.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"


loss_plot <- function(args) {
loss_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"loss.txt",sep="")
df = read.table(loss_file, sep = ";", header=TRUE)

df = df[seq(1, nrow(df), 100), ]

colnames(df) = c("Log-likelihood","KL loss ligands","KL loss receptors")
df$epoch <- 1:nrow(df)
dfm = melt(df,id="epoch")
 
p1 <-
  .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
    geom_point(stroke = 0, color="gray", size=1) +
    geom_smooth(color="red", se=FALSE, size=1) +
    labs(x = "Optimization setp", title = "", y = "Log-likelihood")
  

p2 <-
  .gg.plot(dfm[dfm$variable=="KL loss ligands",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization setp", title = "", y = "KL loss ligands")

p3 <-
  .gg.plot(dfm[dfm$variable=="KL loss receptors",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization setp", title = "", y = "KL loss receptors")

plotlist = list()
plotlist[[1]] = p1 
plotlist[[2]] = p2 
plotlist[[3]] = p3 

stplt <- grid.arrange(grobs=plotlist,ncol=3,
heights = c(1/3, 1/3,1/3))

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"loss_plot.pdf",sep="")
ggsave(f,stplt)
}

loss_plot(args)