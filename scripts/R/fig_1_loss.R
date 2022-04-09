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


loss_plot <- function(args) {
loss_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"loss2.txt",sep="")
df = read.table(loss_file, sep = ",", header=TRUE)
colnames(df) = c("Log-likelihood","KL loss")
df$epoch <- 1:nrow(df)
dfm = melt(df,id="epoch")
 
p1 <-
  .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
    geom_point(stroke = 0, color="gray", size=1) +
    geom_smooth(color="red", se=FALSE, size=1) +
    labs(x = "Optimization setp", title = "", y = "Log-likelihood")
  

p2 <-
  .gg.plot(dfm[dfm$variable=="KL loss",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization setp", title = "", y = "KL loss")

plotlist = list()
plotlist[[1]] = p1 
plotlist[[2]] = p2 

stplt <- grid.arrange(grobs=plotlist,ncol=2,
heights = c(1/2, 1/2))
ggsave("loss_plot.pdf",stplt)
}

loss_plot(args)