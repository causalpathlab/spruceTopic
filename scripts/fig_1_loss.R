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
loss_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"loss2.txt",sep="")
df = read.table(loss_file, sep = ",", header=TRUE)
colnames(df) = c("Log-likelihood","KL loss immune","KL loss non immune")
df$epoch <- 1:nrow(df)
dfm = melt(df,id="epoch")
 
p1 <-
  .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
    geom_point(stroke = 0, color="gray", size=1) +
    geom_smooth(color="red", se=FALSE, size=1) +
    labs(x = "Optimization setp", title = "", y = "Log-likelihood")
  

p2 <-
  .gg.plot(dfm[dfm$variable=="KL loss immune",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization setp", title = "", y = "KL loss immune marker genes")

p3 <-
  .gg.plot(dfm[dfm$variable=="KL loss non immune",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization setp", title = "", y = "KL loss other genes ")

plotlist = list()
plotlist[[1]] = p1 
plotlist[[2]] = p2 
plotlist[[3]] = p3 

stplt <- grid.arrange(grobs=plotlist,ncol=3,
heights = c(1/3, 1/3,1/3))

f = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"loss_plot.pdf",sep="")
ggsave(f,stplt)
}

loss_plot(args)