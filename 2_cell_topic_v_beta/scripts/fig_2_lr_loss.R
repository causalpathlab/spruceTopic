library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(RColorBrewer)
setwd(box::file())
source("Util.R")

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/spruce_topic/"
config = paste(args_home,"/config/",args[1],".yaml",sep="") 
args = read_yaml(config)


loss_plot <- function(args) {
loss_file = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"loss.txt.gz",sep="")
print(loss_file)
df = read.table(loss_file, sep = ";", header=TRUE)

df = df[seq(1, nrow(df), 300), ]

colnames(df) = c("Log-likelihood","KL loss ligands","KL loss receptors")
df$epoch <- 1:nrow(df)
dfm = melt(df,id="epoch")
 
p1 <-
  .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
    geom_point(stroke = 0, color="gray", size=1) +
    geom_smooth(color="red", se=FALSE, size=1) +
    labs(x = "Optimization step", title = "", y = "Log-likelihood")
  

p2 <-
  .gg.plot(dfm[dfm$variable=="KL loss ligands",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization step", title = "", y = "KL loss ligands")

p3 <-
  .gg.plot(dfm[dfm$variable=="KL loss receptors",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "Optimization step", title = "", y = "KL loss receptors")

plotlist = list()
plotlist[[1]] = p1 
plotlist[[2]] = p2 
plotlist[[3]] = p3 

stplt <- grid.arrange(grobs=plotlist,ncol=3,
heights = c(1/3, 1/3,1/3))

f = paste(args_home,args$output,args$interaction_topic$out,args$interaction_topic$model_info,"loss_plot.pdf",sep="")
ggsave(f,stplt)
}

loss_plot(args)