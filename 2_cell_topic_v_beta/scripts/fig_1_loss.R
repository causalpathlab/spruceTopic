library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(RColorBrewer)
setwd(box::file())
source("Util.R")

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/0_augmented_cell_topic_v_beta/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)

loss_plot <- function(args) {
loss_file = paste(args$home,args$experiment,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cell_topic_loss2.txt.gz",sep="")
print(loss_file)
df = read.table(loss_file, sep = ",", header=TRUE)
print(head(df))
colnames(df) = c("Log-likelihood","KL loss","KLB loss")
df$epoch <- 1:nrow(df)
dfm = melt(df,id="epoch")
 
p1 <-
  .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
    geom_point(stroke = 0, color="gray", size=1) +
    geom_smooth(color="red", se=FALSE, size=1) +
    labs(x = "epoch", title = "", y = "Log-likelihood")
  

p2 <-
  .gg.plot(dfm[dfm$variable=="KL loss",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "epoch", title = "", y = "KL loss")

p3 <-
  .gg.plot(dfm[dfm$variable=="KLB loss",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
  geom_smooth(color="red", se=FALSE, size=1) +
  labs(x = "epoch", title = "", y = "KLB loss")


plotlist = list()
plotlist[[1]] = p1 
plotlist[[2]] = p2 
plotlist[[3]] = p3

stplt <- grid.arrange(grobs=plotlist,ncol=3,
heights = c(1/3, 1/3, 1/3))

f = paste(args$home,args$experiment,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_cell_topic_loss_plot.pdf",sep="")
ggsave(f,stplt)
}

loss_plot(args)