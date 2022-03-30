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


struct_plot <- function(args) {

hh_file = paste(args$home,args$output,args$model$mfile,"hh_cell_topic_sample.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)

colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")

df_h_m = melt(df_h,id=c("cell","cluster"))
df_h_m$cluster <- factor(df_h_m$cluster)#,      
      # levels = c("15", "11", "0", "4", "12","9","14","3","6","8","13","7","5","1","10","2"))

colnames(df_h_m) = c("cell", "cluster", "Topic", "hvalue")
p <-
  ggplot(df_h_m, aes(x=cell, y=hvalue,fill=Topic)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_manual(values=colours[1:17]) +
  facet_grid(~ cluster, scales = "free", switch = "x", space = "free")+
  labs(x = "Cells", title = "ETM", y = "Topic Proportion")+
   theme(
    text = element_text(size=50),
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) 
ggsave("hh_struct_plot.pdf",p,width = 45, height = 20)
}

umap_plot <- function(args) {
library(umap)
hh_file = paste(args$home,args$output,args$model["out"],args$model$mfile,"etm_hh_data.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)
hh_umap = umap(df_h[,2:17])
df = data.frame(x = hh_umap$layout[,1],y = hh_umap$layout[,2])

p = ggplot(df) + 
    geom_point(aes(x,y),size=0.000001,color='darkblue')+
    coord_fixed(xlim = c(-10, 10),ylim = c(-10, 10))+
    theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

ggsave("hh_umap_plot.png",p,width = 10, height = 10)
}
