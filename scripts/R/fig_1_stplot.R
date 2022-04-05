library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
source("Util.R")
config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/scmetm.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"


struct_plot <- function(args) {

hh_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"hh_cell_topic_sample.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)

df_h_m = melt(df_h,id=c("cell","cluster"))
df_h_m$cluster <- factor(df_h_m$cluster)

colnames(df_h_m) = c("cell", "cluster", "Topic", "hvalue")

df_h_m$Topic <- gsub("X","",df_h_m$Topic)

p <-
ggplot(df_h_m, aes(x=cell, y=hvalue,fill=Topic)) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_brewer("Latent topic", palette = "Paired")+
  facet_grid(~ cluster, scales = "free", switch = "x", space = "free")+
  labs(x = "Cells", y = "Topic Proportion")+
  theme(
    legend.position = "top",
    legend.justification = "left", 
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(10,10,10,10),
    text = element_text(size=50),
    panel.spacing.x = unit(0.005, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA))+
  guides(fill = guide_legend(nrow = 1))
 

ggsave("hh_struct_plot.pdf",p,width = 45, height = 20)
}

struct_plot(args)