library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(pheatmap)
library(RColorBrewer)
setwd(box::file())

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/sishirsubedi/projects/experiments/spruce_topic/0_augmented_cell_topic_v_beta/"
config = paste(args_home,"config/",args[1],".yaml",sep="") 
args = read_yaml(config)

struct_plot <- function(args) {

hh_file = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_netm_h_topic_sample.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)

df_h_m = melt(df_h,id=c("cell","Topic"))
print(head(df_h_m))
df_h_m$Topic <- factor(df_h_m$Topic)

colnames(df_h_m) = c("cell", "cluster", "Topic", "hvalue")

df_h_m$Topic <- gsub("X","",df_h_m$Topic)
# df_h_m$Topic <- factor(df_h_m$Topic, levels = c("0", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14"))

n <- 25
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_vector = c("#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")

# col_vector <- distinctColorPalette(25)

p <-
ggplot(df_h_m, aes(x=cell, y=hvalue,fill=Topic)) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_manual("Latent topic",values=col_vector)+
  facet_grid(~ cluster, scales = "free", switch = "x", space = "free")+
  labs(x = "Cells", y = "Topic proportion")+
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

f = paste(args_home,args$output,args$cell_topic$out,args$cell_topic$model_info,args$cell_topic$model_id,"_netm_st_plot.pdf",sep="")
ggsave(f,p,width = 60, height = 20,limitsize=F)
}

struct_plot(args)