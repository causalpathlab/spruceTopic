library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(pheatmap)
library(dendextend)
library(RColorBrewer)

config = "/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/tcell.yaml" 
args = read_yaml(config)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"


struct_plot <- function(args) {

hh_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"hh_cell_topic_sample.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)

df_h_m = melt(df_h,id=c("cell","cluster"))
df_h_m$cluster <- factor(df_h_m$cluster)

colnames(df_h_m) = c("cell", "cluster", "Topic", "hvalue")

df_h_m$Topic <- gsub("X","",df_h_m$Topic)
df_h_m$Topic <- factor(df_h_m$Topic, levels = c("0", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14"))

# n <- 15
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = c("#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")[5:20]



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

ggsave("hh_struct_plot.pdf",p,width = 45, height = 20)
}

struct_plot_tclust <- function(args) {

hh_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"hh_cell_topic_sample_CD4CD8.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)

df_h_m = melt(df_h,id=c("cell","cluster"))
df_h_m$cluster <- factor(df_h_m$cluster)

colnames(df_h_m) = c("cell", "cluster", "Topic", "hvalue")

df_h_m$Topic <- gsub("k","",df_h_m$Topic)
df_h_m$Topic <- factor(df_h_m$Topic, levels = c("0", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14"))

# n <- 15
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = c("#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")[5:20]

p <-
ggplot(df_h_m, aes(x=cell, y=hvalue,fill=Topic)) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_manual("Latent topic",values=col_vector)+
  facet_grid(~ cluster, scales = "free", switch = "x", space = "free")+
  scale_y_discrete(expand = c(0, 0))+
  labs(x = "Cells", y = "Topic proportion")+
  theme(

    legend.position = "top",
    legend.justification = "left", 
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(10,10,10,10),

    text = element_text(size=50),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    strip.text.x = element_text(angle=90, color='black'),
    strip.background = element_rect(fill='transparent'),

    panel.background = element_rect(fill='transparent', color=NA),

    plot.background = element_rect(fill='transparent', color=NA))+
  guides(fill = guide_legend(nrow = 1))
 
f = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"cell_topic_struct_plot_cd4cd8.pdf",sep="")
ggsave(f,p,width = 45, height = 20)
}

struct_plot_pbmc <- function(args) {

hh_file = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"hh_cell_topic_sample.tsv",sep="")
df_h = read.table(hh_file, sep = "\t", header=TRUE)

df_h_m = melt(df_h,id=c("cell","cluster"))
df_h_m$cluster <- factor(df_h_m$cluster)

colnames(df_h_m) = c("cell", "cluster", "Topic", "hvalue")

df_h_m$Topic <- gsub("k","",df_h_m$Topic)
df_h_m$Topic <- factor(df_h_m$Topic, levels = c("0", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14"))

# n <- 15
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = c("#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")[5:20]

p <-
ggplot(df_h_m, aes(x=cell, y=hvalue,fill=Topic)) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_manual("Latent topic",values=col_vector)+
  facet_grid(~ cluster, scales = "free", switch = "x", space = "free")+
  scale_y_discrete(expand = c(0, 0))+
  labs(x = "Cells", y = "Topic proportion")+
  theme(

    legend.position = "top",
    legend.justification = "left", 
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(10,10,10,10),

    text = element_text(size=50),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    strip.text.x = element_text(angle=90, color='black'),
    strip.background = element_rect(fill='transparent'),

    panel.background = element_rect(fill='transparent', color=NA),

    plot.background = element_rect(fill='transparent', color=NA))+
  guides(fill = guide_legend(nrow = 1))
 
f = paste(args_home,args$output,args$nbr_model$out,args$nbr_model$mfile,"cell_topic_struct_plot_pbmc.pdf",sep="")
ggsave(f,p,width = 45, height = 20)
}
struct_plot_tclust(args)