library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(RColorBrewer)
library(Polychrome)
library(data.table)
library(ggh4x)
source("Util.R")
library(viridis)


subtype_heatmap <- function(){
file = paste(it_id,"5_subtype_heatmap_r.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

dfm = melt(df,id=c('Cancer', 'subtype', 'cell_topic' ))

dfm$value = ifelse(dfm$value > 3.0,3.0,dfm$value)

p <-
ggplot(dfm, aes(x=Cancer, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "black", high = "yellow")+
  facet_grid(~ cell_topic + subtype, scales = "free", switch = "x", space = "free")+
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      text = element_text(size=8,colour="black"),
        panel.margin=unit(.01, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))
    #   axis.text.y=element_blank(),
    #   axis.ticks.y=element_blank())

ggsave(paste(it_id,"5_subtype_heatmap_r.pdf",sep=""),p,width = 8, height = 6,dpi=1200)

file = paste(it_id,"5_subtype_heatmap_l.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

dfm = melt(df,id=c('Cancer', 'subtype', 'cell_topic' ))

dfm$value = ifelse(dfm$value > 3.0,3.0,dfm$value)

p <-
ggplot(dfm, aes(x=Cancer, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "black", high = "yellow")+
  facet_grid(~ cell_topic + subtype, scales = "free", switch = "x", space = "free")+
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      text = element_text(size=8,colour="black"),
        panel.margin=unit(.01, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))
    #   axis.text.y=element_blank(),
    #   axis.ticks.y=element_blank())

ggsave(paste(it_id,"5_subtype_heatmap_l.pdf",sep=""),p,width = 8, height = 6,dpi=1200)

}

subtype_heatmap()