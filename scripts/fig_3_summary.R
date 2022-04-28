library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(yaml)
library(RColorBrewer)
library(ggh4x)
library(data.table)
setwd(box::file())
source("Util.R")

summary_plot_tcell <- function(args,df,is_colors,ctype) {

plotlist = list()

dfm = df[df$label %like% ctype]


if (ctype=="CD4"){cats = c("Tn","Tm","Temra","Th","Treg")}
else {cats = c("Tm","Tem","Tk","Tex")}


i=1
for (ct in cats){
    p <-
    ggplot(dfm[dfm$label %like% ct ], aes(x=topic, y=cell,fill=as.factor(state))) +
      geom_bar(position="stack",stat="identity",size=0) +
      scale_fill_manual("Interaction topic",values=is_colors)+
      facet_nested(~ label + topic, scales = "free", switch = "x", space = "free")+
      labs(x = "", y = "Cells")+
      # labs(x = "Neighbour topic / Cell type", y = "Cells")+
      theme(
        legend.position = "top",
        legend.justification = "left", 
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(10,10,10,10),
        text = element_text(size=50),
        panel.spacing.x = unit(0.5, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))+
      guides(fill = guide_legend(nrow = 1))

    plotlist[[i]] = p 
    i = i + 1

}

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"summary_plot_",ctype,"v2.pdf",sep="")
stplt <- grid.arrange(grobs=plotlist,ncol=1)
ggsave(f,stplt,width = 60, height = 100,limitsize = FALSE)
dev.off()
}

summary_plot_pbmc <- function(args,df,is_colors) {

    p <-
    ggplot(df, aes(x=topic, y=cell,fill=as.factor(state))) +
      geom_bar(position="stack",stat="identity",size=0) +
      scale_fill_manual("Interaction topic",values=is_colors)+
      facet_nested(~ label + topic, scales = "free", switch = "x", space = "free")+
      labs(x = "", y = "Cells")+
      # labs(x = "Neighbour topic / Cell type", y = "Cells")+
      theme(
        legend.position = "top",
        legend.justification = "left", 
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(10,10,10,10),
        text = element_text(size=50),
        panel.spacing.x = unit(0.5, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))+
      guides(fill = guide_legend(nrow = 1))

f = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"summary_plot.pdf",sep="")
ggsave(f,p,width = 60, height = 20,limitsize = FALSE)
}

args = commandArgs(trailingOnly=TRUE)
args_home ="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/"
config = paste(args_home,"/config/",args[1],".yaml",sep="") 
args = read_yaml(config)

summary_file = paste(args_home,args$output,args$lr_model$out,args$lr_model$mfile,"model_summary.csv",sep="")


df = read.table(summary_file, sep = ",", header=TRUE)
df = data.table(df)
df = df[df$cell>10,]
df$topic <- gsub("hh","",df$topic)


col_vector = c("#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFFF00")

interaction_states = sort(unique(df$state))

is_colors <-setNames(col_vector[0:length(interaction_states)], interaction_states)

if (args[1]=="tcell"){

df = df[df$state != 17,]
df = df[df$state != 22,]

summary_plot_tcell(args,df,is_colors,"CD4") 
summary_plot_tcell(args,df,is_colors,"CD8") 
}

summary_plot_pbmc(args,df,is_colors)