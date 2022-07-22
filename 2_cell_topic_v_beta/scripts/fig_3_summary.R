library(ggplot2, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.6/')
library(gridExtra)
library(reshape)
library(yaml)
library(RColorBrewer)
library(Polychrome)
library(data.table)
library(ggh4x)
source("Util.R")

summary_plot_all <- function(df,f,col) {

  col_vector = c("#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFFF00")
  interaction_states = sort(unique(df$state))
  is_colors <-setNames(col_vector[0:length(interaction_states)], interaction_states)

    plotlist = list()

    cats = unique(df$celltype)
    i=1
    for (ct in cats){
        p <-
        ggplot(df[df$celltype == ct, ], aes(x=topic, y=cell,fill=as.factor(state))) +
          geom_bar(position="stack",stat="identity",size=0) +
          scale_fill_manual("Interaction topic",values=is_colors)+
          facet_nested(~ celltype + topic, scales = "free", switch = "x", space = "free")+
          labs(x = "", y = "Cells")+
          # labs(x = "Neighbour topic / Cell type", y = "Cells")+
          theme(
            legend.position = "top",
            legend.justification = "left", 
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin=margin(10,10,10,10),
            text = element_text(size=12),
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

  stplt <- grid.arrange(grobs=plotlist,ncol=col)
  ggsave(f,stplt,width =10, height = 15)
}


summary_plot_cancer <- function(df,f,col) {

# n <- length(unique(df$interact_topic))
# print(n)
col_vector <- as.vector(kelly.colors(22))[16:22]
# print(col_vector)


    plotlist = list()

    cats = unique(df$celltype)
    i=1
    for (ct in cats){
        p <-
        ggplot(df[df$celltype == ct, ], aes(x=cell_topic, y=ncount,fill=as.factor(interact_topic))) +
          geom_bar(position="stack",stat="identity",size=0) +
          scale_fill_manual("Interaction topic",values=col_vector)+
          facet_nested(~ celltype + cell_topic, scales = "free", switch = "x", space = "free")+
          labs(x ="", y = "")+
          # labs(x = "Neighbour topic / Cell type", y = "Cells")+
          theme(
            legend.position = "none",
            text = element_text(size=12),
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

  stplt <- grid.arrange(grobs=plotlist,ncol=col)
  ggsave(f,stplt,width =20, height = 8)
}
