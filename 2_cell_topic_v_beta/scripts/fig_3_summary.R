library(ggplot2)
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


summary_plot_cancer <- function(df,f,tag,col) {

  col_vector <- as.vector(kelly.colors(22))[16:22]
    plotlist = list()

    cats = unique(df$celltype)
    i=1
    for (ct in cats){
      print(ct)
        p <-
        ggplot(df[df$celltype == ct, ], aes(x=cell_topic, y=ncount,fill=as.factor(interact_topic))) +
          geom_bar(position="stack",stat="identity",size=0) +
          scale_fill_manual("Interaction topic",values=col_vector)+
          facet_nested(~ celltype + cell_topic, scales = "free")+
          labs(x ="", y = "")+
          # labs(x = "Neighbour topic / Cell type", y = "Cells")+
          theme(
            legend.position = "none",
            text = element_text(size=20),
            panel.spacing.x = unit(0.1, "lines"),
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
  ggsave(f,stplt,width =15, height = 12)
}

summary_plot_cancer_v1 <- function(df,f,tag,col) {

  col_vector <- as.vector(kelly.colors(22))[16:22]

        p <- ggplot(df, aes(x=celltype, y=ncount,fill=as.factor(interact_topic))) +
          geom_bar(position="stack",stat="identity",size=0) +
          scale_fill_manual("Interaction topic",values=col_vector)+
          facet_nested(~ celltype,scales = "free", space = "free")+
          labs(x ="", y = "")+
          # labs(x = "Neighbour topic / Cell type", y = "Cells")+
          theme(
            legend.position = "none",
            text = element_text(size=30),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_blank(),
            # axis.ticks.y = element_blank(),
            panel.spacing.x = unit(0.1, "lines"),
            panel.grid = element_blank(),
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA))+
          guides(fill = guide_legend(nrow = 1))


      ggsave(f,p,width =10, height = 8)
}

summary_plot_cancer_v2 <- function(df,f,tag,col) {

# n <- length(unique(df$interact_topic))
# print(n)
col_vector <- as.vector(kelly.colors(22))[11:12]
print(col_vector)


    plotlist = list()

    cats = unique(df$cells) 
    i=1
    for (ct in cats){
      print(ct)
        p <-
        ggplot(df[df$cells == ct, ], aes(x=cell_topic, y=ncount,fill=as.factor(cell_status))) +
          geom_bar(position="stack",stat="identity",size=0) +
          scale_fill_manual("Cell type",values=col_vector)+
          facet_nested(~ cells + cell_topic, scales = "free", switch = "x", space = "free")+
          labs(x ="", y = "")+
          # labs(x = "Neighbour topic / Cell type", y = "Cells")+
          theme(
            legend.direction = "vertical",
            legend.position = "right",
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

summary_boxplot<- function(df_summary,f){
# grouped boxplot
col_vector <- as.vector(kelly.colors(22))[16:22]
df_summary$interact_topic = as.factor(df_summary$interact_topic)
p <- ggplot(df_summary, aes(x=interact_topic, y=ncount, fill=interact_topic)) + 
    geom_boxplot()+
    facet_wrap(~celltype,ncol=9)+
    labs(x ="", y = "")+
    scale_fill_manual("Interaction topic",values=col_vector)+    
      theme(
      legend.position = "none",
      text = element_text(size=25),
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill='grey95'),
      plot.background = element_rect(fill='transparent', color=NA))+
    guides(fill = guide_legend(nrow = 1))

ggsave(f,p,width =10, height = 4)
}