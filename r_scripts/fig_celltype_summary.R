library(ggplot2)
library(gridExtra)
library(reshape)
library(pheatmap)
library(RColorBrewer)
library(Polychrome)
library(randomcoloR)
library(ggh4x)
source("Util.R")

ccv_struct_plot <- function(df,f,tag) {

n <- length(unique(df$cluster_celltype))
print(n)
col_vector <- as.vector(alphabet.colors(26))

col_vector <- c("orange", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941","yellowgreen" , "#7A4900","#FFDBE5",  "#0000A6")

p <-
ggplot(df, aes(x=interact_topic, y=ncount,fill=as.factor(cluster_celltype))) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_manual("Cell type ",values=col_vector)+
  facet_grid(~ interact_topic, scales = "free", switch = "x", space = "free")+
  labs(x = "Interaction topic of cancer cells", y = "Celltype distribution")+
  theme(
    legend.position = "right",
    legend.justification = "left", 
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(10,10,10,10),
    text = element_text(size=75),
    panel.spacing.x = unit(0.005, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank())+
    # panel.background = element_rect(fill='transparent'),
    # plot.background = element_rect(fill='transparent', color=NA))+
  guides(fill = guide_legend(nrow = n))

ggsave(f,p,width = 30, height = 20,limitsize=F)
}

ccv_struct_plot_v2 <- function(df,f,tag) {

n <- length(unique(df$cluster_celltype))
print(n)
col_vector1 <- as.vector(kelly.colors(22))[16:22]

col_vector <- c("orange", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941","yellowgreen" , "#7A4900","#FFDBE5",  "#0000A6")

  df$interact_topic = as.factor(df$interact_topic)
  # cats = unique(df$interact_topic)
  # i=1
  # plotlist = list()

  # for (ct in cats){
  #   print(ct)

p <-
ggplot(df, aes(x=interact_topic, y=ncount,fill=as.factor(cluster_celltype))) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_manual("Cell type ",values=col_vector)+
  facet_nested(~ interact_topic+celltype,scales = "free",
      strip = strip_nested(
      background_x = elem_list_rect(fill=c(col_vector1)))
  ,labeller = labeller(celltype = function(x) {rep("", length(x))})
  )+
          labs(x ="", y = "")+
          # labs(x = "Neighbour topic / Cell type", y = "Cells")+
          theme(
            legend.position = "right",
            text = element_text(size=25),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_blank(),
            # axis.ticks.y = element_blank(),
            panel.spacing.x = unit(0.005, "lines"),
            panel.grid = element_blank(),
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA)
                        )+
          guides(fill = guide_legend(nrow = n))
    
  #   plotlist[[i]] = p 
  #   i = i + 1

  #   }

  # stplt <- grid.arrange(grobs=plotlist,ncol=3)
  ggsave(f,p,width =15, height = 5)
}
# ccv_struct_plot_v2 <- function(df,f,tag) {

# n <- length(unique(df$cluster_celltype))
# print(n)
# col_vector <- as.vector(alphabet.colors(26))

# col_vector <- c("orange", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941","yellowgreen" , "#7A4900","#FFDBE5",  "#0000A6")

#   df$interact_topic = as.factor(df$interact_topic)
#   # cats = unique(df$interact_topic)
#   # i=1
#   # plotlist = list()

#   # for (ct in cats){
#   #   print(ct)

# p <-
# ggplot(df, aes(x=interact_topic, y=ncount,fill=as.factor(cluster_celltype))) +
#   geom_bar(position="stack",stat="identity",size=0) +
#   scale_fill_manual("Cell type ",values=col_vector)+
#   facet_nested(~ interact_topic+celltype,scales = "free", space = "free")+
#           labs(x ="", y = "")+
#           # labs(x = "Neighbour topic / Cell type", y = "Cells")+
#           theme(
#             legend.position = "right",
#             text = element_text(size=25),
#             axis.text.x = element_blank(),
#             axis.ticks.x = element_blank(),
#             # axis.text.y = element_blank(),
#             # axis.ticks.y = element_blank(),
#             panel.spacing.x = unit(0.005, "lines"),
#             panel.grid = element_blank(),
#             panel.background = element_rect(fill='transparent'),
#             plot.background = element_rect(fill='transparent', color=NA),
#             strip.background = element_rect(colour = "black", size = 1)
#             )+
#           guides(fill = guide_legend(nrow = n))
    
#   #   plotlist[[i]] = p 
#   #   i = i + 1

#   #   }

#   # stplt <- grid.arrange(grobs=plotlist,ncol=3)
#   ggsave(f,p,width =20, height = 8)
# }

cancer_nbr_lr_plot <- function(df,f) {


df = as.data.frame(df)
rownames(df) = df$index
df$index = NULL

row_order = row.order(df)

df_t = df
df_t$topic = rownames(df)
df_t = melt(df_t)
colnames(df_t)=c('row','col','weight')
col_order = col.order(df_t,row_order)

df = df[,col_order]
df = df[row_order,]

df = df + 1
df = log10(df)

p1 <- pheatmap(df,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),fontsize_row=6,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=T,show_colnames=T)

ggsave(f,p1)
}