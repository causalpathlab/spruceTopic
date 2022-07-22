library(ggplot2)
library(gridExtra)
library(reshape)
library(pheatmap)
library(RColorBrewer)
library(Polychrome)
library(randomcoloR)
source("Util.R")

ccv_struct_plot <- function(df,f) {

n <- length(unique(df$cluster_celltype))
print(n)
col_vector <- as.vector(alphabet.colors(26))

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