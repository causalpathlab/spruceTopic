library(ggplot2)
library(gridExtra)
library(reshape)
library(pheatmap)
library(RColorBrewer)
library(Polychrome)
library(randomcoloR)
source("Util.R")


ccv_struct_plot <- function(df,f) {

n <- 50
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_vector = c("#F0A3FF", "#0075DC","#808080" ,"#4C005C","#2BCE48","#FFCC99","#993F00","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")

col_vector <- distinctColorPalette(n)

# col_vector <- c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
#           "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
#           "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
#           "springgreen2", "yellowgreen", "palegreen4",
#           "wheat2", "tan", "tan2", "tan3", "brown")

p <-
ggplot(df, aes(x=cluster_celltype, y=ncount,fill=as.factor(state))) +
  geom_bar(position="stack",stat="identity",size=0) +
  scale_fill_manual("Latent topic",values=col_vector)+
  facet_grid(~ cluster_celltype, scales = "free", switch = "x", space = "free")+
  labs(x = "Neighbour cells of cancer cells", y = "Interaction topic distribution")+
  theme(
    legend.position = "top",
    legend.justification = "left", 
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(10,10,10,10),
    text = element_text(size=50),
    panel.spacing.x = unit(0.005, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank())+
    # panel.background = element_rect(fill='transparent'),
    # plot.background = element_rect(fill='transparent', color=NA))+
  guides(fill = guide_legend(nrow = 1))

ggsave(f,p,width = 60, height = 20,limitsize=F)
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