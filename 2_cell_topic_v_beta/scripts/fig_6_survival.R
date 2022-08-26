###survival analysis ###
library("survminer")
require("survival")
library(ggplot2)


file = paste(it_id,"11_survival_analysis_it_2g.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

topics = c('topic2','topic4','topic7','topic10','topic18','topic22','topic24')
formulas = list()
fnames = list()
for ( t in topics)
{
  fnames[[t]] = paste(it_id,"11_survival_analysis_",t,".pdf",sep="")

  formulas[[t]] <- paste0("Surv(overall_time, donor_vital_status) ~ ", t) %>% as.formula()
}

fits <- surv_fit(formulas, data=df)
plots <- ggsurvplot(fits, risk.table = TRUE,pval=TRUE)

findx = 1
for (p in plots){
# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Remember to pass object `p`.
ggsave(fnames[[findx]], plot = p, 
      dpi=300, width = 10, height = 7, units = "in")

findx = findx + 1

}




file = paste(it_id,"11_survival_analysis_ct.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

topics = as.vector(colnames(df)[4:53])
formulas = list()
fnames = list()
for ( t in topics)
{
  fnames[[t]] = paste(it_id,"11_survival_analysis_",t,".pdf",sep="")

  formulas[[t]] <- paste0("Surv(overall_time, donor_vital_status) ~ ", t) %>% as.formula()
}

fits <- surv_fit(formulas, data=df)
plots <- ggsurvplot(fits, risk.table = TRUE,pval=TRUE)

findx = 1
for (p in plots){
# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Remember to pass object `p`.
ggsave(fnames[[findx]], plot = p, 
      dpi=300, width = 10, height = 7, units = "in")

findx = findx + 1

}


#####
library(Polychrome)
library(ggplot2)

file = paste(it_id,"11_survival_analysis_it_mod.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)
df = df[df$donor_vital_status==TRUE,]
df = df[order(df$icgc_donor_id),]

col_vector <- as.vector(kelly.colors(22))[16:22]
p <- ggplot(df, aes(x=icgc_donor_id, y=value,fill=as.factor(variable))) +
geom_bar(position="stack",stat="identity",size=0) +
scale_fill_manual("Interaction topic",values=col_vector)+
# facet_nested(~ celltype + topic, scales = "free", switch = "x", space = "free")+
labs(x = "Patient", y = "score")+
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

ggsave('test_dead_patient.png',p,width =10, height = 15)




file = paste(it_id,"11_survival_analysis_it_mod_1.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

df = df[,c(4,5,6,7,8,9,10)]

source('fig_1_hmap.R')

df = as.data.frame(df)
row_order = row.order(df)
df = df[row_order,]

df_beta_t = df
df_beta_t$topic = rownames(df)
df_beta_t = melt(df_beta_t)
colnames(df_beta_t)=c('row','col','weight')
df_beta_t = as.data.frame(df_beta_t)
col_order = col.order(df_beta_t,rownames(df)) 
df = df[,col_order]

p1 <- pheatmap(df,color = colorRampPalette(c("white", "blue"))(20),fontsize_row=6,fontsize_col=8,show_rownames=F,cluster_columns=F,cluster_rows=F)
ggsave('test4.png',p1,width =10, height = 15)


library(pheatmap)
p1 <- pheatmap(df[,c(4,5,6,7,8,9,10)],color = colorRampPalette(c("white", "blue"))(20),fontsize_row=6,fontsize_col=8,show_rownames=F)
ggsave('test2.png',p1,width =10, height = 15)

