library(fgsea)
library(data.table)
library(ggplot2)
library(pheatmap)
source("Util.R")


run_gsea <- function(df,cutoff){
msigdbr_df <- msigdbr::msigdbr(species = "human", category = "C2",subcategory = "KEGG")
pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

selected_topics = c(3,5,8,11,19,23,25)

fr = data.frame()
for ( i in selected_topics ){
print(paste('processing',i-1))
ranks = as.double(sort(df[i,],descending=FALSE))
names(ranks) = colnames(sort(df[i,],descending=FALSE))
fgseaRes = fgsea(pathways = pathways, stats =ranks ,nperm=10000)
# fgseaRes = fgseaRes[order(pval), ]
fgseaRes[,'topic'] = paste('topic',i-1) 
fgseaRes = fgseaRes[fgseaRes$padj<cutoff,]
fgseaRes = as.data.frame(fgseaRes)
fr <- rbind(fr,fgseaRes)
}
dfm=fr
dfm$padj = -log10(dfm$padj)
dfm = dcast(dfm,pathway~topic,value.var="padj")
dfm[is.na(dfm)] <- 0
# dfm$topic = as.factor(dfm$topic)

rownames(dfm) = dfm$pathway
dfm$pathway = NULL

row_order = row.order(dfm)
dfm = dfm[row_order,]

# p <-   pheatmap(dfm,color =colorRampPalette(c("white", "tan2"))(100),fontsize_row=6,fontsize_col=8,cluster_rows=FALSE,cluster_cols=FALSE,column_names_side = c("top"),angle_col = c("45"))
dfm$pathway = rownames(dfm)

dfm2 = melt(dfm)
p = ggplot(dfm, aes(x = variable, y = pathway)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
#   scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
#   labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
#   theme(legend.key=element_blank(), 
#   axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
#   axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
#   legend.text = element_text(size = 10, face ="bold", colour ="black"), 
#   legend.title = element_text(size = 12, face = "bold"), 
#   panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
#   legend.position = "right") +  
#   scale_fill_manual(values = colours, guide = FALSE) + 
#   scale_y_discrete(limits = rev(levels(pcm$variable))) 

ggsave(paste(it_id,"12_gsea.pdf",sep=""),p,width = 5, height = 6,dpi=600)

}
