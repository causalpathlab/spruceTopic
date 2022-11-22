#### bipartite 
#https://rdrr.io/cran/bipartite/man/plotweb.html
library(bipartite)
library(ggplot2)
library(gridExtra)
library(reshape)

file = paste(it_id,"8_lrnetwork_stringdb.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)

for ( i in c(2,4,7,10,18,22,24)){
print(i)
dfm = df[df$topic==i,]
if (dim(dfm)[1]>100){dfm = dfm[c(1:100),]}
print(max(dfm$score))
print(min(dfm$score))
matrix = table(dfm[,c(1,2)])
edge_colors = colorRampPalette(c("grey80","tan2"))(length(dfm$score))
f=paste(it_id,"8_lrnetwork_string_db",i,".pdf",sep="")
pdf(f)
p <- plotweb(matrix, method="normal", 
labsize=1,text.rot = 90,
col.high='salmon',col.low='cyan',
col.interaction = edge_colors,
bor.col.high='black',bor.col.low='black',bor.col.interaction='NA',
y.width.low=0.05,y.width.high=0.05,
arrow='no',ybig=0.1)
dev.off()
}
