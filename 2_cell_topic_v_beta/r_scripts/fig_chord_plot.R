### chord diagram  ###
plot_chordDiagram <- function(){

file = paste(it_id,"8_chorddata.csv.gz",sep="")
df = read.table(file, sep = ",", header=TRUE)
df = df[order(df$topic),]
df$topic = as.factor(df$topic)

library(reshape)

dfm = melt(df,id=c('topic'))
dfm = dfm[,c('value','variable')]
colnames(dfm) = c('gene','genetype')

dfl = df[,c(1,2)]
library("circlize")
f = paste(it_id,"8_chorddata_ligand.pdf",sep="")
pdf(f)
# df[,c(1,3,4)]
genes<- as.character(dfl[[2]])
othercol = structure(rep("grey", length(genes)), names = genes)
grid_col = c("2" = "#b3446c", "4" = "#dcd300", "7" = "#882d17",
"10"="#8db600","18"="#654522","22"="#e25822","24"="#2b3d26", othercol )




# chordDiagram(df,,grid.col=grid_col,col = grid_col[as.character(df[[1]])],annotationTrack=NULL)


# Plot chord diagram
chordDiagram(dfl,
             grid.col=grid_col,
             col = grid_col[as.character(dfl[[1]])],
             annotationTrack = c("grid"), 
             preAllocateTracks = 1
             )


highlight.sector(dfl$ligand,
                 track.index = 1, col = "white",
                 text = "Ligand", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
highlight.sector(dfl$topic,
                 track.index = 1, col = "white",
                 text = "Interaction topics", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.4)
#   circos.axis(h = "top", labels.cex = 0, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  circos.axis(h = "top",labels=NULL)
}, bg.border = NA)
dev.off()
circos.clear()


dfr = df[,c(1,3)]
library("circlize")
f = paste(it_id,"8_chorddata_receptor.pdf",sep="")
pdf(f)
# df[,c(1,3,4)]
genes<- as.character(dfr[[2]])
othercol = structure(rep("grey", length(genes)), names = genes)
grid_col = c("2" = "#b3446c", "4" = "#dcd300", "7" = "#882d17",
"10"="#8db600","18"="#654522","22"="#e25822","24"="#2b3d26", othercol )




# chordDiagram(df,,grid.col=grid_col,col = grid_col[as.character(df[[1]])],annotationTrack=NULL)


# Plot chord diagram
chordDiagram(dfr,
             grid.col=grid_col,
             col = grid_col[as.character(dfr[[1]])],
             annotationTrack = c("grid"), 
             preAllocateTracks = 1, 

             )
highlight.sector(df$receptor,
                 track.index = 1, col = "white",
                 text = "Receptor", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
highlight.sector(dfr$topic,
                 track.index = 1, col = "white",
                 text = "Interaction topics", cex = 1, text.col = "black", 
                 niceFacing = TRUE, font=2)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.4)
#   circos.axis(h = "top", labels.cex = 0, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  circos.axis(labels=NULL)
}, bg.border = NA)
dev.off()
circos.clear()

}

