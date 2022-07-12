install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")


library(celldex)
library(SingleR)
library(SingleCellExperiment)



experiment_f = 'gse_normal_scanpy_processed.csv.gz'
df = read.table(experiment_f,sep=',',header=TRUE)
rownames(df) = df$index
df$index=NULL
df = t(df)

df <- as.matrix(df) 
sce <- SingleCellExperiment(assays = list(counts = df))

# ref <- ImmGenData()
ref <- load("CHETAH_TME_reference.Rdata") 
assayNames(reference) <- 'logcounts'
assayNames(sce) <- 'logcounts'
pred_chetah <- SingleR(test = sce, ref = reference, assay.type.test=1,labels = reference$celltypes)
pred_chetah = as.data.frame(pred_chetah)
pred_chetah = pred_chetah[,c("labels","pruned.labels")]
write.csv(pred_chetah,'gse_normal_ann_chetah.csv')

ref <- HumanPrimaryCellAtlasData()
pred_hpc <- SingleR(test = sce, ref = ref, assay.type.test=1,labels = ref$label.main)


require(remotes)
install_version("rjson", version = "0.5.0", repos = "http://cran.us.r-project.org")