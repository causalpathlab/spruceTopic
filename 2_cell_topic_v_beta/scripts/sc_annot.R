install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")


library(celldex)
library(SingleR)
ref <- HumanPrimaryCellAtlasData()



experiment_f = '../output/scanpy/gse_normal_scanpy_processed.csv.gz'
df = read.table(experiment_f,sep=',')

pred.hesc <- SingleR(test = df, ref = hpca.se, assay.type.test=1,
    labels = hpca.se$label.main)


require(remotes)
install_version("dplyr", version = "0.5.0", repos = "http://cran.us.r-project.org")