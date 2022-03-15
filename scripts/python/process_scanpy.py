import numpy as np
import pandas as pd
import scanpy as sc


sc.settings.verbosity = 3            
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = '../output/pbmc3k.h5ad'

adata = sc.read_10x_mtx(
    '../input/pbmc3k/filtered_gene_bc_matrices/hg19/', 
    var_names='gene_symbols',                
    cache=True)   

# df = adata.to_df()
# df = df.reset_index()
# df.to_csv("../output/pbmc3k_scanpy_raw_counts.txt.gz",index=False,compression="gzip")
### run through scanpy filtering as suggested by authors
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')

##this section needed to assign cell type
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
clust = [int(x) for x in adata.obs.leiden.values]
cell_type = {
    0:'CD4 T', 1:'CD14 Monocytes', 2:'CD4 Mem',
    3:'B', 4:'CD8 T',
    5:'NK', 6:'FCGR3A Monocytes',
    7:'Dendritic', 8:'Megakaryocytes'}
cell_type_list = [cell_type[x] for x in clust]

df = adata.to_df()

df.to_csv("../output/pbmc3k_scanpy_filtered_counts.txt.gz",index=False,compression="gzip")
