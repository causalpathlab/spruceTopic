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
plt.savefig("../scanpy_pbmc_pca.png")

df = adata.to_df()

df.to_csv("../output/pbmc3k_scanpy_filtered_counts.txt.gz",index=False,compression="gzip")