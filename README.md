## SPRUCE: Single-cell Pairwise Relationships Untangled by Composite Embedding model

<div align="center">
    <img src="images/spruce.png" alt="Logo" width="500" height="500">
</div>

###
This is a project repository for our paper-
* Subedi, S. and Park, Y.P., <a href="https://www.cell.com/iscience/fulltext/S2589-0042(23)00102-5"> Single-cell Pairwise Relationships Untangled by Composite Embedding model, </a> iScience, 2023.

### Summary
In multi-cellular organisms, cell identity and functions are primed and refined through interactions with other surrounding cells. Here, we propose a scalable machine learning method, termed SPURCE, which is designed to systematically ascertain common cell-cell communication patterns embedded in single-cell RNA-seq data. We applied our approach to investigate tumour microenvironments consolidating multiple breast cancer data sets and found seven frequently-observed interaction signatures and underlying gene-gene interaction networks. Our results implicate that a part of tumour heterogeneity, especially within the same subtype, is better understood by differential interaction patterns rather than the static expression of known marker genes.

### Prerequisites

* python - numpy, pandas, scipy, sklearn, annoy, pytorch, igraph, seaborn
* R - celldex, SingleR, SingleCellExperiment, ggplot2, pheatmap, circlize, bipartite

### Dataset
 * <a href="https://pubmed.ncbi.nlm.nih.gov/34493872/">Breast cancer cells </a>
 * <a href="https://pubmed.ncbi.nlm.nih.gov/33763657/">Normal breast cells</a>
 * <a href="https://pubmed.ncbi.nlm.nih.gov/34914499/">Immune cells from breast Cancer </a>

### Installation

* Clone the repo
   ```sh
   git clone https://github.com/causalpathlab/spruceTopic.git
   ```
