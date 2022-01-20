#!/bin/bash

for i in BCL BC ESCA MM PACA RC THCA UCEC
do
    # wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156728/suppl/GSE156728_${i}_10X.CD4.counts.txt.gz"
    wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156728/suppl/GSE156728_${i}_10X.CD8.counts.txt.gz"
done
