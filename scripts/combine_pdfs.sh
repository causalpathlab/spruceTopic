#!/bin/bash

# path="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/output/tcell/model/lr_net/*.log"
# flist=$(ls $path | sed 's/\.log//g')

# for i in ${flist[0]}
# do
#     echo $i
#     pdfunite $i"loss_plot.pdf" $i"top5_genes_hmap.pdf"  $i"top5_lr_pair_topic_hmap.pdf" $i"beta1_bias.pdf" $i"beta2_bias.pdf"  $i"CD4_network.pdf"  $i"CD8_network.pdf"  $i"COMBINE.pdf" 
                
# done

path="/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/output/pbmc/model/lr_net/*.log"
flist=$(ls $path | sed 's/\.log//g')

for i in ${flist[0]}
do
    echo $i
    pdfunite $i"loss_plot.pdf" $i"top5_genes_hmap.pdf"  $i"top5_lr_pair_topic_hmap.pdf" $i"beta1_bias.pdf" $i"beta2_bias.pdf"  $i"network.pdf" $i"COMBINE.pdf" 
                
done
