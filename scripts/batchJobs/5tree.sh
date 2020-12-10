#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=FastTree
#SBATCH --dependency=afterok:41736706
#SBATCH --output=CanadaJayMicrobiome/outputs/%x-%j.out

#script starts here
#----------------------------------
# dependent on 4repSeqs.sh

# create the phylogenetic tree
# fasttree pipeline
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-no-singletons.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree trees/unrooted-tree.qza \
  --o-rooted-tree trees/rooted-tree.qza
