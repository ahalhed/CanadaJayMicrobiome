#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=FastTreePipeline
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# create the phylogenetic tree
# fasttree pipeline
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree ./trees/unrooted-tree.qza \
  --o-rooted-tree ./trees/rooted-tree.qza
