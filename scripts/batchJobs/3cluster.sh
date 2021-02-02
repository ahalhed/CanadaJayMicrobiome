#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=vsearch
#SBATCH --dependency=afterok:44114169
#SBATCH --output=CanadaJayMicrobiome/outputs/%x-%j.out

#script starts here
#----------------------------------
# depends on 2denoise.sh

# Obtaining SILVA reference database 
# need this to do closed reference otu picking
#wget -O "references/silva-138-99-seqs.qza" "https://data.qiime2.org/2020.11/common/silva-138-99-seqs.qza"


qiime vsearch cluster-features-closed-reference \
  --i-table table-dada2.qza \
  --i-sequences rep-seqs-dada2.qza \
  --i-reference-sequences references/silva-138-99-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table table-cr-99.qza \
  --o-clustered-sequences rep-seqs-cr-99.qza \
  --o-unmatched-sequences unmatched-cr-99.qza
