#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=repSeqs
#SBATCH --output=CanadaJayMicrobiome/output/%x-%j.out

#script starts here
#----------------------------------

# Going to work with the DADA2 results
# will use this figure to make sampling depth decision
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file input/jay-met.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

# drop the singletons from the OTU table
qiime feature-table filter-features \
  --i-table table-dada2.qza \
  --p-min-samples 2 \
  --o-filtered-table filtered-table-no-singletons.qza
# filter representative sequences for this group
qiime feature-table filter-seqs \
  --i-data rep-seqs-dada2.qza \
  --i-table filtered-table-no-singletons.qza \
  --o-filtered-data rep-seqs-no-singletons.qza