#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=jays-gg-march13
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# green genes taxonomy
qiime feature-classifier classify-sklearn \
 --i-classifier ./references/gg-13-8-99-nb-classifier.qza \
 --p-n-jobs 16 \
 --i-reads rep-seqs-dada2.qza \
 --o-classification ./taxonomy/GG-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --m-metadata-file ./input/jay-FINAL.tsv \
  --o-visualization dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./taxonomy/GG-table-l2.qza
# Class
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./taxonomy/GG-table-l3.qza
# Order
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./taxonomy/GG-table-l4.qza
# Family
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./taxonomy/GG-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./taxonomy/GG-table-l6.qza
# Species
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./taxonomy/GG-table-l7.qza
