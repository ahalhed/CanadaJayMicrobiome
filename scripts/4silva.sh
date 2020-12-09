#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --ntasks-per-node=16
#SBATCH --time=0-00:20:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=jays-silva-march13
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------


# Classifying taxonomies
qiime feature-classifier classify-sklearn \
  --i-classifier ./references/silva-132-99-nb-classifier.qza \
  --p-n-jobs 16 \
  --i-reads rep-seqs-dada2.qza \
  --o-classification ./taxonomy/SILVA-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --m-metadata-file ./input/jay-FINAL.tsv \
  --o-visualization ./taxonomy/SILVA-dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./taxonomy/SILVA-table-l2.qza
# Class
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./taxonomy/SILVA-table-l3.qza
# Order
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./taxonomy/SILVA-table-l4.qza
# Family
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./taxonomy/SILVA-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./taxonomy/SILVA-table-l6.qza
# Species
qiime taxa collapse \
  --i-table OTU-table-dada2.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./taxonomy/SILVA-table-l7.qza
