qiime feature-classifier classify-sklearn \
  --i-classifier references/silva-132-99-nb-classifier.qza \
  --p-n-jobs 16 \
  --i-reads rep-seqs-no-singletons-cr-99.qza \
  --o-classification taxonomy/SILVA-taxonomy-cr-99.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --m-metadata-file input/jay-met.tsv \
  --o-visualization taxonomy/SILVA-cr99-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-level 2 \
  --o-collapsed-table taxonomy/SILVA-table-l2-cr-99.qza
# Class
qiime taxa collapse \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-level 3 \
  --o-collapsed-table taxonomy/SILVA-table-l3-cr-99.qza
# Order
qiime taxa collapse \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-level 4 \
  --o-collapsed-table taxonomy/SILVA-table-l4-cr-99.qza
# Family
qiime taxa collapse \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-level 5 \
  --o-collapsed-table taxonomy/SILVA-table-l5-cr-99.qza
# Genus
qiime taxa collapse \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-level 6 \
  --o-collapsed-table taxonomy/SILVA-table-l6-cr-99.qza
# Species
qiime taxa collapse \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-level 7 \
  --o-collapsed-table taxonomy/SILVA-table-l7-cr-99.qza

# drop mitochondrial and chloroplast sequences
# these aren't necessarily informative of bacterial diversity
qiime taxa filter-table \
  --i-table filtered-table-no-singletons-cr-99.qza \
  --i-taxonomy taxonomy/SILVA-taxonomy-cr-99.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table filtered-table-no-singletons-mitochondria-chloroplast-cr-99.qza