#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=deicode
#SBATCH --dependency=afterok:41736717
#SBATCH --output=CanadaJayMicrobiome/outputs/%x-%j.out

#script starts here
#----------------------------------
# dependent on 6taxonomy.sh

# filtering out blanks (for H1 & H2)
qiime feature-table filter-samples \
    --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
    --m-metadata-file input/jay-met.tsv \
    --p-exclude-ids 'TRUE' \
    --p-where "[JayID]='BLANK'" \
    --o-filtered-table filtered-table-no-blanks.qza

# Hypothesis 1 - rarefaction
# rarefying to feed into core definition in R (needed nowhere else)
# 344 retains all samples
qiime feature-table rarefy \
    --i-table filtered-table-no-blanks.qza \
    --p-sampling-depth 344 \
    --p-no-with-replacement \
    --o-rarefied-table rarefied-table
# ran gj-core.R here, then produced a list of core and rare features
# filter for core variants
qiime feature-table filter-features \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/coreFeatures.tsv \
  --o-filtered-table filtered-table-core.qza
# filter for rare variants
qiime feature-table filter-features \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/coreFeatures.tsv \
  --p-exclude-ids TRUE \
  --o-filtered-table filtered-table-rare.qza
# Hypothesis 1 & 2 - full dataset
# compute aitchison distance matrix 
qiime deicode rpca \
    --i-table filtered-table-no-blanks.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination.qza \
    --o-distance-matrix aitchison-distance.qza

qiime emperor biplot \
    --i-biplot aitchison-ordination.qza \
    --m-sample-metadata-file input/jay-met.tsv \
    --m-feature-metadata-file taxonomy/SILVA-taxonomy.qza \
    --o-visualization aitchison-biplot.qzv \
    --p-number-of-features 8
# repeat above for core variants
qiime deicode rpca \
    --i-table filtered-table-core.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination-core.qza \
    --o-distance-matrix aitchison-distance-core.qza

qiime emperor biplot \
    --i-biplot aitchison-ordination-core.qza \
    --m-sample-metadata-file input/jay-met.tsv \
    --m-feature-metadata-file taxonomy/SILVA-taxonomy.qza \
    --o-visualization aitchison-biplot-core.qzv \
    --p-number-of-features 8
# repeat above for rare variants
qiime deicode rpca \
    --i-table filtered-table-rare.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination-rare.qza \
    --o-distance-matrix aitchison-distance-rare.qza

qiime emperor biplot \
    --i-biplot aitchison-ordination-rare.qza \
    --m-sample-metadata-file input/jay-met.tsv \
    --m-feature-metadata-file taxonomy/SILVA-taxonomy.qza \
    --o-visualization aitchison-biplot-rare.qzv \
    --p-number-of-features 8
# Hypothesis 3 - spring 2020 samples (nest groups)
# using this i-table since we're filtering, so blanks will be dropped anyways
qiime feature-table filter-samples \
  --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason]='Spring' AND [CollectionYear]='2020'" \
  --o-filtered-table H3-filtered-table.qza

qiime deicode rpca \
    --i-table H3-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H3-aitchison-ordination.qza \
    --o-distance-matrix H3-aitchison-distance.qza
# creates a sample tree
qiime diversity beta-rarefaction \
    --i-table H3-filtered-table.qza \
    --i-phylogeny trees/rooted-tree.qza \
    --p-metric 'aitchison' \
    --p-clustering-method 'nj' \
    --m-metadata-file input/jay-met.tsv \
    --p-sampling-depth 500 \
    --p-iterations 100 \
    --p-correlation-method 'spearman' \
    --p-color-scheme 'PuOr_r' \
    --o-visualization H3-aitchison-beta-rarefaction

# Hypothesis 4 - only samples with origin data
# H4-samples.tsv is a list of sampleid's to keep
qiime feature-table filter-samples \
  --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
  --m-metadata-file CanadaJayMicrobiome/data/H4-samples.tsv \
  --o-filtered-table H4-filtered-table.qza

qiime deicode rpca \
    --i-table H4-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H4-aitchison-ordination.qza \
    --o-distance-matrix H4-aitchison-distance.qza