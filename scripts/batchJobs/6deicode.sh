#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=deicode
#SBATCH --dependency=afterok:44114231
#SBATCH --output=CanadaJayMicrobiome/outputs/%x-%j.out

#script starts here
#----------------------------------
# dependent on 5taxonomy.sh

# filtering out blanks
qiime feature-table filter-samples \
    --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
    --m-metadata-file input/jay-met.tsv \
    --p-exclude-ids 'TRUE' \
    --p-where "[JayID]='BLANK'" \
    --o-filtered-table filtered-table-no-blanks.qza
# rarefying to feed into core definition in R (needed nowhere else)
# 344 retains all samples
qiime feature-table rarefy \
    --i-table filtered-table-no-blanks.qza \
    --p-sampling-depth 344 \
    --p-no-with-replacement \
    --o-rarefied-table rarefied-table
# ran gj-core.R here to produce a list of core and rare features

# Hypothesis 1
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
    --o-visualization H1biplot.qzv \
    --p-number-of-features 8

# Hypothesis 2
# Prediction 2B
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionYear] IN ('2017', '2018')" \
  --o-filtered-table H2filtered-table.qza

qiime deicode rpca \
    --i-table H2filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H2aitchison-ordination.qza \
    --o-distance-matrix H2aitchison-distance.qza


# Hypothesis 3 - spring 2020 samples (nest groups)
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason]='Spring' AND [CollectionYear]='2020'" \
  --o-filtered-table H3filtered-table.qza

qiime deicode rpca \
    --i-table H3filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H3aitchison-ordination.qza \
    --o-distance-matrix H3aitchison-distance.qza

# Hypothesis 4 - only samples with origin data
# H4-samples.tsv is a list of sampleid's to keep
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/H4-samples.tsv \
  --o-filtered-table H4filtered-table.qza

qiime deicode rpca \
    --i-table H4filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot H4aitchison-ordination.qza \
    --o-distance-matrix H4aitchison-distance.qza
