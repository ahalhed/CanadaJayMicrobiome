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
# Prediction 1A - only breeders without supplementation
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[BreedingStatus]='Breeder' AND [FoodSupplement]='N'" \
  --o-filtered-table P1A-filtered-table.qza

# Prediction 1B - only breeders without supplementation and with territory infor
qiime feature-table filter-samples \
  --i-table P1A-filtered-table.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[TerritoryQuality] IN ('H', 'M', 'L')" \
  --o-filtered-table P1B-filtered-table.qza
qiime deicode rpca \
    --i-table P1B-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P1B-aitchison-ordination.qza \
    --o-distance-matrix P1B-aitchison-distance.qza

# Prediction 1C - nest groups
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason]='Spring' AND [CollectionYear]='2020'" \
  --o-filtered-table P1C-filtered-table.qza

qiime deicode rpca \
    --i-table P1C-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P1C-aitchison-ordination.qza \
    --o-distance-matrix P1C-aitchison-distance.qza

# Hypothesis 2
# Prediction 2A - host associated factors (all samples)
qiime deicode rpca \
    --i-table filtered-table-no-blanks.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P2A-aitchison-ordination.qza \
    --o-distance-matrix P2A-aitchison-distance.qza


# Hypothesis 3
# Prediction 3A/B - food caching (winter and spring only)
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason] IN ('Winter', 'Spring')" \
  --o-filtered-table P3AB-filtered-table.qza

qiime deicode rpca \
    --i-table P3AB-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P3AB-aitchison-ordination.qza \
    --o-distance-matrix P3AB-aitchison-distance.qza

# Prediction 3C/D - food supplementation
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionYear] IN ('2017', '2018') AND [BreedingStatus]='Breeder'" \
  --o-filtered-table P3CD-filtered-table.qza

qiime deicode rpca \
    --i-table P3CD-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P3CD-aitchison-ordination.qza \
    --o-distance-matrix P3CD-aitchison-distance.qza

# Hypothesis 4 - parental care
# all predictions
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason]='Spring' AND [CollectionYear]='2020'" \
  --o-filtered-table P4ABCD-filtered-table.qza

qiime deicode rpca \
    --i-table P4ABCD-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P4ABCD-aitchison-ordination.qza \
    --o-distance-matrix P4ABCD-aitchison-distance.qza

# Hypothesis 5 - only samples with origin data
# H5-samples.tsv is a list of sampleid's to keep
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file CanadaJayMicrobiome/data/H5-samples.tsv \
  --o-filtered-table P5AB-filtered-table.qza

qiime deicode rpca \
    --i-table P5AB-filtered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot P5AB-aitchison-ordination.qza \
    --o-distance-matrix P5AB-aitchison-distance.qza
