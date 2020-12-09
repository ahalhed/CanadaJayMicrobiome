#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=deicode
#SBATCH --output=CanadaJayMicrobiome/output/%x-%j.out

#script starts here
#----------------------------------
# dependent on 6taxonomy.sh

# compositional data analysis
# compute aitchison distance matrix 
# full dataset (H1, H2)
qiime deicode rpca \
    --i-table filtered-table-no-singletons-mitochondria-chloroplast.qza \
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


# will test these two sections of code once have complete sequence set
# Hypothesis 3 - spring 2020 samples (nest groups)
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

# Hypothesis 4 - only samples with origin data
# H4-samples.tsv is a list of sampleid's to keep (will complete list once have all metadata)
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