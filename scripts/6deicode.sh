#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=deicode
#SBATCH --output=output/%x-%j.out

#script starts here
#----------------------------------
# dependent on 5taxonomy.sh

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