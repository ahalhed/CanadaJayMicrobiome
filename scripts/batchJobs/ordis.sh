#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --ntasks-per-node=16
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=ordis
#SBATCH --output=CanadaJayMicrobiome/outputs/%x-%j.out

#script starts here
#----------------------------------

# dada2 tax (for ordination comparison)
qiime feature-classifier classify-sklearn \
  --i-classifier references/silva-138-99-nb-classifier.qza \
  --p-n-jobs 16 \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy/SILVA-taxonomy-dada2.qza
# filtering dada2 table (for ordination comparison)
qiime feature-table filter-samples \
    --i-table table-dada2.qza \
    --m-metadata-file input/jay-met.tsv \
    --p-exclude-ids 'TRUE' \
    --p-where "[JayID]='BLANK'" \
    --i-taxonomy taxonomy/SILVA-taxonomy-dada2.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table filtered-table-dada2.qza

# full ordinations (for comparison purposes)
qiime deicode rpca \
    --i-table filtered-table-no-blanks.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination-cr.qza \
    --o-distance-matrix aitchison-distance-cr.qza

qiime deicode rpca \
    --i-table filtered-table-dada2.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 2 \
    --o-biplot aitchison-ordination-dn.qza \
    --o-distance-matrix aitchison-distance-dn.qza