# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2020.11

# script starts here
# ________________________________________

# working in a SHARCNET folder for the jays
cd /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome
# request interactive session for testing
# salloc --mem-per-cpu=4G --account=def-cottenie --time=0-01:00:00
# data transfer to graham (not github - for large files)
# rsync -avz --no-g --no-p <source> <destination>
# load miniconda
module load nixpkgs/16.09 miniconda3
# Activate QIIME2
conda activate qiime2-2020.11

# 6 - only breeders without supplementation and with territory quality information
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

# Close QIIME2
conda deactivate