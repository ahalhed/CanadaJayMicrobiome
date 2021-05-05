# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2021.2 (end of this file only)

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

# filtering
qiime feature-table filter-samples \
  --i-table filtered-table-no-blanks.qza \
  --m-metadata-file input/jay-met.tsv \
  --p-where "[CollectionSeason] IN ('Fall')" \
  --o-filtered-table 1B-filtered-table-fall.qza

qiime feature-table filter-features \
  --i-table 1B-filtered-table-fall.qza \
  --m-metadata-file CanadaJayMicrobiome/data/coreFeatures.tsv \
  --o-filtered-table 1B-filtered-table.qza

# local for 1B
# Activate QIIME2
conda activate qiime2-2021.2
# Prediction 1B - The most common microbiota will putatively function in food preservation.
# this run locally in QIIME2-2021.2 (issues with picrust 2installation on cluster in 2020.11)
qiime picrust2 full-pipeline \
   --i-table 1B-filtered-table.qza \
   --i-seq rep-seqs-cr-99.qza \
   --output-dir q2-picrust2_output \
   --p-placement-tool sepp \
   --p-threads 8 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose
# 2 of 2301 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.
# make visualization (pathway abundance like an OTU table)
qiime feature-table summarize \
   --i-table q2-picrust2_output/pathway_abundance.qza \
   --o-visualization q2-picrust2_output/pathway_abundance.qzv
# differential abundance testing - using aldex2 in r

conda deactivate