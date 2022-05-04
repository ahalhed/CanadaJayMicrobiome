# Author: Alicia Halhed
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2021.2 (this file only)

# script starts here
# ________________________________________

# local for 1B
# Activate QIIME2
conda activate qiime2-2021.2

# Prediction 1B - There will be seassonally associated functions due to the need for food preservation in the fall.
# this run locally in QIIME2-2021.2 (issues with picrust 2installation on cluster in 2020.11)
qiime picrust2 full-pipeline \
   --i-table filtered-table-no-blanks.qza \
   --i-seq rep-seqs-cr-99.qza \
   --output-dir 1B-picrust2_output \
   --p-placement-tool sepp \
   --p-threads 8 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose
# 2 of 2301 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.
# make visualization (pathway abundance like an OTU table)
qiime feature-table summarize \
   --i-table 1B-picrust2_output/pathway_abundance.qza \
   --o-visualization 1B-picrust2_output/pathway_abundance.qzv
# differential abundance testing - using aldex2 in r

conda deactivate