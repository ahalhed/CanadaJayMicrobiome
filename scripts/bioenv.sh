#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-06:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=jays-bioenv
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------
# Species: Canada (Grey) Jay
# Sample: Oral Swabs
# 16S rRNA
# working with qiime2-2020.2

# script starts here
# ________________________________________

# working in a SHARCNET folder for the jays
cd /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome
# request interactive session for testing
# salloc --mem-per-cpu=4G --account=def-cottenie --time=0-01:00:00
# data transfer to graham (not github - for large files)
# rsync -avz --no-g --no-p <source> <destination>
#submit job with conda environment activated

qiime diversity bioenv \
  --i-distance-matrix aitchison-distance.qza \
  --m-metadata-file input/jay-met.tsv \
  --o-visualization bioenv-aitchison.qzv
