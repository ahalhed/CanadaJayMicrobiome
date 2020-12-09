#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=ImportDemux
#SBATCH --output=./CanadaJayMicrobiome/output/%x-%j.out

#script starts here
#----------------------------------

# import the sequence data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path rawFiles \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

# demultiplex the sequences
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv
 
# import the metadata
qiime metadata tabulate \
  --m-input-file input/jay-met.tsv \
  --o-visualization tabulated-metadata.qzv