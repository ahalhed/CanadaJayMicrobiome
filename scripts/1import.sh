#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=import
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# import the sequence data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./rawFiles \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ./demux-paired-end.qza
