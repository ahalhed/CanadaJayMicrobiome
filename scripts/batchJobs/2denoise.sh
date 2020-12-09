#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=02:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=DADA2
#SBATCH --dependency=afterok:
#SBATCH --output=CanadaJayMicrobiome/output/%x-%j.out

#script starts here
#----------------------------------
# depends on 1import.sh
# going to trim at 250 bp
# this took just over two hours to run
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 250 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv