#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:02:00
#SBATCH --mem-per-cpu 1G
#SBATCH --job-name=4Parental
#SBATCH --output=./outputs/%x-%j.out


# cd /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/
module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2 #StdEnv/2016.4  

# run R script
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/CanadaJayMicrobiome/scripts/4Parent.R
