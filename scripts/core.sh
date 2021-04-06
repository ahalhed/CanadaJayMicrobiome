#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:12:00
#SBATCH --mem-per-cpu 1G
#SBATCH --job-name=core
#SBATCH --output=outputs/%x-%j.out


# cd /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/
module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2

# run R script
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/CanadaJayMicrobiome/scripts/core.R
# run by season
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/CanadaJayMicrobiome/scripts/SeasonalCores/coreF17.R
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/CanadaJayMicrobiome/scripts/SeasonalCores/coreF18.R
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/CanadaJayMicrobiome/scripts/SeasonalCores/coreS20.R
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/CanadaJayMicrobiome/scripts/SeasonalCores/coreF20.R