# for Hypothesis 1 - related to seasonal differences

print("Working directory, packages, functions")
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(qiime2R)
library(phyloseq)
library(tidyverse)

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)


print("Should we do analyses by year/season?")
# confirmation that year/season was the right way to go
print("Collection Year Ordination")
(orY <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionYear))
anova(orY)
print("Collection Season Ordination")
(orS <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionSeason))
anova(orS)
print("Collection Year*Season Ordination")
(orSY <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionYear * CollectionSeason))
anova(orSY)
# clean up
rm(gj_meta, gj_ps, orS, orSY, orY)