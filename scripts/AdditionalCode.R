# additional code that don't necessarily fit into one hypothesis
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/gj_meta-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

print("counting")
print("breeding status")
gj_meta %>%
  group_by(BreedingStatus) %>%
  count()
print("breeding status and juvenile status")
gj_meta %>%
  group_by(BreedingStatus, JuvenileStatus) %>%
  count()
print("breeding status, juvenile status, and collection season")
gj_meta %>%
  group_by(BreedingStatus, JuvenileStatus, CollectionSeason) %>%
  count()
print("breeding status, juvenile status, collection season, and collection year")
gj_meta %>%
  group_by(BreedingStatus, JuvenileStatus, CollectionSeason, CollectionYear) %>%
  count()
print("Breeders by territory, year, season")
gj_meta %>% filter(BreedingStatus == "Breeder") %>%
  group_by(Territory, CollectionYear, CollectionSeason) %>%
  count()
print("Non-breeders by territory, year, season")
gj_meta %>% filter(BreedingStatus == "Non-breeder") %>%
  group_by(Territory, CollectionYear, CollectionSeason) %>%
  count()

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