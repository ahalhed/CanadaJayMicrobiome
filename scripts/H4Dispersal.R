# Testing Hypothesis 4 - host dispersal
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
## Script set up
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(geosphere)
library(vegan)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the ordination aitchison matrix (only samples with origins)
ordiAitchison <- read_qza("H4-aitchison-ordination.qza")
ordiAitchison$data$Species
## Calculate distance between origin and sampling locations
# this uses the Haversine to calculate half circle distance
oriDist <- distm(gj_meta %>% select(LongitudeSamplingDD, LatitudeSamplingDD), 
      gj_meta %>% select(LongitudeOriginDD, LatitudeOriginDD), 
      fun = distHaversine) %>% as.data.frame
# add extract labels to distances
rownames(oriDist) <- rownames(gj_meta) 
colnames(oriDist) <- rownames(gj_meta)
# retain only calculation of distance between origin and same sample location
oriDistCol <- oriDist %>% rownames_to_column(var = "SampleID") %>%
  # pivot to longer (proper formatting for plotting), keep only same sample distances
  pivot_longer(-SampleID, names_to = "SampleID2", values_to = "DistanceFromOrigin") %>%
  .[which(.$SampleID == .$SampleID2),] %>% select(-SampleID2) %>%
  # drop NAs from the Distance From Origin Column
  subset(!is.na(DistanceFromOrigin)) 

## Combine data from 3 above sources to plots
# select columns of interest from meta (no location or extraction info)
ordiDF <- gj_meta %>%
  rownames_to_column(var = "SampleID") %>%
  # distance from origin and principle component axes
  full_join(oriDistCol) %>% full_join(ordiAitchison$data$Vectors)


## Make an ordination plot
pdf("CanadaJayMicrobiome/plots/H4DistanceFromOrigin.pdf")
ggplot(ordiDF, aes(x = PC1, y = PC2, color = DistanceFromOrigin, size = AgeAtCollection)) +
  geom_point() + labs(size = "Age at Collection", shape = "Distance from Origin (m)") +
  coord_fixed(ylim = c(-0.3, 1), xlim = c(-0.3, 1)) + scale_color_viridis_c()
dev.off()


ggplot(ordiDF, aes(x=PC1, size= DistanceFromOrigin, y = AgeAtCollection)) + 
  geom_point() + labs(size = "Distance From Origin", y = "Age at Collection")

# permanova
dmAitchison <- read_qza("H4-aitchison-distance.qza")$data
adonis2(dmAitchison ~ DistanceFromOrigin + JayID + AgeAtCollection + CollectionYear,
        data = ordiDF)