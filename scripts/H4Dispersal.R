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
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-singletons-mitochondria-chloroplast.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met2.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the ordination aitchison matrix
ordiAitchison <- read_qza("aitchison-ordination.qza")
# read in the aitchison distance matrix
dmAitchison <- read_qza("aitchison-distance.qza")

## Calculate distance between origin and sampling locations
# this uses the Haversine to calculate half circle distance
# longitude is labeled x, latitude labeled y
oriDist <- distm(gj_meta %>% select(XsampDD, YsampDD), 
      gj_meta %>% select(XoriginDD, YoriginDD), 
      fun = distHaversine) %>% as.data.frame
# add extract labels to distances
rownames(oriDist) <- rownames(gj_meta) 
colnames(oriDist) <- rownames(gj_meta)
# retain only calculation of distance between origin and same sample location
oriDistCol <- oriDist %>% rownames_to_column(var = "SampleID") %>%
  # drop G6, because it is the blank (not needed)
  .[ which(.$SampleID != "G6"), ] %>%
  # pivot to longer (proper formatting for plotting), keep only same sample distances
  pivot_longer(-SampleID, names_to = "SampleID2", values_to = "DistanceFromOrigin") %>%
  .[which(.$SampleID == .$SampleID2),] %>% select(-SampleID2) %>%
  # drop NAs from the Distance From Origin Column
  subset(!is.na(DistanceFromOrigin)) 
# filter the distance matrix to only contain samples of interest
distAitch <- dmAitchison$data %>% as.matrix %>% as.data.frame %>%
  rownames_to_column(var = "SampleID") %>%
  pivot_longer(-SampleID, names_to = "SampleID2", values_to = "AitchisonDistance") %>%
  # remove same-sample comparisons
  .[-which(.$SampleID == .$SampleID2),] %>%
  # select only the samples that we have harvesine distance from origin for
  .[which(.$SampleID %in% oriDistCol$SampleID & .$SampleID2 %in% oriDistCol$SampleID),] %>%
  # PIVOT wide
  pivot_wider(id_cols = SampleID, names_from = SampleID2, values_from = AitchisonDistance) %>%
  remove_rownames %>% column_to_rownames(var = "SampleID") %>%
  # replace NA with 0 and convert back to distance matrix
  mutate_all(~replace(., is.na(.), 0)) %>% as.matrix %>% as.dist()
View(distAitch)
## Combine data from 3 above sources to plots
# select columns of interest from meta (no location or extraction info)
ordiDF <- gj_meta %>% select(1:5, 7:19, 28:33) %>%
  rownames_to_column(var = "SampleID") %>%
  # distance from origin
  full_join(oriDistCol)
ordiDF
capscale(dmAitchison$data ~ Prop_Spruce_On_Territory + CollectionYear + AgeAtCollection,
         data = ordiDF, comm = otu_table(gj_ps), na.action = na.exclude)

## Make an ordination plot
pdf("CanadaJayMicrobiome/plots/H4DistanceFromOrigin.pdf")
ggplot(ordiDF, aes(x = PC1, y = PC2, color = DistanceFromOrigin, size = AgeAtCollection)) +
  geom_point() + labs(size = "Age at Collection", shape = "Distance from Origin (m)")
dev.off()
