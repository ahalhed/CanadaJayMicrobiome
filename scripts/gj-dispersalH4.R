# Testing Hypothesis 4 - host dispersal

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
                         taxonomy = "taxonomy/silva-taxonomy.qza",
                         metadata = "input/jay-met2.tsv")
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the aitchison distance matrix
ordiAitchison <- read_qza("aitchison-ordination.qza")

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
  pivot_longer(-SampleID, names_to = "SampleID2", values_to = "DistanceFromOrigin") %>%
  .[which(.$SampleID == .$SampleID2),] %>% select(-SampleID2)

## Combine data from 3 above sources to plots
# select columns of interest from meta (no location or extraction info)
ordiDF <- gj_meta %>% select(1:5, 7:19, 28:33) %>%
  rownames_to_column(var = "SampleID") %>%
  # combine metadata with ordination vectors and distance from origin
  full_join(ordiAitchison$data$Vectors) %>% full_join(oriDistCol)

## Make an ordination plot
pdf("CanadaJayMicrobiome/plots/DistanceFromOrigin.pdf")
ggplot(ordiDF, aes(x = PC1, y = PC2, color = DistanceFromOrigin, size = AgeAtCollection)) +
  geom_point()
dev.off()
