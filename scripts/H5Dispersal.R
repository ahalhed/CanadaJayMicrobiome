# Testing Hypothesis 5 - host dispersal
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

## Calculate distance between origin and sampling locations
# this uses the Haversine to calculate half circle distance
oriDist <- function(samp) {
  # samp is the sample data
  df1 <- distm(samp %>% select(LongitudeSamplingDD, LatitudeSamplingDD),
               samp %>% select(LongitudeOriginDD, LatitudeOriginDD),
               fun = distHaversine) %>% as.data.frame
  # add extract labels to distances
  rownames(df1) <- rownames(samp)
  colnames(df1) <- rownames(samp)
  # retain only calculation of distance between origin and same sample location
  df2 <- df1 %>% rownames_to_column(var = "SampleID") %>%
    # pivot to longer (proper formatting for plotting), keep only same sample distances
    pivot_longer(-SampleID, names_to = "SampleID2", values_to = "DistanceFromOrigin") %>%
    .[which(.$SampleID == .$SampleID2),] %>% select(-SampleID2) %>%
    # drop NAs from the Distance From Origin Column
    subset(!is.na(DistanceFromOrigin)) 
  return(df2)
}


## Prediction 5A
# Non-breeders who have travelled further from their origin will have more diverse microbiomes.
# considering the number of OTUs observed to be diversity here
# non-breeders are not established in a territory and are likely travelling around more
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P5A-filtered-table.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

# convert OTU table to presence/absences
OTUs <- otu_table(gj_ps) %>% t %>% as.data.frame()
OTUs[which(OTUs>0, arr.ind=TRUE)] <- 1

oriDF <- gj_meta %>% mutate(OTUs = colSums(OTUs),
                   DistanceFromOrigin = oriDist(gj_meta)$DistanceFromOrigin)

# Make an ordination plot
pdf("CanadaJayMicrobiome/plots/P5A.pdf", width = 10)
ggplot(oriDF, aes(y = OTUs, x = DistanceFromOrigin, shape = as.factor(CollectionYear))) +
  labs(x = "Distance From Origin (m)",
       y = "Number of OTUs Observed",
       shape = "Collection Year") +
  geom_point() +
  scale_color_viridis_d()
dev.off()

# stats on linear model
# linear regression on all samples
print("All years") # can't do year by year b/c of low sample number
lm(OTUs~CollectionYear+DistanceFromOrigin, data = oriDF) %>%
  summary
lm(OTUs~DistanceFromOrigin, data = oriDF) %>%
  summary

# clean up
rm(gj_meta, gj_ps, OTUs, oriDF)

## Prediction 5B
# Breeders established in a specific territory closer to their natal territory will have less diverse microbial communities.
gj_ps <- qza_to_phyloseq(features = "P5B-filtered-table.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

# convert OTU table to presence/absences
OTUs <- otu_table(gj_ps) %>% t %>% as.data.frame()
OTUs[which(OTUs>0, arr.ind=TRUE)] <- 1

oriDF <- gj_meta %>% mutate(OTUs = colSums(OTUs),
                            DistanceFromOrigin = oriDist(gj_meta)$DistanceFromOrigin)

# all fall samples
pdf("CanadaJayMicrobiome/plots/P5B.pdf", width = 10)
ggplot(oriDF, aes(x = DistanceFromOrigin, y = OTUs,
                      shape = as.factor(CollectionYear))) +
  geom_point() + #log scale?
  labs(y = "Number of OTUs", shape = "Collection Year",
       x = "Distance from Origin (m)")
dev.off()

# stats on linear models
# linear regression on all samples
print("All years")
lm(OTUs~CollectionYear+DistanceFromOrigin, data = oriDF) %>%
  summary
lm(OTUs~DistanceFromOrigin, data = oriDF) %>%
  summary
# 2017
print("2017 Only")
lm(OTUs~DistanceFromOrigin, 
   data=filter(oriDF,CollectionYear == 2017)) %>%
  summary
# 2018
print("2018 Only")
lm(OTUs~DistanceFromOrigin, 
   data=filter(oriDF,CollectionYear == 2018)) %>%
  summary
# 2020
print("2020 Only")
lm(OTUs~DistanceFromOrigin, 
   data=filter(oriDF,CollectionYear == 2020)) %>%
  summary

# clean up
rm(gj_meta, gj_ps, OTUs, oriDF)
