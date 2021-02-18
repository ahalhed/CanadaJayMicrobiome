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
pdf("CanadaJayMicrobiome/plots/P5A.pdf")
ggplot(oriDF, aes(y = OTUs, x = DistanceFromOrigin,
                  shape = as.factor(CollectionYear))) +
  labs(x = "Distance From Origin (m)",
       y = "Number of OTUs Observed",
       shape = "Collection Year") +
  geom_point() +
  scale_color_viridis_d()
dev.off()

# stats on linear model
# linear regression on all samples
print("All years") # can't do year by year b/c of low sample number
lm(OTUs~DistanceFromOrigin, data = oriDF) %>%
  summary

# permanova
dmAitchison <- read_qza("P5A-aitchison-distance.qza")$data
adonis2(dmAitchison ~ DistanceFromOrigin + CollectionYear,
        data = oriDF)

# clean up
rm(gj_meta, gj_ps, OTUs, oriDF, dmAitchison)

## Prediction 5B
# Breeders established in a specific territory will have microbial communities unique from their natal territory.
mixed <- distm(gj_meta[rownames(gj_meta) %in% oriDistCol$SampleID,] %>%
                 select(LongitudeSamplingDD, LatitudeSamplingDD), 
               gj_meta[rownames(gj_meta) %in% oriDistCol$SampleID,] %>%
                 select(LongitudeOriginDD, LatitudeOriginDD), 
               fun = distHaversine) %>% as.data.frame
# add extract labels to distances
rownames(mixed) <- rownames(gj_meta[rownames(gj_meta) %in% oriDistCol$SampleID,])
colnames(mixed) <- rownames(gj_meta[rownames(gj_meta) %in% oriDistCol$SampleID,])
# pivot to longer (proper formatting for plotting)
mixedLong <- mixed %>% rownames_to_column(var = "SampleID") %>%
  pivot_longer(-SampleID, names_to = "SampleID2", values_to = "DistanceFromOrigin") %>%
  .[which(.$SampleID != .$SampleID2),]
# collect only metadata of interest
meta_sub <- gj_meta %>%
  select(CollectionYear, CollectionSeason, Territory, JayID) %>%
  rownames_to_column(var = "SampleID") %>%
  .[.$SampleID %in% oriDistCol$SampleID,]
# samples collected within the same year and season only
# sample data
mixedMeta <- left_join(mixedLong, meta_sub) %>%
  left_join(meta_sub, suffix = c("", "2"), by = c("SampleID2" = "SampleID")) %>%
  .[.$CollectionYear==.$CollectionYear2,] %>%
  .[.$CollectionSeason==.$CollectionSeason2,] %>%
  .[.$CollectionSeason=="Fall",] %>% # omitting spring b/c of nest thing
  select(JayID, JayID2, DistanceFromOrigin, CollectionYear,SampleID,SampleID2)
# distance data
AitchisonLong <- dmAitchison %>% 
  as.matrix %>% as.data.frame %>% 
  rownames_to_column(var = "SampleID") %>%
  # pivot to longer (proper formatting for plotting)
  pivot_longer(-SampleID, names_to = "SampleID2", values_to = "AitchisonDistance")
# combined
plotDFall <- left_join(AitchisonLong, mixedMeta, by = c("SampleID", "SampleID2")) %>%
  .[!is.na(.$JayID),]

# all distances within the same year and season only
pdf("CanadaJayMicrobiome/plots/H4B.pdf")
ggplot(plotDFall, aes(x = DistanceFromOrigin, y = AitchisonDistance,
                      colour = CollectionYear)) +
  geom_point() + scale_color_viridis_c() +
  geom_smooth(method='lm', formula= y~x, color="black") +
  #facet_grid(~CollectionYear) +
  labs(y = "Aitchison Distance",
       x = "Distance from Origin of Focal Individual (m)")
dev.off()


# 2020
print("2020 Only")
lm(AitchisonDistance~DistanceFromOrigin, 
   data=filter(plotDFall,CollectionYear == 2020)) %>%
  summary