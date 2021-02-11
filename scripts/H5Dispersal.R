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
gj_ps <- qza_to_phyloseq(features = "H4filtered-table.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the ordination aitchison matrix (only samples with origins)
ordiAitchison <- read_qza("H4aitchison-ordination.qza")

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

# select columns of interest from meta (no location or extraction info)
ordiDF <- gj_meta %>%
  rownames_to_column(var = "SampleID") %>%
  # distance from origin and principle component axes
  full_join(oriDistCol) %>% full_join(ordiAitchison$data$Vectors)

## Prediction 4A
# Make an ordination plot
pdf("CanadaJayMicrobiome/plots/H4pc.pdf")
ggplot(ordiDF, aes(x = PC1, y = PC2, size = DistanceFromOrigin, color = AgeAtCollection)) +
  geom_point() + labs(color = "Age at Collection", size = "Distance from Origin (m)") +
  coord_fixed(ylim = c(-0.35, 1), xlim = c(-0.35, 1)) +
  #stat_ellipse() +
  scale_color_viridis_c()
dev.off()

# permanova
dmAitchison <- read_qza("H4aitchison-distance.qza")$data
adonis2(dmAitchison ~ DistanceFromOrigin + JayID + AgeAtCollection + CollectionYear,
        data = ordiDF)

## Prediction 4B
# Pairing off samples from the same collection year. 
# Determining least distance from origin of one sample to current location of another
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

# stats on linear model
# linear regression on all samples
print("All years")
lm(AitchisonDistance~DistanceFromOrigin, data = plotDFall) %>%
  summary
# 2017
print("2017 Only")
lm(AitchisonDistance~DistanceFromOrigin, 
   data=filter(plotDFall,CollectionYear == 2017)) %>%
  summary
# 2018
print("2018 Only")
lm(AitchisonDistance~DistanceFromOrigin, 
   data=filter(plotDFall,CollectionYear == 2018)) %>%
  summary
# 2020
print("2020 Only")
lm(AitchisonDistance~DistanceFromOrigin, 
   data=filter(plotDFall,CollectionYear == 2020)) %>%
  summary