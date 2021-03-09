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


print("Prediction 5A")
# Breeders established in a specific territory closer to their natal territory will have less diverse microbial communities.
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

# save plot
pdf("CanadaJayMicrobiome/plots/P5A.pdf", width = 10)
ggplot(oriDF, aes(x = DistanceFromOrigin, y = OTUs, shape = as.factor(CollectionYear))) +
  geom_point() + scale_x_log10() + geom_smooth(method=lm, se=FALSE, color = "black") +
  labs(shape = "Collection Year", y = "Number of OTUs Observed",
       x = "Distance from Origin (m, log10)")
dev.off()

# stats on linear models
# linear regression on all samples
print("All years")
Yall <- (lm(OTUs~DistanceFromOrigin, data = oriDF))
summary(Yall)
anova(Yall)
rm(Yall)
# 2017
print("2017 Only")
(Y2017 <- lm(OTUs~DistanceFromOrigin, 
   data=filter(oriDF,CollectionYear == 2017)))
summary(Y2017)
anova(Y2017)
rm(Y2017)
# 2018
print("2018 Only")
(Y2018 <- lm(OTUs~DistanceFromOrigin, 
             data=filter(oriDF,CollectionYear == 2018)))
summary(Y2018)
anova(Y2018)
rm(Y2018)
# 2020
print("2020 Only")
(Y2020 <- lm(OTUs~DistanceFromOrigin, 
             data=filter(oriDF,CollectionYear == 2020)))
summary(Y2020)
anova(Y2020)
rm(Y2020)

# clean up
rm(gj_meta, gj_ps, oriDF, OTUs, oriDist)
