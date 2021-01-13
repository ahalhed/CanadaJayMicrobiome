# for Hypothesis 1 - related to environmental/spatial differences
# this script is for the environmental analysis
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(vegan)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the aitchison ordination
ordiAitchison <- read_qza("aitchison-ordination.qza")
# join aitchison vectors with metadata
# select columns of interest from meta (no location or extraction info)
gj_aitch_V <- gj_meta %>% select(1:5, 7:17, 24:27) %>%
  rownames_to_column(var = "SampleID") %>% # row names need to be a column to join
  left_join(ordiAitchison$data$Vectors,.) %>%
  remove_rownames()
# read in aitchison distance matrix
# aitchison distances
dmAitchison <- read_qza("aitchison-distance.qza")

## Territory Figure
# Corresponds to figure 3-1 from proposal. 
# Version 1 - plot PC1/PC2 by territory
pdf("CanadaJayMicrobiome/plots/H1pc.pdf", width = 10)
# need to sort out NA thing
ggplot(gj_aitch_V, aes(y=PC2, x=PC1, shape = as.factor(CollectionYear), group = JayID)) + #, group = JayID
  geom_point() + #geom_line()
  labs(shape = "Collection Year")
dev.off()
# considering also averaging the samples by territory (so fewer panels)
aitchSum <- gj_aitch_V %>% group_by(CollectionYear, CollectionSeason, Territory) %>%
  # averages all individuals collected on the same territory in same year/season
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3))
# Version 2 - Cleveland dot plot by territory could be a different way to look at these
# see section 13.7 in R Graphics cookbook for more information
pdf("CanadaJayMicrobiome/plots/H1cleveland.pdf")
# need to sort out NA season
ggplot(aitchSum, aes(x=PC1, y=Territory, size = CollectionYear, shape = CollectionSeason)) + 
  geom_point() + labs(size = "Collection Year", shape = "Collection Season")
dev.off()


# Vegan-based analysis starts here
# distance based RDA using aitchison distance matrix
# should throw in collection season when there are more samples
# might switch to the phyloseq implementation
gj_cap <- capscale(dmAitchison$data ~ ProportionSpruceOnTerritory + MeanTempC + CollectionSeason + CollectionYear,
                   data = gj_meta, comm = otu_table(gj_ps), na.action = na.exclude)
# look at summaries
summary(gj_cap)
# simple biplot
pdf("CanadaJayMicrobiome/plots/H1envBiplot.pdf")
plot(gj_cap, main = "Aitchison Distance-based RDA")
dev.off()

gj_meta %>% mutate(across(everything(), ~replace_na(.x, '')))
# PERMANOVA
adonis2(dmAitchison$data ~ TerritoryQuality + MeanTempC + CollectionSeason + CollectionYear + JayID,
        data = gj_meta, na.action = na.exclude)
