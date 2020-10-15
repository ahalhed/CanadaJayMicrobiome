# for Hypothesis 1 - related to environmental/spatial differences
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
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-singletons-mitochondria-chloroplast.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met2.tsv")
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the aitchison distance matrix
ordiAitchison <- read_qza("aitchison-ordination.qza")
# join aitchison vectors with metadata
# select columns of interest from meta (no location or extraction info)
gj_aitch_V <- gj_meta %>% select(1:5, 7:19, 28:33) %>%
  rownames_to_column(var = "SampleID") %>% # row names need to be a column to join
  left_join(ordiAitchison$data$Vectors,.) %>%
  # drop G6, because it is the blank (not needed)
  .[ which(.$SampleID != "G6"), ] %>%
  remove_rownames()

## Territory Figure
# Corresponds to figure 3-1 from proposal. 

# plot PC1/PC2 by territory
pdf("CanadaJayMicrobiome/plots/H1pc.pdf")
ggplot(gj_aitch_V, aes(y=PC2, x=PC1, shape = as.factor(CollectionYear), group = JayID)) + 
  geom_point() + geom_line() + # points are samples
  facet_grid(~Territory) + labs(shape = "Collection Year")
dev.off()

# considering also averaging the samples by territory (so fewer panels)
aitchSum <- gj_aitch_V %>% group_by(CollectionYear, CollectionSeason, Territory) %>%
  # averages all individuals collected on the same territory in same year/season
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3))
# Cleveland dot plot by territory could be a different way to look at these
# see section 13.7 in R Graphics cookbook for more information
pdf("CanadaJayMicrobiome/plots/H1cleveland.pdf")
ggplot(aitchSum, aes(x=PC1, y=Territory, size = as.factor(CollectionYear), shape = CollectionSeason)) + 
  geom_point() + labs(size = "Collection Year", shape = "Collection Season")
dev.off()

# add PCNM spatial analysis here