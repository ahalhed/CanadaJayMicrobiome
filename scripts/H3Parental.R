setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(vegan)
library(phyloseq)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-singletons-mitochondria-chloroplast-cr-99.qza", 
                         tree = "trees/rooted-tree-cr-99.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy-cr-99.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame") %>%
  mutate_all(na_if,"")
rownames(gj_meta) <- sample_names(gj_ps)

# select only offspring with parent information
offspring <- gj_meta %>% rownames_to_column('sampleid') %>%
  subset(!is.na(SampledBreederMale) | !is.na(SampledBreederFemale))
# retain parents identified in nestlings rows
parents <- gj_meta %>% rownames_to_column('sampleid') %>%
  subset(JayID %in% unique(offspring$SampledBreederMale) | JayID %in% unique(offspring$SampledBreederFemale)) %>%
  subset(CollectionYear %in% unique(offspring$CollectionYear) & CollectionSeason %in% unique(offspring$CollectionSeason))
# join parent and offspring together
offPar <- full_join(offspring, parents) %>%
  remove_rownames() %>% column_to_rownames('sampleid')
# subset the phyloseq object accordingly
# this needs final testing since we don't yet have sequence data for all samples
offParPS <- prune_samples(sample_names(gj_ps) %in% rownames(offPar), gj_ps)
# clean up
rm(gj_ps, gj_meta, parents, offspring)

# read in Aitchison ordination data
ordiAitchison <- read_qza("aitchison-ordination-cr-99-wmc.qza")$data$Vectors %>%
  # subset to include only samples in phyloseq object
  # this needs final testing since we don't yet have sequence data for all samples
  filter(SampleID %in% sample_names(offParPS)) %>%
  # combined with filtered metadata
  # this needs final testing since we don't yet have sequence data for all samples
  full_join(rownames_to_column(offPar, var = "SampleID"))
# plot ordination data
# territory is a proxy for nest group
ggplot(ordiAitchison, aes(x = PC1, y = PC2, group = Territory, linetype = Territory)) +
  geom_point() + stat_ellipse(type = "norm")
#working on the appearance of this one