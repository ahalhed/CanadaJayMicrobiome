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
# metadata will all samples (temp fix)
#jay_ALL <- readxl::read_excel("~/OneDrive - University of Guelph/Alicia's Thesis/Grey Jays/metadata/jay-ALL.xlsx")
# drop rows juveniles without parent(s) sampled from same season
# jay_ALL[-1,] %>% 
offspring <- gj_meta %>% rownames_to_column('sampleid') %>%
  subset(!is.na(Sampled_Parent_Male) | !is.na(Sampled_Parent_Female)) %>% 
  subset(JayID != 'BLANK')
# retain parents identified in nestlings rows
parents <- gj_meta %>% rownames_to_column('sampleid') %>%
  subset(JayID %in% unique(offspring$Sampled_Parent_Male) | JayID %in% unique(offspring$Sampled_Parent_Female)) %>%
  subset(CollectionYear %in% unique(offspring$CollectionYear) & CollectionSeason %in% unique(offspring$CollectionSeason))
# join parent and offspring together
offPar <- full_join(offspring, parents) %>%
  select(1, 3:5, 8, 11:20, 25, 27, 29:34) %>% # dropping irrelevant columns
  remove_rownames() %>% column_to_rownames('sampleid')
# subset the phyloseq object accordingly
# this needs final testing since we don't yet have sequence data for all samples
offParPS <- prune_samples(sample_names(gj_ps) %in% rownames(offPar), gj_ps)
# clean up
rm(gj_ps, gj_meta, parents, offspring)

# read in Aitchison distance matrix data
dmAitchison <- read_qza("aitchison-distance.qza")$data %>% as.matrix %>%
  # subset to include only samples in phyloseq object
  # this needs final testing since we don't yet have sequence data for all samples
  .[sample_names(offParPS), sample_names(offParPS)] %>% 
  as.dist() # converting back to distance matrix
# plot on dbRDA
# or could just use aitchison ordination from Q2