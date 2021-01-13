# for Hypothesis 2 - feeding differences
# this script is for the diet supplementation plot
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(vegan)
library(ggsignif)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the aitchison distance matrix
dmAitchison <- read_qza("aitchison-distance.qza")
# getting ready for boxplot
# combine the aitchison distance data with metadata data
dm_meta <- dmAitchison$data %>% as.matrix %>% as.data.frame %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "AitchisonDistance") %>%
  # remove same-sample comparisons
  .[-which(.$Sample1 == .$Sample2),] %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample1")) %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample2"), by = "Sample2") %>%
  # add shared food information
  mutate(host = ifelse(.$JayID.x == .$JayID.y, "Same Bird", "Different Bird"),
         group = ifelse(.$`FoodSupplement.x` == "Y" & .$`FoodSupplement.y` == "Y", "Yes",
                        ifelse(.$`FoodSupplement.x` == "N" & .$`FoodSupplement.y` == "N",
                               "No", "Either"))) %>%
  # select only the variables of interest for the boxplot
  select(1:4, ExtractID.y, host, group)

# save boxplot to PDF
pdf("CanadaJayMicrobiome/plots/H2foodBox.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = group)) +
  geom_boxplot() + labs(x = "Food Supplementation")
dev.off()

# read in the aitchison ordination (probs won't keep this)
ordiAitchison <- read_qza("aitchison-ordination.qza")
# combine aitchison vectors with environmental data
aitch <- gj_meta %>% rownames_to_column(var = "SampleID") %>%
  full_join(ordiAitchison$data$Vectors)

# save ordination plot to PDF
pdf("CanadaJayMicrobiome/plots/H2foodOrdi.pdf", width = 9)
# ellipse by territory
# shape by food
# too few samples to get much meaningful from this
ggplot(aitch, aes(x=PC1, y=PC2, shape = FoodSupplement, linetype = FoodSupplement)) + 
  geom_point() + # points are samples
  stat_ellipse(type = "norm") + 
  coord_fixed(ylim = c(-0.6, 0.4), xlim = c(-0.6, 0.4)) +
  labs(shape = "Food Supplementation", linetype = "Food Supplementation")
# turn off graphics device
dev.off()

# PERMANOVA
adonis2(dmAitchison$data ~ FoodSupplement + JayID + AgeAtCollection + CollectionYear,
        data = gj_meta)
