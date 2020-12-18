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
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame") %>%
  mutate_all(na_if,"")
rownames(gj_meta) <- sample_names(gj_ps)

# select only offspring with parent information
offspring <- gj_meta %>% rownames_to_column('sampleid') %>%
  subset(CollectionYear == 2020) %>% subset(CollectionSeason == "Spring")
# retain parents identified in nestlings rows
parents <- gj_meta %>% rownames_to_column('sampleid') %>%
  subset(JayID %in% unique(offspring$SampledBreederMale) | JayID %in% unique(offspring$SampledBreederFemale)) %>%
  subset(CollectionYear == 2020) %>% subset(CollectionSeason == "Spring")
# join parent and offspring together
offPar <- full_join(offspring, parents) %>%
  # dropping offspring where parent was only sampled in past seasons
  #subset(Territory != "SundayCreek") %>%
  remove_rownames() %>% column_to_rownames('sampleid')
# subset the phyloseq object accordingly
# this needs final testing since we don't yet have sequence data for all samples
offParPS <- prune_samples(sample_names(gj_ps) %in% rownames(offPar), gj_ps)
# clean up
rm(gj_ps, parents, offspring)

# read in Aitchison ordination data
ordiAitchison <- read_qza("aitchison-ordination.qza")$data$Vectors %>%
  # subset to include only samples in phyloseq object
  # this needs final testing since we don't yet have sequence data for all samples
  filter(SampleID %in% sample_names(offParPS)) %>%
  # combined with filtered metadata
  # this needs final testing since we don't yet have sequence data for all samples
  full_join(rownames_to_column(offPar, var = "SampleID"))

pdf("CanadaJayMicrobiome/plots/H3nestOrdi.pdf", width = 9)
# plot ordination data
# territory is a proxy for nest group
ggplot(ordiAitchison, aes(x = PC1, y = PC2, colour = Territory)) + # , group = Territory, linetype = Territory
  geom_point() + stat_ellipse(type = "norm") +
  scale_color_viridis_d()
# working on the appearance of this one
dev.off()

dmAitchison <- read_qza("H3-aitchison-distance.qza")$data
dm_meta <- dmAitchison %>% as.matrix %>% as.data.frame %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "AitchisonDistance") %>%
  # remove same-sample comparisons
  .[-which(.$Sample1 == .$Sample2),] %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample1")) %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample2"), by = "Sample2") %>%
  # remove across territory comparisons
  .[-which(.$Territory.x == .$Territory.y),]
# boxplot
pdf("CanadaJayMicrobiome/plots/H3parentalBox.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Territory.x)) +
  geom_boxplot() + labs(x = "Territory")
dev.off()
# permanova
adonis2(dmAitchison ~ Territory + JuvenileStatus + BreedingStatus,
        data = offPar)
