# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(qiime2R)
library(vegan)
library(phyloseq)
library(ggsignif)
library(tidyverse)
# set plotting theme
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

# boxplots - data prep (to long)
dmAitchison <- read_qza("H3-aitchison-distance.qza")$data
dm_all <- dmAitchison %>% as.matrix %>% as.data.frame %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "AitchisonDistance") %>%
  # remove same-sample comparisons
  .[-which(.$Sample1 == .$Sample2),] %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample1")) %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample2"), by = "Sample2")

# within territories (3A)
# remove within territory comparisons
dm_between <- dm_all[-which(dm_all$Territory.x == dm_all$Territory.y),] %>%
  mutate(Group = "Between")
# only within territory groups
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  mutate(Group = .$Territory.x)
# put back together
dm_meta <- rbind(dm_between, dm_within) %>% select(Group, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/H3A.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Territory", y = "Aitchison Distance") +
  scale_x_discrete(limits = c("Between", "SWAir", "DaviesBog", "CamLkRd",
                              "Mile36", "Arowhon", "NorthBog", "WolfHowl", "BatLake"))
dev.off()
# clean up
rm(dm_between, dm_within, dm_meta)
# need to find a test for determining if the within group variation is < than between
# permanova
adonis2(dmAitchison ~ Territory + JuvenileStatus + BreedingStatus, data = offPar)

# dominant (3B)
# only one dominant juvenile non-breeder
# breeder (any) to non-breeder (dominant) in same territory
dm_dj <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$JuvenileStatus.y == "DominantJuvenile" &
                        .$BreedingStatus.y == "Non-breeder"),] %>%
  select(BreedingStatus.x, BreedingStatus.y, Territory.x, Territory.y,
         JuvenileStatus.x, JuvenileStatus.y, everything()) %>%
  mutate(Group = "Dominant")
# non-breeder (not dominant juvenile) to breeder (any) in same territory
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$JuvenileStatus.y != "DominantJuvenile" &
            .$BreedingStatus.y == "Non-breeder"),] %>%
  select(BreedingStatus.x, BreedingStatus.y, Territory.x, Territory.y,
         JuvenileStatus.x, JuvenileStatus.y, everything()) %>%
  mutate(Group = "Not Dominant")
# put back together
dm_meta <- rbind(dm_dj, dm_within) %>% select(Group, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/H3B.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Juvenile Status of Non-Breeder", y = "Aitchison Distance") +
  geom_signif(comparisons = list(c("Dominant", "Not Dominant")), 
              map_signif_level=TRUE)
dev.off()
# clean up
rm(dm_dj, dm_within, dm_meta)

# own offspring vs other offspring (3C)
# breeders with non-breeders on same territory (does not include between breeders)
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.y == "Non-breeder"),] %>% # y's will be non-breeders
  select(BreedingStatus.x, BreedingStatus.y, Territory.x, Territory.y,
         JuvenileStatus.x, JuvenileStatus.y, everything()) %>%
  mutate(Group = "Same Territory")
# breeders with non-breeders on different territory (does not include between breeders)
dm_between <- dm_all[which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.y == "Non-breeder"),] %>% # y's will be non-breeders
  select(BreedingStatus.x, BreedingStatus.y, Territory.x, Territory.y,
         JuvenileStatus.x, JuvenileStatus.y, everything()) %>%
  mutate(Group = "Different Territory")
# put back together
dm_meta <- rbind(dm_between, dm_within) %>% select(Group, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/H3C.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Territory of Non-Breeder", y = "Aitchison Distance") +
  geom_signif(comparisons = list(c("Different Territory", "Same Territory")), 
              map_signif_level=TRUE)
dev.off()
# clean up
rm(dm_dj, dm_within, dm_meta)