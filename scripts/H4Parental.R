print("Set up - working directory, packages, data")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(car)
library(qiime2R)
library(vegan)
library(phyloseq)
library(ggsignif)
library(tidyverse)
# set plotting theme
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
gj_meta <- as(sample_data(gj_ps), "data.frame") %>%
  mutate_all(na_if,"")
rownames(gj_meta) <- sample_names(gj_ps)

# boxplots - data prep (to long)
dmAitchison <- read_qza("H3aitchison-distance.qza")$data
# need to remove duplicate comparisons
dm_all <- dmAitchison %>% as.matrix %>% as.data.frame %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "AitchisonDistance") %>%
  # filter out duplicated comparisons (taling "lower" part of dm df)
  mutate(Sample1 = gsub("G", "", as.character(factor(.$Sample1))) %>% as.numeric(), 
         Sample2 = gsub("G", "", as.character(factor(.$Sample2))) %>% as.numeric()) %>%
  .[as.numeric(.$Sample1) > as.numeric(.$Sample2), ] %>%
  mutate(Sample1 = paste0("G", as.character(Sample1)), # Fixing sample namanes
         Sample2 = paste0("G", as.character(Sample2))) %>% # joing with sample data
  left_join(., rownames_to_column(gj_meta, var = "Sample1")) %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample2"), by = "Sample2")

print("Within territories (3A)")
# remove within territory comparisons
dm_between <- dm_all[-which(dm_all$Territory.x == dm_all$Territory.y),] %>%
  mutate(Territory = "Between", Group = "Between")
# only within territory groups
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  mutate(Territory = .$Territory.x, Group = "Within")
# put back together
dm_meta <- rbind(dm_between, dm_within) %>% select(Group, Territory, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/H3A.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Territory)) +
  geom_boxplot() + labs(x = "Territory", y = "Aitchison Distance") +
  scale_x_discrete(limits = c("Between", "SWAir", "DaviesBog", "CamLkRd",
                              "Mile36", "Arowhon", "NorthBog", "WolfHowl", "BatLake"))
dev.off()
# need to find a test for determining if the within group variation is < than between
# levene's test (car package)
print("Levene's test")
print("Within or between")
leveneTest(AitchisonDistance ~ Group, data = dm_meta)
print("Individual territories")
leveneTest(AitchisonDistance ~ Territory.x*Territory.y, data = dm_meta)
# clean up
rm(dm_between, dm_within, dm_meta)

print("Dominant Juveniles (3B)")
# only one dominant juvenile non-breeder
# breeder (any) to non-breeder (dominant) in same territory
dm_dj <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$JuvenileStatus.x == "DominantJuvenile" &
                        .$BreedingStatus.x == "Non-breeder"),] %>%
  select(BreedingStatus.x, BreedingStatus.y, Territory.x, Territory.y,
         JuvenileStatus.x, JuvenileStatus.y, everything()) %>%
  mutate(Group = "Dominant")
# non-breeder (not dominant juvenile) to breeder (any) in same territory
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$JuvenileStatus.x != "DominantJuvenile" &
            .$BreedingStatus.x == "Non-breeder"),] %>%
  select(BreedingStatus.x, BreedingStatus.y, Territory.x, Territory.y,
         JuvenileStatus.x, JuvenileStatus.y, everything()) %>%
  mutate(Group = "Not Dominant")
# put back together
dm_meta <- rbind(dm_dj, dm_within) %>% select(Group, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/H3B.pdf")
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Juvenile Status of Non-Breeder", y = "Aitchison Distance") +
  geom_signif(comparisons = list(c("Dominant", "Not Dominant")), 
              map_signif_level=TRUE, test = wilcox.test)
dev.off()
# wilcox-test
wilcox.test(dm_dj$AitchisonDistance, dm_within$AitchisonDistance)
# clean up
rm(dm_dj, dm_within, dm_meta)


print("Own offspring vs other offspring (3C)")
# breeders with non-breeders on same territory (does not include between breeders)
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.x == "Non-breeder"),] %>% # y's will be non-breeders
  mutate(Group = "Non-breeder (Same Territory)")
# breeders with breeders on same territory (does not include between non-breeders)
dm_withinB <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x == .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.x == "Breeder"),] %>%
  mutate(Group = "Breeder (Same Territory)")
# breeders with non-breeders on different territory (does not include between breeders)
dm_between <- dm_all[which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.x == "Non-breeder"),] %>% # y's will be non-breeders
  mutate(Group = "Non-breeder (Different Territory)")
# breeders with breeders on different territory (does not include between non-breeders)
dm_breeders <- dm_all[which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x == .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.x == "Breeder"),] %>%
  mutate(Group = "Breeder (Different Territory)")
# put back together
dm_meta <- rbind(dm_between, dm_within) %>%
  rbind(., dm_breeders) %>% rbind(., dm_withinB) %>%
  select(Group, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/H3C.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Focal Breeder Compared to...", y = "Aitchison Distance") +
  geom_signif(comparisons = list(c("Breeder (Different Territory)", "Breeder (Same Territory)"),
                                 c("Non-breeder (Different Territory)", "Non-breeder (Same Territory)")), 
              map_signif_level=TRUE, test = wilcox.test)
dev.off()
# Mann-Whitney test
print("Non-breeder Mann-Whitney")
wilcox.test(dm_between$AitchisonDistance, dm_within$AitchisonDistance)
print("Breeder Mann-Whitney")
wilcox.test(dm_breeders$AitchisonDistance, dm_withinB$AitchisonDistance)
# clean up
rm(dm_between, dm_within, dm_meta, dm_breeders, dm_withinB)
