print("Set up - working directory, packages, functions, data")
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

#functions for analysis
longDM <- function(dm, metric, samp){
  # dm is a distance matrix of interest
  # metric is the name of the distance metric in matrix
  # samp is the data frame containing the sample data
  df1 <- dm %>% as.matrix %>% as.data.frame %>%
    rownames_to_column(var = "Sample1")
  df2 <- df1 %>% pivot_longer(-Sample1, names_to = "Sample2", values_to = metric)
  # filter out duplicated comparisons (taling "lower" part of dm df)
  df3 <- df2 %>%
    mutate(Sample1 = gsub("G", "",as.character(factor(.$Sample1))) %>% as.numeric(), 
           Sample2 = gsub("G", "", as.character(factor(.$Sample2))) %>% as.numeric()) %>%
    .[as.numeric(.$Sample1) > as.numeric(.$Sample2), ]
  df4 <- df3 %>% mutate(Sample1 = paste0("G", as.character(Sample1)), # Fixing sample names
                        Sample2 = paste0("G", as.character(Sample2)))
  df5 <- left_join(df4, rownames_to_column(samp, var = "Sample1")) %>% # joining with sample data
    left_join(., rownames_to_column(samp, var = "Sample2"), by = "Sample2")
  return(df5)
}
# in process - can't find phy within sample_names for whatever reason
share <- function(ps, b, nb) {
  # ps is the phyloseq object with all samples
  # b is the sample ID for the breeder
  # nb is the sample ID for the non-breeder
  phy1 <- subset_samples(ps, sample_names(ps) == b | sample_names(ps) == nb)
  phy2 <- filter_taxa(phy1, function(x) sum(x >= 1) == (2), TRUE)
  # breeder only
  breeder <- subset_samples(phy1, sample_names(phy1) == b)
  breeder2 <- filter_taxa(breeder, function(x) sum(x >= 1) == (1), TRUE)
  # non-breeder only
  nbreed <- subset_samples(phy1, sample_names(phy1) == nb)
  nbreed2 <- filter_taxa(nbreed, function(x) sum(x >= 1) == (1), TRUE)
  # get the number of OTUs
  shared <- otu_table(phy2) %>% ncol
  Bonly <- otu_table(breeder2) %>% ncol
  NBonly <- otu_table(nbreed2) %>% ncol
  df <- data.frame(Bonly, shared, NBonly)
  return(df)
}
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P4ABCD-filtered-table.qza",
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
dmAitchison <- read_qza("P4ABCD-aitchison-distance.qza")$data
# need to remove duplicate comparisons
dm_all <- longDM(dmAitchison, "AitchisonDistance", gj_meta)

print("Dominant Juveniles (4A)")
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
pdf("CanadaJayMicrobiome/plots/P4A.pdf")
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Juvenile Status of Non-Breeder", y = "Aitchison Distance") +
  geom_signif(comparisons = list(c("Dominant", "Not Dominant")), 
              map_signif_level=TRUE, test = wilcox.test)
dev.off()
# wilcox-test
wilcox.test(dm_dj$AitchisonDistance, dm_within$AitchisonDistance)
# clean up
rm(dm_dj, dm_within, dm_meta)

print("Own offspring vs other offspring (4B)")
# breeders with non-breeders on same territory (does not include between breeders)
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.x == "Non-breeder"),] %>% # y's will be non-breeders
  mutate(Group = "Same Territory")
# breeders with non-breeders on different territory (does not include between breeders)
dm_between <- dm_all[which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$BreedingStatus.x != .$BreedingStatus.y),] %>%
  .[which(.$BreedingStatus.x == "Non-breeder"),] %>% # y's will be non-breeders
  mutate(Group = "Different Territory")
# put back together
dm_meta <- rbind(dm_between, dm_within) %>%
  select(Group, everything())
# save figure
pdf("CanadaJayMicrobiome/plots/P4B.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + labs(x = "Non-Breeder Location", y = "Aitchison Distance") +
  geom_signif(comparisons = list(c("Different Territory", "Same Territory")), 
              map_signif_level=TRUE, test = wilcox.test)
dev.off()
# Mann-Whitney test
print("Non-breeder Mann-Whitney")
wilcox.test(dm_between$AitchisonDistance, dm_within$AitchisonDistance)
# clean up
rm(dm_between, dm_within, dm_meta, dm_all, dmAitchison)

print("Counts of Genera (4C+D)")
# looking at the number of OTUs (not reads) that occur per sample
# could do reads using sample_sums
OTUs <- otu_table(gj_ps) %>% t %>% as.data.frame()
OTUs[which(OTUs>0, arr.ind=TRUE)] <- 1

print("Juvenile Diversity (4C)")
# count the number of genera that occur
plot4C <- colSums(OTUs) %>%
  as.data.frame() %>%
  rename("NumberOfOTUs" = ".") %>%
  # join with sample data
  merge(gj_meta, by=0, all=TRUE)

pdf("CanadaJayMicrobiome/plots/P4C.pdf")
# make a scatter plot
ggplot(plot4C, aes(y = NumberOfOTUs, x = BreedingStatus)) +
  geom_violin() + geom_jitter(width = 0.05, height = 0) +
  labs(y = "Number of OTUs", x = "Breeding Status")
dev.off()

# group summary
aggregate(plot4C$NumberOfOTUs, list(plot4C$BreedingStatus), summary)
# are the means different?
br <- plot4C[plot4C$BreedingStatus=="Breeder",] %>% .$NumberOfOTUs
nb <- plot4C[plot4C$BreedingStatus=="Non-breeder",] %>% .$NumberOfOTUs
wilcox.test(br,nb)

#clean up
rm(plot4C, nb, br, gentab, OTUs)

print("Shared Microbiota (4D)")
OTUs <- otu_table(gj_ps) %>% t %>% as.data.frame()
# get sample ID's
br <- gj_meta[which(gj_meta$BreedingStatus == "Breeder"),] %>% rownames()
nb <- gj_meta[which(gj_meta$BreedingStatus == "Non-breeder"),] %>% rownames()
# replace counts with rownames
rep <- which(OTUs>0, arr.ind=TRUE)
OTUs[rep] <- rownames(rep)
OTUs[OTUs == 0] <- NA

# get genera for sample ID's
otu_br <- OTUs[, colnames(OTUs) %in% br] %>% as.data.frame()
otu_nb <- OTUs[, colnames(OTUs) %in% nb] %>% as.data.frame()

# need to function or loop this
phy1 <- subset_samples(gj_ps, sample_names(gj_ps) == "G31" | sample_names(gj_ps) == "G48")
phy2 <- filter_taxa(phy1, function(x) sum(x >= 1) == (2), TRUE)
# breeder only
breeder <- subset_samples(phy1, sample_names(phy1) == "G31")
breeder2 <- filter_taxa(breeder, function(x) sum(x >= 1) == (1), TRUE)
# non-breeder only
nbreed <- subset_samples(phy1, sample_names(phy1) == "G48")
nbreed2 <- filter_taxa(nbreed, function(x) sum(x >= 1) == (1), TRUE)
# get the number of OTUs
otu_table(phy2) %>% ncol
otu_table(breeder2) %>% ncol
otu_table(nbreed2) %>% ncol

data.frame(Bonly, shared, NBonly)

# section clean up
rm(br, nb, rep, otu_br, otu_nb, OTUs)