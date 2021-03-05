print("Set up - working directory, packages, functions, data")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(qiime2R)
library(vegan)
library(phyloseq)
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
share <- function(otu, b, nb) {
  # otu is the otu table with all samples where columns are OTUs
  # b is the sample ID for the breeder
  # nb is the sample ID for the non-breeder
  phy1 <- otu[which(rownames(otu) %in% c(b, nb)), ]
  phy2 <- phy1[,colSums(phy1 >= 1) >=1]
  # breeder only
  breeder <- otu[which(rownames(otu) %in% b), ]
  breeder2 <- breeder[,colSums(breeder >= 1) >=1]
  # non-breeder only
  nbreed <- otu[which(rownames(otu) %in% nb), ]
  nbreed2 <- nbreed[,colSums(nbreed>= 1) >=1]
  # get the number of OTUs
  shared <- ncol(phy2) %>% as.numeric
  Bonly <- ncol(breeder2) %>% as.numeric
  NBonly <- ncol(nbreed2) %>% as.numeric
  df <- data.frame(b, nb, Bonly, shared, NBonly)
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
# clr OTUs - may not need unless calculating distances here
# OTUclr <- zCompositions::cmultRepl(otu_table(gj_ps), label=0, method="CZM") %>% CoDaSeq::codaSeq.clr(.)

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
  geom_boxplot() + labs(x = "Juvenile Status of Non-Breeder", y = "Aitchison Distance")
dev.off()

# ANOSIM - need to arrange this so the distances are only between breeders and non-breeders
gj_meta$Group <- ifelse(gj_meta$BreedingStatus == "Breeder", "Breeder",
       ifelse(gj_meta$JuvenileStatus == "DominantJuvenile",
              "DominantJuvenile", "Nestling"))

# need to filter out breeder-breeder and non-breeder-non-breeder comparisons (?)
# should also add Territory,
with(gj_meta, anosim(dmAitchison, interaction(Group,Group))) %>% summary

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
  geom_boxplot() +
  geom_dotplot(binaxis = "y", binwidth = 0.1,
               stackdir = "center", fill = NA) +
  labs(x = "Non-Breeder Location", y = "Aitchison Distance")
dev.off()

# ANOSIM - need to arrange this so the distances are only between breeders and non-breeders
# not sure interaction(Territory, BreedingStatus, BreedingStatus) is doing the right thing
with(gj_meta, anosim(dmAitchison, interaction(Territory, BreedingStatus, BreedingStatus))) %>% summary

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
  geom_violin() +
  geom_dotplot(binaxis = "y", binwidth = 10,
               stackdir = "center", fill = NA) +
  labs(y = "Number of OTUs", x = "Breeding Status")
dev.off()

print("Two sample t-test")
# get data
br <- plot4C[plot4C$BreedingStatus=="Breeder",] %>% select(NumberOfOTUs, Territory)
nb <- plot4C[plot4C$BreedingStatus=="Non-breeder",] %>% select(NumberOfOTUs, Territory)
# null and alternative hypotheses
print("H0 - equal mean OTUs (u1-u2=0)")
print("HA - breeder mean OTUs greater than non-breeder mean OTUs (u1-u2>0)")
# significance (going to 0.05)
# group sd & mean
aggregate(plot4C$NumberOfOTUs, by = list(plot4C$BreedingStatus), FUN=sd)
aggregate(plot4C$NumberOfOTUs, by = list(plot4C$BreedingStatus), FUN=mean)
# two sample t-test
t.test(br$NumberOfOTUs, nb$NumberOfOTUs, alternative = "greater", var.equal = T)

#clean up
rm(plot4C, nb, br, OTUs)

print("Shared Microbiota (4D)")
# otu and sample data
OTUs <- otu_table(gj_ps) %>% as.matrix %>% as.data.frame()
br <- gj_meta %>% rownames_to_column(var = "b") %>%
  filter(BreedingStatus == "Breeder") %>%
  select(b, JuvenileStatus, Territory)
nb <- gj_meta %>% rownames_to_column(var = "nb") %>%
  filter(BreedingStatus == "Non-breeder") %>%
  select(nb, JuvenileStatus, Territory)
# get grouped samples
otu_br <- OTUs[rownames(OTUs) %in% br$b,]
otu_nb <- OTUs[rownames(OTUs) %in% nb$nb,]

# get shared OTUs
pairedNames <- as.list(outer(rownames(otu_br), rownames(otu_nb), paste))
shareOTUs <- mapply(share, word(pairedNames), word(pairedNames,2),
                    MoreArgs = list(otu = OTUs)) %>% t %>%
  as.data.frame %>% remove_rownames() %>%
  mutate(b = as.character(b), nb = as.character(nb))

# combine shareOTUs with sample data
OTUsamples <- shareOTUs %>% left_join(br) %>%
  left_join(nb, by = "nb", suffix = c(".b", ".nb")) %>%
  mutate(Territory = ifelse(Territory.b == Territory.nb, "Within", "Between"))

# for plotting
plot4D <- OTUsamples %>%
  select(Territory, Bonly, NBonly, shared) %>%
  pivot_longer(-Territory, names_to = "Sharing", values_to = "NumberOfOTUs") %>%
  filter(Territory == "Within") %>%
  mutate(NumberOfOTUs = as.numeric(NumberOfOTUs))

# make plot
pdf("CanadaJayMicrobiome/plots/P4D.pdf")
ggplot(plot4D, aes(y = as.numeric(NumberOfOTUs), x = Sharing)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("Breeder Only", "Non-breeder Only", "Both")) +
  labs(x = "Sample(s)", y = "Number of OTUs Present") +
  ggtitle("Within Territory OTUs Shared Between Individuals")
dev.off()

print("Two sample t-test")
# get data
br <- plot4D[plot4D$Sharing=="Bonly",]
nb <- plot4D[plot4D$Sharing=="NBonly",]
sh <- plot4D[plot4D$Sharing=="shared",]
# null and alternative hypotheses
print("H0 - equal mean OTUs (u1-u2=0)")
print("HA1 - more shared OTUs than unique non-breeder OTUs (u1-u2 > u0)")
print("HA2 - more unique breeder OTUs than unique non-breeder OTUs (u1-u3 > u0)")
# significance (going to 0.05)
# group sd & mean
aggregate(plot4D$NumberOfOTUs, by = list(plot4D$Sharing), FUN=sd)
aggregate(plot4D$NumberOfOTUs, by = list(plot4D$Sharing), FUN=mean)
# two sample t-test
print("Shared & Non-breeder")
t.test(sh$NumberOfOTUs, nb$NumberOfOTUs, alternative = "greater", var.equal = T)
print("Breeder & Non-breeder")
t.test(br$NumberOfOTUs, nb$NumberOfOTUs, alternative = "greater", var.equal = T)

# section clean up
rm(br, nb, sh, otu_br, otu_nb, OTUs, plot4D, shareOTUs, pairedNames, OTUsamples,
   gj_meta, gj_ps, longDM, share)
