print("Set up - working directory, packages, functions, data")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(boot)
library(qiime2R)
library(vegan)
library(phyloseq)
library(ggpubr)
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
# for bootstrapping
fc <- function(d, i){
  # d is the data frame
  # i is the iteration number
  d2 <- d[i,]
  return(mean(d2$diff))
}

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "4-filtered-table.qza",
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
dmAitchison <- read_qza("4A-aitchison-distance.qza")$data

# need to remove duplicate comparisons
dm_all <- longDM(dmAitchison, "AitchisonDistance", gj_meta)

print("Own parent vs other offspring (4A)")
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
dm_meta <- rbind(dm_between, dm_within)
plot4A <- aggregate(dm_meta$AitchisonDistance,
          list(interaction(dm_meta$Group, dm_meta$JayID.x)), mean) %>%
  mutate(AitchisonDistance = x, Group = word(.$'Group.1', 1, sep=fixed('.')),
         NonBreeder = word(.$'Group.1', 2, sep=fixed('.'))) %>%
  select(AitchisonDistance, Group, NonBreeder)
# make figure
figA <- ggplot(plot4A, aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() +
  geom_line(alpha = 0.3, aes(group = NonBreeder)) +
  labs(x = "Territory of Breeder", y = "Mean Aitchison Distance")

# paired t-test
# step 0 - check assumptions
# step 1 - set null and alternative hypotheses
print("H0: d = d0")
print("HA: d < 0")
# step 2 - alpha = 0.05
# step 3 - test statistic
test4A <- plot4A %>%
  pivot_wider(-c(Group), values_from = AitchisonDistance, names_from = Group) %>%
  na.omit() %>% mutate(diff = `Same Territory` - `Different Territory` )
summary(test4A)
SD <- sd(test4A$diff)
samp <- nrow(test4A)
Tstat <- (-0.1062)/(SD/sqrt(samp))
# step 5 - p value
1-pt(Tstat, df=samp)
# step 6 is p<a? no
# step 7 check
t.test(test4A$`Same Territory`, test4A$`Different Territory`,
       paired = T, alternative = "less")

# bootstrapping
boot(test4A, fc, R = 999)

# clean up
rm(dm_between, dm_within, dm_meta, dm_all, dmAitchison,
   plot4A, samp, SD, Tstat, test4A, fc)

print("Shared Microbiota (4B)")
# looking at the number of OTUs (not reads) that occur per sample
# could do reads using sample_sums
OTUs <- otu_table(gj_ps) %>% t %>% as.data.frame()
OTUs[which(OTUs>0, arr.ind=TRUE)] <- 1

# otu and sample data
OTUs <- otu_table(gj_ps) %>% as.matrix %>% as.data.frame()
br <- gj_meta %>% rownames_to_column(var = "b") %>%
  filter(BreedingStatus == "Breeder") %>%
  select(b, JuvenileStatus, Territory)
nb <- gj_meta %>% rownames_to_column(var = "nb") %>%
  filter(BreedingStatus == "Non-breeder") %>%
  filter(JuvenileStatus != "DominantJuvenile") %>%
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
plot4B <- OTUsamples %>% rownames_to_column(var = "pair") %>%
  select(Territory, pair, Bonly, NBonly, shared, `JuvenileStatus.nb`) %>%
  pivot_longer(-c(Territory, pair, `JuvenileStatus.nb`), names_to = "Sharing", values_to = "NumberOfOTUs") %>%
  filter(Territory == "Within") %>%
  mutate(NumberOfOTUs = as.numeric(NumberOfOTUs))

# make percentage, average, line graph
# make plot
figB <- ggplot(plot4B, aes(y = as.numeric(NumberOfOTUs), x = Sharing)) +
  geom_boxplot() +
  geom_line(stat = "smooth", method=loess, alpha = 0.25,
            se=FALSE, color = "black", aes(group = pair)) +
  scale_x_discrete(labels = c("Breeder Only", "Non-breeder Only", "Both")) +
  labs(x = "Sample(s)", y = "Number of OTUs Present")

print("Two sample t-test")
# get data
br <- plot4B[plot4B$Sharing=="Bonly",]
nb <- plot4B[plot4B$Sharing=="NBonly",]
sh <- plot4B[plot4B$Sharing=="shared",]
# null and alternative hypotheses
print("H0 - equal mean OTUs (u1-u2=0)")
print("HA1 - more shared OTUs than unique non-breeder OTUs (u1-u2 > u0)")
print("HA2 - more unique breeder OTUs than unique non-breeder OTUs (u1-u3 > u0)")
# significance (going to 0.05)
# group sd & mean
aggregate(plot4B$NumberOfOTUs, by = list(plot4B$Sharing), FUN=sd)
aggregate(plot4B$NumberOfOTUs, by = list(plot4B$Sharing), FUN=mean)
# two sample t-test
print("Shared & Non-breeder")
t.test(sh$NumberOfOTUs, nb$NumberOfOTUs, alternative = "greater", var.equal = T)
print("Breeder & Non-breeder")
t.test(br$NumberOfOTUs, nb$NumberOfOTUs, alternative = "greater", var.equal = T)

# section clean up
rm(br, nb, sh, otu_br, otu_nb, OTUs, plot4B, shareOTUs, pairedNames, OTUsamples,
   gj_meta, gj_ps, longDM, share)

# arrange figures
pdf("CanadaJayMicrobiome/plots/H4.pdf", width = 10)
ggarrange(figA, figB, labels = c("A", "B"))
dev.off()

# top outlier is a dominant juvenile
print("Dominant Juvenile")
plot4B %>% filter(JuvenileStatus.nb == "DominantJuvenile") %>%
  group_by(Sharing) %>% summarize(count = nrow(.)/3,
                                  mean = mean(NumberOfOTUs),
                                  max = max(NumberOfOTUs),
                                  min = min(NumberOfOTUs))
print("Other Juveniles")
plot4B %>% filter(JuvenileStatus.nb != "DominantJuvenile") %>%
  group_by(Sharing) %>% summarize(count = nrow(.)/3,
                                  mean = mean(NumberOfOTUs),
                                  max = max(NumberOfOTUs),
                                  min = min(NumberOfOTUs))
