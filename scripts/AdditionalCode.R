# additional code that don't necessarily fit into one hypothesis
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

print("counting")
print("breeding status")
gj_meta %>%
  group_by(BreedingStatus) %>%
  count()
print("breeding status and juvenile status")
gj_meta %>%
  group_by(BreedingStatus, JuvenileStatus) %>%
  count()
print("breeding status, juvenile status, and collection season")
gj_meta %>%
  group_by(BreedingStatus, JuvenileStatus, CollectionSeason) %>%
  count()
print("breeding status, juvenile status, collection season, and collection year")
gj_meta %>%
  group_by(BreedingStatus, JuvenileStatus, CollectionSeason, CollectionYear) %>%
  count()
print("Breeders by territory, year, season")
gj_meta %>% filter(BreedingStatus == "Breeder") %>%
  group_by(Territory, CollectionYear, CollectionSeason) %>%
  count()
print("Non-breeders by territory, year, season")
gj_meta %>% filter(BreedingStatus == "Non-breeder") %>%
  group_by(Territory, CollectionYear, CollectionSeason) %>%
  count()

print("Should we do analyses by year/season?")
# confirmation that year/season was the right way to go
print("Collection Year Ordination")
(orY <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionYear))
anova(orY)
print("Collection Season Ordination")
(orS <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionSeason))
anova(orS)
print("Collection Year*Season Ordination")
(orSY <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionYear * CollectionSeason))
anova(orSY)
# clean up
rm(gj_meta, gj_ps, orS, orSY, orY)

print("Are there any core differences in each season?")
F17 <- read.csv("CanadaJayMicrobiome/data/coreJayF17.csv") %>%
  .[which(.$fill == "Core"),]
F18 <- read.csv("CanadaJayMicrobiome/data/coreJayF18.csv") %>%
  .[which(.$fill == "Core"),]
F20 <- read.csv("CanadaJayMicrobiome/data/coreJayF20.csv") %>%
  .[which(.$fill == "Core"),]
S20 <- read.csv("CanadaJayMicrobiome/data/coreJayS20.csv") %>%
  .[which(.$fill == "Core"),]
# replace missing genera
F17$Genus[is.na(F17$Genus)] <- paste0("Unassigned",F17$Family[is.na(F17$Genus)])
F18$Genus[is.na(F18$Genus)] <- paste0("Unassigned",F18$Family[is.na(F18$Genus)])
S20$Genus[which(S20$Genus == " uncultured")] <- paste("Uncultured",S20$Family[which(S20$Genus == " uncultured")])
S20$Genus[is.na(S20$Genus)] <- paste("Unassigned",S20$Family[is.na(S20$Genus)])
F17$Genus <- str_replace_all(F17$Genus, "^ ", "") %>%
  str_replace_all(., "  ", " ")
F18$Genus <- str_replace_all(F18$Genus, "^ ", "") %>%
  str_replace_all(., "  ", " ")
F20$Genus <- str_replace_all(F20$Genus, "^ ", "") %>%
  str_replace_all(., "  ", " ")
S20$Genus <- str_replace_all(S20$Genus, "^ ", "") %>%
  str_replace_all(., "  ", " ")
# look at the common genera
comG <- function(dat) {
  # data is a dataframe
  df1 <- table(dat$Genus) %>%
    as.data.frame %>%
    rename(Genus = Var1)
  df1$SY <- deparse(substitute(dat))
  df2 <- aggregate(list(MeanOcc = dat$otu_occ),
                   by=list(Genus=dat$Genus), FUN=mean)
  df3 <- aggregate(list(MeanRel = dat$otu_rel),
                   by=list(Genus=dat$Genus), FUN=mean)
  df4 <- left_join(df1, df2) %>% left_join(df3)
  return(df4)
}

# combine data frames
all <- rbind(comG(F17), comG(F18)) %>%
  rbind(., comG(F20)) %>% rbind(., comG(S20))
# put on figure
f <- ggplot(all, aes(y = MeanOcc, x = MeanRel, shape = SY,
                size = Freq, colour = Genus)) +
  scale_color_viridis_d() + scale_x_log10() +
  geom_point() + guides(color = FALSE) +
  labs(shape = "Season/Year", size = "Number of OTUs",
       x = "Mean Relative Abundance", y = "Mean Occupancy")
# anova
a <- rbind(F17, F18) %>% rbind(., F20) %>% rbind(., S20) %>%
  aov(otu_occ ~ Genus, data = .)
summary(a)
aov(otu_occ ~ SY, data = all) %>% summary
aov(otu_occ ~ Genus, data = all) %>% summary
