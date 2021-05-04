# for Hypothesis 1 - related to seasonal differences

print("Working directory, packages, functions")
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(qiime2R)
library(phyloseq)
library(ALDEx2)
library(tidyverse)

aldAn <- function(met, year, group, lev, v = FALSE) {
  # met is the sample metadata
  # year is the year(s) of interest
  # group is the column to select
  # lev is a level(s) of group to select (optional)
  # v is whether or not to view aldex output
  df1 <- met[which(met$CollectionYear %in% year),]
  if(missing(lev)){
    df <- df1
  } else {
    df <- df1[which(df1[,group] %in% lev),]
  }
  conds <- df %>%
    select(group) %>% t %>% as.data.frame %>%
    .[ , order(names(.))] %>% unlist
  orderPathAbun <- pathAbunI %>% select(names(conds))
  Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)
  # volcano plots
  if(missing(lev)){
    pa <- paste0("CanadaJayMicrobiome/plots/volcanoes/", year, group, ".pdf")
  } else {
    pa <- paste0("CanadaJayMicrobiome/plots/volcanoes/", year, group, lev, ".pdf")
  }
  pdf(pa)
  par(mfrow=c(1,2))
  aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
             ylab="Difference")
  aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
             ylab="Difference")
  dev.off()
  # clean up
  rm(Ald, conds, orderPathAbun)
}


print("prediction 1A - by season?")
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

print("prediction 1B - functional?")
print("read in sample data")
meta <- read_q2metadata("input/jay-met.tsv") %>%
  .[which(.$CollectionSeason == "Fall"),] %>%
  remove_rownames() %>% column_to_rownames(var = "SampleID")
pathAbun <- read_qza("q2-picrust2_output/pathway_abundance.qza")$data
pathAbunI <- pathAbun %>% as.data.frame %>% mutate_all(~as.integer(.))

# year - 2017 vs 2018 (sig)
# can only have 2 levels
aldAn(meta, c(2017, 2018), "CollectionYear")

# year - 2017 vs 2020 (sig)
aldAn(meta, c(2017, 2020), "CollectionYear")

# year - 2018 vs 2020 (sig)
aldAn(meta,  c(2018, 2020), "CollectionYear")

# breeder vs non-breeder  (no sig)
aldAn(meta, 2017, "BreedingStatus", c("Breeder", "Non-breeder"))
aldAn(meta, 2018, "BreedingStatus", c("Breeder", "Non-breeder"))
aldAn(meta, 2020, "BreedingStatus", c("Breeder", "Non-breeder"))

# sex (no sig)
aldAn(meta, 2017, "Sex", c("M", "F"))
aldAn(meta, 2018, "Sex", c("M", "F"))
aldAn(meta, 2020, "Sex", c("M", "F"))

# territory quality - h vs l (no sig)
aldAn(meta, 2017, "TerritoryQuality", c("H", "L"))
aldAn(meta, 2018, "TerritoryQuality", c("H", "L"))
aldAn(meta, 2020, "TerritoryQuality", c("H", "L"))
# territory quality - h vs m (no sig)
aldAn(meta, 2017, "TerritoryQuality", c("M", "H"))
aldAn(meta, 2018, "TerritoryQuality", c("M", "H"))
aldAn(meta, 2020, "TerritoryQuality", c("M", "H"))

# territory quality - l vs m (no sig)
aldAn(meta,  2017, "TerritoryQuality", c("L", "M"))
aldAn(meta,  2018, "TerritoryQuality", c("L", "M"))
aldAn(meta,  2020, "TerritoryQuality", c("L", "M"))