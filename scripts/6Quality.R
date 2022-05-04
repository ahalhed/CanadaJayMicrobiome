print("Working directory, packages, functions")
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(vegan)
library(zCompositions)
# devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(tidyverse)

# subset full metadata by season (placed inside for loop)
met_filter <- function(meta, season, year, te=FALSE) {
  # meta is the full metadata table
  # season is the sampling season of interest
  # year is the sampling year of interest
  # te indicates whether to select territory data or not
  if(missing(te)){
    met <- rownames_to_column(meta, var = "SampleID")
    df1 <- subset(met, CollectionSeason == season) %>%
      subset(., CollectionYear == year)
    df2 <- df1 %>%
      select("SampleID", "Sex", "BirthYear", "CollectionSeason",
             "CollectionYear", "JuvenileStatus") # make sure this is dplyr select
    df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
    # select columns with more than one level
    df4 <- df3[sapply(df3, function(x) length(unique(x))>1)]
    return(df4)
  } else {
    met <- rownames_to_column(meta, var = "SampleID")
    df1 <- subset(met, CollectionSeason == season) %>%
      subset(., CollectionYear == year)
    df5 <- df1 %>% select("SampleID", "ProportionSpruceOnTerritory", "TerritoryQuality",
                          "CollectionSeason", "CollectionYear") # make sure this is dplyr select
    df6 <- column_to_rownames(remove_rownames(df5), var = "SampleID") %>%
      .[sapply(., function(x) length(unique(x))>1)]
    return(df6)
  }
}
# community object
comm_obj <- function(n, c) {
  # subset the OTUs (c is OTU table being subset)
  # n is sample data with rownames to be subset
  comm <- c %>%
    subset(., rownames(.) %in% rownames(n)) %>%
    .[ , colSums(.)>0 ]
  return(comm)
}
# filtering matrix
dmFilter <- function(data, dm) {
  # dm is the distance matrix
  # data is the community object whose row names in the DM labels
  # filtering based on the row names present in data
  mat <- dm %>% as.matrix %>%
    .[rownames(.) %in% rownames(data),
      colnames(.) %in% rownames(data)]
  d <- as.dist(mat)
  return(d)
}


print("6 - territory quality")
print("Read in the Data")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "6-filtered-table.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))

# based on the meta function from the microbiome package
print("Extract the metadata")
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta[is.na(gj_meta)] <- "" # putting blanks instead of NA to fix complete case issue

print("Accessing the metadata by season/year")
# loop to create individual season data frames
for (YeaR in unique(gj_meta$CollectionYear)) {
  for (SeaSon in unique(gj_meta$CollectionSeason)) {
    met <- met_filter(gj_meta, SeaSon, YeaR, te=TRUE)
    assign(paste0("sea",SeaSon,YeaR),met)
    rm(met)
  }
}
# make a list of the season data frames generated from the loop
sea_list <- list(seaFall2017 = seaFall2017, seaFall2018 = seaFall2018,
                 seaFall2020 = seaFall2020, seaSpring2020 = seaSpring2020)
# clean up individual data frames, now that the list is there
rm(SeaSon, YeaR, seaFall2017, seaFall2018, seaFall2020, seaSpring2020,
   # not enough samples in these ones
   seaSpring2017, seaSpring2018, seaWinter2017, seaWinter2018, seaWinter2020)

print("CLR transformation")
# rows are OTUs, then transposed to OTUs as column
# impute the OTU table
OTUclr <- cmultRepl(otu_table(gj_ps), label=0, method="CZM") %>% # all OTUs
  codaSeq.clr # compute the CLR values
print("Build the community objects (OTU table) for season")
commFull <- lapply(sea_list, comm_obj, c=OTUclr)

# read in aitchison distance matrix
dmAitchison <- read_qza("6-aitchison-distance.qza")$data
# filter distance matrix by year/season
dmYear <- lapply(commFull, dmFilter, dmAitchison)
# non transformed OTUs
# working on looping the capscale below
# might switch to the phyloseq implementation
for (i in names(dmYear)) {
  cap <- capscale(dmYear[[i]] ~ ProportionSpruceOnTerritory + TerritoryQuality,
                  data = sea_list[[i]], comm = commFull[[i]])
  assign(paste0("cap",i),cap)
  rm(i, cap)
}

# make a list of the capscale data generated from the loop
cap_list <- list(capseaFall2017 = capseaFall2017, capseaFall2018 = capseaFall2018, 
                 capseaFall2020 = capseaFall2020, capseaSpring2020 = capseaSpring2020)
# clean up individual data frames, now that the list is there
rm(capseaFall2017, capseaFall2018, capseaFall2020, capseaSpring2020)

# look at summaries
cap_list
lapply(cap_list, summary)
lapply(cap_list, anova, step=200, perm.max=1000)

# simple biplot
pdf("CanadaJayMicrobiome/plots/P5B.pdf", width = 17.5, height = 9)
lapply(cap_list, plot, main = "Aitchison Distance-based RDA")
dev.off()

#clean up
rm(commFull, dmYear, cap_list, gj_meta, gj_ps, OTUclr, sea_list, dmAitchison,
   comm_obj, dmFilter, met_filter)
