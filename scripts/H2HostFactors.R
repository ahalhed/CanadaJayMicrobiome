# for Hypothesis 2 - related to host differences
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(zCompositions)
# devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(geosphere)
library(vegan)
library(tidyverse)

theme_set(theme_bw())

print("Initiate functions for analysis")
# community object
comm_obj <- function(XY, c) {
  # subset the OTUs (c is OTU table being subset)
  # XY is the location metadata
  comm <- c %>%
    subset(., rownames(.) %in% rownames(XY)) %>%
    .[ , colSums(.)>0 ]
  return(comm)
}
# subset full metadata by season (placed inside for loop)
met_filter <- function(meta, season, year, status) {
  # meta is the full metadata table
  # season is the collection season
  # year is the collection year
  # status is the breeding status of the individuals (optional)
  met <- rownames_to_column(meta, var = "SampleID")
  df1 <- subset(met, CollectionSeason == season) %>%
    subset(., CollectionYear == year) %>%
    # may need to fiddle with the env data been uses
    # not including JayID b/c it fully constrains models
    select("SampleID", "Sex", "AgeAtCollection", "BirthYear", "CollectionSeason", 
           "CollectionYear", "BreedingStatus") # make sure this is dplyr select
  df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
  # select columns with more than one level
  df4 <- df2[sapply(df2, function(x) length(unique(x))>1)]
  return(df4)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))

# based on the meta function from the microbiome package
print("Extract the metadata")
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta[is.na(gj_meta)] <- "" # putting blanks instead of NA to fix complete case issue
# adding more specific age values (ensuring consistent calculation)
gj_meta$AgeAtCollection <- ifelse(gj_meta$CollectionSeason == "Fall",
                                  0.5 + gj_meta$CollectionYear - gj_meta$BirthYear,
                                  ifelse(gj_meta$CollectionSeason == "Winter",
                                         0.75 + gj_meta$CollectionYear - gj_meta$BirthYear,
                                         gj_meta$CollectionYear - gj_meta$BirthYear))
# making birth year categorical
gj_meta$BirthYear <- as.character(gj_meta$BirthYear)

# example in https://github.com/ggloor/CoDaSeq/blob/master/Intro_tiger_ladybug.Rmd
print("CLR transformation")
# rows are OTUs, then transposed to OTUs as column
# impute the OTU table
OTUclr <- cmultRepl(otu_table(gj_ps), label=0, method="CZM") %>% # all OTUs
  codaSeq.clr # compute the CLR values

print("Accessing the metadata by season/year")
# loop to create individual season data frames
for (YeaR in unique(gj_meta$CollectionYear)) {
  for (SeaSon in unique(gj_meta$CollectionSeason)) {
  met <- met_filter(gj_meta, SeaSon, YeaR)
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
   seaFall2016, seaFall2019, seaSpring2016, seaSpring2017, seaSpring2018, seaSpring2019,
   seaWinter2016, seaWinter2017, seaWinter2018, seaWinter2019, seaWinter2020)


## community object
print("Build the community object (OTU table) for season")
commFull <- lapply(sea_list, comm_obj, c=OTUclr)

# cleanup
rm(OTUclr, comm_obj, met_filter)

print("Prediction 2A - Host associated factors")
# create a tiny anonymous function to include formula syntax in call
RDAs <- mapply(function(x,data) rda(x~., data), 
                  commFull, sea_list, SIMPLIFY=FALSE) # Reduced model
# anova
lapply(RDAs, anova)

print("Spring 2020 - individual variables")
br <- rda(commFull[["seaSpring2020"]] ~ BreedingStatus, sea_list[["seaSpring2020"]])
anova(br)

age <- rda(commFull[["seaSpring2020"]] ~ AgeAtCollection, sea_list[["seaSpring2020"]])
anova(age)

#cleanup
rm(gj_meta, gj_ps, RDAs, sea_list, commFull, age, br)
