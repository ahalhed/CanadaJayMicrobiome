# for Hypothesis 1 - related to environmental/spatial differences
# this script is for the PCNM spatial analysis
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
library(geosphere) # using to compute greater circle distance
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
           "CollectionYear", "BreedingStatus", "JuvenileStatus") # make sure this is dplyr select
  df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
  # select columns with more than one level
  df4 <- df2[sapply(df2, function(x) length(unique(x))>1)]
  return(df4)
}
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

print("Accessing the XY metadata by season")
# loop to create individual season/year data frames
for (YeaR in unique(gj_meta$CollectionYear)) {
  for (SeaSon in unique(gj_meta$CollectionSeason)) {
    met <- rownames_to_column(gj_meta, var = "SampleID") %>% 
      rename("Latitude"="LatitudeSamplingDD", "Longitude"="LongitudeSamplingDD")
    df1 <- subset(met, CollectionSeason == SeaSon) %>%
      subset(., CollectionYear == YeaR) %>%
      select("SampleID", "Longitude", "Latitude") #dplyr select
    df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
    assign(paste0("xy",SeaSon, YeaR),df2)
    rm(met, df1, df2)
  }
}
# make a list of the xy data frames generated from the loop
XY_list <- list(xyFall2017 = xyFall2017, xyFall2018 = xyFall2018, 
                xyFall2020 = xyFall2020, xySpring2020 = xySpring2020)

# clean up individual data frames, now that the list is there
rm(SeaSon, YeaR, xyFall2017, xyFall2018, xyFall2020, xySpring2020,
   # not in list b/c not enough samples for this analysis
   xyFall2016, xyFall2019, xySpring2016, xySpring2017, xySpring2018, xySpring2019,
   xyWinter2016, xyWinter2017, xyWinter2018, xyWinter2019, xyWinter2020)

print("Computing Haversine Distances")
# using Haversine distance to get distance between sampling locations in meters
# reasonable alternative to doing euclidean distances
dist_list <- lapply(XY_list, function(x) distm(x, x, fun = distHaversine))

## community object
print("Build the community object (OTU table) for season")
commFull <- lapply(sea_list, comm_obj, c=OTUclr)
# non-breeders only
nb_list <- lapply(sea_list, function(x) x[which(x$BreedingStatus == "Non-breeder"),])
commJ <- lapply(nb_list, comm_obj, c=OTUclr)

# unweighted PCNM
print("Unweighted PCNM - for use with all OTU tables")
pcnm_list <- lapply(dist_list, pcnm)
lapply(pcnm_list, function(x) x$vectors)
# cleanup
rm(OTUclr, comm_obj, met_filter, XY_list)

# analysis really starts here
print("Acessing PCNM scores")
scores_list <- lapply(pcnm_list, scores)
# test with RDA
print("Testing with RDA (full model)")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, sea_list, SIMPLIFY=FALSE)
abFrac # Full model

print("Prediction 2A - Host associated factors")
# create a tiny anonymous function to include formula syntax in call
abFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commFull, sea_list, SIMPLIFY=FALSE) # Reduced model

step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                   abFrac0, abFrac, SIMPLIFY=FALSE)
step.env # an rda model, with the final model predictor variables
# focus on the environment that is the host
print("Summary of environmental selection process")
lapply(step.env, function(x) x$anova)
print("ANOVA on full environmental selection")
lapply(step.env, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/P2AStepEnv.pdf", width = 10)
# make plot
lapply(step.env, plot)
dev.off()

#cleanup
# remove objects to be replaced/no longer needed
rm(abFrac, abFrac0, step.env, dist_list, gj_meta, gj_ps, scores_list, pcnm_list)

# distance based RDA using aitchison distance matrix
# read in aitchison distance matrix
dmAitchison <- read_qza("P2A-aitchison-distance.qza")$data
# filter dm by year/season
dmYear <- lapply(commFull, dmFilter, dmAitchison)

# might switch to the phyloseq implementation
for (i in names(dmYear)) {
  cap <- capscale(dmYear[[i]] ~ AgeAtCollection + Sex + BirthYear + BreedingStatus + JuvenileStatus,
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

# simple biplot
pdf("CanadaJayMicrobiome/plots/P2AenvBiplot.pdf", width = 12)
lapply(cap_list, plot, main = "Aitchison Distance-based RDA")
dev.off()

# PERMANOVA
for (i in names(dmYear)) {
  p <- adonis2(dmYear[[i]] ~ BreedingStatus + AgeAtCollection + JuvenileStatus + BirthYear + Sex,
          data = sea_list[[i]])
  print(i)
  print(p)
  rm(i, p)
}

# clean up
rm(sea_list,dmYear, cap_list, commFull)

# looking just within non-breeders for juvenile status
print("non-breeders")
dmYear <- lapply(commJ, dmFilter, dmAitchison)

# might switch to the phyloseq implementation
for (i in names(dmYear)) {
  cap <- capscale(dmYear[[i]] ~ JuvenileStatus,
                  data = nb_list[[i]], comm = commJ[[i]])
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
summary(cap_list[[1]])
#summary(cap_list[[2]]) # number 2 acting up
summary(cap_list[[3]])
summary(cap_list[[4]])
#lapply(cap_list, plot, main = "Aitchison Distance-based RDA")

print("F17 only has one level")
#adonis2(dmYear[[1]] ~ JuvenileStatus, data = nb_list[[1]])
adonis2(dmYear[[2]] ~ JuvenileStatus, data = nb_list[[2]])
print("F20 only has one level")
#adonis2(dmYear[[3]] ~ JuvenileStatus, data = nb_list[[3]])
adonis2(dmYear[[4]] ~ JuvenileStatus, data = nb_list[[4]])

# clean up
rm(dmAitchison, sea_list, cap_list)