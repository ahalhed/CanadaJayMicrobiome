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
# maximum distance
max_dist <- function(dm) {
  # dm is a distance matrix containing the distances between
  # sampling locations
  df1 <- as.data.frame(as.matrix(dm))
  # find max in each column
  summ <-  df1 %>% summarise_if(is.numeric, ~ max(., na.rm=TRUE))
  m <- max(summ)
  return(m)
}
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
met_filter <- function(meta, season, year) {
  # meta is the full metadata table
  met <- rownames_to_column(meta, var = "SampleID")
  df1 <- subset(met, CollectionSeason == season) %>%
    subset(., CollectionYear == year) %>%
    # may need to fiddle with the env data been uses
    # not including JayID b/c it fully constrains models
    select("SampleID", "Sex", "BirthYear", "CollectionSeason", 
           "CollectionYear", "BreedingStatus", "JuvenileStatus") # make sure this is dplyr select
  df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
  # select columns with more than one level
  df4 <- df2[sapply(df2, function(x) length(unique(x))>1)]
  return(df4)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))

# based on the meta function from the microbiome package
print("Extract the metadata")
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta[is.na(gj_meta)] <- "" # putting blanks instead of NA to fix complete case issue

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

# need to fix something with sex in above loop
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
print("Maximum Distance (m) by Collection Season")
lapply(dist_list, max_dist)


## community objects
# subset the samples from the core microbiome
print("Build the community object (OTU table) for season")
commFull <- lapply(XY_list, comm_obj, c=OTUclr)

# unweighted PCNM
print("Unweighted PCNM - for use with all OTU tables")
pcnm_list <- lapply(dist_list, pcnm)
lapply(pcnm_list, function(x) x$vectors)

# put ordisurfs here
# ordisurf(XY_list[["xyFall2017"]], scores(pcnm_list[["xyFall2017"]], choi=5), bubble = 4, col = "black", main = "PCNM 5")
# ordisurf(XY_list[["xySpring2020"]], scores(pcnm_list[["xySpring2020"]], choi=8), bubble = 4, col = "black", main = "PCNM 8")
# cleanup
rm(OTUclr, comm_obj, max_dist, met_filter, XY_list)

# analysis really starts here
print("Acessing PCNM scores")
scores_list <- lapply(pcnm_list, scores)
# analysis for all OTUs
print("Analysis for All OTUs")
print("Variance partitioning")
vp_mod1_list <- mapply(varpart, commFull, scores_list, data=sea_list, 
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)
vp_mod1_list
# plot the partitioning
pdf(file = "CanadaJayMicrobiome/plots/AdditionalFigures/H1VPmod1.pdf")
# make plot
# plotted in numerical order by month
lapply(vp_mod1_list, plot)
dev.off()
# clean up
rm(vp_mod1_list, dist_list)

# test with RDA
print("Testing with RDA (full model)")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, sea_list, SIMPLIFY=FALSE)

abFrac # Full model

lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a])")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commFull, pcnm_list, sea_list, SIMPLIFY=FALSE)
aFrac
# anova
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# forward selection for parsimonious model
print("Forward selection for parsimonious model")

print("Prediction 1A")
# spatial variables
print("Spatial variables")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commFull, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

# summary of selection process
print("Summary of spatial selection process")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection")
lapply(step.space, anova)

# save plot


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
pdf(file = "CanadaJayMicrobiome/plots/H1StepEnv.pdf")
# make plot
lapply(step.env, plot)
dev.off()

#cleanup
# remove objects to be replaced/no longer needed
rm(vdist,pbcd, commFull, sea_list, abFrac, aFrac, abFrac0,
   step.env, pcnm_df, bcFrac, bcFrac0, step.space)


# distance based RDA using aitchison distance matrix

# read in aitchison distance matrix
dmAitchison <- read_qza("aitchison-distance.qza")
# might switch to the phyloseq implementation
gj_cap <- capscale(dmAitchison$data ~ ProportionSpruceOnTerritory + MeanTempC + CollectionSeason + CollectionYear,
                   data = gj_meta, comm = otu_table(gj_ps), na.action = na.exclude)
# look at summaries
summary(gj_cap)
# simple biplot
pdf("CanadaJayMicrobiome/plots/H1envBiplot.pdf")
plot(gj_cap, main = "Aitchison Distance-based RDA")
dev.off()

# PERMANOVA
adonis2(dmAitchison$data ~ TerritoryQuality + MeanTempC + CollectionSeason + CollectionYear + JayID,
        data = gj_meta, na.action = na.exclude)
# clean up
rm(gj_meta, gj_ps, dmAitchison, gj_cap)
