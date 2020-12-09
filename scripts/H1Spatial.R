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
met_filter <- function(meta, season) {
  # meta is the full metadata table
  met <- rownames_to_column(meta, var = "SampleID")
  df1 <- subset(met, CollectionSeason == season) %>%
    # may need to fiddle with the env data been uses
    select(1, 4, 5, 12, 16) # make sure this is dplyr select
  df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
  # replace NAs with blanks (so as to retain the columns)
  df3 <- df2 %>% mutate_all(~replace(., is.na(.), ''))
  # select columns with more than one level
  df4 <- df3[sapply(df3, function(x) length(unique(x))>1)]
  return(df4)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-singletons-mitochondria-chloroplast.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))

# based on the meta function from the microbiome package
print("Extract the metadata")
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

# example in https://github.com/ggloor/CoDaSeq/blob/master/Intro_tiger_ladybug.Rmd
print("CLR transformation")
# rows are OTUs, then transposed to OTUs as column
# impute the OTU table
OTUclr <- cmultRepl(otu_table(gj_ps), label=0, method="CZM") %>% # all OTUs
  codaSeq.clr # compute the CLR values

## Core and rare divide
print("Separating core microbiome")
# extract the core identified OTUs from the occ-abun results (coreJay)
cOTU <- read.csv("CanadaJayMicrobiome/data/coreJay.csv") %>% 
  .[which(.$fill == "core"),]
# make the new data frames
print("Subset the OTU table to find core and rare OTUs")
OTU_core <- OTUclr[, cOTU$Feature.ID]
# make sure this is dplyr select
OTU_rare <- select(as.data.frame(OTUclr), -one_of(cOTU$Feature.ID))
print("Accessing the metadata by season/year")
# loop to create individual season data frames
for (SeaSon in unique(gj_meta$CollectionSeason)) {
  met <- met_filter(gj_meta, SeaSon)
  assign(paste0("sea",SeaSon),met)
  rm(met)
}

# clean up from loop (don't need blanks in list)
rm(seaBLANK, SeaSon)
# make a list of the season data frames generated from the loop
sea_list <- do.call("list",
                    # searching the global environment for the pattern
                    mget(grep("sea", names(.GlobalEnv), value=TRUE)))
# clean up individual data frames, now that the list is there
rm(seaFall)
print("Accessing the XY metadata by season")
# loop to create individual season/year data frames
for (SeaSon in unique(gj_meta$CollectionSeason)) {
  met <- rownames_to_column(gj_meta, var = "SampleID") %>% 
    rename("Latitude"="YsampDD", "Longitude"="XsampDD")
  df1 <- subset(met, CollectionSeason == SeaSon) %>%
    select("SampleID", "Longitude", "Latitude") #dplyr select
  df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
  assign(paste0("xy",SeaSon),df2)
  rm(met, df1, df2)
}
# we aren't interested in the blanks, so we will remove that (clean up)
rm(xyBLANK, SeaSon)
# make a list of the xy data frames generated from the loop
XY_list <- do.call("list",
                   # searching the global environment for the pattern
                   mget(grep("xy", names(.GlobalEnv), value=TRUE)))
# clean up individual data frames, now that the list is there
rm(xyFall)

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
commCore <- lapply(XY_list, comm_obj, c=OTU_core)
commRare <- lapply(XY_list, comm_obj, c=OTU_rare)

# Remove objects we're done with
print("Removing pobjects that are no longer needed")
rm(gj_meta, cOTU, OTUclr, gj_ps, OTU_core, OTU_rare)
# leaving functions, just in case

## Analysis time!
# unweighted PCNM
print("Unweighted PCNM - for use with all OTU tables")
pcnm_list <- lapply(dist_list, pcnm)
# need to sort out how to generalize printing  label with the vectors
for (month in pcnm_list) {
  print(month$vectors)
}
print("Acessing PCNM scores")
scores_list <- lapply(pcnm_list, scores)

# core OTUs
print("Analysis for Core OTUs")
print("Variance partitioning - Core OTUs")
# for some reason 2016 is the second option in sea_list, so dropping that
vp_mod1_list <- mapply(varpart, commCore, scores_list, data=sea_list,
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)

# plot the partitioning
pdf(file = "CanadaJayMicrobiome/plots/core_vp_mod1.pdf")
# make plot
# plotted in numerical order by season
lapply(vp_mod1_list, plot)
dev.off()
#remove vp object, to repeat with new OTU table
rm(vp_mod1_list)

# test with RDA
print("Testing with RDA (full model) - core OTUS")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commCore, sea_list, SIMPLIFY=FALSE)
abFrac # Full model
lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a]) - core OTUS")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commCore, pcnm_list, sea_list, SIMPLIFY=FALSE)
aFrac
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# forward selection for parsimonious model
print("Forward selection for parsimonious model - core OTUs")
# env variables
print("Environmental variables - core OTUs")
# create a tiny anonymous function to include formula syntax in call
abFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commCore, sea_list, SIMPLIFY=FALSE) # Reduced model
step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                   abFrac0, abFrac, SIMPLIFY=FALSE)
step.env # an rda model, with the final model predictor variables

print("Summary of environmental selection process - core OTUs")
lapply(step.env, function(x) x$anova)
print("ANOVA on full environmental selection - core OTUs")

lapply(step.env, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/core_step_env.pdf")
# make plot
lapply(step.env, plot)
dev.off()

# spatial variables
print("Spatial variables - core OTU")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commCore, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commCore, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

print("Summary of spatial selection process - core OTU")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection - core OTU")
lapply(step.space, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/core_step_space.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Aitchison dissimilarities - core OTUs")
vdist <- lapply(commCore, dist) # euclidean on clr = aitchison
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, sea_list, SIMPLIFY=FALSE)
pbcd

#cleanup
# remove objects to be replaced/no longer needed
rm(vdist,pbcd, commCore)
rm(abFrac, aFrac,abFrac0, step.env, pcnm_df, bcFrac, bcFrac0, step.space)


# Rare OTUs
print("Analysis for Rare OTUs")
print("Variance partitioning - Rare OTUs")
vp_mod1_list <- mapply(varpart, commRare, scores_list, data=sea_list, 
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)
vp_mod1_list
# plot the partitioning
pdf(file = "CanadaJayMicrobiome/plots/rare_vp_mod1.pdf")
# make plot
# plotted in numerical order by month
lapply(vp_mod1_list, plot)
dev.off()
#remove vp object, to repeat with new OTU table
rm(vp_mod1_list)

# test with RDA
print("Testing with RDA (full model) - rare OTUS")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commRare, sea_list, SIMPLIFY=FALSE)
abFrac # Full model
# anova
lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a]) - rare OTUS")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commRare, pcnm_list, sea_list, SIMPLIFY=FALSE)
aFrac
# anova
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# forward selection for parsimonious model
print("Forward selection for parsimonious model - rare OTUs")
# env variables
print("Environmental variables - rare OTUs")
# create a tiny anonymous function to include formula syntax in call
abFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commRare, sea_list, SIMPLIFY=FALSE) # Reduced model

step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                   abFrac0, abFrac, SIMPLIFY=FALSE)
step.env # an rda model, with the final model predictor variables

print("Summary of environmental selection process - rare OTUs")
lapply(step.env, function(x) x$anova)
print("ANOVA on full environmental selection - rare OTUs")
lapply(step.env, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/rare_step_env.pdf")
# make plot
lapply(step.env, plot)
dev.off()

# spatial variables
print("Spatial variables - rare OTU")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commRare, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commRare, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

# summary of selection process
print("Summary of spatial selection process - rare OTU")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection - rare OTU")
lapply(step.space, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/rare_step_space.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Bray-Curtis dissimilarities - rare OTUs")
vdist <- lapply(commRare, dist) # eclidean dist on clr = aitchison
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, sea_list, SIMPLIFY=FALSE)
pbcd

#cleanup
# remove objects to be replaced
rm(vdist, pbcd, commRare)
rm(abFrac, aFrac,abFrac0, step.env, pcnm_df, bcFrac, bcFrac0, step.space)

# analysis for all OTUs
print("Analysis for All OTUs")
print("Variance partitioning - All OTUs")
vp_mod1_list <- mapply(varpart, commFull, scores_list, data=sea_list, 
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)
vp_mod1_list
# plot the partitioning
pdf(file = "CanadaJayMicrobiome/plots/vp_mod1M.pdf")
# make plot
# plotted in numerical order by month
lapply(vp_mod1_list, plot)
dev.off()
#remove vp object, done with it now
rm(vp_mod1_list)

# test with RDA
print("Testing with RDA (full model) - all OTUS")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, sea_list, SIMPLIFY=FALSE)

abFrac # Full model

lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a]) - all OTUS")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commFull, pcnm_list, sea_list, SIMPLIFY=FALSE)
aFrac
# anova
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# forward selection for parsimonious model
print("Forward selection for parsimonious model - all OTUs")

# env variables
print("Environmental variables - all OTUs")
# create a tiny anonymous function to include formula syntax in call
abFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commFull, sea_list, SIMPLIFY=FALSE) # Reduced model

step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                   abFrac0, abFrac, SIMPLIFY=FALSE)

step.env # an rda model, with the final model predictor variables

print("Summary of environmental selection process - all OTUs")
lapply(step.env, function(x) x$anova)
print("ANOVA on full environmental selection - all OTUs")
lapply(step.env, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/step_envM.pdf")
# make plot
lapply(step.env, plot)
dev.off()

# spatial variables
print("Spatial variables - all OTU")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commFull, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

# summary of selection process
print("Summary of spatial selection process - all OTU")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection - all OTU")
lapply(step.space, anova)

# save plot
pdf(file = "CanadaJayMicrobiome/plots/step_spaceM.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Bray-Curtis dissimilarities - all OTUs")
vdist <- lapply(commFull, dist) # euclidean dist on clr = aitchison
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, sea_list, SIMPLIFY=FALSE)
pbcd

#cleanup
# remove objects to be replaced/no longer needed
rm(vdist,pbcd, commFull)
rm(abFrac, aFrac,abFrac0, step.env, pcnm_df, bcFrac, bcFrac0, step.space)

# I have removed the variation decomposition with parsimonious variables, 
# since it was frequently failing and would likely cause issues.
