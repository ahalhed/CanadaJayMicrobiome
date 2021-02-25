# for Hypothesis 1 - related to territorial differences

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
library(geosphere) # using to compute greater circle distance
library(ggpubr)
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
comm_obj <- function(n, c) {
  # subset the OTUs (c is OTU table being subset)
  # n is sample data with rownames to be subset
  comm <- c %>%
    subset(., rownames(.) %in% rownames(n)) %>%
    .[ , colSums(.)>0 ]
  return(comm)
}
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

# make long data frame with sample data and distances
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

# ANOSIM repititions
anoRep <- function(samp, season, year, dis) {
  # samp is df with sample data
  # season is the collection season of interest
  # year is the collection year of interest
  # dis is a distance matrix
  print(paste(season, year))
  SampD <- samp[which(samp$CollectionSeason == season & samp$CollectionYear == year),]
  SampN <- SampD %>% rownames
  dmSamp <- dis %>% as.matrix %>%
    .[which(rownames(.) %in% SampN), which(colnames(.) %in% SampN)] %>%
    as.dist
  ano1A <- with(SampD, anosim(dmSamp, Territory))
  summary(ano1A)
  rm(SampD, SampN, dmSamp, ano1A)
}

# get the data
print("Read in the Data for A and B")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "P1AB-filtered-table.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))

# based on the meta function from the microbiome package
print("Extract the metadata")
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta[is.na(gj_meta)] <- "" # putting blanks instead of NA to fix complete case issue

print("Prediction 1A - Territorial distribution")
# boxplots - data prep (to long)
dmAitchison <- read_qza("P1A-aitchison-distance.qza")$data
# need to remove duplicate comparisons
dm_all <- longDM(dmAitchison, "AitchisonDistance", gj_meta)

# remove within territory comparisons
dm_between <- dm_all[-which(dm_all$Territory.x == dm_all$Territory.y),] %>%
  .[which(.$CollectionYear.x == .$CollectionYear.y),] %>%
  .[which(.$CollectionSeason.x == .$CollectionSeason.y),] %>%
  mutate(Territory = "Between", Group = "Between",
         CollectionYear = CollectionYear.x, CollectionSeason = CollectionSeason.x)
# only within territory groups
dm_within <- dm_all[-which(dm_all$Territory.x != dm_all$Territory.y),] %>%
  .[which(.$CollectionYear.x == .$CollectionYear.y),] %>%
  .[which(.$CollectionSeason.x == .$CollectionSeason.y),] %>%
  mutate(Territory = .$Territory.x, Group = "Within",
         CollectionYear = CollectionYear.x, CollectionSeason = CollectionSeason.x)
# put back together
dm_meta <- rbind(dm_between, dm_within) %>%
  select(Group, Territory, CollectionYear, everything())
# separate plot for each season/year
F17 <- dm_meta[which(dm_meta$CollectionSeason == "Fall" & dm_meta$CollectionYear == 2017),] %>%
  ggplot(aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + rremove("xylab")
F18 <- dm_meta[which(dm_meta$CollectionSeason == "Fall" & dm_meta$CollectionYear == 2018),] %>%
  ggplot(aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + rremove("xylab")
F20 <- dm_meta[which(dm_meta$CollectionSeason == "Fall" & dm_meta$CollectionYear == 2020),] %>%
  ggplot(aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + rremove("xylab")
S20 <- dm_meta[which(dm_meta$CollectionSeason == "Spring" & dm_meta$CollectionYear == 2020),] %>%
  ggplot(aes(y = AitchisonDistance, x = Group)) +
  geom_boxplot() + rremove("xylab")

# save figure
pdf("CanadaJayMicrobiome/plots/P1A.pdf", width = 12, height = 6)
fig <- ggarrange(F17, F18, F20, S20, nrow = 1, vjust = 2, font.label = list(size = 10),
                 labels = c("Fall 2017", "Fall 2018", "Fall 2020", "Spring 2020"))
annotate_figure(fig, bottom = text_grob("Territory Group"),
                left = text_grob("Aitchison Distance", rot = 90))
dev.off()

# ANOSIM
with(gj_meta, anosim(dmAitchison, Territory)) %>% summary
anoRep(gj_meta, "Fall", 2017, dmAitchison)
anoRep(gj_meta, "Fall", 2018, dmAitchison)
anoRep(gj_meta, "Spring", 2020, dmAitchison)
anoRep(gj_meta, "Fall", 2020, dmAitchison)

# clean up
rm(dm_all, dm_between, dm_meta, dm_within, dmAitchison, fig,
   F17, F18, F20, S20)

print("Prediction 1B - Spatial distribution")
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
# try using group_split instead of loop
# need to fix something with sex in above loop
# make a list of the season data frames generated from the loop
sea_list <- list(seaFall2017 = seaFall2017, seaFall2018 = seaFall2018,
                 seaFall2020 = seaFall2020, seaSpring2020 = seaSpring2020)
# clean up individual data frames, now that the list is there
rm(SeaSon, YeaR, seaFall2017, seaFall2018, seaFall2020, seaSpring2020,
   # not enough samples in these ones
   seaSpring2017, seaSpring2018, seaWinter2017, seaWinter2018, seaWinter2020)

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
   xySpring2017, xySpring2018, xyWinter2017, xyWinter2018, xyWinter2020)

print("Computing Haversine Distances")
# using Haversine distance to get distance between sampling locations in meters
# reasonable alternative to doing euclidean distances
dist_list <- lapply(XY_list, function(x) distm(x, x, fun = distHaversine))
print("Maximum Distance (m) by Collection Season")
lapply(dist_list, max_dist)

print("Build the community objects (OTU table) for season")
commFull <- lapply(XY_list, comm_obj, c=OTUclr)

print("Unweighted PCNM - for use with all OTU tables")
pcnm_list <- lapply(dist_list, pcnm)
lapply(pcnm_list, function(x) x$vectors)

# cleanup
rm(OTUclr, max_dist, XY_list)

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
pdf(file = "CanadaJayMicrobiome/plots/AdditionalFigures/P1BVPmod1.pdf")
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
pdf(file = "CanadaJayMicrobiome/plots/P1BStepSpace.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Bray-Curtis dissimilarities")
vdist <- lapply(commFull, dist) # euclidean dist on clr = aitchison
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, sea_list, SIMPLIFY=FALSE)
pbcd

# clean up
rm(gj_meta, abFrac, aFrac, bcFrac, bcFrac0, pbcd, pcnm_df, commFull,
   gj_ps, pcnm_list, scores_list, step.space, vdist, sea_list)

print("Prediction 1C - territory quality")
print("Read in the Data")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "P1C-filtered-table.qza",
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
dmAitchison <- read_qza("P1C-aitchison-distance.qza")$data
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
pdf("CanadaJayMicrobiome/plots/P1CterBiplot.pdf", width = 17.5, height = 9)
lapply(cap_list, plot, main = "Aitchison Distance-based RDA")
dev.off()

#clean up
rm(commFull, dmYear, cap_list, gj_meta, gj_ps, OTUclr, sea_list, dmAitchison,
   comm_obj, dmFilter, met_filter)
