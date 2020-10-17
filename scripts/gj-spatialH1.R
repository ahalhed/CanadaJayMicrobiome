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
library(tidyverse)

theme_set(theme_bw())

print("Initiate functions for analysis")
# maximum distance
max_dist <- function(dm) {
  df1 <- as.data.frame(as.matrix(dm))
  # message: `summarise_each_()` is deprecated as of dplyr 0.7.0. (use across)
  summ <- summarise_each(df1, ~ max(df1, na.rm=TRUE))
  m <- apply(summ, 1, max)
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

# get the data
print("Read in the Data")
print("Building phyloseq object")
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-singletons-mitochondria-chloroplast.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         # will eventually need to remove 2 once all samples are there
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met2.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.))

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
  for (Year in unique(gj_meta$CollectionYear)) {
    met <- rownames_to_column(gj_meta, var = "SampleID")
    df1 <- subset(met, CollectionSeason == SeaSon) %>%
      subset(CollectionYear == Year) %>%
      select(1, 3:6, 8:10, 13:20, 29:34) # make sure this is dplyr select
    df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
    assign(paste0("sea",SeaSon,Year),df2)
    rm(met, df1, df2)
  }
}
# clean up from loop (don't need blanks in list)
rm(seaBLANK2016, seaBLANK2017, seaBLANK2018, SeaSon, Year)
# make a list of the season data frames generated from the loop
sea_list <- do.call("list",
                   # searching the global environment for the pattern
                   mget(grep("sea", names(.GlobalEnv), value=TRUE)))
# clean up individual data frames, now that the list is there
rm(seaFall2016, seaFall2017, seaFall2018)

print("Accessing the XY metadata by season/year")
# loop to create individual month data frames
for (SeaSon in unique(gj_meta$CollectionSeason)) {
  for (Year in unique(gj_meta$CollectionYear)) {
    met <- rownames_to_column(gj_meta, var = "SampleID") %>% 
      rename("Latitude"="YsampDD", "Longitude"="XsampDD")
    df1 <- subset(met, CollectionSeason == SeaSon) %>%
      subset(CollectionYear == Year) %>%
      select("SampleID", "Latitude", "Longitude") #dplyr select
    df2 <- column_to_rownames(remove_rownames(df1), var = "SampleID")
    assign(paste0("xy",SeaSon,Year),df2)
    rm(met, df1, df2)
  }
}
# we aren't interested in the blanks, so we will remove those (clean up)
rm(xyBLANK2016, xyBLANK2017, xyBLANK2018, SeaSon, Year)
# make a list of the xy data frames generated from the loop
XY_list <- do.call("list",
                   # searching the global environment for the pattern
                   mget(grep("xy", names(.GlobalEnv), value=TRUE)))
# clean up individual data frames, now that the list is there
rm(xyFall2016,xyFall2017,xyFall2018)

print("Computing Euclidean Distances")
dist_list <- lapply(XY_list, dist)
print("Maximum Euclidean Distance by Month")
lapply(dist_list, max_dist)

## community objects
# subset the samples from the core microbiome
print("Build the community object (OTU table) for season/year")
commFull <- lapply(XY_list, comm_obj, c=OTUclr)
commCore <- lapply(XY_list, comm_obj, c=OTU_core)
commRare <- lapply(XY_list, comm_obj, c=OTU_rare)

# Remove objects we're done with
print("Removing pobjects that are no longer needed")
rm(gj_meta, cOTU, OTUclr, gj_ps, OTU_core, OTU_rare)
# leaving functions, just in case

