# for Hypothesis 3 - feeding differences
# this script is for the diet supplementation plot
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ALDEx2)
library(tidyverse)
# read in data
pathAbun <- read_qza("q2-picrust2_output/pathway_abundance.qza")$data
pathAbunI <- pathAbun %>% as.data.frame %>% mutate_all(~as.integer(.))
meta <- read_q2metadata("input/jay-met.tsv") %>%
  .[which(.$JayID != "BLANK"),] %>%
  remove_rownames() %>% column_to_rownames(var = "SampleID")

# season - spring vs fall 2020 (sig)
# within a year accounts between year variation
conds <- meta[which(meta$CollectionSeason %in% c("Fall", "Spring") &
                      meta$CollectionYear == 2020),] %>%
  select(CollectionSeason) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
                   include.sample.summary=FALSE, denom="all", verbose=FALSE)
# volcano plots
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")
# clean up
rm(Ald, conds, orderPathAbun)

# year - 2017 vs 2018 (sig)
# can only have 2 levels (2017 vs 2018 as example)
conds <- meta[which(meta$CollectionYear %in% c(2017, 2018) &
                      meta$CollectionSeason == "Fall"),] %>%
  select(CollectionYear) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
             include.sample.summary=FALSE, denom="all", verbose=FALSE)
# volcano plots
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")
# clean up
rm(Ald, conds, orderPathAbun)

# sex (no sig)
conds <- meta[which(meta$BreedingStatus == "Breeder"),] %>%
  select(Sex) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
             include.sample.summary=FALSE, denom="all", verbose=FALSE)
# volcano plots
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")
# clean up
rm(Ald, conds, orderPathAbun)

# breeder vs non-breeder  (no sig)
conds <- meta[which(meta$BreedingStatus != "Unknown"),] %>%
  select(BreedingStatus) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
             include.sample.summary=FALSE, denom="all", verbose=FALSE)
# volcano plots
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")
# clean up
rm(Ald, conds, orderPathAbun)

# territory quality (no sig)
conds <- meta[which(meta$BreedingStatus == "Breeder" & meta$TerritoryQuality %in% c("L", "H")),] %>% # can only have 2 levels
  select(TerritoryQuality) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
             include.sample.summary=FALSE, denom="all", verbose=FALSE)
# volcano plots
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")
# clean up
rm(Ald, conds, orderPathAbun)
