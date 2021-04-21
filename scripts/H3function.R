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

aldAn <- function(met, season, year, group, lev, v = FALSE) {
  # met is the sample metadata
  # season is the season(s) of interest
  # year is the year(s) of interest
  # group is the column to select
  # lev is a level(s) of group to select (optional)
  # v is whether or not to view aldex output
  df1 <- met[which(met$CollectionSeason %in% season &
                    met$CollectionYear %in% year),]
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
  # do we view the data?
  if(!missing(v)){
    View(Ald)
  }
  # volcano plots
  par(mfrow=c(1,2))
  aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
             ylab="Difference")
  aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
             ylab="Difference")
  # clean up
  rm(Ald, conds, orderPathAbun)
}

# season - spring vs fall 2020 (sig)
# within a year accounts between year variation
aldAn(meta, c("Fall", "Spring"), 2020, "CollectionSeason")

# year - 2017 vs 2018 (sig)
# can only have 2 levels
aldAn(meta, "Fall", c(2017, 2018), "CollectionYear")

# year - 2017 vs 2020 (sig)
aldAn(meta, "Fall", c(2017, 2020), "CollectionYear")

# year - 2018 vs 2020 (sig)
aldAn(meta, "Fall", c(2018, 2020), "CollectionYear")

# breeder vs non-breeder  (no sig)
aldAn(meta, "Fall", 2017, "BreedingStatus", c("Breeder", "Non-breeder"))
aldAn(meta, "Fall", 2018, "BreedingStatus", c("Breeder", "Non-breeder"))
aldAn(meta, "Fall", 2020, "BreedingStatus", c("Breeder", "Non-breeder"))
aldAn(meta, "Spring", 2020, "BreedingStatus", c("Breeder", "Non-breeder"))

# sex (no sig)
aldAn(meta, "Fall", 2017, "Sex", c("M", "F"))
aldAn(meta, "Fall", 2018, "Sex", c("M", "F"))
aldAn(meta, "Fall", 2020, "Sex", c("M", "F"))
aldAn(meta, "Spring", 2020, "Sex", c("M", "F"))

# territory quality - h vs l (no sig)
aldAn(meta, "Fall", 2017, "TerritoryQuality", c("L", "H"))
aldAn(meta, "Fall", 2018, "TerritoryQuality", c("L", "H"))
aldAn(meta, "Fall", 2020, "TerritoryQuality", c("L", "H"))
aldAn(meta, "Spring", 2020, "TerritoryQuality", c("L", "H"))

# territory quality - h vs m (no sig in fall, sig in spring)
aldAn(meta, "Fall", 2017, "TerritoryQuality", c("M", "H"))
aldAn(meta, "Fall", 2018, "TerritoryQuality", c("M", "H"))
aldAn(meta, "Fall", 2020, "TerritoryQuality", c("M", "H"))
aldAn(meta, "Spring", 2020, "TerritoryQuality", c("M", "H"))

# territory quality - l vs m (no sig)
aldAn(meta, "Fall", 2017, "TerritoryQuality", c("M", "L"))
aldAn(meta, "Fall", 2018, "TerritoryQuality", c("M", "L"))
aldAn(meta, "Fall", 2020, "TerritoryQuality", c("M", "L"))
aldAn(meta, "Spring", 2020, "TerritoryQuality", c("M", "L"))

# working on this
# multiple levels
mm <- model.matrix(~ CollectionSeason*CollectionYear, meta)
orderPathAbun <- pathAbunI %>% select(row.names(mm))
x <- aldex.clr(orderPathAbun, mm, mc.samples=128, denom="all")
glm.test <- aldex.glm(x)
summary(glm.test)
