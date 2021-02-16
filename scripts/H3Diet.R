# for Hypothesis 2 - feeding differences
# this script is for the diet supplementation plot
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(vegan)
library(car)
library(lubridate)
library(tidyverse)

theme_set(theme_bw())

# generate functions for analysis
cacheGroup <- function(Data, Group, Date1, Date2) {
  # Data is a dataframe with all weather data
  # Group is the label assigned to the output
  # Date1 and Date2 are the start/end dates of the group
  df <- Data %>% mutate(Group = Group) %>% # labelling groups
    subset(`Date/Time` > Date1 & `Date/Time` < Date2) %>%
    # assign a number to each date in order
    .[order(.$`Date/Time`), ] %>% mutate(Number = rep(1:(0.5*nrow(. )), each=2))
  return(df)
}

eventCount <- function(column, station, event){
  # station is the data from the weather station with freeze thaw & season information
  # column is the data frame column containing the sampling dates of interest
  # event is the weather event (freeze-thaw, warm, frozen)
  df1 <- station %>% filter(`Date/Time` %in% seq(column - days(14), column - days(1), by="day"))
  df2 <- df1 %>% select(FreezeThaw) %>% count(vars = FreezeThaw)
  df3 <- df2 %>% filter(vars == event) %>% .$n
  return(df3)
}

print("Prediction 3A + B - Data for freeze thaw")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3AB-filtered-table.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         tree = "trees/rooted_tree.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

# read in weather data
weather <- read_csv("CanadaJayMicrobiome/data/2020.csv") %>% 
  rbind(read_csv("CanadaJayMicrobiome/data/2019.csv")) %>% 
  rbind(read_csv("CanadaJayMicrobiome/data/2018.csv")) %>% 
  rbind(read_csv("CanadaJayMicrobiome/data/2017.csv")) %>% 
  rbind(read_csv("CanadaJayMicrobiome/data/2016.csv"))
# add new freeze-thaw column
weather$FreezeThaw <- ifelse(weather$`Max Temp (°C)` > -1.9 & weather$`Min Temp (°C)` < -1.9, "Freeze-Thaw",
                             ifelse(weather$`Max Temp (°C)` > -1.9 & weather$`Min Temp (°C)` > -1.9, "Warm", "Frozen"))
# add new season type column
weather$SeasonType <- ifelse(as.numeric(weather$Month) > 9, "Caching",
                             ifelse(as.numeric(weather$Month) < 4, "Pre-breeding",
                                    ifelse(as.numeric(weather$Month) < 6, "Breeding",
                                           "Summer")))
# fomatting columns
# caching period = October-December
# pre-breeding period = January-February
longWeather <-  weather %>%
  select(`Min Temp (°C)`, `Max Temp (°C)`,`Date/Time`, FreezeThaw, SeasonType) %>%
  rename('Daily Minimum' = `Min Temp (°C)`, 
         'Daily Maximum' = `Max Temp (°C)`) %>%
  # long formatting (for making plot)
  pivot_longer(-c(`Date/Time`, FreezeThaw, SeasonType), names_to = "Type", values_to = "Temperature")
# select only caching period and pre-breeding period by annual season
# 2016-2017
weather1 <- cacheGroup(longWeather, "2016-2017", "2016-09-01", "2017-08-31")
# 2017-2018
weather2 <- cacheGroup(longWeather, "2017-2018", "2017-09-01", "2018-08-31")
# 2018-2019
weather3 <- cacheGroup(longWeather, "2018-2019", "2018-09-01", "2019-08-31")
# 2019-2020
weather4 <- cacheGroup(longWeather, "2019-2020", "2019-09-01", "2020-08-31")
# 2020 Caching
weather5 <-cacheGroup(longWeather, "Fall 2020", "2020-09-01", "2021-01-01")
# combine groups
weatherCombo <- rbind(weather1, weather2) %>% rbind(., weather3) %>%
  rbind(., weather4) %>% rbind(., weather5)
# clean up
rm(longWeather, weather1, weather2, weather3, weather4, weather5)

print("Prediction 3A - Freeze thaw variation")
# select most relevant sample data
# get only those who did not receive food supplementation
dates <- gj_meta[gj_meta$FoodSupplement == "N",] %>%
  rownames_to_column(var = "sampleID") %>%
  select(JayID, CollectionDate, CollectionDay, sampleID,
         CollectionMonth, CollectionSeason, CollectionYear,
         TerritoryQuality, BreedingStatus) %>%
  # modify the dates to match weather dates
  mutate(CollectionDate = dmy(.$CollectionDate)) %>%
  .[!is.na(.$CollectionDate),] %>% # remove any missing dates
  remove_rownames() %>% column_to_rownames(var = "sampleID")
# merge sample data with weather data
metaWeather <- weatherCombo %>%
  rename("CollectionDate" = "Date/Time") %>%
  select(CollectionDate, SeasonType, Group) %>% 
  unique %>% right_join(dates)

# count the number of freeze-thaw days in the 14 days prior to collection
metaWeather$FreezeThaw <- sapply(metaWeather$CollectionDate, eventCount, 
                                 station = weather, event = 'Freeze-Thaw')
# count the number of warm days in the 14 days prior to collection
metaWeather$Warm <- sapply(metaWeather$CollectionDate, eventCount, 
                           station = weather, event = 'Warm')
metaWeather$Warm[metaWeather$Warm== "integer(0)"] <- 0
# count the number of cold days in the 14 days prior to collection
metaWeather$Frozen <- sapply(metaWeather$CollectionDate, eventCount, 
                             station = weather, event = 'Frozen')
metaWeather$Frozen[metaWeather$Frozen== "integer(0)"] <- 0

# read in the aitchison distance matrix
dmAitchison <- read_qza("P3AB-aitchison-distance.qza")$data %>%
  as.matrix() %>% # remove row without date information
  .[rownames(.) %in% rownames(dates), colnames(.) %in% rownames(dates)]

# PERMANOVA
adonis2(dmAitchison ~ FreezeThaw + SeasonType + Group,
        data = metaWeather, na.action = na.exclude)

# clean up
rm(cacheGroup, eventCount, dates, metaWeather, weatherCombo, weather, dmAitchison)

print("Prediction 3B - Shared OTUs")
#going to put heatmap?

# clean up - removing 3A+B data
rm(gj_meta, gj_ps)

print("Prediction 3C - Food supplementation")
# Read in 2017-2018 distance matrix
dmAitchisonB <- read_qza("H2aitchison-distance.qza")$data 
# combine the aitchison distance data with metadata
dm_meta <- dmAitchisonB %>% as.matrix %>% as.data.frame %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "AitchisonDistance") %>%
  # filter out duplicated comparisons (taling "lower" part of dm df)
  mutate(Sample1 = gsub("G", "", as.character(factor(.$Sample1))) %>% as.numeric(), 
         Sample2 = gsub("G", "", as.character(factor(.$Sample2))) %>% as.numeric()) %>%
  .[as.numeric(.$Sample1) > as.numeric(.$Sample2), ] %>%
  mutate(Sample1 = paste0("G", as.character(Sample1)), # Fixing sample namanes
         Sample2 = paste0("G", as.character(Sample2))) %>% # joing with sample data
  left_join(., rownames_to_column(gj_meta, var = "Sample1")) %>%
  left_join(., rownames_to_column(gj_meta, var = "Sample2"), by = "Sample2") %>%
  # add shared food information
  mutate(host = ifelse(.$JayID.x == .$JayID.y, "Same Bird", "Different Bird"),
         group = ifelse(.$`FoodSupplement.x` == "Y" & .$`FoodSupplement.y` == "Y", "Yes",
                        ifelse(.$`FoodSupplement.x` == "N" & .$`FoodSupplement.y` == "N",
                               "No", "Either"))) %>%
  # select only the variables of interest for the boxplot
  select(1:4, ExtractID.y, host, group)

pdf("CanadaJayMicrobiome/plots/H2BBox.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = group)) +
  geom_boxplot() + labs(x = "Food Supplementation")
dev.off()

dm_meta[dm_meta$group=="No",] %>% .$AitchisonDistance
dm_meta[dm_meta$group=="Yes",] %>% .$AitchisonDistance
# Mann-Whitney test
print("2b Mann-Whitney")
wilcox.test(dm_meta[dm_meta$group=="No",] %>% .$AitchisonDistance,
            dm_meta[dm_meta$group=="Yes",] %>% .$AitchisonDistance)
print("2b t-test")
t.test(dm_meta[dm_meta$group=="No",] %>% .$AitchisonDistance,
       dm_meta[dm_meta$group=="Yes",] %>% .$AitchisonDistance)
# read in the aitchison ordination (probs won't keep this)
ordiAitchison <- read_qza("H2aitchison-ordination.qza")
# combine aitchison vectors with environmental data
aitch <- gj_meta %>% rownames_to_column(var = "SampleID") %>%
  right_join(ordiAitchison$data$Vectors)
# save ordination plot to PDF
pdf("CanadaJayMicrobiome/plots/H2BOrdi.pdf", width = 8)
# ellipse by territory
# shape by food
# too few samples to get much meaningful from this
ggplot(aitch, aes(x=PC1, y=PC2, shape = FoodSupplement, linetype = FoodSupplement)) + 
  geom_point() + # points are samples
  stat_ellipse(type = "norm") + 
  coord_fixed(ylim = c(-0.8, 0.8), xlim = c(-0.8, 0.8)) +
  labs(shape = "Food Supplementation", linetype = "Food Supplementation")
# turn off graphics device
dev.off()
