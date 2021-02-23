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

print("Prediction 3A + B - Weather data for freeze thaw")
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

# join with sample data
samp <- read_q2metadata("input/jay-met.tsv") %>%
  rename("sampleID"="SampleID")
# select most relevant sample data
# double check how this is running
# get only those who did not receive food supplementation
dates <- samp[samp$FoodSupplement == "N",] %>%
  select(JayID, CollectionDate, CollectionDay, sampleID,
         CollectionMonth, CollectionSeason, CollectionYear,
         TerritoryQuality, BreedingStatus) %>%
  # modify the dates to match weather dates
  mutate(CollectionDate = dmy(.$CollectionDate)) %>%
  .[!is.na(.$CollectionDate),] # remove any missing dates
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
#clean up
rm(cacheGroup, eventCount, weatherCombo, weather)

print("Prediction 3A - Freeze thaw variation")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3A-filtered-table.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# double check the combination here
# read in the aitchison distance matrix
dmAitchison <- read_qza("P3A-aitchison-distance.qza")$data %>%
  as.matrix() %>% # remove row without date information
  .[rownames(.) %in% rownames(gj_meta), colnames(.) %in% rownames(gj_meta)]

# get data frame with sames collected on the same date
plot3A <- longDM(dmAitchison, "AitchisonDistance", gj_meta) %>%
  .[which(.$CollectionDate.x == .$CollectionDate.y),] %>%
  mutate(CollectionDate = dmy(CollectionDate.x),
         Year = CollectionYear.x) %>%
  select(CollectionDate, AitchisonDistance, Year, Sample1, Sample2) %>%
  left_join(metaWeather %>% select(CollectionDate, FreezeThaw, CollectionSeason),
            by="CollectionDate") %>%
  unique() %>% # remove duplicate rows
  # add in the breeding status information
  left_join(metaWeather %>% select(sampleID, BreedingStatus),
          by = c("Sample1" = "sampleID")) %>%
  left_join(metaWeather %>% select(sampleID, BreedingStatus),
            by = c("Sample2" = "sampleID"), suffix = c("1", "2"))
plot3A$Group <- ifelse(plot3A$BreedingStatus1 == plot3A$BreedingStatus2 &
                         plot3A$BreedingStatus1 == "Breeder",
                       "Within Breeders",
                       ifelse(plot3A$BreedingStatus1 == plot3A$BreedingStatus2 &
                           plot3A$BreedingStatus1 == "Non-breeder",
                         "Within Non-Breeders", "Between Breeder & Non-Breeder"))


# make figure
pdf("CanadaJayMicrobiome/plots/H3ABox.pdf", width = 9)
ggplot(plot3A, aes(y = AitchisonDistance, fill = Group, x = factor(FreezeThaw))) +
  geom_boxplot() + scale_fill_viridis_d() +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  labs(x = "Number of Freeze Thaw Events (14 days prior to sampling)",
       y = "Aitchison Distance",
       fill = "Within Territory Comparison")
dev.off()

# PERMANOVA
adonis2(dmAitchison ~ FreezeThaw,
        data = metaWeather[which(metaWeather$sampleID %in% rownames(gj_meta)),],
        na.action = na.exclude)
print("Paired t-test")
tDist <- OTUsamples %>% filter(Territory == "Within") %>%
  mutate(Bonly = as.double(Bonly), NBonly = as.double(NBonly), shared = as.double(shared),
         diffShared = as.numeric(shared) - as.numeric(NBonly),
         diffUnique = as.numeric(Bonly) - as.numeric(NBonly)) %>%
  select(Bonly, NBonly, shared, diffUnique, diffShared) %>% na.omit()
# Step 0 - check assumptions (see Radziwill p 345-6)
# Step 1 - null and alternative hypotheses
print("H0 - equal mean distance (d = d0)")
print("HA - mean distance greater with more ft event than fewer OTUs (d > d0)")
# Step 2 - significance (going to 0.05)
# Step 3 - test statistic
summary(tDist)
sdA <- sd(tDist$diffShared)
meanA <- mean(tDist$diffShared)
meanN <- 0
n <- as.numeric(nrow(tDist))
t <- (meanA-meanN)/(sdA/sqrt(n))
# Step 4 - draw a figure (gonna do this later)
# Step 5 - find the p-value
print("p-value")
1 - pt(t, df=n-1)
#Step 6 - is p-val < a?
print("alpha")
qt(0.975, df=n-1)
#Step 7 - CI and R check
print("t-test")
t.test(tDist$shared, tDist$NBonly,
       paired = TRUE, alternative = "greater")
# clean up (double check items here)
rm(cacheGroup, eventCount, dates, gj_meta, weatherCombo, weather, dmAitchison,
   plot3A, gj_ps)

print("Prediction 3B - Shared microbiota")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3B-filtered-table.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

# extract taxonomy (was thinking I would collapse by this)
tax <- as(tax_table(gj_ps), "matrix") %>% as.data.frame %>%
  rownames_to_column(var = "OTU")
# scatter - x is freeze thaw events, y is % total microbiota shared
otu_df <- as(otu_table(gj_ps), "matrix") %>%
  as.data.frame %>% rownames_to_column(var = "sampleID") %>%
  pivot_longer(-sampleID, names_to = "OTU", values_to = "Count") %>%
  select(sampleID, OTU, Count) %>%
  .[which(.$Count > 0),] # select only those present

# combine samples by shared OTUs
plot3B <- full_join(otu_df, otu_df, by = "OTU") %>%
  .[-which(.$sampleID.x == .$sampleID.y),] %>%
  select(sampleID.x, sampleID.y) %>%
  group_by(sampleID.x, sampleID.y) %>% # group by samples
  rename(Sample1 = sampleID.x, Sample2 = sampleID.y) %>%
  count %>% # count the number of times the groups occur
  left_join(plot3A) %>% # add sample data for plotting
  na.omit() %>% ungroup # tidy data for plotting

pdf("CanadaJayMicrobiome/plots/H3B.pdf", width = 9)
ggplot(plot3B, aes(x = FreezeThaw, y = n)) +
  geom_point() + facet_grid(CollectionSeason~Year) +
  labs(x = "Number of Freeze Thaw Events (14 Days prior to sampling)",
       y = "Number of Shared Taxa")
dev.off()

# clean up - removing 3A+B data
rm(metaWeather, gj_ps, plot3A, plot3B,
   otu_df)

print("Prediction 3C + D - Data for Supplementation")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3CD-filtered-table.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

print("Prediction 3CD - Food supplementation")
# Read in 2017-2018 distance matrix
dmAitchisonB <- read_qza("P3CD-aitchison-distance.qza")$data 
# combine the aitchison distance data with metadata
dm_meta <- longDM(dmAitchisonB, "AitchisonDistance", gj_meta) %>%
  # add shared food information
  mutate(host = ifelse(.$JayID.x == .$JayID.y, "Same Bird", "Different Bird"),
         group = ifelse(.$`FoodSupplement.x` == "Y" & .$`FoodSupplement.y` == "Y", "Yes",
                        ifelse(.$`FoodSupplement.x` == "N" & .$`FoodSupplement.y` == "N",
                               "No", "Between")),
         Year = word(ExtractID.x, 2, sep = "-")) %>%
  # select only the variables of interest for the boxplot
  select(1:4, ExtractID.y, host, group, Year)

pdf("CanadaJayMicrobiome/plots/H3CDBox.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = group)) +
  geom_boxplot() + labs(x = "Food Supplementation") +
  facet_grid(~Year)
dev.off()

dm_meta[dm_meta$group=="No",] %>% .$AitchisonDistance
dm_meta[dm_meta$group=="Yes",] %>% .$AitchisonDistance
# Mann-Whitney test
print("Mann-Whitney")
wilcox.test(dm_meta[dm_meta$group=="No",] %>% .$AitchisonDistance,
            dm_meta[dm_meta$group=="Yes",] %>% .$AitchisonDistance)
print("t-test")
t.test(dm_meta[dm_meta$group=="No",] %>% .$AitchisonDistance,
       dm_meta[dm_meta$group=="Yes",] %>% .$AitchisonDistance)

