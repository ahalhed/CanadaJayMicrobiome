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

print("Prediction 3A - Weather data for freeze thaw")
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
         TerritoryQuality, BreedingStatus, Territory) %>%
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

print("Prediction 3A - FT number of microbiota")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3A-filtered-table.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta <- rownames_to_column(gj_meta, var = "sampleID")

# scatter - x is freeze thaw events, y is % total microbiota shared
otu_df <- as(otu_table(gj_ps), "matrix") %>%
  as.data.frame %>% rownames_to_column(var = "sampleID") %>%
  pivot_longer(-sampleID, names_to = "OTU", values_to = "Count") %>%
  select(sampleID, OTU, Count) %>%
  .[which(.$Count > 0),] # select only those present

# summarize to number of OTUs (not counts per)
plot3A <- otu_df %>% mutate(Count = 1) %>%
  select(sampleID, Count) %>%
  group_by(sampleID) %>%
  count %>% ungroup %>% # count the number of times the groups occur
  full_join(gj_meta, by = "sampleID") %>% # add sample data for plotting
  left_join(metaWeather %>% select(sampleID, FreezeThaw))

pdf("CanadaJayMicrobiome/plots/P3A.pdf", width = 9)
ggplot(plot3A, aes(x = FreezeThaw, y = n)) +
  geom_jitter() + geom_abline() +
  labs(x = "Number of Freeze Thaw Events (14 Days prior to sampling)",
       y = "Number of OTUs")
dev.off()

# linear model
(lm3A <- lm(n~FreezeThaw+BreedingStatus, data = plot3A))
summary(lm3A)
anova(lm3A)
# clean up - removing 3A+B data
rm(metaWeather, gj_ps, plot3A,
   otu_df, lm3A, gj_meta)

print("Prediction 3B+C - Data for Supplementation")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3BC-filtered-table.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta <- gj_meta %>% rownames_to_column(var = "sampleID")

print("Prediction 3B - Food supplementation (OTUs)")
# scatter - x is freeze thaw events, y is % total microbiota shared
otu_df <- as(otu_table(gj_ps), "matrix") %>%
  as.data.frame %>% rownames_to_column(var = "sampleID") %>%
  pivot_longer(-sampleID, names_to = "OTU", values_to = "Count") %>%
  select(sampleID, OTU, Count) %>%
  .[which(.$Count > 0),] # select only those present

# summarize to number of OTUs (not counts per)
plot3B <- otu_df %>% mutate(Count = 1) %>%
  select(sampleID, Count) %>%
  group_by(sampleID) %>%
  count %>% ungroup %>% # count the number of OTUs per group
  full_join(gj_meta, by = "sampleID") %>% # add sample data for plotting
  filter(FoodSupplement != "U")

pdf("CanadaJayMicrobiome/plots/P3B.pdf", width = 9)
ggplot(plot3B, aes(x = FoodSupplement, y = n)) +
  geom_boxplot() + facet_grid(~CollectionYear) +
  labs(x = "Food Supplementation Group",
       y = "Number of OTUs")
dev.off()

print("one sample t-test")
# step 0 - check assumptions
# step 1 - set null and alternative hypotheses
print("H0 - Means are not different (u=0)")
print("HA - yes FS mean is less than no FS mean (u>u0)")
# step 2 - choose alpha (0.05)
# step 3 - calculate test statistic
print("all samples")
print("test statistic")
yes <- plot3B %>% filter(FoodSupplement == "Y") %>% select(n)
no <- plot3B %>% filter(FoodSupplement == "N") %>% select(n)
print("summary of yes")
summary(yes)
print("summary of no")
summary(no)
sdA <- sd(yes$n)
meanA <- mean(yes$n)
meanN <- mean(no$n)
n <- as.numeric(nrow(yes))
t <- (meanA-meanN)/(sdA/sqrt(n))
# Step 4 - draw a figure (gonna do this later)
# Step 5 - find the p-value
print("p-value")
pt(t,df=n-1)
# step 6 - calculate with R
t.test(yes, mu=meanN, alternative = "less")
# clean up
rm(yes, no, meanA, meanN, n, sdA, t)
print("2017 only")
yes <- plot3B %>% filter(FoodSupplement == "Y", CollectionYear == 2017) %>% select(n)
no <- plot3B %>% filter(FoodSupplement == "N", CollectionYear == 2017) %>% select(n)
meanN <- mean(no$n)
t.test(yes, mu=meanN, alternative = "less")
rm(yes, no, meanN)
print("2018 only")
yes <- plot3B %>% filter(FoodSupplement == "Y", CollectionYear == 2018) %>% select(n)
no <- plot3B %>% filter(FoodSupplement == "N", CollectionYear == 2018) %>% select(n)
meanN <- mean(no$n)
t.test(yes, mu=meanN, alternative = "less")
rm(yes, no, meanN,plot3B)

print("Prediction 3C - Food supplementation (distances)")
# Read in data for 3C
dmAitchisonB <- read_qza("P3C-aitchison-distance.qza")$data 
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# combine the aitchison distance data with metadata
dm_meta <- longDM(dmAitchisonB, "AitchisonDistance", gj_meta) %>%
  filter(CollectionYear.x == CollectionYear.y) %>%
  # add shared food information
  mutate(host = ifelse(.$JayID.x == .$JayID.y, "Same Bird", "Different Bird"),
         group = ifelse(.$`FoodSupplement.x` == "Y" & .$`FoodSupplement.y` == "Y", "Yes",
                        ifelse(.$`FoodSupplement.x` == "N" & .$`FoodSupplement.y` == "N",
                               "No", "Between")))

pdf("CanadaJayMicrobiome/plots/P3D.pdf", width = 9)
ggplot(dm_meta, aes(y = AitchisonDistance, x = group)) +
  geom_boxplot() + labs(x = "Food Supplementation", y = "Aitchison Distance") +
  facet_grid(~CollectionYear.x)
dev.off()

# ANOSIM
anoMet <- gj_meta[which(gj_meta$BreedingStatus == "Breeder"),]
print("Both years")
SampN <- anoMet %>% rownames
dmSamp <- dmAitchisonB %>% as.matrix %>%
  .[which(rownames(.) %in% SampN), which(colnames(.) %in% SampN)] %>%
  as.dist
with(anoMet, anosim(dmSamp, FoodSupplement)) %>% summary
rm(SampN, dmSamp)

print("2017 Only")
SampD <- anoMet[which(anoMet$CollectionYear == 2017),] 
SampN <- SampD %>% rownames
dmSamp <- dmAitchisonB %>% as.matrix %>%
  .[which(rownames(.) %in% SampN), which(colnames(.) %in% SampN)] %>%
  as.dist
with(SampD, anosim(dmSamp, FoodSupplement)) %>% summary()
rm(SampD, SampN, dmSamp)

print("2018 Only")
SampD <- anoMet[which(anoMet$CollectionYear == 2018),] 
SampN <- SampD %>% rownames
dmSamp <- dmAitchisonB %>% as.matrix %>%
  .[which(rownames(.) %in% SampN), which(colnames(.) %in% SampN)] %>%
  as.dist
with(SampD, anosim(dmSamp, FoodSupplement)) %>% summary()
rm(SampD, SampN, dmSamp)

# clean up
rm(dm_meta, gj_meta, gj_ps, otu_df, dmAitchisonB, longDM, anoMet)
