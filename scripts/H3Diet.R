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

# ANOSIM repititions
anoRep <- function(samp, breed, dis, year=NULL) {
  # samp is df with sample data
  # breed is the breeding status of interest
  # year is the collection year of interest
  # dis is a distance matrix
  if(missing(year)) {
    print(breed)
    SampD <- samp[which(samp$BreedingStatus == breed),]
    SampN <- SampD %>% rownames
    dmSamp <- dis %>% as.matrix %>%
      .[which(rownames(.) %in% SampN), which(colnames(.) %in% SampN)] %>%
      as.dist
    ano1A <- with(SampD, anosim(dmSamp, Territory))
    summary(ano1A)
  } else{
    print(year)
    SampD <- samp[which(samp$CollectionYear == year),]
    SampN <- SampD %>% rownames
    dmSamp <- dis %>% as.matrix %>%
      .[which(rownames(.) %in% SampN), which(colnames(.) %in% SampN)] %>%
      as.dist
    ano1A <- with(SampD, anosim(dmSamp, Territory))
    summary(ano1A)
  }
  rm(SampD, SampN, dmSamp, ano1A)
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
  left_join(metaWeather %>% select(sampleID, BreedingStatus, Territory),
          by = c("Sample1" = "sampleID")) %>%
  left_join(metaWeather %>% select(sampleID, BreedingStatus, Territory),
            by = c("Sample2" = "sampleID"), suffix = c("1", "2"))
plot3A$Group <- ifelse(plot3A$BreedingStatus1 == plot3A$BreedingStatus2 &
                         plot3A$BreedingStatus1 == "Breeder",
                       "Within Breeders",
                       ifelse(plot3A$BreedingStatus1 == plot3A$BreedingStatus2 &
                           plot3A$BreedingStatus1 == "Non-breeder",
                         "Within Non-Breeders", "Between Breeder & Non-Breeder"))
plot3A$Territory <- ifelse(plot3A$Territory1 == plot3A$Territory2,
                           "Same Territory", "Different Territory")
plot3A <- plot3A %>% select(-c(Territory1, Territory2, BreedingStatus1, BreedingStatus2))

# make figure
pdf("CanadaJayMicrobiome/plots/P3A.pdf", width = 9)
ggplot(plot3A, aes(y = AitchisonDistance, fill = Group, x = factor(FreezeThaw))) +
  geom_boxplot() + scale_fill_viridis_d() +
  facet_grid(~Territory) +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  labs(x = "Number of Freeze Thaw Events (14 days prior to sampling)",
       y = "Aitchison Distance",
       fill = "Within Territory Comparison")
dev.off()

# ANOSIM
anoA <- with(gj_meta, anosim(dmAitchison, interaction(Territory, BreedingStatus)))
summary(anoA)
anoRep(gj_meta, "Breeder", dmAitchison)
anoRep(gj_meta, "Non-breeder", dmAitchison)
# clean up (double check items here)
rm(dates, gj_meta, dmAitchison, samp, gj_ps, anoA)

print("Prediction 3B - FT number of microbiota")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "P3B-filtered-table.qza",
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
plot3B <- otu_df %>% mutate(Count = 1) %>%
  select(sampleID, Count) %>%
  group_by(sampleID) %>%
  count %>% ungroup %>% # count the number of times the groups occur
  full_join(gj_meta, by = "sampleID") %>% # add sample data for plotting
  right_join(plot3A %>% select(Sample1, FreezeThaw),
            by = c("sampleID" = "Sample1"))

pdf("CanadaJayMicrobiome/plots/P3B.pdf", width = 9)
ggplot(plot3B, aes(x = FreezeThaw, y = n)) +
  geom_jitter() + #geom_abline() +
  labs(x = "Number of Freeze Thaw Events (14 Days prior to sampling)",
       y = "Number of OTUs")
dev.off()

# linear model
(lm3B <- lm(n~FreezeThaw, data = plot3B))
summary(lm3B)
anova(lm3B)
# clean up - removing 3A+B data
rm(metaWeather, gj_ps, plot3A, plot3B,
   otu_df, lm3B, gj_meta)

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
gj_meta <- gj_meta %>% rownames_to_column(var = "sampleID")

print("Prediction 3C - Food supplementation (OTUs)")
# scatter - x is freeze thaw events, y is % total microbiota shared
otu_df <- as(otu_table(gj_ps), "matrix") %>%
  as.data.frame %>% rownames_to_column(var = "sampleID") %>%
  pivot_longer(-sampleID, names_to = "OTU", values_to = "Count") %>%
  select(sampleID, OTU, Count) %>%
  .[which(.$Count > 0),] # select only those present

# summarize to number of OTUs (not counts per)
plot3C <- otu_df %>% mutate(Count = 1) %>%
  select(sampleID, Count) %>%
  group_by(sampleID) %>%
  count %>% ungroup %>% # count the number of OTUs per group
  full_join(gj_meta, by = "sampleID") %>% # add sample data for plotting
  filter(FoodSupplement != "U")

pdf("CanadaJayMicrobiome/plots/P3C.pdf", width = 9)
ggplot(plot3C, aes(x = FoodSupplement, y = n)) +
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
print("all")
print("test statistic")
yes <- plot3C %>% filter(FoodSupplement == "Y") %>% select(n)
no <- plot3C %>% filter(FoodSupplement == "N") %>% select(n)
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
yes <- plot3C %>% filter(FoodSupplement == "Y", CollectionYear == 2017) %>% select(n)
no <- plot3C %>% filter(FoodSupplement == "N", CollectionYear == 2017) %>% select(n)
meanN <- mean(no$n)
t.test(yes, mu=meanN, alternative = "less")
rm(yes, no, meanN)
print("2018 only")
yes <- plot3C %>% filter(FoodSupplement == "Y", CollectionYear == 2018) %>% select(n)
no <- plot3C %>% filter(FoodSupplement == "N", CollectionYear == 2018) %>% select(n)
meanN <- mean(no$n)
t.test(yes, mu=meanN, alternative = "less")
rm(yes, no, meanN,plot3C)

print("Prediction 3D - Food supplementation (distances)")
# Read in data for 3D
dmAitchisonB <- read_qza("P3D-aitchison-distance.qza")$data 
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
anoA <- with(gj_meta, anosim(dmAitchisonB, FoodSupplement))
summary(anoA)
anoRep(gj_meta, "Breeder", dmAitchisonB, 2017)
anoRep(gj_meta, "Breeder", dmAitchisonB, 2018)

# clean up
rm(dm_meta, gj_meta, gj_ps, otu_df, dmAitchisonB, longDM, anoA, anoRep)
rm(longWeather, weather)