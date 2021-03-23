# additional figures that don't necessarily fit into one hypothesis
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggmap)
library(phyloseq)
library(vegan)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza",
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
# read in the aitchison ordination (need to generate this ordi)
ordiAitchison <- read_qza("aitchison-ordination.qza")
# join aitchison vectors with metadata
# select columns of interest from meta (no location or extraction info)
gj_aitch_V <- gj_meta %>% select(1:5, 7:17, 24:27) %>%
  rownames_to_column(var = "SampleID") %>% # row names need to be a column to join
  left_join(ordiAitchison$data$Vectors,.) %>%
  remove_rownames()
# read in core data
coreTable <- read.csv("CanadaJayMicrobiome/data/coreJay.csv")

# save core plot
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/coreSites.pdf")
ggplot(coreTable, aes(y = otu_occ, x = otu_rel, color = fill)) + 
  geom_point() +
  # log transform the x axis, set discrete viridis colour scheme
  scale_x_log10() + scale_colour_viridis_d() + 
  # add axis labels
  labs(x = "Mean Relative Abundance of Each OTU (log10)", 
       y = "Occupancy (Proportion of Samples)",
       color = "OTU Type")
dev.off()
# clean up
rm(coreTable)

# counts
# plot b/nb counts
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/countsSampleType.pdf")
ggplot(gj_meta, aes(x = interaction(CollectionSeason, CollectionYear, sep = " "),
                fill = interaction(BreedingStatus, JuvenileStatus))) +
  geom_bar() +
  # reorder x axis and angle of labels
  scale_x_discrete(limits = c("Fall 2016", "Fall 2017", "Fall 2018", "Winter 2019", 
                              "Spring 2019", "Winter 2020", "Spring 2020", "Fall 2020")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  # changing the interaction legend
  scale_fill_viridis_d(labels = c("Breeder (Dominant Juvenile)", "Non-Breeder (Dominant Juvenile)", "Breeder (Ejectee)",
                                  "Breeder (Unknown)", "Non-breeder (Unknown)", "Unknown")) +
  # adjust axis and legend names
  labs(x = "Collection Season and Year", y = "Number of Individuals", 
       fill = "Breeding Status (Juvenile Status)") +
  ggtitle("Seasonal Distribution of Collected Canada Jay Oral Microbiome Samples",
          "All Individuals by Breeding and Juvenile Status")
gj_meta %>% filter(BreedingStatus == "Breeder") %>%
  ggplot(aes(x = interaction(CollectionSeason, CollectionYear, sep = " "),
             fill = Territory)) + geom_bar() +
  scale_fill_viridis_d() +
  # reorder x axis labels
  scale_x_discrete(limits = c("Fall 2016", "Fall 2017", "Fall 2018", 
                              "Winter 2020", "Spring 2020", "Fall 2020")) +
  # adjust axis and legend names
  labs(x = "Collection Season and Year", y = "Number of Individuals") +
  ggtitle("Seasonal Distribution of Collected Canada Jay Oral Microbiome Samples",
          "Breeding Individuals Only by Territory")
gj_meta %>% filter(BreedingStatus == "Non-breeder") %>%
  ggplot(aes(x = interaction(CollectionSeason, CollectionYear, sep = " "),
             fill = Territory)) + geom_bar() +
  scale_fill_viridis_d() +
  # reorder x axis and angle of labels
  scale_x_discrete(limits = c("Fall 2017", "Fall 2018", "Winter 2019", "Spring 2019", "Spring 2020", "Fall 2020")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  facet_grid(~JuvenileStatus) +
  # adjust axis and legend names
  labs(x = "Collection Season and Year", y = "Number of Individuals") +
  ggtitle("Seasonal Distribution of Collected Canada Jay Oral Microbiome Samples",
          "Non-breeders Only by Territory")
dev.off()

# breeding status by age plot
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/ageBreedingStatus.pdf")
ggplot(jay, aes(x = BreedingStatus, y = as.numeric(AgeAtCollection))) +
  geom_boxplot() +
  labs(y = "Age At Collection", x = "Breeding Status")
dev.off()

# figure 2
# map of sampling locations
# these locations account for ALL sample locations (2016-2020)
map_gj <- get_map(
  location = c(left = -79.5, bottom = 45.2, right = -77.9, top = 46.2),
  source = "osm",
  force = TRUE) # adding force = TRUE to get_map to force a re-rendering of the map
# export map
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/mapSamples.pdf")
ggmap(map_gj) + 
  geom_count(data = gj_meta, 
             aes(y = LatitudeSamplingDD, x = LongitudeSamplingDD, 
                 shape = as.character(CollectionYear))) + 
  theme(legend.position = "bottom", legend.box = "vertical") +
  labs(shape = "Collection Year", size = "Number of Samples",
       title = "Map of Canada Jay Sampling Locations",
       subtitle = "Algonquin Park, Ontario (2016-2020)")
dev.off()


