print("Figures")
# additional figures that don't necessarily fit into one hypothesis
# run locally (on equivalent local directory)
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load required packages
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
ordiAitchison <- read_qza("aitchison-ordination-cr.qza")
ordiAitchisonDN <- read_qza("aitchison-ordination-dn.qza")
# join aitchison vectors with metadata
# select columns of interest from meta (no location or extraction info)
gj_aitch_V <- gj_meta %>% select(1:5, 7:17, 24:27) %>%
  rownames_to_column(var = "SampleID") %>% # row names need to be a column to join
  left_join(ordiAitchison$data$Vectors,.) %>%
  remove_rownames()
gj_aitch_VDN <- gj_meta %>% select(1:5, 7:17, 24:27) %>%
  rownames_to_column(var = "SampleID") %>% # row names need to be a column to join
  left_join(ordiAitchisonDN$data$Vectors,.) %>%
  remove_rownames()

print("Counts")
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

print("Breeding status by age")
# breeding status by age plot
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/ageBreedingStatus.pdf")
ggplot(gj_meta, aes(x = BreedingStatus, y = as.numeric(AgeAtCollection))) +
  geom_boxplot() +
  labs(y = "Age At Collection", x = "Breeding Status")
dev.off()

# figure 2
print("Map of sampling locations")
# these locations account for ALL sample locations (2016-2020)
map_gj <- get_map(
  location = c(left = -79.5, bottom = 45.2, right = -77.9, top = 46.2),
  source = "osm", color = "bw",
  force = TRUE) # adding force = TRUE to get_map to force a re-rendering of the map

# add season/year column
gj_meta$SeasonYear <- paste(gj_meta$CollectionSeason, gj_meta$CollectionYear)
gj_meta2 <- gj_meta %>%
  count(LatitudeSamplingDD, LongitudeSamplingDD, SeasonYear) %>%
  inner_join(., gj_meta)
#filter for sample size >2
tt <- table(gj_meta2$SeasonYear)
gj_meta3 <- subset(gj_meta2, SeasonYear %in% names(tt[tt > 2]))
# export map (increase width/height)
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/mapSamples.pdf")
ggmap(map_gj) + #facet_grid(~SeasonYear) +
  geom_count(data = gj_meta3, 
             aes(y = LatitudeSamplingDD, x = LongitudeSamplingDD)) + #, color = SeasonYear
  #geom_text(hjust = 0, aes(label = 1))+
  theme(legend.position = "bottom", legend.box = "vertical") +
  #scale_color_viridis_d() +
  ggsn::scalebar(x.min = -79.4, x.max = -78, 
                 y.min = 45.3, y.max = 46.1,
                 location = "bottomleft",
                 dist = 10, dist_unit = "km",
                 transform = TRUE, model = "WGS84") +
  labs(shape = "Collection Year", size = "Number of Samples") +
  theme(text = element_text(size = 20)) +
  labs(y = "Lattitude", x = "Longitude")
dev.off()

# samples
samples <- ggmap(map_gj) + #facet_grid(~SeasonYear) +
  geom_jitter(data = gj_meta3, 
             aes(y = LatitudeSamplingDD, x = LongitudeSamplingDD)) +
  ggsn::scalebar(x.min = -79.4, x.max = -78, 
                 y.min = 45.3, y.max = 46.1,
                 location = "bottomleft",
                 dist = 10, dist_unit = "km",
                 transform = TRUE, model = "WGS84",
                 st.size = 2, border.size = 0.1)
# canada
map_ca <- get_map(
  location = c(left = -141, bottom = 42, right = -53, top = 83),
  source = "osm", color = "bw",
  force = TRUE)

Canada <- ggmap(map_ca) + 
  annotate("text", x = -78.3791, y = 45.8372, label = "o", colour = "black") +
  ggsn::scalebar(x.min = -139, x.max = -53.1, 
                 y.min = 44, y.max = 82.9,
                 location = "bottomleft",
                 dist = 1000, dist_unit = "km",
                 transform = TRUE, model = "WGS84",
                 st.size = 2, border.size = 0.1) +
  theme(text = element_text(size = 20)) +
  labs(y = "Lattitude", x = "Longitude")
# inset APP onto canada
library(grid)
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/mapInsetSamples.pdf")
Canada +
  inset(grob = ggplotGrob(samples + theme_inset()),
    xmin = -141, xmax = -90, ymin = 55, ymax = 86.9) +
  annotate("segment", x = -141, xend = -79, y = 74.4, yend = 45.5,
           colour = "black") +
  annotate("segment", x = -90, xend = -79, y = 74.4, yend = 46,
           colour = "black") + theme(text = element_text(size = 20)) 
dev.off()
# clean up
rm(map_gj)

# ordinations
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/PCA.pdf")
ggplot(gj_aitch_V, aes(PC1, PC2, shape = as.factor(CollectionYear))) +
  geom_point() + labs(shape = "Collection Year") +
  ggtitle("Closed Reference Clustering")
ggplot(gj_aitch_VDN, aes(PC1, PC2, shape = as.factor(CollectionYear))) +
  geom_point() + labs(shape = "Collection Year") +
  ggtitle("DADA2 De Novo Clustering")
dev.off()
