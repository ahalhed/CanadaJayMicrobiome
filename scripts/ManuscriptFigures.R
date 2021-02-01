# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")

# addition plots scripts
# not associated with a particular hypothesis
library(phyloseq)
library(qiime2R)
library(ggmap)
library(tidyverse)

theme_set(theme_bw())

## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "filtered-table-no-blanks.qza", 
                         tree = "trees/rooted-tree.qza", 
                         taxonomy = "taxonomy/SILVA-taxonomy.qza",
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.), tax_table(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)

# figure 1
# created in PPT (jay hyp/pred)

# figure 2
# map of sampling locations
# plot the locations (test run)
ggplot(gj_meta, aes(y = LatitudeSamplingDD, x = LongitudeSamplingDD,
                 shape=as.character(CollectionYear))) +
  geom_count() +  #facet_grid(CollectionYear~.) +
  labs(y = "Latitude", x = "Longitude", shape = "Collection Year", size = "Number of Samples")


# these locations account for ALL sample locations (2016-2020)
map_gj <- get_map(
  location = c(left = -79.5, bottom = 45.2, right = -77.9, top = 46.2),
  source = "osm",
  force = TRUE) # adding force = TRUE to get_map to force a re-rendering of the map
# Let's look at our map
ggmap(map_gj)


# Let's just throw everything onto the map to see what we've got
ggmap(map_gj) + 
  # switch count to jitter to get around being right on top of one another
  geom_count(data = gj_meta, 
             aes(y = LatitudeSamplingDD, x = LongitudeSamplingDD, 
                 shape = as.character(CollectionYear))) + 
  #facet_grid(CollectionYear~.) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  labs(shape = "Collection Year", size = "Number of Samples",
       title = "Map of Canada Jay Sampling Locations",
       subtitle = "Algonquin Park, Ontario (2016-2020)")

