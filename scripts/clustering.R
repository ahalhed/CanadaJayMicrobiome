setwd("~/OneDrive - University of Guelph/Alicia's Thesis/GreyJayMicrobiome")
library(qiime2R)
library(tidyverse)
# read in the data
dist <- read_qza("aitchison-distance-cr.qza")
met <- read_q2metadata("input/jay-met.tsv") %>% .[which(.$JayID != "BLANK"),]
# clustering
hc <- hclust(dist$data)
plot(hc)

# function for relabelling
labClust <- function(d, m, c) {
  # d is a distance matrix
  # m is the metadata df
  # c is the metadata column
  nM <- m[order(m$SampleID),] %>% select(c) 
  named <- usedist::dist_setNames(d, nM[,1])
  hc <- hclust(named)
  plot(hc)
}

labClust(dist$data, met, 'Sex')
labClust(dist$data, met, 'AgeAtCollection')
labClust(dist$data, met, 'BirthYear')
labClust(dist$data, met, 'TerritoryQuality')
labClust(dist$data, met, 'Territory')
labClust(dist$data, met, 'CollectionYear')
