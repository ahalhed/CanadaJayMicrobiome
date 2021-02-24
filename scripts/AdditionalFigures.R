# additional figures that don't necessarily fit into one hypothesis
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
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

# save plot
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

## PCA plot (closed reference)
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/PCA-cr.pdf", width = 10)
# need to sort out NA thing
ggplot(gj_aitch_V, aes(y=PC2, x=PC1, shape = as.factor(CollectionYear), group = JayID)) + #, group = JayID
  geom_point() + #geom_line()
  labs(shape = "Collection Year")
dev.off()
## PCA plot (de novo) - need to do this