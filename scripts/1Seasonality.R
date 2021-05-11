# for Hypothesis 1 - related to seasonal differences

print("Working directory, packages, functions")
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(qiime2R)
library(phyloseq)
library(ALDEx2)
library(tidyverse)
# set theme 
theme_set(theme_bw())
aldAn <- function(met, year, group, lev, v = FALSE) {
  # met is the sample metadata
  # year is the year(s) of interest
  # group is the column to select
  # lev is a level(s) of group to select (optional)
  # v is whether or not to view aldex output
  df1 <- met[which(met$CollectionYear %in% year),]
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
  # volcano plots
  if(missing(lev)){
    pa <- paste0("CanadaJayMicrobiome/plots/volcanoes/", year, group, ".pdf")
  } else {
    pa <- paste0("CanadaJayMicrobiome/plots/volcanoes/", year, group, lev, ".pdf")
  }
  pdf(pa)
  par(mfrow=c(1,2))
  aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
             ylab="Difference")
  aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
             ylab="Difference")
  dev.off()
  # clean up
  rm(Ald, conds, orderPathAbun)
}


print("prediction 1A - by season?")
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
# confirmation that year/season was the right way to go
print("Collection Year Ordination")
(orY <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionYear))
anova(orY)
print("Collection Season Ordination")
(orS <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionSeason))
anova(orS)
print("Collection Year*Season Ordination")
(orSY <- ordinate(gj_ps, method = "RDA", formula = . ~ CollectionYear * CollectionSeason))
anova(orSY)
# clean up
rm(gj_meta, gj_ps, orS, orSY, orY)

print("prediction 1B - functional?")
print("Core figure")
# read in core data
coreTable <- read.csv("CanadaJayMicrobiome/data/coreJay.csv")
coreTable$samples <- round(coreTable$otu_occ*88)

corePlot <- ggplot(coreTable, aes(y = otu_occ, x = otu_rel, color = fill)) + 
  geom_point() +
  # log transform the x axis, set discrete viridis colour scheme
  scale_x_log10() + scale_colour_viridis_d() + 
  # add axis labels
  labs(x = "Mean Relative Abundance of Each OTU (log10)", 
       y = "Occupancy (Proportion of Samples)",
       color = "OTU Type")
# save core plot
pdf("CanadaJayMicrobiome/plots/P1B.pdf")
corePlot    # without labels
corePlot +  # with text labels
  geom_text(data=coreTable[which(coreTable$fill == "Core"),],
            aes(y = otu_occ, x = otu_rel, label= Genus),
            color='black', size=2.5,
            position=position_jitter(width=0.01,height=0.01))
ggplot(coreTable, aes(y = FallOcc, x = FallRel, color = fill)) + 
  geom_point() + scale_colour_viridis_d() + 
  # add axis labels
  labs(x = "Mean Relative Abundance of Each OTU", 
       y = "Occupancy (Proportion of Fall Samples)",
       color = "OTU Type")
ggplot(coreTable, aes(y = WSOcc, x = WSRel, color = fill)) + 
  geom_point() + scale_colour_viridis_d() + 
  # add axis labels
  labs(x = "Mean Relative Abundance of Each OTU", 
       y = "Occupancy (Proportion of Winter/Spring Samples)",
       color = "OTU Type")
dev.off()
# clean up
rm(coreTable, corePlot)

print("differential abundance")
meta <- read_q2metadata("input/jay-met.tsv") %>%
  .[which(.$JayID != "BLANK"),] %>%
  remove_rownames() %>% column_to_rownames(var = "SampleID")
pathAbun <- read_qza("q2-picrust2_output/pathway_abundance.qza")$data
pathAbunI <- pathAbun %>% as.data.frame %>% mutate_all(~as.integer(.))

#can only have 2 levels (fall)
fall <- meta[which(meta$CollectionSeason == "Fall"),]
# year - 2017 vs 2018 
aldAn(fall, c(2017, 2018), "CollectionYear")

# year - 2017 vs 2020 
aldAn(fall, c(2020, 2017), "CollectionYear")

# year - 2018 vs 2020 (fall)
aldAn(fall,  c(2018, 2020), "CollectionYear")

# Season
springFall <- meta[which(meta$CollectionSeason %in% c("Spring", "Fall")),]
conds <- springFall %>%
  select(CollectionSeason) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
             include.sample.summary=FALSE, denom="all", verbose=FALSE)

pdf("CanadaJayMicrobiome/plots/volcanoes/seasons.pdf")
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(Ald, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")
dev.off()
# clean up
rm(Ald, conds, orderPathAbun, springFall)