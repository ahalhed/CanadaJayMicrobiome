# for Hypothesis 1 - related to seasonal differences

print("Working directory, packages, functions")
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
# load packages
library(qiime2R)
library(phyloseq)
library(zCompositions)
library(CoDaSeq)
library(ALDEx2)
library(tidyverse)
# set theme 
theme_set(theme_bw())

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

# impute the OTU table
OTUclr <- cmultRepl(otu_table(gj_ps), label=0, method="CZM") %>% # all OTUs
  codaSeq.clr # compute the CLR values
# update phyloseq object
up_ps <- phyloseq(otu_table(OTUclr, taxa_are_rows = F), sample_data(gj_ps), tax_table(gj_ps))
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
rm(gj_meta, gj_ps, up_ps, orS, orSY, orY)

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
pdf("CanadaJayMicrobiome/plots/AdditionalFigures/core.pdf")
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

print("prediction 1B - functional?")
print("differential abundance")
meta <- read_q2metadata("input/jay-met.tsv") %>%
  .[which(.$JayID != "BLANK"),] %>%
  remove_rownames() %>% column_to_rownames(var = "SampleID")
pathAbun <- read_qza("1B-picrust2_output/pathway_abundance.qza")$data
pathAbunI <- pathAbun %>% as.data.frame %>% mutate_all(~as.integer(.))

# Season (main analysis)
springFall <- meta[which(meta$CollectionSeason %in% c("Spring", "Fall")),]
conds <- springFall %>%
  select(CollectionSeason) %>% t %>% as.data.frame %>%
  .[ , order(names(.))] %>% unlist
orderPathAbun <- pathAbunI %>% select(names(conds))
# log-ratio transformation and statistical testing 
Ald <- aldex(orderPathAbun, conds, test="t", effect=TRUE,
             include.sample.summary=FALSE, denom="all", verbose=FALSE)
# save results to file
write.csv(Ald, file = "CanadaJayMicrobiome/data/diffAbun.csv")
# check the "explaining the outputs" section of the vignette
# https://rdrr.io/bioc/ALDEx2/f/vignettes/ALDEx2_vignette.Rmd

pdf("CanadaJayMicrobiome/plots/1B.pdf")
# create aldex plot
par(mfrow=c(1,2))
aldex.plot(Ald, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", cutoff.pval=0.05)
aldex.plot(Ald, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", cutoff.pval=0.05)
dev.off()

# save significantly associated pathways to file
# we.eBH is the Expected Benjamini-Hochberg corrected P value of Welch's t test
# associated with spring (positive diff)
spring <- Ald[which(Ald$we.eBH<0.05 & Ald$effect>0),]
write.csv(spring, file = "CanadaJayMicrobiome/data/diffAbunSpring.csv")
# associated with fall
fall <- Ald[which(Ald$we.eBH<0.05 & Ald$effect<0),]
write.csv(fall, file = "CanadaJayMicrobiome/data/diffAbunFall.csv")
