# for Hypothesis 3 - feeding differences
# this script is for the diet supplementation plot
# to load R on interactive graham
# module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/GreyJayMicrobiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
library(tidyverse)

theme_set(theme_bw())

print("Prediction 3 - Food supplementation (OTUs)")
## Load in the required data
# build the phyloseq object
gj_ps <- qza_to_phyloseq(features = "3-filtered-table.qza",
                         # q2 types line causes issues (so removed in the tsv file input here)
                         metadata = "input/jay-met.tsv") %>%
  # transposing the OTU table into the format expected by vegan (OTUs as columns)
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# extract the metadata from the phyloseq object
gj_meta <- as(sample_data(gj_ps), "data.frame")
rownames(gj_meta) <- sample_names(gj_ps)
gj_meta <- gj_meta %>% rownames_to_column(var = "sampleID")
# scatter - x is freeze thaw events, y is % total microbiota shared
otu_df <- as(otu_table(gj_ps), "matrix") %>%
  as.data.frame %>% rownames_to_column(var = "sampleID") %>%
  pivot_longer(-sampleID, names_to = "OTU", values_to = "Count") %>%
  select(sampleID, OTU, Count) %>%
  .[which(.$Count > 0),] # select only those present

# summarize to number of OTUs (not counts per)
plot3 <- otu_df %>% mutate(Count = 1) %>%
  select(sampleID, Count) %>%
  group_by(sampleID) %>%
  count %>% ungroup %>% # count the number of OTUs per group
  full_join(gj_meta, by = "sampleID") %>% # add sample data for plotting
  filter(FoodSupplement != "U")

pdf("CanadaJayMicrobiome/plots/H3.pdf", width = 9)
ggplot(plot3, aes(x = FoodSupplement, y = n)) +
  geom_boxplot() + facet_grid(~CollectionYear) +
  geom_dotplot(binaxis = "y", binwidth = 1, stackdir = "center", fill = NA) +
  labs(x = "Food Supplementation Group",
       y = "Number of OTUs") + theme(text = element_text(size = 20)) 
dev.off()

print("two sample t-test")
# step 0 - check assumptions
# step 1 - set null and alternative hypotheses
print("H0 - Means are not different (u1-u2=D0)")
print("HA - yes FS mean is less than no FS mean (u1-u2>D0)")
# step 2 - choose alpha (0.05)
# step 3 - calculate test statistic
print("all samples")
yes <- plot3C %>% filter(FoodSupplement == "Y") %>% select(n)
no <- plot3C %>% filter(FoodSupplement == "N") %>% select(n)
var.test(yes$n, no$n)
print("summary of yes")
summary(yes)
print("summary of no")
summary(no)
# step 6 - calculate with R
t.test(yes, no, alternative = "less")
# clean up
rm(yes, no)
print("2017 only")
yes <- plot3C %>% filter(FoodSupplement == "Y", CollectionYear == 2017) %>% select(n)
no <- plot3C %>% filter(FoodSupplement == "N", CollectionYear == 2017) %>% select(n)
var.test(yes$n, no$n)
t.test(yes, no, alternative = "less")
rm(yes, no)
print("2018 only")
yes <- plot3C %>% filter(FoodSupplement == "Y", CollectionYear == 2018) %>% select(n)
no <- plot3C %>% filter(FoodSupplement == "N", CollectionYear == 2018) %>% select(n)
var.test(yes$n, no$n)
t.test(yes, no, alternative = "less")
rm(yes, no, plot3C)

# full clean up
rm(dates, gj_meta, gj_ps, otu_df, samp)
